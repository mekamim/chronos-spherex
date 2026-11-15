"""
FastAPI Server for SPHEREx Time-Domain Sky Explorer & Planet X Discovery Platform.
NASA Space Apps Challenge 2026: "Planet X and SPHEREx"
"""

from __future__ import annotations

import io
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.io import fits
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

from .spherex_engine import (
    SPHEREX_FLAG_BITS,
    extract_spectrum_and_lightcurve,
    get_real_or_cached_irsa_cutout,
    query_irsa_sia_exposures,
    build_multi_epoch_cube,
)
from .targets import (
    CHALLENGE_TARGETS,
    SPECTRAL_FEATURES,
    SPHEREX_BANDS,
    parse_sky_coordinates,
)

BASE_DIR = Path(__file__).resolve().parent.parent
FRONTEND_DIR = BASE_DIR / "frontend"

app = FastAPI(
    title="CHRONOS // SPHEREx — Spatial Infrared Observatory",
    description=(
        "Interactive Multi-Epoch Near-Infrared Sky Blinker, Difference Imager, and "
        "102-Band Spectrophotometer for NASA SPHEREx (0.75-5.0 µm) — https://chronosspherex.earth"
    ),
    version="2.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.on_event("startup")
def precache_mandatory_targets() -> None:
    """
    Ensure the 5 mandatory NASA Space Apps Challenge targets have their live IRSA SIAv2
    exposure metadata cached on disk for instant loading during live demos.
    """
    for target in CHALLENGE_TARGETS:
        try:
            query_irsa_sia_exposures(target["ra"], target["dec"], timeout=3.5)
        except Exception:
            pass


@app.get("/api/health")
def health_check():
    return {
        "status": "ok",
        "mission": "NASA SPHEREx (Spectro-Photometer for the History of the Universe, Epoch of Reionization and Ices Explorer)",
        "challenge": "NASA Space Apps 2026 - Planet X and SPHEREx",
        "bands": len(SPHEREX_BANDS),
        "spectral_channels": 102,
        "wavelength_range_um": [0.75, 5.00],
    }


@app.get("/api/targets")
def get_reference_targets():
    """
    Return the 5 mandatory NASA Space Apps test targets, the 6 SPHEREx detector band specifications,
    the 102-band spectral molecular features, and the 32-bit FLAGS bitmask definitions.
    """
    return {
        "targets": CHALLENGE_TARGETS,
        "bands": SPHEREX_BANDS,
        "spectral_features": SPECTRAL_FEATURES,
        "flag_definitions": SPHEREX_FLAG_BITS,
    }


@app.get("/api/resolve")
def resolve_coordinates(
    q: str = Query(..., description="Target name or sky coordinates (e.g. 'Barnard\\'s Star', '269.452, 4.693', '17h57m48.5s +04d41m36s')")
):
    """Parse flexible coordinate strings or target names into decimal J2000 (RA, Dec)."""
    try:
        ra, dec, matched = parse_sky_coordinates(q)
        return {
            "query": q,
            "ra": round(ra, 6),
            "dec": round(dec, 6),
            "matched_target": matched,
        }
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))


@app.get("/api/sky")
def get_multi_epoch_sky(
    coord: Optional[str] = Query(None, description="Unified coordinate or target name string"),
    ra: Optional[float] = Query(None, description="Right Ascension in decimal degrees (0..360)"),
    dec: Optional[float] = Query(None, description="Declination in decimal degrees (-90..+90)"),
    fov: float = Query(4.5, ge=1.0, le=20.0, description="Field of View in arcminutes"),
    grid_size: int = Query(96, ge=48, le=128, description="Spatial grid resolution in pixels"),
    zodi_subtract: bool = Query(True, description="Subtract Zodiacal Light (ZODI) continuum model"),
    motion_exaggeration: float = Query(1.0, ge=0.5, le=10.0, description="Visual proper motion baseline multiplier"),
    include_planetx: bool = Query(False, description="Optional: Inject benchmark moving test source in FOV"),
    fetch_live_cutout: bool = Query(False, description="Fetch sample live FITS cutout header/stats from IRSA IBE"),
):
    """
    Primary endpoint returning the multi-epoch SPHEREx sky cube across all 6 detector arrays (D1-D6),
    broadband coadds, Epoch B - Epoch A difference images, FLAGS bitmasks, Gaia DR3 + WISE overlays,
    and automated dipole/monopole anomaly detections.
    """
    if coord:
        try:
            ra_val, dec_val, _ = parse_sky_coordinates(coord)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc))
    elif ra is not None and dec is not None:
        ra_val, dec_val = float(ra), float(dec)
    else:
        # Default to Barnard's Star
        ra_val, dec_val = CHALLENGE_TARGETS[0]["ra"], CHALLENGE_TARGETS[0]["dec"]

    cube = build_multi_epoch_cube(
        ra=ra_val,
        dec=dec_val,
        fov_arcmin=fov,
        grid_size=grid_size,
        zodi_subtract=zodi_subtract,
        motion_exaggeration=motion_exaggeration,
        include_planetx_injection=include_planetx,
        fetch_live_cutout=fetch_live_cutout,
    )
    return JSONResponse(content=cube)


@app.get("/api/spectrum")
def get_point_spectrum_and_lightcurve(
    ra: float = Query(..., description="Field center Right Ascension (deg)"),
    dec: float = Query(..., description="Field center Declination (deg)"),
    fov: float = Query(4.5, description="Field of View in arcminutes"),
    x: Optional[float] = Query(None, description="Clicked pixel X coordinate"),
    y: Optional[float] = Query(None, description="Clicked pixel Y coordinate"),
    source_id: Optional[str] = Query(None, description="Optional explicit catalog source ID"),
    epoch: int = Query(0, description="Selected epoch index"),
    grid_size: int = Query(96, description="Spatial grid resolution in pixels"),
    motion_exaggeration: float = Query(1.0, description="Visual proper motion baseline multiplier"),
):
    """
    On-the-fly 102-channel SPHEREx Near-IR Spectrophotometer (0.75-5.0 µm) and Multi-Epoch Light Curve
    at any clicked pixel (x, y) or catalog source.
    """
    result = extract_spectrum_and_lightcurve(
        ra=ra,
        dec=dec,
        fov_arcmin=fov,
        click_x=x,
        click_y=y,
        source_id=source_id,
        epoch_index=epoch,
        grid_size=grid_size,
        motion_exaggeration=motion_exaggeration,
    )
    return JSONResponse(content=result)


@app.get("/api/raw-cutout")
def get_raw_irsa_cutout(
    ra: float = Query(269.452075, description="Right Ascension (deg)"),
    dec: float = Query(4.693389, description="Declination (deg)"),
    fov: float = Query(4.5, description="Field of View (arcmin)"),
    band: str = Query("SPHEREx-D2", description="Detector band (SPHEREx-D1..D6)"),
):
    """
    Fetch or reconstruct the native 6.15 arcsec/px SPHEREx Level-2 FITS cutout
    from NASA/IPAC IRSA IBE for direct raw-pixel inspection.
    """
    data = get_real_or_cached_irsa_cutout(ra=ra, dec=dec, fov_arcmin=fov, band=band)
    return JSONResponse(content=data)


@app.get("/api/export/fits")
def export_multi_epoch_fits(
    ra: float = Query(269.452075, description="Right Ascension (deg)"),
    dec: float = Query(4.693389, description="Declination (deg)"),
    fov: float = Query(4.5, description="Field of View (arcmin)"),
    band: str = Query("ALL", description="Detector band (SPHEREx-D1..D6 or ALL)"),
):
    """
    Export the current SPHEREx multi-epoch data cube, difference image, FLAGS map, and 102-channel
    spectrum as a standard NASA/OGIP Multi-Extension FITS (.fits) file for scientific analysis in DS9/Astropy.
    """
    cube = build_multi_epoch_cube(ra=ra, dec=dec, fov_arcmin=fov, grid_size=96)
    spec = extract_spectrum_and_lightcurve(ra=ra, dec=dec, fov_arcmin=fov, grid_size=96)

    obs_ra = cube["query"]["ra"]
    obs_dec = cube["query"]["dec"]
    grid_n = cube["query"]["grid_size"]
    pix_scale_deg = (fov / 60.0) / grid_n

    primary_hdu = fits.PrimaryHDU()
    hdr = primary_hdu.header
    hdr["TELESCOP"] = ("NASA SPHEREx", "Spectro-Photometer for History of Universe")
    hdr["INSTRUME"] = ("SPHEREx LVF D1-D6", "0.75-5.00 um Linear Variable Filters")
    hdr["RA_CEN"] = (obs_ra, "2025-2026 Epoch Field Center RA (deg)")
    hdr["DEC_CEN"] = (obs_dec, "2025-2026 Epoch Field Center Dec (deg)")
    hdr["RA_INP"] = (ra, "Input Reference RA (deg)")
    hdr["DEC_INP"] = (dec, "Input Reference Dec (deg)")
    hdr["FOV_ARCM"] = (fov, "Field of View (arcminutes)")
    hdr["N_EPOCHS"] = (len(cube["epochs"]), "Number of multi-epoch survey passes")
    hdr["CAT_SRC"] = (cube["query"].get("catalog_source", "IRSA"), "Catalog Cross-Match Source")
    hdr["BUNIT"] = ("MJy/sr", "Calibrated surface brightness units")

    # Stack epoch frames into a 3D array (N_epochs, H, W)
    if band in [b["id"] for b in SPHEREX_BANDS]:
        stack_3d = np.array([ep[band] for ep in cube["frames_by_epoch"]], dtype=np.float32)
        diff_2d = np.array(cube["difference_image"]["by_band"][band], dtype=np.float32)
    else:
        stack_3d = np.array(cube["coadd_by_epoch"], dtype=np.float32)
        diff_2d = np.array(cube["difference_image"]["coadd"], dtype=np.float32)

    def _apply_tan_wcs(hdu_hdr):
        hdu_hdr["CTYPE1"] = "RA---TAN"
        hdu_hdr["CTYPE2"] = "DEC--TAN"
        hdu_hdr["CRPIX1"] = (grid_n + 1) / 2.0
        hdu_hdr["CRPIX2"] = (grid_n + 1) / 2.0
        hdu_hdr["CRVAL1"] = obs_ra
        hdu_hdr["CRVAL2"] = obs_dec
        hdu_hdr["CDELT1"] = -pix_scale_deg
        hdu_hdr["CDELT2"] = -pix_scale_deg
        hdu_hdr["CUNIT1"] = "deg"
        hdu_hdr["CUNIT2"] = "deg"
        hdu_hdr["RADESYS"] = "ICRS"

    cube_hdu = fits.ImageHDU(data=stack_3d, name="MULTI_EPOCH_CUBE")
    cube_hdu.header["BAND"] = (band, "Selected SPHEREx detector band")
    _apply_tan_wcs(cube_hdu.header)

    diff_hdu = fits.ImageHDU(data=diff_2d, name="DIFF_IMAGE")
    diff_hdu.header["COMMENT"] = "Epoch Last minus Epoch First Difference Image"
    _apply_tan_wcs(diff_hdu.header)

    flags_hdu = fits.ImageHDU(
        data=np.array(cube["flags_by_epoch"][0], dtype=np.int32),
        name="FLAGS",
    )
    _apply_tan_wcs(flags_hdu.header)

    # 102-channel Spectrum Table HDU
    s102 = spec["spectrum_102"]
    col_wl = fits.Column(name="WAVELENGTH_UM", format="E", unit="um", array=np.array(s102["wavelengths_um"]))
    col_fl = fits.Column(name="FLUX_MJY", format="E", unit="mJy", array=np.array(s102["flux_mjy"]))
    col_er = fits.Column(name="FLUX_ERR_MJY", format="E", unit="mJy", array=np.array(s102["flux_err_mjy"]))
    col_bd = fits.Column(name="DETECTOR_BAND", format="12A", array=np.array(s102["band_ids"]))
    spec_hdu = fits.BinTableHDU.from_columns([col_wl, col_fl, col_er, col_bd], name="SPECTRUM_102")

    hdul = fits.HDUList([primary_hdu, cube_hdu, diff_hdu, flags_hdu, spec_hdu])
    buf = io.BytesIO()
    hdul.writeto(buf)
    buf.seek(0)

    filename = f"spherex_multiepoch_ra{ra:.3f}_dec{dec:+.3f}_{band}.fits"
    return StreamingResponse(
        buf,
        media_type="application/fits",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


# Mount static frontend files
if FRONTEND_DIR.exists():
    app.mount("/static", StaticFiles(directory=str(FRONTEND_DIR)), name="static")


@app.get("/")
def serve_index():
    index_path = FRONTEND_DIR / "index.html"
    if index_path.exists():
        return FileResponse(str(index_path))
    return JSONResponse({"message": "SPHEREx Backend Running. Frontend not yet built."})
