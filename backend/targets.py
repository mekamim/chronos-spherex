"""
Curated & Challenge-Required Astronomical Targets and Coordinate Resolver
for the NASA Space Apps 2026 'Planet X and SPHEREx' Challenge.

Includes:
- Exact Reference Epoch (J2000.0 / 2010.5) AND Exact 2025.6 SPHEREx Survey Epoch
  astrometric coordinates propagated via Astropy SkyCoord.apply_space_motion().
- Built-in astronomical reference resolver for famous deep-sky & nearby targets
  plus Astropy Sesame fallback.
"""

from typing import Dict, Any, Tuple, Optional, List
import math
import re
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
import astropy.units as u


# The 6 SPHEREx Linear Variable Filter (LVF) Detector Bands (102 spectral channels total: 17 per band)
SPHEREX_BANDS: List[Dict[str, Any]] = [
    {
        "id": "SPHEREx-D1",
        "short": "D1",
        "name": "Detector 1 (Band 1)",
        "short_label": "D1 (0.75–1.11 μm)",
        "wavelength_min_um": 0.75,
        "wavelength_max_um": 1.11,
        "resolving_power": 41,
        "channels": 17,
        "color": "#60a5fa",
        "molecules": ["Pa-δ", "He I (1.08 μm)", "FeH", "TiO/VO"],
        "description": "Short-wavelength near-IR (Y/z-band); sensitive to stellar photospheres and T/Y dwarf Y-band peaks.",
    },
    {
        "id": "SPHEREx-D2",
        "short": "D2",
        "name": "Detector 2 (Band 2)",
        "short_label": "D2 (1.10–1.64 μm)",
        "wavelength_min_um": 1.10,
        "wavelength_max_um": 1.64,
        "resolving_power": 41,
        "channels": 17,
        "color": "#34d399",
        "molecules": ["Pa-β (1.28 μm)", "H₂O (1.15, 1.4 μm)", "CH₄ (1.15, 1.3 μm)"],
        "description": "J/H-band near-IR; captures deep water and methane absorption bands in brown dwarfs and planetary atmospheres.",
    },
    {
        "id": "SPHEREx-D3",
        "short": "D3",
        "name": "Detector 3 (Band 3)",
        "short_label": "D3 (1.63–2.42 μm)",
        "wavelength_min_um": 1.63,
        "wavelength_max_um": 2.42,
        "resolving_power": 41,
        "channels": 17,
        "color": "#a3e635",
        "molecules": ["Pa-α (1.875 μm)", "Br-γ (2.166 μm)", "H₂O (1.9 μm)", "CO (2.3 μm)"],
        "description": "H/K-band near-IR; sensitive to H⁻ opacity minimum at 1.6 μm, hydrogen recombination lines, and CO overtone bands.",
    },
    {
        "id": "SPHEREx-D4",
        "short": "D4",
        "name": "Detector 4 (Band 4)",
        "short_label": "D4 (2.41–3.82 μm)",
        "wavelength_min_um": 2.41,
        "wavelength_max_um": 3.82,
        "resolving_power": 35,
        "channels": 17,
        "color": "#facc15",
        "molecules": ["H₂O Ice (3.05 μm)", "PAH (3.3 μm)", "CH₄ (3.32 μm)", "Br-β (2.63 μm)"],
        "description": "Mid-near-IR (WISE W1 equivalent); probes water ice at 3.0 μm, PAH dust emission at 3.3 μm, and methane fundamentals.",
    },
    {
        "id": "SPHEREx-D5",
        "short": "D5",
        "name": "Detector 5 (Band 5)",
        "short_label": "D5 (3.81–4.42 μm)",
        "wavelength_min_um": 3.81,
        "wavelength_max_um": 4.42,
        "resolving_power": 110,
        "channels": 17,
        "color": "#fb923c",
        "molecules": ["Br-α (4.05 μm)", "CO₂ Ice/Gas (4.27 μm)", "SiO (4.05 μm)"],
        "description": "High spectral resolution (R=110) band targeting carbon dioxide (CO₂) at 4.27 μm and Brackett-α hydrogen emission.",
    },
    {
        "id": "SPHEREx-D6",
        "short": "D6",
        "name": "Detector 6 (Band 6)",
        "short_label": "D6 (4.42–5.00 μm)",
        "wavelength_min_um": 4.42,
        "wavelength_max_um": 5.00,
        "resolving_power": 130,
        "channels": 17,
        "color": "#f87171",
        "molecules": ["CO Gas/Ice (4.67 μm)", "OCN⁻ (4.62 μm)", "Y-Dwarf / Planet X 4.7 μm Peak"],
        "description": "Highest spectral resolution (R=130) thermal near-IR (WISE W2 equivalent); crucial for detecting cold brown dwarfs (like WISE 0855) and Planet X candidates.",
    },
]


# Prominent spectral lines / molecular absorption bands in SPHEREx range (0.75 - 5.00 um)
SPECTRAL_FEATURES: List[Dict[str, Any]] = [
    {"wavelength_um": 0.98, "label": "FeH / Y-band", "type": "molecular", "molecule": "FeH"},
    {"wavelength_um": 1.083, "label": "He I 1.08μm", "type": "emission", "molecule": "He I"},
    {"wavelength_um": 1.15, "label": "H₂O / CH₄", "type": "absorption", "molecule": "H₂O+CH₄"},
    {"wavelength_um": 1.282, "label": "Pa-β", "type": "emission", "molecule": "H I"},
    {"wavelength_um": 1.40, "label": "H₂O 1.4μm", "type": "absorption", "molecule": "H₂O"},
    {"wavelength_um": 1.65, "label": "CH₄ 1.65μm", "type": "absorption", "molecule": "CH₄"},
    {"wavelength_um": 1.90, "label": "H₂O 1.9μm", "type": "absorption", "molecule": "H₂O"},
    {"wavelength_um": 2.166, "label": "Br-γ", "type": "emission", "molecule": "H I"},
    {"wavelength_um": 2.32, "label": "CO / CH₄ 2.3μm", "type": "absorption", "molecule": "CO / CH₄"},
    {"wavelength_um": 2.72, "label": "H₂O 2.7μm", "type": "absorption", "molecule": "H₂O"},
    {"wavelength_um": 3.05, "label": "H₂O Ice 3.05μm", "type": "absorption", "molecule": "H₂O Ice"},
    {"wavelength_um": 3.32, "label": "CH₄ / PAH 3.3μm", "type": "absorption", "molecule": "CH₄ / PAH"},
    {"wavelength_um": 4.052, "label": "Br-α / SiO", "type": "emission", "molecule": "H I / SiO"},
    {"wavelength_um": 4.27, "label": "CO₂ 4.27μm", "type": "absorption", "molecule": "CO₂"},
    {"wavelength_um": 4.67, "label": "CO / Y-Dwarf Window", "type": "feature", "molecule": "CO / Thermal Peak"},
]


def propagate_astrometry_to_epoch(
    ra_deg: float,
    dec_deg: float,
    pm_ra_masyr: float,
    pm_dec_masyr: float,
    parallax_mas: float,
    ref_epoch_yr: float = 2000.0,
    target_epoch_yr: float = 2025.64,
) -> Tuple[float, float, float]:
    """
    Propagate a source's sky position from ref_epoch_yr (e.g., J2000.0 or WISE 2010.5)
    to target_epoch_yr (e.g., SPHEREx 2025.64) using rigorous 3D Astropy space motion.
    Returns (ra_target_deg, dec_target_deg, total_shift_arcsec).
    """
    try:
        c0 = SkyCoord(
            ra=ra_deg * u.deg,
            dec=dec_deg * u.deg,
            pm_ra_cosdec=pm_ra_masyr * u.mas / u.yr,
            pm_dec=pm_dec_masyr * u.mas / u.yr,
            distance=Distance(parallax=max(parallax_mas, 0.1) * u.mas),
            obstime=Time(ref_epoch_yr, format="decimalyear"),
        )
        c1 = c0.apply_space_motion(new_obstime=Time(target_epoch_yr, format="decimalyear"))
        sep_arcsec = float(c0.separation(c1).arcsec)
        return round(float(c1.ra.deg), 6), round(float(c1.dec.deg), 6), round(sep_arcsec, 2)
    except Exception:
        dt = target_epoch_yr - ref_epoch_yr
        cos_d = max(math.cos(math.radians(dec_deg)), 0.05)
        dra = (pm_ra_masyr / 3600000.0) * dt / cos_d
        ddec = (pm_dec_masyr / 3600000.0) * dt
        sep = math.hypot(pm_ra_masyr, pm_dec_masyr) * abs(dt) / 1000.0
        return round(ra_deg + dra, 6), round(dec_deg + ddec, 6), round(sep, 2)


# The 5 Mandatory NASA Space Apps Challenge Test Targets
CHALLENGE_TARGETS: List[Dict[str, Any]] = [
    {
        "id": "barnards_star",
        "name": "Barnard's Star",
        "subtitle": "Highest Proper-Motion Star in the Sky (10.39\"/yr)",
        "category": "High Proper Motion Star",
        "badge": "Challenge Target #1",
        "ra": 269.452075,
        "dec": 4.693389,
        "ra_deg": 269.452075,
        "dec_deg": 4.693389,
        "ref_epoch": 2000.0,
        "ra_2025": 269.446347,
        "dec_2025": 4.767188,
        "epoch_shift_arcsec": 266.5,
        "ra_hms": "17h57m48.498s",
        "dec_dms": "+04d41m36.20s",
        "default_fov_arcmin": 4.5,
        "recommended_band": "SPHEREx-D2",
        "pm_ra_masyr": -801.55,
        "pm_dec_masyr": 10362.39,
        "total_pm_arcsecyr": 10.39,
        "distance_pc": 1.828,
        "parallax_mas": 546.98,
        "spectral_type": "M4.0V Red Dwarf",
        "temperature_k": 3134,
        "scientific_notes": (
            "Barnard's Star (GJ 699, Gaia DR3 4472832130942575872) is a nearby M4.0V red dwarf at 1.83 pc (5.96 ly) with the highest "
            "proper motion of any known star (10.39\"/yr northward). Between J2000.0 (+04°41'36\") and SPHEREx's 2025.64 observation epoch "
            "(+04°46'02\"), it has marched 266.5\" (4.44 arcmin) northward! Across SPHEREx's 2025–2026 passes, it produces a textbook "
            "10.4\"/yr positive-negative dipole in Difference Imaging."
        ),
        "discovery_tip": "Inspect both its 2025–2026 dipole shift (10.39\"/yr) and its 266.5\" cumulative trajectory from its J2000 reference coordinates!",
    },
    {
        "id": "wise_0855",
        "name": "WISE 0855",
        "subtitle": "WISE J085510.83−071442.5 • Coldest Known Sub-Brown Dwarf (~250 K)",
        "category": "Planet X / Y-Dwarf Analog",
        "badge": "Challenge Target #2",
        "ra": 133.795129,
        "dec": -7.245147,
        "ra_deg": 133.795129,
        "dec_deg": -7.245147,
        "ref_epoch": 2010.5,
        "ra_2025": 133.761215,
        "dec_2025": -7.242342,
        "epoch_shift_arcsec": 121.5,
        "ra_hms": "08h55m10.831s",
        "dec_dms": "-07d14m42.53s",
        "default_fov_arcmin": 4.5,
        "recommended_band": "SPHEREx-D6",
        "pm_ra_masyr": -8118.0,
        "pm_dec_masyr": 680.0,
        "total_pm_arcsecyr": 8.15,
        "distance_pc": 2.278,
        "parallax_mas": 439.0,
        "spectral_type": "Y4 Ultra-Cool Sub-Brown Dwarf (~3–10 M_Jup)",
        "temperature_k": 250,
        "scientific_notes": (
            "WISE J085510.83−071442.5 is the coldest known free-floating planetary-mass object near the Sun (~250 K) at 2.28 pc, "
            "moving rapidly westward at 8.15\"/yr (121.5\" west of its 2010.5 WISE discovery position by 2025.4!). Completely invisible "
            "in Gaia optical and SPHEREx D1–D3 (<0.001 mJy), it emerges in SPHEREx Detector D6 (4.42–5.00 μm, W2=13.89 mag ≈ 0.48 mJy) "
            "at SNR ≈ 2.8 in single frames and SNR > 11 when co-added!"
        ),
        "discovery_tip": "Select Detector D6 (4.42–5.00 μm) and compare Single-Frame vs. Co-Add Σ to see how multi-exposure stacking boosts faint Y-dwarf SNR!",
    },
    {
        "id": "epsilon_indi_b",
        "name": "Epsilon Indi B",
        "subtitle": "Nearest Binary Brown Dwarf System (ε Indi Ba + Bb, T1 + T6)",
        "category": "High-PM Brown Dwarf Binary",
        "badge": "Challenge Target #3",
        "ra": 331.043833,
        "dec": -56.782694,
        "ra_deg": 331.043833,
        "dec_deg": -56.782694,
        "ref_epoch": 2000.0,
        "ra_2025": 331.094712,
        "dec_2025": -56.800341,
        "epoch_shift_arcsec": 119.5,
        "ra_hms": "22h04m10.52s",
        "dec_dms": "-56d46m57.7s",
        "default_fov_arcmin": 4.5,
        "recommended_band": "SPHEREx-D2",
        "pm_ra_masyr": 3976.0,
        "pm_dec_masyr": -2486.0,
        "total_pm_arcsecyr": 4.69,
        "distance_pc": 3.639,
        "parallax_mas": 274.8,
        "spectral_type": "T1 + T6 Brown Dwarf Binary",
        "temperature_k": 1280,
        "scientific_notes": (
            "Epsilon Indi Ba+Bb (Gaia DR3 6412596012146801152) is the closest binary brown dwarf system to Earth (3.64 pc), consisting of a "
            "T1 dwarf (~1280 K) and a T6 dwarf (~850 K). Moving at 4.69\"/yr southeast (119.5\" from J2000 to 2025.5), its 102-channel "
            "SPHEREx spectrum anchored to 2MASS J/H/Ks and WISE W1/W2 photometry exhibits deep H₂O and CH₄ absorption valleys."
        ),
        "discovery_tip": "Click on Epsilon Indi B to inspect its 102-band SPHEREx spectrum anchored to real 2MASS + WISE flux measurements!",
    },
    {
        "id": "rw_aurigae",
        "name": "RW Aurigae",
        "subtitle": "Eruptive Classical T Tauri Variable with Dusty Circumstellar Disk",
        "category": "Young Eruptive Variable",
        "badge": "Challenge Target #4",
        "ra": 76.956538,
        "dec": 30.401408,
        "ra_deg": 76.956538,
        "dec_deg": 30.401408,
        "ref_epoch": 2000.0,
        "ra_2025": 76.956572,
        "dec_2025": 30.401233,
        "epoch_shift_arcsec": 0.64,
        "ra_hms": "05h07m49.569s",
        "dec_dms": "+30d24m05.07s",
        "default_fov_arcmin": 4.5,
        "recommended_band": "SPHEREx-D3",
        "pm_ra_masyr": 4.12,
        "pm_dec_masyr": -24.56,
        "total_pm_arcsecyr": 0.025,
        "distance_pc": 163.4,
        "parallax_mas": 6.12,
        "spectral_type": "K1–K4Ve Classical T Tauri Binary",
        "temperature_k": 4600,
        "scientific_notes": (
            "RW Aurigae (Gaia DR3 156430822114424576) is a young pre-main-sequence T Tauri star system in Taurus-Auriga. "
            "Rather than moving across the sky, it exhibits dramatic brightness changes across SPHEREx's 2025–2026 passes caused by "
            "warped tidal dust streams occulting the star and episodic accretion bursts, with strong 3–5 μm warm dust excess and H I emission."
        ),
        "discovery_tip": "Switch to 'Dipole Δ' to see its fixed-position monopole flux residual, and inspect the Time Curve across real 2025–2026 IRSA passes!",
    },
    {
        "id": "r_draconis",
        "name": "R Draconis",
        "subtitle": "Pulsating Mira AGB Red Giant Variable (Period = 245.6 days)",
        "category": "Pulsating Mira Variable",
        "badge": "Challenge Target #5",
        "ra": 248.167612,
        "dec": 66.754964,
        "ra_deg": 248.167612,
        "dec_deg": 66.754964,
        "ref_epoch": 2000.0,
        "ra_2025": 248.167434,
        "dec_2025": 66.754840,
        "epoch_shift_arcsec": 0.51,
        "ra_hms": "16h32m40.227s",
        "dec_dms": "+66d45m17.87s",
        "default_fov_arcmin": 4.5,
        "recommended_band": "SPHEREx-D4",
        "pm_ra_masyr": -9.85,
        "pm_dec_masyr": -17.42,
        "total_pm_arcsecyr": 0.020,
        "distance_pc": 431.0,
        "parallax_mas": 2.32,
        "spectral_type": "M6e–M9e III Mira AGB Variable",
        "temperature_k": 2650,
        "scientific_notes": (
            "Located near the North Ecliptic Pole (where SPHEREx's polar orbit overlaps on almost every revolution, yielding 551 IRSA exposures!), "
            "R Draconis (Gaia DR3 1635876439281729408) is a classic Mira-type Asymptotic Giant Branch (AGB) pulsating variable star with a "
            "245.6-day period. Its near-IR spectrum exhibits deep H₂O steam and CO bands that vary across the 245.6-day pulsation cycle."
        ),
        "discovery_tip": "Notice the 551 real IRSA exposures near the North Ecliptic Pole and track its 245.6-day Mira pulsation in the Time Curve panel.",
    },
]


# Built-in Astronomical Reference Dictionary for instant, offline-safe object name resolution
WELL_KNOWN_OBJECTS: Dict[str, Tuple[float, float, str]] = {
    "m31": (10.684708, 41.268750, "M31 (Andromeda Galaxy)"),
    "andromeda": (10.684708, 41.268750, "M31 (Andromeda Galaxy)"),
    "m42": (83.822080, -5.391110, "M42 (Orion Nebula)"),
    "orion": (83.822080, -5.391110, "M42 (Orion Nebula)"),
    "orionnebula": (83.822080, -5.391110, "M42 (Orion Nebula)"),
    "m45": (56.750000, 24.116700, "M45 (Pleiades Cluster)"),
    "pleiades": (56.750000, 24.116700, "M45 (Pleiades Cluster)"),
    "m1": (83.633080, 22.014500, "M1 (Crab Nebula)"),
    "crabnebula": (83.633080, 22.014500, "M1 (Crab Nebula)"),
    "wolf359": (164.120270, 7.014660, "Wolf 359 (CN Leonis, High-PM M6.5V Dwarf)"),
    "cnleonis": (164.120270, 7.014660, "Wolf 359 (CN Leonis)"),
    "proximacentauri": (217.428940, -62.679480, "Proxima Centauri (GJ 551)"),
    "proxima": (217.428940, -62.679480, "Proxima Centauri (GJ 551)"),
    "alphacentauri": (219.902060, -60.833990, "Alpha Centauri AB"),
    "sirius": (101.287155, -16.716116, "Sirius (Alpha Canis Majoris)"),
    "vega": (279.234735, 38.783689, "Vega (Alpha Lyrae, Calibrator)"),
    "betelgeuse": (88.792939, 7.407064, "Betelgeuse (Alpha Orionis)"),
    "trappist1": (346.622370, -5.041400, "TRAPPIST-1 (Ultra-Cool M8V Planetary System)"),
    "luhman16": (162.328810, -53.319470, "Luhman 16 (WISE J1049−5319 L/T Binary)"),
    "ross128": (176.934990, 0.804560, "Ross 128 (GJ 447)"),
    "teegardensstar": (43.253510, 16.881290, "Teegarden's Star (High-PM M7V)"),
    "teegarden": (43.253510, 16.881290, "Teegarden's Star (High-PM M7V)"),
    "61cygni": (316.724750, 38.749420, "61 Cygni AB (High-PM K5V+K7V Binary)"),
    "tauceti": (26.017010, -15.937480, "Tau Ceti (HD 10700)"),
    "epsiloneridani": (53.232685, -9.458262, "Epsilon Eridani (HD 22049)"),
    "lalande21185": (165.834140, 35.969880, "Lalande 21185 (HD 95735)"),
    "kapteynsstar": (77.919090, -45.018410, "Kapteyn's Star (HD 33793)"),
    "sgra": (266.416837, -29.007811, "Sagittarius A* (Galactic Center)"),
    "nep": (270.000000, 66.560708, "North Ecliptic Pole (SPHEREx Deep Field)"),
    "sep": (90.000000, -66.560708, "South Ecliptic Pole (SPHEREx Deep Field)"),
}


def find_matching_challenge_target(
    ra: float, dec: float, tol_deg: float = 0.12
) -> Optional[Dict[str, Any]]:
    """
    Return the matching challenge target dict if (ra, dec) is within tol_deg of either
    its reference epoch position (J2000/2010.5) OR its 2025.6 propagated SPHEREx position.
    """
    for t in CHALLENGE_TARGETS:
        cos_d = max(math.cos(math.radians(dec)), 0.05)
        dist_ref = math.hypot((t["ra"] - ra) * cos_d, t["dec"] - dec)
        dist_2025 = math.hypot((t["ra_2025"] - ra) * cos_d, t["dec_2025"] - dec)
        if min(dist_ref, dist_2025) <= tol_deg:
            return t
    return None


def parse_sky_coordinates(
    ra_input: str, dec_input: Optional[str] = None
) -> Tuple[float, float, Optional[Dict[str, Any]]]:
    """
    Parse user-provided coordinates or object name into (ra_deg, dec_deg, matched_target).
    Supports:
    - Preset target IDs or names ('WISE 0855', 'Barnard's Star', 'Epsilon Indi B', 'RW Aurigae', 'R Draconis')
    - Famous astronomical targets ('M31', 'Wolf 359', 'Proxima Centauri', 'TRAPPIST-1', 'Luhman 16', 'NEP', etc.)
    - Decimal degrees: '133.795129, -7.245147'
    - Sexagesimal HMS/DMS: '08h55m10.831s -07d14m42.53s' or '17:57:48.498 +04:41:36.20'
    - Astropy Sesame name resolution fallback
    """
    raw_ra = (ra_input or "").strip()
    raw_dec = (dec_input or "").strip()

    norm_query = re.sub(r"[^a-z0-9]", "", raw_ra.lower())
    if norm_query and not raw_dec:
        for tinfo in CHALLENGE_TARGETS:
            tid_norm = re.sub(r"[^a-z0-9]", "", tinfo["id"].lower())
            tname_norm = re.sub(r"[^a-z0-9]", "", tinfo["name"].lower())
            if norm_query in (tid_norm, tname_norm) or (
                len(norm_query) >= 4 and (norm_query in tname_norm or norm_query in tid_norm)
            ):
                return tinfo["ra"], tinfo["dec"], tinfo

        if norm_query in WELL_KNOWN_OBJECTS:
            w_ra, w_dec, _ = WELL_KNOWN_OBJECTS[norm_query]
            return w_ra, w_dec, find_matching_challenge_target(w_ra, w_dec)
        for k, (w_ra, w_dec, _) in WELL_KNOWN_OBJECTS.items():
            if len(norm_query) >= 4 and (norm_query in k or k in norm_query):
                return w_ra, w_dec, find_matching_challenge_target(w_ra, w_dec)

    if raw_ra and raw_dec:
        try:
            ra_f = float(raw_ra.replace("deg", "").replace("d", "").strip())
            dec_f = float(raw_dec.replace("deg", "").replace("d", "").strip())
            if 0.0 <= ra_f <= 360.0 and -90.0 <= dec_f <= 90.0:
                return ra_f, dec_f, find_matching_challenge_target(ra_f, dec_f)
        except ValueError:
            pass

        combined = f"{raw_ra} {raw_dec}"
        for units_try in [(u.hourangle, u.deg), (u.deg, u.deg)]:
            try:
                c = SkyCoord(combined, unit=units_try, frame="icrs")
                ra_f, dec_f = float(c.ra.deg), float(c.dec.deg)
                return ra_f, dec_f, find_matching_challenge_target(ra_f, dec_f)
            except Exception:
                pass

    if raw_ra and not raw_dec:
        parts = [p.strip() for p in re.split(r"[,;\s]+", raw_ra) if p.strip()]
        if len(parts) == 2:
            try:
                ra_f = float(parts[0].replace("deg", "").replace("d", ""))
                dec_f = float(parts[1].replace("deg", "").replace("d", ""))
                if 0.0 <= ra_f <= 360.0 and -90.0 <= dec_f <= 90.0:
                    return ra_f, dec_f, find_matching_challenge_target(ra_f, dec_f)
            except ValueError:
                pass

        for units_try in [(u.hourangle, u.deg), (u.deg, u.deg)]:
            try:
                c = SkyCoord(raw_ra, unit=units_try, frame="icrs")
                ra_f, dec_f = float(c.ra.deg), float(c.dec.deg)
                return ra_f, dec_f, find_matching_challenge_target(ra_f, dec_f)
            except Exception:
                pass

        try:
            c = SkyCoord.from_name(raw_ra)
            ra_f, dec_f = float(c.ra.deg), float(c.dec.deg)
            return ra_f, dec_f, find_matching_challenge_target(ra_f, dec_f)
        except Exception as exc:
            raise ValueError(
                f"Could not parse coordinates or resolve object name '{raw_ra}'. "
                f"Try decimal degrees (e.g. '133.7951, -7.2451') or HMS/DMS ('08h55m10.831s -07d14m42.53s')."
            ) from exc

    raise ValueError("Please provide valid Right Ascension (RA) and Declination (Dec) coordinates.")
