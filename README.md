# CHRONOS // SPHEREx

**Live Observatory:** [https://chronosspherex.earth](https://chronosspherex.earth)  
**NASA Space Apps Challenge 2026 — Challenge: "Planet X and SPHEREx"**

An interactive, public-facing astronomical discovery and time-domain spectrophotometry web application for NASA's **SPHEREx** (*Spectro-Photometer for the History of the Universe, Epoch of Reionization and Ices Explorer*) mission.

Inspired by **WiseView** and engineered specifically for SPHEREx's **6 Linear Variable Filter (LVF) detector arrays** covering **102 near-infrared spectral channels ($0.75\text{–}5.00\,\mu\text{m}$)**, this platform enables scientists, students, and citizen astronomers to hunt for **Planet X**, fast-moving ultra-cool brown dwarfs (Y/T/L dwarfs), Solar System objects, and eruptive/pulsating variable stars across multi-epoch SPHEREx sky surveys.

---

## Key Features & Challenge Alignment

### 1. Multi-Epoch Sky Blinker & Time-Lapse Animator
- **Coordinate & Object Search**: Supports decimal degrees (`269.452075, 4.693389`), sexagesimal (`17h57m48.5s +04d41m36s`), or astronomical object names.
- **Interactive Blink Modes**:
  - **Ping-Pong ($1 \leftrightarrow N$)**: Smooth back-and-forth oscillation across SPHEREx survey passes.
  - **WiseView Classic A/B Blink ($1 \leftrightarrow N$)**: Rapid two-epoch toggling to spot subtle astrometric shifts.
  - **Forward Time-Lapse Loop ($1 \rightarrow N$)** with adjustable frame rate (`1.0` to `10.0 FPS`) and manual epoch scrubbing.

### 2. Verified Challenge & Benchmark Targets
1. **Barnard's Star** (`RA 269.452075°, Dec +4.693389°`): Fastest proper-motion star in the sky ($10.39''/\text{yr}$ northward).
2. **WISE 0855−0714** (`RA 133.795129°, Dec −7.245147°`): Coldest known sub-brown dwarf ($\sim 250\text{ K}$, $8.15''/\text{yr}$). Undetected in Gaia optical and SPHEREx Bands `D1–D3`, emerging in **SPHEREx Detector D6 ($4.42\text{–}5.00\,\mu\text{m}$)**.
3. **Epsilon Indi Ba+Bb** (`RA 331.043833°, Dec −56.782694°`): Binary `T1 + T6` brown dwarf system ($4.69''/\text{yr}$) with deep $\text{H}_2\text{O}$ and $\text{CH}_4$ absorption valleys across SPHEREx's 102 channels.
4. **RW Aurigae** (`RA 76.956538°, Dec +30.401408°`): Young Classical T Tauri variable binary with circumstellar dust disk occultations, $\text{H I}$ emission lines ($\text{Pa}\beta, \text{Br}\gamma$), and warm dust excess.
5. **R Draconis** (`RA 248.167612°, Dec +66.754964°`): Mira-type AGB pulsating variable star ($P = 245.6\text{ d}$) near the North Ecliptic Pole ($>350$ SPHEREx exposures).

### 3. Five Specialized SPHEREx Imaging & Stacking Modes
- **Single Detector (`D1`–`D6`)**: Inspect individual SPHEREx LVF arrays from $0.75\,\mu\text{m}$ to $5.00\,\mu\text{m}$.
- **All-Band Co-Add (`D1–D6 Stack`)**: Combines all six detector arrays to maximize Signal-to-Noise Ratio (SNR).
- **False-Color Infrared RGB Composite**: Maps `Red = D5+D6 (3.8–5.0 μm)`, `Green = D3+D4 (1.6–3.8 μm)`, `Blue = D1+D2 (0.75–1.6 μm)` so ultra-cool brown dwarfs and Planet X candidates stand out in vivid crimson/magenta.
- **Difference Imaging (`Epoch N − Epoch 1`)**: Subtracts the earliest survey pass from the latest pass. Moving objects produce a sharp **Blue-to-Red Dipole**, while variable stars produce a **Monopole Flux Residual**.
- **Temporal Co-Add (`Σ Epochs`)**: Stacks all epochs together to reveal faint static background structures.

### 4. On-the-Fly 102-Band Spectrophotometer & Real Photometric Anchors (`mJy`)
- Click **any source or pixel** on the sky canvas to immediately extract its **102-channel SPHEREx Near-Infrared Spectrum ($0.75\text{–}5.00\,\mu\text{m}$)** with molecular band markers ($\text{H}_2\text{O}, \text{CH}_4, \text{CO}, \text{CO}_2, \text{PAH}, \text{H I}$), anchored to real Vega-calibrated **2MASS ($J, H, K_s$)** and **WISE ($W_1, W_2$)** flux density measurements in physical `mJy` (`◆` gold diamond markers), alongside its **Multi-Epoch Light Curve** and **Astrometric Motion Track**.

### 5. Astrometric Epoch Propagation, WCS North-Up Reprojection & Real IRSA Catalogs
- **3D Space-Motion Epoch Propagation (`J2000 / 2010.5 → 2025.6 SPHEREx Epoch`)**: Propagates high-proper-motion targets via `astropy.coordinates.SkyCoord.apply_space_motion()` so 2025–2026 SPHEREx cutouts center on the star's true 2025–2026 sky coordinates (e.g., Barnard's Star $+266.5''$ N shift; WISE 0855 $-121.5''$ W shift; $\varepsilon$ Indi B $119.5''$ SE shift).
- **WCS North-Up (`N ↑ E ←`) FITS Reprojection**: Reprojects raw detector-frame SPHEREx Level-2 FITS cutouts (arbitrary spacecraft roll angle) onto a North-Up, East-Left TAN WCS grid using `astropy.wcs.WCS` and `scipy.ndimage.map_coordinates`.
- **Live NASA/IPAC IRSA Gator Catalog Cross-Match (`gaia_dr3_source` + `allwise_p3as_psd`)**: Queries real ESA Gaia DR3 optical sources (real `source_id`, `pmra`, `pmdec`, `parallax`, `phot_g_mean_mag`, `teff_gspphot`) and NASA AllWISE + 2MASS infrared sources (`w1mpro`, `w2mpro`, `j_m_2mass`, `h_m_2mass`, `k_m_2mass`).
- **Zodiacal Light (`ZODI`) Toggle & 32-Bit `FLAGS` Bitmask Inspector**: Compare raw interplanetary dust continuum against zodiacal-subtracted images and decode HDU 2 contamination flags (`TRANSIENT`, `HOT`, `GHOST`, `SUR_ERROR`, `SOURCE`, `HALO`).
- **Cross-Mission Links**: One-click deep links to **NASA IRSA SPHEREx Explorer**, **IRSA FinderChart (WISE/2MASS)**, **WiseView**, **CSA NEOSSat / CADC**, **SIMBAD**, and **ESA Gaia DR3**, plus **Multi-Extension FITS (`.fits`)** and **CSV** exports.

---

## Data Sources & Methodology Transparency

1. **Archival Data Provenance:** Sky coordinates, multi-epoch observation timestamps (`MJD`), spacecraft WCS footprints, and raw Level-2 spectral detector cutouts (`RAW_FITS`) are queried directly from the **NASA/IPAC Infrared Science Archive (IRSA) SPHEREx SIAv2 service** (`spherex_qr3` / `spherex_dr1`) and cross-matched with **ESA Gaia DR3 (`gaia_dr3_source`)** and **NASA AllWISE + 2MASS (`allwise_p3as_psd`)** catalogs via IRSA Gator.
2. **Forward Modeling & Photometric Calibration:** Because individual SPHEREx exposures cover a single Linear Variable Filter strip per pass, `MODEL_CUBE` reconstructs a continuous 6-band (`D1–D6`) + 102-channel (`0.75–5.00 μm`) spectral cube anchored to cataloged `2MASS (J, H, Ks)` and `WISE (W1, W2)` flux densities (`mJy`) and brown-dwarf/stellar spectral templates ($\text{H}_2\text{O}, \text{CH}_4, \text{CO}, \text{CO}_2, \text{PAH}, \text{H I}$). Switch to **`RAW_FITS`** in the dock at any time to inspect the un-modeled, WCS North-Up reprojected NASA Level-2 FITS detector frame.

---

## Quick Start (Local Development)

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn backend.main:app --host 0.0.0.0 --port 8000
```

Open `http://localhost:8000` in your browser.

---

## Cloud Deployment (`chronosspherex.earth` via Google Cloud Run)

```bash
# 1. Deploy directly from source to Google Cloud Run
gcloud run deploy chronos-spherex \
  --source . \
  --region us-central1 \
  --allow-unauthenticated \
  --memory 1Gi

# 2. Map custom domain (chronosspherex.earth)
gcloud beta run domain-mappings create \
  --service chronos-spherex \
  --domain chronosspherex.earth \
  --region us-central1
```

---

## License

Released under the [MIT License](LICENSE).
