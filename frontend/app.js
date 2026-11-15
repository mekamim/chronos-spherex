/**
 * SPHEREx // CHRONOS — Observatory Viewport & Spectral Analysis Client
 */

const state = {
  targets: [],
  bands: [],
  spectralFeatures: [],
  flagDefinitions: [],

  ra: 269.452075,
  dec: 4.693389,
  fov: 4.5,
  motionExaggeration: 1.0,
  zodiSubtract: true,

  // SINGLE_BAND | SPLIT_CURTAIN | DIFFERENCE | RGB_COMPOSITE
  viewMode: "SINGLE_BAND",
  activeBand: "SPHEREx-D1",
  activeEpoch: 0,

  splitRatio: 0.5,
  isDraggingSplit: false,

  isPlaying: true,
  blinkPattern: "PINGPONG",
  blinkDirection: 1,
  fps: 2.6,
  timerId: null,

  stretch: "asinh",
  colormap: "inferno",
  contrast: 1.45,
  bias: 0.04,

  showGaiaOverlay: true,
  showDetectionsOverlay: true,
  showFlagsOverlay: false,
  showWcsGrid: true,
  showMolecularLines: true,
  lightCurveMode: "FLUX",

  skyCube: null,
  spectrumData: null,
  rawCutout: null,
  selectedSource: null,

  apertureX: 47.5,
  apertureY: 47.5,
  hoverX: null,
  hoverY: null,
  hoveredSpecChannel: null,

  // Animation progress for the 102-band spectrum curve (0..1)
  specAnimProgress: 1.0,
  specAnimFrame: null,
};

const offscreenA = document.createElement("canvas");
const ctxA = offscreenA.getContext("2d");
const offscreenB = document.createElement("canvas");
const ctxB = offscreenB.getContext("2d");

// Ambient Wavelength Aura & Colormaps

const AURA_BY_BAND = {
  "SPHEREx-D1": "rgba(96, 165, 250, 0.20)",
  "SPHEREx-D2": "rgba(52, 211, 153, 0.19)",
  "SPHEREx-D3": "rgba(163, 230, 53, 0.17)",
  "SPHEREx-D4": "rgba(250, 204, 21, 0.18)",
  "SPHEREx-D5": "rgba(251, 146, 60, 0.20)",
  "SPHEREx-D6": "rgba(248, 113, 113, 0.24)",
};

function updateAmbientAura() {
  const aura = document.getElementById("ambientAura");
  if (!aura) return;
  let color = AURA_BY_BAND[state.activeBand] || "rgba(56, 189, 248, 0.18)";
  if (state.viewMode === "DIFFERENCE") {
    color = "rgba(245, 176, 65, 0.22)";
  } else if (state.viewMode === "RGB_COMPOSITE") {
    color = "rgba(192, 132, 252, 0.22)";
  }
  document.documentElement.style.setProperty("--aura-color", color);
}

function clamp01(t) {
  return Math.max(0.0, Math.min(1.0, t));
}

function interpolateStops(stops, t) {
  const x = clamp01(t) * (stops.length - 1);
  const i = Math.floor(x);
  const f = x - i;
  if (i >= stops.length - 1) return stops[stops.length - 1];
  const c0 = stops[i];
  const c1 = stops[i + 1];
  return [
    Math.round(c0[0] + (c1[0] - c0[0]) * f),
    Math.round(c0[1] + (c1[1] - c0[1]) * f),
    Math.round(c0[2] + (c1[2] - c0[2]) * f),
  ];
}

const COLORMAP_STOPS = {
  inferno: [
    [2, 3, 10],
    [32, 11, 68],
    [102, 21, 104],
    [176, 50, 72],
    [234, 110, 30],
    [249, 192, 44],
    [252, 254, 170],
  ],
  viridis: [
    [10, 7, 38],
    [58, 82, 139],
    [33, 145, 140],
    [94, 201, 98],
    [253, 231, 37],
  ],
  magma: [
    [3, 3, 14],
    [59, 15, 112],
    [140, 41, 129],
    [222, 73, 104],
    [254, 159, 109],
    [252, 253, 191],
  ],
  coolwarm: [
    [32, 82, 218],
    [88, 146, 240],
    [7, 9, 16],
    [240, 108, 72],
    [228, 36, 44],
  ],
};

function mapValueToRGB(t, cmap) {
  const v = clamp01(t);
  if (cmap === "grayscale") {
    const g = Math.round(v * 255);
    return [g, g, g];
  }
  if (cmap === "inverted") {
    const g = Math.round((1.0 - v) * 246);
    return [g, g, g];
  }
  const stops = COLORMAP_STOPS[cmap] || COLORMAP_STOPS.inferno;
  return interpolateStops(stops, v);
}

function applyStretch(normVal, stretchMode, contrast, bias) {
  let x = (normVal - bias) * contrast;
  x = Math.max(0.0, x);
  if (stretchMode === "asinh") {
    return clamp01(Math.asinh(x * 6.0) / Math.asinh(6.0));
  }
  if (stretchMode === "sqrt") {
    return clamp01(Math.sqrt(x));
  }
  if (stretchMode === "log") {
    return clamp01(Math.log10(1.0 + x * 19.0) / Math.log10(20.0));
  }
  return clamp01(x);
}

// URL Sync & Utilities

function degToHMS(raDeg) {
  const totalSec = (((raDeg % 360) + 360) % 360) * 240.0;
  const h = Math.floor(totalSec / 3600);
  const m = Math.floor((totalSec % 3600) / 60);
  const s = (totalSec % 60).toFixed(1);
  return `${String(h).padStart(2, "0")}h${String(m).padStart(2, "0")}m${String(s).padStart(4, "0")}s`;
}

function degToDMS(decDeg) {
  const sign = decDeg >= 0 ? "+" : "-";
  const absDec = Math.abs(decDeg);
  const d = Math.floor(absDec);
  const m = Math.floor((absDec - d) * 60);
  const s = Math.round(((absDec - d) * 60 - m) * 60);
  return `${sign}${String(d).padStart(2, "0")}d${String(m).padStart(2, "0")}m${String(s).padStart(2, "0")}s`;
}

function readUrlParameters() {
  const params = new URLSearchParams(window.location.search);
  if (params.has("ra") && params.has("dec")) {
    const r = parseFloat(params.get("ra"));
    const d = parseFloat(params.get("dec"));
    if (!isNaN(r) && !isNaN(d)) {
      state.ra = r;
      state.dec = d;
      const raEl = document.getElementById("raInput");
      const decEl = document.getElementById("decInput");
      if (raEl) raEl.value = r.toFixed(6);
      if (decEl) decEl.value = (d >= 0 ? "+" : "") + d.toFixed(6);
    }
  }
  if (params.has("fov")) {
    const f = parseFloat(params.get("fov"));
    if (!isNaN(f)) {
      state.fov = f;
      document.getElementById("fovSelect").value = String(f);
    }
  }
  if (params.has("band")) state.activeBand = params.get("band");
  if (params.has("mode")) state.viewMode = params.get("mode");
}

function syncUrlParameters() {
  const url = new URL(window.location.href);
  url.searchParams.set("ra", state.ra.toFixed(6));
  url.searchParams.set("dec", state.dec.toFixed(6));
  url.searchParams.set("fov", state.fov.toFixed(1));
  url.searchParams.set("band", state.activeBand);
  url.searchParams.set("mode", state.viewMode);
  window.history.replaceState({}, "", url.toString());
}

// Target Presets Dropdown & Detector Dock

function renderTargetSelectorStrip() {
  const container = document.getElementById("presetsListContainer");
  if (!container) return;
  container.innerHTML = "";

  state.targets.forEach((t, idx) => {
    const isMatched = Math.hypot(t.ra - state.ra, t.dec - state.dec) < 0.08;
    const item = document.createElement("div");
    item.className = `preset-item ${isMatched ? "active" : ""}`;
    item.innerHTML = `
      <div>
        <div style="display:flex;align-items:center;gap:8px;">
          <span style="font-family:var(--font-mono);font-size:10.5px;color:var(--accent-gold);">0${idx + 1}</span>
          <span style="font-weight:600;font-size:13px;color:var(--text-primary);">${t.name}</span>
          <span style="font-size:11px;color:var(--text-secondary);">· ${t.category}</span>
        </div>
        <div style="font-family:var(--font-mono);font-size:10.5px;color:var(--text-tertiary);margin-top:2px;">
          RA ${t.ra.toFixed(4)}° (${t.ra_hms}) · Dec ${t.dec >= 0 ? "+" : ""}${t.dec.toFixed(4)}° (${t.dec_dms})
        </div>
      </div>
      <div style="text-align:right;font-family:var(--font-mono);font-size:11px;color:var(--accent-cyan);">
        μ = ${t.total_pm_arcsecyr}"/yr
      </div>
    `;
    item.addEventListener("click", () => {
      document.getElementById("presetsDropdown").classList.add("hidden");
      selectTargetWithTransition(t);
    });
    container.appendChild(item);
  });
}

function selectTargetWithTransition(t) {
  const portal = document.getElementById("telescopePortal");
  const edCol = document.getElementById("editorialColumn");

  portal.classList.add("slewing");
  edCol.style.animation = "none";
  void edCol.offsetWidth;
  edCol.style.animation = "fadeSlideUp 0.65s var(--ease-spring)";

  state.ra = t.ra;
  state.dec = t.dec;
  state.fov = t.default_fov_arcmin || 4.5;
  document.getElementById("raInput").value = t.ra.toFixed(6);
  document.getElementById("decInput").value = (t.dec >= 0 ? "+" : "") + t.dec.toFixed(6);
  document.getElementById("fovSelect").value = String(state.fov);

  if (t.recommended_band && t.recommended_band !== "ALL") {
    state.activeBand = t.recommended_band;
    state.viewMode = "SINGLE_BAND";
  } else {
    state.viewMode = "SINGLE_BAND";
    state.activeBand = "SPHEREx-D2";
  }

  updateViewModeTabs();
  updateAmbientAura();
  loadSkyCube().finally(() => {
    setTimeout(() => portal.classList.remove("slewing"), 220);
  });
}

function renderDetectorBandStrip() {
  const container = document.getElementById("detectorBandStrip");
  container.innerHTML = "";

  // 1. All-Band Stack (SNR Boost) pill
  const allBtn = document.createElement("button");
  allBtn.type = "button";
  allBtn.className = `dock-band-btn ${state.activeBand === "ALL" ? "active" : ""}`;
  allBtn.style.setProperty("--b-color", "#38efdf");
  allBtn.title = "Stack all 6 detector arrays (D1–D6, 0.75–5.00 μm) to maximize Signal-to-Noise Ratio";
  allBtn.innerHTML = `<span class="dot"></span><span>ALL</span>`;
  allBtn.addEventListener("click", () => {
    state.activeBand = "ALL";
    updateAmbientAura();
    renderDetectorBandStrip();
    renderSkyCanvas();
    syncUrlParameters();
  });
  container.appendChild(allBtn);

  // 2. Individual D1..D6 detector pills
  state.bands.forEach((b) => {
    const isActive = state.activeBand === b.id;
    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = `dock-band-btn ${isActive ? "active" : ""}`;
    btn.style.setProperty("--b-color", b.color);
    btn.title = `${b.name}: ${b.wavelength_min_um}–${b.wavelength_max_um} μm (R=${b.resolving_power})`;
    btn.innerHTML = `<span class="dot"></span><span>${b.short}</span>`;
    btn.addEventListener("click", () => {
      state.activeBand = b.id;
      if (state.viewMode === "RGB_COMPOSITE") {
        state.viewMode = "SINGLE_BAND";
        updateViewModeTabs();
      }
      if (state.viewMode === "RAW_FITS") {
        loadRawCutout();
      }
      updateAmbientAura();
      renderDetectorBandStrip();
      renderSkyCanvas();
      syncUrlParameters();
    });
    container.appendChild(btn);
  });
}

function updateViewModeTabs() {
  document.querySelectorAll(".dock-mode-btn").forEach((btn) => {
    const m = btn.getAttribute("data-mode");
    btn.classList.toggle("active", m === state.viewMode);
  });
}

function renderEpochTimelineButtons() {
  const container = document.getElementById("epochTimelineButtons");
  container.innerHTML = "";
  if (!state.skyCube || !state.skyCube.epochs) return;

  state.skyCube.epochs.forEach((ep, idx) => {
    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = `epoch-pill ${idx === state.activeEpoch ? "active" : ""}`;
    btn.textContent = `E${idx + 1}`;
    btn.title = `${ep.label} (MJD ${ep.mjd})`;
    btn.addEventListener("click", () => {
      state.activeEpoch = idx;
      renderEpochTimelineButtons();
      renderSkyCanvas();
      renderLightCurveCanvas();
    });
    container.appendChild(btn);
  });
}

// Sky Viewport Renderer (Time-Lapse, Split Curtain, Dipole, RGB)

function renderEpochToOffscreen(ctx, N, epIdx, modeOverride = null) {
  const cube = state.skyCube;
  const mode = modeOverride || state.viewMode;
  const imgData = ctx.createImageData(N, N);
  const data = imgData.data;

  if (mode === "RGB_COMPOSITE") {
    const f = cube.frames_by_epoch[epIdx];
    const d1 = f["SPHEREx-D1"];
    const d2 = f["SPHEREx-D2"];
    const d3 = f["SPHEREx-D3"];
    const d4 = f["SPHEREx-D4"];
    const d5 = f["SPHEREx-D5"];
    const d6 = f["SPHEREx-D6"];
    const maxVal = 8.0;
    for (let y = 0; y < N; y++) {
      for (let x = 0; x < N; x++) {
        const rRaw = 0.4 * d5[y][x] + 0.6 * d6[y][x];
        const gRaw = 0.5 * d3[y][x] + 0.5 * d4[y][x];
        const bRaw = 0.5 * d1[y][x] + 0.5 * d2[y][x];
        const idx = (y * N + x) * 4;
        data[idx] = Math.round(applyStretch(rRaw / maxVal, state.stretch, state.contrast, state.bias) * 255);
        data[idx + 1] = Math.round(applyStretch(gRaw / maxVal, state.stretch, state.contrast, state.bias) * 255);
        data[idx + 2] = Math.round(applyStretch(bRaw / maxVal, state.stretch, state.contrast, state.bias) * 255);
        data[idx + 3] = 255;
      }
    }
  } else if (mode === "DIFFERENCE") {
    const mat = cube.difference_image.by_band[state.activeBand] || cube.difference_image.coadd;
    let maxAbs = 1.0;
    for (let y = 0; y < N; y++) {
      for (let x = 0; x < N; x++) {
        const a = Math.abs(mat[y][x]);
        if (a > maxAbs) maxAbs = a;
      }
    }
    const scale = maxAbs * 0.62;
    for (let y = 0; y < N; y++) {
      for (let x = 0; x < N; x++) {
        const v = mat[y][x] / scale;
        const signed = Math.asinh(v * state.contrast * 2.5) / Math.asinh(2.5);
        const t = clamp01(0.5 + 0.5 * signed);
        const rgb = mapValueToRGB(t, "coolwarm");
        const idx = (y * N + x) * 4;
        data[idx] = rgb[0];
        data[idx + 1] = rgb[1];
        data[idx + 2] = rgb[2];
        data[idx + 3] = 255;
      }
    }
  } else {
    let mat;
    if (mode === "TEMPORAL_COADD") {
      mat = cube.temporal_coadd.by_band[state.activeBand] || cube.temporal_coadd.coadd;
    } else if (state.activeBand === "ALL") {
      mat = cube.coadd_by_epoch[epIdx];
    } else {
      mat = cube.frames_by_epoch[epIdx][state.activeBand] || cube.coadd_by_epoch[epIdx];
    }
    const maxVal = 8.5;
    for (let y = 0; y < N; y++) {
      for (let x = 0; x < N; x++) {
        const norm = mat[y][x] / maxVal;
        const stretched = applyStretch(norm, state.stretch, state.contrast, state.bias);
        const rgb = mapValueToRGB(stretched, state.colormap);
        const idx = (y * N + x) * 4;
        data[idx] = rgb[0];
        data[idx + 1] = rgb[1];
        data[idx + 2] = rgb[2];
        data[idx + 3] = 255;
      }
    }
  }
  ctx.putImageData(imgData, 0, 0);
}

function renderSkyCanvas() {
  const cube = state.skyCube;
  if (!cube) return;

  const canvas = document.getElementById("skyCanvas");
  const ctx = canvas.getContext("2d");
  const W = canvas.width;
  const H = canvas.height;
  const N = cube.query.grid_size;

  offscreenA.width = N;
  offscreenA.height = N;
  offscreenB.width = N;
  offscreenB.height = N;

  const epIdx = Math.max(0, Math.min(state.activeEpoch, cube.epochs.length - 1));
  const lastEpIdx = cube.epochs.length - 1;

  ctx.save();
  ctx.clearRect(0, 0, W, H);
  ctx.imageSmoothingEnabled = true;
  ctx.imageSmoothingQuality = "high";

  if (state.viewMode === "RAW_FITS" && state.rawCutout && state.rawCutout.image_mjy_sr) {
    const cutN = state.rawCutout.cutout_size;
    offscreenA.width = cutN;
    offscreenA.height = cutN;
    const imgData = ctxA.createImageData(cutN, cutN);
    const data = imgData.data;
    const rawMat = state.zodiSubtract
      ? state.rawCutout.image_mjy_sr
      : state.rawCutout.raw_image_mjy_sr || state.rawCutout.image_mjy_sr;

    let minV = Infinity;
    let maxV = -Infinity;
    for (let y = 0; y < cutN; y++) {
      for (let x = 0; x < cutN; x++) {
        const v = rawMat[y][x];
        if (v < minV) minV = v;
        if (v > maxV) maxV = v;
      }
    }
    const span = Math.max(1e-3, maxV - minV);
    for (let y = 0; y < cutN; y++) {
      for (let x = 0; x < cutN; x++) {
        const norm = (rawMat[y][x] - minV) / span;
        const stretched = applyStretch(norm, state.stretch, state.contrast, state.bias);
        const rgb = mapValueToRGB(stretched, state.colormap);
        const idx = (y * cutN + x) * 4;
        data[idx] = rgb[0];
        data[idx + 1] = rgb[1];
        data[idx + 2] = rgb[2];
        data[idx + 3] = 255;
      }
    }
    ctxA.putImageData(imgData, 0, 0);
    // Preserve crisp native 6.15"/px detector pixels
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(offscreenA, 0, 0, W, H);
  } else if (state.viewMode === "SPLIT_CURTAIN") {
    renderEpochToOffscreen(ctxA, N, 0, "SINGLE_BAND");
    renderEpochToOffscreen(ctxB, N, lastEpIdx, "SINGLE_BAND");

    const splitX = Math.round(state.splitRatio * W);

    ctx.save();
    ctx.beginPath();
    ctx.rect(0, 0, splitX, H);
    ctx.clip();
    ctx.drawImage(offscreenA, 0, 0, W, H);
    ctx.restore();

    ctx.save();
    ctx.beginPath();
    ctx.rect(splitX, 0, W - splitX, H);
    ctx.clip();
    ctx.drawImage(offscreenB, 0, 0, W, H);
    ctx.restore();

    // Sleek Curtain Line & Pill Handle
    ctx.strokeStyle = "rgba(255, 255, 255, 0.75)";
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(splitX, 0);
    ctx.lineTo(splitX, H);
    ctx.stroke();

    ctx.fillStyle = "rgba(10, 12, 18, 0.9)";
    ctx.strokeStyle = "rgba(255, 255, 255, 0.8)";
    ctx.beginPath();
    ctx.roundRect(splitX - 15, H / 2 - 14, 30, 28, 14);
    ctx.fill();
    ctx.stroke();

    ctx.fillStyle = "#fff";
    ctx.font = "600 10px 'JetBrains Mono', monospace";
    ctx.fillText("↔", splitX - 5, H / 2 + 4);
  } else {
    renderEpochToOffscreen(ctxA, N, epIdx, state.viewMode);
    ctx.drawImage(offscreenA, 0, 0, W, H);
  }
  ctx.restore();

  const scalePx = W / N;
  const ep = cube.epochs[epIdx];
  const ep0 = cube.epochs[0];
  const epLast = cube.epochs[lastEpIdx];
  const isTimeLapseMode = state.viewMode === "SINGLE_BAND" || state.viewMode === "RGB_COMPOSITE";
  const epochFraction =
    isTimeLapseMode && epLast.decimal_year > ep0.decimal_year
      ? (ep.decimal_year - ep0.decimal_year) / (epLast.decimal_year - ep0.decimal_year)
      : 0.0;

  // 1. Subtle WCS Grid
  if (state.showWcsGrid) {
    ctx.save();
    ctx.strokeStyle = "rgba(255, 255, 255, 0.06)";
    ctx.lineWidth = 1;
    for (let i = 1; i < 4; i++) {
      const p = (W / 4) * i;
      ctx.beginPath();
      ctx.moveTo(p, 0);
      ctx.lineTo(p, H);
      ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(0, p);
      ctx.lineTo(W, p);
      ctx.stroke();
    }
    ctx.restore();
  }

  // 2. Contamination FLAGS Bitmask
  if (state.showFlagsOverlay) {
    if (state.viewMode === "RAW_FITS" && state.rawCutout && state.rawCutout.flags) {
      const cutN = state.rawCutout.cutout_size;
      const cutScale = W / cutN;
      ctx.save();
      for (let y = 0; y < cutN; y++) {
        for (let x = 0; x < cutN; x++) {
          const mask = state.rawCutout.flags[y][x];
          if (mask & ((1 << 0) | (1 << 10) | (1 << 22))) {
            ctx.strokeStyle = "#ff5e4d";
            ctx.lineWidth = 1.4;
            ctx.strokeRect(x * cutScale, y * cutScale, cutScale, cutScale);
          }
        }
      }
      ctx.restore();
    } else {
      const flagsMat = cube.flags_by_epoch[epIdx];
      ctx.save();
      for (let y = 0; y < N; y++) {
        for (let x = 0; x < N; x++) {
          const mask = flagsMat[y][x];
          if (mask & ((1 << 0) | (1 << 10) | (1 << 22))) {
            ctx.strokeStyle = "#ff5e4d";
            ctx.lineWidth = 1.4;
            ctx.strokeRect(x * scalePx, y * scalePx, scalePx, scalePx);
          }
        }
      }
      ctx.restore();
    }
  }

  // 3. Gaia DR3 Counterparts & Proper Motion Vectors
  if (state.showGaiaOverlay && cube.catalog_sources) {
    ctx.save();
    cube.catalog_sources.forEach((s) => {
      const cx = (s.x_px + 0.5) * scalePx;
      const cy = (s.y_px + 0.5) * scalePx;

      if (s.g_mag < 21.5) {
        ctx.strokeStyle = "rgba(52, 211, 153, 0.55)";
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.arc(cx, cy, 6.5, 0, Math.PI * 2);
        ctx.stroke();
      }

      if (s.total_pm_arcsec_yr > 0.25) {
        const vecScale = 4.5 * state.motionExaggeration;
        const ex = cx + s.pm_dx_px_yr * scalePx * vecScale;
        const ey = cy + s.pm_dy_px_yr * scalePx * vecScale;

        ctx.strokeStyle = "#34d399";
        ctx.fillStyle = "#34d399";
        ctx.lineWidth = 1.6;
        ctx.beginPath();
        ctx.moveTo(cx, cy);
        ctx.lineTo(ex, ey);
        ctx.stroke();

        const angle = Math.atan2(ey - cy, ex - cx);
        ctx.beginPath();
        ctx.moveTo(ex, ey);
        ctx.lineTo(ex - 6 * Math.cos(angle - 0.45), ey - 6 * Math.sin(angle - 0.45));
        ctx.lineTo(ex - 6 * Math.cos(angle + 0.45), ey - 6 * Math.sin(angle + 0.45));
        ctx.closePath();
        ctx.fill();
      }
    });
    ctx.restore();
  }

  // 4. Minimalist Corner Reticles for Moving Dipoles & Variables (Dynamic Astrometric Tracking)
  if (state.showDetectionsOverlay && cube.detections) {
    ctx.save();
    ctx.font = "500 10px 'JetBrains Mono', monospace";
    cube.detections.forEach((d) => {
      const shiftX = isTimeLapseMode ? (d.dx_px || 0) * epochFraction : 0;
      const shiftY = isTimeLapseMode ? (d.dy_px || 0) * epochFraction : 0;
      const x = (d.x_px + shiftX + 0.5) * scalePx;
      const y = (d.y_px + shiftY + 0.5) * scalePx;
      const isMoving = d.signature.includes("DIPOLE");
      const color = isMoving ? "#f5b041" : "#ff5e4d";

      ctx.strokeStyle = color;
      ctx.lineWidth = 1.4;
      const r = 14;
      const b = 4;
      ctx.beginPath();
      ctx.moveTo(x - r, y - r + b); ctx.lineTo(x - r, y - r); ctx.lineTo(x - r + b, y - r);
      ctx.moveTo(x + r - b, y - r); ctx.lineTo(x + r, y - r); ctx.lineTo(x + r, y - r + b);
      ctx.moveTo(x + r, y + r - b); ctx.lineTo(x + r, y + r); ctx.lineTo(x + r - b, y + r);
      ctx.moveTo(x - r + b, y + r); ctx.lineTo(x - r, y + r); ctx.lineTo(x - r, y + r - b);
      ctx.stroke();
    });
    ctx.restore();
  }

  // 5. Active Spectrophotometer Target Ring (with Dynamic Moving-Star Astrometric Lock)
  ctx.save();
  let ax = (state.apertureX + 0.5) * scalePx;
  let ay = (state.apertureY + 0.5) * scalePx;
  const sel = state.selectedSource;

  if (
    sel &&
    sel.total_pm_arcsec_yr > 0.15 &&
    isTimeLapseMode &&
    Math.hypot(sel.x_px - state.apertureX, sel.y_px - state.apertureY) < 4.5
  ) {
    const dtYr = (ep.decimal_year - ep0.decimal_year) * state.motionExaggeration;
    const refX = (sel.x_px + 0.5) * scalePx;
    const refY = (sel.y_px + 0.5) * scalePx;
    ax = (sel.x_px + (sel.pm_dx_px_yr || 0) * dtYr + 0.5) * scalePx;
    ay = (sel.y_px + (sel.pm_dy_px_yr || 0) * dtYr + 0.5) * scalePx;

    // Subtle trajectory trail from Epoch 1 to Current Epoch
    if (Math.hypot(ax - refX, ay - refY) > 0.8) {
      ctx.strokeStyle = "rgba(56, 239, 223, 0.55)";
      ctx.setLineDash([2, 3]);
      ctx.lineWidth = 1.3;
      ctx.beginPath();
      ctx.moveTo(refX, refY);
      ctx.lineTo(ax, ay);
      ctx.stroke();
      ctx.setLineDash([]);

      ctx.fillStyle = "rgba(56, 239, 223, 0.7)";
      ctx.beginPath();
      ctx.arc(refX, refY, 2.2, 0, Math.PI * 2);
      ctx.fill();

      const realShiftArcsec = sel.total_pm_arcsec_yr * Math.abs(ep.decimal_year - ep0.decimal_year);
      ctx.fillStyle = "#38efdf";
      ctx.font = "600 9.5px 'JetBrains Mono', monospace";
      ctx.fillText(`Δθ ${realShiftArcsec.toFixed(2)}"`, ax + 15, ay - 6);
    }
  }

  ctx.strokeStyle = "#38efdf";
  ctx.lineWidth = 1.6;
  ctx.beginPath();
  ctx.arc(ax, ay, 11, 0, Math.PI * 2);
  ctx.stroke();
  ctx.restore();

  // 6. Precision Astronomical Angular Scale Bar (Bottom-Left of Portal)
  ctx.save();
  const barArcmin = state.fov <= 3.0 ? 0.5 : 1.0;
  const barLabel = barArcmin === 0.5 ? `30"` : `1' (60")`;
  const barPx = (barArcmin / state.fov) * W;
  const bx0 = 18;
  const by0 = H - 48;

  ctx.strokeStyle = "rgba(248, 250, 252, 0.65)";
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(bx0, by0 - 3);
  ctx.lineTo(bx0, by0 + 3);
  ctx.moveTo(bx0, by0);
  ctx.lineTo(bx0 + barPx, by0);
  ctx.moveTo(bx0 + barPx, by0 - 3);
  ctx.lineTo(bx0 + barPx, by0 + 3);
  ctx.stroke();

  ctx.fillStyle = "rgba(248, 250, 252, 0.72)";
  ctx.font = "500 9.5px 'JetBrains Mono', monospace";
  ctx.fillText(barLabel, bx0 + 4, by0 - 6);
  ctx.restore();

  // Update Portal Pills & Mode Legend Pill
  const bandObj = state.bands.find((b) => b.id === state.activeBand) || state.bands[0];
  const legendPill = document.getElementById("dipoleLegendPill");

  if (state.viewMode === "DIFFERENCE") {
    document.getElementById("hudEpochTitle").textContent = "Dipole Difference (E6 − E1)";
    document.getElementById("hudBandLabel").textContent = bandObj.short_label;
    if (legendPill) {
      legendPill.innerHTML = `<span style="color:#60a5fa;">● 2025 E1 (Old)</span><span style="color:var(--text-tertiary);">───</span><span style="color:#f87171;">● 2026 E6 (New)</span>`;
      legendPill.classList.remove("hidden");
    }
  } else if (state.viewMode === "RAW_FITS") {
    const cutInfo = state.rawCutout;
    document.getElementById("hudEpochTitle").textContent = cutInfo
      ? `RAW L2 FITS · 6.15"/px (${cutInfo.cutout_size}×${cutInfo.cutout_size})`
      : "Loading Raw IRSA FITS...";
    document.getElementById("hudBandLabel").textContent = bandObj.short_label;
    if (legendPill && cutInfo) {
      const paTxt = cutInfo.detector_pa_deg ? ` · WCS N↑ (Roll ${cutInfo.detector_pa_deg}°)` : "";
      legendPill.innerHTML = `<span style="color:var(--accent-cyan);">IRSA L2 FITS</span><span>·</span><span>${cutInfo.obs_id}${paTxt}</span>`;
      legendPill.classList.remove("hidden");
    }
  } else if (state.viewMode === "SPLIT_CURTAIN") {
    document.getElementById("hudEpochTitle").textContent = "Epoch 1 ◄ | ► Epoch 6";
    document.getElementById("hudBandLabel").textContent = bandObj.short_label;
    if (legendPill) legendPill.classList.add("hidden");
  } else if (state.viewMode === "RGB_COMPOSITE") {
    document.getElementById("hudEpochTitle").textContent = `Epoch ${epIdx + 1} · ${ep.iso_date.slice(0, 7)}`;
    document.getElementById("hudBandLabel").textContent = "Thermal RGB (D6/D4/D1)";
    if (legendPill) legendPill.classList.add("hidden");
  } else {
    document.getElementById("hudEpochTitle").textContent = `Epoch ${epIdx + 1} · ${ep.iso_date.slice(0, 7)}`;
    document.getElementById("hudBandLabel").textContent = bandObj.short_label;
    if (legendPill) legendPill.classList.add("hidden");
  }
  document.getElementById("hudFovLabel").textContent = `${state.fov.toFixed(1)}'`;
}

function renderOpticalLoupe(mainCanvas) {
  const loupe = document.getElementById("loupeCanvas");
  if (!loupe || !state.skyCube || state.hoverX === null) return;
  const lctx = loupe.getContext("2d");
  const LW = loupe.width;
  const LH = loupe.height;
  const N = state.skyCube.query.grid_size;
  const scalePx = mainCanvas.width / N;

  const focusX = state.hoverX * scalePx;
  const focusY = state.hoverY * scalePx;
  const zoom = 3.2;
  const srcW = LW / zoom;
  const srcH = LH / zoom;

  lctx.clearRect(0, 0, LW, LH);
  lctx.save();
  lctx.imageSmoothingEnabled = false;
  lctx.drawImage(mainCanvas, focusX - srcW / 2, focusY - srcH / 2, srcW, srcH, 0, 0, LW, LH);

  lctx.strokeStyle = "rgba(56, 239, 223, 0.55)";
  lctx.lineWidth = 1;
  lctx.beginPath();
  lctx.moveTo(LW / 2 - 10, LH / 2); lctx.lineTo(LW / 2 + 10, LH / 2);
  lctx.moveTo(LW / 2, LH / 2 - 10); lctx.lineTo(LW / 2, LH / 2 + 10);
  lctx.stroke();
  lctx.restore();
}

// 102-Band Near-IR Spectrum Renderer

function triggerSpectrumAnimation() {
  if (state.specAnimFrame) cancelAnimationFrame(state.specAnimFrame);
  const start = performance.now();
  const duration = 520;

  function step(now) {
    const t = Math.min(1.0, (now - start) / duration);
    // Cubic ease-out
    state.specAnimProgress = 1.0 - Math.pow(1.0 - t, 3);
    renderSpectrumCanvas();
    if (t < 1.0) {
      state.specAnimFrame = requestAnimationFrame(step);
    }
  }
  state.specAnimFrame = requestAnimationFrame(step);
}

function renderSpectrumCanvas() {
  const specData = state.spectrumData;
  if (!specData || !specData.spectrum_102) return;

  const canvas = document.getElementById("spectrumCanvas");
  const ctx = canvas.getContext("2d");
  const W = canvas.width;
  const H = canvas.height;
  ctx.clearRect(0, 0, W, H);

  const padL = 36;
  const padR = 12;
  const padT = 20;
  const padB = 26;
  const plotW = W - padL - padR;
  const plotH = H - padT - padB;

  const wl = specData.spectrum_102.wavelengths_um;
  const flux = specData.spectrum_102.flux_mjy;
  const bIds = specData.spectrum_102.band_ids;

  const minWl = 0.72;
  const maxWl = 5.05;
  const maxFlux = Math.max(...flux) * 1.18 || 10.0;
  const animP = state.specAnimProgress;

  const xOfWl = (w) => padL + ((w - minWl) / (maxWl - minWl)) * plotW;
  const yOfFlux = (f) => padT + plotH - ((f * animP) / maxFlux) * plotH;

  // Band labels at top
  state.bands.forEach((b) => {
    const x0 = xOfWl(b.wavelength_min_um);
    const x1 = xOfWl(b.wavelength_max_um);
    ctx.fillStyle = `${b.color}0d`;
    ctx.fillRect(x0, padT, x1 - x0, plotH);
    ctx.fillStyle = b.color;
    ctx.font = "600 9px 'JetBrains Mono', monospace";
    ctx.fillText(b.short, (x0 + x1) / 2 - 6, padT + 10);
  });

  // Molecular Absorption Lines
  if (state.showMolecularLines && specData.spectral_features) {
    ctx.save();
    specData.spectral_features.forEach((feat) => {
      const x = xOfWl(feat.wavelength_um);
      if (x < padL || x > W - padR) return;
      ctx.strokeStyle = "rgba(255, 255, 255, 0.08)";
      ctx.setLineDash([2, 3]);
      ctx.beginPath();
      ctx.moveTo(x, padT + 14);
      ctx.lineTo(x, padT + plotH);
      ctx.stroke();
    });
    ctx.restore();
  }

  // Glowing Area Fill Under 102-Band Curve
  const grad = ctx.createLinearGradient(0, padT, 0, padT + plotH);
  grad.addColorStop(0, "rgba(56, 239, 223, 0.25)");
  grad.addColorStop(1, "rgba(56, 239, 223, 0.0)");

  ctx.beginPath();
  ctx.moveTo(xOfWl(wl[0]), padT + plotH);
  for (let i = 0; i < wl.length; i++) {
    ctx.lineTo(xOfWl(wl[i]), yOfFlux(flux[i]));
  }
  ctx.lineTo(xOfWl(wl[wl.length - 1]), padT + plotH);
  ctx.closePath();
  ctx.fillStyle = grad;
  ctx.fill();

  // Smooth Connecting Curve
  ctx.strokeStyle = "rgba(248, 250, 252, 0.65)";
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  for (let i = 0; i < wl.length; i++) {
    const x = xOfWl(wl[i]);
    const y = yOfFlux(flux[i]);
    if (i === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  }
  ctx.stroke();

  const bandColorMap = {};
  state.bands.forEach((b) => (bandColorMap[b.id] = b.color));

  for (let i = 0; i < wl.length; i++) {
    const x = xOfWl(wl[i]);
    const y = yOfFlux(flux[i]);
    const isHov = state.hoveredSpecChannel === i;
    ctx.fillStyle = isHov ? "#ffffff" : bandColorMap[bIds[i]] || "#38efdf";
    ctx.beginPath();
    ctx.arc(x, y, isHov ? 4.2 : 2.0, 0, Math.PI * 2);
    ctx.fill();
  }

  // Overlay Real 2MASS (J, H, Ks) & WISE (W1, W2) Photometric Anchor Diamonds
  const anchors = specData.spectrum_102.photometric_anchors || [];
  if (anchors.length > 0) {
    ctx.save();
    anchors.forEach((a) => {
      const ax = xOfWl(a.wavelength_um);
      const ay = yOfFlux(a.flux_mjy);
      if (ax < padL || ax > W - padR || ay < padT || ay > padT + plotH) return;
      ctx.fillStyle = "#f5b041";
      ctx.strokeStyle = "#06080e";
      ctx.lineWidth = 1.2;
      ctx.beginPath();
      ctx.moveTo(ax, ay - 4.2);
      ctx.lineTo(ax + 4.2, ay);
      ctx.lineTo(ax, ay + 4.2);
      ctx.lineTo(ax - 4.2, ay);
      ctx.closePath();
      ctx.fill();
      ctx.stroke();
    });
    ctx.restore();
  }

  // Wavelength axis
  ctx.fillStyle = "rgba(248, 250, 252, 0.45)";
  ctx.font = "9.5px 'JetBrains Mono', monospace";
  [1.0, 2.0, 3.0, 4.0, 5.0].forEach((wTick) => {
    const x = xOfWl(wTick);
    ctx.fillText(`${wTick.toFixed(1)}μm`, x - 14, H - 7);
  });
}

// Light Curve & Archive Renderers

function renderLightCurveCanvas() {
  const specData = state.spectrumData;
  if (!specData || !specData.light_curve) return;

  const canvas = document.getElementById("lightCurveCanvas");
  const ctx = canvas.getContext("2d");
  const W = canvas.width;
  const H = canvas.height;
  ctx.clearRect(0, 0, W, H);

  const lc = specData.light_curve;
  if (!lc.length) return;

  const padL = 38;
  const padR = 14;
  const padT = 20;
  const padB = 26;
  const plotW = W - padL - padR;
  const plotH = H - padT - padB;
  const xOfIdx = (i) => padL + (i / Math.max(1, lc.length - 1)) * plotW;

  const activeX = xOfIdx(Math.min(state.activeEpoch, lc.length - 1));
  ctx.strokeStyle = "rgba(56, 239, 223, 0.45)";
  ctx.setLineDash([3, 3]);
  ctx.beginPath();
  ctx.moveTo(activeX, padT);
  ctx.lineTo(activeX, padT + plotH);
  ctx.stroke();
  ctx.setLineDash([]);

  if (state.lightCurveMode === "FLUX") {
    const keys = [
      { id: "SPHEREx-D1", label: "D1", color: "#60a5fa" },
      { id: "SPHEREx-D2", label: "D2", color: "#34d399" },
      { id: "SPHEREx-D4", label: "D4", color: "#facc15" },
      { id: "SPHEREx-D6", label: "D6", color: "#f87171" },
    ];
    let maxF = 1.0;
    let minF = Infinity;
    lc.forEach((pt) => {
      keys.forEach((k) => {
        const v = pt.band_fluxes_mjy[k.id] || 0;
        if (v > maxF) maxF = v;
        if (v < minF) minF = v;
      });
    });
    minF = Math.max(0, minF * 0.8);
    maxF = maxF * 1.15;
    if (maxF - minF < 1e-3) maxF = minF + 1.0;

    const yOfF = (f) => padT + plotH - ((f - minF) / (maxF - minF)) * plotH;

    keys.forEach((k, idx) => {
      ctx.strokeStyle = k.color;
      ctx.lineWidth = 1.8;
      ctx.beginPath();
      lc.forEach((pt, i) => {
        const x = xOfIdx(i);
        const y = yOfF(pt.band_fluxes_mjy[k.id]);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      });
      ctx.stroke();

      lc.forEach((pt, i) => {
        const x = xOfIdx(i);
        const y = yOfF(pt.band_fluxes_mjy[k.id]);
        ctx.fillStyle = k.color;
        ctx.beginPath();
        ctx.arc(x, y, i === state.activeEpoch ? 4.2 : 2.6, 0, Math.PI * 2);
        ctx.fill();
      });

      ctx.fillStyle = k.color;
      ctx.font = "600 9.5px 'JetBrains Mono', monospace";
      ctx.fillText(k.label, padL + idx * 45, 12);
    });
  } else {
    let maxS = 0.5;
    lc.forEach((pt) => {
      if (pt.total_offset_arcsec > maxS) maxS = pt.total_offset_arcsec;
    });
    maxS *= 1.2;
    const yOfS = (s) => padT + plotH - (s / maxS) * plotH;

    ctx.strokeStyle = "#f5b041";
    ctx.lineWidth = 2.0;
    ctx.beginPath();
    lc.forEach((pt, i) => {
      const x = xOfIdx(i);
      const y = yOfS(pt.total_offset_arcsec);
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    });
    ctx.stroke();

    lc.forEach((pt, i) => {
      const x = xOfIdx(i);
      const y = yOfS(pt.total_offset_arcsec);
      ctx.fillStyle = "#f5b041";
      ctx.beginPath();
      ctx.arc(x, y, i === state.activeEpoch ? 4.4 : 2.8, 0, Math.PI * 2);
      ctx.fill();
    });
  }

  ctx.fillStyle = "rgba(248, 250, 252, 0.45)";
  ctx.font = "9.5px 'JetBrains Mono', monospace";
  lc.forEach((pt, i) => {
    const x = xOfIdx(i);
    ctx.fillText(pt.iso_date.slice(2, 7), x - 14, H - 7);
  });
}

function renderArchivePane() {
  const cube = state.skyCube;
  if (!cube) return;

  const detBox = document.getElementById("detectionsListContainer");
  detBox.innerHTML = (cube.detections || [])
    .map(
      (d) => `
      <div data-x="${d.x_px}" data-y="${d.y_px}" data-id="${d.source_id}" class="det-card" style="padding:8px 11px;border-radius:12px;background:rgba(255,255,255,0.04);border:1px solid var(--glass-border);margin-bottom:6px;cursor:pointer;display:flex;justify-content:space-between;align-items:center;">
        <div>
          <div style="font-size:12px;font-weight:600;color:var(--text-primary);">${d.name}</div>
          <div style="font-family:var(--font-mono);font-size:10px;color:var(--accent-gold);">${d.signature}</div>
        </div>
        <div style="font-family:var(--font-mono);font-size:11px;color:var(--accent-cyan);">μ=${d.proper_motion_arcsec_yr}"/yr</div>
      </div>
    `
    )
    .join("");

  detBox.querySelectorAll(".det-card").forEach((el) => {
    el.addEventListener("click", () => {
      state.apertureX = parseFloat(el.getAttribute("data-x"));
      state.apertureY = parseFloat(el.getAttribute("data-y"));
      renderSkyCanvas();
      loadPointSpectrum(state.apertureX, state.apertureY, el.getAttribute("data-id"));
    });
  });

  const links = cube.external_links || {};
  const archiveItems = [
    { label: "IRSA SPHEREx", url: links.irsa_spherex_explorer },
    { label: "WISE FinderChart", url: links.irsa_finderchart_wise },
    { label: "WiseView Blinker", url: links.wiseview_blinker },
    { label: "CSA NEOSSat", url: links.csa_neossat_cadc },
  ];
  document.getElementById("externalArchiveLinks").innerHTML = archiveItems
    .map(
      (a) => `
      <a href="${a.url}" target="_blank" rel="noopener" style="text-decoration:none;padding:7px 10px;border-radius:10px;background:rgba(255,255,255,0.03);border:1px solid var(--glass-border);color:var(--text-secondary);font-family:var(--font-mono);font-size:10.5px;display:flex;justify-content:space-between;">
        <span>${a.label}</span><span>↗</span>
      </a>
    `
    )
    .join("");

  const expBox = document.getElementById("irsaExposuresContainer");
  expBox.innerHTML = (cube.sample_irsa_exposures || [])
    .slice(0, 8)
    .map((e) => {
      const cutUrl = `${e.access_url}?center=${state.ra.toFixed(5)},${state.dec.toFixed(5)}d&size=${(state.fov / 60).toFixed(4)}`;
      return `
        <div style="padding:5px 8px;border-bottom:1px solid var(--glass-border);display:flex;justify-content:space-between;font-family:var(--font-mono);font-size:10px;color:var(--text-secondary);">
          <span><strong style="color:var(--text-primary);">${e.band}</strong> · ${e.iso_date.slice(0, 10)}</span>
          <a href="${cutUrl}" target="_blank" rel="noopener" style="color:var(--accent-cyan);text-decoration:none;">FITS ↗</a>
        </div>
      `;
    })
    .join("");
}

// API Client & Signature Reveal Controls

async function loadSkyCube() {
  const badge = document.getElementById("irsaLiveStatusBadge");
  badge.textContent = "Syncing IRSA...";

  const params = new URLSearchParams({
    ra: String(state.ra),
    dec: String(state.dec),
    fov: String(state.fov),
    zodi_subtract: String(state.zodiSubtract),
    motion_exaggeration: String(state.motionExaggeration),
  });

  try {
    const resp = await fetch(`/api/sky?${params.toString()}`);
    const cube = await resp.json();
    state.skyCube = cube;
    if (state.activeEpoch >= cube.epochs.length) state.activeEpoch = 0;

    badge.textContent = `${cube.query.total_irsa_exposures} IRSA Frames`;

    if (cube.matched_target) {
      const mt = cube.matched_target;
      const tIdx = state.targets.findIndex((x) => x.id === mt.id);
      const prop = cube.astrometric_propagation;
      const shiftSuffix =
        prop && prop.cumulative_shift_arcsec > 1.0
          ? ` · ${prop.ref_epoch}→${prop.obs_epoch} Δ=${prop.cumulative_shift_arcsec}"`
          : "";
      document.getElementById("storyBadge").textContent =
        `NASA DEMO TARGET 0${tIdx + 1} · ${mt.category.toUpperCase()}${shiftSuffix}`;
      document.getElementById("storyTitle").textContent = mt.name;
      document.getElementById("storyDescription").textContent = mt.scientific_notes;

      if (mt.id === "wise_0855") {
        document.getElementById("revealBtnLabel").textContent = "Reveal at 4.7 μm (Band D6)";
      } else if (mt.id === "barnards_star") {
        document.getElementById("revealBtnLabel").textContent = "Isolate Motion Dipole (Δ)";
      } else if (mt.id === "epsilon_indi_b") {
        document.getElementById("revealBtnLabel").textContent = "Composite Thermal RGB";
      } else {
        document.getElementById("revealBtnLabel").textContent = "Reveal Variability Curve";
      }
    } else {
      const nDet = cube.detections ? cube.detections.length : 0;
      document.getElementById("storyBadge").textContent =
        `DEEP SKY SECTOR · ${degToHMS(state.ra)} ${degToDMS(state.dec)}`;
      document.getElementById("storyTitle").textContent =
        `Sector J${state.ra.toFixed(2)}${state.dec >= 0 ? "+" : ""}${state.dec.toFixed(2)}`;
      document.getElementById("storyDescription").textContent =
        `Synchronized ${cube.query.total_irsa_exposures} SPHEREx calibrated frames across ${cube.epochs.length} epochs (FOV ${state.fov}'). Automated difference scanner isolated ${nDet} moving dipole / variable candidates in this field. Click any reticle on the sky portal to extract its 102-channel spectrum.`;
      document.getElementById("revealBtnLabel").textContent = "Isolate Sector Dipoles (Δ)";
    }

    if (cube.catalog_sources && cube.catalog_sources.length > 0) {
      state.apertureX = cube.catalog_sources[0].x_px;
      state.apertureY = cube.catalog_sources[0].y_px;
      state.selectedSource = cube.catalog_sources[0];
    }

    renderTargetSelectorStrip();
    renderDetectorBandStrip();
    renderEpochTimelineButtons();
    renderSkyCanvas();
    renderArchivePane();
    syncUrlParameters();

    if (state.viewMode === "RAW_FITS") {
      loadRawCutout();
    }

    await loadPointSpectrum(state.apertureX, state.apertureY);
  } catch (err) {
    console.error("Error loading sky cube:", err);
  }
}

async function loadRawCutout() {
  const bandParam = state.activeBand === "ALL" ? "SPHEREx-D2" : state.activeBand;
  const params = new URLSearchParams({
    ra: String(state.ra),
    dec: String(state.dec),
    fov: String(state.fov),
    band: bandParam,
  });
  try {
    const resp = await fetch(`/api/raw-cutout?${params.toString()}`);
    if (!resp.ok) return;
    state.rawCutout = await resp.json();
    renderSkyCanvas();
  } catch (err) {
    console.error("Failed to load raw IRSA cutout:", err);
  }
}

async function loadPointSpectrum(xPx, yPx, sourceId = null) {
  const params = new URLSearchParams({
    ra: String(state.ra),
    dec: String(state.dec),
    fov: String(state.fov),
    x: String(xPx),
    y: String(yPx),
    epoch: String(state.activeEpoch),
    motion_exaggeration: String(state.motionExaggeration),
  });
  if (sourceId) params.set("source_id", sourceId);

  try {
    const resp = await fetch(`/api/spectrum?${params.toString()}`);
    const spec = await resp.json();
    state.spectrumData = spec;

    const src = spec.source;
    state.selectedSource = src;
    document.getElementById("specSourceName").textContent = src.name;
    renderSkyCanvas();

    const optBadge = document.getElementById("opticalStatusBadge");
    if (src.g_mag > 22.0 && src.obj_class !== "ZODI_BACKGROUND") {
      optBadge.className = "insp-target-badge blindspot";
      optBadge.textContent = "Optical Blind Spot · IR Only";
    } else {
      optBadge.className = "insp-target-badge";
      optBadge.textContent = src.spectral_type;
    }

    document.getElementById("diagProperMotion").textContent =
      `${src.total_pm_arcsec_yr.toFixed(2)}"/yr`;
    document.getElementById("diagTemp").textContent = `${Math.round(src.temp_k)} K`;
    document.getElementById("diagIrRatio").textContent =
      String(spec.diagnostics.ir_to_optical_ratio_d6_d1);

    triggerSpectrumAnimation();
    renderLightCurveCanvas();
  } catch (err) {
    console.error("Failed to load point spectrum:", err);
  }
}

function revealTargetDiscoverySignature() {
  if (!state.skyCube) return;
  const mt = state.skyCube.matched_target;
  const tid = mt ? mt.id : "";

  if (tid === "barnards_star") {
    state.viewMode = "DIFFERENCE";
    state.activeBand = "SPHEREx-D2";
  } else if (tid === "wise_0855") {
    state.viewMode = "SINGLE_BAND";
    state.activeBand = "SPHEREx-D6";
  } else if (tid === "epsilon_indi_b") {
    state.viewMode = "RGB_COMPOSITE";
    state.activeBand = "SPHEREx-D2";
  } else if (tid === "rw_aurigae" || tid === "r_draconis") {
    state.viewMode = "DIFFERENCE";
    state.activeBand = "SPHEREx-D3";
    // Also switch inspector card to Time Curve
    document.querySelector('[data-card="CARD_LIGHTCURVE"]').click();
  } else {
    state.viewMode = "DIFFERENCE";
  }

  updateAmbientAura();
  updateViewModeTabs();
  renderDetectorBandStrip();
  renderSkyCanvas();
}

// Playback Loop & Event Bindings

function stepBlinkEpoch() {
  if (!state.skyCube || !state.skyCube.epochs) return;
  if (state.viewMode === "DIFFERENCE" || state.viewMode === "SPLIT_CURTAIN") return;

  const n = state.skyCube.epochs.length;
  if (n <= 1) return;

  let next = state.activeEpoch + state.blinkDirection;
  if (next >= n) {
    state.blinkDirection = -1;
    next = n - 2;
  } else if (next < 0) {
    state.blinkDirection = 1;
    next = 1;
  }
  state.activeEpoch = Math.max(0, Math.min(n - 1, next));

  renderEpochTimelineButtons();
  renderSkyCanvas();
  renderLightCurveCanvas();
}

function restartBlinkTimer() {
  if (state.timerId) {
    clearInterval(state.timerId);
    state.timerId = null;
  }
  if (state.isPlaying) {
    state.timerId = setInterval(stepBlinkEpoch, Math.round(1000 / state.fps));
    document.getElementById("playPauseIcon").textContent = "❚❚";
  } else {
    document.getElementById("playPauseIcon").textContent = "▶";
  }
}

function bindEvents() {
  // Magic Reveal Button
  document.getElementById("revealSignatureBtn").addEventListener("click", revealTargetDiscoverySignature);

  // Play / Pause
  document.getElementById("playPauseBtn").addEventListener("click", () => {
    state.isPlaying = !state.isPlaying;
    restartBlinkTimer();
  });

  // View Modes in Bottom Dock
  document.querySelectorAll(".dock-mode-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      state.viewMode = btn.getAttribute("data-mode");
      if (state.viewMode === "RAW_FITS") {
        loadRawCutout();
      }
      updateAmbientAura();
      updateViewModeTabs();
      renderSkyCanvas();
      syncUrlParameters();
    });
  });

  // Inspector Card Switcher (Spectrum | Time Curve | IRSA Archive)
  document.querySelectorAll("[data-card]").forEach((btn) => {
    btn.addEventListener("click", () => {
      const targetId = btn.getAttribute("data-card");
      document.querySelectorAll("[data-card]").forEach((b) => b.classList.toggle("active", b === btn));
      ["CARD_SPECTRUM", "CARD_LIGHTCURVE", "CARD_ARCHIVE"].forEach((id) => {
        document.getElementById(id).style.display = id === targetId ? "block" : "none";
      });
      if (targetId === "CARD_SPECTRUM") triggerSpectrumAnimation();
      if (targetId === "CARD_LIGHTCURVE") renderLightCurveCanvas();
    });
  });

  // Light Curve Mode Toggle
  const fluxBtn = document.getElementById("lcModeFluxBtn");
  const motionBtn = document.getElementById("lcModeMotionBtn");
  fluxBtn.addEventListener("click", () => {
    state.lightCurveMode = "FLUX";
    fluxBtn.classList.add("active");
    motionBtn.classList.remove("active");
    renderLightCurveCanvas();
  });
  motionBtn.addEventListener("click", () => {
    state.lightCurveMode = "MOTION";
    motionBtn.classList.add("active");
    fluxBtn.classList.remove("active");
    renderLightCurveCanvas();
  });

  // Layers Popover Toggle
  const layersPopover = document.getElementById("layersPopover");
  document.getElementById("toggleLayersPopoverBtn").addEventListener("click", (e) => {
    e.stopPropagation();
    layersPopover.classList.toggle("hidden");
  });
  window.addEventListener("click", (e) => {
    if (!layersPopover.contains(e.target)) {
      layersPopover.classList.add("hidden");
    }
  });

  const bindSwitch = (id, key, reload = false) => {
    const el = document.getElementById(id);
    el.addEventListener("click", () => {
      state[key] = !state[key];
      el.classList.toggle("on", state[key]);
      if (reload) loadSkyCube();
      else {
        renderSkyCanvas();
        renderSpectrumCanvas();
      }
    });
  };
  bindSwitch("swZodi", "zodiSubtract", true);
  bindSwitch("swGaia", "showGaiaOverlay", false);
  bindSwitch("swAnomalies", "showDetectionsOverlay", false);
  bindSwitch("swFlags", "showFlagsOverlay", false);
  bindSwitch("swGrid", "showWcsGrid", false);
  bindSwitch("swLines", "showMolecularLines", false);

  document.getElementById("colormapSelect").addEventListener("change", (e) => {
    state.colormap = e.target.value;
    renderSkyCanvas();
  });

  document.getElementById("motionScaleSelect").addEventListener("change", (e) => {
    state.motionExaggeration = parseFloat(e.target.value) || 1.0;
    loadSkyCube();
  });

  // Sky Canvas Interaction (Split Curtain + Smooth Loupe + Click Aperture)
  const skyCanvas = document.getElementById("skyCanvas");
  const portalLoupe = document.getElementById("portalLoupe");

  skyCanvas.addEventListener("mousedown", (e) => {
    const rect = skyCanvas.getBoundingClientRect();
    const relX = (e.clientX - rect.left) / rect.width;
    if (state.viewMode === "SPLIT_CURTAIN" && Math.abs(relX - state.splitRatio) < 0.07) {
      state.isDraggingSplit = true;
    }
  });
  window.addEventListener("mouseup", () => {
    state.isDraggingSplit = false;
  });

  skyCanvas.addEventListener("mousemove", (e) => {
    if (!state.skyCube) return;
    const rect = skyCanvas.getBoundingClientRect();
    const relX = clamp01((e.clientX - rect.left) / rect.width);
    const relY = clamp01((e.clientY - rect.top) / rect.height);

    if (state.isDraggingSplit && state.viewMode === "SPLIT_CURTAIN") {
      state.splitRatio = clamp01(Math.max(0.08, Math.min(0.92, relX)));
      renderSkyCanvas();
      return;
    }

    const N = state.skyCube.query.grid_size;
    state.hoverX = relX * N;
    state.hoverY = relY * N;

    const cx = (N - 1) / 2.0;
    const cy = (N - 1) / 2.0;
    const scaleArcsec = state.skyCube.query.pixel_scale_arcsec;
    const cosDec = Math.max(Math.cos((state.dec * Math.PI) / 180.0), 0.05);

    const curRa = state.ra - ((state.hoverX - cx) * (scaleArcsec / 3600.0)) / cosDec;
    const curDec = state.dec - (state.hoverY - cy) * (scaleArcsec / 3600.0);

    const ix = Math.max(0, Math.min(N - 1, Math.floor(state.hoverX)));
    const iy = Math.max(0, Math.min(N - 1, Math.floor(state.hoverY)));
    const epIdx = Math.max(0, Math.min(state.activeEpoch, state.skyCube.epochs.length - 1));
    const mat =
      state.activeBand === "ALL"
        ? state.skyCube.coadd_by_epoch[epIdx]
        : state.skyCube.frames_by_epoch[epIdx][state.activeBand] || state.skyCube.coadd_by_epoch[epIdx];
    const fluxVal = mat ? mat[iy][ix] : 0.0;
    const flagMask = state.skyCube.flags_by_epoch[epIdx][iy][ix];
    const activeFlags = [];
    for (const f of state.flagDefinitions) {
      if (flagMask & f.mask) activeFlags.push(f.name);
    }
    const flagStr = activeFlags.length ? activeFlags.join("|") : "CLEAN";

    document.getElementById("cursorCoordReadout").textContent =
      `RA ${curRa.toFixed(4)}° · Dec ${curDec >= 0 ? "+" : ""}${curDec.toFixed(4)}° · ${fluxVal.toFixed(2)} MJy/sr · FLAGS: ${flagStr}`;

    portalLoupe.classList.add("visible");
    renderOpticalLoupe(skyCanvas);
  });

  skyCanvas.addEventListener("mouseleave", () => {
    state.hoverX = null;
    state.hoverY = null;
    portalLoupe.classList.remove("visible");
  });

  skyCanvas.addEventListener("click", (e) => {
    if (!state.skyCube || state.isDraggingSplit) return;
    const rect = skyCanvas.getBoundingClientRect();
    const relX = clamp01((e.clientX - rect.left) / rect.width);
    const relY = clamp01((e.clientY - rect.top) / rect.height);
    if (state.viewMode === "SPLIT_CURTAIN") state.splitRatio = relX;

    const N = state.skyCube.query.grid_size;
    state.apertureX = relX * N;
    state.apertureY = relY * N;
    renderSkyCanvas();
    loadPointSpectrum(state.apertureX, state.apertureY);
  });

  // 102-Channel Spectrum Hover
  const specCanvas = document.getElementById("spectrumCanvas");
  specCanvas.addEventListener("mousemove", (e) => {
    if (!state.spectrumData || !state.spectrumData.spectrum_102) return;
    const rect = specCanvas.getBoundingClientRect();
    const relX = clamp01((e.clientX - rect.left) / rect.width);
    const targetWl = 0.72 + relX * (5.05 - 0.72);
    const wl = state.spectrumData.spectrum_102.wavelengths_um;
    let bestIdx = 0;
    let bestDist = 999;
    for (let i = 0; i < wl.length; i++) {
      const d = Math.abs(wl[i] - targetWl);
      if (d < bestDist) {
        bestDist = d;
        bestIdx = i;
      }
    }
    state.hoveredSpecChannel = bestIdx;
    const s = state.spectrumData.spectrum_102;
    document.getElementById("spectrumHoverReadout").innerHTML = `
      <span>Ch ${bestIdx + 1}/102 · <strong style="color:var(--text-primary);">${s.band_ids[bestIdx]}</strong> · ${s.wavelengths_um[bestIdx].toFixed(3)} μm</span>
      <span style="color:var(--accent-cyan);">${s.flux_mjy[bestIdx].toFixed(2)} mJy (R=${s.resolving_power[bestIdx]})</span>
    `;
    renderSpectrumCanvas();
  });
  specCanvas.addEventListener("mouseleave", () => {
    state.hoveredSpecChannel = null;
    renderSpectrumCanvas();
  });

  // NASA Demo Presets Dropdown & Universal Sky Search Bar
  const presetsDropdown = document.getElementById("presetsDropdown");
  const togglePresetsBtn = document.getElementById("togglePresetsBtn");

  togglePresetsBtn.addEventListener("click", (e) => {
    e.stopPropagation();
    presetsDropdown.classList.toggle("hidden");
  });

  window.addEventListener("click", (e) => {
    if (!presetsDropdown.contains(e.target) && e.target !== togglePresetsBtn) {
      presetsDropdown.classList.add("hidden");
    }
  });

  window.addEventListener("keydown", (e) => {
    if (e.key === "Escape") {
      presetsDropdown.classList.add("hidden");
      document.getElementById("shareModal").classList.add("hidden");
      document.getElementById("judgeModal").classList.add("hidden");
    }
  });

  document.getElementById("fovSelect").addEventListener("change", (e) => {
    state.fov = parseFloat(e.target.value) || 4.5;
    loadSkyCube();
  });

  document.getElementById("randomSectorBtn").addEventListener("click", () => {
    const portal = document.getElementById("telescopePortal");
    portal.classList.add("slewing");
    // Pick an unexplored Ecliptic / Outer Solar System sky sector
    const randRa = +(Math.random() * 360.0).toFixed(5);
    const randDec = +((Math.random() - 0.5) * 110.0).toFixed(5);
    state.ra = randRa;
    state.dec = randDec;
    document.getElementById("raInput").value = randRa.toFixed(5);
    document.getElementById("decInput").value = (randDec >= 0 ? "+" : "") + randDec.toFixed(5);
    loadSkyCube().finally(() => {
      setTimeout(() => portal.classList.remove("slewing"), 220);
    });
  });

  document.getElementById("skySearchForm").addEventListener("submit", async (e) => {
    e.preventDefault();
    const rawRa = document.getElementById("raInput").value.trim();
    const rawDec = document.getElementById("decInput").value.trim();
    state.fov = parseFloat(document.getElementById("fovSelect").value) || 4.5;

    // If user typed a target name in RA field (e.g. "WISE 0855"), resolve by name alone;
    // otherwise combine RA + Dec
    const isNameQuery = /[a-zA-Z]{3,}/.test(rawRa) && !/^\d+h/i.test(rawRa);
    const queryStr = isNameQuery || !rawDec ? rawRa : `${rawRa} ${rawDec}`;

    const portal = document.getElementById("telescopePortal");
    portal.classList.add("slewing");
    try {
      const resp = await fetch(`/api/resolve?q=${encodeURIComponent(queryStr)}`);
      if (!resp.ok) return;
      const res = await resp.json();
      state.ra = res.ra;
      state.dec = res.dec;
      document.getElementById("raInput").value = res.ra.toFixed(6);
      document.getElementById("decInput").value = (res.dec >= 0 ? "+" : "") + res.dec.toFixed(6);
      presetsDropdown.classList.add("hidden");
      await loadSkyCube();
    } catch (err) {
      console.error(err);
    } finally {
      setTimeout(() => portal.classList.remove("slewing"), 220);
    }
  });

  document.getElementById("shareUrlBtn").addEventListener("click", () => {
    syncUrlParameters();
    const url = window.location.href;
    document.getElementById("shareDirectUrlInput").value = url;
    document.getElementById("shareIframeInput").value =
      `<iframe src="${url}" width="100%" height="780" frameborder="0" allowfullscreen></iframe>`;
    document.getElementById("shareModal").classList.remove("hidden");
  });
  document.getElementById("closeShareModalBtn").addEventListener("click", () => {
    document.getElementById("shareModal").classList.add("hidden");
  });
  document.getElementById("copyDirectUrlBtn").addEventListener("click", () => {
    const inp = document.getElementById("shareDirectUrlInput");
    inp.select();
    navigator.clipboard.writeText(inp.value);
    document.getElementById("copyDirectUrlBtn").textContent = "Copied";
    setTimeout(() => (document.getElementById("copyDirectUrlBtn").textContent = "Copy"), 1400);
  });

  document.getElementById("judgeBriefBtn").addEventListener("click", () => {
    document.getElementById("judgeModal").classList.remove("hidden");
  });
  document.getElementById("closeJudgeModalBtn").addEventListener("click", () => {
    document.getElementById("judgeModal").classList.add("hidden");
  });

  document.getElementById("exportFitsBtn").addEventListener("click", () => {
    window.open(`/api/export/fits?ra=${state.ra}&dec=${state.dec}&fov=${state.fov}&band=${state.activeBand}`, "_blank");
  });

  document.getElementById("downloadSpectrumCsvBtn").addEventListener("click", () => {
    if (!state.spectrumData || !state.spectrumData.spectrum_102) return;
    const s = state.spectrumData.spectrum_102;
    let csv = "channel,detector_band,wavelength_um,flux_mjy,flux_err_mjy,resolving_power\n";
    for (let i = 0; i < s.wavelengths_um.length; i++) {
      csv += `${i + 1},${s.band_ids[i]},${s.wavelengths_um[i]},${s.flux_mjy[i]},${s.flux_err_mjy[i]},${s.resolving_power[i]}\n`;
    }
    const blob = new Blob([csv], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `spherex_102spectrum_ra${state.ra.toFixed(3)}_dec${state.dec.toFixed(3)}.csv`;
    a.click();
  });
}

async function initApp() {
  readUrlParameters();
  bindEvents();
  try {
    const resp = await fetch("/api/targets");
    const meta = await resp.json();
    state.targets = meta.targets || [];
    state.bands = meta.bands || [];
    state.spectralFeatures = meta.spectral_features || [];
    state.flagDefinitions = meta.flag_definitions || [];

    updateAmbientAura();
    updateViewModeTabs();
    await loadSkyCube();
    restartBlinkTimer();
  } catch (err) {
    console.error("Init error:", err);
  }
}

window.addEventListener("DOMContentLoaded", initApp);
