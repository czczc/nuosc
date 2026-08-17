// L/E view: oscillation probabilities vs log10(L/E), NuFast engine.
// In vacuum every P is an exact function of L/E, so each channel collapses onto a
// single curve (x = log L/E, y = P) — both dm^2 frequencies visible at once, with
// the real experiments pinned at their L/E. The third axis carries the energy E
// across the active experiment's window: the surface is the header channel
// computed IN MATTER at the shared rho, drawn oscillogram-style (palette colors)
// and clipped to the experiment's reach (L up to twice the baseline), so it is
// extruded flat along E at rho = 0 (the L/E degeneracy) and twists otherwise.
import * as THREE from 'three';
import { SceneBase, PALETTES, textSprite } from '../three/SceneBase.js';
import { probabilityMatter } from '../engines/nufast.js';
import {
  engineParams, PRESETS, eRangeOf, lRangeOf, eUnitOf, eStepOf, lStepOf, fmtE, CHANNELS, channelLabel, pLabel,
} from '../engines/constants.js';
import { theme } from '../theme.js';
import { plot2d, legend } from './plot2d.js';

const SX = 10, SY = 2.5, DE = 6;
const LX0 = Math.log10(30), LX1 = 5;    // log10(L/E [km/GeV]) axis span
const NCURVE = 8192;                    // vacuum-curve samples (fast wiggles at large L/E)
const NX = 128, NZ = 32;                // matter surface grid
const LOOP_S = 10;                      // ~10 s per marker sweep, as in the other views
const ZF = DE / 2 + 0.5;                // front plane: vacuum curves + markers
const CH_COLOR = { mue: '#e05545', mumu: '#4488ee', mutau: '#44aa55', ee: '#e0a030' };
const CH_KEY = { mue: 'cMue', mumu: 'cMumu', mutau: 'cMutau', ee: 'cEe' };
const palette = PALETTES.rainbow; // surface colors, as on the oscillogram page

const xw = (lx) => -SX / 2 + SX * (lx - LX0) / (LX1 - LX0);
// the experiment's live L/E point (from the shared sliders) as a 0..1 axis fraction
function fracOfStore(store) {
  const lx = Math.log10(Math.max(1e-9, store.L / store.E));
  return Math.max(0, Math.min(1, (lx - LX0) / (LX1 - LX0)));
}
const zOfE = (E, E0, E1) => DE / 2 - ((E - E0) / (E1 - E0)) * DE;
const fmtX = (x) => (x >= 9999.5 ? `${(x / 1000).toFixed(1)}k` : `${Math.round(x)}`);

// vacuum P as a function of x = L/E alone: NuFast with rho = 0, L = x, E = ±1
function pVac(ep, ch, x, anti) {
  const { a, b } = CHANNELS[ch];
  const P = probabilityMatter(ep.s12sq, ep.s13sq, ep.s23sq, ep.delta, ep.dm21, ep.dm31, x, anti ? -1 : 1, 0, 1);
  return Math.max(0, Math.min(1, P[a][b]));
}

// <sin^2 D31> averaging weight: 0 below an atmospheric phase of 15 rad, 1 above
// 25, smoothstepped in between — a hard threshold leaves a visible step in the
// curve where the averaging switches on
function avgWeight(ep, x) {
  const t = (1.267 * Math.abs(ep.dm31) * x - 15) / 10;
  if (t <= 0) return 0;
  if (t >= 1) return 1;
  return t * t * (3 - 2 * t);
}
const atmPeriod = (ep) => Math.PI / (1.267 * Math.abs(ep.dm31)); // one sin^2 period in x

// optional <sin^2 D31> averaging for the vacuum curves (detectors can't resolve
// the fast wiggles at large L/E): mean over one atmospheric period in x,
// cross-faded in over the transition window
function pVacAvg(ep, ch, x, anti, avg) {
  const w = avg ? avgWeight(ep, x) : 0;
  if (w === 0) return pVac(ep, ch, x, anti);
  const period = atmPeriod(ep);
  let s = 0;
  const K = 9;
  for (let k = 0; k < K; k++) s += pVac(ep, ch, x + period * ((k + 0.5) / K - 0.5), anti);
  const pa = s / K;
  return w === 1 ? pa : (1 - w) * pVac(ep, ch, x, anti) + w * pa;
}

// matter P at a real (L, E) point, shared rho
function pMat(ep, ch, L, E, rho, anti) {
  const { a, b } = CHANNELS[ch];
  const P = probabilityMatter(ep.s12sq, ep.s13sq, ep.s23sq, ep.delta, ep.dm21, ep.dm31, L, anti ? -E : E, rho * ep.Ye, 1);
  return Math.max(0, Math.min(1, P[a][b]));
}

// matter P vs x = L/E at fixed E, ALWAYS averaged over one atmospheric period at
// large L/E (cross-faded like the curves) — the surface shows what a detector
// with finite resolution would see. K = 5 (vs 9 for the vacuum curves): the
// surface rebuilds live during the ρ animation, and 5 evenly spaced samples
// already average a sinusoid exactly.
function pMatAvg(ep, ch, x, E, rho, anti) {
  const w = avgWeight(ep, x);
  if (w === 0) return pMat(ep, ch, x * E, E, rho, anti);
  const period = atmPeriod(ep);
  let s = 0;
  const K = 5;
  for (let k = 0; k < K; k++) s += pMat(ep, ch, (x + period * ((k + 0.5) / K - 0.5)) * E, E, rho, anti);
  const pa = s / K;
  return w === 1 ? pa : (1 - w) * pMat(ep, ch, x * E, E, rho, anti) + w * pa;
}

// companion memos: curves keyed by the physics inputs (marker excluded)
let compCache = null;
let comp2Cache = null;

// built-in experiments at their L/E (flux peak); reactor markers in the ee color
const MARKS = Object.entries(PRESETS)
  .map(([name, p]) => ({ name, lx: Math.log10(p.L / p.Epeak), reactor: !!p.anti }))
  .sort((m1, m2) => m1.lx - m2.lx);

export default {
  id: 'loe',
  label: 'Worldline (L/E)',
  extras: [
    { key: 'cMue', type: 'checkbox', label: 'νμ→νe', lock: (s) => s.channel === 'mue' },
    { key: 'cMumu', type: 'checkbox', label: 'νμ→νμ', lock: (s) => s.channel === 'mumu' },
    { key: 'cMutau', type: 'checkbox', label: 'νμ→ντ', lock: (s) => s.channel === 'mutau' },
    { key: 'cEe', type: 'checkbox', label: 'νe→νe', lock: (s) => s.channel === 'ee' },
    { key: 'surf', type: 'checkbox', label: 'matter surface (exp. range)' },
    { key: 'exps', type: 'checkbox', label: 'experiment markers' },
    {
      key: 'avg', type: 'checkbox', label: '⟨sin²Δ₃₁⟩ averaged',
      title: 'Detector resolution can\'t follow the fast Δ₃₁ wiggles at large L/E, so they are '
        + 'cross-faded into their average as the phase grows past ~15–25 rad.',
    },
    {
      // the marker drives a shared slider; the white L/E cross-section always
      // tracks the experiment's live L/E point (E and L move it, ρ morphs the
      // surface through the vacuum collapse at ρ = 0, δCP reshapes the curves)
      key: 'marker', type: 'marker', label: 'animate', step: 0.002,
      select: {
        key: 'anim',
        options: [
          { value: 'L', label: 'L' },
          { value: 'E', label: 'E' },
          { value: 'dcp', label: 'δCP' },
          { value: 'rho', label: 'ρ' },
        ],
      },
    },
  ],

  create(container, store) {
    const base = new SceneBase(container, { camPos: [-7, 4.5, 9], target: [0, 1, 0] });

    const grid = new THREE.GridHelper(SX, 10, theme().grid1, theme().grid2);
    grid.position.y = -0.01;
    base.scene.add(grid);

    // decade marks on the floor + front-edge labels
    for (const d of [2, 3, 4]) {
      const x = xw(d);
      base.scene.add(new THREE.Line(
        new THREE.BufferGeometry().setFromPoints(
          [new THREE.Vector3(x, 0, ZF), new THREE.Vector3(x, 0, -DE / 2)]),
        new THREE.LineBasicMaterial({ color: theme().axis, transparent: true, opacity: 0.5 })));
      const s = textSprite(`10${{ 2: '²', 3: '³', 4: '⁴' }[d]}`, 0.8);
      s.position.set(x, 0.05, ZF + 0.45);
      base.scene.add(s);
    }
    const xLabel = textSprite('L/E [km/GeV]  (30 – 10⁵, log)');
    xLabel.position.set(0, -0.45, ZF + 1.3);
    base.scene.add(xLabel);
    const yLabel = textSprite('P');
    yLabel.position.set(-SX / 2 - 0.6, SY, ZF);
    base.scene.add(yLabel);
    let eLabel = null, lastEKey = null;

    // vacuum curves: thin line per channel; the header channel drawn as a tube
    const lines = {};
    for (const ch of Object.keys(CHANNELS)) {
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(NCURVE * 3), 3));
      const line = new THREE.Line(geo, new THREE.LineBasicMaterial({
        color: new THREE.Color(CH_COLOR[ch]), transparent: true, opacity: 0.85,
      }));
      base.scene.add(line);
      lines[ch] = line;
    }
    const focusTube = new THREE.Mesh(new THREE.BufferGeometry(), new THREE.MeshBasicMaterial());
    base.scene.add(focusTube);

    // matter surface (indexed grid), oscillogram-style: vertex colors from the palette
    const sheetGeo = new THREE.BufferGeometry();
    sheetGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(NX * NZ * 3), 3));
    sheetGeo.setAttribute('color', new THREE.BufferAttribute(new Float32Array(NX * NZ * 3), 3));
    {
      const idx = [];
      for (let j = 0; j < NZ - 1; j++) for (let i = 0; i < NX - 1; i++) {
        const a = j * NX + i, b = a + 1, c = a + NX, d = c + 1;
        idx.push(a, b, c, b, d, c);
      }
      sheetGeo.setIndex(idx);
    }
    const sheetMesh = new THREE.Mesh(sheetGeo, new THREE.MeshLambertMaterial({
      vertexColors: true, side: THREE.DoubleSide,
    }));
    base.scene.add(sheetMesh);

    // experiment markers: static dashed verticals + labels at their true L/E
    const markGroup = new THREE.Group();
    base.scene.add(markGroup);
    MARKS.forEach((m, i) => {
      const x = xw(m.lx);
      const col = m.reactor ? CH_COLOR.ee : theme().axis;
      const lineGeo = new THREE.BufferGeometry().setFromPoints(
        [new THREE.Vector3(x, 0, ZF), new THREE.Vector3(x, SY, ZF)]);
      const line = new THREE.LineSegments(lineGeo, new THREE.LineDashedMaterial({
        color: new THREE.Color(col), dashSize: 0.12, gapSize: 0.1, transparent: true, opacity: 0.7,
      }));
      line.computeLineDistances();
      markGroup.add(line);
      const s = textSprite(m.name, 0.8, m.reactor ? CH_COLOR.ee : undefined);
      s.position.set(x, SY + 0.28 + (i % 4) * 0.3, ZF);
      markGroup.add(s);
    });

    // marker slice: vertical line + dot on the focus curve, plus the fixed-L/E
    // cross-section across the sheet (flat in vacuum, wiggly in matter)
    const dropGeo = new THREE.BufferGeometry();
    dropGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(2 * 3), 3));
    const dropLine = new THREE.Line(dropGeo, new THREE.LineBasicMaterial({ color: theme().hi }));
    base.scene.add(dropLine);
    const dot = new THREE.Mesh(new THREE.SphereGeometry(0.07, 12, 12), new THREE.MeshBasicMaterial({ color: theme().hi }));
    base.scene.add(dot);
    const sliceGeo = new THREE.BufferGeometry();
    sliceGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(NZ * 3), 3));
    const sliceLine = new THREE.Line(sliceGeo, new THREE.LineBasicMaterial({ color: theme().hi }));
    base.scene.add(sliceLine);

    let cache = null; // { sheetP, sheetLxMax, E0, E1, focusYs, focus, showSurf }
    let lastCurveKey = null;
    const focusYs = new Float32Array(NCURVE); // persists across rho-only rebuilds

    function update() {
      const ep = engineParams(store);
      const anti = store.anti;
      const focus = store.channel;
      const rho = store.rho;
      const vs = store.views.loe; // checkboxes only — never play/marker (tick-only)
      const [E0, E1] = eRangeOf(store.basePreset);

      // vacuum curves don't depend on rho — skip them while the ρ animation only
      // morphs the sheet (same memo pattern as the oscillogram surface)
      const curveKey = JSON.stringify([ep, anti, focus, vs.avg, vs.cMue, vs.cMumu, vs.cMutau, vs.cEe]);
      if (curveKey !== lastCurveKey) {
        lastCurveKey = curveKey;
        // the focus channel becomes a tube (GL lines are 1px)
        const focusPts = [];
        for (const ch of Object.keys(CHANNELS)) {
          const isFocus = ch === focus;
          const show = isFocus || vs[CH_KEY[ch]];
          lines[ch].visible = show && !isFocus;
          if (!show) continue;
          const pos = lines[ch].geometry.attributes.position;
          for (let i = 0; i < NCURVE; i++) {
            const lx = LX0 + (i / (NCURVE - 1)) * (LX1 - LX0);
            const p = pVacAvg(ep, ch, 10 ** lx, anti, vs.avg);
            if (isFocus) {
              focusYs[i] = p;
              if (i % 4 === 0) focusPts.push(new THREE.Vector3(xw(lx), p * SY, ZF));
            } else {
              pos.setXYZ(i, xw(lx), p * SY, ZF);
            }
          }
          if (!isFocus) pos.needsUpdate = true;
        }
        focusTube.geometry.dispose();
        focusTube.geometry = new THREE.TubeGeometry(new THREE.CatmullRomCurve3(focusPts), 1024, 0.025, 6, false);
        focusTube.material.color.set(CH_COLOR[focus]);
      }

      // matter surface for the focus channel over (log L/E, E) at the shared rho,
      // clipped to the experiment's reach: row j (energy E) spans L/E only up to
      // L = Lmax (twice the baseline, the app's standard L span) — a slanted band
      const sheetP = new Float32Array(NX * NZ);
      const sheetLxMax = new Float32Array(NZ);
      sheetMesh.visible = vs.surf;
      if (vs.surf) {
        const Lmax = lRangeOf(store.basePreset)[1];
        const pos = sheetGeo.attributes.position;
        const col = sheetGeo.attributes.color;
        let maxP = 1e-6;
        for (let j = 0; j < NZ; j++) {
          const E = E0 + (j / (NZ - 1)) * (E1 - E0);
          const z = zOfE(E, E0, E1);
          const lxMax = Math.max(LX0 + 0.05, Math.min(LX1, Math.log10(Lmax / E)));
          sheetLxMax[j] = lxMax;
          for (let i = 0; i < NX; i++) {
            const lx = LX0 + (i / (NX - 1)) * (lxMax - LX0);
            const p = pMatAvg(ep, focus, 10 ** lx, E, rho, anti);
            sheetP[j * NX + i] = p;
            if (p > maxP) maxP = p;
            pos.setXYZ(j * NX + i, xw(lx), p * SY, z);
          }
        }
        for (let k = 0; k < NX * NZ; k++) {
          const [r, g, b] = palette(sheetP[k] / maxP);
          col.setXYZ(k, r, g, b);
        }
        pos.needsUpdate = true;
        col.needsUpdate = true;
        sheetGeo.computeVertexNormals();
        sheetGeo.computeBoundingSphere();
      }

      markGroup.visible = vs.exps;

      const { unit, scale } = eUnitOf(store.basePreset);
      const eKey = `E [${unit}]  (${eRangeOf(store.basePreset).map((v) => +(v * scale).toFixed(4)).join(' - ')})`;
      if (eKey !== lastEKey) {
        if (eLabel) base.scene.remove(eLabel);
        eLabel = textSprite(eKey);
        eLabel.position.set(-SX / 2 - 1.6, -0.45, 0);
        base.scene.add(eLabel);
        lastEKey = eKey;
      }

      cache = { sheetP, sheetLxMax, E0, E1, focusYs, focus, showSurf: vs.surf };
    }

    let lastMarker = store.views.loe.marker;

    function tick(dt) {
      const vs = store.views.loe;
      if (vs.play) vs.marker = (vs.marker + dt / LOOP_S) % 1;
      // the marker drives the shared slider of the chosen variable (rounded to
      // its step, like the other views); manual slider drags stay untouched
      if (vs.marker !== lastMarker) {
        if (vs.anim === 'rho') {
          store.rho = Math.round((vs.marker * 5) / 0.05) * 0.05;
        } else if (vs.anim === 'dcp') {
          store.dcp = Math.round(vs.marker * 360);
        } else if (vs.anim === 'L') {
          const Lmax = lRangeOf(store.basePreset)[1];
          const s = lStepOf(store.basePreset);
          store.L = Math.round((vs.marker * Lmax) / s) * s;
        } else {
          const [E0, E1] = eRangeOf(store.basePreset);
          const s = eStepOf(store.basePreset);
          store.E = Math.round((E0 + vs.marker * (E1 - E0)) / s) * s;
        }
      }
      lastMarker = vs.marker;
      if (!cache) return;
      const frac = fracOfStore(store); // cross-section at the live L/E point
      const fi = frac * (NCURVE - 1), i0 = Math.floor(fi), t = fi - i0;
      const pv = cache.focusYs[i0] * (1 - t) + cache.focusYs[Math.min(NCURVE - 1, i0 + 1)] * t;
      const x = xw(LX0 + frac * (LX1 - LX0));
      const dp = dropGeo.attributes.position;
      dp.setXYZ(0, x, 0, ZF);
      dp.setXYZ(1, x, SY, ZF);
      dp.needsUpdate = true;
      dropGeo.computeBoundingSphere();
      dot.position.set(x, pv * SY, ZF);
      // cross-section across the surface: rows are individually parameterized
      // (slanted band), so include only rows the marker's L/E actually reaches
      let count = 0;
      if (cache.showSurf) {
        const lxm = LX0 + frac * (LX1 - LX0);
        const sp = sliceGeo.attributes.position;
        for (let j = 0; j < NZ; j++) {
          const lxMax = cache.sheetLxMax[j];
          if (lxm > lxMax) continue;
          const fi = ((lxm - LX0) / (lxMax - LX0)) * (NX - 1);
          const g0 = Math.floor(fi), gt = fi - g0, g1 = Math.min(NX - 1, g0 + 1);
          const p = cache.sheetP[j * NX + g0] * (1 - gt) + cache.sheetP[j * NX + g1] * gt;
          sp.setXYZ(count++, x, p * SY + 0.01, zOfE(cache.E0 + (j / (NZ - 1)) * (cache.E1 - cache.E0), cache.E0, cache.E1));
        }
        sp.needsUpdate = true;
        sliceGeo.setDrawRange(0, count);
        sliceGeo.computeBoundingSphere();
      }
      sliceLine.visible = cache.showSurf && count >= 2;
    }

    function probe(event) {
      if (!cache) return null;
      const ep = engineParams(store);
      const hit = cache.showSurf ? base.raycast(event, sheetMesh) : null;
      if (hit) {
        const lx = LX0 + (hit.point.x / SX + 0.5) * (LX1 - LX0);
        const E = cache.E0 + (0.5 - hit.point.z / DE) * (cache.E1 - cache.E0);
        const x = 10 ** lx;
        const pm = pMatAvg(ep, cache.focus, x, E, store.rho, store.anti);
        const pv = pVac(ep, cache.focus, x, store.anti);
        return `L/E ${fmtX(x)} km/GeV · E ${fmtE(E, store.basePreset)} · L ${fmtX(x * E)} km` +
          ` · P ${pm.toFixed(4)} · vac ${pv.toFixed(4)}`;
      }
      // no hover -> live status at the experiment's current L/E point
      const frac = fracOfStore(store);
      const x = 10 ** (LX0 + frac * (LX1 - LX0));
      let out = `L/E ${fmtX(x)} km/GeV · ρ ${store.rho.toFixed(2)} · ${pLabel(cache.focus, store.anti)} vac ${pVac(ep, cache.focus, x, store.anti).toFixed(4)}`;
      if (cache.showSurf) {
        const lxm = LX0 + frac * (LX1 - LX0);
        let lo = 1, hi = 0, any = false;
        for (let j = 0; j < NZ; j++) {
          const lxMax = cache.sheetLxMax[j];
          if (lxm > lxMax) continue;
          const g0 = Math.round(((lxm - LX0) / (lxMax - LX0)) * (NX - 1));
          const p = cache.sheetP[j * NX + g0];
          if (p < lo) lo = p;
          if (p > hi) hi = p;
          any = true;
        }
        if (any) out += ` · matter over E: ${lo.toFixed(4)} – ${hi.toFixed(4)}`;
      }
      return out;
    }

    return { base, update, tick, probe, dispose: () => base.dispose() };
  },

  // front view: the classic 2D L/E plot (vacuum curves). Unresolved fast wiggles
  // are drawn as their per-column min–max envelope (a plain polyline aliases into
  // jagged moiré); curves are memoized so per-frame marker repaints stay cheap.
  companion: {
    title: () => 'P vs L/E (vacuum)',
    markerDriven: true,
    draw(canvas, store) {
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.clientWidth * dpr, h = canvas.clientHeight * dpr;
      canvas.width = w; canvas.height = h;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = theme().canvas; ctx.fillRect(0, 0, w, h);

      const ep = engineParams(store);
      const vs = store.views.loe;
      const focus = store.channel;
      const shown = Object.keys(CHANNELS).filter((ch) => ch === focus || vs[CH_KEY[ch]]);
      const NCOL = 380, KSUB = 5; // envelope columns x subsamples
      const key = JSON.stringify([ep, store.anti, vs.avg, shown]);
      if (compCache?.key !== key) {
        const env = {};
        for (const ch of shown) {
          const lo = new Float32Array(NCOL), hi = new Float32Array(NCOL), mid = new Float32Array(NCOL);
          for (let i = 0; i < NCOL; i++) {
            let a = 1, b = 0, s = 0;
            for (let k = 0; k < KSUB; k++) {
              const lx = LX0 + (LX1 - LX0) * (i + (k + 0.5) / KSUB) / NCOL;
              const p = pVacAvg(ep, ch, 10 ** lx, store.anti, vs.avg);
              if (p < a) a = p;
              if (p > b) b = p;
              s += p;
            }
            lo[i] = a; hi[i] = b; mid[i] = s / KSUB;
          }
          env[ch] = { lo, hi, mid };
        }
        compCache = { key, env };
      }

      const P = plot2d(ctx, w, h, dpr, {
        x: [LX0, LX1], y: [0, 1], xTitle: 'log₁₀(L/E [km/GeV])', yTitle: 'P', xTicks: 4,
      });
      const XC = (i) => P.X(LX0 + (LX1 - LX0) * (i + 0.5) / NCOL);
      const items = [];
      for (const ch of shown) {
        const { lo, hi, mid } = compCache.env[ch];
        items.push({ label: channelLabel(ch), color: CH_COLOR[ch] });
        ctx.fillStyle = CH_COLOR[ch];
        ctx.globalAlpha = 0.6;
        ctx.beginPath();
        for (let i = 0; i < NCOL; i++) ctx.lineTo(XC(i), P.Y(hi[i]));
        for (let i = NCOL - 1; i >= 0; i--) ctx.lineTo(XC(i), P.Y(lo[i]));
        ctx.closePath();
        ctx.fill();
        ctx.globalAlpha = 1;
        ctx.strokeStyle = CH_COLOR[ch];
        ctx.lineWidth = (ch === focus ? 2.4 : 1.2) * dpr;
        ctx.beginPath();
        for (let i = 0; i < NCOL; i++) ctx.lineTo(XC(i), P.Y(mid[i]));
        ctx.stroke();
      }
      const mx = P.X(LX0 + (LX1 - LX0) * fracOfStore(store));
      ctx.strokeStyle = theme().hiCss;
      ctx.lineWidth = dpr;
      ctx.beginPath(); ctx.moveTo(mx, P.mt); ctx.lineTo(mx, P.mt + P.ph); ctx.stroke();
      legend(ctx, dpr, P, items);
    },
  },

  // the spectrum the experiment actually measures: P vs E for the header channel
  // at the shared L and ρ, matter effect included, with the vacuum curve dashed
  // for reference — the gap between the two is the matter effect at each energy
  companion2: {
    // "ρ" is avoided here: the panel caption is uppercased by CSS, turning it into
    // a capital rho that reads like a Latin P
    title: (store) => `P vs E at L = ${store.L} km, ${store.rho.toFixed(2)} g/cm³`,
    markerDriven: true,
    draw(canvas, store) {
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.clientWidth * dpr, h = canvas.clientHeight * dpr;
      canvas.width = w; canvas.height = h;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = theme().canvas; ctx.fillRect(0, 0, w, h);

      const ep = engineParams(store);
      const focus = store.channel;
      const [E0, E1] = eRangeOf(store.basePreset);
      const { unit, scale } = eUnitOf(store.basePreset);
      const NPT = 300;
      const key = JSON.stringify([ep, store.anti, focus, store.L, store.rho, E0, E1]);
      if (comp2Cache?.key !== key) {
        const mat = new Float32Array(NPT), vac = new Float32Array(NPT);
        let pmax = 1e-6;
        for (let i = 0; i < NPT; i++) {
          const E = E0 + (E1 - E0) * i / (NPT - 1);
          mat[i] = pMat(ep, focus, store.L, E, store.rho, store.anti);
          vac[i] = pVac(ep, focus, store.L / E, store.anti);
          pmax = Math.max(pmax, mat[i], vac[i]);
        }
        comp2Cache = { key, mat, vac, pmax };
      }

      // autoscaled y: small-probability channels (numu->nue) fill the panel
      const P = plot2d(ctx, w, h, dpr, {
        x: [E0 * scale, E1 * scale], y: [0, comp2Cache.pmax * 1.08],
        xTitle: `E [${unit}]`, yTitle: pLabel(focus, store.anti),
      });
      const XE = (i) => P.X((E0 + (E1 - E0) * i / (NPT - 1)) * scale);
      // vacuum reference (dashed)
      ctx.strokeStyle = CH_COLOR[focus];
      ctx.globalAlpha = 0.55;
      ctx.lineWidth = 1.2 * dpr;
      ctx.setLineDash([5 * dpr, 4 * dpr]);
      ctx.beginPath();
      for (let i = 0; i < NPT; i++) ctx.lineTo(XE(i), P.Y(comp2Cache.vac[i]));
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.globalAlpha = 1;
      // matter (solid)
      ctx.lineWidth = 1.8 * dpr;
      ctx.beginPath();
      for (let i = 0; i < NPT; i++) ctx.lineTo(XE(i), P.Y(comp2Cache.mat[i]));
      ctx.stroke();
      legend(ctx, dpr, P, [
        { label: 'matter', color: CH_COLOR[focus] },
        { label: 'vacuum (dashed)', color: theme().ink },
      ]);
      // white line at the shared E (the spectrum position the L/E marker sits on)
      const mx = P.X(Math.max(E0, Math.min(E1, store.E)) * scale);
      ctx.strokeStyle = theme().hiCss;
      ctx.lineWidth = dpr;
      ctx.beginPath(); ctx.moveTo(mx, P.mt); ctx.lineTo(mx, P.mt + P.ph); ctx.stroke();
    },
  },
};
