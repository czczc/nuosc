// Oscillogram surface: P(numu->nue) over (E x second-axis), NuFast engine.
// Reference implementation of the view contract (see views/index.js).
import * as THREE from 'three';
import { SceneBase, PALETTES, textSprite } from '../three/SceneBase.js';
import { probabilityMatter } from '../engines/nufast.js';
import { engineParams, PRESETS, eRangeOf, lRangeOf } from '../engines/constants.js';
import { theme } from '../theme.js';
import { plot2d, legend } from './plot2d.js';

const SX = 8, SZ = 8, SY = 2.2, N = 161;
const RANGES = { dcp: [0, 360], rho: [0, 5] };
const AXIS_LABEL = { L: 'L [km]', dcp: 'dCP [deg]', rho: 'rho [g/cm3]' };
const LOOP_S = 10; // ~10 s per sweep of the marker, as in phasors

function pMuE(ep, E, L, rho, dcp, anti) {
  const P = probabilityMatter(ep.s12sq, ep.s13sq, ep.s23sq, dcp * Math.PI / 180, ep.dm21, ep.dm31, L, anti ? -E : E, rho * ep.Ye, 1);
  return Math.max(0, Math.min(1, P[1][0]));
}


export default {
  id: 'oscillogram',
  label: 'Oscillogram',
  extras: [
    {
      key: 'axis2', type: 'select', label: '2nd axis',
      options: [
        { value: 'L', label: (store) => `L [0–${lRangeOf(store.basePreset)[1]} km]` },
        { value: 'dcp', label: 'δCP [0–360°]' },
        { value: 'rho', label: 'ρ [0–5 g/cm³]' },
      ],
    },
    {
      key: 'palette', type: 'select', label: 'palette',
      options: [
        { value: 'rainbow', label: 'rainbow' },
        { value: 'viridis', label: 'viridis' },
        { value: 'coolwarm', label: 'cool–warm' },
        { value: 'grayscale', label: 'grayscale' },
      ],
    },
    {
      key: 'marker', type: 'marker', label: 'animate', step: 0.002,
      select: {
        key: 'anim',
        options: [
          { value: 'L', label: 'L' },
          { value: 'E', label: 'E' },
          { value: 'dcp', label: 'δCP' },
        ],
      },
    },
  ],

  create(container, store) {
    const base = new SceneBase(container, { camPos: [-7, 5, 8] });

    const geo = new THREE.PlaneGeometry(SX, SZ, N - 1, N - 1);
    geo.rotateX(-Math.PI / 2);
    const colors = new Float32Array(geo.attributes.position.count * 3);
    geo.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    const surface = new THREE.Mesh(geo, new THREE.MeshLambertMaterial({ vertexColors: true, side: THREE.DoubleSide }));
    base.scene.add(surface);

    const sliceGeo = new THREE.BufferGeometry();
    sliceGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(N * 3), 3));
    const slice = new THREE.Line(sliceGeo, new THREE.LineBasicMaterial({ color: theme().hi }));
    base.scene.add(slice);

    // marker cross-section along z at the swept E (visible only when animating E)
    const eSliceGeo = new THREE.BufferGeometry();
    eSliceGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(N * 3), 3));
    const eSlice = new THREE.Line(eSliceGeo, new THREE.LineBasicMaterial({ color: theme().hi }));
    eSlice.visible = false;
    base.scene.add(eSlice);

    const grid = new THREE.GridHelper(Math.max(SX, SZ), 10, theme().grid1, theme().grid2);
    grid.position.y = -0.01;
    base.scene.add(grid);

    let xLabel = null, lastXKey = null;
    let zLabel = null, lastZKey = null;

    let maxP = 1, lastHeld = null, lastSurfKey = null;

    function update() {
      const ep = engineParams(store);
      const [E_MIN, E_MAX] = eRangeOf(store.basePreset);
      const axis2 = store.views.oscillogram.axis2;
      const [a2min, a2max] = axis2 === 'L' ? lRangeOf(store.basePreset) : RANGES[axis2];
      const held = { dcp: store.dcp, L: store.L, rho: store.rho };
      const anti = store.anti;

      // the surface doesn't depend on the swept variable's held value — skip the
      // rebuild when only that slider moved (e.g. play sweeping L with z axis = L)
      const paletteName = store.views.oscillogram.palette;
      const palette = PALETTES[paletteName] ?? PALETTES.rainbow;
      const surfKey = JSON.stringify([ep, E_MIN, E_MAX, axis2, a2min, a2max, anti, paletteName, { ...held, [axis2]: 0 }]);
      if (surfKey !== lastSurfKey) {
        lastSurfKey = surfKey;
        const pos = geo.attributes.position;
        const P = new Float32Array(pos.count);
        maxP = 1e-6;
        for (let i = 0; i < pos.count; i++) {
          const E = E_MIN + (pos.getX(i) / SX + 0.5) * (E_MAX - E_MIN);
          const v = a2min + (0.5 - pos.getZ(i) / SZ) * (a2max - a2min);
          held[axis2] = v;
          P[i] = pMuE(ep, E, held.L, held.rho, held.dcp, anti);
          if (P[i] > maxP) maxP = P[i];
        }
        for (let i = 0; i < pos.count; i++) {
          pos.setY(i, (P[i] / maxP) * SY);
          const [r, g, b] = palette(P[i] / maxP);
          colors[3 * i] = r; colors[3 * i + 1] = g; colors[3 * i + 2] = b;
        }
        pos.needsUpdate = true;
        geo.attributes.color.needsUpdate = true;
        geo.computeVertexNormals();
      }

      // slice: the legacy 2D curve at the axis-matching slider value
      const sv = { dcp: store.dcp, L: store.L, rho: store.rho }[axis2];
      const zSlice = (0.5 - (sv - a2min) / (a2max - a2min)) * SZ;
      const sp = sliceGeo.attributes.position;
      held[axis2] = sv;
      for (let i = 0; i < N; i++) {
        const E = E_MIN + (i / (N - 1)) * (E_MAX - E_MIN);
        const p = pMuE(ep, E, held.L, held.rho, held.dcp, anti);
        sp.setXYZ(i, (i / (N - 1) - 0.5) * SX, (p / maxP) * SY + 0.02, zSlice);
      }
      sp.needsUpdate = true;

      const xKey = `E [GeV]  (${E_MIN} - ${E_MAX})`;
      if (xKey !== lastXKey) {
        if (xLabel) base.scene.remove(xLabel);
        xLabel = textSprite(xKey);
        xLabel.position.set(0, 0.1, SZ / 2 + 0.8);
        base.scene.add(xLabel);
        lastXKey = xKey;
      }
      const zKey = `${AXIS_LABEL[axis2]}  (${a2min} - ${a2max})`;
      if (zKey !== lastZKey) {
        if (zLabel) base.scene.remove(zLabel);
        zLabel = textSprite(zKey);
        zLabel.position.set(-SX / 2 - 1.2, 0.1, 0);
        base.scene.add(zLabel);
        lastZKey = zKey;
      }
      lastHeld = { ...held, axis2, a2min, a2max, E_MIN, E_MAX, ep, anti };
    }

    // play/marker animation: the marker drives the shared slider of the chosen
    // variable. E sweeps the cross-section line (drawn at store.E), L or dCP moves
    // the white slice when it is the z axis and morphs the surface when held.
    let lastMarker = store.views.oscillogram.marker;
    let eApplied = null;

    function tick(dt) {
      const vs = store.views.oscillogram;
      if (vs.play) vs.marker = (vs.marker + dt / LOOP_S) % 1;
      if (vs.marker !== lastMarker) {
        const [v0, v1] = vs.anim === 'L' ? lRangeOf(store.basePreset)
          : vs.anim === 'E' ? eRangeOf(store.basePreset) : RANGES.dcp;
        const v = v0 + vs.marker * (v1 - v0);
        // rounded to the shared sliders' steps
        if (vs.anim === 'L') store.L = Math.round(v / 5) * 5;
        else if (vs.anim === 'E') store.E = Math.round(v * 100) / 100;
        else store.dcp = Math.round(v);
      }
      lastMarker = vs.marker;
      eSlice.visible = vs.anim === 'E';
      if (vs.anim === 'E' && lastHeld && (!eApplied || eApplied.held !== lastHeld || eApplied.E !== store.E)) {
        eApplied = { held: lastHeld, E: store.E };
        const x = ((store.E - lastHeld.E_MIN) / (lastHeld.E_MAX - lastHeld.E_MIN) - 0.5) * SX;
        const held = { dcp: lastHeld.dcp, L: lastHeld.L, rho: lastHeld.rho };
        const sp = eSliceGeo.attributes.position;
        for (let i = 0; i < N; i++) {
          const z = (i / (N - 1) - 0.5) * SZ;
          held[lastHeld.axis2] = lastHeld.a2min + (0.5 - z / SZ) * (lastHeld.a2max - lastHeld.a2min);
          const p = pMuE(lastHeld.ep, store.E, held.L, held.rho, held.dcp, lastHeld.anti);
          sp.setXYZ(i, x, (p / maxP) * SY + 0.02, z);
        }
        sp.needsUpdate = true;
      }
    }

    function probe(event) {
      if (!lastHeld) return null;
      const hit = base.raycast(event, surface);
      if (!hit) {
        // no hover -> live status at the shared sliders (updates as the animation drives them)
        const p = pMuE(lastHeld.ep, store.E, store.L, store.rho, store.dcp, lastHeld.anti);
        return `E ${store.E.toFixed(2)} GeV · L ${store.L.toFixed(0)} km · δCP ${store.dcp.toFixed(0)}° · P ${p.toFixed(4)}`;
      }
      const E = lastHeld.E_MIN + (hit.point.x / SX + 0.5) * (lastHeld.E_MAX - lastHeld.E_MIN);
      const v = lastHeld.a2min + (0.5 - hit.point.z / SZ) * (lastHeld.a2max - lastHeld.a2min);
      const held = { dcp: lastHeld.dcp, L: lastHeld.L, rho: lastHeld.rho };
      held[lastHeld.axis2] = v;
      const p = pMuE(lastHeld.ep, E, held.L, held.rho, held.dcp, lastHeld.anti);
      return `E ${E.toFixed(2)} GeV · ${AXIS_LABEL[lastHeld.axis2]} ${v.toFixed(1)} · P ${p.toFixed(4)} · maxP ${maxP.toFixed(4)}`;
    }

    return { base, update, tick, probe, dispose: () => base.dispose() };
  },

  companion: {
    title: (store) => `Spectrum at L = ${store.L} km`,
    markerDriven: true,
    draw(canvas, store) {
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.clientWidth * dpr, h = canvas.clientHeight * dpr;
      canvas.width = w; canvas.height = h;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = theme().canvas; ctx.fillRect(0, 0, w, h);
      const ep = engineParams(store);
      const [E_MIN, E_MAX] = eRangeOf(store.basePreset);
      const NPT = 240;
      const series = [
        { anti: false, color: '#e0604f', label: 'ν' },
        { anti: true, color: '#5588e8', label: 'ν̅' },
      ];
      let pmax = 1e-6;
      const data = series.map((s) => {
        const ys = new Float32Array(NPT);
        for (let i = 0; i < NPT; i++) {
          const E = E_MIN + (E_MAX - E_MIN) * i / (NPT - 1);
          ys[i] = pMuE(ep, E, store.L, store.rho, store.dcp, s.anti);
          if (ys[i] > pmax) pmax = ys[i];
        }
        return ys;
      });
      const P = plot2d(ctx, w, h, dpr, { x: [E_MIN, E_MAX], y: [0, pmax * 1.08], xTitle: 'E [GeV]', yTitle: 'P(νμ→νe)' });

      // experiment beam: dashed line at the flux peak
      const beam = PRESETS[store.basePreset];
      if (beam) {
        const xp = P.X(beam.Epeak);
        ctx.strokeStyle = theme().beam;
        ctx.lineWidth = dpr;
        ctx.setLineDash([4 * dpr, 4 * dpr]);
        ctx.beginPath(); ctx.moveTo(xp, P.mt); ctx.lineTo(xp, P.mt + P.ph); ctx.stroke();
        ctx.setLineDash([]);
        ctx.font = `${10 * dpr}px ui-monospace, monospace`;
        ctx.fillStyle = theme().beamText;
        ctx.fillText(`peak ${beam.Epeak} GeV`, xp + 4 * dpr, P.mt + P.ph - 6 * dpr);
      }
      data.forEach((ys, si) => {
        ctx.strokeStyle = series[si].color;
        ctx.lineWidth = 1.8 * dpr;
        ctx.beginPath();
        for (let i = 0; i < NPT; i++) {
          const x = P.X(E_MIN + (E_MAX - E_MIN) * i / (NPT - 1));
          const y = P.Y(ys[i]);
          if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        }
        ctx.stroke();
      });
      if (store.views.oscillogram.anim === 'E') {
        const mx = P.X(store.E);
        ctx.strokeStyle = theme().hiCss;
        ctx.lineWidth = dpr;
        ctx.beginPath(); ctx.moveTo(mx, P.mt); ctx.lineTo(mx, P.mt + P.ph); ctx.stroke();
      }
      legend(ctx, dpr, P, series);
    },
  },

  // second panel: P vs dCP at the shared energy
  companion2: {
    title: (store) => `P vs δCP at E = ${store.E.toFixed(2)} GeV`,
    markerDriven: true,
    draw(canvas, store) {
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.clientWidth * dpr, h = canvas.clientHeight * dpr;
      canvas.width = w; canvas.height = h;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = theme().canvas; ctx.fillRect(0, 0, w, h);
      const ep = engineParams(store);
      const E = store.E;
      const NPT = 240;
      const series = [
        { anti: false, color: '#e0604f', label: 'ν' },
        { anti: true, color: '#5588e8', label: 'ν̅' },
      ];
      let pmax = 1e-6;
      const data = series.map((s) => {
        const ys = new Float32Array(NPT);
        for (let i = 0; i < NPT; i++) {
          const dcp = 360 * i / (NPT - 1);
          ys[i] = pMuE(ep, E, store.L, store.rho, dcp, s.anti);
          if (ys[i] > pmax) pmax = ys[i];
        }
        return ys;
      });
      const P = plot2d(ctx, w, h, dpr, { x: [0, 360], y: [0, pmax * 1.08], xTitle: 'δCP [deg]', yTitle: 'P(νμ→νe)' });

      // dashed line at the current dCP slider value
      const xp = P.X(store.dcp);
      ctx.strokeStyle = theme().beam;
      ctx.lineWidth = dpr;
      ctx.setLineDash([4 * dpr, 4 * dpr]);
      ctx.beginPath(); ctx.moveTo(xp, P.mt); ctx.lineTo(xp, P.mt + P.ph); ctx.stroke();
      ctx.setLineDash([]);

      data.forEach((ys, si) => {
        ctx.strokeStyle = series[si].color;
        ctx.lineWidth = 1.8 * dpr;
        ctx.beginPath();
        for (let i = 0; i < NPT; i++) {
          const x = P.X(360 * i / (NPT - 1));
          const y = P.Y(ys[i]);
          if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        }
        ctx.stroke();
      });
      legend(ctx, dpr, P, series);
    },
  },
};
