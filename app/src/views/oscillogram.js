// Oscillogram surface: P(numu->nue) over (E x second-axis), NuFast engine.
// Reference implementation of the view contract (see views/index.js).
import * as THREE from 'three';
import { SceneBase, viridis, textSprite } from '../three/SceneBase.js';
import { probabilityMatter } from '../engines/nufast.js';
import { engineParams, PRESETS, eRangeOf, lRangeOf } from '../engines/constants.js';
import { theme } from '../theme.js';
import { plot2d, legend } from './plot2d.js';

const SX = 8, SZ = 8, SY = 2.2, N = 161;
const RANGES = { dcp: [0, 360], rho: [0, 15] };
const AXIS_LABEL = { L: 'L [km]', dcp: 'dCP [deg]', rho: 'rho [g/cm3]' };

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
        { value: 'L', label: 'L [0 – 2×baseline km]' },
        { value: 'dcp', label: 'δCP [0–360°]' },
        { value: 'rho', label: 'ρ [0–15 g/cm³]' },
      ],
    },
  ],

  create(container, store) {
    const base = new SceneBase(container, { camPos: [7, 5, 8], ortho: store.ortho });

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

    const grid = new THREE.GridHelper(Math.max(SX, SZ), 10, theme().grid1, theme().grid2);
    grid.position.y = -0.01;
    base.scene.add(grid);

    let xLabel = null, lastXKey = null;
    let zLabel = null, lastZKey = null;

    let maxP = 1, lastHeld = null;

    function update() {
      const ep = engineParams(store);
      const [E_MIN, E_MAX] = eRangeOf(store.basePreset);
      const axis2 = store.views.oscillogram.axis2;
      const [a2min, a2max] = axis2 === 'L' ? lRangeOf(store.basePreset) : RANGES[axis2];
      const held = { dcp: store.dcp, L: store.L, rho: store.rho };
      const anti = store.anti;

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
        const [r, g, b] = viridis(P[i] / maxP);
        colors[3 * i] = r; colors[3 * i + 1] = g; colors[3 * i + 2] = b;
      }
      pos.needsUpdate = true;
      geo.attributes.color.needsUpdate = true;
      geo.computeVertexNormals();

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

    function probe(event) {
      if (!lastHeld) return null;
      const hit = base.raycast(event, surface);
      if (!hit) return null;
      const E = lastHeld.E_MIN + (hit.point.x / SX + 0.5) * (lastHeld.E_MAX - lastHeld.E_MIN);
      const v = lastHeld.a2min + (0.5 - hit.point.z / SZ) * (lastHeld.a2max - lastHeld.a2min);
      const held = { dcp: lastHeld.dcp, L: lastHeld.L, rho: lastHeld.rho };
      held[lastHeld.axis2] = v;
      const p = pMuE(lastHeld.ep, E, held.L, held.rho, held.dcp, lastHeld.anti);
      return `E ${E.toFixed(2)} GeV · ${AXIS_LABEL[lastHeld.axis2]} ${v.toFixed(1)} · P ${p.toFixed(4)} · maxP ${maxP.toFixed(4)}`;
    }

    return { base, update, probe, dispose: () => base.dispose() };
  },

  companion: {
    title: (store) => `Spectrum at L = ${store.L} km`,
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
      legend(ctx, dpr, P, series);
    },
  },
};
