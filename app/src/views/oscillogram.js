// Oscillogram surface: P(numu->nue) over (E x second-axis), NuFast engine.
// Reference implementation of the view contract (see views/index.js).
import * as THREE from 'three';
import { SceneBase, viridis, textSprite } from '../three/SceneBase.js';
import { probabilityMatter } from '../engines/nufast.js';
import { engineParams } from '../engines/constants.js';

const SX = 8, SZ = 8, SY = 2.2, N = 161;
const E_MIN = 0.2, E_MAX = 6.0;
const RANGES = { L: [0, 5000], dcp: [0, 360], rho: [0, 15] };
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
        { value: 'L', label: 'L [0–5000 km]' },
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
    const slice = new THREE.Line(sliceGeo, new THREE.LineBasicMaterial({ color: 0xffffff }));
    base.scene.add(slice);

    const grid = new THREE.GridHelper(Math.max(SX, SZ), 10, 0x334455, 0x223344);
    grid.position.y = -0.01;
    base.scene.add(grid);

    const xLabel = textSprite(`E [GeV]  (${E_MIN} - ${E_MAX})`);
    xLabel.position.set(0, 0.1, SZ / 2 + 0.8);
    base.scene.add(xLabel);
    let zLabel = null;

    let maxP = 1, lastAxis2 = null, lastHeld = null;

    function update() {
      const ep = engineParams(store);
      const axis2 = store.views.oscillogram.axis2;
      const [a2min, a2max] = RANGES[axis2];
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

      if (axis2 !== lastAxis2) {
        if (zLabel) base.scene.remove(zLabel);
        zLabel = textSprite(`${AXIS_LABEL[axis2]}  (${a2min} - ${a2max})`);
        zLabel.position.set(-SX / 2 - 1.2, 0.1, 0);
        base.scene.add(zLabel);
        lastAxis2 = axis2;
      }
      lastHeld = { ...held, axis2, a2min, a2max, ep, anti };
    }

    function probe(event) {
      if (!lastHeld) return null;
      const hit = base.raycast(event, surface);
      if (!hit) return null;
      const E = E_MIN + (hit.point.x / SX + 0.5) * (E_MAX - E_MIN);
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
      ctx.fillStyle = '#0b0e13'; ctx.fillRect(0, 0, w, h);
      const ep = engineParams(store);
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
      ctx.strokeStyle = 'rgba(140,155,170,0.15)';
      ctx.lineWidth = dpr;
      for (let g = 1; g < 5; g++) { ctx.beginPath(); ctx.moveTo(0, h * g / 5); ctx.lineTo(w, h * g / 5); ctx.stroke(); }
      data.forEach((ys, si) => {
        ctx.strokeStyle = series[si].color;
        ctx.lineWidth = 1.8 * dpr;
        ctx.beginPath();
        for (let i = 0; i < NPT; i++) {
          const x = w * i / (NPT - 1);
          const y = h - 6 * dpr - (h - 20 * dpr) * ys[i] / pmax;
          if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        }
        ctx.stroke();
      });
      ctx.font = `${10 * dpr}px ui-monospace, monospace`;
      ctx.fillStyle = '#9aa7b5';
      ctx.fillText(`P(νμ→νe) vs E [${E_MIN}–${E_MAX} GeV] · max ${pmax.toFixed(3)}`, 8 * dpr, 14 * dpr);
      ctx.fillStyle = '#e0604f'; ctx.fillText('ν', w - 30 * dpr, 14 * dpr);
      ctx.fillStyle = '#5588e8'; ctx.fillText('ν̅', w - 16 * dpr, 14 * dpr);
    },
  },
};
