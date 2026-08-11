// State sphere: νe–νμ subspace Bloch projection of the exact 3-flavor state,
// Jacobi amplitude engine. Ported from prototypes/state-sphere.html (ticket #5).
// Bloch coords map to world as (bx, bz, by) -> (x, y, z): νe pole up, νμ pole down.
import * as THREE from 'three';
import { SceneBase, viridis, textSprite } from '../three/SceneBase.js';
import { C, hamiltonian, eigH, amp, prob } from '../engines/jacobi.js';
import { engineParams, DEG } from '../engines/constants.js';

const R = 2.5;                 // sphere radius in world units
const NSAMP = 1024;            // trajectory buffer size (L sweep; E/dcp use 512 via setDrawRange)
const PERIOD_S = 10;           // play loop: marker sweeps the full range over ~10 s
const SWEEP_LABEL = { L: 'L [km]', E: 'E [GeV]', dcp: 'δCP [°]' };

// Sweep context from the store: range [lo, hi] and stateAt(sv) -> {eig, L}.
// sweep=L needs one eigH total; sweep=E/dcp rebuilds hamiltonian+eigH per sample.
function makeSweep(store) {
  const ep = engineParams(store);
  const vs = store.views.sphere;
  const sweep = vs.sweep;
  const rho = store.rho, anti = store.anti, Lfixed = store.L, Efixed = vs.E;
  if (sweep === 'L') {
    const eig = eigH(hamiltonian(ep, Efixed, rho, anti));
    return { sweep, lo: 0, hi: vs.Lmax, stateAt: (sv) => ({ eig, L: sv }) };
  }
  if (sweep === 'E') {
    return { sweep, lo: 0.2, hi: 6, stateAt: (sv) => ({ eig: eigH(hamiltonian(ep, sv, rho, anti)), L: Lfixed }) };
  }
  return {
    sweep, lo: 0, hi: 360,
    stateAt: (sv) => ({ eig: eigH(hamiltonian({ ...ep, delta: sv * DEG }, Efixed, rho, anti)), L: Lfixed }),
  };
}

// Bloch vector of the νe–νμ subspace, initial state νμ (index 1):
// bx = 2 Re(a_e conj(a_mu)), by = 2 Im(a_e conj(a_mu)), bz = |a_e|^2 - |a_mu|^2
function bloch(eig, L) {
  const ae = amp(eig, 1, 0, L), am = amp(eig, 1, 1, L);
  const cr = C.mul(ae, C.conj(am));
  return [2 * cr[0], 2 * cr[1], C.abs2(ae) - C.abs2(am)];
}

export default {
  id: 'sphere',
  label: 'State sphere',
  note: 'νe–νμ subspace projection of the 3-flavor state (8-D full space). Vector length < 1 = leakage into ντ; a pure 2-flavor oscillation would trace a surface circle.',
  extras: [
    {
      key: 'sweep', type: 'select', label: 'sweep',
      options: [
        { value: 'L', label: 'trajectory over L' },
        { value: 'E', label: 'over E (state at detector)' },
        { value: 'dcp', label: 'over δCP' },
      ],
    },
    { key: 'E', type: 'range', label: 'E [GeV]', min: 0.2, max: 6, step: 0.01 },
    { key: 'Lmax', type: 'range', label: 'L max [km]', min: 500, max: 20000, step: 100 },
    { key: 'play', type: 'checkbox', label: 'play' },
    { key: 'marker', type: 'range', label: 'marker', min: 0, max: 1, step: 0.002 },
  ],

  create(container, store) {
    const base = new SceneBase(container, { camPos: [5, 3.5, 6], ortho: store.ortho });

    // unit sphere: faint solid shell + subtle wireframe + equator ring
    base.scene.add(new THREE.Mesh(
      new THREE.SphereGeometry(R, 48, 32),
      new THREE.MeshLambertMaterial({ color: 0x334455, transparent: true, opacity: 0.10, depthWrite: false })));
    base.scene.add(new THREE.Mesh(
      new THREE.SphereGeometry(R, 24, 16),
      new THREE.MeshBasicMaterial({ color: 0x334455, wireframe: true, transparent: true, opacity: 0.35 })));
    {
      const pts = [];
      for (let i = 0; i <= 128; i++) {
        const a = (i / 128) * 2 * Math.PI;
        pts.push(new THREE.Vector3(R * Math.cos(a), 0, R * Math.sin(a)));
      }
      base.scene.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts),
        new THREE.LineBasicMaterial({ color: 0x556677, transparent: true, opacity: 0.6 })));
    }

    const nLabel = textSprite('νe');
    nLabel.position.set(0, R + 0.4, 0);
    base.scene.add(nLabel);
    const sLabel = textSprite('νμ');
    sLabel.position.set(0, -R - 0.4, 0);
    base.scene.add(sLabel);

    // trajectory line, colored viridis along the sweep fraction
    const trajGeo = new THREE.BufferGeometry();
    trajGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(NSAMP * 3), 3));
    trajGeo.setAttribute('color', new THREE.BufferAttribute(new Float32Array(NSAMP * 3), 3));
    const traj = new THREE.Line(trajGeo, new THREE.LineBasicMaterial({ vertexColors: true }));
    base.scene.add(traj);

    // animated state vector: arrow from origin + marker sphere at tip
    const arrow = new THREE.ArrowHelper(new THREE.Vector3(0, -1, 0), new THREE.Vector3(0, 0, 0), R, 0xff8800, 0.3, 0.15);
    base.scene.add(arrow);
    const tip = new THREE.Mesh(new THREE.SphereGeometry(0.08, 16, 12),
      new THREE.MeshBasicMaterial({ color: 0xff8800 }));
    base.scene.add(tip);

    // cache built by update(); tick/probe read the marker through it
    let markerState = null;

    function update() {
      const { sweep, lo, hi, stateAt } = makeSweep(store);
      const n = sweep === 'L' ? NSAMP : 512;
      const pos = trajGeo.attributes.position, col = trajGeo.attributes.color;
      for (let i = 0; i < n; i++) {
        const t = i / (n - 1);
        const st = stateAt(lo + t * (hi - lo));
        const b = bloch(st.eig, st.L);
        pos.setXYZ(i, R * b[0], R * b[2], R * b[1]);
        const [cr, cg, cb] = viridis(t);
        col.setXYZ(i, cr, cg, cb);
      }
      trajGeo.setDrawRange(0, n);
      pos.needsUpdate = true;
      col.needsUpdate = true;
      trajGeo.computeBoundingSphere();

      markerState = (frac) => {
        const sv = lo + frac * (hi - lo);
        const st = stateAt(sv);
        const ae = amp(st.eig, 1, 0, st.L), am = amp(st.eig, 1, 1, st.L), at = amp(st.eig, 1, 2, st.L);
        const cr = C.mul(ae, C.conj(am));
        return {
          sweep, sv,
          b: [2 * cr[0], 2 * cr[1], C.abs2(ae) - C.abs2(am)],
          Pe: C.abs2(ae), Pm: C.abs2(am), Pt: C.abs2(at),
        };
      };
    }

    function tick(dt) {
      const vs = store.views.sphere;
      if (vs.play) vs.marker = (vs.marker + dt / PERIOD_S) % 1;
      if (!markerState) return;
      const st = markerState(vs.marker);
      const v = new THREE.Vector3(R * st.b[0], R * st.b[2], R * st.b[1]);
      const len = v.length();
      if (len > 1e-6) arrow.setDirection(v.clone().normalize());
      arrow.setLength(Math.max(len, 1e-4), 0.3, 0.15);
      tip.position.copy(v);
    }

    // readout at the CURRENT marker (no raycast)
    function probe() {
      if (!markerState) return null;
      const st = markerState(store.views.sphere.marker);
      const name = st.sweep === 'L' ? `L ${Math.round(st.sv)} km`
        : st.sweep === 'E' ? `E ${st.sv.toFixed(2)} GeV`
        : `δCP ${Math.round(st.sv)}°`;
      const bn = Math.hypot(st.b[0], st.b[1], st.b[2]);
      return `${name} · Pe ${st.Pe.toFixed(4)} Pμ ${st.Pm.toFixed(4)} Pτ ${st.Pt.toFixed(4)} · |b| ${bn.toFixed(4)}`;
    }

    return { base, update, tick, probe, dispose: () => base.dispose() };
  },

  companion: {
    title: (store) => `P vs ${SWEEP_LABEL[store.views.sphere.sweep]}`,
    markerDriven: true,
    draw(canvas, store) {
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.clientWidth * dpr, h = canvas.clientHeight * dpr;
      canvas.width = w; canvas.height = h;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = '#0b0e13'; ctx.fillRect(0, 0, w, h);

      const { sweep, lo, hi, stateAt } = makeSweep(store);
      const NPT = 200;
      const series = [
        { beta: 0, color: '#e05545', label: 'Pe' },
        { beta: 1, color: '#4488ee', label: 'Pμ' },
        { beta: 2, color: '#44aa55', label: 'Pτ' },
      ];
      const ys = series.map(() => new Float32Array(NPT));
      for (let i = 0; i < NPT; i++) {
        const st = stateAt(lo + (hi - lo) * i / (NPT - 1));
        for (let s = 0; s < 3; s++) ys[s][i] = prob(st.eig, 1, series[s].beta, st.L);
      }

      ctx.strokeStyle = 'rgba(140,155,170,0.15)';
      ctx.lineWidth = dpr;
      for (let g = 1; g < 5; g++) { ctx.beginPath(); ctx.moveTo(0, h * g / 5); ctx.lineTo(w, h * g / 5); ctx.stroke(); }

      const yPix = (p) => h - 6 * dpr - (h - 20 * dpr) * p; // y axis fixed to [0, 1]
      ys.forEach((yv, si) => {
        ctx.strokeStyle = series[si].color;
        ctx.lineWidth = 1.8 * dpr;
        ctx.beginPath();
        for (let i = 0; i < NPT; i++) {
          const x = w * i / (NPT - 1);
          if (i === 0) ctx.moveTo(x, yPix(yv[i])); else ctx.lineTo(x, yPix(yv[i]));
        }
        ctx.stroke();
      });

      // marker (markerDriven: CompanionPanel's rAF repaints this)
      const mx = w * store.views.sphere.marker;
      ctx.strokeStyle = 'rgba(255,255,255,0.8)';
      ctx.lineWidth = dpr;
      ctx.beginPath(); ctx.moveTo(mx, 0); ctx.lineTo(mx, h); ctx.stroke();

      ctx.font = `${10 * dpr}px ui-monospace, monospace`;
      ctx.fillStyle = '#9aa7b5';
      ctx.fillText(`P(νμ→νx) vs ${SWEEP_LABEL[sweep]} [${lo}–${hi}]`, 8 * dpr, 14 * dpr);
      series.forEach((s, si) => {
        ctx.fillStyle = s.color;
        ctx.fillText(s.label, w - (3 - si) * 26 * dpr, 14 * dpr);
      });
    },
  },
};
