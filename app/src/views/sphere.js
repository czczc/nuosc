// State sphere: νe–νμ subspace Bloch projection of the exact 3-flavor state,
// Jacobi amplitude engine. Ported from prototypes/state-sphere.html (ticket #5).
// Bloch coords map to world as (bx, bz, by) -> (x, y, z): νe pole up, νμ pole down.
import * as THREE from 'three';
import { SceneBase, viridis, textSprite } from '../three/SceneBase.js';
import { C, hamiltonian, eigH, amp, prob } from '../engines/jacobi.js';
import { engineParams, DEG, eRangeOf, lRangeOf } from '../engines/constants.js';
import { theme } from '../theme.js';
import { plot2d, legend } from './plot2d.js';

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
  const rho = store.rho, anti = store.anti, Lfixed = store.L, Efixed = store.E;
  if (sweep === 'L') {
    const eig = eigH(hamiltonian(ep, Efixed, rho, anti));
    return { sweep, lo: 0, hi: lRangeOf(store.basePreset)[1], stateAt: (sv) => ({ eig, L: sv }) };
  }
  if (sweep === 'E') {
    const [lo, hi] = eRangeOf(store.basePreset);
    return { sweep, lo, hi, stateAt: (sv) => ({ eig: eigH(hamiltonian(ep, sv, rho, anti)), L: Lfixed }) };
  }
  return {
    sweep, lo: 0, hi: 360,
    stateAt: (sv) => ({ eig: eigH(hamiltonian({ ...ep, delta: sv * DEG }, Efixed, rho, anti)), L: Lfixed }),
  };
}

// Bloch vectors of the νe–νμ and ντ–νμ subspaces, initial state νμ (index 1):
// b = ( 2 Re(a_t conj(a_mu)), 2 Im(a_t conj(a_mu)), |a_t|^2 - |a_mu|^2 ), t = e or tau
function blochPair(eig, L) {
  const ae = amp(eig, 1, 0, L), am = amp(eig, 1, 1, L), at = amp(eig, 1, 2, L);
  const crE = C.mul(ae, C.conj(am)), crT = C.mul(at, C.conj(am));
  return {
    bE: [2 * crE[0], 2 * crE[1], C.abs2(ae) - C.abs2(am)],
    bT: [2 * crT[0], 2 * crT[1], C.abs2(at) - C.abs2(am)],
    Pe: C.abs2(ae), Pm: C.abs2(am), Pt: C.abs2(at),
  };
}

// flavor colors, matching the companion plot's Pe/Pμ/Pτ curves
const FLAV = {
  e: { color: 0xe05545, css: '#e05545', tint: [0.88, 0.33, 0.27] },
  tau: { color: 0x44aa55, css: '#44aa55', tint: [0.27, 0.67, 0.33] },
  mu: { css: '#4488ee' },
};

export default {
  id: 'sphere',
  label: 'Statesphere',
  note: 'νe–νμ and/or ντ–νμ subspace projection of the 3-flavor state (8-D full space); arrow colors match the 2D panel (red νe, green ντ). Vector length < 1 = leakage into the third flavor; a pure 2-flavor oscillation would trace a surface circle.',
  extras: [
    {
      key: 'pole', type: 'select', label: 'north pole',
      options: [
        { value: 'both', label: 'νe & ντ' },
        { value: 'e', label: 'νe' },
        { value: 'tau', label: 'ντ' },
      ],
    },
    {
      key: 'marker', type: 'marker', label: 'animate', step: 0.002,
      select: {
        key: 'sweep',
        options: [
          { value: 'L', label: 'L' },
          { value: 'E', label: 'E' },
          { value: 'dcp', label: 'δCP' },
        ],
      },
    },
  ],

  create(container, store) {
    const base = new SceneBase(container, { camPos: [5, 3.5, 6] });

    // unit sphere: faint solid shell + subtle wireframe + equator ring
    base.scene.add(new THREE.Mesh(
      new THREE.SphereGeometry(R, 48, 32),
      new THREE.MeshLambertMaterial({ color: theme().grid1, transparent: true, opacity: 0.10, depthWrite: false })));
    base.scene.add(new THREE.Mesh(
      new THREE.SphereGeometry(R, 24, 16),
      new THREE.MeshBasicMaterial({ color: theme().grid1, wireframe: true, transparent: true, opacity: 0.35 })));
    {
      const pts = [];
      for (let i = 0; i <= 128; i++) {
        const a = (i / 128) * 2 * Math.PI;
        pts.push(new THREE.Vector3(R * Math.cos(a), 0, R * Math.sin(a)));
      }
      base.scene.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts),
        new THREE.LineBasicMaterial({ color: theme().axis, transparent: true, opacity: 0.6 })));
    }

    const poleLabels = [];
    let lastPole = null;
    const sLabel = textSprite('νμ', 1, FLAV.mu.css);
    sLabel.position.set(0, -R - 0.4, 0);
    base.scene.add(sLabel);

    // one trajectory line + one mesh arrow per subspace (GL lines are stuck at
    // 1px, so the arrow is a cylinder + cone, oriented per frame)
    const HEAD = 0.15;
    const yUp = new THREE.Vector3(0, 1, 0);

    function mkTraj() {
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(NSAMP * 3), 3));
      geo.setAttribute('color', new THREE.BufferAttribute(new Float32Array(NSAMP * 3), 3));
      const line = new THREE.Line(geo, new THREE.LineBasicMaterial({ vertexColors: true }));
      base.scene.add(line);
      return { geo, line };
    }

    function mkArrow(color) {
      const mat = new THREE.MeshBasicMaterial({ color });
      const shaft = new THREE.Mesh(new THREE.CylinderGeometry(0.02, 0.02, 1, 12), mat);
      const head = new THREE.Mesh(new THREE.ConeGeometry(0.08, HEAD, 16), mat);
      const group = new THREE.Group();
      group.add(shaft, head);
      base.scene.add(group);
      const v = new THREE.Vector3();
      return {
        group,
        set(b) {
          v.set(R * b[0], R * b[2], R * b[1]);
          const len = v.length();
          const shaftLen = Math.max(len - HEAD, 1e-4);
          shaft.scale.set(1, shaftLen, 1);
          shaft.position.y = shaftLen / 2;
          head.position.y = shaftLen + HEAD / 2;
          if (len > 1e-6) group.quaternion.setFromUnitVectors(yUp, v.normalize());
        },
      };
    }

    const sets = {
      e: { traj: mkTraj(), arrow: mkArrow(FLAV.e.color) },
      tau: { traj: mkTraj(), arrow: mkArrow(FLAV.tau.color) },
    };

    // cache built by update(); tick/probe read the marker through it
    let markerState = null;

    function update() {
      const { sweep, lo, hi, stateAt } = makeSweep(store);
      const pole = store.views.sphere.pole;

      if (pole !== lastPole) {
        for (const l of poleLabels) base.scene.remove(l);
        poleLabels.length = 0;
        const addLabel = (text, y, color) => {
          const s = textSprite(text, 1, color);
          s.position.set(0, y, 0);
          base.scene.add(s);
          poleLabels.push(s);
        };
        if (pole === 'both') {
          addLabel([['νe', FLAV.e.css], [' / ', null], ['ντ', FLAV.tau.css]], R + 0.4);
        } else if (pole === 'tau') {
          addLabel('ντ', R + 0.4, FLAV.tau.css);
        } else {
          addLabel('νe', R + 0.4, FLAV.e.css);
        }
        lastPole = pole;
      }

      const showE = pole !== 'tau', showT = pole !== 'e';
      sets.e.traj.line.visible = sets.e.arrow.group.visible = showE;
      sets.tau.traj.line.visible = sets.tau.arrow.group.visible = showT;

      const n = sweep === 'L' ? NSAMP : 512;
      const posE = sets.e.traj.geo.attributes.position, colE = sets.e.traj.geo.attributes.color;
      const posT = sets.tau.traj.geo.attributes.position, colT = sets.tau.traj.geo.attributes.color;
      for (let i = 0; i < n; i++) {
        const t = i / (n - 1);
        const st = stateAt(lo + t * (hi - lo));
        const { bE, bT } = blochPair(st.eig, st.L);
        posE.setXYZ(i, R * bE[0], R * bE[2], R * bE[1]);
        posT.setXYZ(i, R * bT[0], R * bT[2], R * bT[1]);
        if (pole === 'both') {
          // flavor-tinted gradients (dim -> bright along the sweep) so the two
          // trajectories stay distinguishable; single mode keeps viridis
          const k = 0.35 + 0.65 * t;
          colE.setXYZ(i, FLAV.e.tint[0] * k, FLAV.e.tint[1] * k, FLAV.e.tint[2] * k);
          colT.setXYZ(i, FLAV.tau.tint[0] * k, FLAV.tau.tint[1] * k, FLAV.tau.tint[2] * k);
        } else {
          const [cr, cg, cb] = viridis(t);
          colE.setXYZ(i, cr, cg, cb);
          colT.setXYZ(i, cr, cg, cb);
        }
      }
      for (const s of Object.values(sets)) {
        s.traj.geo.setDrawRange(0, n);
        s.traj.geo.attributes.position.needsUpdate = true;
        s.traj.geo.attributes.color.needsUpdate = true;
        s.traj.geo.computeBoundingSphere();
      }

      markerState = (frac) => {
        const sv = lo + frac * (hi - lo);
        const st = stateAt(sv);
        return { sweep, sv, ...blochPair(st.eig, st.L) };
      };
    }

    function tick(dt) {
      const vs = store.views.sphere;
      if (vs.play) vs.marker = (vs.marker + dt / PERIOD_S) % 1;
      if (!markerState) return;
      const st = markerState(vs.marker);
      sets.e.arrow.set(st.bE);
      sets.tau.arrow.set(st.bT);
    }

    // readout at the CURRENT marker (no raycast)
    function probe() {
      if (!markerState) return null;
      const pole = store.views.sphere.pole;
      const st = markerState(store.views.sphere.marker);
      const name = st.sweep === 'L' ? `L ${Math.round(st.sv)} km`
        : st.sweep === 'E' ? `E ${st.sv.toFixed(2)} GeV`
        : `δCP ${Math.round(st.sv)}°`;
      const bnE = Math.hypot(st.bE[0], st.bE[1], st.bE[2]);
      const bnT = Math.hypot(st.bT[0], st.bT[1], st.bT[2]);
      const btxt = pole === 'both' ? `|b(e)| ${bnE.toFixed(4)} |b(τ)| ${bnT.toFixed(4)}`
        : `|b| ${(pole === 'tau' ? bnT : bnE).toFixed(4)}`;
      return `${name} · Pe ${st.Pe.toFixed(4)} Pμ ${st.Pm.toFixed(4)} Pτ ${st.Pt.toFixed(4)} · ${btxt}`;
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
      ctx.fillStyle = theme().canvas; ctx.fillRect(0, 0, w, h);

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

      const P = plot2d(ctx, w, h, dpr, { x: [lo, hi], y: [0, 1], xTitle: SWEEP_LABEL[sweep], yTitle: 'P(νμ→νx)' });
      ys.forEach((yv, si) => {
        ctx.strokeStyle = series[si].color;
        ctx.lineWidth = 1.8 * dpr;
        ctx.beginPath();
        for (let i = 0; i < NPT; i++) {
          const x = P.X(lo + (hi - lo) * i / (NPT - 1));
          if (i === 0) ctx.moveTo(x, P.Y(yv[i])); else ctx.lineTo(x, P.Y(yv[i]));
        }
        ctx.stroke();
      });

      // marker (markerDriven: CompanionPanel's rAF repaints this)
      const mx = P.X(lo + (hi - lo) * store.views.sphere.marker);
      ctx.strokeStyle = theme().hiCss;
      ctx.lineWidth = dpr;
      ctx.beginPath(); ctx.moveTo(mx, P.mt); ctx.lineTo(mx, P.mt + P.ph); ctx.stroke();

      legend(ctx, dpr, P, series);
    },
  },
};
