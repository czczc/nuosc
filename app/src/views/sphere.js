// State sphere: two-flavor subspace Bloch projections of the exact 3-flavor state,
// Jacobi amplitude engine. Ported from prototypes/state-sphere.html (ticket #5).
// The initial flavor (from the shared channel: νμ for beam channels, νe for
// reactor survival) sits at the south pole; the two other flavors are the north
// poles. Bloch coords map to world as (bx, bz, by) -> (x, y, z).
import * as THREE from 'three';
import { SceneBase, viridis, textSprite } from '../three/SceneBase.js';
import { C, hamiltonian, eigH, amp, prob } from '../engines/jacobi.js';
import { engineParams, DEG, eRangeOf, lRangeOf, CHANNELS, eUnitOf, fmtE } from '../engines/constants.js';
import { theme } from '../theme.js';
import { plot2d, legend } from './plot2d.js';

const R = 2.5;                 // sphere radius in world units
const NSAMP = 1024;            // trajectory buffer size (L sweep; E/dcp use 512 via setDrawRange)
const PERIOD_S = 10;           // play loop: marker sweeps the full range over ~10 s
const SWEEP_LABEL = { L: 'L [km]', E: 'E [GeV]', dcp: 'δCP [°]' };
// E sweeps display in the experiment's unit (MeV for reactors)
const sweepLabel = (sweep, store) => (sweep === 'E' ? `E [${eUnitOf(store.basePreset).unit}]` : SWEEP_LABEL[sweep]);

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

// Bloch vectors of the (other flavor)–(initial flavor) subspaces, initial flavor ai:
// b = ( 2 Re(a_o conj(a_i)), 2 Im(a_o conj(a_i)), |a_o|^2 - |a_i|^2 ), o = the two others
function blochPair(eig, L, ai, o1, o2) {
  const amps = [amp(eig, ai, 0, L), amp(eig, ai, 1, L), amp(eig, ai, 2, L)];
  const aI = amps[ai];
  const crA = C.mul(amps[o1], C.conj(aI)), crB = C.mul(amps[o2], C.conj(aI));
  return {
    bA: [2 * crA[0], 2 * crA[1], C.abs2(amps[o1]) - C.abs2(aI)],
    bB: [2 * crB[0], 2 * crB[1], C.abs2(amps[o2]) - C.abs2(aI)],
    P: [C.abs2(amps[0]), C.abs2(amps[1]), C.abs2(amps[2])],
  };
}

// per-flavor colors (indexed 0 = e, 1 = mu, 2 = tau), matching the companion curves
const FLAV = [
  { name: 'νe', color: 0xe05545, css: '#e05545', tint: [0.88, 0.33, 0.27] },
  { name: 'νμ', color: 0x4488ee, css: '#4488ee', tint: [0.27, 0.53, 0.93] },
  { name: 'ντ', color: 0x44aa55, css: '#44aa55', tint: [0.27, 0.67, 0.33] },
];
// north-pole flavors for each initial flavor: the two others, in display order
const OTHERS = { 0: [1, 2], 1: [0, 2] };
const initOf = (store) => CHANNELS[store.channel].a;

export default {
  id: 'sphere',
  label: 'Statesphere',
  note: 'Two-flavor subspace projections of the 3-flavor state (8-D full space): the initial flavor sits at the south pole, the two other flavors at the north; arrow colors match the 2D panel. Vector length < 1 = leakage into the third flavor; a pure 2-flavor oscillation would trace a surface circle.',
  extras: [
    {
      key: 'pole', type: 'select', label: 'north pole',
      options: [
        { value: 'both', label: (s) => OTHERS[initOf(s)].map((f) => FLAV[f].name).join(' & ') },
        { value: 'A', label: (s) => FLAV[OTHERS[initOf(s)][0]].name },
        { value: 'B', label: (s) => FLAV[OTHERS[initOf(s)][1]].name },
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
    let sLabel = null; // south pole = initial flavor, rebuilt on channel change

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
        mat, // recolored when the initial flavor (channel) changes
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

    // A/B = the two north-pole subspaces (flavor identity follows the channel)
    const sets = {
      A: { traj: mkTraj(), arrow: mkArrow(FLAV[0].color) },
      B: { traj: mkTraj(), arrow: mkArrow(FLAV[2].color) },
    };

    // cache built by update(); tick/probe read the marker through it
    let markerState = null;

    function update() {
      const { sweep, lo, hi, stateAt } = makeSweep(store);
      const pole = store.views.sphere.pole;
      const ai = initOf(store);
      const [o1, o2] = OTHERS[ai];

      const poleKey = `${pole}-${ai}`;
      if (poleKey !== lastPole) {
        for (const l of poleLabels) base.scene.remove(l);
        poleLabels.length = 0;
        const addLabel = (text, y, color) => {
          const s = textSprite(text, 1, color);
          s.position.set(0, y, 0);
          base.scene.add(s);
          poleLabels.push(s);
        };
        if (pole === 'both') {
          addLabel([[FLAV[o1].name, FLAV[o1].css], [' / ', null], [FLAV[o2].name, FLAV[o2].css]], R + 0.4);
        } else if (pole === 'B') {
          addLabel(FLAV[o2].name, R + 0.4, FLAV[o2].css);
        } else {
          addLabel(FLAV[o1].name, R + 0.4, FLAV[o1].css);
        }
        if (sLabel) base.scene.remove(sLabel);
        sLabel = textSprite(FLAV[ai].name, 1, FLAV[ai].css);
        sLabel.position.set(0, -R - 0.4, 0);
        base.scene.add(sLabel);
        sets.A.arrow.mat.color.setHex(FLAV[o1].color);
        sets.B.arrow.mat.color.setHex(FLAV[o2].color);
        lastPole = poleKey;
      }

      const showA = pole !== 'B', showB = pole !== 'A';
      sets.A.traj.line.visible = sets.A.arrow.group.visible = showA;
      sets.B.traj.line.visible = sets.B.arrow.group.visible = showB;

      const n = sweep === 'L' ? NSAMP : 512;
      const posA = sets.A.traj.geo.attributes.position, colA = sets.A.traj.geo.attributes.color;
      const posB = sets.B.traj.geo.attributes.position, colB = sets.B.traj.geo.attributes.color;
      for (let i = 0; i < n; i++) {
        const t = i / (n - 1);
        const st = stateAt(lo + t * (hi - lo));
        const { bA, bB } = blochPair(st.eig, st.L, ai, o1, o2);
        posA.setXYZ(i, R * bA[0], R * bA[2], R * bA[1]);
        posB.setXYZ(i, R * bB[0], R * bB[2], R * bB[1]);
        if (pole === 'both') {
          // flavor-tinted gradients (dim -> bright along the sweep) so the two
          // trajectories stay distinguishable; single mode keeps viridis
          const k = 0.35 + 0.65 * t;
          colA.setXYZ(i, FLAV[o1].tint[0] * k, FLAV[o1].tint[1] * k, FLAV[o1].tint[2] * k);
          colB.setXYZ(i, FLAV[o2].tint[0] * k, FLAV[o2].tint[1] * k, FLAV[o2].tint[2] * k);
        } else {
          const [cr, cg, cb] = viridis(t);
          colA.setXYZ(i, cr, cg, cb);
          colB.setXYZ(i, cr, cg, cb);
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
        return { sweep, sv, o1, o2, ...blochPair(st.eig, st.L, ai, o1, o2) };
      };
    }

    function tick(dt) {
      const vs = store.views.sphere;
      if (vs.play) vs.marker = (vs.marker + dt / PERIOD_S) % 1;
      if (!markerState) return;
      const st = markerState(vs.marker);
      sets.A.arrow.set(st.bA);
      sets.B.arrow.set(st.bB);
    }

    // readout at the CURRENT marker (no raycast)
    function probe() {
      if (!markerState) return null;
      const pole = store.views.sphere.pole;
      const st = markerState(store.views.sphere.marker);
      const name = st.sweep === 'L' ? `L ${Math.round(st.sv)} km`
        : st.sweep === 'E' ? `E ${fmtE(st.sv, store.basePreset)}`
        : `δCP ${Math.round(st.sv)}°`;
      const bnA = Math.hypot(st.bA[0], st.bA[1], st.bA[2]);
      const bnB = Math.hypot(st.bB[0], st.bB[1], st.bB[2]);
      const btxt = pole === 'both' ? `|b(${FLAV[st.o1].name})| ${bnA.toFixed(4)} |b(${FLAV[st.o2].name})| ${bnB.toFixed(4)}`
        : `|b| ${(pole === 'B' ? bnB : bnA).toFixed(4)}`;
      return `${name} · Pe ${st.P[0].toFixed(4)} Pμ ${st.P[1].toFixed(4)} Pτ ${st.P[2].toFixed(4)} · ${btxt}`;
    }

    return { base, update, tick, probe, dispose: () => base.dispose() };
  },

  companion: {
    title: (store) => `P vs ${sweepLabel(store.views.sphere.sweep, store)}`,
    markerDriven: true,
    draw(canvas, store) {
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.clientWidth * dpr, h = canvas.clientHeight * dpr;
      canvas.width = w; canvas.height = h;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = theme().canvas; ctx.fillRect(0, 0, w, h);

      const { sweep, lo, hi, stateAt } = makeSweep(store);
      const ai = initOf(store);
      const scale = sweep === 'E' ? eUnitOf(store.basePreset).scale : 1; // x axis in display units
      const NPT = 200;
      const series = [
        { beta: 0, color: FLAV[0].css, label: 'Pe' },
        { beta: 1, color: FLAV[1].css, label: 'Pμ' },
        { beta: 2, color: FLAV[2].css, label: 'Pτ' },
      ];
      const ys = series.map(() => new Float32Array(NPT));
      for (let i = 0; i < NPT; i++) {
        const st = stateAt(lo + (hi - lo) * i / (NPT - 1));
        for (let s = 0; s < 3; s++) ys[s][i] = prob(st.eig, ai, series[s].beta, st.L);
      }

      const P = plot2d(ctx, w, h, dpr, {
        x: [lo * scale, hi * scale], y: [0, 1],
        xTitle: sweepLabel(sweep, store), yTitle: `P(${FLAV[ai].name}→νx)`,
      });
      ys.forEach((yv, si) => {
        ctx.strokeStyle = series[si].color;
        ctx.lineWidth = 1.8 * dpr;
        ctx.beginPath();
        for (let i = 0; i < NPT; i++) {
          const x = P.X((lo + (hi - lo) * i / (NPT - 1)) * scale);
          if (i === 0) ctx.moveTo(x, P.Y(yv[i])); else ctx.lineTo(x, P.Y(yv[i]));
        }
        ctx.stroke();
      });

      // marker (markerDriven: CompanionPanel's rAF repaints this)
      const mx = P.X((lo + (hi - lo) * store.views.sphere.marker) * scale);
      ctx.strokeStyle = theme().hiCss;
      ctx.lineWidth = dpr;
      ctx.beginPath(); ctx.moveTo(mx, P.mt); ctx.lineTo(mx, P.mt + P.ph); ctx.stroke();

      legend(ctx, dpr, P, series);
    },
  },
};
