// Eigenstate phasors: the three matter-eigenstate contributions to the numu->nue amplitude,
// A(x) = sum_i c_i e^{-i lam_i L}, drawn as a resultant 3D curve (x = swept variable, transverse
// (z, y) = complex plane), with head-to-tail arrows at a marker cross-section and the floor
// probability curve P = |A|^2. Ported from prototypes/phasor-waves.html (validated, ticket #6).
import * as THREE from 'three';
import { SceneBase, textSprite } from '../three/SceneBase.js';
import { C, hamiltonian, eigH, phasorTerms } from '../engines/jacobi.js';
import { engineParams, DEG } from '../engines/constants.js';

const SX = 10, N = 1024, NSWEEP = 512;   // N samples for x = L (one eigH); NSWEEP for x = E / dcp (eigH per sample)
const AMP_SCALE = 2.5 / 0.35;            // |A| = 0.35 -> radius 2.5
const FLOOR_Y = -3, P_SCALE = 8;         // floor curve: z = P * P_SCALE at y = FLOOR_Y
const COL = [0xe0a030, 0x30b0d0, 0xb070e0];
const AXIS_LABEL = { L: 'L [km]', E: 'E [GeV]', dcp: 'δCP [deg]' };
const LOOP_S = 10;                       // ~10 s per sweep of the marker

function xOfFrac(f) { return -SX / 2 + SX * f; }

function xRange(xaxis, Lmax) {
  if (xaxis === 'E') return [0.2, 6];
  if (xaxis === 'dcp') return [0, 360];
  return [0, Lmax];
}

function fmtVal(xaxis, xv) {
  if (xaxis === 'E') return `${xv.toFixed(2)} GeV`;
  if (xaxis === 'dcp') return `${xv.toFixed(0)}°`;
  return `${xv.toFixed(0)} km`;
}

// Matter-eigenstate terms of the numu->nue amplitude, sorted by eigenvalue ascending
// (stable arm color identity) with min(lam) subtracted (overall phase removal; |A|^2 unchanged).
function decompose(ep, E, rho, anti) {
  const terms = phasorTerms(eigH(hamiltonian(ep, E, rho, anti)), 1, 0);
  terms.sort((a, b) => a.lam - b.lam);
  const lmin = terms[0].lam;
  return terms.map((t) => ({ c: t.c, lam: t.lam - lmin }));
}

// Terms at swept value xv; in the E and dcp modes the swept quantity overrides the held one.
function decomposeAt(ep, xaxis, xv, E, rho, anti) {
  if (xaxis === 'E') return decompose(ep, xv, rho, anti);
  if (xaxis === 'dcp') return decompose({ ...ep, delta: xv * DEG }, E, rho, anti);
  return decompose(ep, E, rho, anti);
}

// Running phasor sums at one L: [0, s1, s2, s3 = A], each a complex [re, im].
function partials(terms, L) {
  const pts = [[0, 0]];
  let s = [0, 0];
  for (const t of terms) {
    s = C.add(s, C.mul(t.c, C.expi(-t.lam * L)));
    pts.push(s);
  }
  return pts;
}

export default {
  id: 'phasors',
  label: 'Phasors',
  note: 'Phasors = matter-eigenstate contributions to the νe appearance amplitude; overall phase removed; arm lengths |cᵢ| depend on E and ρ.',
  extras: [
    {
      key: 'xaxis', type: 'select', label: 'x axis',
      options: [
        { value: 'L', label: 'x = L' },
        { value: 'E', label: 'x = E (spectrum)' },
        { value: 'dcp', label: 'x = δCP' },
      ],
    },
    { key: 'E', type: 'range', label: 'E [GeV]', min: 0.2, max: 6, step: 0.01 },
    { key: 'Lmax', type: 'range', label: 'L max [km]', min: 500, max: 10000, step: 100 },
    { key: 'play', type: 'checkbox', label: 'play' },
    { key: 'marker', type: 'range', label: 'marker', min: 0, max: 1, step: 0.002 },
  ],

  create(container, store) {
    const base = new SceneBase(container, { camPos: [8, 5, 10], ortho: store.ortho });

    function makeLine(n, color, opacity = 1) {
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.BufferAttribute(new Float32Array(n * 3), 3));
      const line = new THREE.Line(g, new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity }));
      base.scene.add(line);
      return line;
    }

    // x-axis (|A| = 0 reference) and floor grid
    const axisLine = makeLine(2, 0x556677, 1);
    axisLine.geometry.attributes.position.set([-SX / 2, 0, 0, SX / 2, 0, 0]);
    axisLine.geometry.attributes.position.needsUpdate = true;
    const grid = new THREE.GridHelper(SX + 2, 12, 0x334455, 0x223344);
    grid.position.y = FLOOR_Y - 0.01;
    base.scene.add(grid);

    // curves: partial sums s1, s2 and the white resultant A(x); floor P(x)
    const s1Line = makeLine(N, COL[0], 0.5);
    const s2Line = makeLine(N, COL[1], 0.5);
    const ALine = makeLine(N, 0xffffff, 1);
    const floorLine = makeLine(N, 0x66ccff, 1);

    // cross-section disk (complex plane) at the marker
    const diskGeo = new THREE.CircleGeometry(1, 48);
    diskGeo.rotateY(Math.PI / 2); // normal along +x -> disk lies in the (z, y) complex plane
    const disk = new THREE.Mesh(diskGeo, new THREE.MeshBasicMaterial({
      color: 0x8899bb, transparent: true, opacity: 0.12, side: THREE.DoubleSide, depthWrite: false,
    }));
    base.scene.add(disk);

    // head-to-tail phasor arrows + white resultant arrow
    const arrows = COL.map((c) => {
      const a = new THREE.ArrowHelper(new THREE.Vector3(0, 1, 0), new THREE.Vector3(), 1, c, 0.14, 0.08);
      base.scene.add(a);
      return a;
    });
    const resArrow = new THREE.ArrowHelper(new THREE.Vector3(0, 1, 0), new THREE.Vector3(), 1, 0xffffff, 0.14, 0.08);
    base.scene.add(resArrow);

    // drop-line from the cross-section to the floor probability curve, plus a dot on it
    const dropLine = makeLine(3, 0xffffff, 0.4);
    const floorDot = new THREE.Mesh(new THREE.SphereGeometry(0.07, 12, 12), new THREE.MeshBasicMaterial({ color: 0xffffff }));
    base.scene.add(floorDot);

    // static labels
    const imLabel = textSprite('Im', 0.4); imLabel.position.set(-SX / 2 - 0.5, 2.2, 0); base.scene.add(imLabel);
    const reLabel = textSprite('Re', 0.4); reLabel.position.set(-SX / 2 - 0.5, 0.15, 2.2); base.scene.add(reLabel);
    const pLabel = textSprite('P(νμ→νe)', 0.45); pLabel.position.set(-SX / 2 - 1.5, FLOOR_Y + 0.1, 1.2); base.scene.add(pLabel);
    let xLabel = null, lastLabelKey = null;

    let armState = null; // (frac) -> { xv, terms, pts, P }; cached by update() for tick/probe
    let last = null;     // marker snapshot for probe

    function update() {
      const ep = engineParams(store);
      const { xaxis, E, Lmax } = store.views.phasors; // NOT play/marker (tick-only, per contract)
      const rho = store.rho, anti = store.anti, Lfix = store.L;
      const [x0, x1] = xRange(xaxis, Lmax);
      const M = xaxis === 'L' ? N : NSWEEP;
      const d0 = xaxis === 'L' ? decompose(ep, E, rho, anti) : null; // one eigH when x = L

      const p1 = s1Line.geometry.attributes.position;
      const p2 = s2Line.geometry.attributes.position;
      const pA = ALine.geometry.attributes.position;
      const pF = floorLine.geometry.attributes.position;
      for (let i = 0; i < M; i++) {
        const f = i / (M - 1);
        const xv = x0 + f * (x1 - x0);
        const x = xOfFrac(f);
        const terms = d0 ?? decomposeAt(ep, xaxis, xv, E, rho, anti);
        const pts = partials(terms, xaxis === 'L' ? xv : Lfix); // fixed shared L anchors the phase when x != L
        p1.setXYZ(i, x, AMP_SCALE * pts[1][1], AMP_SCALE * pts[1][0]);
        p2.setXYZ(i, x, AMP_SCALE * pts[2][1], AMP_SCALE * pts[2][0]);
        pA.setXYZ(i, x, AMP_SCALE * pts[3][1], AMP_SCALE * pts[3][0]);
        pF.setXYZ(i, x, FLOOR_Y, C.abs2(pts[3]) * P_SCALE);
      }
      for (const line of [s1Line, s2Line, ALine, floorLine]) {
        line.geometry.setDrawRange(0, M);
        line.geometry.attributes.position.needsUpdate = true;
      }

      const key = `${AXIS_LABEL[xaxis]}  (${x0} - ${x1})`;
      if (key !== lastLabelKey) {
        if (xLabel) base.scene.remove(xLabel);
        xLabel = textSprite(key);
        xLabel.position.set(0, FLOOR_Y - 0.5, 2.5);
        base.scene.add(xLabel);
        lastLabelKey = key;
      }

      armState = (frac) => {
        const xv = x0 + frac * (x1 - x0);
        const terms = d0 ?? decomposeAt(ep, xaxis, xv, E, rho, anti);
        const pts = partials(terms, xaxis === 'L' ? xv : Lfix);
        return { xv, terms, pts, P: C.abs2(pts[3]) };
      };
    }

    function applyMarker(frac) {
      const { xv, terms, pts, P } = armState(frac);
      const x = xOfFrac(frac);

      disk.position.x = x;
      const armSum = terms.reduce((t, z) => t + Math.sqrt(C.abs2(z.c)), 0);
      disk.scale.setScalar(Math.max(0.3, AMP_SCALE * armSum));

      const v = pts.map((p) => new THREE.Vector3(x, AMP_SCALE * p[1], AMP_SCALE * p[0]));
      for (let i = 0; i < 3; i++) {
        const seg = v[i + 1].clone().sub(v[i]), len = seg.length();
        arrows[i].visible = len > 1e-4;
        arrows[i].position.copy(v[i]);
        if (len > 1e-4) {
          arrows[i].setDirection(seg.normalize());
          arrows[i].setLength(len, Math.min(0.14, 0.4 * len), Math.min(0.08, 0.25 * len));
        }
      }
      const res = v[3].clone().sub(v[0]), rlen = res.length();
      resArrow.visible = rlen > 1e-4;
      resArrow.position.copy(v[0]);
      if (rlen > 1e-4) {
        resArrow.setDirection(res.normalize());
        resArrow.setLength(rlen, Math.min(0.14, 0.4 * rlen), Math.min(0.08, 0.25 * rlen));
      }

      const dp = dropLine.geometry.attributes.position;
      dp.setXYZ(0, x, 0, 0);
      dp.setXYZ(1, x, FLOOR_Y, 0);
      dp.setXYZ(2, x, FLOOR_Y, P * P_SCALE);
      dp.needsUpdate = true;
      floorDot.position.set(x, FLOOR_Y, P * P_SCALE);

      last = { xv, terms, P };
    }

    function tick(dt) {
      if (!armState) return;
      const vs = store.views.phasors;
      if (vs.play) vs.marker = (vs.marker + dt / LOOP_S) % 1;
      applyMarker(vs.marker);
    }

    function probe() {
      if (!last) return null;
      const xaxis = store.views.phasors.xaxis;
      const name = xaxis === 'dcp' ? 'δCP' : xaxis;
      const mags = last.terms.map((t) => Math.sqrt(C.abs2(t.c)).toFixed(3));
      return `${name} ${fmtVal(xaxis, last.xv)} · |c1| ${mags[0]} · |c2| ${mags[1]} · |c3| ${mags[2]} · P ${last.P.toFixed(4)}`;
    }

    return { base, update, tick, probe, dispose: () => base.dispose() };
  },

  companion: {
    title: (store) => `P(νμ→νe) vs ${AXIS_LABEL[store.views.phasors.xaxis]}`,
    markerDriven: true,
    draw(canvas, store) {
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.clientWidth * dpr, h = canvas.clientHeight * dpr;
      canvas.width = w; canvas.height = h;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = '#0b0e13'; ctx.fillRect(0, 0, w, h);

      const ep = engineParams(store);
      const { xaxis, E, Lmax, marker } = store.views.phasors; // marker OK here: markerDriven
      const [x0, x1] = xRange(xaxis, Lmax);
      const NPT = 240;
      const d0 = xaxis === 'L' ? decompose(ep, E, store.rho, store.anti) : null;
      const ys = new Float32Array(NPT);
      let pmax = 1e-6;
      for (let i = 0; i < NPT; i++) {
        const xv = x0 + (x1 - x0) * i / (NPT - 1);
        const terms = d0 ?? decomposeAt(ep, xaxis, xv, E, store.rho, store.anti);
        const pts = partials(terms, xaxis === 'L' ? xv : store.L);
        ys[i] = C.abs2(pts[3]);
        if (ys[i] > pmax) pmax = ys[i];
      }

      ctx.strokeStyle = 'rgba(140,155,170,0.15)';
      ctx.lineWidth = dpr;
      for (let g = 1; g < 5; g++) { ctx.beginPath(); ctx.moveTo(0, h * g / 5); ctx.lineTo(w, h * g / 5); ctx.stroke(); }

      ctx.strokeStyle = '#66ccff';
      ctx.lineWidth = 1.8 * dpr;
      ctx.beginPath();
      for (let i = 0; i < NPT; i++) {
        const x = w * i / (NPT - 1);
        const y = h - 6 * dpr - (h - 20 * dpr) * ys[i] / pmax;
        if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.stroke();

      const mx = w * marker;
      ctx.strokeStyle = 'rgba(255,255,255,0.8)';
      ctx.lineWidth = dpr;
      ctx.beginPath(); ctx.moveTo(mx, 0); ctx.lineTo(mx, h); ctx.stroke();

      ctx.font = `${10 * dpr}px ui-monospace, monospace`;
      ctx.fillStyle = '#9aa7b5';
      ctx.fillText(`${AXIS_LABEL[xaxis]} [${x0}–${x1}] · max ${pmax.toFixed(3)}`, 8 * dpr, 14 * dpr);
    },
  },
};
