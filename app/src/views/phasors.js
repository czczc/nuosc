// Eigenstate phasors: the three matter-eigenstate contributions to the shared channel's
// transition amplitude, A(x) = sum_i c_i e^{-i lam_i L}, drawn as a resultant 3D curve
// (x = swept variable, transverse (z, y) = complex plane), with head-to-tail arrows at a
// marker cross-section and the floor probability curve P = |A|^2.
// Ported from prototypes/phasor-waves.html (validated, ticket #6).
import * as THREE from 'three';
import { SceneBase, textSprite } from '../three/SceneBase.js';
import { C, hamiltonian, eigH, phasorTerms } from '../engines/jacobi.js';
import { engineParams, DEG, eRangeOf, lRangeOf, CHANNELS, pLabel, eUnitOf, fmtE } from '../engines/constants.js';
import { theme } from '../theme.js';
import { plot2d } from './plot2d.js';

const SX = 10, N = 1024, NSWEEP = 512;   // N samples for x = L (one eigH); NSWEEP for x = E / dcp (eigH per sample)
// per-channel amplitude scale (typical |A| -> radius 2.5; the numu->nue appearance
// amplitude is ~3x smaller than the others)
const AMP_SCALE = { mue: 2.5 / 0.35, mumu: 2.5 / 1.0, mutau: 2.5 / 1.0, ee: 2.5 / 1.0 };
const FLOOR_Y = -3, P_DEPTH = 4;         // floor curve at y = FLOOR_Y, toward -z (top view: up-screen), autoscaled so pmax spans P_DEPTH
const COL = [0xe0a030, 0x4cb85e, 0xb070e0]; // arm 2 is green, not blue: the floor P curve is already blue
const AXIS_LABEL = { L: 'L [km]', E: 'E [GeV]', dcp: 'δCP [deg]' };
const LOOP_S = 10;                       // ~10 s per sweep of the marker

function xOfFrac(f) { return -SX / 2 + SX * f; }

function xRange(xaxis, preset) {
  if (xaxis === 'E') return eRangeOf(preset);
  if (xaxis === 'dcp') return [0, 360];
  return lRangeOf(preset);
}

function fmtVal(xaxis, xv, preset) {
  if (xaxis === 'E') return fmtE(xv, preset);
  if (xaxis === 'dcp') return `${xv.toFixed(0)}°`;
  return `${xv.toFixed(0)} km`;
}

// Matter-eigenstate terms of the (alpha)->(beta) amplitude, sorted by eigenvalue ascending
// (stable arm color identity) with min(lam) subtracted (overall phase removal; |A|^2 unchanged).
function decompose(ep, E, rho, anti, alpha, beta) {
  const terms = phasorTerms(eigH(hamiltonian(ep, E, rho, anti)), alpha, beta);
  terms.sort((a, b) => a.lam - b.lam);
  const lmin = terms[0].lam;
  return terms.map((t) => ({ c: t.c, lam: t.lam - lmin }));
}

// Terms at swept value xv; in the E and dcp modes the swept quantity overrides the held one.
function decomposeAt(ep, xaxis, xv, E, rho, anti, alpha, beta) {
  if (xaxis === 'E') return decompose(ep, xv, rho, anti, alpha, beta);
  if (xaxis === 'dcp') return decompose({ ...ep, delta: xv * DEG }, E, rho, anti, alpha, beta);
  return decompose(ep, E, rho, anti, alpha, beta);
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
  note: 'Phasors = matter-eigenstate contributions to the transition amplitude of the channel selected in the header; overall phase removed; arm lengths |cᵢ| depend on E and ρ.',
  extras: [
    {
      key: 'marker', type: 'marker', label: 'animate', step: 0.002,
      select: {
        key: 'xaxis',
        options: [
          { value: 'L', label: 'L' },
          { value: 'E', label: 'E' },
          { value: 'dcp', label: 'δCP' },
        ],
      },
    },
  ],

  create(container, store) {
    const base = new SceneBase(container, { camPos: [8, 5, 10] });

    function makeLine(n, color, opacity = 1) {
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.BufferAttribute(new Float32Array(n * 3), 3));
      const line = new THREE.Line(g, new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity }));
      base.scene.add(line);
      return line;
    }

    // x-axis (|A| = 0 reference) and floor grid
    const axisLine = makeLine(2, theme().axis, 1);
    axisLine.geometry.attributes.position.set([-SX / 2, 0, 0, SX / 2, 0, 0]);
    axisLine.geometry.attributes.position.needsUpdate = true;
    const grid = new THREE.GridHelper(SX + 2, 12, theme().grid1, theme().grid2);
    grid.position.y = FLOOR_Y - 0.01;
    base.scene.add(grid);

    // curves: partial sums s1, s2 and the white resultant A(x); floor P(x)
    const s1Line = makeLine(N, COL[0], 0.5);
    const s2Line = makeLine(N, COL[1], 0.5);
    const ALine = makeLine(N, theme().hi, 1);
    // floor P curve as a tube mesh (GL lines are 1px); matches the biprob rings' on-screen thickness
    const floorTube = new THREE.Mesh(new THREE.BufferGeometry(), new THREE.MeshBasicMaterial({ color: theme().pCurve }));
    base.scene.add(floorTube);

    // head-to-tail phasor arrows + white resultant arrow
    const arrows = COL.map((c) => {
      const a = new THREE.ArrowHelper(new THREE.Vector3(0, 1, 0), new THREE.Vector3(), 1, c, 0.14, 0.08);
      base.scene.add(a);
      return a;
    });
    const resArrow = new THREE.ArrowHelper(new THREE.Vector3(0, 1, 0), new THREE.Vector3(), 1, theme().hi, 0.14, 0.08);
    base.scene.add(resArrow);

    // drop-line from the cross-section to the floor probability curve, plus a dot on it
    const dropLine = makeLine(3, theme().hi, 0.4);
    const floorDot = new THREE.Mesh(new THREE.SphereGeometry(0.07, 12, 12), new THREE.MeshBasicMaterial({ color: theme().hi }));
    base.scene.add(floorDot);

    // static labels
    const imLabel = textSprite('Im', 0.8); imLabel.position.set(-SX / 2 - 0.5, 2.2, 0); base.scene.add(imLabel);
    const reLabel = textSprite('Re', 0.8); reLabel.position.set(-SX / 2 - 0.5, 0.15, 2.2); base.scene.add(reLabel);
    let pSprite = null, lastPLabelKey = null;
    let xLabel = null, lastLabelKey = null;

    let armState = null; // (frac) -> { xv, terms, pts, P }; cached by update() for tick/probe
    let pScale = P_DEPTH; // world units per unit P, set by update() from the sweep's pmax
    let ampScale = AMP_SCALE.mue; // channel-dependent |A| -> world radius, set by update()
    let last = null;     // marker snapshot for probe

    function update() {
      const ep = engineParams(store);
      const { xaxis } = store.views.phasors; // NOT play/marker (tick-only, per contract)
      const channel = store.channel;
      const { a: alpha, b: beta } = CHANNELS[channel];
      ampScale = AMP_SCALE[channel];
      const E = store.E;
      const rho = store.rho, anti = store.anti, Lfix = store.L;
      const [x0, x1] = xRange(xaxis, store.basePreset);
      const M = xaxis === 'L' ? N : NSWEEP;
      const d0 = xaxis === 'L' ? decompose(ep, E, rho, anti, alpha, beta) : null; // one eigH when x = L

      const p1 = s1Line.geometry.attributes.position;
      const p2 = s2Line.geometry.attributes.position;
      const pA = ALine.geometry.attributes.position;
      const Ps = new Float32Array(M);
      for (let i = 0; i < M; i++) {
        const f = i / (M - 1);
        const xv = x0 + f * (x1 - x0);
        const x = xOfFrac(f);
        const terms = d0 ?? decomposeAt(ep, xaxis, xv, E, rho, anti, alpha, beta);
        const pts = partials(terms, xaxis === 'L' ? xv : Lfix); // fixed shared L anchors the phase when x != L
        p1.setXYZ(i, x, ampScale * pts[1][1], ampScale * pts[1][0]);
        p2.setXYZ(i, x, ampScale * pts[2][1], ampScale * pts[2][0]);
        pA.setXYZ(i, x, ampScale * pts[3][1], ampScale * pts[3][0]);
        Ps[i] = C.abs2(pts[3]);
      }
      // autoscale the floor curve like the 2D panel (pmax + 8% headroom fills P_DEPTH)
      let pmax = 1e-6;
      for (let i = 0; i < M; i++) if (Ps[i] > pmax) pmax = Ps[i];
      pScale = P_DEPTH / (pmax * 1.08);
      const fpts = [];
      for (let i = 0; i < M; i++) fpts.push(new THREE.Vector3(xOfFrac(i / (M - 1)), FLOOR_Y, -Ps[i] * pScale));
      floorTube.geometry.dispose();
      floorTube.geometry = new THREE.TubeGeometry(new THREE.CatmullRomCurve3(fpts), 512, 0.022, 8, false);
      for (const line of [s1Line, s2Line, ALine]) {
        line.geometry.setDrawRange(0, M);
        line.geometry.attributes.position.needsUpdate = true;
      }

      const { unit, scale } = eUnitOf(store.basePreset);
      const key = xaxis === 'E'
        ? `E [${unit}]  (${+(x0 * scale).toFixed(4)} - ${+(x1 * scale).toFixed(4)})`
        : `${AXIS_LABEL[xaxis]}  (${x0} - ${x1})`;
      if (key !== lastLabelKey) {
        if (xLabel) base.scene.remove(xLabel);
        xLabel = textSprite(key);
        xLabel.position.set(0, FLOOR_Y - 0.5, 2.5);
        base.scene.add(xLabel);
        lastLabelKey = key;
      }
      const chLabel = pLabel(channel, anti);
      if (chLabel !== lastPLabelKey) {
        if (pSprite) base.scene.remove(pSprite);
        pSprite = textSprite(chLabel, 0.9);
        pSprite.position.set(-SX / 2 - 1.5, FLOOR_Y + 0.1, -1.2);
        base.scene.add(pSprite);
        lastPLabelKey = chLabel;
      }

      armState = (frac) => {
        const xv = x0 + frac * (x1 - x0);
        const terms = d0 ?? decomposeAt(ep, xaxis, xv, E, rho, anti, alpha, beta);
        const pts = partials(terms, xaxis === 'L' ? xv : Lfix);
        return { xv, terms, pts, P: C.abs2(pts[3]) };
      };
    }

    function applyMarker(frac) {
      const { xv, terms, pts, P } = armState(frac);
      const x = xOfFrac(frac);

      const v = pts.map((p) => new THREE.Vector3(x, ampScale * p[1], ampScale * p[0]));
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
      dp.setXYZ(2, x, FLOOR_Y, -P * pScale);
      dp.needsUpdate = true;
      floorDot.position.set(x, FLOOR_Y, -P * pScale);

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
      return `${name} ${fmtVal(xaxis, last.xv, store.basePreset)} · |c1| ${mags[0]} · |c2| ${mags[1]} · |c3| ${mags[2]} · P ${last.P.toFixed(4)}`;
    }

    return { base, update, tick, probe, dispose: () => base.dispose() };
  },

  companion: {
    title: (store) => `${pLabel(store.channel, store.anti)} vs ${store.views.phasors.xaxis === 'E'
      ? `E [${eUnitOf(store.basePreset).unit}]` : AXIS_LABEL[store.views.phasors.xaxis]}`,
    markerDriven: true,
    draw(canvas, store) {
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.clientWidth * dpr, h = canvas.clientHeight * dpr;
      canvas.width = w; canvas.height = h;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = theme().canvas; ctx.fillRect(0, 0, w, h);

      const ep = engineParams(store);
      const { xaxis, marker } = store.views.phasors; // marker OK here: markerDriven
      const { a: alpha, b: beta } = CHANNELS[store.channel];
      const E = store.E;
      const [x0, x1] = xRange(xaxis, store.basePreset);
      const scale = xaxis === 'E' ? eUnitOf(store.basePreset).scale : 1; // x axis in display units
      const NPT = 240;
      const d0 = xaxis === 'L' ? decompose(ep, E, store.rho, store.anti, alpha, beta) : null;
      const ys = new Float32Array(NPT);
      let pmax = 1e-6;
      for (let i = 0; i < NPT; i++) {
        const xv = x0 + (x1 - x0) * i / (NPT - 1);
        const terms = d0 ?? decomposeAt(ep, xaxis, xv, E, store.rho, store.anti, alpha, beta);
        const pts = partials(terms, xaxis === 'L' ? xv : store.L);
        ys[i] = C.abs2(pts[3]);
        if (ys[i] > pmax) pmax = ys[i];
      }

      const P = plot2d(ctx, w, h, dpr, {
        x: [x0 * scale, x1 * scale], y: [0, pmax * 1.08],
        xTitle: xaxis === 'E' ? `E [${eUnitOf(store.basePreset).unit}]` : AXIS_LABEL[xaxis],
        yTitle: pLabel(store.channel, store.anti),
      });

      ctx.strokeStyle = theme().pCurveCss;
      ctx.lineWidth = 1.8 * dpr;
      ctx.beginPath();
      for (let i = 0; i < NPT; i++) {
        const x = P.X((x0 + (x1 - x0) * i / (NPT - 1)) * scale);
        const y = P.Y(ys[i]);
        if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.stroke();

      const mx = P.X((x0 + (x1 - x0) * marker) * scale);
      ctx.strokeStyle = theme().hiCss;
      ctx.lineWidth = dpr;
      ctx.beginPath(); ctx.moveTo(mx, P.mt); ctx.lineTo(mx, P.mt + P.ph); ctx.stroke();
    },
  },
};
