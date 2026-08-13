// Biprobability (Minakata-Nunokawa): P(numu->nue) vs P(antinumu->antinue), delta-CP rings
// stacked over E for both mass orderings. NuFast engine. Ported from the validated prototype.
import * as THREE from 'three';
import { SceneBase, textSprite } from '../three/SceneBase.js';
import { probabilityMatter } from '../engines/nufast.js';
import { engineParams, eRangeOf, lRangeOf } from '../engines/constants.js';
import { theme } from '../theme.js';
import { plot2d } from './plot2d.js';

const SP = 8, SZ = 8;
const E_FLOOR = 0.5, NE = 100, NDCP = 72, NR = 25;
// preset E window clipped to the readability floor (low-E rings blow up the shared scale)
const eSpan = (preset) => {
  const [e0, e1] = eRangeOf(preset);
  return [Math.max(E_FLOOR, e0), e1];
};
const zOfE = (E, E_MAX) => -(E / E_MAX) * SZ; // one-sided: z = 0 is E = 0; E grows toward -z as in oscillogram
const ORD_DEF = { NO: { color: 0xe07040, sign: 1 }, IO: { color: 0x4090e0, sign: -1 } };
const LOOP_S = 10; // ~10 s per sweep of the marker, as in oscillogram

function pair(ep, dm31, E, L, rho, dcpRad) {
  const args = [ep.s12sq, ep.s13sq, ep.s23sq, dcpRad, ep.dm21, dm31, L, E, rho * ep.Ye, 1];
  const Pn = probabilityMatter(...args);
  args[7] = -E;
  const Pa = probabilityMatter(...args);
  return [Math.max(0, Pn[1][0]), Math.max(0, Pa[1][0])];
}

export default {
  id: 'biprob',
  label: 'Biprobability',
  extras: [
    { key: 'showNO', type: 'checkbox', label: 'Normal ordering' },
    { key: 'showIO', type: 'checkbox', label: 'Inverted ordering' },
    { key: 'showSurf', type: 'checkbox', label: 'tube surface' },
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
  note: 'Ring = the ellipse traced as δCP sweeps 0–360° at one energy; the shared δCP slider moves the markers. E < 0.5 GeV excluded from the tube for readability.',

  create(container, store) {
    // side from +x so E (which grows toward -z) increases to the right on screen
    const base = new SceneBase(container, { camPos: [15, 11, -17], target: [4, 4, -4], snapPlanes: { side: [1, 0, 0] } });

    const axMat = new THREE.LineBasicMaterial({ color: theme().axis });
    for (const pts of [
      [[0, 0, 0], [SP + 0.6, 0, 0]], [[0, 0, 0], [0, SP + 0.6, 0]], [[0, 0, 0], [0, 0, -SZ - 0.6]],
    ]) base.scene.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts.map((p) => new THREE.Vector3(...p))), axMat));
    const xl = textSprite('P(νμ→νe)'); xl.position.set(SP + 1.6, 0, 0); base.scene.add(xl);
    const yl = textSprite('P̄(ν̄μ→ν̄e)'); yl.position.set(0, SP + 1, 0); base.scene.add(yl);
    let zl = null, lastZKey = null;

    const ORD = {};
    for (const [name, def] of Object.entries(ORD_DEF)) {
      const o = { ...def };
      const rgeo = new THREE.BufferGeometry();
      rgeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(NR * NDCP * 2 * 3), 3));
      o.rings = new THREE.LineSegments(rgeo, new THREE.LineBasicMaterial({ color: def.color, transparent: true, opacity: 0.55 }));
      base.scene.add(o.rings);
      // translucent tube surface: the dCP ellipses extruded along E (from the prototype)
      const tgeo = new THREE.BufferGeometry();
      tgeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(NE * NDCP * 3), 3));
      const idx = [];
      for (let i = 0; i < NE - 1; i++) for (let j = 0; j < NDCP; j++) {
        const a = i * NDCP + j, b = i * NDCP + (j + 1) % NDCP, c = (i + 1) * NDCP + j, d = (i + 1) * NDCP + (j + 1) % NDCP;
        idx.push(a, b, c, b, d, c);
      }
      tgeo.setIndex(idx);
      o.tube = new THREE.Mesh(tgeo, new THREE.MeshBasicMaterial({ color: def.color, transparent: true, opacity: 0.15, side: THREE.DoubleSide, depthWrite: false }));
      base.scene.add(o.tube);
      o.slice = new THREE.Mesh(new THREE.BufferGeometry(), new THREE.MeshBasicMaterial({ color: def.color }));
      base.scene.add(o.slice);
      o.marker = new THREE.Mesh(new THREE.SphereGeometry(0.16, 16, 16), new THREE.MeshBasicMaterial({ color: def.color }));
      base.scene.add(o.marker);
      ORD[name] = o;
    }
    const connector = new THREE.Line(
      new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(), new THREE.Vector3()]),
      new THREE.LineBasicMaterial({ color: theme().hi })
    );
    base.scene.add(connector);

    let readout = '';
    let maxP = 1e-6, lastGridKey = null;

    function update() {
      const ep = engineParams(store);
      const dm31mag = Math.abs(ep.dm31);
      const vs = store.views.biprob;
      const L = store.L, rho = store.rho, Es = store.E, dcpM = store.dcp * Math.PI / 180;
      const [E_MIN, E_MAX] = eSpan(store.basePreset);

      const zKey = `E [GeV] (${E_MIN} - ${E_MAX})`;
      if (zKey !== lastZKey) {
        if (zl) base.scene.remove(zl);
        zl = textSprite(zKey);
        zl.position.set(0, -0.5, -SZ - 1.6);
        base.scene.add(zl);
        lastZKey = zKey;
      }

      // the ring stack doesn't depend on E/dCP — skip the NE x NDCP grids when
      // only those sliders moved (e.g. play sweeping E)
      const gridKey = JSON.stringify([ep, L, rho, E_MIN, E_MAX]);
      if (gridKey !== lastGridKey) {
        lastGridKey = gridKey;

        // pass 1: grids + shared scale
        maxP = 1e-6;
        for (const o of Object.values(ORD)) {
          o.P = new Float32Array(NE * NDCP);
          o.Pb = new Float32Array(NE * NDCP);
          for (let i = 0; i < NE; i++) {
            const E = E_MIN + (i / (NE - 1)) * (E_MAX - E_MIN);
            for (let j = 0; j < NDCP; j++) {
              const [p, pb] = pair(ep, o.sign * dm31mag, E, L, rho, (j / NDCP) * 2 * Math.PI);
              const k = i * NDCP + j;
              o.P[k] = p; o.Pb[k] = pb;
              if (p > maxP) maxP = p;
              if (pb > maxP) maxP = pb;
            }
          }
        }

        for (const o of Object.values(ORD)) {
          const rp = o.rings.geometry.attributes.position;
          let v = 0;
          for (let r = 0; r < NR; r++) {
            const i = Math.round((r / (NR - 1)) * (NE - 1));
            const z = zOfE(E_MIN + (i / (NE - 1)) * (E_MAX - E_MIN), E_MAX);
            for (let j = 0; j < NDCP; j++) {
              const k1 = i * NDCP + j, k2 = i * NDCP + (j + 1) % NDCP;
              rp.setXYZ(v++, (o.P[k1] / maxP) * SP, (o.Pb[k1] / maxP) * SP, z);
              rp.setXYZ(v++, (o.P[k2] / maxP) * SP, (o.Pb[k2] / maxP) * SP, z);
            }
          }
          rp.needsUpdate = true;

          const tp = o.tube.geometry.attributes.position;
          for (let i = 0; i < NE; i++) {
            const z = zOfE(E_MIN + (i / (NE - 1)) * (E_MAX - E_MIN), E_MAX);
            for (let j = 0; j < NDCP; j++) {
              const k = i * NDCP + j;
              tp.setXYZ(k, (o.P[k] / maxP) * SP, (o.Pb[k] / maxP) * SP, z);
            }
          }
          tp.needsUpdate = true;
        }
      }

      for (const [name, o] of Object.entries(ORD)) {
        const vis = name === 'NO' ? vs.showNO : vs.showIO;
        o.rings.visible = o.slice.visible = o.marker.visible = vis;
        o.tube.visible = vis && vs.showSurf;

        const pts = [];
        for (let j = 0; j < NDCP; j++) {
          const [p, pb] = pair(ep, o.sign * dm31mag, Es, L, rho, (j / NDCP) * 2 * Math.PI);
          pts.push(new THREE.Vector3((p / maxP) * SP, (pb / maxP) * SP, zOfE(Es, E_MAX)));
        }
        o.slice.geometry.dispose();
        o.slice.geometry = new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pts, true), 144, 0.045, 8, true);

        [o.pM, o.pbM] = pair(ep, o.sign * dm31mag, Es, L, rho, dcpM);
        o.marker.position.set((o.pM / maxP) * SP, (o.pbM / maxP) * SP, zOfE(Es, E_MAX));
      }
      connector.visible = vs.showNO && vs.showIO;
      connector.geometry.setFromPoints([ORD.NO.marker.position.clone(), ORD.IO.marker.position.clone()]);
      readout =
        `E ${Es} GeV · δCP ${store.dcp}° · NO P ${ORD.NO.pM.toFixed(4)} P̄ ${ORD.NO.pbM.toFixed(4)}` +
        ` · IO P ${ORD.IO.pM.toFixed(4)} P̄ ${ORD.IO.pbM.toFixed(4)}` +
        ` · ΔP ${(ORD.NO.pM - ORD.IO.pM).toFixed(4)}`;
    }

    // play/marker animation: the marker drives the shared slider of the chosen variable
    let lastMarker = store.views.biprob.marker;

    function tick(dt) {
      const vs = store.views.biprob;
      if (vs.play) vs.marker = (vs.marker + dt / LOOP_S) % 1;
      if (vs.marker !== lastMarker) {
        const [v0, v1] = vs.anim === 'L' ? lRangeOf(store.basePreset)
          : vs.anim === 'E' ? eRangeOf(store.basePreset) : [0, 360];
        const v = v0 + vs.marker * (v1 - v0);
        // rounded to the shared sliders' steps
        if (vs.anim === 'L') store.L = Math.round(v / 5) * 5;
        else if (vs.anim === 'E') store.E = Math.round(v * 100) / 100;
        else store.dcp = Math.round(v);
      }
      lastMarker = vs.marker;
    }

    return { base, update, tick, probe: () => readout, dispose: () => base.dispose() };
  },

  companion: {
    title: (store) => `Biprobability at E = ${store.E} GeV`,
    height: 280,
    draw(canvas, store) {
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.clientWidth * dpr, h = canvas.clientHeight * dpr;
      canvas.width = w; canvas.height = h;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = theme().canvas; ctx.fillRect(0, 0, w, h);

      const ep = engineParams(store);
      const dm31mag = Math.abs(ep.dm31);
      const vs = store.views.biprob;
      const Es = store.E, L = store.L, rho = store.rho, dcpM = store.dcp * Math.PI / 180;

      const orders = [];
      if (vs.showNO) orders.push({ sign: 1, color: '#e07040' });
      if (vs.showIO) orders.push({ sign: -1, color: '#4090e0' });
      let maxP = 1e-6;
      for (const o of orders) {
        o.pts = [];
        for (let j = 0; j <= NDCP; j++) {
          const [p, pb] = pair(ep, o.sign * dm31mag, Es, L, rho, ((j % NDCP) / NDCP) * 2 * Math.PI);
          o.pts.push([p, pb]);
          maxP = Math.max(maxP, p, pb);
        }
        [o.pM, o.pbM] = pair(ep, o.sign * dm31mag, Es, L, rho, dcpM);
      }
      const m = maxP * 1.05;
      const P = plot2d(ctx, w, h, dpr, { x: [0, m], y: [0, m], xTitle: 'P(νμ→νe)', yTitle: 'P̄(ν̄μ→ν̄e)', square: true, xTicks: 3, yTicks: 3 });
      const X = (p) => P.X(p), Y = (p) => P.Y(p);

      ctx.strokeStyle = theme().plotFrame;
      ctx.lineWidth = dpr;
      ctx.setLineDash([4 * dpr, 4 * dpr]);
      ctx.beginPath(); ctx.moveTo(X(0), Y(0)); ctx.lineTo(X(m), Y(m)); ctx.stroke();
      ctx.setLineDash([]);

      for (const o of orders) {
        ctx.strokeStyle = o.color;
        ctx.lineWidth = 1.8 * dpr;
        ctx.beginPath();
        o.pts.forEach(([p, pb], i) => { if (i === 0) ctx.moveTo(X(p), Y(pb)); else ctx.lineTo(X(p), Y(pb)); });
        ctx.stroke();
        ctx.fillStyle = o.color;
        ctx.beginPath(); ctx.arc(X(o.pM), Y(o.pbM), 4 * dpr, 0, 2 * Math.PI); ctx.fill();
      }
      ctx.font = `${10 * dpr}px ui-monospace, monospace`;
      ctx.fillStyle = theme().ink;
      ctx.fillText('δCP sweep', P.ml + 5 * dpr, P.mt + 13 * dpr);
    },
  },
};
