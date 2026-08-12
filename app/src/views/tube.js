// Flavor tube: composition of an initial numu along the baseline, Jacobi engine.
// Port of prototypes/flavor-ribbon.html — tube mode (cross-section = pie of the three
// flavor fractions, constant radius = unitarity, nue sector anchored at angle 0) plus
// stacked-bands mode (stack top edge = unitarity check), with a pie-disk marker (tube)
// or sphere + vertical line (bands) riding the marker L.
import * as THREE from 'three';
import { SceneBase, textSprite } from '../three/SceneBase.js';
import { hamiltonian, eigH, prob } from '../engines/jacobi.js';
import { engineParams, eRangeOf, lRangeOf } from '../engines/constants.js';
import { theme } from '../theme.js';
import { plot2d, legend } from './plot2d.js';

// world extents: x = L over [-SX/2, SX/2], y = stacked flavor fraction (P=1 -> SY), z = extrusion
const SX = 10, SY = 2.5, HD = 0.15; // HD = half depth of the band extrusion
const NS = 512;                     // length samples
const ASEG = 96, SPS = ASEG / 3, TAU = 2 * Math.PI; // angular segments; SPS per sector
const TUBE_R = 1.1, TUBE_CY = SY / 2;
const DISK_R = 1.4; // slightly proud of the tube

const FLAVOR_HEX = ['#e05545', '#4488ee', '#44aa55']; // e, mu, tau
const SECTOR_RGB = FLAVOR_HEX.map((c) => new THREE.Color(c));

// vertex angles of one cross-section ring from cumulative probabilities; sector
// boundaries land exactly on vertices so they stay crisp and migrate smoothly along L
function fillRingAngles(out, cumE, cumMu, cumTau) {
  const b = [0, TAU * cumE, TAU * cumMu, TAU * cumTau];
  for (let j = 0; j <= ASEG; j++) {
    const s = Math.min(2, Math.floor(j / SPS));
    out[j] = b[s] + ((j - s * SPS) / SPS) * (b[s + 1] - b[s]);
  }
}

export default {
  id: 'tube',
  label: 'Flavor tube',
  extras: [
    {
      key: 'mode', type: 'select', label: 'display',
      options: [
        { value: 'tube', label: 'flavor tube' },
        { value: 'bands', label: 'stacked bands' },
      ],
    },
    { key: 'E', type: 'range', label: 'E [GeV]', min: (s) => eRangeOf(s.basePreset)[0], max: (s) => eRangeOf(s.basePreset)[1], step: 0.01 },
    { key: 'Lmax', type: 'range', label: 'L max [km]', min: 100, max: (s) => lRangeOf(s.basePreset)[1], step: 5 },
    { key: 'play', type: 'checkbox', label: 'play' },
    { key: 'marker', type: 'range', label: 'marker', min: 0, max: 1, step: 0.002 },
  ],

  create(container, store) {
    const base = new SceneBase(container, { camPos: [-6, 4.5, 8.5], target: [0, TUBE_CY, 0], ortho: store.ortho });

    const grid = new THREE.GridHelper(SX, 10, theme().grid1, theme().grid2);
    grid.position.y = -0.01;
    base.scene.add(grid);

    let xLabel = null, labeledLmax = null;
    const yLabel = textSprite('flavor fraction');
    yLabel.position.set(-SX / 2 - 1.4, SY / 2, 0);
    base.scene.add(yLabel);

    // ---------- stacked band meshes ----------
    // Each band = front face + back face + top strip (disjoint vertex groups: front 2*NS,
    // back 2*NS, top 2*NS). Per sample i: bottom vertex 2i, top vertex 2i+1 (front & back);
    // top strip: back edge 2i, front edge 2i+1.
    const bands = [];
    for (const hex of FLAVOR_HEX) {
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(6 * NS * 3), 3));
      const idx = [];
      const quad = (a, b, c, d) => idx.push(a, b, d, b, c, d);
      const OB = 2 * NS, OT = 4 * NS;
      for (let i = 0; i < NS - 1; i++) {
        quad(2 * i, 2 * i + 2, 2 * i + 3, 2 * i + 1);                     // front (+z out)
        quad(OB + 2 * i + 1, OB + 2 * i + 3, OB + 2 * i + 2, OB + 2 * i); // back (-z out)
        quad(OT + 2 * i, OT + 2 * i + 2, OT + 2 * i + 3, OT + 2 * i + 1); // top (+y out)
      }
      geo.setIndex(idx);
      const mesh = new THREE.Mesh(geo, new THREE.MeshLambertMaterial({ color: new THREE.Color(hex), side: THREE.DoubleSide }));
      base.scene.add(mesh);
      bands.push(mesh);
    }

    function setBand(geo, cum0, cum1) {
      const p = geo.attributes.position.array;
      const OB = 2 * NS * 3, OT = 4 * NS * 3;
      for (let i = 0; i < NS; i++) {
        const x = -SX / 2 + (i / (NS - 1)) * SX;
        const y0 = cum0[i] * SY, y1 = cum1[i] * SY;
        let o = 6 * i;                                     // front: bottom, top
        p[o] = x; p[o + 1] = y0; p[o + 2] = HD;
        p[o + 3] = x; p[o + 4] = y1; p[o + 5] = HD;
        o = OB + 6 * i;                                    // back: bottom, top
        p[o] = x; p[o + 1] = y0; p[o + 2] = -HD;
        p[o + 3] = x; p[o + 4] = y1; p[o + 5] = -HD;
        o = OT + 6 * i;                                    // top strip: back, front
        p[o] = x; p[o + 1] = y1; p[o + 2] = -HD;
        p[o + 3] = x; p[o + 4] = y1; p[o + 5] = HD;
      }
      geo.attributes.position.needsUpdate = true;
      geo.computeVertexNormals();
      geo.computeBoundingSphere();
    }

    // ---------- flavor tube mesh ----------
    // Tube along x, CONSTANT radius (constant radius = unitarity). Non-indexed
    // (duplicated vertices) with per-vertex colors so sector boundaries are crisp.
    const TUBE_NV = (NS - 1) * ASEG * 6; // 2 triangles per quad, duplicated vertices
    const tubeGeo = new THREE.BufferGeometry();
    tubeGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(TUBE_NV * 3), 3));
    tubeGeo.setAttribute('normal', new THREE.BufferAttribute(new Float32Array(TUBE_NV * 3), 3));
    tubeGeo.setAttribute('color', new THREE.BufferAttribute(new Float32Array(TUBE_NV * 3), 3));
    const tubeMesh = new THREE.Mesh(tubeGeo, new THREE.MeshLambertMaterial({ vertexColors: true, side: THREE.DoubleSide }));
    base.scene.add(tubeMesh);

    const ringPrev = new Float64Array(ASEG + 1), ringCur = new Float64Array(ASEG + 1);
    function setTube(Pe, cumMu, cumTau) {
      const p = tubeGeo.attributes.position.array;
      const n = tubeGeo.attributes.normal.array;
      const c = tubeGeo.attributes.color.array;
      let rA = ringPrev, rB = ringCur, o = 0;
      fillRingAngles(rA, Pe[0], cumMu[0], cumTau[0]);
      const put = (x, th, col) => {
        const cy = Math.cos(th), sz = Math.sin(th);
        p[o] = x; p[o + 1] = TUBE_CY + TUBE_R * cy; p[o + 2] = TUBE_R * sz;
        n[o] = 0; n[o + 1] = cy; n[o + 2] = sz;
        c[o] = col.r; c[o + 1] = col.g; c[o + 2] = col.b;
        o += 3;
      };
      for (let i = 0; i < NS - 1; i++) {
        const x0 = -SX / 2 + (i / (NS - 1)) * SX;
        const x1 = -SX / 2 + ((i + 1) / (NS - 1)) * SX;
        fillRingAngles(rB, Pe[i + 1], cumMu[i + 1], cumTau[i + 1]);
        for (let j = 0; j < ASEG; j++) {
          const col = SECTOR_RGB[Math.floor(j / SPS)];
          put(x0, rA[j], col); put(x1, rB[j + 1], col); put(x1, rB[j], col);
          put(x0, rA[j], col); put(x0, rA[j + 1], col); put(x1, rB[j + 1], col);
        }
        [rA, rB] = [rB, rA];
      }
      tubeGeo.attributes.position.needsUpdate = true;
      tubeGeo.attributes.normal.needsUpdate = true;
      tubeGeo.attributes.color.needsUpdate = true;
      tubeGeo.computeBoundingSphere();
    }

    // ---------- marker pie disk (tube mode): instantaneous composition at marker L ----------
    const diskGeo = new THREE.BufferGeometry();
    diskGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(ASEG * 9), 3));
    diskGeo.setAttribute('color', new THREE.BufferAttribute(new Float32Array(ASEG * 9), 3));
    const diskMesh = new THREE.Mesh(diskGeo, new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide }));
    diskMesh.position.set(0, TUBE_CY, 0);
    base.scene.add(diskMesh);
    const diskAngles = new Float64Array(ASEG + 1);
    function setDisk(pe, pm, pt) {
      fillRingAngles(diskAngles, pe, pe + pm, pe + pm + pt);
      const p = diskGeo.attributes.position.array, c = diskGeo.attributes.color.array;
      let o = 0;
      for (let j = 0; j < ASEG; j++) {
        const col = SECTOR_RGB[Math.floor(j / SPS)];
        const a0 = diskAngles[j], a1 = diskAngles[j + 1];
        for (const [cy, sz] of [[0, 0], [Math.cos(a0), Math.sin(a0)], [Math.cos(a1), Math.sin(a1)]]) {
          p[o] = 0; p[o + 1] = DISK_R * cy; p[o + 2] = DISK_R * sz;
          c[o] = col.r; c[o + 1] = col.g; c[o + 2] = col.b;
          o += 3;
        }
      }
      diskGeo.attributes.position.needsUpdate = true;
      diskGeo.attributes.color.needsUpdate = true;
      diskGeo.computeBoundingSphere();
    }

    // ---------- marker (bands mode): sphere on top of the stack + vertical line ----------
    const markerSphere = new THREE.Mesh(new THREE.SphereGeometry(0.08, 16, 12),
      new THREE.MeshLambertMaterial({ color: theme().hi, emissive: 0x666666 }));
    base.scene.add(markerSphere);
    const markerLineGeo = new THREE.BufferGeometry();
    markerLineGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(2 * 3), 3));
    const markerLine = new THREE.Line(markerLineGeo, new THREE.LineBasicMaterial({ color: theme().hi }));
    base.scene.add(markerLine);

    let eig = null, lastLmax = 5000;
    const Pe = new Float64Array(NS), cumMu = new Float64Array(NS), cumTau = new Float64Array(NS);
    const zero = new Float64Array(NS);

    function update() {
      const st = store.views.tube; // reads mode/E/Lmax only — never play/marker (tick-only)
      const ep = engineParams(store);
      lastLmax = st.Lmax;
      eig = eigH(hamiltonian(ep, st.E, store.rho, store.anti));
      for (let i = 0; i < NS; i++) {
        const L = (i / (NS - 1)) * lastLmax;
        const pe = prob(eig, 1, 0, L);
        const pm = prob(eig, 1, 1, L);
        const pt = prob(eig, 1, 2, L);
        Pe[i] = pe; cumMu[i] = pe + pm; cumTau[i] = pe + pm + pt; // top edge = unitarity check
      }
      setBand(bands[0].geometry, zero, Pe);
      setBand(bands[1].geometry, Pe, cumMu);
      setBand(bands[2].geometry, cumMu, cumTau);
      setTube(Pe, cumMu, cumTau);

      if (lastLmax !== labeledLmax) {
        if (xLabel) base.scene.remove(xLabel);
        xLabel = textSprite(`L [km]  (0 - ${lastLmax})`);
        xLabel.position.set(0, -0.4, HD + 1.2);
        base.scene.add(xLabel);
        labeledLmax = lastLmax;
      }

      const tube = st.mode === 'tube';
      tubeMesh.visible = tube;
      diskMesh.visible = tube;
      for (const m of bands) m.visible = !tube;
      markerSphere.visible = !tube;
      markerLine.visible = !tube;
      yLabel.visible = !tube; // the "flavor fraction" y axis only applies to the stacked bands
    }

    function tick(dt) {
      const st = store.views.tube;
      if (st.play) st.marker = (st.marker + dt / 10) % 1; // ~10 s per traversal
      if (!eig) return;
      const L = st.marker * st.Lmax;
      const pe = prob(eig, 1, 0, L), pm = prob(eig, 1, 1, L), pt = prob(eig, 1, 2, L);
      const sum = pe + pm + pt;
      const x = -SX / 2 + st.marker * SX;
      diskMesh.position.x = x;
      setDisk(pe, pm, pt);
      markerSphere.position.set(x, sum * SY + 0.06, 0);
      const lp = markerLineGeo.attributes.position.array;
      lp[0] = x; lp[1] = 0; lp[2] = HD + 0.01;
      lp[3] = x; lp[4] = sum * SY; lp[5] = HD + 0.01;
      markerLineGeo.attributes.position.needsUpdate = true;
      markerLineGeo.computeBoundingSphere();
    }

    function probe(event) {
      if (!eig) return null;
      let hit = null;
      if (store.views.tube.mode === 'tube') {
        hit = base.raycast(event, tubeMesh);
      } else {
        for (const m of bands) {
          const h = base.raycast(event, m);
          if (h && (!hit || h.distance < hit.distance)) hit = h;
        }
      }
      if (!hit) return null;
      const L = Math.max(0, Math.min(1, hit.point.x / SX + 0.5)) * lastLmax;
      const pe = prob(eig, 1, 0, L), pm = prob(eig, 1, 1, L), pt = prob(eig, 1, 2, L);
      return `L ${L.toFixed(0)} km · Pe ${pe.toFixed(4)} Pμ ${pm.toFixed(4)} Pτ ${pt.toFixed(4)}`;
    }

    return { base, update, tick, probe, dispose: () => base.dispose() };
  },

  companion: {
    title: (store) => `Flavor fractions vs L at E = ${store.views.tube.E} GeV`,
    markerDriven: true, // repainted per frame by CompanionPanel's rAF, so draw may read marker
    draw(canvas, store) {
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.clientWidth * dpr, h = canvas.clientHeight * dpr;
      canvas.width = w; canvas.height = h;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = theme().canvas; ctx.fillRect(0, 0, w, h);

      const st = store.views.tube;
      const ep = engineParams(store);
      const eig = eigH(hamiltonian(ep, st.E, store.rho, store.anti));
      const NPT = 240;
      const cums = [new Float64Array(NPT), new Float64Array(NPT), new Float64Array(NPT)];
      for (let i = 0; i < NPT; i++) {
        const L = (i / (NPT - 1)) * st.Lmax;
        const pe = prob(eig, 1, 0, L), pm = prob(eig, 1, 1, L), pt = prob(eig, 1, 2, L);
        cums[0][i] = pe; cums[1][i] = pe + pm; cums[2][i] = pe + pm + pt;
      }

      const P = plot2d(ctx, w, h, dpr, { x: [0, st.Lmax], y: [0, 1], xTitle: 'L [km]', yTitle: 'P(νμ→να)' });
      const xOf = (i) => P.X((i / (NPT - 1)) * st.Lmax);
      for (let k = 0; k < 3; k++) {
        const lo = k === 0 ? null : cums[k - 1], hi = cums[k];
        ctx.fillStyle = FLAVOR_HEX[k];
        ctx.beginPath();
        for (let i = 0; i < NPT; i++) {
          const y = P.Y(hi[i]);
          if (i === 0) ctx.moveTo(xOf(0), y); else ctx.lineTo(xOf(i), y);
        }
        for (let i = NPT - 1; i >= 0; i--) ctx.lineTo(xOf(i), P.Y(lo ? lo[i] : 0));
        ctx.closePath();
        ctx.fill();
      }
      P.frame(); // the stacked fills cover the frame edges

      // white vertical marker line at the marker fraction (allowed: markerDriven)
      const mx = P.X(st.marker * st.Lmax);
      ctx.strokeStyle = theme().hiCss;
      ctx.lineWidth = dpr;
      ctx.beginPath(); ctx.moveTo(mx, P.mt); ctx.lineTo(mx, P.mt + P.ph); ctx.stroke();

      legend(ctx, dpr, P, [
        { label: 'νe', color: FLAVOR_HEX[0] },
        { label: 'νμ', color: FLAVOR_HEX[1] },
        { label: 'ντ', color: FLAVOR_HEX[2] },
      ]);
    },
  },
};
