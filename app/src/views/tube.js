// Flavor tube: composition of an initial numu along the baseline, Jacobi engine.
// Port of prototypes/flavor-ribbon.html — tube mode (cross-section = pie of the three
// flavor fractions, constant radius = unitarity, nue sector anchored at angle 0) plus
// stacked-bands mode (stack top edge = unitarity check), with a pie-disk marker (tube)
// or sphere + vertical line (bands) riding the marker L.
import * as THREE from 'three';
import { SceneBase, textSprite } from '../three/SceneBase.js';
import { hamiltonian, eigH, prob } from '../engines/jacobi.js';
import { engineParams, lRangeOf, eRangeOf } from '../engines/constants.js';
import { theme } from '../theme.js';
import { plot2d, legend } from './plot2d.js';

// world extents: y = stacked flavor fraction (P=1 -> SY); tube mode runs x = L over
// [-SX/2, SX/2] and slides along z = E (toward -z, biprob convention); bands mode
// (the "cube") matches the oscillogram's axes instead: x = E over the beam window,
// z = L toward -z over [DE/2, -DE/2], with the stacked slice plane sliding along x
const SX = 10, SY = 2.5;
const DE = 6;
const xOfE = (E, E0, E1) => -SX / 2 + ((E - E0) / (E1 - E0)) * SX;
const zOfE = (E, E0, E1) => DE / 2 - ((E - E0) / (E1 - E0)) * DE;
const zOfL = (L, Lmax) => DE / 2 - (L / Lmax) * DE;
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
  label: 'Flavortube',
  extras: [
    {
      key: 'mode', type: 'select', label: 'display',
      options: [
        { value: 'tube', label: 'flavor tube' },
        { value: 'bands', label: 'stacked bands' },
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
    const base = new SceneBase(container, { camPos: [-6, 4.5, 8.5], target: [0, TUBE_CY, 0] });

    const grid = new THREE.GridHelper(SX, 10, theme().grid1, theme().grid2);
    grid.position.y = -0.01;
    base.scene.add(grid);

    let xLabel = null, labeledLmax = null;
    let eLabel = null, labeledE = null;
    const yLabel = textSprite('flavor fraction');
    yLabel.position.set(-SX / 2 - 1.4, SY / 2, -DE / 2 - 0.4); // back-left: clear of the L label in the top view
    base.scene.add(yLabel);

    // ---------- slice bands: zero-thickness translucent stacked plane at the shared E ----------
    // Per L sample i: bottom vertex 2i, top vertex 2i+1, in the local x = 0 plane; the
    // meshes slide along x to the cube's position for the shared E.
    const bands = [];
    for (const hex of FLAVOR_HEX) {
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(2 * NS * 3), 3));
      const idx = [];
      for (let i = 0; i < NS - 1; i++) idx.push(2 * i, 2 * i + 2, 2 * i + 3, 2 * i, 2 * i + 3, 2 * i + 1);
      geo.setIndex(idx);
      const mesh = new THREE.Mesh(geo, new THREE.MeshBasicMaterial({
        color: new THREE.Color(hex), transparent: true, opacity: 0.5, side: THREE.DoubleSide, depthWrite: false,
      }));
      base.scene.add(mesh);
      bands.push(mesh);
    }

    function setBand(geo, cum0, cum1) {
      const p = geo.attributes.position.array;
      for (let i = 0; i < NS; i++) {
        const z = DE / 2 - (i / (NS - 1)) * DE; // L toward -z
        const o = 6 * i; // bottom, top
        p[o] = 0; p[o + 1] = cum0[i] * SY; p[o + 2] = z;
        p[o + 3] = 0; p[o + 4] = cum1[i] * SY; p[o + 5] = z;
      }
      geo.attributes.position.needsUpdate = true;
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

    // ---------- bands-mode "cube": flavor boundary surfaces over the full (L, E) window ----------
    // red surface = top of the nue share, blue = top of the numu share; the gap from the blue
    // surface up to the box top (fraction 1, unitarity) is the nutau share.
    const NL = 128, NE = 48;
    function mkSurf(hex) {
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(NL * NE * 3), 3));
      const idx = [];
      for (let i = 0; i < NE - 1; i++) for (let j = 0; j < NL - 1; j++) {
        const a = i * NL + j, b = a + 1, c = a + NL, d = c + 1;
        idx.push(a, b, c, b, d, c);
      }
      geo.setIndex(idx);
      const mesh = new THREE.Mesh(geo, new THREE.MeshLambertMaterial({
        color: new THREE.Color(hex), transparent: true, opacity: 0.55, side: THREE.DoubleSide, depthWrite: false,
      }));
      base.scene.add(mesh);
      return mesh;
    }
    const surfE = mkSurf(FLAVOR_HEX[0]);
    const surfMu = mkSurf(FLAVOR_HEX[1]);
    const boxEdges = new THREE.LineSegments(new THREE.EdgesGeometry(new THREE.BoxGeometry(SX, SY, DE)),
      new THREE.LineBasicMaterial({ color: theme().grid1 }));
    boxEdges.position.y = SY / 2;
    base.scene.add(boxEdges);
    let lastSurfKey = null;

    let eig = null, lastLmax = 5000;
    const Pe = new Float64Array(NS), cumMu = new Float64Array(NS), cumTau = new Float64Array(NS);
    const zero = new Float64Array(NS);

    function update() {
      const st = store.views.tube; // reads mode only — never play/marker (tick-only)
      const ep = engineParams(store);
      lastLmax = lRangeOf(store.basePreset)[1]; // tube always spans the experiment's 0 - 2L
      eig = eigH(hamiltonian(ep, store.E, store.rho, store.anti));
      for (let i = 0; i < NS; i++) {
        const L = (i / (NS - 1)) * lastLmax;
        const pe = prob(eig, 1, 0, L);
        const pm = prob(eig, 1, 1, L);
        const pt = prob(eig, 1, 2, L);
        Pe[i] = pe; cumMu[i] = pe + pm; cumTau[i] = pe + pm + pt; // top edge = unitarity check
      }
      const tube = st.mode === 'tube';
      if (tube) {
        setTube(Pe, cumMu, cumTau);
        // the tube rides the E axis: it slides toward -z as the shared E grows
        const [E0, E1] = eRangeOf(store.basePreset);
        tubeMesh.position.z = zOfE(store.E, E0, E1);
        diskMesh.position.z = tubeMesh.position.z;
      } else {
        // slice slab at the shared E, parked on the cube's E axis
        setBand(bands[0].geometry, zero, Pe);
        setBand(bands[1].geometry, Pe, cumMu);
        setBand(bands[2].geometry, cumMu, cumTau);
        const [E0, E1] = eRangeOf(store.basePreset);
        for (const m of bands) m.position.set(xOfE(store.E, E0, E1), 0, 0);

        // boundary surfaces over the full (L, E) window; independent of store.E, so
        // memoized — an E animation only moves the slice through a fixed cube
        const surfKey = JSON.stringify([ep, store.rho, store.anti, lastLmax, E0, E1]);
        if (surfKey !== lastSurfKey) {
          const pE = surfE.geometry.attributes.position;
          const pM = surfMu.geometry.attributes.position;
          for (let i = 0; i < NE; i++) {
            const Ei = E0 + (i / (NE - 1)) * (E1 - E0);
            const x = xOfE(Ei, E0, E1);
            const eg = eigH(hamiltonian(ep, Ei, store.rho, store.anti));
            for (let j = 0; j < NL; j++) {
              const L = (j / (NL - 1)) * lastLmax;
              const z = zOfL(L, lastLmax);
              const pe = prob(eg, 1, 0, L), pm = prob(eg, 1, 1, L);
              pE.setXYZ(i * NL + j, x, pe * SY, z);
              pM.setXYZ(i * NL + j, x, (pe + pm) * SY, z);
            }
          }
          for (const s of [surfE, surfMu]) {
            s.geometry.attributes.position.needsUpdate = true;
            s.geometry.computeVertexNormals();
            s.geometry.computeBoundingSphere();
          }
          lastSurfKey = surfKey;
        }
      }

      if (lastLmax !== labeledLmax) {
        if (xLabel) base.scene.remove(xLabel);
        xLabel = textSprite(`L [km]  (0 - ${lastLmax})`);
        base.scene.add(xLabel);
        labeledLmax = lastLmax;
      }
      // axis labels: tube mode has L along x (front) and E along z (left);
      // bands mode follows the oscillogram's layout — E along x (front), L along z (left)
      if (tube) xLabel.position.set(0, -0.4, DE / 2 + 1.2);
      else xLabel.position.set(-SX / 2 - 1.6, -0.4, 0);
      const eKey = `E [GeV]  (${eRangeOf(store.basePreset).join(' - ')})`;
      if (eKey !== labeledE) {
        if (eLabel) base.scene.remove(eLabel);
        eLabel = textSprite(eKey);
        base.scene.add(eLabel);
        labeledE = eKey;
      }
      if (tube) eLabel.position.set(-SX / 2 - 1.6, -0.4, 0);
      else eLabel.position.set(0, -0.4, DE / 2 + 1.0);

      tubeMesh.visible = tube;
      diskMesh.visible = tube;
      for (const m of bands) m.visible = !tube;
      markerSphere.visible = !tube;
      markerLine.visible = !tube;
      surfE.visible = surfMu.visible = boxEdges.visible = !tube;
      yLabel.visible = !tube; // the "flavor fraction" y axis only applies to the stacked bands
    }

    let lastMarker = null, lastML = 0;
    function tick(dt) {
      const st = store.views.tube;
      if (st.play) st.marker = (st.marker + dt / 10) % 1; // ~10 s per sweep
      if (!eig) return;
      let frac = st.marker; // anim = L: the marker rides the tube
      if (st.anim !== 'L') {
        // anim = E / δCP: the marker drives the shared slider (tube morphs) and the
        // disk parks at the experiment baseline, showing the composition at the detector
        if (st.marker !== lastMarker) {
          if (st.anim === 'E') {
            const [v0, v1] = eRangeOf(store.basePreset);
            store.E = Math.round((v0 + st.marker * (v1 - v0)) * 100) / 100;
          } else {
            store.dcp = Math.round(st.marker * 360);
          }
        }
        frac = Math.max(0, Math.min(1, store.L / lastLmax));
      }
      lastMarker = st.marker;
      const L = frac * lastLmax;
      const pe = prob(eig, 1, 0, L), pm = prob(eig, 1, 1, L), pt = prob(eig, 1, 2, L);
      const sum = pe + pm + pt;
      lastML = L; // marker L for the no-hover status chip
      diskMesh.position.x = -SX / 2 + frac * SX; // tube-mode marker: L along x
      setDisk(pe, pm, pt);
      const [E0, E1] = eRangeOf(store.basePreset);
      const slabX = xOfE(store.E, E0, E1); // bands-mode marker rides the slice plane, L along z
      const zM = zOfL(frac * lastLmax, lastLmax);
      markerSphere.position.set(slabX, sum * SY + 0.06, zM);
      const lp = markerLineGeo.attributes.position.array;
      lp[0] = slabX + 0.02; lp[1] = 0; lp[2] = zM;
      lp[3] = slabX + 0.02; lp[4] = sum * SY; lp[5] = zM;
      markerLineGeo.attributes.position.needsUpdate = true;
      markerLineGeo.computeBoundingSphere();
    }

    function probe(event) {
      if (!eig) return null;
      const tube = store.views.tube.mode === 'tube';
      let hit = null;
      if (tube) {
        hit = base.raycast(event, tubeMesh);
      } else {
        for (const m of [...bands, surfE, surfMu]) {
          const h = base.raycast(event, m);
          if (h && (!hit || h.distance < hit.distance)) hit = h;
        }
      }
      // no hover -> live status at the animation marker (the chip re-polls while playing)
      let L = lastML, Eh = store.E, eg = eig;
      if (hit && tube) {
        L = Math.max(0, Math.min(1, hit.point.x / SX + 0.5)) * lastLmax;
      } else if (hit) {
        L = Math.max(0, Math.min(1, 0.5 - hit.point.z / DE)) * lastLmax;
        // slab hit -> the shared E; surface hit -> decode E from the hit's x
        if (!bands.includes(hit.object)) {
          const [E0, E1] = eRangeOf(store.basePreset);
          Eh = Math.max(E0, Math.min(E1, E0 + (hit.point.x / SX + 0.5) * (E1 - E0)));
          eg = eigH(hamiltonian(engineParams(store), Eh, store.rho, store.anti));
        }
      }
      const pe = prob(eg, 1, 0, L), pm = prob(eg, 1, 1, L), pt = prob(eg, 1, 2, L);
      return `L ${L.toFixed(0)} km · E ${Eh.toFixed(2)} GeV · Pe ${pe.toFixed(4)} Pμ ${pm.toFixed(4)} Pτ ${pt.toFixed(4)}`;
    }

    return { base, update, tick, probe, dispose: () => base.dispose() };
  },

  companion: {
    title: (store) => `Flavor fractions vs L at E = ${store.E} GeV`,
    markerDriven: true, // repainted per frame by CompanionPanel's rAF, so draw may read marker
    draw(canvas, store) {
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.clientWidth * dpr, h = canvas.clientHeight * dpr;
      canvas.width = w; canvas.height = h;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = theme().canvas; ctx.fillRect(0, 0, w, h);

      const st = store.views.tube;
      const ep = engineParams(store);
      const Lmax = lRangeOf(store.basePreset)[1];
      const eig = eigH(hamiltonian(ep, store.E, store.rho, store.anti));
      const NPT = 240;
      const cums = [new Float64Array(NPT), new Float64Array(NPT), new Float64Array(NPT)];
      for (let i = 0; i < NPT; i++) {
        const L = (i / (NPT - 1)) * Lmax;
        const pe = prob(eig, 1, 0, L), pm = prob(eig, 1, 1, L), pt = prob(eig, 1, 2, L);
        cums[0][i] = pe; cums[1][i] = pe + pm; cums[2][i] = pe + pm + pt;
      }

      const P = plot2d(ctx, w, h, dpr, { x: [0, Lmax], y: [0, 1], xTitle: 'L [km]', yTitle: 'P(νμ→να)' });
      const xOf = (i) => P.X((i / (NPT - 1)) * Lmax);
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

      // white vertical marker line (allowed: markerDriven); tracks the 3D disk:
      // the marker fraction in L mode, the experiment baseline in E/δCP mode
      const mx = P.X(st.anim === 'L' ? st.marker * Lmax : Math.min(store.L, Lmax));
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

  companion2: {
    title: (store) => `Flavor fractions vs E at L = ${store.L} km`,
    markerDriven: true, // repainted per frame so the white line can track an E animation
    draw(canvas, store) {
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.clientWidth * dpr, h = canvas.clientHeight * dpr;
      canvas.width = w; canvas.height = h;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = theme().canvas; ctx.fillRect(0, 0, w, h);

      const ep = engineParams(store);
      const [E0, E1] = eRangeOf(store.basePreset);
      const NPT = 160; // one eigH per point
      const cums = [new Float64Array(NPT), new Float64Array(NPT), new Float64Array(NPT)];
      for (let i = 0; i < NPT; i++) {
        const E = E0 + (i / (NPT - 1)) * (E1 - E0);
        const eig = eigH(hamiltonian(ep, E, store.rho, store.anti));
        const pe = prob(eig, 1, 0, store.L), pm = prob(eig, 1, 1, store.L), pt = prob(eig, 1, 2, store.L);
        cums[0][i] = pe; cums[1][i] = pe + pm; cums[2][i] = pe + pm + pt;
      }

      const P = plot2d(ctx, w, h, dpr, { x: [E0, E1], y: [0, 1], xTitle: 'E [GeV]', yTitle: 'P(νμ→να)' });
      const xOf = (i) => P.X(E0 + (i / (NPT - 1)) * (E1 - E0));
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

      // white vertical marker line at the shared E (the tube/slice position)
      const mx = P.X(Math.max(E0, Math.min(E1, store.E)));
      ctx.strokeStyle = theme().hiCss;
      ctx.lineWidth = dpr;
      ctx.beginPath(); ctx.moveTo(mx, P.mt); ctx.lineTo(mx, P.mt + P.ph); ctx.stroke();
    },
  },
};
