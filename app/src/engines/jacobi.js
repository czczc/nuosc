// Exact constant-density 3-flavor amplitude engine (complex-Hermitian Jacobi diagonalization).
// Provides complex amplitudes and matter-eigenstate decompositions, which NuFast (probabilities
// only) cannot — used by the state-sphere and phasor views.
// H (flavor basis, per km) = KM_PER_EVSQ_OVER_GEV/E * U diag(0,dm21,dm31) U^dag + diag(+-Vkm,0,0).
// Antineutrino: delta -> -delta, V -> -V. Flavor indices: 0=e, 1=mu, 2=tau.
// Cross-validated against the NuFast port to 5.8e-11 (tests/engines.test.js).
import { KM_PER_EVSQ_OVER_GEV, YerhoE2a } from './constants.js';

export const C = {
  add: (a, b) => [a[0] + b[0], a[1] + b[1]],
  sub: (a, b) => [a[0] - b[0], a[1] - b[1]],
  mul: (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]],
  conj: (a) => [a[0], -a[1]],
  expi: (t) => [Math.cos(t), Math.sin(t)],
  abs2: (a) => a[0] * a[0] + a[1] * a[1],
};

export function pmns(th12, th13, th23, dcp) {
  const s12 = Math.sin(th12), c12 = Math.cos(th12);
  const s13 = Math.sin(th13), c13 = Math.cos(th13);
  const s23 = Math.sin(th23), c23 = Math.cos(th23);
  const eid = C.expi(dcp), emid = C.expi(-dcp);
  const s13eid = [s13 * eid[0], s13 * eid[1]];
  return [
    [[c12 * c13, 0], [s12 * c13, 0], [s13 * emid[0], s13 * emid[1]]],
    [C.sub([-s12 * c23, 0], C.mul([c12 * s23, 0], s13eid)),
     C.sub([c12 * c23, 0], C.mul([s12 * s23, 0], s13eid)),
     [s23 * c13, 0]],
    [C.sub([s12 * s23, 0], C.mul([c12 * c23, 0], s13eid)),
     C.sub([-c12 * s23, 0], C.mul([s12 * c23, 0], s13eid)),
     [c23 * c13, 0]],
  ];
}

// ep: output of engineParams(); E [GeV], rho [g/cm^3].
export function hamiltonian(ep, E, rho, anti = false) {
  const d = anti ? -ep.delta : ep.delta;
  const U = pmns(ep.th12, ep.th13, ep.th23, d);
  const k = KM_PER_EVSQ_OVER_GEV / E;
  const m = [0, k * ep.dm21, k * ep.dm31];
  const H = [[[0, 0], [0, 0], [0, 0]], [[0, 0], [0, 0], [0, 0]], [[0, 0], [0, 0], [0, 0]]];
  for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) {
    let s = [0, 0];
    for (let i = 0; i < 3; i++) s = C.add(s, C.mul([m[i], 0], C.mul(U[a][i], C.conj(U[b][i]))));
    H[a][b] = s;
  }
  const V = (anti ? -1 : 1) * YerhoE2a * KM_PER_EVSQ_OVER_GEV * ep.Ye * rho;
  H[0][0] = [H[0][0][0] + V, 0];
  return H;
}

// Jacobi eigensolver for a complex Hermitian 3x3. Returns {lam: [3], W: columns = eigenvectors}.
export function eigH(Hin) {
  const A = Hin.map((r) => r.map((z) => [z[0], z[1]]));
  const W = [[[1, 0], [0, 0], [0, 0]], [[0, 0], [1, 0], [0, 0]], [[0, 0], [0, 0], [1, 0]]];
  for (let sweep = 0; sweep < 20; sweep++) {
    let off = 0;
    for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) off += C.abs2(A[p][q]);
    if (off < 1e-26) break;
    for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
      const h = A[p][q], r = Math.hypot(h[0], h[1]);
      if (r < 1e-18) continue;
      const phi = Math.atan2(h[1], h[0]);
      const th = 0.5 * Math.atan2(2 * r, A[q][q][0] - A[p][p][0]);
      const c = Math.cos(th), s = Math.sin(th);
      const spq = [s * Math.cos(phi), s * Math.sin(phi)];
      for (const M of [A, W]) for (let i = 0; i < 3; i++) {
        const mp = M[i][p], mq = M[i][q];
        M[i][p] = C.sub([c * mp[0], c * mp[1]], C.mul(mq, C.conj(spq)));
        M[i][q] = C.add(C.mul(mp, spq), [c * mq[0], c * mq[1]]);
      }
      for (let j = 0; j < 3; j++) {
        const ap = A[p][j], aq = A[q][j];
        A[p][j] = C.sub([c * ap[0], c * ap[1]], C.mul(spq, aq));
        A[q][j] = C.add(C.mul(C.conj(spq), ap), [c * aq[0], c * aq[1]]);
      }
    }
  }
  return { lam: [A[0][0][0], A[1][1][0], A[2][2][0]], W };
}

// Amplitude alpha -> beta after L km.
export function amp(eig, alpha, beta, L) {
  let s = [0, 0];
  for (let i = 0; i < 3; i++)
    s = C.add(s, C.mul(eig.W[beta][i], C.mul(C.expi(-eig.lam[i] * L), C.conj(eig.W[alpha][i]))));
  return s;
}

export function prob(eig, alpha, beta, L) {
  return C.abs2(amp(eig, alpha, beta, L));
}

// Phasor decomposition of the alpha -> beta amplitude: A(L) = sum_i c_i e^{-i lam_i L}.
export function phasorTerms(eig, alpha, beta) {
  return [0, 1, 2].map((i) => ({
    c: C.mul(eig.W[beta][i], C.conj(eig.W[alpha][i])),
    lam: eig.lam[i],
  }));
}
