// NuFast-LBL, JavaScript port of NuFast_LBL.cpp (Peter Denton & Stephen Parke, arXiv:2405.02400, MIT).
// Computes all nine oscillation probabilities in constant-density matter with real arithmetic only.
// Cross-validated against the complex-Jacobi engine to 5.8e-11 (N_Newton=2) / 1.4e-7 (N_Newton=1);
// see tests/engines.test.js and github.com/czczc/nuglass/issues/8.
import { eVsqkm_to_GeV_over4, YerhoE2a } from './constants.js';

// E > 0 for neutrinos, E < 0 for antineutrinos (upstream convention).
// Returns P[alpha][beta] with flavor indices 0=e, 1=mu, 2=tau.
export function probabilityMatter(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rhoYe, N_Newton = 1) {
  const c13sq = 1 - s13sq;
  const Ue2sq_tmp = c13sq * s12sq;
  const Ue3sq_tmp = s13sq;
  const Um3sq_tmp = c13sq * s23sq;
  let Ut2sq = s13sq * s12sq * s23sq;
  let Um2sq_tmp = (1 - s12sq) * (1 - s23sq);
  const Jrr = Math.sqrt(Um2sq_tmp * Ut2sq);
  const sind = Math.sin(delta), cosd = Math.cos(delta);
  Um2sq_tmp += Ut2sq - 2 * Jrr * cosd;
  const Jmatter_tmp = 8 * Jrr * c13sq * sind;
  const Dmsqee = Dmsq31 - s12sq * Dmsq21;
  const A_tmp = Dmsq21 + Dmsq31;
  const See = A_tmp - Dmsq21 * Ue2sq_tmp - Dmsq31 * Ue3sq_tmp;
  const Tmm_tmp = Dmsq21 * Dmsq31;
  const Tee = Tmm_tmp * (1 - Ue3sq_tmp - Ue2sq_tmp);

  const Amatter = rhoYe * E * YerhoE2a;
  const C = Amatter * Tee;
  const A = A_tmp + Amatter;
  const B = Tmm_tmp + Amatter * See;

  let lambda3;
  if (N_Newton < 0) {
    const rootAsqB = Math.sqrt(A * A - 3 * B);
    let ss0 = Math.acos((A * A * A - 4.5 * A * B + 13.5 * C) / (rootAsqB * rootAsqB * rootAsqB));
    if (Dmsq31 < 0) ss0 += 2 * Math.PI;
    lambda3 = (A + 2 * rootAsqB * Math.cos(ss0 / 3)) / 3;
  } else {
    const xmat = Amatter / Dmsqee;
    const t = 1 - xmat;
    lambda3 = Dmsq31 + 0.5 * Dmsqee * (xmat - 1 + Math.sqrt(t * t + 4 * s13sq * xmat));
    for (let j = 0; j < N_Newton; j++)
      lambda3 = (lambda3 * lambda3 * (lambda3 + lambda3 - A) + C) / (lambda3 * (2 * (lambda3 - A) + lambda3) + B);
  }

  const tmp = A - lambda3;
  const Dlambda21 = Math.sqrt(tmp * tmp - 4 * C / lambda3);
  const lambda2 = 0.5 * (A - lambda3 + Dlambda21);
  const Dlambda32 = lambda3 - lambda2;
  const Dlambda31 = Dlambda32 + Dlambda21;

  const PiDlambdaInv = 1 / (Dlambda31 * Dlambda32 * Dlambda21);
  const Xp3 = PiDlambdaInv * Dlambda21;
  const Xp2 = -PiDlambdaInv * Dlambda31;

  const Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3;
  const Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2;

  const Smm = A - Dmsq21 * Um2sq_tmp - Dmsq31 * Um3sq_tmp;
  const Tmm = Tmm_tmp * (1 - Um3sq_tmp - Um2sq_tmp) + Amatter * (See + Smm - A);

  const Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3;
  const Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2;

  const Jmatter = Jmatter_tmp * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21) * PiDlambdaInv;

  const Ue1sq = 1 - Ue3sq - Ue2sq;
  const Um1sq = 1 - Um3sq - Um2sq;
  const Ut3sq = 1 - Um3sq - Ue3sq;
  Ut2sq = 1 - Um2sq - Ue2sq;
  const Ut1sq = 1 - Um1sq - Ue1sq;

  const Lover4E = eVsqkm_to_GeV_over4 * L / E;
  const D21 = Dlambda21 * Lover4E;
  const D32 = Dlambda32 * Lover4E;
  const sinD21 = Math.sin(D21);
  const sinD31 = Math.sin(D32 + D21);
  const sinD32 = Math.sin(D32);
  const triple_sin = sinD21 * sinD31 * sinD32;
  const sinsqD21_2 = 2 * sinD21 * sinD21;
  const sinsqD31_2 = 2 * sinD31 * sinD31;
  const sinsqD32_2 = 2 * sinD32 * sinD32;

  const Pme_TC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2
               + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2
               + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
  const Pme_TV = -Jmatter * triple_sin;
  const Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2 + Um3sq * Um1sq * sinsqD31_2 + Um3sq * Um2sq * sinsqD32_2);
  const Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2 + Ue3sq * Ue1sq * sinsqD31_2 + Ue3sq * Ue2sq * sinsqD32_2);

  const P = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  P[0][0] = Pee;
  P[0][1] = Pme_TC - Pme_TV;
  P[0][2] = 1 - Pee - P[0][1];
  P[1][0] = Pme_TC + Pme_TV;
  P[1][1] = Pmm;
  P[1][2] = 1 - P[1][0] - Pmm;
  P[2][0] = 1 - Pee - P[1][0];
  P[2][1] = 1 - P[0][1] - Pmm;
  P[2][2] = 1 - P[0][2] - P[1][2];
  return P;
}

// Convenience wrapper over engineParams()-style objects.
// ep: output of engineParams(); anti flips via negative E.
export function probs(ep, E, L, rho, anti = false) {
  return probabilityMatter(ep.s12sq, ep.s13sq, ep.s23sq, ep.delta, ep.dm21, ep.dm31, L, anti ? -E : E, rho * ep.Ye, 1);
}
