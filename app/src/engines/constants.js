// Shared physics constants and defaults.
// Unit constants are FULL PRECISION, adopted from NuFast-LBL (Denton & Parke, arXiv:2405.02400).
// Do NOT use the legacy html/script.js matter constant (0.76e-4) — it is half the standard value
// (documented on github.com/czczc/nuosc/issues/8).

export const eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4; // (Delta m^2 [eV^2] * L [km] / E [GeV]) -> phase/4
export const KM_PER_EVSQ_OVER_GEV = 1e-6 / 1.97327e-7 / 2;      // Jacobi engine phase factor: H entries per km
export const YerhoE2a = 1.52588e-4;                              // a = 2*sqrt2*GF*Ne*E [eV^2 per (g/cc * GeV)]

export const DEG = Math.PI / 180;

// PDG 2023 defaults (matching the legacy app and the validated prototypes)
export const DEFAULTS = {
  th12: 33.82,   // degrees
  th13: 8.61,
  th23: 48.3,
  dm21: 7.39,    // 1e-5 eV^2
  dm31: 2.525,   // 1e-3 eV^2 (magnitude; sign from mass ordering)
  dcp: 217,      // degrees
  Ye: 0.5,
};

export const PRESETS = {
  DUNE: { L: 1300, rho: 2.85 },
  NOvA: { L: 810, rho: 2.84 },
  T2K: { L: 295, rho: 2.6 },
  JUNO: { L: 52.5, rho: 2.45 },
};

// Convert UI parameter state to engine-ready values.
export function engineParams(p) {
  return {
    s12sq: Math.sin(p.th12 * DEG) ** 2,
    s13sq: Math.sin(p.th13 * DEG) ** 2,
    s23sq: Math.sin(p.th23 * DEG) ** 2,
    th12: p.th12 * DEG,
    th13: p.th13 * DEG,
    th23: p.th23 * DEG,
    delta: p.dcp * DEG,
    dm21: p.dm21 * 1e-5,
    dm31: (p.normalOrdering ? 1 : -1) * p.dm31 * 1e-3,
    Ye: p.Ye ?? DEFAULTS.Ye,
  };
}
