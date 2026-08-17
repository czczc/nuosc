// Shared physics constants and defaults.
// Unit constants are FULL PRECISION, adopted from NuFast-LBL (Denton & Parke, arXiv:2405.02400).
// Do NOT use the legacy/html/script.js matter constant (0.76e-4) — it is half the standard value
// (documented on github.com/czczc/nuglass/issues/8).

export const eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4; // (Delta m^2 [eV^2] * L [km] / E [GeV]) -> phase/4
export const KM_PER_EVSQ_OVER_GEV = 1e-6 / 1.97327e-7 / 2;      // Jacobi engine phase factor: H entries per km
export const YerhoE2a = 1.52588e-4;                              // a = 2*sqrt2*GF*Ne*E [eV^2 per (g/cc * GeV)]

export const DEG = Math.PI / 180;

// NuFit 6.1 (2025) global-fit best fit, normal ordering, IC24 with SK atmospheric
// data (www.nu-fit.org, arXiv:2410.05380). IO values differ slightly; the app
// approximates IO by flipping the sign of dm31.
export const DEFAULTS = {
  th12: 33.76,   // degrees
  th13: 8.62,
  th23: 43.29,
  dm21: 7.537,   // 1e-5 eV^2
  dm31: 2.511,   // 1e-3 eV^2 (magnitude; sign from mass ordering)
  dcp: 212,      // degrees
  Ye: 0.5,
};

// Oscillation channels selectable from the header rail. a/b = initial/final
// flavor indices into the engines' 3x3 probability matrix (0 = e, 1 = mu, 2 = tau).
// Labels use plain nu; antineutrinos stay on the shared particle toggle.
export const CHANNELS = {
  mue: { a: 1, b: 0, nu: ['νμ', 'νe'] },
  mumu: { a: 1, b: 1, nu: ['νμ', 'νμ'] },
  mutau: { a: 1, b: 2, nu: ['νμ', 'ντ'] },
  ee: { a: 0, b: 0, nu: ['νe', 'νe'] },
};
export function channelLabel(ch) {
  return CHANNELS[ch].nu.join('→');
}
// 'P(νμ→νe)', with bars (and P̄) when antineutrinos are selected
export function pLabel(ch, anti = false) {
  const [a, b] = CHANNELS[ch].nu;
  const bar = (s) => s.replace('ν', 'ν̄');
  return anti ? `P̄(${bar(a)}→${bar(b)})` : `P(${a}→${b})`;
}

// Experiment presets. Erange = flux span [GeV], Epeak = flux peak [GeV] (reactor
// windows are MeV-scale; see eUnitOf). channels = which oscillation channels the
// experiment measures; anti = the source is antineutrinos (checks the shared
// particle toggle on selection).
export const PRESETS = {
  DUNE: { L: 1300, rho: 2.85, Erange: [0.5, 6], Epeak: 2.5, channels: ['mue', 'mumu', 'mutau'] },
  NOvA: { L: 810, rho: 2.84, Erange: [1, 3], Epeak: 2.0, channels: ['mue', 'mumu'] },
  T2K: { L: 295, rho: 2.6, Erange: [0.2, 1.5], Epeak: 0.6, channels: ['mue', 'mumu'] },
  JUNO: { L: 52.5, rho: 2.6, Erange: [0.0018, 0.009], Epeak: 0.004, channels: ['ee'], anti: true },
  KamLAND: { L: 180, rho: 2.6, Erange: [0.0018, 0.009], Epeak: 0.004, channels: ['ee'], anti: true },
  'Daya Bay': { L: 1.66, rho: 2.6, Erange: [0.0018, 0.009], Epeak: 0.004, channels: ['ee'], anti: true },
};

// User-defined experiments (src/experiments.js) are registered here so the range
// helpers and applyPreset resolve them exactly like built-ins. Entries share the
// PRESETS shape plus an optional `params` object of oscillation-parameter overrides.
export const USER_PRESETS = {};
export function presetOf(name) {
  return PRESETS[name] ?? USER_PRESETS[name];
}

// Active baseline span [km]: 0 to twice the experiment's L, or the full default span.
export function lRangeOf(preset) {
  const p = presetOf(preset);
  return p ? [0, 2 * p.L] : [0, 5000];
}

// Active energy span [GeV]: the selected experiment's beam window, or the full default span.
export const E_RANGE_DEFAULT = [0.2, 6];
export function eRangeOf(preset) {
  return presetOf(preset)?.Erange ?? E_RANGE_DEFAULT;
}

// Channels an experiment measures; user experiments without a declaration
// measure everything (github.com/czczc/nuglass/issues/15 adds the declaration).
export function channelsOf(preset) {
  return presetOf(preset)?.channels ?? Object.keys(CHANNELS);
}

// Display unit for the experiment's energies: reactor windows are MeV-scale.
// E is stored in GeV everywhere; only labels and slider steps switch.
export function eUnitOf(preset) {
  return eRangeOf(preset)[1] < 0.05 ? { unit: 'MeV', scale: 1000 } : { unit: 'GeV', scale: 1 };
}
export function eStepOf(preset) {
  return 0.01 / eUnitOf(preset).scale; // 0.01 in display units
}
export function fmtE(E, preset) {
  const { unit, scale } = eUnitOf(preset);
  return `${(E * scale).toFixed(2)} ${unit}`;
}
// Slider/animation step for L [km], scaled to the experiment's span (Daya Bay's
// 1.66 km baseline needs finer steps than DUNE's 1300 km).
export function lStepOf(preset) {
  const Lmax = lRangeOf(preset)[1];
  return Lmax >= 1000 ? 5 : Lmax >= 100 ? 0.5 : 0.01;
}

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
