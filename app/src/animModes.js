// Shared animation modes for the per-view "animate" marker row.
// Every view offers the same list: the four built-ins (L, E, δCP, ρ — L and E
// spans follow the active experiment) plus user-defined modes, each sweeping one
// shared parameter over a custom range. User modes are persisted in localStorage
// and referenced as 'u:<name>' in store.views[id].anim.
import { reactive } from 'vue';
import { store } from './store.js';
import { engineParams, lRangeOf, eRangeOf, lStepOf, eStepOf, eUnitOf, fmtE } from './engines/constants.js';

const KEY = 'nuglass-anim-modes';

// Parameters a mode can sweep: shared store keys. min/max/step may be
// (store) => value (preset-dependent L and E spans); ranges and steps match the
// shared sliders in ControlsCard. E is stored in GeV (reactor MeV is display-only).
export const ANIM_PARAMS = {
  L: {
    short: 'L', unit: 'km', axis: () => 'L [km]',
    min: 0, max: (s) => lRangeOf(s.basePreset)[1], step: (s) => lStepOf(s.basePreset),
  },
  E: {
    short: 'E', unit: 'GeV', axis: (s) => `E [${eUnitOf(s.basePreset).unit}]`,
    min: (s) => eRangeOf(s.basePreset)[0], max: (s) => eRangeOf(s.basePreset)[1], step: (s) => eStepOf(s.basePreset),
  },
  dcp: { short: 'δCP', unit: '°', axis: () => 'δCP [°]', min: 0, max: 360, step: 1 },
  rho: { short: 'ρ', unit: 'g/cm³', axis: () => 'ρ [g/cm³]', min: 0, max: 5, step: 0.05 },
  th12: { short: 'θ₁₂', unit: '°', axis: () => 'θ₁₂ [°]', min: 25, max: 45, step: 0.01 },
  th13: { short: 'θ₁₃', unit: '°', axis: () => 'θ₁₃ [°]', min: 7, max: 10, step: 0.01 },
  th23: { short: 'θ₂₃', unit: '°', axis: () => 'θ₂₃ [°]', min: 38, max: 54, step: 0.1 },
  dm21: { short: 'Δm²₂₁', unit: '10⁻⁵ eV²', axis: () => 'Δm²₂₁ [10⁻⁵ eV²]', min: 6, max: 9, step: 0.01 },
  dm31: { short: 'Δm²₃₁', unit: '10⁻³ eV²', axis: () => '|Δm²₃₁| [10⁻³ eV²]', min: 2.3, max: 2.7, step: 0.001 },
};

export const BUILTIN_MODES = ['L', 'E', 'dcp', 'rho'];

const lim = (v, s) => (typeof v === 'function' ? v(s) : v);

// user modes: name -> { param, min, max }
export const userModes = reactive({});

function persist() {
  try { localStorage.setItem(KEY, JSON.stringify(userModes)); } catch { /* private mode */ }
}

export function loadUserModes() {
  try {
    const obj = JSON.parse(localStorage.getItem(KEY) ?? '{}');
    for (const [name, d] of Object.entries(obj)) {
      if (ANIM_PARAMS[d?.param] && Number.isFinite(d.min) && Number.isFinite(d.max)) {
        userModes[name] = { param: d.param, min: d.min, max: d.max };
      }
    }
  } catch { /* corrupted entry -> start empty */ }
}

export function saveUserMode(name, def) {
  userModes[name] = { ...def };
  persist();
}

export function deleteUserMode(name) {
  delete userModes[name];
  persist();
  // any view still animating the deleted mode falls back to the default
  const gone = `u:${name}`;
  for (const vs of Object.values(store.views)) if (vs.anim === gone) vs.anim = 'L';
}

// Options for the animate select, identical in every view.
export function animOptions() {
  return [
    ...BUILTIN_MODES.map((m) => ({ value: m, label: ANIM_PARAMS[m].short })),
    ...Object.keys(userModes).map((n) => ({ value: `u:${n}`, label: n })),
  ];
}

// Mode value -> resolved sweep: the swept store key, its [lo, hi] span, slider
// step and labels. Unknown values (a deleted user mode) fall back to L.
export function modeDef(value, s = store) {
  let param = value, lo, hi;
  if (typeof value === 'string' && value.startsWith('u:')) {
    const d = userModes[value.slice(2)];
    if (d) { ({ param } = d); lo = d.min; hi = d.max; } else param = 'L';
  }
  let p = ANIM_PARAMS[param];
  if (!p) { param = 'L'; p = ANIM_PARAMS.L; }
  return {
    param,
    lo: lo ?? lim(p.min, s),
    hi: hi ?? lim(p.max, s),
    step: lim(p.step, s),
    short: p.short,
    axis: p.axis(s),
  };
}

// Tick helper for the slider-driven views: write the marker's swept value into
// the shared store key, rounded to that slider's step.
export function driveShared(s, value, frac) {
  const { param, lo, hi, step } = modeDef(value, s);
  s[param] = Math.round((lo + frac * (hi - lo)) / step) * step;
}

// 'name value' readout for the probe chips, e.g. 'θ₁₃ 8.62°', 'E 2.50 GeV'.
export function fmtSweep(value, v, s = store) {
  const { param, step, short } = modeDef(value, s);
  if (param === 'E') return `E ${fmtE(v, s.basePreset)}`;
  const prec = (String(step).split('.')[1] || '').length;
  const unit = ANIM_PARAMS[param].unit;
  return `${short} ${v.toFixed(prec)}${unit === '°' ? '°' : ` ${unit}`}`;
}

// Engine params at swept value xv: for oscillation parameters (δCP, θij, Δm²)
// the swept one is overridden; L, E, ρ are propagation inputs the caller passes
// to the hamiltonian itself, so the shared params are returned unchanged.
export function epAt(value, xv, s = store) {
  const { param } = modeDef(value, s);
  if (param === 'L' || param === 'E' || param === 'rho') return engineParams(s);
  const { th12, th13, th23, dm21, dm31, dcp, Ye, normalOrdering } = s;
  return engineParams({ th12, th13, th23, dm21, dm31, dcp, Ye, normalOrdering, [param]: xv });
}
