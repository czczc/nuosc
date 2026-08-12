import { reactive } from 'vue';
import { DEFAULTS, PRESETS } from './engines/constants.js';

const storedTheme = typeof localStorage !== 'undefined' ? localStorage.getItem('nuosc-theme') : null;

export const store = reactive({
  theme: storedTheme === 'light' ? 'light' : 'dark',
  view: 'oscillogram',
  faq: null, // null = closed; 'engines' or a view id = FAQ page open, scrolled to that section
  preset: 'DUNE',
  basePreset: 'DUNE', // last experiment clicked; keeps L spans stable while tweaking into "custom"

  // shared physics controls
  L: PRESETS.DUNE.L,
  rho: PRESETS.DUNE.rho,
  dcp: DEFAULTS.dcp,
  anti: false,
  normalOrdering: true,

  // "all parameters" expansion (PDG 2023 defaults, legacy app ranges)
  th12: DEFAULTS.th12,
  th13: DEFAULTS.th13,
  th23: DEFAULTS.th23,
  dm21: DEFAULTS.dm21, // 1e-5 eV^2
  dm31: DEFAULTS.dm31, // 1e-3 eV^2, magnitude
  Ye: DEFAULTS.Ye,
  showAllParams: false,

  // display
  ortho: true,

  // per-view state (marker is a 0..1 fraction of the swept range)
  views: {
    oscillogram: { axis2: 'L' },
    tube: { mode: 'tube', play: true, marker: 0.25, Lmax: 2 * PRESETS.DUNE.L, E: 2.5 },
    sphere: { sweep: 'L', play: true, marker: 0, Lmax: 2 * PRESETS.DUNE.L, E: 2.5 },
    phasors: { xaxis: 'L', play: true, marker: 0, Lmax: 2 * PRESETS.DUNE.L, E: 2.5 },
    biprob: { Eslice: 2.5, showNO: true, showIO: true },
  },
});

export function applyPreset(name) {
  const p = PRESETS[name];
  if (!p) return;
  store.preset = name;
  store.basePreset = name;
  store.L = p.L;
  store.rho = p.rho;
  // snap every per-view energy to the experiment's beam peak
  store.views.tube.E = p.Epeak;
  store.views.sphere.E = p.Epeak;
  store.views.phasors.E = p.Epeak;
  store.views.biprob.Eslice = p.Epeak;
  // snap the L sweeps to the experiment's span (0 - 2L)
  store.views.tube.Lmax = 2 * p.L;
  store.views.sphere.Lmax = 2 * p.L;
  store.views.phasors.Lmax = 2 * p.L;
}

export function setTheme(t) {
  store.theme = t;
  document.documentElement.dataset.theme = t;
  try { localStorage.setItem('nuosc-theme', t); } catch { /* private mode */ }
}

// Any manual change to L or rho makes the preset "custom".
export function markCustom() {
  store.preset = 'custom';
}
