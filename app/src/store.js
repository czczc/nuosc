import { reactive } from 'vue';
import { DEFAULTS, PRESETS } from './engines/constants.js';

export const store = reactive({
  view: 'oscillogram',
  preset: 'DUNE',

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
    tube: { mode: 'tube', play: true, marker: 0.25, Lmax: 5000, E: 2.5 },
    sphere: { sweep: 'L', play: true, marker: 0, Lmax: 10000, E: 2.5 },
    phasors: { xaxis: 'L', play: true, marker: 0, Lmax: 5000, E: 2.5 },
    biprob: { Eslice: 2.5, showNO: true, showIO: true },
  },
});

export function applyPreset(name) {
  const p = PRESETS[name];
  if (!p) return;
  store.preset = name;
  store.L = p.L;
  store.rho = p.rho;
}

// Any manual change to L or rho makes the preset "custom".
export function markCustom() {
  store.preset = 'custom';
}
