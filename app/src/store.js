import { reactive } from 'vue';
import { DEFAULTS, PRESETS, presetOf, channelsOf } from './engines/constants.js';

const storedTheme = typeof localStorage !== 'undefined' ? localStorage.getItem('nuglass-theme') : null;

export const store = reactive({
  theme: storedTheme === 'light' ? 'light' : 'dark',
  view: 'oscillogram',
  faq: null, // null = closed; 'engines' or a view id = FAQ page open, scrolled to that section
  exps: false, // "my experiments" page open
  preset: 'DUNE',
  basePreset: 'DUNE', // last experiment clicked; keeps L spans stable while tweaking into "custom"
  channel: 'mue', // oscillation channel (header rail), a CHANNELS key

  // shared physics controls
  E: PRESETS.DUNE.Epeak,
  L: PRESETS.DUNE.L,
  rho: PRESETS.DUNE.rho,
  dcp: DEFAULTS.dcp,
  anti: false,
  normalOrdering: true,

  // "all parameters" expansion (NuFit 6.1 defaults, legacy app ranges)
  th12: DEFAULTS.th12,
  th13: DEFAULTS.th13,
  th23: DEFAULTS.th23,
  dm21: DEFAULTS.dm21, // 1e-5 eV^2
  dm31: DEFAULTS.dm31, // 1e-3 eV^2, magnitude
  Ye: DEFAULTS.Ye,
  showAllParams: false,

  // per-view state (marker is a 0..1 fraction of the swept range)
  views: {
    oscillogram: { axis2: 'L', anim: 'dcp', palette: 'rainbow', play: true, marker: 0.5 },
    tube: { mode: 'tube', anim: 'L', play: true, marker: 0.25 },
    sphere: { sweep: 'L', pole: 'both', play: true, marker: 0 },
    phasors: { xaxis: 'L', play: true, marker: 0 },
    biprob: { showNO: true, showIO: true, showSurf: false, anim: 'E', play: true, marker: 0.5 },
    loe: {
      cMue: true, cMumu: true, cMutau: false, cEe: true, surf: true, exps: true, avg: false,
      anim: 'dcp', play: true, marker: 0.5,
    },
  },
});

export function applyPreset(name) {
  const p = presetOf(name); // built-in or user-defined experiment
  if (!p) return;
  store.preset = name;
  store.basePreset = name;
  store.L = p.L;
  store.rho = p.rho;
  // snap the shared energy to the experiment's beam peak
  store.E = p.Epeak;
  // oscillation parameters: the experiment's own overrides (user experiments)
  // or the global-fit defaults (built-ins)
  const o = p.params ?? DEFAULTS;
  store.th12 = o.th12;
  store.th13 = o.th13;
  store.th23 = o.th23;
  store.dm21 = o.dm21;
  store.dm31 = o.dm31;
  store.dcp = o.dcp;
  store.Ye = o.Ye ?? DEFAULTS.Ye;
  // the particle toggle follows the source: reactors are antineutrinos, beams reset to ν
  store.anti = p.anti ?? false;
  // keep the channel if this experiment measures it, else its first channel
  const chs = channelsOf(name);
  if (!chs.includes(store.channel)) store.channel = chs[0];
}

// Header rail: switch channel; if the active experiment doesn't measure it,
// fall back to the first built-in that does.
export function setChannel(ch) {
  store.channel = ch;
  if (!channelsOf(store.basePreset).includes(ch)) {
    const fallback = Object.keys(PRESETS).find((n) => PRESETS[n].channels.includes(ch));
    if (fallback) applyPreset(fallback);
  }
}

export function setTheme(t) {
  store.theme = t;
  document.documentElement.dataset.theme = t;
  try { localStorage.setItem('nuglass-theme', t); } catch { /* private mode */ }
}

// Any manual change to L or rho makes the preset "custom".
export function markCustom() {
  store.preset = 'custom';
}
