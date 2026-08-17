<script setup>
import { computed, reactive, ref } from 'vue';
import { store, applyPreset } from '../store.js';
import { router } from '../router.js';
import { DEFAULTS, PRESETS, CHANNELS, channelLabel } from '../engines/constants.js';
import { userExps, saveUserExp, deleteUserExp } from '../experiments.js';

// form starts from DUNE's setup and the NuFit 6.1 global fit; all channels checked
const blank = () => ({
  name: '',
  L: PRESETS.DUNE.L, rho: PRESETS.DUNE.rho,
  Emin: PRESETS.DUNE.Erange[0], Emax: PRESETS.DUNE.Erange[1], Epeak: PRESETS.DUNE.Epeak,
  channels: Object.keys(CHANNELS), anti: false,
  th12: DEFAULTS.th12, th13: DEFAULTS.th13, th23: DEFAULTS.th23,
  dm21: DEFAULTS.dm21, dm31: DEFAULTS.dm31, dcp: DEFAULTS.dcp, Ye: DEFAULTS.Ye,
});
const f = reactive(blank());
const err = ref('');
const savedFlash = ref('');

const expRows = [
  { key: 'L', label: 'baseline L [km]' },
  { key: 'rho', label: 'density ρ [g/cm³]' },
  { key: 'Emin', label: 'E min [GeV]' },
  { key: 'Emax', label: 'E max [GeV]' },
  { key: 'Epeak', label: 'E peak [GeV]' },
];
const oscRows = [
  { key: 'th12', label: 'θ₁₂ [°]' },
  { key: 'th13', label: 'θ₁₃ [°]' },
  { key: 'th23', label: 'θ₂₃ [°]' },
  { key: 'dm21', label: 'Δm²₂₁ [10⁻⁵ eV²]' },
  { key: 'dm31', label: '|Δm²₃₁| [10⁻³ eV²]' },
  { key: 'dcp', label: 'δCP [°]' },
  { key: 'Ye', label: 'electron fraction Yₑ' },
];
const NUM_KEYS = [...expRows, ...oscRows].map((r) => r.key);

const names = computed(() => Object.keys(userExps).sort());

function validate() {
  const name = f.name.trim();
  if (!name) return 'give the experiment a name';
  if (PRESETS[name] || name === 'custom') return `"${name}" is a built-in name; pick another`;
  for (const k of NUM_KEYS) if (!Number.isFinite(f[k])) return 'all fields must be numbers';
  if (f.L <= 0) return 'L must be positive';
  if (!(f.Emin > 0 && f.Emax > f.Emin)) return 'need 0 < E min < E max';
  if (f.Epeak < f.Emin || f.Epeak > f.Emax) return 'E peak must lie inside the E range';
  if (!f.channels.length) return 'pick at least one channel the experiment measures';
  return '';
}

function save() {
  err.value = validate();
  if (err.value) return;
  const name = f.name.trim();
  const { name: _, ...def } = f;
  def.channels = [...f.channels]; // detach from the form's reactive array
  saveUserExp(name, def);
  savedFlash.value = name;
  setTimeout(() => { savedFlash.value = ''; }, 1500);
}

function load(name) {
  applyPreset(name);
  router.push('/' + store.view);
}

function edit(name) {
  // blank() first: experiments saved before channels/anti existed keep the defaults
  Object.assign(f, blank(), userExps[name], { name });
  f.channels = [...f.channels]; // detach from the saved definition's array
  err.value = '';
}

// summary line for a saved experiment's channels
function chanSummary(d) {
  const chs = d.channels ?? Object.keys(CHANNELS);
  return chs.length === Object.keys(CHANNELS).length ? 'all channels' : chs.map(channelLabel).join(' ');
}

function del(name) {
  deleteUserExp(name);
  if (store.basePreset === name) applyPreset('DUNE'); // don't leave the app on a deleted preset
}

function resetForm() {
  Object.assign(f, blank());
  err.value = '';
}
</script>

<template>
  <div class="exp-page">
    <article class="doc">
      <h2>My experiments</h2>
      <p class="sub">
        Define your own setup (baseline, matter density, energy window, which oscillation channels it
        measures, and optionally your own oscillation parameters) and save it under a name. Saved
        experiments are kept in this browser's storage and can be loaded here any time; loading one
        makes it the active experiment everywhere (sliders, sweep ranges, the reset ↺ button); switching
        the channel rail to a channel it doesn't measure falls back to a built-in experiment. The form
        starts from DUNE's numbers and the NuFit 6.1 global fit.
      </p>

      <section>
        <h3>{{ userExps[f.name.trim()] ? `edit “${f.name.trim()}”` : 'new experiment' }}</h3>
        <div class="grid">
          <label for="exp-name">name</label>
          <input id="exp-name" v-model="f.name" type="text" placeholder="e.g. ESSnuSB" />
          <template v-for="r in expRows" :key="r.key">
            <label :for="'exp-' + r.key">{{ r.label }}</label>
            <input :id="'exp-' + r.key" v-model.number="f[r.key]" type="number" step="any" />
          </template>
        </div>
        <div class="chanrow">
          <label>channels measured</label>
          <span class="checks">
            <label v-for="(c, id) in CHANNELS" :key="id" class="check">
              <input v-model="f.channels" type="checkbox" :value="id" /> {{ channelLabel(id) }}
            </label>
            <label class="check">
              <input v-model="f.anti" type="checkbox" /> ν̄ source (reactor-like: loads with the particle toggle on ν̄)
            </label>
          </span>
        </div>
        <h4>oscillation parameters</h4>
        <div class="grid">
          <template v-for="r in oscRows" :key="r.key">
            <label :for="'exp-' + r.key">{{ r.label }}</label>
            <input :id="'exp-' + r.key" v-model.number="f[r.key]" type="number" step="any" />
          </template>
        </div>
        <p v-if="err" class="err">{{ err }}</p>
        <p class="btnrow">
          <button class="primary" @click="save">{{ savedFlash ? 'saved ✓' : 'save' }}</button>
          <button @click="resetForm">reset form to DUNE defaults</button>
        </p>
      </section>

      <section>
        <h3>saved experiments</h3>
        <p v-if="!names.length" class="sub">Nothing saved yet; fill in the form above and hit save.</p>
        <div v-for="n in names" :key="n" class="exprow">
          <div class="meta">
            <strong>{{ n }}</strong>
            <span class="sum">
              L {{ userExps[n].L }} km · ρ {{ userExps[n].rho }} g/cm³ ·
              E {{ userExps[n].Emin }}–{{ userExps[n].Emax }} GeV (peak {{ userExps[n].Epeak }}) ·
              δCP {{ userExps[n].dcp }}° ·
              {{ chanSummary(userExps[n]) }}<template v-if="userExps[n].anti"> · ν̄ source</template>
            </span>
          </div>
          <span class="acts">
            <button class="primary" @click="load(n)">load</button>
            <button @click="edit(n)">edit</button>
            <button class="danger" @click="del(n)">delete</button>
          </span>
        </div>
      </section>
    </article>
  </div>
</template>

<style scoped>
.exp-page {
  flex: 1;
  min-width: 0;
  overflow-y: auto;
  background: var(--bg);
}
.doc {
  max-width: 820px;
  margin: 0 auto;
  padding: 34px 28px 80px;
  font-size: 15px;
  line-height: 1.6;
}
h2 { font-size: 26px; margin: 0 0 6px; }
.sub { color: var(--muted); margin: 0 0 10px; }
section {
  margin-top: 26px;
  padding-top: 14px;
  border-top: 1px solid var(--border);
}
h3 {
  font-family: var(--font-mono);
  font-size: 17px;
  color: var(--accent);
  margin: 0 0 10px;
}
h4 {
  font-family: var(--font-mono);
  font-size: 13px;
  color: var(--muted);
  margin: 18px 0 8px;
}
.grid {
  display: grid;
  grid-template-columns: max-content 1fr max-content 1fr;
  gap: 8px 12px;
  align-items: center;
}
.grid label {
  font-family: var(--font-mono);
  font-size: 13px;
  color: var(--muted);
  text-align: right;
}
.grid input {
  min-width: 0;
  background: var(--surface-3);
  color: var(--text);
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 5px 8px;
  font-family: var(--font-mono);
  font-size: 13px;
}
.grid input:focus { outline: none; border-color: var(--accent); }
.chanrow {
  display: flex;
  gap: 12px;
  align-items: baseline;
  margin-top: 10px;
}
.chanrow > label {
  font-family: var(--font-mono);
  font-size: 13px;
  color: var(--muted);
  flex: none;
}
.checks {
  display: flex;
  flex-wrap: wrap;
  gap: 6px 16px;
}
.check {
  display: inline-flex;
  align-items: center;
  gap: 5px;
  font-family: var(--font-mono);
  font-size: 13px;
  cursor: pointer;
}
.check input { accent-color: var(--accent); margin: 0; }
.err { color: #e05545; font-family: var(--font-mono); font-size: 13px; }
.btnrow { display: flex; gap: 10px; margin-top: 14px; }
button {
  background: var(--surface-3);
  color: var(--text);
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 5px 12px;
  font-family: var(--font-mono);
  font-size: 13px;
  cursor: pointer;
}
button:hover { border-color: var(--accent); }
button.primary { background: var(--accent); border-color: var(--accent); color: var(--accent-ink); }
button.danger:hover { border-color: #e05545; color: #e05545; }
.exprow {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 12px;
  padding: 9px 0;
  border-bottom: 1px solid var(--border);
}
.meta { display: flex; flex-direction: column; min-width: 0; }
.meta strong { font-family: var(--font-mono); }
.sum { color: var(--muted); font-size: 13px; }
.acts { display: flex; gap: 8px; flex: none; }
@media (max-width: 700px) {
  .grid { grid-template-columns: max-content 1fr; }
  .exprow { flex-direction: column; align-items: flex-start; }
}
</style>
