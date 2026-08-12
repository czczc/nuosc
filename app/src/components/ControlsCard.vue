<script setup>
import { computed } from 'vue';
import { store, markCustom } from '../store.js';
import { PRESETS, lRangeOf } from '../engines/constants.js';

const beam = computed(() => PRESETS[store.basePreset] ?? null);

const shared = [
  { key: 'dcp', label: 'δCP [°]', min: 0, max: 360, step: 1 },
  { key: 'L', label: 'L [km]', min: 0, max: (s) => lRangeOf(s.basePreset)[1], step: 5, custom: true },
  { key: 'rho', label: 'ρ [g/cm³]', min: 0, max: 15, step: 0.05, custom: true },
];
const rare = [
  { key: 'th12', label: 'θ₁₂ [°]', min: 25, max: 45, step: 0.01 },
  { key: 'th13', label: 'θ₁₃ [°]', min: 7, max: 10, step: 0.01 },
  { key: 'th23', label: 'θ₂₃ [°]', min: 38, max: 54, step: 0.1 },
  { key: 'dm21', label: 'Δm²₂₁ [10⁻⁵]', min: 6, max: 9, step: 0.01 },
  { key: 'dm31', label: '|Δm²₃₁| [10⁻³]', min: 2.3, max: 2.7, step: 0.001 },
];
function fmt(v, step) {
  const prec = (String(step).split('.')[1] || '').length;
  return Number(v).toFixed(prec);
}

// min/max may be (store) => value, e.g. preset-dependent L span
function lim(v) {
  return typeof v === 'function' ? v(store) : v;
}
</script>

<template>
  <div class="card">
    <h3>Controls</h3>
    <div v-if="beam" class="note">
      {{ store.basePreset }} beam: {{ beam.Erange[0] }}–{{ beam.Erange[1] }} GeV · peak {{ beam.Epeak }} GeV
    </div>
    <div v-for="p in shared" :key="p.key" class="ctl-row">
      <label :for="p.key">{{ p.label }}</label>
      <input :id="p.key" v-model.number="store[p.key]" type="range" :min="lim(p.min)" :max="lim(p.max)" :step="p.step"
        @input="p.custom && markCustom()" />
      <span class="val">{{ fmt(store[p.key], p.step) }}</span>
    </div>
    <div class="ctl-row">
      <label>particle</label>
      <span class="seg">
        <button :class="{ on: !store.anti }" @click="store.anti = false">ν</button>
        <button :class="{ on: store.anti }" @click="store.anti = true">ν̅</button>
      </span>
      <label style="flex: none; margin-left: 8px;">ordering</label>
      <span class="seg">
        <button :class="{ on: store.normalOrdering }" @click="store.normalOrdering = true">NO</button>
        <button :class="{ on: !store.normalOrdering }" @click="store.normalOrdering = false">IO</button>
      </span>
    </div>
    <div class="ctl-row" style="justify-content: flex-end;">
      <button class="linkish" @click="store.showAllParams = !store.showAllParams">
        all parameters {{ store.showAllParams ? '▴' : '▾' }}
      </button>
    </div>
    <template v-if="store.showAllParams">
      <div v-for="p in rare" :key="p.key" class="ctl-row">
        <label :for="p.key">{{ p.label }}</label>
        <input :id="p.key" v-model.number="store[p.key]" type="range" :min="p.min" :max="p.max" :step="p.step" />
        <span class="val">{{ fmt(store[p.key], p.step) }}</span>
      </div>
    </template>
  </div>
</template>
