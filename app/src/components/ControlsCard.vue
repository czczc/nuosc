<script setup>
import { computed } from 'vue';
import { store, markCustom } from '../store.js';
import { presetOf, lRangeOf, eRangeOf, eUnitOf, eStepOf, lStepOf } from '../engines/constants.js';

const beam = computed(() => presetOf(store.basePreset) ?? null);
const eu = computed(() => eUnitOf(store.basePreset)); // GeV or MeV display

const shared = [
  { key: 'dcp', label: 'δCP [°]', min: 0, max: 360, step: 1 },
  {
    key: 'E', label: (s) => `E [${eUnitOf(s.basePreset).unit}]`,
    min: (s) => eRangeOf(s.basePreset)[0], max: (s) => eRangeOf(s.basePreset)[1],
    step: (s) => eStepOf(s.basePreset),
    // E is stored in GeV; reactor windows display in MeV
    disp: (v, s) => fmt(v * eUnitOf(s.basePreset).scale, 0.01),
  },
  { key: 'L', label: 'L [km]', min: 0, max: (s) => lRangeOf(s.basePreset)[1], step: (s) => lStepOf(s.basePreset), custom: true },
  { key: 'rho', label: 'ρ [g/cm³]', min: 0, max: 5, step: 0.05, custom: true },
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

// min/max/step/label may be (store) => value, e.g. preset-dependent L span
function lim(v) {
  return typeof v === 'function' ? v(store) : v;
}
// energies in the experiment's display unit (avoids float dust like 1.7999999...)
function eDisp(v) {
  return +(v * eu.value.scale).toFixed(4);
}
</script>

<template>
  <div class="card">
    <h3>Controls</h3>
    <div v-if="beam" class="note">
      {{ store.basePreset }} flux: {{ eDisp(beam.Erange[0]) }}–{{ eDisp(beam.Erange[1]) }} {{ eu.unit }} ·
      peak {{ eDisp(beam.Epeak) }} {{ eu.unit }}
    </div>
    <div v-for="p in shared" :key="p.key" class="ctl-row">
      <label :for="p.key">{{ lim(p.label) }}</label>
      <input :id="p.key" v-model.number="store[p.key]" type="range" :min="lim(p.min)" :max="lim(p.max)"
        :step="lim(p.step)" @input="p.custom && markCustom()" />
      <span class="val">{{ p.disp ? p.disp(store[p.key], store) : fmt(store[p.key], lim(p.step)) }}</span>
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
