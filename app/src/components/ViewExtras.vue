<script setup>
import { store } from '../store.js';

const props = defineProps({ viewDef: { type: Object, required: true } });
const vs = store.views[props.viewDef.id];

function fmt(v, step) {
  const prec = (String(step ?? 1).split('.')[1] || '').length;
  return Number(v).toFixed(prec);
}

// range extras may give min/max as (store) => value, e.g. preset-dependent E spans
function lim(v) {
  return typeof v === 'function' ? v(store) : v;
}
</script>

<template>
  <div v-if="viewDef.extras.length || viewDef.note" class="card">
    <h3>{{ viewDef.label }}</h3>
    <template v-for="e in viewDef.extras" :key="e.key">
      <div v-if="e.type === 'select'" class="ctl-row">
        <label :for="'x-' + e.key">{{ e.label }}</label>
        <select :id="'x-' + e.key" v-model="vs[e.key]">
          <option v-for="o in e.options" :key="o.value" :value="o.value">{{ o.label }}</option>
        </select>
      </div>
      <div v-else-if="e.type === 'checkbox'" class="ctl-row">
        <label :for="'x-' + e.key">{{ e.label }}</label>
        <input :id="'x-' + e.key" v-model="vs[e.key]" type="checkbox" />
      </div>
      <div v-else-if="e.type === 'range'" class="ctl-row">
        <label :for="'x-' + e.key">{{ e.label }}</label>
        <input :id="'x-' + e.key" v-model.number="vs[e.key]" type="range" :min="lim(e.min)" :max="lim(e.max)" :step="e.step" />
        <span class="val">{{ fmt(vs[e.key], e.step) }}</span>
      </div>
    </template>
    <p v-if="viewDef.note" class="note">{{ viewDef.note }}</p>
  </div>
</template>

<style scoped>
select {
  background: var(--surface-3);
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 3px 6px;
  flex: 2;
  min-width: 0;
}
</style>
