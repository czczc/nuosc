<script setup>
import { store, applyPreset } from '../store.js';
import { router } from '../router.js';

const props = defineProps({ viewDef: { type: Object, required: true } });
const vs = store.views[props.viewDef.id];

function fmt(v, step) {
  const prec = (String(step ?? 1).split('.')[1] || '').length;
  return Number(v).toFixed(prec);
}

// range extras may give min/max as (store) => value, e.g. preset-dependent E spans;
// select options may likewise give label as (store) => string
function lim(v) {
  return typeof v === 'function' ? v(store) : v;
}
</script>

<template>
  <div v-if="viewDef.extras.length || viewDef.note" class="card">
    <h3>
      {{ viewDef.label }}
      <button class="faqlink" title="About this view (FAQ)" @click="router.push('/faq/' + viewDef.id)">?</button>
    </h3>
    <template v-for="e in viewDef.extras" :key="e.key">
      <div v-if="e.type === 'select'" class="ctl-row">
        <label :for="'x-' + e.key">{{ e.label }}</label>
        <select :id="'x-' + e.key" v-model="vs[e.key]">
          <option v-for="o in e.options" :key="o.value" :value="o.value">{{ lim(o.label) }}</option>
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
      <div v-else-if="e.type === 'marker'" class="ctl-row">
        <label :for="'x-' + e.key">{{ e.label }}</label>
        <select v-if="e.select" v-model="vs[e.select.key]" class="mini" :aria-label="e.label + ' variable'">
          <option v-for="o in e.select.options" :key="o.value" :value="o.value">{{ o.label }}</option>
        </select>
        <button class="playbtn" :title="vs.play ? 'pause' : 'play'" @click="vs.play = !vs.play">
          {{ vs.play ? '❚❚' : '▶' }}
        </button>
        <button class="playbtn reset" title="reset to experiment defaults" @click="applyPreset(store.basePreset)">
          ↺
        </button>
        <input :id="'x-' + e.key" v-model.number="vs[e.key]" type="range" min="0" max="1" :step="e.step" />
      </div>
    </template>
    <p v-if="viewDef.note" class="note">{{ viewDef.note }}</p>
  </div>
</template>

<style scoped>
h3 { display: flex; align-items: center; justify-content: space-between; }
.faqlink {
  width: 16px; height: 16px;
  border: 1px solid var(--border); border-radius: 50%;
  background: none; color: var(--muted);
  font-size: 10px; line-height: 1;
  padding: 0; cursor: pointer;
}
.faqlink:hover { color: var(--accent); border-color: var(--accent); }
.playbtn {
  flex: 0 0 24px; height: 22px;
  border: 1px solid var(--border); border-radius: 6px;
  background: var(--surface-3); color: var(--accent);
  font-size: 9px; line-height: 1;
  padding: 0; cursor: pointer;
}
.playbtn:hover { border-color: var(--accent); }
.playbtn.reset { font-size: 13px; }
select.mini { flex: 0 0 62px; }
select {
  background: var(--surface-3);
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 3px 6px;
  flex: 2;
  min-width: 0;
}
</style>
