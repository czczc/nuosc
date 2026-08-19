<script setup>
import { computed, reactive, ref } from 'vue';
import { store, applyPreset } from '../store.js';
import { router } from '../router.js';
import { ANIM_PARAMS, animOptions, saveUserMode, deleteUserMode } from '../animModes.js';

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

// ---- animate select: shared modes (4 built-ins + user modes) in every view ----
const animOpts = computed(() => animOptions());
const creating = ref(false);
const draft = reactive({ name: '', param: 'th13', min: 7, max: 10 });

function setDraftParam(param) {
  draft.param = param;
  draft.min = lim(ANIM_PARAMS[param].min);
  draft.max = lim(ANIM_PARAMS[param].max);
}

function onAnimChange(key, ev) {
  const v = ev.target.value;
  if (v === '__new__') {
    ev.target.value = vs[key]; // keep the current mode; the form takes over
    setDraftParam(draft.param);
    draft.name = '';
    creating.value = true;
  } else {
    vs[key] = v;
    creating.value = false;
  }
}

function saveDraft(key) {
  const min = +draft.min, max = +draft.max;
  if (!Number.isFinite(min) || !Number.isFinite(max) || min >= max) return;
  const name = draft.name.trim() || `${ANIM_PARAMS[draft.param].short} ${min}–${max}`;
  saveUserMode(name, { param: draft.param, min, max });
  vs[key] = `u:${name}`;
  creating.value = false;
}

// delete the selected user mode (views animating it fall back to L)
function deleteCurrent(key) {
  deleteUserMode(vs[key].slice(2));
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
      <div v-else-if="e.type === 'checkbox'" class="ctl-row" :title="e.title || null">
        <label :for="'x-' + e.key">{{ e.label }}</label>
        <!-- lock: (store) => bool — checkbox forced on and disabled (e.g. the header channel) -->
        <input :id="'x-' + e.key" type="checkbox" :checked="vs[e.key] || (e.lock ? lim(e.lock) : false)"
          :disabled="e.lock ? lim(e.lock) : false" @change="vs[e.key] = $event.target.checked" />
      </div>
      <div v-else-if="e.type === 'range'" class="ctl-row">
        <label :for="'x-' + e.key">{{ e.label }}</label>
        <input :id="'x-' + e.key" v-model.number="vs[e.key]" type="range" :min="lim(e.min)" :max="lim(e.max)" :step="e.step" />
        <span class="val">{{ fmt(vs[e.key], e.step) }}</span>
      </div>
      <template v-else-if="e.type === 'marker'">
        <div class="ctl-row">
          <label :for="'x-' + e.key">{{ e.label }}</label>
          <select class="mini" :value="vs[e.select.key]" :aria-label="e.label + ' variable'"
            @change="onAnimChange(e.select.key, $event)">
            <option v-for="o in animOpts" :key="o.value" :value="o.value">{{ o.label }}</option>
            <option value="__new__">+ new…</option>
          </select>
          <button v-if="String(vs[e.select.key]).startsWith('u:')" class="playbtn del"
            title="delete this animation mode" @click="deleteCurrent(e.select.key)">
            ✕
          </button>
          <button class="playbtn" :title="vs.play ? 'pause' : 'play'" @click="vs.play = !vs.play">
            {{ vs.play ? '❚❚' : '▶' }}
          </button>
          <button class="playbtn reset" title="reset to experiment defaults" @click="applyPreset(store.basePreset)">
            ↺
          </button>
          <input :id="'x-' + e.key" v-model.number="vs[e.key]" type="range" min="0" max="1" :step="e.step" />
        </div>
        <!-- "+ new" mode form: one parameter + a custom range, saved for all views -->
        <div v-if="creating" class="newmode">
          <div class="ctl-row">
            <label>parameter</label>
            <select :value="draft.param" @change="setDraftParam($event.target.value)">
              <option v-for="(p, k) in ANIM_PARAMS" :key="k" :value="k">{{ p.short }} [{{ p.unit }}]</option>
            </select>
          </div>
          <div class="ctl-row">
            <label>range</label>
            <input v-model.number="draft.min" type="number" step="any" aria-label="range min" />
            <span class="dash">–</span>
            <input v-model.number="draft.max" type="number" step="any" aria-label="range max" />
          </div>
          <div class="ctl-row">
            <label>name</label>
            <input v-model="draft.name" type="text"
              :placeholder="`${ANIM_PARAMS[draft.param].short} ${draft.min}–${draft.max}`" />
          </div>
          <div class="ctl-row btns">
            <button class="linkish" @click="creating = false">cancel</button>
            <button class="linkish save" :disabled="!(+draft.min < +draft.max)" @click="saveDraft(e.select.key)">
              save
            </button>
          </div>
        </div>
      </template>
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
.playbtn.del { color: var(--muted); font-size: 10px; }
.playbtn.del:hover { color: var(--accent); }
select.mini { flex: 0 0 62px; }
select {
  background: var(--surface-3);
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 3px 6px;
  flex: 2;
  min-width: 0;
}
.newmode {
  margin: 4px 0 6px;
  padding: 6px 8px;
  border: 1px solid var(--border);
  border-radius: 8px;
}
.newmode input[type='number'], .newmode input[type='text'] {
  flex: 1;
  min-width: 0;
  background: var(--surface-3);
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 3px 6px;
  color: inherit;
  font: inherit;
  font-size: 11px;
}
.newmode .dash { flex: none; color: var(--muted); }
.newmode .btns { justify-content: flex-end; }
.newmode .save:disabled { opacity: 0.4; cursor: default; }
</style>
