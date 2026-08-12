<script setup>
import { onMounted, onBeforeUnmount, ref, watch, watchEffect } from 'vue';
import { store } from '../store.js';

const props = defineProps({ viewDef: { type: Object, required: true } });
const host = ref(null);
const chip = ref('');

let instance = null;
let stopUpdate = null;
let stopOrtho = null;
let chipTimer = null;
let lastEvent = null;

function onMove(ev) {
  lastEvent = ev;
  refreshChip();
}

function snapTo(plane) {
  instance?.base?.snapTo(plane);
}

function refreshChip() {
  if (!instance?.probe) return;
  const s = instance.probe(lastEvent ?? new MouseEvent('mousemove'));
  if (s) chip.value = s;
}

onMounted(() => {
  instance = props.viewDef.create(host.value, store);
  stopUpdate = watchEffect(() => instance.update());
  stopOrtho = watch(() => store.ortho, (v) => instance.base?.setOrtho(v));
  if (instance.tick) instance.base.onTick = (dt) => instance.tick(dt);
  chipTimer = setInterval(refreshChip, 200);
  refreshChip();
});

onBeforeUnmount(() => {
  clearInterval(chipTimer);
  stopUpdate?.();
  stopOrtho?.();
  instance?.dispose();
  instance = null;
});
</script>

<template>
  <div ref="host" class="stage" @pointermove="onMove">
    <span class="seg proj" aria-label="Projection">
      <button v-for="p in ['yx', 'zx', 'zy']" :key="p" @click="snapTo(p)">
        {{ p[0] }}-{{ p[1] }}
      </button>
    </span>
    <div v-if="chip" class="hoverchip">{{ chip }}</div>
  </div>
</template>

<style scoped>
.stage { position: relative; flex: 1; min-width: 0; min-height: 320px; background: var(--stage-bg); }
.stage :deep(canvas) { display: block; }
.proj { position: absolute; left: 10px; top: 10px; z-index: 5; }
.proj button { background: var(--chip-bg); }
.hoverchip {
  position: absolute; right: 10px; bottom: 10px; z-index: 5;
  font-family: var(--font-mono); font-size: 11px;
  font-variant-numeric: tabular-nums;
  background: var(--chip-bg); color: var(--chip-ink);
  padding: 4px 8px; border-radius: 5px;
  pointer-events: none;
  max-width: 90%;
}
</style>
