<script setup>
import { computed, onMounted, onBeforeUnmount, ref, watchEffect } from 'vue';
import { store } from '../store.js';

const props = defineProps({ viewDef: { type: Object, required: true } });
const canvas = ref(null);
const title = computed(() => props.viewDef.companion.title(store));

let stopDraw = null;
let raf = 0;
let ro = null;

function drawNow() {
  if (canvas.value) props.viewDef.companion.draw(canvas.value, store);
}

onMounted(() => {
  if (props.viewDef.companion.markerDriven) {
    // marker moves every animation frame during play — repaint on rAF, read store non-reactively
    const loop = () => { drawNow(); raf = requestAnimationFrame(loop); };
    raf = requestAnimationFrame(loop);
  } else {
    stopDraw = watchEffect(drawNow);
  }
  ro = new ResizeObserver(drawNow);
  ro.observe(canvas.value);
});

onBeforeUnmount(() => {
  cancelAnimationFrame(raf);
  stopDraw?.();
  ro?.disconnect();
});
</script>

<template>
  <div class="card companion">
    <h3>{{ title }}</h3>
    <div class="frame">
      <canvas ref="canvas"></canvas>
    </div>
  </div>
</template>

<style scoped>
.companion .frame {
  position: relative;
  height: 180px;
  border-radius: 6px;
  overflow: hidden;
  background: #0b0e13;
}
canvas { position: absolute; inset: 0; width: 100%; height: 100%; }
</style>
