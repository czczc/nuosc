<script setup>
import { computed, onMounted, onBeforeUnmount, ref, watchEffect } from 'vue';
import { store } from '../store.js';

const props = defineProps({ companion: { type: Object, required: true } });
const canvas = ref(null);
const title = computed(() => props.companion.title(store));

let stopDraw = null;
let raf = 0;
let ro = null;

function drawNow() {
  if (canvas.value) props.companion.draw(canvas.value, store);
}

onMounted(() => {
  if (props.companion.markerDriven) {
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
  height: 150px;
  border-radius: 6px;
  overflow: hidden;
  background: var(--stage-bg);
}
canvas { position: absolute; inset: 0; width: 100%; height: 100%; }
</style>
