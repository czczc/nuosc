<script setup>
import { store, applyPreset } from './store.js';
import { PRESETS } from './engines/constants.js';
import { VIEWS, VIEW_MAP } from './views/index.js';
import StageHost from './components/StageHost.vue';
import ControlsCard from './components/ControlsCard.vue';
import ViewExtras from './components/ViewExtras.vue';
import CompanionPanel from './components/CompanionPanel.vue';

const presetNames = Object.keys(PRESETS);
</script>

<template>
  <div class="app">
    <header>
      <div class="wordmark">nuosc<span>3D NEUTRINO OSCILLATION</span></div>
      <nav class="seg views" aria-label="View">
        <button v-for="v in VIEWS" :key="v.id" :class="{ on: store.view === v.id }" @click="store.view = v.id">
          {{ v.label }}
        </button>
      </nav>
      <div class="right">
        <span v-for="p in presetNames" :key="p" class="chip" :class="{ on: store.preset === p }" role="button"
          tabindex="0" @click="applyPreset(p)" @keydown.enter="applyPreset(p)">{{ p }}</span>
        <span class="chip" :class="{ on: store.preset === 'custom' }">custom</span>
        <span class="seg">
          <button :class="{ on: store.ortho }" @click="store.ortho = true">ortho</button>
          <button :class="{ on: !store.ortho }" @click="store.ortho = false">persp</button>
        </span>
      </div>
    </header>

    <main>
      <StageHost :key="store.view" :view-def="VIEW_MAP[store.view]" />
      <aside>
        <ControlsCard />
        <ViewExtras :key="store.view" :view-def="VIEW_MAP[store.view]" />
        <CompanionPanel :key="store.view" :view-def="VIEW_MAP[store.view]" />
        <div class="colophon">
          engines: NuFast-LBL (<a href="https://arxiv.org/abs/2405.02400" target="_blank" rel="noopener">arXiv:2405.02400</a>)
          + exact 3-flavor · PDG 2023 · cross-validated to 1e-7
        </div>
      </aside>
    </main>
  </div>
</template>

<style scoped>
.app { height: 100%; display: flex; flex-direction: column; }
header {
  height: 48px; flex: none;
  display: flex; align-items: center; gap: 14px;
  background: var(--surface-2);
  border-bottom: 1px solid var(--border);
  padding: 0 14px;
}
.wordmark { font-family: var(--font-mono); font-weight: 700; font-size: 15px; }
.wordmark span {
  display: block; font-weight: 400; font-size: 8px;
  letter-spacing: 0.09em; color: var(--muted);
}
.right { margin-left: auto; display: flex; gap: 8px; align-items: center; }
main { flex: 1; display: flex; min-height: 0; }
aside {
  width: 340px; flex: none;
  background: var(--surface-2);
  border-left: 1px solid var(--border);
  overflow-y: auto;
  display: flex; flex-direction: column;
}
.colophon {
  margin-top: auto;
  padding: 10px 12px;
  font-size: 10.5px;
  color: var(--muted);
  border-top: 1px solid var(--border);
}
.colophon a { color: var(--muted); }
@media (max-width: 900px) {
  main { flex-direction: column; overflow-y: auto; }
  aside { width: 100%; border-left: none; border-top: 1px solid var(--border); }
}
</style>
