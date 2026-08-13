<script setup>
import { store, applyPreset, setTheme } from './store.js';
import { PRESETS } from './engines/constants.js';
import { PALETTES } from './three/SceneBase.js';
import { VIEWS, VIEW_MAP } from './views/index.js';
import StageHost from './components/StageHost.vue';
import ControlsCard from './components/ControlsCard.vue';
import ViewExtras from './components/ViewExtras.vue';
import CompanionPanel from './components/CompanionPanel.vue';
import FaqPage from './components/FaqPage.vue';

const presetNames = Object.keys(PRESETS);
const PALETTE_NAMES = Object.keys(PALETTES);

const companions = (v) => [v.companion, v.companion2].filter(Boolean);
</script>

<template>
  <div class="app">
    <header>
      <img class="logo" src="./assets/logo.svg" alt="nuglass" />
      <nav class="tabs" aria-label="View">
        <button v-for="v in VIEWS" :key="v.id" :class="{ on: store.view === v.id && !store.faq }"
          @click="store.view = v.id; store.faq = null">
          {{ v.label }}
        </button>
      </nav>
      <div class="right">
        <div class="group" role="group" aria-label="Experiment">
          <span v-for="p in presetNames" :key="p" class="chip" :class="{ on: store.preset === p }" role="button"
            tabindex="0" @click="applyPreset(p)" @keydown.enter="applyPreset(p)">{{ p }}</span>
          <span class="chip" :class="{ on: store.preset === 'custom' }">custom</span>
        </div>
        <div class="group" role="group" aria-label="Display">
          <button class="faqbtn" :class="{ on: !!store.faq }" @click="store.faq = store.faq ? null : 'engines'">
            FAQ
          </button>
          <select v-model="store.palette" class="palette" aria-label="Color palette" title="color palette">
            <option v-for="p in PALETTE_NAMES" :key="p" :value="p">{{ p }}</option>
          </select>
          <button class="themebtn" :title="store.theme === 'dark' ? 'switch to light mode' : 'switch to dark mode'"
            @click="setTheme(store.theme === 'dark' ? 'light' : 'dark')">
            {{ store.theme === 'dark' ? '☀' : '☾' }}
          </button>
        </div>
      </div>
    </header>

    <main>
      <FaqPage v-if="store.faq" />
      <template v-else>
      <StageHost :key="store.view + '-' + store.theme" :view-def="VIEW_MAP[store.view]" />
      <aside>
        <ControlsCard />
        <ViewExtras :key="store.view" :view-def="VIEW_MAP[store.view]" />
        <CompanionPanel v-for="(c, i) in companions(VIEW_MAP[store.view])" :key="store.view + '-' + i"
          :companion="c" />
      </aside>
      </template>
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
.tabs { display: flex; }
.tabs button {
  position: relative;
  padding: 4px 12px 5px;
  font-family: var(--font-mono);
  font-size: 14px;
  color: var(--muted);
  background: none;
  border: none;
  cursor: pointer;
}
.tabs button:hover { color: var(--text); }
.group {
  display: flex; align-items: center; gap: 2px;
  padding: 3px;
  border: 1px solid var(--border);
  border-radius: 999px;
  background: var(--surface-3);
}
.group .chip { border-color: transparent; background: transparent; }
.group .chip.on { background: var(--accent); border-color: var(--accent); }
.faqbtn {
  font-family: var(--font-mono);
  font-size: 11.5px;
  color: var(--muted);
  background: none;
  border: none;
  border-radius: 999px;
  padding: 3px 12px;
  cursor: pointer;
}
.faqbtn:hover { color: var(--text); }
.faqbtn.on { color: var(--accent-ink); background: var(--accent); font-weight: 700; }
.themebtn {
  width: 28px; height: 28px;
  border: none; border-radius: 50%;
  background: none; color: var(--muted);
  font-size: 17px; line-height: 1;
  padding: 0; cursor: pointer;
}
.themebtn:hover { color: var(--accent); }
.palette {
  font-family: var(--font-mono);
  font-size: 11px;
  color: var(--muted);
  background: transparent;
  border: none;
  border-radius: 999px;
  padding: 3px 8px;
  cursor: pointer;
}
.palette:hover { color: var(--text); }
.tabs button.on { color: var(--accent); font-weight: 600; }
.tabs button.on::after {
  content: '';
  position: absolute;
  left: 12px; right: 12px; bottom: 0;
  height: 2px;
  background: var(--accent);
}
.logo { height: 38px; display: block; }
.right { margin-left: auto; display: flex; gap: 8px; align-items: center; }
main { flex: 1; display: flex; min-height: 0; }
aside {
  width: 340px; flex: none;
  background: var(--surface-2);
  border-left: 1px solid var(--border);
  overflow-y: auto;
  display: flex; flex-direction: column;
}
@media (max-width: 900px) {
  main { flex-direction: column; overflow-y: auto; }
  aside { width: 100%; border-left: none; border-top: 1px solid var(--border); }
}
</style>
