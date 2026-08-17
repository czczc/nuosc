<script setup>
import { computed, ref } from 'vue';
import { store, applyPreset, setTheme, setChannel } from './store.js';
import { router } from './router.js';
import { PRESETS, CHANNELS, channelLabel, channelsOf } from './engines/constants.js';
import { VIEWS, VIEW_MAP } from './views/index.js';
import StageHost from './components/StageHost.vue';
import ControlsCard from './components/ControlsCard.vue';
import ViewExtras from './components/ViewExtras.vue';
import CompanionPanel from './components/CompanionPanel.vue';
import FaqPage from './components/FaqPage.vue';
import ExperimentsPage from './components/ExperimentsPage.vue';
import logoLight from './assets/logo-light.png';
import logoDark from './assets/logo-dark.png';

const channelIds = Object.keys(CHANNELS);
// experiment chips: only presets that measure the active channel
const presetNames = computed(() =>
  Object.keys(PRESETS).filter((p) => PRESETS[p].channels.includes(store.channel)));

const companions = (v) => [v.companion, v.companion2].filter(Boolean);
const vs = computed(() => store.views[store.view]); // current view's play/marker state
// a loaded user experiment gets its own chip after the built-ins
const userChip = computed(() => (PRESETS[store.basePreset] ? null : store.basePreset));
// tweaked sliders restyle the active chip instead of showing a separate "custom" chip
const isCustom = computed(() => store.preset === 'custom');
const chipTitle = (p) => (isCustom.value && store.basePreset === p
  ? `parameters modified — click to reset to ${p} defaults` : p);
const logoSrc = computed(() => (store.theme === 'light' ? logoLight : logoDark));

// channel rail: chips the active experiment doesn't measure are dimmed, and
// their tooltip warns about the experiment switch the click would trigger
const chanDim = (c) => !channelsOf(store.basePreset).includes(c);
function chanTitle(c) {
  if (!chanDim(c)) return null;
  const fb = Object.keys(PRESETS).find((n) => PRESETS[n].channels.includes(c));
  if (!fb) return null;
  const p = PRESETS[fb];
  return `not measured at ${store.basePreset} — switches to ${fb} (${p.anti ? 'reactor ν̄, ' : ''}L = ${p.L} km)`;
}
// when a channel click forces an experiment switch, pulse the new chip
const flashing = ref('');
function pickChannel(c) {
  const before = store.basePreset;
  setChannel(c);
  if (store.basePreset !== before) {
    flashing.value = ''; // retrigger the CSS animation on repeat switches
    requestAnimationFrame(() => { flashing.value = store.basePreset; });
    setTimeout(() => { flashing.value = ''; }, 2000);
  }
}
</script>

<template>
  <div class="app">
    <header>
      <img class="logo" :src="logoSrc" alt="NuGlass" />
      <nav class="tabs" aria-label="View">
        <button v-for="v in VIEWS" :key="v.id" :class="{ on: store.view === v.id && !store.faq && !store.exps }"
          @click="router.push('/' + v.id)">
          {{ v.label }}
        </button>
      </nav>
      <div v-if="!store.faq && !store.exps" class="group" role="group" aria-label="Animation">
        <button class="navbtn play" :title="vs.play ? 'pause' : 'play'" @click="vs.play = !vs.play">
          {{ vs.play ? '❚❚' : '▶' }}
        </button>
        <button class="navbtn" title="reset to experiment defaults" @click="applyPreset(store.basePreset)">
          ↺
        </button>
      </div>
      <div class="right">
        <div class="group" role="group" aria-label="Channel">
          <span v-for="c in channelIds" :key="c" class="chip" :class="{ on: store.channel === c, dim: chanDim(c) }"
            role="button" tabindex="0" :title="chanTitle(c)" @click="pickChannel(c)"
            @keydown.enter="pickChannel(c)">{{ channelLabel(c) }}</span>
        </div>
        <span class="vsep" aria-hidden="true"></span>
        <div class="group" role="group" aria-label="Experiment">
          <span v-for="p in [...presetNames, ...(userChip ? [userChip] : [])]" :key="p" class="chip"
            :class="{ on: store.basePreset === p, custom: isCustom && store.basePreset === p, flash: flashing === p }"
            role="button" tabindex="0" :title="chipTitle(p)" @click="applyPreset(p)" @keydown.enter="applyPreset(p)">
            {{ p }}{{ isCustom && store.basePreset === p ? '*' : '' }}</span>
          <span class="chip" :class="{ on: store.exps }" role="button" tabindex="0" title="define your own experiment"
            @click="router.push(store.exps ? '/' + store.view : '/experiments')"
            @keydown.enter="router.push(store.exps ? '/' + store.view : '/experiments')">+</span>
        </div>
        <div class="group" role="group" aria-label="Display">
          <button class="themebtn" :title="store.theme === 'dark' ? 'switch to light mode' : 'switch to dark mode'"
            @click="setTheme(store.theme === 'dark' ? 'light' : 'dark')">
            {{ store.theme === 'dark' ? '☀' : '☾' }}
          </button>
          <button class="faqbtn" :class="{ on: !!store.faq }" title="FAQ"
            @click="router.push(store.faq ? '/' + store.view : '/faq')">
            ?
          </button>
        </div>
      </div>
    </header>

    <main>
      <ExperimentsPage v-if="store.exps" />
      <FaqPage v-else-if="store.faq" />
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
/* active experiment with tweaked parameters: hollow + dashed, until reset/reapplied */
.group .chip.on.custom {
  background: transparent;
  border: 1px dashed var(--accent);
  color: var(--accent);
  font-weight: 700;
}
/* channels the active experiment doesn't measure: dimmed; the tooltip explains
   that clicking switches experiments */
.group .chip.dim { opacity: 0.4; }
.group .chip.dim:hover { opacity: 0.8; }
/* the experiment chip a channel click just switched to: expanding ring pulse */
.group .chip.flash { animation: chip-flash 0.9s ease-out 2; }
@keyframes chip-flash {
  0% { box-shadow: 0 0 0 0 var(--accent); }
  100% { box-shadow: 0 0 0 9px transparent; }
}
.navbtn {
  width: 26px; height: 24px;
  border: none; border-radius: 999px;
  background: none; color: var(--accent);
  font-size: 15px; line-height: 1;
  padding: 0; cursor: pointer;
}
.navbtn.play { font-size: 10px; }
.navbtn:hover { background: var(--surface-2); }
.faqbtn {
  width: 28px; height: 28px;
  border: 1px solid var(--border); border-radius: 50%;
  background: none; color: var(--muted);
  font-family: var(--font-mono);
  font-size: 15px; line-height: 1;
  padding: 0; cursor: pointer;
}
.faqbtn:hover { color: var(--accent); border-color: var(--accent); }
.faqbtn.on { color: var(--accent-ink); background: var(--accent); border-color: var(--accent); font-weight: 700; }
.themebtn {
  width: 28px; height: 28px;
  border: none; border-radius: 50%;
  background: none; color: var(--muted);
  font-size: 17px; line-height: 1;
  padding: 0; cursor: pointer;
}
.themebtn:hover { color: var(--accent); }
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
/* divider between the channel and experiment groups (the pill borders alone are
   too subtle in dark mode to read as separate controls) */
.vsep { width: 1px; height: 22px; flex: none; background: var(--muted); opacity: 0.4; }
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
  header { height: auto; flex-wrap: wrap; gap: 6px 10px; padding: 6px 10px; }
  .logo { height: 30px; }
  .tabs {
    order: 3;
    flex: 1 1 100%;
    overflow-x: auto;
    scrollbar-width: none;
  }
  .tabs::-webkit-scrollbar { display: none; }
  .tabs button { padding: 4px 8px 5px; font-size: 13px; white-space: nowrap; }
  .right { flex-wrap: wrap; justify-content: flex-end; gap: 6px; }
  .group .chip { padding: 3px 7px; }
}
</style>
