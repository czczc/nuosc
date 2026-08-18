// Hash-based routes (GitHub Pages friendly): #/<view> for the views,
// #/faq/<section?> for the FAQ page. The URL is the source of truth for
// view/faq navigation — a guard syncs it into the reactive store, which the
// components keep reading as before.
import { createRouter, createWebHashHistory } from 'vue-router';
import { store, applyPreset } from './store.js';
import { VIEW_MAP } from './views/index.js';

const VIEW_IDS = Object.keys(VIEW_MAP).join('|');
// App.vue renders from the store (no <router-view>); routes carry a null component
const none = { render: () => null };

export const router = createRouter({
  history: createWebHashHistory(),
  routes: [
    { path: '/', redirect: `/${store.view}` },
    { path: '/faq/:section?', name: 'faq', component: none },
    { path: '/experiments', name: 'experiments', component: none },
    { path: `/:view(${VIEW_IDS})`, name: 'view', component: none },
    { path: '/:pathMatch(.*)*', redirect: '/' },
  ],
});

router.beforeEach((to) => {
  if (to.name === 'faq') {
    store.faq = to.params.section || 'engines';
    store.exps = false;
  } else if (to.name === 'experiments') {
    store.exps = true;
    store.faq = null;
  } else if (to.name === 'view') {
    // switching to a different view starts it from the experiment's initial state
    // (animations leave the shared sliders wherever they stopped)
    if (to.params.view !== store.view) applyPreset(store.basePreset);
    store.view = to.params.view;
    store.faq = null;
    store.exps = false;
  }
});
