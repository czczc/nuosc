import { createApp } from 'vue';
import App from './App.vue';
import { store, setTheme } from './store.js';
import { router } from './router.js';
import { loadUserExps } from './experiments.js';
import './style.css';

loadUserExps();
setTheme(store.theme);
const app = createApp(App).use(router);
// wait for the initial route so a deep link doesn't first mount the default view
router.isReady().then(() => app.mount('#app'));
