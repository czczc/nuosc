import { createApp } from 'vue';
import App from './App.vue';
import { store, setTheme } from './store.js';
import './style.css';

setTheme(store.theme);
createApp(App).mount('#app');
