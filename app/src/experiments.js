// User-defined experiments: a reactive name -> definition map persisted in
// localStorage and mirrored into USER_PRESETS so the preset machinery
// (applyPreset, lRangeOf/eRangeOf, the reset button) resolves them like built-ins.
// Definition shape: { L, rho, Emin, Emax, Epeak, th12, th13, th23, dm21, dm31, dcp, Ye }.
import { reactive } from 'vue';
import { USER_PRESETS } from './engines/constants.js';

const KEY = 'nuglass-experiments';

export const userExps = reactive({});

function sync() {
  for (const k of Object.keys(USER_PRESETS)) delete USER_PRESETS[k];
  for (const [name, d] of Object.entries(userExps)) {
    USER_PRESETS[name] = {
      L: d.L, rho: d.rho, Erange: [d.Emin, d.Emax], Epeak: d.Epeak,
      params: { th12: d.th12, th13: d.th13, th23: d.th23, dm21: d.dm21, dm31: d.dm31, dcp: d.dcp, Ye: d.Ye },
    };
  }
}

function persist() {
  try { localStorage.setItem(KEY, JSON.stringify(userExps)); } catch { /* private mode */ }
}

export function loadUserExps() {
  try {
    const obj = JSON.parse(localStorage.getItem(KEY) ?? '{}');
    if (obj && typeof obj === 'object') Object.assign(userExps, obj);
  } catch { /* corrupted entry -> start empty */ }
  sync();
}

export function saveUserExp(name, def) {
  userExps[name] = { ...def };
  persist();
  sync();
}

export function deleteUserExp(name) {
  delete userExps[name];
  persist();
  sync();
}
