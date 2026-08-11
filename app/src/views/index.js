// View contract — every view module default-exports:
// {
//   id: string, label: string,
//   extras: [                       // per-view controls, rendered by ViewExtras.vue,
//     { key, type: 'select', label, options: [{value,label}] }       // bound to store.views[id][key]
//     { key, type: 'checkbox', label }
//     { key, type: 'range', label, min, max, step, unit? }
//   ],
//   create(container, store) -> {
//     base: SceneBase,              // StageHost drives base.setOrtho(store.ortho)
//     update(),                     // heavy rebuild; runs inside watchEffect — every reactive read
//                                   // it makes becomes a dependency. Must NOT read play/marker.
//     tick?(dt),                    // cheap per-frame work (play/marker animation); reads store freely
//     probe?(event) -> string|null, // hover readout for the stage chip
//     dispose(),
//   },
//   companion: {                    // the live 2D projection panel (present in ALL views)
//     title(store) -> string,
//     draw(canvas, store),          // canvas 2D; runs inside watchEffect (same dependency rule);
//                                   // marker-driven repaints are triggered by CompanionPanel's rAF
//                                   // only if markerDriven is true
//     markerDriven?: boolean,
//   },
// }
import oscillogram from './oscillogram.js';
import tube from './tube.js';
import sphere from './sphere.js';
import phasors from './phasors.js';
import biprob from './biprob.js';

export const VIEWS = [oscillogram, tube, sphere, phasors, biprob];
export const VIEW_MAP = Object.fromEntries(VIEWS.map((v) => [v.id, v]));
