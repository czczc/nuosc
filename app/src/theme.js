// Palette for the 3D scenes and canvas plots. DOM colors live in style.css as
// CSS variables switched by [data-theme]; this module serves the code-drawn side.
// Reading theme() inside a watchEffect-driven draw makes it repaint on toggle;
// the 3D stages are re-created on toggle (StageHost is keyed on store.theme).
import { store } from './store.js';

const THEMES = {
  dark: {
    stage: 0x0b0e13, canvas: '#0b0e13',
    grid1: 0x334455, grid2: 0x223344,
    axis: 0x667788,
    hi: 0xffffff, hiCss: '#ffffff',
    pCurve: 0x66ccff, pCurveCss: '#66ccff',
    ink: '#9aa7b5',
    label: '#9aa7b5', // 3D scene text sprites
    plotGrid: 'rgba(140,155,170,0.15)', plotFrame: 'rgba(140,155,170,0.5)',
    legendBg: 'rgba(11,14,19,0.75)',
    beam: 'rgba(78,205,180,0.55)', beamText: 'rgba(78,205,180,0.8)',
  },
  light: {
    stage: 0xf4f6f9, canvas: '#fdfdfe',
    grid1: 0xc3ccd6, grid2: 0xdde3ea,
    axis: 0x7c8b9a,
    hi: 0x1d242d, hiCss: '#1d242d',
    pCurve: 0x1a7fc4, pCurveCss: '#1a7fc4',
    ink: '#5f6b78',
    label: '#14181d', // 3D scene text sprites: near-black for contrast against the pale stage
    plotGrid: 'rgba(70,90,110,0.15)', plotFrame: 'rgba(70,90,110,0.55)',
    legendBg: 'rgba(253,253,254,0.8)',
    beam: 'rgba(13,125,104,0.55)', beamText: 'rgba(13,125,104,0.85)',
  },
};

export function theme() {
  return THEMES[store.theme] ?? THEMES.dark;
}
