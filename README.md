# NuGlass

A fast and interactive 3D visualization of neutrino oscillations, built with Vue 3 and three.js.

**Live site: https://czczc.github.io/nuglass/**

The app offers six linked views of the same physics (oscillogram, biprobability, statesphere, phasors, flavortube, and worldline), each paired with live 2D plots, switchable between oscillation channels (νμ→νe, νμ→νμ, νμ→ντ, νe→νe), with experiment presets spanning accelerator and reactor scales (DUNE, NOvA, T2K, JUNO, KamLAND, Daya Bay). All calculations are performed live in the browser by two independent engines: a JavaScript port of NuFast-LBL ([Denton & Parke, arXiv:2405.02400](https://arxiv.org/abs/2405.02400)) and an exact amplitude engine that tracks the full quantum state, cross-checked against each other to better than one part in ten million. See the in-app FAQ page for the physics details.

## Local setup

Requires Node.js 20+.

```
cd app
npm install
npm run dev     # dev server at http://localhost:5173/nuglass/
npm test        # engine unit tests
npm run build   # production build to app/dist
```

Pushes to `main` that touch `app/` are automatically tested, built, and deployed to GitHub Pages via `.github/workflows/deploy.yml`.

## How to cite

If you use NuGlass in a paper or talk, please cite it as:

```bibtex
@misc{nuglass,
  author       = {Zhang, Chao},
  title        = {NuGlass: a fast and interactive 3D visualization of neutrino oscillations},
  year         = {2026},
  howpublished = {\url{https://czczc.github.io/nuglass/}},
  note         = {Source code: \url{https://github.com/czczc/nuglass}}
}
```

## Repository layout

- `app/`: the web app (Vite + Vue 3 + three.js)
- `legacy/`: the original tools this app replaced (a 2D Chart.js page, a ROOT plotting macro, and a wrapper around [Brett Viren's nuosc](https://github.com/brettviren/nuosc) code). Kept for reference, no longer maintained.
