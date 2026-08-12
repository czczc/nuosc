# nuosc

Interactive 3D visualization of long-baseline neutrino oscillation (νμ → νe), built with Vue 3 and three.js.

**Live site: https://czczc.github.io/nuosc/**

The app offers five linked views of the same physics — an oscillogram surface, a bi-probability plot, a flavor-state sphere, an amplitude phasor diagram, and a flavor-composition tube — each with a companion 2D plot and experiment presets (DUNE, NOvA, T2K). Probabilities come from two independent engines: a JavaScript port of NuFast-LBL ([Denton & Parke, arXiv:2405.02400](https://arxiv.org/abs/2405.02400)) and an exact complex-Jacobi amplitude engine, cross-validated against each other. See the in-app FAQ page for the physics details.

## Local setup

Requires Node.js 20+.

```
cd app
npm install
npm run dev     # dev server at http://localhost:5173/nuosc/
npm test        # engine unit tests
npm run build   # production build to app/dist
```

Pushes to `main` that touch `app/` are automatically tested, built, and deployed to GitHub Pages via `.github/workflows/deploy.yml`.

## Repository layout

- `app/` — the web app (Vite + Vue 3 + three.js)
- `legacy/` — the original tools this app replaced: a 2D Chart.js page, a ROOT plotting macro, and a wrapper around [Brett Viren's nuosc](https://github.com/brettviren/nuosc) code. Kept for reference, no longer maintained.
