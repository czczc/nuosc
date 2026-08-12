# nuglass — 3D Neutrino Oscillation

Interactive five-view 3D visualization of neutrino oscillation with matter effects.
Built per the locked design spec ([issue #9](https://github.com/czczc/nuglass/issues/9)),
charted by the wayfinder map ([issue #1](https://github.com/czczc/nuglass/issues/1)).

## Views

| View | Shows | Engine |
|---|---|---|
| Oscillogram | P(νμ→νe) surface over E × (L / δCP / ρ) | NuFast |
| Flavor tube | flavor composition flowing along the baseline | Jacobi |
| State sphere | νe–νμ Bloch projection of the propagating state | Jacobi |
| Phasors | matter-eigenstate amplitudes interfering | Jacobi |
| Biprobability | Minakata–Nunokawa P vs P̄ rings, NO + IO | NuFast |

Every view pairs with a live 2D companion panel (its natural projection).

## Physics

- `src/engines/nufast.js` — JS port of [NuFast-LBL](https://github.com/PeterDenton/NuFast-LBL)
  (Denton & Parke, [arXiv:2405.02400](https://arxiv.org/abs/2405.02400), MIT): all nine
  probabilities in constant-density matter, N_Newton=1.
- `src/engines/jacobi.js` — exact complex-Hermitian 3×3 diagonalization: amplitudes,
  phases, and matter-eigenstate decompositions.
- Cross-validated against each other to ≤5e-7 in `tests/engines.test.js` (run `npm test`).

## Develop

```sh
npm install
npm run dev     # local dev server
npm test        # engine cross-validation
npm run build   # production build (served at /nuglass/ on GitHub Pages)
```
