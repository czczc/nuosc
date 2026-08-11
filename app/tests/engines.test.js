// Cross-validation of the two physics engines (the wayfinder map's build-time gate, issue #8).
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { probabilityMatter } from '../src/engines/nufast.js';
import { hamiltonian, eigH, prob } from '../src/engines/jacobi.js';
import { engineParams, DEFAULTS } from '../src/engines/constants.js';

const base = { ...DEFAULTS, normalOrdering: true };

test('unitarity of the Jacobi engine', () => {
  const ep = engineParams(base);
  const eig = eigH(hamiltonian(ep, 2.5, 3.0));
  for (const L of [0, 500, 1300, 5000]) {
    const s = prob(eig, 1, 0, L) + prob(eig, 1, 1, L) + prob(eig, 1, 2, L);
    assert.ok(Math.abs(s - 1) < 1e-12, `sum=${s} at L=${L}`);
  }
});

test('two-flavor vacuum limit matches the analytic formula', () => {
  const ep = engineParams({ ...base, th12: 0, th13: 0, dm21: 0 });
  const eig = eigH(hamiltonian(ep, 2.5, 0));
  const got = prob(eig, 1, 1, 1300);
  const want = 1 - Math.sin(2 * ep.th23) ** 2 * Math.sin(1.2669 * ep.dm31 * 1300 / 2.5) ** 2;
  assert.ok(Math.abs(got - want) < 1e-4, `got=${got} want=${want}`);
});

test('NuFast (N_Newton=1) agrees with Jacobi across the grid to 5e-7', () => {
  let worst = 0;
  for (const E of [0.3, 0.5, 1.0, 2.5, 5.0])
    for (const L of [295, 810, 1300, 3000])
      for (const rho of [0, 2.85, 8.0])
        for (const dcp of [0, 90, 217, 270])
          for (const anti of [false, true])
            for (const no of [true, false]) {
              const ep = engineParams({ ...base, dcp, normalOrdering: no });
              const eig = eigH(hamiltonian(ep, E, rho, anti));
              const P = probabilityMatter(ep.s12sq, ep.s13sq, ep.s23sq, ep.delta, ep.dm21, ep.dm31, L, anti ? -E : E, rho * ep.Ye, 1);
              for (let a = 0; a < 3; a++)
                for (let b = 0; b < 3; b++)
                  worst = Math.max(worst, Math.abs(P[a][b] - prob(eig, a, b, L)));
            }
  assert.ok(worst < 5e-7, `max |dP| = ${worst}`);
});

test('DUNE reference point', () => {
  const ep = engineParams(base);
  const P = probabilityMatter(ep.s12sq, ep.s13sq, ep.s23sq, ep.delta, ep.dm21, ep.dm31, 1300, 2.5, 1.5, 1);
  assert.ok(Math.abs(P[1][0] - 0.08437) < 2e-4, `P(mu->e)=${P[1][0]}`);
  const Pb = probabilityMatter(ep.s12sq, ep.s13sq, ep.s23sq, ep.delta, ep.dm21, ep.dm31, 1300, -2.5, 1.5, 1);
  assert.ok(Math.abs(Pb[1][0] - 0.0223) < 2e-3, `Pbar(mu->e)=${Pb[1][0]}`);
});
