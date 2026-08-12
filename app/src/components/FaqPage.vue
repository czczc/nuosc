<script setup>
import { onMounted } from 'vue';
import { store } from '../store.js';

function openView(id) {
  store.view = id;
  store.faq = null;
}

onMounted(() => {
  if (typeof store.faq === 'string') {
    document.getElementById('faq-' + store.faq)?.scrollIntoView({ block: 'start' });
  }
});
</script>

<template>
  <div class="faq-page">
    <article class="doc">
      <h2>FAQ</h2>
      <p class="sub">
        nuosc visualizes three-flavor neutrino oscillation in constant-density matter — five 3D views,
        each paired with a live 2D companion plot. This page describes the physics engines and what each
        view is plotting.
      </p>

      <section id="faq-engines">
        <h3>The physics engines</h3>
        <p>
          Everything on screen is computed live from the same physics: the 3×3 flavor Hamiltonian in
          constant-density matter, with the standard matter potential
          <span class="k">a = 2√2 G_F N_e E</span>. Oscillation parameters default to PDG 2023
          (θ₁₂ = 33.82°, θ₁₃ = 8.61°, θ₂₃ = 48.3°, Δm²₂₁ = 7.39×10⁻⁵ eV², |Δm²₃₁| = 2.525×10⁻³ eV²,
          δCP = 217°), all adjustable under “all parameters”. Antineutrinos are computed by flipping the
          sign of E, which conjugates both the matter potential and the CP phase.
        </p>
        <p>Two engines share the work:</p>
        <p>
          <strong>NuFast-LBL</strong> — a JavaScript port of Denton &amp; Parke’s algorithm
          (<a href="https://arxiv.org/abs/2405.02400" target="_blank" rel="noopener">arXiv:2405.02400</a>,
          <a href="https://github.com/PeterDenton/NuFast-LBL" target="_blank" rel="noopener">reference code</a>).
          It gets the matter eigenvalues from one Newton refinement of the characteristic polynomial and
          produces all nine oscillation probabilities in pure real arithmetic, at roughly 100 ns per point.
          It powers the grid-heavy views: the oscillogram surface and the biprobability ring stack.
        </p>
        <p>
          <strong>Exact amplitude engine</strong> — builds the complex Hamiltonian and diagonalizes it with
          a Jacobi eigensolver, giving full transition <em>amplitudes</em> and phases, not just
          probabilities. The state sphere and phasor views need this: their content (coherences, phasor
          arms) cannot be reconstructed from probabilities alone.
        </p>
        <p>
          The two engines are cross-validated against each other to better than 10⁻⁷ in probability over an
          (E, L) grid (10⁻¹⁰-level with a second Newton iteration).
        </p>
        <p>
          <strong>Experiment presets</strong> — DUNE (1300 km, ρ = 2.85 g/cm³), NOvA (810 km, 2.84) and
          T2K (295 km, 2.6) each carry their beam window and flux-peak energy. Selecting one snaps every
          energy control to the peak, restricts the energy axes to the beam window, and scopes the baseline
          range to 0 – 2L, so the first oscillation maximum sits mid-plot. Dragging any slider afterwards
          shows the “custom” chip but keeps you inside that experiment’s ranges.
        </p>
      </section>

      <section id="faq-oscillogram">
        <h3>Oscillogram</h3>
        <p>
          A surface of the appearance probability <span class="k">P(νμ→νe)</span> over energy E (x axis)
          × a second axis of your choice — baseline L, CP phase δCP, or matter density ρ. Height and
          viridis color both encode the probability. The white line is the slice at the current slider
          value of the second-axis quantity, and it is exactly the curve shown in the companion spectrum.
        </p>
        <p>
          How to read it: the ridges are oscillation maxima (at fixed L they march to lower E as
          Δm²₃₁L/4E passes odd multiples of π/2); along the δCP axis you see the CP interference shift the
          peaks; along ρ you see the matter resonance enhance ν (normal ordering) or ν̄ (inverted).
          The companion plot overlays ν and ν̄ spectra at the current L, with the experiment’s flux peak
          marked by a dashed line.
        </p>
        <p><button class="linkish" @click="openView('oscillogram')">open this view →</button></p>
      </section>

      <section id="faq-biprob">
        <h3>Biprobability</h3>
        <p>
          The Minakata–Nunokawa biprobability plot
          (<a href="https://arxiv.org/abs/hep-ph/0108085" target="_blank" rel="noopener">arXiv:hep-ph/0108085</a>):
          neutrino appearance <span class="k">P(νμ→νe)</span> on x against antineutrino appearance
          <span class="k">P̄(ν̄μ→ν̄e)</span> on y, at fixed E and L. As δCP sweeps 0–360° the point traces
          an ellipse. The 3D view stacks these ellipses along energy, for both mass orderings
          (orange = normal, blue = inverted); the spheres mark the current δCP and the white connector
          shows the ordering difference at the slice energy.
        </p>
        <p>
          How to read it: displacement of the point from the <span class="k">P = P̄</span> diagonal is
          direct CP violation; the separation between the orange and blue ellipses is the mass-ordering
          signal (driven by the matter effect). Where the ellipses overlap, an experiment at that energy
          cannot distinguish the orderings from these channels alone.
        </p>
        <p><button class="linkish" @click="openView('biprob')">open this view →</button></p>
      </section>

      <section id="faq-sphere">
        <h3>State sphere</h3>
        <p>
          A Bloch-sphere picture of the evolving flavor state, projected onto the νe–νμ subspace: the north
          pole is pure νe, the south pole pure νμ, and the transverse components are the quantum coherences
          between them. The trajectory shows the state of an initial νμ evolving along L (or, in the other
          sweep modes, the state arriving at the detector as a function of E or δCP).
        </p>
        <p>
          How to read it: a genuine two-flavor oscillation would trace a circle on the sphere’s surface.
          The full three-flavor state lives in an 8-dimensional space, so this 3D shadow can shrink:
          a vector shorter than 1 means probability has leaked into ντ. The companion plot shows all three
          probabilities P(νμ→νe/μ/τ) along the sweep with the marker synchronized to the 3D animation.
        </p>
        <p><button class="linkish" @click="openView('sphere')">open this view →</button></p>
      </section>

      <section id="faq-phasors">
        <h3>Phasors</h3>
        <p>
          The νe appearance amplitude decomposed into its three matter-eigenstate contributions,
          <span class="k">A(L) = Σᵢ cᵢ e^(−iλᵢL)</span>, where the coefficients cᵢ come from the mixing
          matrix in matter and the λᵢ are the matter eigenvalues. The transverse plane is the complex
          plane; the colored curves are the running partial sums, and the white curve is the resultant
          amplitude A. At the marker, the arrows show the head-to-tail phasor sum at that point of the
          sweep. The floor curve is <span class="k">|A|² = P(νμ→νe)</span>.
        </p>
        <p>
          How to read it: oscillation is literally interference here — the probability peaks where the
          three arms align and dips where they cancel. δCP rotates the relative phases of the arms
          (changing where they align), while E and ρ change the arm lengths through the matter mixing.
        </p>
        <p><button class="linkish" @click="openView('phasors')">open this view →</button></p>
      </section>

      <section id="faq-tube">
        <h3>Flavor tube</h3>
        <p>
          The flavor content of an initial νμ along the baseline. Each cross-section of the tube is a pie
          of the three flavor probabilities P(νμ→νe/μ/τ) at that L; the radius is constant because the
          fractions always sum to 1 (unitarity). The disk rides the marker and shows the instantaneous
          composition. The “stacked bands” mode is the unrolled version — a stacked area plot extruded
          into 3D, whose flat top edge is the same unitarity check.
        </p>
        <p>
          How to read it: the dominant blue (νμ) sector thins where the oscillation transfers probability
          to green (ντ, the main partner at atmospheric frequencies) with a small red (νe) appearance
          wedge — the wedge an appearance experiment actually measures. The companion plot is the matching
          2D stacked area chart.
        </p>
        <p><button class="linkish" @click="openView('tube')">open this view →</button></p>
      </section>
    </article>
  </div>
</template>

<style scoped>
.faq-page {
  flex: 1;
  min-width: 0;
  overflow-y: auto;
  background: var(--bg);
}
.doc {
  max-width: 820px;
  margin: 0 auto;
  padding: 34px 28px 80px;
  font-size: 16px;
  line-height: 1.7;
}
h2 { font-size: 26px; margin: 0 0 6px; }
.sub { color: var(--muted); margin: 0 0 10px; }
section {
  margin-top: 26px;
  padding-top: 14px;
  border-top: 1px solid var(--border);
}
h3 {
  font-family: var(--font-mono);
  font-size: 17px;
  color: var(--accent);
  margin: 0 0 10px;
}
p { margin: 10px 0; }
.k {
  font-family: var(--font-mono);
  font-size: 14.5px;
  background: var(--surface-2);
  border-radius: 4px;
  padding: 0 5px;
  white-space: nowrap;
}
.linkish { font-size: 14px; }
a { color: var(--accent); }
</style>
