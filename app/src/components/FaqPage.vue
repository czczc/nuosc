<script setup>
import { onMounted, ref } from 'vue';
import { store } from '../store.js';
import { router } from '../router.js';

const BIBTEX = `@misc{nuglass,
  author       = {Zhang, Chao},
  title        = {NuGlass: a fast and interactive 3D visualization of neutrino oscillations},
  year         = {2026},
  howpublished = {\\url{https://czczc.github.io/nuglass/}},
  note         = {Source code: \\url{https://github.com/czczc/nuglass}}
}`;

const copied = ref(false);

function copyBibtex() {
  navigator.clipboard.writeText(BIBTEX);
  copied.value = true;
  setTimeout(() => { copied.value = false; }, 1500);
}

function openView(id) {
  router.push('/' + id);
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
        NuGlass is a fast, interactive, and animated 3D visualization of neutrino oscillations: the strange way neutrinos change
        from one type to another as they travel. Six 3D views are provided, each paired with live 2D plots.
        This page explains what each view shows and how the numbers behind it are computed.
      </p>

      <section id="faq-engines">
        <h3>How the numbers are computed</h3>
        <p>
          All calculations are performed live in your web browser from the equations of three-flavor neutrino oscillations, including the extra effect neutrinos feel when they travel
          through the Earth (the “matter effect”). The default parameter values come from the
          world combined fit of all oscillation experiments
          (<a href="http://www.nu-fit.org" target="_blank" rel="noopener">NuFit 6.1, 2025</a>), and can be adjusted under “all parameters”. Antineutrinos are computed the same way with the relevant signs flipped.
        </p>
        <p>
          Two calculators share the work. <strong>NuFast-LBL</strong>
          (<a href="https://arxiv.org/abs/2405.02400" target="_blank" rel="noopener">arXiv:2405.02400</a>,
          <a href="https://github.com/PeterDenton/NuFast-LBL" target="_blank" rel="noopener">reference code</a>)
          is a very fast method by Denton &amp; Parke, quick enough to recompute a whole surface of
          ~26,000 points every time you drag a slider. It powers the oscillogram, biprobability and
          worldline views.
          The <strong>exact amplitude engine</strong> tracks the full quantum state of the neutrino, not
          just the final probabilities; the statesphere and phasors views need that extra information.
          The two calculators are cross-checked against each other to better than one part in ten million.
        </p>
        <p>
          <strong>Experiment presets</strong>: DUNE, NOvA and T2K each set the travel
          distance, rock density and beam energy range of a real accelerator experiment; the reactor
          experiments JUNO, KamLAND and Daya Bay (offered for the νe→νe channel) do the same at MeV
          energies. Picking one sets every parameter to that experiment’s values.
        </p>
      </section>

      <section id="faq-oscillogram">
        <h3>Oscillogram</h3>
        <p>
          A landscape of probability: the surface’s height and color show the chance that a neutrino flavor has turned into another flavor (or itself) by the time it is detected. Left-to-right is the
          neutrino’s energy E; the depth axis is your choice of a second variable: the travel distance L,
          the CP phase δCP (the parameter that can make neutrinos and antineutrinos behave differently),
          or the density ρ of the rock the beam passes through.
        </p>
        <p>
          The white line is a slice through the landscape at the current slider setting, and it is
          exactly the curve drawn in the first 2D panel below (probability vs energy, for neutrinos and
          antineutrinos, with the experiment’s energy peak dashed). The second 2D panel shows how the
          probability at the current energy changes as δCP goes around its full circle.
        </p>
        <p>
          Press ▶ on the “animate” row to set the picture in motion, sweeping the distance, the energy,
          or δCP; the matching slider in the controls panel moves along with it, so you can pause
          anywhere and continue by hand. The front / top / side buttons swing the camera to flat views:
          “top” looks straight down and turns the surface into the classic color-map oscillogram that
          experiments publish.
        </p>
        <p><button class="linkish" @click="openView('oscillogram')">open this view →</button></p>
      </section>

      <section id="faq-biprob">
        <h3>Biprobability</h3>
        <p>
          Two chances plotted against each other: the probability that a neutrino becomes a
          different flavor (horizontal) versus the same probability for antineutrinos (vertical), at
          one energy and distance. As δCP sweeps through its full circle, the point traces out an
          ellipse. The 3D view stacks those ellipses along energy, for the two possible orderings of the
          neutrino masses (orange = normal, blue = inverted); the small spheres mark the current δCP,
          and the connecting line shows the ordering difference.
          (After Minakata &amp; Nunokawa,
          <a href="https://arxiv.org/abs/hep-ph/0108085" target="_blank" rel="noopener">arXiv:hep-ph/0108085</a>.)
        </p>
        <p>
          How to read it: if neutrinos and antineutrinos behaved identically, the point would sit on the
          gray diagonal; its distance from that diagonal is CP violation. The gap between
          the orange and blue rings is what lets experiments tell the two mass orderings apart; where
          the rings overlap, that energy alone cannot decide.
        </p>
        <p><button class="linkish" @click="openView('biprob')">open this view →</button></p>
      </section>

      <section id="faq-sphere">
        <h3>Statesphere</h3>
        <p>
          The neutrino’s quantum state drawn as an arrow inside a globe. The north pole means “certainly
          an electron neutrino”, the south pole “certainly a muon neutrino”, and everywhere in between
          is a quantum mixture of the two. As the neutrino travels, the arrow swings around and traces
          the colored path; press ▶ to watch it move.
        </p>
        <p>
          How to read it: a textbook two-type oscillation would keep the arrow’s tip on the globe’s
          surface. Here the tip dips inside: a shorter arrow means part of the probability has leaked
          into the third type, which this picture cannot show directly. The 2D panel
          plots all three probabilities along the sweep, with its marker synced to the 3D arrow.
        </p>
        <p>
          For readers who want the exact construction: this globe is the
          <a href="https://en.wikipedia.org/wiki/Bloch_sphere" target="_blank" rel="noopener">Bloch
          sphere</a> of a quantum state, applied to the two flavors on the poles. Writing
          a<sub>e</sub> and a<sub>μ</sub> for the quantum amplitudes of the flavor states νe and νμ, the arrow is
          <span class="k">b = ( 2&hairsp;Re(a<sub>e</sub>a<sub>μ</sub>*), 2&hairsp;Im(a<sub>e</sub>a<sub>μ</sub>*), |a<sub>e</sub>|² − |a<sub>μ</sub>|² )</span>. The height toward the north pole is the probability difference
          P<sub>e</sub> − P<sub>μ</sub>; the two horizontal components are the real and imaginary parts
          of the interference term (the “coherence”) between the flavors. |b| = 1 means a pure two-flavor state;
          anything shorter means probability sits in the third flavor. For the ντ-vs-νμ pole choice,
          replace a<sub>e</sub> with a<sub>τ</sub>.
        </p>
        <p><button class="linkish" @click="openView('sphere')">open this view →</button></p>
      </section>

      <section id="faq-phasors">
        <h3>Phasors</h3>
        <p>
          Why oscillation happens, drawn as directly as possible. Quantum mechanics adds up the three
          ways a neutrino can become another neutrino as three arrows in a plane: one for each
          neutrino mass, each rotating at its own speed as the neutrino travels. The colored curves show
          the arrows chained head to tail, the white curve is their total amplitude, and the curve on the floor is
          the resulting probability (amplitude squared).
        </p>
        <p>
          How to read it: the probability peaks where the arrows line up and vanishes where they cancel.
          Oscillation is interference, the same phenomenon as overlapping ripples on a pond. Changing
          δCP rotates the arrows relative to each other (moving where they align), while energy and rock
          density change the arrows’ lengths.
        </p>
        <p><button class="linkish" @click="openView('phasors')">open this view →</button></p>
      </section>

      <section id="faq-tube">
        <h3>Flavortube</h3>
        <p>
          Follow a neutrino along its journey and watch its identity mix. Each slice of the tube is
          a pie chart of the three possibilities at that distance: still a muon neutrino (blue), turned
          into a tau neutrino (green), or turned into an electron neutrino (red). The tube’s cross-section area
          never changes, because the three chances always add up to 100%. The disk rides the animation
          marker and shows the mix at one spot. The tube also sits on an energy axis: drag the E slider
          (or animate over E) and the whole tube glides along it while its pattern stretches; lower
          energies flip flavor faster.
          Switch the display to “stacked bands” and the same information unrolls into a box.
        </p>
        <p><button class="linkish" @click="openView('tube')">open this view →</button></p>
      </section>

      <section id="faq-loe">
        <h3>Worldline</h3>
        <p>
          The other views sweep distance or energy separately; this one uses the variable oscillation
          actually depends on: L/E (distance traveled divided by energy) is proportional to the proper time the neutrino itself
          experiences, so this axis is the neutrino's own worldline, with every experiment pinned
          somewhere along it. In vacuum, every oscillation probability is a function of the single
          combination L/E, so each channel collapses onto one
          universal curve. Both oscillation frequencies
          are visible at once: the fast atmospheric wiggle (first dip near L/E ≈ 500 km/GeV) and the
          slow solar valley (near L/E ≈ 15,000 km/GeV) with the fast wiggle riding on top of it.
          The thick
          curve is the channel selected in the header; the checkboxes overlay the others.
        </p>
        <p>
          The third axis shows what breaks this elegant collapse: matter. Traveling through rock adds
          an effect that grows with energy E itself, not with L/E, so two neutrinos with the same L/E
          but different energies no longer oscillate identically. The colored surface (drawn like the
          Oscillogram view) spreads the selected channel across the experiment's energy window, computed with the matter effect at the
          current density ρ: at ρ = 0 the surface is perfectly flat along the energy axis (the front
          view reproduces the vacuum curve). The line riding the surface is the cut through the
          experiment's energy window at the fixed baseline L, and it is the same spectrum drawn in the second 2D panel.
        </p>
        <p><button class="linkish" @click="openView('loe')">open this view →</button></p>
      </section>

      <section id="faq-experiments">
        <h3>Custom experiments</h3>
        <p>
          The experiment chips in the top bar are presets: each sets the distance to the
          detector, the rock density along the way, and the source's energy window. The chip row
          follows the oscillation channel selected next to the view tabs — beam channels
          (νμ→νe, νμ→νμ, νμ→ντ) offer DUNE / NOvA / T2K, while νe→νe offers the reactor
          experiments JUNO / KamLAND / Daya Bay (which also switch the particle toggle to
          antineutrinos, since reactors emit ν̄e). The <span class="k">+</span>
          chip opens a page where you can build your own: name it, set the baseline and beam, declare
          which channels it measures and whether the source is antineutrinos, and, if you want to
          explore, override the oscillation parameters themselves. The inputs switch between GeV·km
          and MeV·m for reactor-scale setups, and the form can be reset to DUNE or JUNO defaults.
          Saved experiments live in your browser's localStorage (nothing is uploaded), which persisits across visits. Loading one makes it the active experiment everywhere.
        </p>
      </section>

      <section id="faq-cite">
        <h3>How to cite</h3>
        <p>
          If you use NuGlass in a paper or talk, please cite it as:
        </p>
        <pre class="bib">{{ BIBTEX }}</pre>
        <p><button class="linkish" @click="copyBibtex">{{ copied ? 'copied ✓' : 'copy BibTeX' }}</button></p>
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
.bib {
  font-family: var(--font-mono);
  font-size: 13px;
  line-height: 1.5;
  background: var(--surface-2);
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 12px 14px;
  overflow-x: auto;
}
a { color: var(--accent); }
</style>
