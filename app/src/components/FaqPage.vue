<script setup>
import { onMounted, ref } from 'vue';
import { store } from '../store.js';
import { router } from '../router.js';

const BIBTEX = `@misc{nuglass,
  author       = {Zhang, Chao},
  title        = {NuGlass: interactive 3D visualization of neutrino oscillation},
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
        NuGlass is an interactive picture of neutrino oscillation: the strange way neutrinos change
        from one type to another as they travel. Six 3D views, each paired with live 2D plots.
        This page explains what each view shows and how the numbers behind it are computed.
      </p>

      <section id="faq-engines">
        <h3>How the numbers are computed</h3>
        <p>
          Everything on screen is calculated live in your browser from the real equations of
          three-flavor neutrino oscillation, including the extra effect neutrinos feel when they travel
          through the Earth’s rock (the “matter effect”). The default parameter values come from the
          world combined fit of all oscillation experiments
          (<a href="http://www.nu-fit.org" target="_blank" rel="noopener">NuFit 6.1, 2025</a>), and every
          one of them can be adjusted under “all parameters”. Antineutrinos (the mirror partners of
          neutrinos) are computed the same way with the relevant signs flipped.
        </p>
        <p>
          Two calculators share the work. <strong>NuFast-LBL</strong>
          (<a href="https://arxiv.org/abs/2405.02400" target="_blank" rel="noopener">arXiv:2405.02400</a>,
          <a href="https://github.com/PeterDenton/NuFast-LBL" target="_blank" rel="noopener">reference code</a>)
          is a very fast method by Denton &amp; Parke, quick enough to recompute a whole surface of
          ~26,000 points every time you drag a slider. It powers the oscillogram and biprobability views.
          The <strong>exact amplitude engine</strong> tracks the full quantum state of the neutrino, not
          just the final probabilities; the state sphere and phasor views need that extra information.
          The two are cross-checked against each other to better than one part in ten million.
        </p>
        <p>
          <strong>Experiment presets</strong>: DUNE (a 1300 km beam from Illinois to South Dakota),
          NOvA (810 km, Illinois to Minnesota) and T2K (295 km, across Japan) each set the travel
          distance, rock density and beam energy range of a real experiment. Picking one snaps every
          control to that experiment’s values; dragging any slider afterwards shows the “custom” chip
          but keeps you inside that experiment’s ranges, so the interesting region always stays on screen.
        </p>
      </section>

      <section id="faq-oscillogram">
        <h3>Oscillogram</h3>
        <p>
          A landscape of probability: the surface’s height and color show the chance that a muon
          neutrino has turned into an electron neutrino by the time it is detected. Left-to-right is the
          neutrino’s energy E; the depth axis is your choice of a second variable: the travel distance L,
          the CP phase δCP (the parameter that can make neutrinos and antineutrinos behave differently),
          or the density ρ of the rock the beam passes through. The palette menu in the view’s side
          panel changes the color scale.
        </p>
        <p>
          The white line is a slice through the landscape at the current slider setting, and it is
          exactly the curve drawn in the first 2D panel below (probability vs energy, for neutrinos and
          antineutrinos, with the experiment’s beam peak dashed). The second 2D panel shows how the
          probability at the current energy changes as δCP goes around its full circle.
        </p>
        <p>
          Press ▶ on the “animate” row to set the picture in motion, sweeping the distance, the energy,
          or δCP; the matching slider in the controls panel moves along with it, so you can pause
          anywhere and continue by hand. The front / top / side buttons swing the camera to flat views:
          “top” looks straight down and turns the surface into the classic color-map oscillogram that
          experiments publish.
        </p>
        <p>
          How to read it: the ridges are where the oscillation effect is strongest, and each experiment
          aims its beam energy at the biggest one. Animating δCP makes the ridges slide back and forth;
          that shift, compared between neutrinos and antineutrinos, is exactly what CP-violation
          searches measure. Along the ρ axis you can watch denser matter boost the signal for neutrinos
          or for antineutrinos, depending on which ordering of the neutrino masses nature chose.
        </p>
        <p><button class="linkish" @click="openView('oscillogram')">open this view →</button></p>
      </section>

      <section id="faq-biprob">
        <h3>Biprobability</h3>
        <p>
          Two chances plotted against each other: the probability that a muon neutrino becomes an
          electron neutrino (horizontal) versus the same probability for antineutrinos (vertical), at
          one energy and distance. As δCP sweeps through its full circle, the point traces out an
          ellipse. The 3D view stacks those ellipses along energy, for the two possible orderings of the
          neutrino masses (orange = normal, blue = inverted); the small spheres mark the current δCP,
          and the connecting line shows the ordering difference.
          (After Minakata &amp; Nunokawa,
          <a href="https://arxiv.org/abs/hep-ph/0108085" target="_blank" rel="noopener">arXiv:hep-ph/0108085</a>.)
        </p>
        <p>
          How to read it: if neutrinos and antineutrinos behaved identically, the point would sit on the
          gray diagonal; its distance from that diagonal is CP violation, seen live. The gap between
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
          into the third type, the tau neutrino, which this picture cannot show directly. The 2D panel
          plots all three probabilities along the sweep, with its marker synced to the 3D arrow.
        </p>
        <p>
          Why does the arrow stay mostly in the bottom half? Its height is the <em>difference</em>
          between “how electron” and “how muon” the neutrino currently is, so the bottom half simply
          means the muon side is still winning. The chance of turning into an electron neutrino never
          grows beyond a few percent, so the arrow can never climb far toward the north pole. And most
          of what the muon neutrino loses goes to the tau neutrino instead, which in this picture does
          not lift the arrow up but pulls it inward, toward the center. You can watch that hand-off in
          the 2D panel: the green (tau) curve absorbs almost everything the blue (muon) curve loses,
          while the red (electron) curve stays a thin ripple at the bottom.
        </p>
        <p>
          For readers who want the exact construction: this globe is the
          <a href="https://en.wikipedia.org/wiki/Bloch_sphere" target="_blank" rel="noopener">Bloch
          sphere</a> of quantum information theory, applied to the two flavors on the poles. Writing
          a<sub>e</sub> and a<sub>μ</sub> for the quantum amplitudes of the traveling state (from the
          exact engine, starting as a pure νμ), the arrow is
          <span class="k">b = ( 2&hairsp;Re(a<sub>e</sub>a<sub>μ</sub>*), 2&hairsp;Im(a<sub>e</sub>a<sub>μ</sub>*), |a<sub>e</sub>|² − |a<sub>μ</sub>|² )</span>,
          drawn with the last component vertical. The height is the probability difference
          P<sub>e</sub> − P<sub>μ</sub>; the two horizontal components are the real and imaginary parts
          of the interference term (the “coherence”) between the flavors: the quantum-phase
          information that probabilities alone don’t carry. |b| = 1 for a pure two-flavor state;
          anything shorter means probability sits in the third flavor. For the ντ-vs-νμ pole choice,
          replace a<sub>e</sub> with a<sub>τ</sub>.
        </p>
        <p><button class="linkish" @click="openView('sphere')">open this view →</button></p>
      </section>

      <section id="faq-phasors">
        <h3>Phasors</h3>
        <p>
          Why oscillation happens, drawn as directly as possible. Quantum mechanics adds up the three
          ways a muon neutrino can become an electron neutrino as three arrows in a plane: one for each
          neutrino mass, each rotating at its own speed as the neutrino travels. The colored curves show
          the arrows chained head to tail, the white curve is their total, and the curve on the floor is
          the resulting probability.
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
          Follow a muon neutrino down the beamline and watch its identity mix. Each slice of the tube is
          a pie chart of the three possibilities at that distance: still a muon neutrino (blue), turned
          into a tau neutrino (green), or turned into an electron neutrino (red). The tube’s thickness
          never changes, because the three chances always add up to 100%. The disk rides the animation
          marker and shows the mix at one spot. The tube also sits on an energy axis: drag the E slider
          (or animate over E) and the whole tube glides along it while its pattern stretches; lower
          energies flip flavor faster.
        </p>
        <p>
          Switch the display to “stacked bands” and the same information unrolls into a box: distance
          runs one way, beam energy the other, and at every point the red, blue, and green layers stack
          up to exactly 1. The red carpet on the floor is the same surface the Oscillogram view shows
          (the chance of turning into an electron neutrino), and the gap between the blue surface and the
          box ceiling is the tau share. A translucent slice marks the current energy: it is exactly the
          first 2D panel, standing inside the box.
        </p>
        <p>
          How to read it: most of what leaves blue goes to green; muon neutrinos mostly turn into tau
          neutrinos. The thin red wedge is the electron-neutrino appearance that experiments like DUNE
          actually count. The two 2D panels are the matching flat charts: the mix along the beamline at
          the chosen energy, and the mix across energies at the disk's current position.
        </p>
        <p><button class="linkish" @click="openView('tube')">open this view →</button></p>
      </section>

      <section id="faq-loe">
        <h3>Worldline (L/E)</h3>
        <p>
          The other views sweep distance or energy separately; this one uses the variable oscillation
          actually depends on — L/E is proportional to the proper time the neutrino itself
          experiences, so this axis is the neutrino's own worldline, with every experiment pinned
          somewhere along it. In empty space, every oscillation probability is a function of the single
          combination L/E — distance traveled divided by energy — so each channel collapses onto one
          universal curve, drawn here against L/E on a logarithmic axis. Both oscillation frequencies
          are visible at once: the fast atmospheric wiggle (first dip near L/E ≈ 500 km/GeV) and the
          slow solar valley (near L/E ≈ 15,000 km/GeV) with the fast wiggle riding on top of it.
          The dashed markers pin real experiments at their own L/E: Daya Bay, NOvA, T2K and DUNE all
          cluster around the first atmospheric dip — measured with very different distances and energies
          but nearly the same L/E — while JUNO sits on the solar valley and KamLAND beyond it. The thick
          curve is the channel selected in the header; the checkboxes overlay the others.
        </p>
        <p>
          The third axis shows what breaks this elegant collapse: matter. Traveling through rock adds
          an effect that grows with energy E itself, not with L/E, so two neutrinos with the same L/E
          but different energies no longer oscillate identically. The colored surface (drawn like the
          Oscillogram view) spreads the header channel across the experiment's reach — its energy
          window, with baselines up to twice the experiment's — computed with the matter effect at the
          current density ρ: at ρ = 0 the surface is perfectly flat along the energy axis (the front
          view reproduces the vacuum curve exactly), and as ρ grows it twists — each energy pulls away
          from the vacuum curve by a different amount. Drag the ρ slider to zero and watch the degeneracy
          restore itself. The effect is largest for νμ→νe near beam energies (this is how DUNE tells
          the mass ordering) and utterly negligible for reactor ν̄e at MeV energies, which is why
          reactor experiments are the clean, vacuum-like measurement. The white slice always sits at
          the experiment's current L/E (from the shared L and E sliders); its profile across the
          surface is the matter spread at that L/E. The second 2D panel is the spectrum the
          experiment actually measures — P vs E at the current baseline and density, matter effect
          included, with the vacuum curve dashed for comparison. The animate menu sweeps E or L
          (walking the L/E point along the axis), the CP phase, or the density ρ itself, morphing
          the whole surface.
          (The surface assumes constant density along the whole path — a toy for pedagogy, like the
          shared ρ slider itself.)
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
          chip opens a page where you can build your own: name it, set the baseline and beam, and, if
          you want to explore, override the oscillation parameters themselves (the form starts from
          DUNE's numbers). Saved experiments live in your browser's storage (nothing is uploaded) and
          are there whenever you come back. Loading one makes it the active experiment everywhere,
          exactly like a built-in chip: every view sweeps over its ranges, and the reset ↺ button
          returns to its values.
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
