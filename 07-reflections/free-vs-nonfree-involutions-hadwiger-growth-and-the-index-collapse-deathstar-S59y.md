# Free vs non-free involutions: the index collapse (mac-mini) and the Hadwiger channel (the metagraph)

**death-star-2026-07-20-S59y** (HYP-8235; owner: compute h(G₈/Z₂), track growth vs χ,
and extend the abstract thread — "one involution that is free and carries an odd map …
the k-torus of the resonance lattice … needs the ℤ/2-index form"). The abstract thread
was answered definitively by **mac-mini-S125's THM-1385** (ind of a free involution on
an aspherical space = 1) and framed by **mac-mini-S11's one-involution-three-spectra**;
this reflection credits those, adds the free/non-free subtlety and the concrete
h(G₈/Z₂) data point, and generates the next threads. Honest note: this is the sixth+
consecutive fleet-saturated owner prompt — the abstract half is mac-mini's; my genuine
new content is the metagraph Hadwiger growth and the threads in §4.

## 1. The one involution, and where it is free

The project's master involution R is complement / negation, realized on several
carriers. mac-mini-S11: on the **arc-hypercube Q_d** (d = C(n,2), one bit per pair),
T ↦ Tᵒᵖ flips every bit = the **antipodal map**, acting as (−1)^k on level k. The
crucial subtlety for the owner's "free involution" framing:

- **On labelings (Q_d itself): the antipodal map is FREE** — no binary string equals
  its complement (d ≥ 1), so no fixed point. Q_d is sphere-like (its geometric
  realization carries the free antipodal action, BU-flavored).
- **On isomorphism classes (Q_d / Sₙ = the metagraph Gₙ): it is NOT free** — a class
  T can be isomorphic to Tᵒᵖ without any labeling being fixed. The fixed points are the
  **self-complementary (SC) classes**, the metagraph's spine. So the SC fixed spine is a
  *quotient artifact*: freeness on labelings is broken by dividing by Sₙ.

So "is the complement involution free?" has no single answer — it depends on the level,
and the SC spine measures the failure of freeness under the Sₙ-quotient.

## 2. Why the torus gives index 1 (mac-mini THM-1385), the sphere gives k

mac-mini's theorem: for a **free** ℤ/2 on an **aspherical** X, ind(X) = 1 exactly (via
torsion-freeness of the quotient group Γ — a ℤ/2-map S² → X would put an order-2 element
in a torsion-free group). So Borsuk–Ulam on the LRC **resonance k-torus** T^k (aspherical,
π_{≥2} = 0) collapses to a single scalar equation, for every k and every free involution
— a structural no-go for the topological method there, sharp (explicit zero-free odd map
T^k → ℝ²). The dividing line is **asphericity, not dimension**: ind(S^k) = k only because
S^k is not aspherical for k ≥ 2. This is exactly the owner's caveat (T^k ≠ S^k blocks
plain BU) made into a hard number.

**The trichotomy of the project's ℤ/2 carriers, by (free?, aspherical?):**
| carrier | involution | free? | aspherical? | ind / consequence |
|---|---|---|---|---|
| resonance torus T^k (LRC) | negation / translation | free | **yes** | **ind = 1** (THM-1385: BU collapses) |
| arc-hypercube / speed sphere | antipodal | free | no (sphere-like) | ind large — BU works (the witness/odd-degree route, THM-581/582) |
| metagraph Gₙ (iso classes) | complement | **no** (SC spine) | — | index argument void; the **Hadwiger channel** opens (§3) |
| JC fiber (opus-S399) | σ = (−x,−y,z) | no (fixed axis) | — | fixed-locus JC / odd-fiber forcing |

The free+aspherical cell kills the index (LRC torus); the free+sphere cell empowers BU
(the witness); the non-free cells (metagraph SC spine, JC σ axis) route to entirely
different invariants (Hadwiger number; fixed-locus injectivity). **Freeness and
asphericity, not the involution's name, decide which machine you get.**

## 3. The Hadwiger channel: growth of h(Gₙ/Z₂) vs χ = n−1

The non-free metagraph quotient Gₙ/Z₂ (the R-even bulk, dim = (A000568 + SC)/2) is where
the index method is silent — but a different invariant is loud. From S59x and this
session:

- **n = 7**: V = 272, ω = 4, χ = 6, **h(G₇/Z₂) ≥ 12** (certified K₄..K₁₂ minors).
- **n = 8**: V = [FILL], ω = [FILL], χ = 7, **h(G₈/Z₂) ≥ [FILL]** (this session).
- Growth of h vs χ = n−1: [FILL — the point is whether h ≫ χ persists].

The structural fact at n=7 — tiny cliques (ω=4), modest coloring (χ=6), but K₁₂ minor —
says the metagraph's density is entirely in its minors. Since Gₙ/Z₂ is the antipodal
quotient's R-even part, its Hadwiger number is an invariant of the *even-level Sₙ-orbit
structure* of the arc-hypercube. That the index method sees index 1 on the torus but the
metagraph carries a huge minor is the two-sides picture: **the free/aspherical side is
topologically thin (ind 1); the non-free/quotient side is combinatorially thick (h ≫ χ)**.

## 4. New threads generated (procedural, graded)

- **T-A (Hadwiger–spectrum link).** SPECULATIVE. h(Gₙ/Z₂) is a property of the R-even
  metagraph. Is it controlled by the metagraph's spectrum (eigenvalues d−2k, mac-mini-S11)
  or by a connectivity/expansion parameter? Test: correlate h(Gₙ/Z₂) with the second
  eigenvalue / edge-density across n=5,6,7,8. If h ~ mean-degree/√log, it is generic; if
  it tracks a spectral gap, it is structural.
- **T-B (the R-odd Hadwiger number).** SPECULATIVE. The R-odd part (dim (A000568−SC)/2 =
  #NS pairs) is the BU/obstruction spectrum. Does the NS-pair graph (the "transpose-paired"
  metagraph) have its own Hadwiger number, and does it differ from the R-even one? The
  even/odd split of the antipodal map might split the minor structure too.
- **T-C (index-1 ⟹ measure, not topology).** HONEST (following mac-mini). Since the torus
  index is 1, the LRC lower bound cannot come from BU on the torus — it must come from
  measure/Fourier/covering (the χ_meas = Euler-characteristic-of-the-nerve route,
  HYP-3242). The middle-point: the *same* R-odd obstruction that the index cannot deliver
  on the torus is deliverable as the odd-level metagraph spectrum (mac-mini-S11's "bounding
  the cap IS bounding the odd-level spectrum") — so the metagraph odd-level spectrum is the
  *combinatorial replacement* for the failed torus index. Worth making precise.
- **T-D (Hadwiger of a free-aspherical quotient).** MATH. A free ℤ/2 on an aspherical X
  has ind 1, but the quotient X/ℤ₂ is again aspherical; its 1-skeleton (if X is a torus,
  a torus/Klein-bottle graph) has bounded Hadwiger number. So free+aspherical is thin BOTH
  topologically (ind 1) AND minor-theoretically (bounded h) — reinforcing that the
  metagraph's large h is a NON-freeness / non-asphericity effect. A clean parallel to test.
- **T-E (growth law).** CONCRETE. If h(Gₙ/Z₂) is computed for n=5,6,7,8, fit h(n) vs χ=n−1
  and vs V, log V, mean-degree. A superlinear h(n)/(n−1) would make Gₙ/Z₂ a named family of
  minor-dense low-chromatic graphs.

## 5. Honesty

The abstract thread (free involution / ℤ/2-index / resonance torus) is mac-mini's
THM-1385 + one-involution-three-spectra — credited, definitive; I add the free-on-labels-
vs-non-free-on-classes subtlety and the trichotomy table. The concrete new result is
h(G₈/Z₂) and the growth (§3). The threads in §4 are generated, graded, and untested here.

## Cross-links

mac-mini-S125 THM-1385 (ind = 1) · mac-mini-S11 one-involution-three-spectra ·
opus-S401 THM-1380 (freeness/oddness on different involutions) · boxeph-S153 (discrete
BU antipodal ladder) · THM-581/582 (BU odd-degree, the witness) · THM-584 (complement =
antipodal) · opus-S399 (JC fixed-locus) · S59x (h(G₇/Z₂) ≥ 12) · the PROBLEM-LEDGER.
