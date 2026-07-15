---
source: klein-2026-07-15-S313 (cont.4) — LRC(14) session applying the S313 lenses (truncation
  grammar, corona onset, Kakeya tightness, rigidity torsors, Milgram); synthesizes THM-604,
  THM-663, THM-863, THM-867/868/869/870, s703, opus-S316's truncation-grammar thread
status: CORE PROVED (LEM-020's six identities + mid-band law, 8/8 exact checks); the synthesis
  is the content; one honest reversal recorded
tags: [lrc14, coverage-spectrum, bonferroni, corona, kakeya, mid-band, roots-of-unity,
  truncation-grammar, second-moment, rigidity]
---

# The coverage spectrum: one grammar, four instruments

**The object.** Every question in the LRC(14) covering case factors through ONE object I had
not seen named: the multiplicity spectrum μ_k(x) of the thirteen 1/7-arcs. Uncovered mass,
every Bonferroni level, every d-fold resonance sum, the pair-overlap vein — all are LINEAR
FUNCTIONALS of μ (LEM-020). The depth ladder does not probe thirteen speeds; it probes one
14-entry vector of masses.

**The grammar, fourth instrument.** S313 built the truncation grammar three times: Moser/D
(polygonal shadows of Pascal), the Landau corona (shadow of the E8 shells), Bonferroni was
conjectured to be the fourth. It is, verbatim: B_D = the rank-D truncation, deficit = the
binomial tail Σ C(k−1,D)μ_k, EXACT iff max multiplicity ≤ D, and the first failure is the
cheapest (D+1)-fold stack — the corona onset. Landau's onset was a duel of two champions;
Bonferroni's is a champion STACK (a (D+1)-fold resonance); THM-863's ρ-minimizing pair (1,12)
is the d = 2 onset. One law: **truncations fail first at the cheapest deep coincidence, and
the coincidence is always a small cast of extreme characters.**

**Kakeya tightness transfers exactly.** K(A₅) = 15 was proved by "needles pairwise meet ≤ 1
⟹ Bonferroni-2 tight" (THM-870); the optimal witness has every point in ≤ 2 needles. The
LRC twin: a covering x with max multiplicity 2 is FORCED to the spectrum (μ₁, μ₂) = (1/7, 6/7)
— and the tight AP at its clock witnesses x = 1/7, 1/13, 1/14 realizes exactly this. The
extremal instance that killed five fleet routes is, in spectrum language, just "the unique
minimal-multiplicity double cover" — the same object as the icosahedral K₆ witness. Tightness
everywhere = minimal multiplicity + saturated budget (the 6/7 slack identity; note the SAME
6/7 is the base of THM-604's independent limit (6/7)^13 — budget and baseline coincide).

**An honest reversal (recorded because the correction is the insight).** I predicted deep
multiplicity concentrates at low-order roots of unity (worry-set intuition). Backwards: for
the AP, low-order x IS the covering band (max multiplicity 2, perfectly spread), and the
13-deep stacks live at HIGH-order/small x — the single-cluster far tail. The worry set is not
where arcs pile up; it is where they interlock perfectly. Resonance = rigidity = perfect
spreading, not clumping. (The TV-distance of μ to Binom(13, 1/7) is a clean resonance meter:
0.565 at x = 1/7 vs 0.284 at a golden-ratio-like x.) This retro-explains the route
architecture: the mid-band is where truncation certificates are EXACT (finite checks belong
there); the far tail is where the spectrum degenerates to one stacked block (explicit
transfers belong there). THM-710 + Lemma A + THM-667 were forced moves, spectrally.

**Roots of unity, properly.** The covering x's found are clock fractions (q = 7, 13, 14) —
the cyclotomic bottom of s703's solvability inversion. In spectrum language the inversion
reads: *solvable/cyclotomic x's have rigid, low-entropy spectra pinned by the sum rules;
generic x's have near-binomial spectra with slack.* Tightness is a rigidity phenomenon of μ,
and the Milgram lens says what kind: μ at the clock witnesses is the discrete measure of a
perfect 2-cover of Z/7-type — a discriminant-form-flavored object (the 14 = 2·7 clock), which
is where any deeper Gauss-sum structure should be sought (OPEN).

**What this buys next (handoffs).** (1) The maxC = 2 covering-x characterization — if the
clock witnesses are the ONLY minimal-multiplicity covers, the covering case's tight locus is
spectrally enumerable. (2) The d = 3 pincer: c₃ = 2/105 (THM-863) + S₃ = ΣC(k,3)μ_k squeezes
covering spectra between moment constraints — a potential new uniform floor with no new
analysis, just linear algebra on μ. (3) Lean: the whole frame is exact rational identities,
decide-shaped, and slots beside LRCDiscreteBonferroni.
