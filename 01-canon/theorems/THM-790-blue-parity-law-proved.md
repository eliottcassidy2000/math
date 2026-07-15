---
id: THM-790
title: THE LEG LAW and the blue parity law PROVED — every d=m line's axis drop is Δx = 8(e₁ − e_n) (the difference of the two LEG out-degrees); the full flip acts on centered scores as reversal-plus-leg-defect; closed-form drop spectra for the whole layer AND the blue sublayer at every n; the parity law is a one-line corollary; the half-tiling model counts blue objects
status: PROVED (every identity verified on ALL tilings n = 4..7 — the leg law, the score action, both generating functions, all six lemmas) + n=8 CERTIFIED (subH-augmented classifier: 6,880/6,880 classes; every prediction exact incl. the blue GF histogram {0:640, 16:960, 32:384, 48:64}; codex-S9's interim 6,874-bucket caveat is superseded by the completed run)
source: opus-2026-07-14-S305 (owner directive: prove the blue parity law, check n=8, find predictive formulas and recursive structure via the half/quarter tiling models)
depends_on:
  - THM-787   # the flow study whose laws this proves
related: [THM-785, THM-787, THM-791, HYP-6855, everything-is-the-triangle, the geometric-alignment frame]
verification: 04-computation/blue_parity_law_proof_referee_opus_S305.py (six lemmas, all blue tilings n=4..7)
  + 04-computation/leg_law_referee_opus_S305.py (the leg law + score action + GFs, ALL tilings n=4..7)
  + 04-computation/metagraph_flow_n8_check_opus_S305.py (n=8)
---

# THM-790 — the leg law, and the blue parity law as its corollary

**Setup.** Base path n → … → 1; tiles (x,y), x−y ≥ 2, bit 1 = arc x→y; the
d=m line joins t to its full flip t̄; axis x = Σ_v d_v², d_v = 2s_v − (n−1).
Let **e₁(t)** = the out-tile-degree of vertex 1 (tiles (x,1) with bit 0) and
**e_n(t)** = the out-tile-degree of vertex n (tiles (n,y) with bit 1) — the
two LEGS of the staircase (vertex 1's column, vertex n's row); the tile (n,1)
they share is the APEX.

## (1) THE LEG LAW (PROVED — every tiling, not only blue)

> **Δx = x(t̄) − x(t) = 8·(e₁(t) − e_n(t)).**

*Proof.* With a_v = (U+L)_v − 2e_v one has d_v(t̄) = d_v + 2a_v, so
Δx = 4Σ_v a_v(d_v + a_v). The bookkeeping s_v = [v≥2] + e_v and
(U+L)_v = (v−2)⁺ + (n−1−v)⁺ give d_v + a_v = 2[v≥2] − (n−1) + (U+L)_v,
which is 0 for EVERY v = 2, …, n−1 (interior: 2 − (n−1) + (n−3); the boundary
cases v = 2, n−1 check the same way), while
d_1 + a_1 = −1 and d_n + a_n = +1 — identically, for every tiling. Hence
Δx = 4(a_n − a_1) = 8(e₁ − e_n), using a_1 = (n−2) − 2e₁, a_n = (n−2) − 2e_n. ∎

**The interior is silent.** Only the legs carry axis current; the drop-profile
invariant of THM-790's first draft collapses: its support is ALWAYS {1, n}.

## (2) The flip's score action: reversal plus leg defect (PROVED)

> **d_v(t̄) = −d_v(t) for every interior vertex v = 2,…,n−1;
> d_1(t̄) = −d_1 − 2; d_n(t̄) = −d_n + 2.**

The d=m line IS score-reversal with a ±2 defect at the two legs. (This is why
"lines with equal score multisets" cluster at complement-symmetric classes.)

## (3) The drop spectra in CLOSED FORM (PROVED; the predictive formulas)

The n−3 free column tiles, n−3 free row tiles, the apex, and the
m − (2n−5) interior tiles factor the layer's drop generating function:

> **Σ_{all tilings} z^{Δx/8} = 2^{m−2n+5} · (1+z)^{n−3} (1+z^{−1})^{n−3} (z + z^{−1}).**

The apex (the σ-fixed corner tile (n,1)) contributes the (z + z^{−1}) factor —
the unit current generator. For the blue sublayer, the σ-pairing locks the
legs: each of the n−3 column tiles mirrors to a row tile with the SAME bit,
contributing ±1 to e₁ − e_n, and the apex contributes its own ±1 — so
e₁ − e_n is a sum of n−2 independent signs (in particular e₁ + e_n = n−2,
the reflection identity at v = 1):

> **Σ_{blue tilings} z^{Δx/8} = 2^{(m+f)/2 − (n−2)} · (z + z^{−1})^{n−2},
> f = ⌊(n−1)/2⌋.**

Verified exactly: n=7 gives 16·(z+1/z)^5 ⟹ blue-line histogram
{8: 160, 24: 80, 40: 16} ✓; the all-tilings GF matches at n = 4..7 (TRUE).
The apex parity makes Δx/8 ≡ n−2 (mod 2) on blue — level blue lines exist iff
n is even, which is the parity law again from the generating-function side.

## (4) THE BLUE PARITY LAW (PROVED — now a one-line corollary)

Blue lines have Δx = 8(2e₁ − (n−2)):
**n odd ⟹ Δx ≡ 8 (mod 16) — blue lines never level; n even ⟹ Δx ≡ 0 (mod 16).**
Max |Δx| = 8(n−2), attained iff one leg is fully out and the other fully in —
uniquely including the transitive line, whose entire drop 8(n−2) is carried by
the legs ("the strict ordering's transitivity leaves through the legs").

## (5) The node laws (PROVED)

- Blue (grid-symmetric) tilings have antisymmetric centered scores
  (d_{n+1−v} = −d_v): **every pure-blue or mixed class has a
  complement-symmetric score multiset; a class without one is PURE BLACK** —
  a per-node exclusion computable from the score sequence at every n.
- Blue lines never touch pure black (THM-787(1)); the transitive node's whole
  d=m connectivity is ONE blue pipe of drop 8(n−2), landing mixed in the
  certified `n=4..8` atlases.  General mixed landing remains open.

## (6) The half-tiling model and the recursion (PROVED / named)

σ has f = ⌊(n−1)/2⌋ fixed tiles (the anti-diagonal — the staircase's
hypotenuse), so **blue tilings = free tilings of the HALF STAIRCASE:
2^{(m+f)/2}; blue lines = 2^{(m+f)/2−1}** = 1, 2, 8, 32, 256, 2048 (n = 3..8).
The recursive reading (everything-is-the-triangle): the d=m layer's axis flow
factors as APEX × LEGS with silent interior — precisely the triangle's
recursive decomposition (overlap = interior, wiring = legs, apex = apex); the
half model quotients the hypotenuse; the QUARTER model (quotient the half
domain by its residual reflection through the apex) is the named next descent
— its fixed cells sit on the staircase's short diagonal, and its count
formula 2^{(half+fixed′)/2} is the next prediction to verify.

## The predictive-formula ledger (what can be written down at EVERY n)

1. The transitive node: fiber 1, one line, blue, drop 8(n−2). Its landing node
   is mixed in the certified n = 4..8 atlases; a general mixed-landing proof is
   not supplied by the lemmas here (named gap, codex-S9).
2. Blue tilings/lines: 2^{(m+f)/2}, 2^{(m+f)/2−1}.
3. The drop spectra in closed form (section 3): the layer GF and the blue GF.
4. Non-pure-black ⟹ complement-symmetric score multiset — a per-node exclusion
   from the score sequence alone.
5. All axis drops are ≡ 0 (mod 8) — by the leg law directly, and consistently
   with THM-785's identity x = n(n²−1)/3 − 8·C3 (codex): the axis is an affine
   image of the 3-cycle count, so the leg law is equally a C3-flow law:
   ΔC3 = −(e₁ − e_n) per line. THM-785's blue |ΔC3| histogram prediction
   {0: 640, 2: 960, 4: 384, 6: 64} at n=8 is the GF's histogram again — the
   two frames agree exactly.

## (7) The n = 8 certification (COMPLETE)

Predictions logged before the run: 6,880 classes; 4,096 blue tilings / 2,048
blue lines; blue |Δx| ∈ {0, 16, 32, 48}, GF 64·(z+1/z)^6 ⟹ line histogram
{0: 640, 16: 960, 32: 384, 48: 64}; transitive pipe 168 → 120. OUTCOME (the
subH-augmented classifier; an earlier 6,874-bucket interim was honestly flagged
by codex-S9 and is superseded): **6,880/6,880 classes; every prediction exact**,
including the full histogram, parity, max, blue-avoids-pureblack, and the
transitive pipe landing on a mixed class. Bonus n=8 tables (pure blue 3 at
x ∈ {168, 152, 152}; mixed 173; pure black 6,704 peaking at x = 40; the black
sea: 992,684 pureblack–pureblack lines): see
05-knowledge/results/metagraph_flow_n8_check_opus_S305.out.
