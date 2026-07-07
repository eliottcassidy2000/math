---
source: klein-2026-07-07-S157 (HYP-4811)
status: EXACT classification + a proved positivity mechanism. The owner's tiling-analogy
  directive, taken on its voltage-lift face: E[maxgap] record families ARE defective c-lifts
  (the SC-blowup analog); the systematic class bottoms out EXACTLY at the known record
  (12907/65520, +0.0057 above T*) — the mean sidecar survives the lift class. Plus the k=9
  Hunter floor rescued by 7-ADIC DESCENT (THM-638's G vanishes iff the pair's 7-adic
  valuations differ), verified. Proof-first session per owner directive (no Lean).
tags:
  - lonely-runner
  - LRC14
  - tiling-model
  - voltage-lift
  - tournament-analysis
  - density-floor
  - seven-adic
---

# The voltage-lift reading: the mean-record families are defective lifts

**klein-2026-07-07-S157.** Owner: the tiling model views a tournament from one of its
Hamiltonian paths, and that view reveals the iso-class symmetry forced by binary relations
themselves; apply the same analysis to lonely runners; proofs before formalization. Three
agents took the path/step/order-cell faces of this directive concurrently (kps-S62
step-gauge, mac-mini-S43 palindrome/wall-count, opus-S136 order-cells — the movie's
Hamiltonian path is the circular order σ(x), the gaps are the tiles, every pair relation is
a partial sum of gaps, x↦1−x is T↦Tᵒᵖ). This session takes the face nobody claimed: the
**voltage lift** (double round-robin / SC blowup, THM-378).

## 1. The dictionary entry

A **c-lift** of a base set A is c·A: its orbit at x is A's orbit at cx — a c-sheeted cover
in time, the exact analog of sheet-doubling a base tournament. A **defective lift**
`F = c·A ∪ B` breaks sheets with a small defect set B. The two record-holding families of
the mean program are precisely this shape: death-star's prim-sat `2·{1..12}∪{13}` and
monad's record `2·{1..11}∪{11,13}` (odd defects bisecting the even AP's doubled gaps —
"parity interlacing" = imperfect 2-cover). The tournament lesson (voltage classes measured
by triangle parities) transfers as: what a defect buys is set by its RESIDUE-depth against
the lift, not its size — same shape as THM-638's residue-graded pair masses.

## 2. The classification (exact): the record IS the class minimum

Systematic sweep of `F = c·{1..a} ∪ B`, `c ∈ {2,3}`, `|F| = 13`, all defect sets B in range
(8726 primitive-deduped families), numeric screen + exact rational confirmation of leaders
(death-star integrator):

| family | exact E[maxgap] | vs T* = 56291/294294 |
|---|---|---|
| `2·{1..11} ∪ {11,13}` (monad's record) | **12907/65520 = 0.196993** | **+0.005719** |
| `3·{1..11} ∪ {13,23}` | 33463369/169080912 = 0.197913 | +0.006639 |
| `3·{1..11} ∪ {10,26}` | 0.200904 | +0.009630 |
| `2·{1..12} ∪ {13}` (prim-sat) | 145091/720720 = 0.201314 | +0.010039 |

**The known record is exactly the minimum of the class** — no c ∈ {2,3} defect pattern goes
lower; the mean sidecar's +0.0057 margin survives the systematic lift attack. Notably the
runner-up is a **3-lift within 0.001 of the 2-lift record** (opus-S134's "trisection worse"
holds for their variants but the right 3-lift defects {13,23} nearly tie — the mechanism is
defect placement at the right cover depths, not the parity of c). And the **tail never
budges**: every leader has `μ_{1/7} ≈ 0.50–0.61 ≈ 10× m_P` — the tail/mean divergence
(death-star's capstone) in its sharpest form. HONEST SCOPE: c ≥ 4, mixed lifts
`c₁A₁ ∪ c₂A₂`, and non-AP bases are unexplored; within-class exactness only.

## 3. k=9 positivity rescued by 7-adic descent (proved + verified)

MISTAKE-122 left the bare k=9 Hunter-endpoint floor at exactly 0. THM-638 rescues it
arithmetically:

> **The pair mass exceeds θ² iff the two differences share their 7-adic valuation** —
> `G₊ > 0 ⟺ v₇(d_a) = v₇(d_b)` (the gcd strips common 7s; unequal valuations leave a
> reduced part ≡ 0 mod 7, the resonant row).

**Proposition (k=9).** For any 9-set E with top-anchored differences D (|D| = 8): if two
elements of D share a 7-adic valuation, then
`μ_{1/7}(E) ≥ meas W_top ≥ G₊(q_a,q_b)/(49·q_aq_b) > 0` ((q_a,q_b) the reduced pair; Hunter
tree through that edge; base `1 − 8/7 + 7·θ² = 0` exactly). Eight pairwise-distinct
valuations force `max D ≥ 7⁷ = 823543` — so positivity holds for every 9-set of diameter
`< 7⁷`, with an explicit per-shape constant. Verified on four probes (predicted floors
5/98 ≈ 0.051 ≤ measured W ≈ 0.27–0.32). Non-uniform (decays with the reduced product) —
a structural positivity, not the 0.562 leg bar; its value is the MECHANISM: Hunter floors
recurse through 7-adic levels, the arithmetic face of the apex-7.

The same remark strengthens k=8 below diameter 7⁶: some pair shares a valuation, so the
floor strictly exceeds 6/49 with a per-shape G-bonus.

## 4. What the analogy has bought so far, honestly

The Hamiltonian-path/tiling directive, one week in: (i) my Hunter tree "is" the path on 8
events — and the law makes specific trees exactly computable (`m(1,j) = 1/(7j)`: the star
and path masses of the AP endpoint are `(1/7)(H₇−1)` in closed form); (ii) gaps-as-tiles =
the order-cell decomposition (opus-S136's engine); (iii) reversal x↦1−x = complement
symmetry (already silently used in every mirror-dedupe); (iv) voltage lifts = the record
mechanism, now exactly classified; (v) steps/palindromes = kps/mac-mini's lane. The frame
is generating exact statements, not just metaphor — but no leg of the witness floor has
closed because of it yet.

## Files
`04-computation/lrc14_lift_defect_mean_class_klein_S157.py` (+.out; 8726 families, exact
leaders), `05-knowledge/results/lrc14_k9_sevenadic_check_klein_S157.out` (proposition
probes). Pointers: THM-638 (the law + C1 resonance), MISTAKE-122 (the bare-0), monad
HYP-4787/4807 (T*, the record), death-star-S1 (tail/mean capstone), opus-S134/S136,
kps-S62, mac-mini-S43 (the other faces), THM-378 (the tournament voltage-lift theorem this
reading transports).
