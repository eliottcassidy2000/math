---
source: kind-pasteur-2026-07-07-S74
status: VERIFIED (adversarially hardened) + one clean REFUTATION + a classical base case. A
  second-layer naive/coprime correction on monad-S11's degree-3 tail. Owner: rethink past work
  through coprimality at the LRC14 bleeding edge. Not a proof of R2.
tags:
  - lonely-runner
  - LRC14
  - coprime
  - gcd
  - triple-atom
  - density-floor
  - majorization
---

# The gcd-graded triple count is cumulative, not per-pattern

**kind-pasteur-2026-07-07-S74 (HYP-5187).** The owner's lens — *think in coprime; we were
wrong where we assumed a naive relationship that was really coprime* — landed the day the
fleet turned my S72 per-q reduction into its correct (coprime) form. monad-explorer-S11 found
the mechanism: the **triple atom law** `I(0,p,q;θ) = θ²·gcd(p,q)/q = θ²/q'`, where
`q' = q/gcd` is the **reduced** largest difference — *only the reduced (coprime) max-difference
matters*. My naive S72 "residue-spread" was refuted three ways in a day; monad's `gcd/q` is the
truth, and `Σ₃(AP)` is a **Pillai** (A018804, a sum of gcd's) convolution.

The broad coprime re-audit of past mistakes is already done by the fleet this session
(monad-S12's "eight past corrections, one lens": the dilation artifact = missing `gcd=1`;
MISTAKE-116 escapes = CRT freedom through fresh primes; Farey cells = determinant-1 coprime
adjacency; residue blindness = the floor is coprime to every finite quotient; …) and
mac-mini-S55's `W_q` dilation law `gcd(c,q)=1 ⟹ W_q(cE)=W_q(E)`. I do **not** repeat it. This
note contributes the one piece the degree-3 frontier still had backwards — **a second layer of
the same naive/coprime trap, one level up.**

## The load-bearing lemma, and its naive twin

monad-S11's route to "AP maximizes the degree-3 term" is the **layer-cake dominance**

> `M_m(E) := #{triples of E with reduced max-diff q' ≤ m} ≤ M_m(AP)` for every `m ≤ 7`,

which by Abel summation (the weights `1/q'` are **decreasing**) gives `Σ₃(E) ≤ Σ₃(AP)`. monad
stated it and tested it lightly (random + a few structured shapes).

There is an *obvious stronger* claim sitting right next to it — the naive one:

> **(per-pattern)** for **each** primitive coprime shape `π = {0,i,i+j}` (`gcd(i,j)=1`), the AP
> maximizes the number of affine copies `copies(π,E) = #{(a,g≥1): {a,a+gi,a+g(i+j)} ⊆ E}`.

Per-pattern would be cleaner (it implies the cumulative bound term by term) and it *looks*
right — the AP is the "most commensurate" set. **It is false.** Exhaustively (k=8: 2024
violations; k=13: 117): e.g. the primitive near-affine family `{0,2,4,…,22,25}`
(spread-AP with one bump — one of the very families that refuted my S72 per-q claim) has
**4 copies of `{0,2,5}` versus the AP's 3**, and `{0,3,…,33,37}` over-realizes `{0,3,7}`.
A non-AP can beat the AP at an *individual* coprime pattern.

But the **cumulative** `M_m(E) ≤ M_m(AP)` **survives** — 0 violations across the exact
adversary class that killed the naive per-q claim (near-affine bumps at every scale,
Fibonacci/Lucas `μ>M` separators, GW, prim-sat, parity record), 24 000 random shapes, and a
hill-climb that *actively maximizes* `M_7` and `Σ₃` — which returns **exactly the AP**. The
mechanism is a majorization: a non-AP that over-realizes one coprime pattern pays for it by
**depleting the lower-complexity patterns even more**. `{0,2,4,…,22,25}` steals one `{0,2,5}`
but loses six 3-APs; the cumulative-from-the-bottom count stays dominated.

## Why this is the coprime story, stated correctly

The extremality of the AP on the gcd-graded triple sum is a **cumulative (layer-cake /
Hardy–Littlewood majorization) property, not a term-by-term one.** `{0,…,k-1}` is
*coprime-complete* — a complete residue system mod every `q ≤ k`, so it realizes the whole
spread of reduced patterns **and** front-loads the low-`q'` (low-complexity) ones. The correct
invariant reads the *distribution* of reduced max-diffs, dominated from below; the naive
per-shape reading is exactly the kind of "assumed a naive relation" the owner is warning
against — and here it fails while its coprime, cumulative refinement holds.

The two facts fit together through the atom law itself: the weight `1/q'` it puts on a triple
is **decreasing in the reduced denominator**, so Abel summation reads precisely the *cumulative*
dominance and is blind to the per-pattern spikes. The coprime weight and the cumulative
dominance are the same coin — the frontier just had the near side (per-pattern) facing up.

## The base case is classical (and clean)

The bottom rung is provable and known. `M_2(E)` = the number of 3-term APs (reduced max-diff
`= 2` ⟺ the pattern `{0,1,2}`), and

> the AP maximizes the 3-term-AP count over all `k`-element sets, at
> `M_2(AP_k) = ⌊(k-1)²/4⌋` (12 at k=8, 36 at k=13),

via the midpoint/additive-energy identity `#3-APs(E) = #{(x,y)∈E², x<y, (x+y)/2 ∈ E}` (verified
exactly; hill-climb confirms AP-maximality: best non-AP 11<12, 31<36). This rung is both
per-pattern and cumulative — the base of the layer cake — and grounds any induction on `m`.

## Honest scope

- **Verified, hardened:** cumulative `M_m(E) ≤ M_m(AP)` and `Σ₃(E) ≤ Σ₃(AP)` (θ²-order,
  exact rational) against the per-q killer class + hill-climb, k=8 and k=13. This *upgrades*
  monad's lightly-tested lemma past the MISTAKE-102 weak-census gap.
- **Refuted (new correction):** the stronger per-pattern dominance — 2024 / 117 violations.
  Filing it so no one spends a session trying to prove the clean-looking false statement.
- **Proved (classical):** the m=2 base, `⌊(k-1)²/4⌋`.
- **Open:** cumulative dominance for `3 ≤ m ≤ 7` in general (a finite gcd-graded counting
  statement — the coprime-completeness of the AP; the right target for a real proof).
- **Does NOT prove R2 / (A').** Even `Σ₃(E) ≤ Σ₃(AP)` is only one degree-3 slice; monad-S11's
  honest negative stands — the full deficit has massive `Σ₃/Σ₄` cancellation, so the tail is a
  dichotomy (spread side Bonferroni-tight, structured side finite arithmetic). This nails the
  structured-side object as cumulative, not per-pattern.

## Ledger

- New: per-pattern gcd-graded dominance is FALSE (2024/117 exact violations) while cumulative
  `M_m` dominance holds (0, hardened) — the second-layer naive/coprime split on the degree-3
  tail. m=2 base pinned to the classical `⌊(k-1)²/4⌋`.
- Files: `lrc_coprime_completeness_kps_S73.py`, `lrc_coprime_perpattern_kps_S73.py` (+outs).
- Builds on: monad-S11 (atom law `θ²·gcd/q`, Pillai `Σ₃`, the cumulative lemma stated),
  monad-S12 + mac-mini-S55 (the broad coprime audit — cited, not duplicated), opus-S144
  (near-affine adversaries, V_r invariant frame), kps-S72 (the refuted naive per-q claim this
  corrects), kps-S67 (σ-grading: the cumulative count is the σ-odd-computable coefficient of a
  σ-even functional).
- Does NOT prove LRC(14), R2, or (A'). It fixes which coprime statement is true.
