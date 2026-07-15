---
id: THM-841
title: THE VIOLATION LADDER (the Benjamini–Hochberg lens on interval cores) — for λ < 1/(k+1) the FULL ORDER-STATISTIC structure of the danger indicators is exactly computable on Farey gaps. Let W(t) = #{j ≤ k : ‖jt‖ < λ} (the violation count = the density route's arc-occupancy variable). Then W is supported on the Farey-gap ends: inside the gap (a/i, b/j) of F_k, a point at distance s from a/i violates EXACTLY the multiples mi ≤ k with s < λ/(mi) (a divisor chain; nothing else — no third-fraction intrusion below 1/(k+1)), and symmetrically from b/j. Hence every ladder measure m_r(λ) = |{t : W(t) < r}| is an explicit piecewise-linear sum over gaps with breakpoints on the divisor-refined grid {λ : λ = (1/(ij))·(i j)/(m i + m' j)-type collisions} — and the density route's factorial moments S_ρ = E[C(W,ρ)] (LEM-006's exact inputs) are CLOSED-FORM rationals on interval cores: S_ρ(λ) = Σ_gaps Σ-over-divisor-chain-pairs of explicit interval lengths. r = 1 recovers THM-826 (m_1 = the Farey profile); the OU/BH synthesis: the violation process along t is a two-sided divisor-chain birth process on each gap
status: PROVED (support lemma = two-line corollary of THM-826 Lemmas 1-2) + REFEREED (m_r for r ≤ 4 and S_ρ for ρ ≤ 3 exact vs direct at k ≤ 10 over rational λ grids, zero mismatches; e.g. S_1 = 3/2, S_2 = 277/240 at k=6, λ=1/8). Original claim text: the support lemma (violations ⊆ the two end-fractions' divisor chains, λ < 1/(k+1)) is a two-line corollary of THM-826's Lemmas 1-2; m_r and S_ρ formulas refereed exactly for r ≤ 4, ρ ≤ 3, k ≤ 10 on a rational λ grid (script); the general-λ breakpoint taxonomy stated, full write-up of the S_ρ closed form deferred to the referee output tables
source: kind-pasteur-2026-07-15-S128 (cont.12; owner: BH/dobriban + LEM-006 merge; comprehensive-by-formula mandate)
depends_on:
  - THM-826   # the r=1 rung + the two lemmas whose corollary is the support statement
related:
  - LEM-006 (the factorial-moment E[W] Bonferroni ladder — its S_j integrals are now closed-form on interval cores)
  - THM-819 (the primitive harmonic law = m_1 at λ=1/(k+2)), THM-833 (the OU drift; W is the spatial side of the same danger bookkeeping)
  - BH/FDR reading (owner's dobriban/BH pointer): the speeds are k simultaneous tests, ‖jt‖ the p-values, W the rejection count; the ladder measures are the FDR operating characteristics of the interval core; the Farey gaps make the joint law EXACT — a rare case of fully explicit multiple-testing dependence
  - A139250 toothpick lead (owner): the divisor-chain two-sided growth on gaps is toothpick-like self-similar accretion; the 2-adic structure of the chains (m doubling) echoes THM-466's 2-adic H digits — recorded as a lead, not a claim
---

# THM-841 — the violation ladder

**Support lemma.** For λ < 1/(k+1) and t in the Farey gap (a/i, b/j): every violated speed is a
multiple of i (if the violation is at the left end) or of j (right end); speed mi is violated iff
s = t − a/i < λ/(mi). (Proof: THM-826 Lemma 1 nests every arc in its primitive fraction's arc;
Lemma 2 says only the two end fractions' arcs reach into the gap; the multiples of i with centers
AT a/i are the imprimitive copies — the divisor chain.) So W(t) = #{m : mi ≤ k, s < λ/(mi)} +
#{m' : m'j ≤ k, s' < λ/(m'j)} with s' the distance to b/j — a two-sided staircase in t.

**Ladder measures.** |{t ∈ gap : W ≥ r}| is a sum of explicit interval lengths cut by the thresholds
λ/(mi), λ/(m'j); each m_r(λ) is piecewise linear with breakpoints where thresholds collide with each
other or the gap length. m_1 = THM-826. The factorial moments S_ρ(λ) = ∫ C(W,ρ) dt are rational
polynomials in λ per segment — LEM-006's density-route inputs, exact on interval cores.

## Codex-S14 correction/addendum: the apparent divisor refinement collapses

The preceding potential-breakpoint description is correct but not sharp.  If `i,j` are the endpoint
denominators of a gap of `F_k`, then `i+j>k`, so

`min(floor(k/i),floor(k/j))=1`.

Thus at most one endpoint has a nontrivial multiple staircase.  More strongly, the collision of its
`m`-th copy with the other endpoint's primitive copy occurs at

`lambda_m = m/(j+mi)`

(after interchanging `i,j` if necessary).  For `m>=2` and `mi<=k`,

`m(k+1-i) = m(k+1)-mi >= (m-1)k+m > k >= j`,

so `lambda_m>1/(k+1)`.  Individual thresholds also cannot hit the opposite gap boundary in this
window.  Consequently the **only actual cross-collision in the scope of THM-841 is the primitive
one** `lambda=1/(i+j)`.  There is no divisor-refined or dyadic breakpoint tree: each gap's
`lambda`-parameter graph is a path with at most one kink.

This gives the deferred all-order closed forms.  Let `G_k` be the circular Farey gaps, represented by
their endpoint denominators `(i,j)`, put `phi(1)=1`, and set

`O_k(lambda) = sum_((i,j) in G_k) [lambda(i+j)-1]_+/(ij)`.

Then, throughout `0<lambda<1/(k+1)`, THM-826 gives `m_1`, while

`m_2(lambda) = 1 - lambda sum_(d<=floor(k/2)) phi(d)/d - O_k(lambda)`,

and for every `r>=3`,

`m_r(lambda) = 1 - (2lambda/r) sum_(d<=floor(k/r)) phi(d)/d`.

In particular every rung `r>=3` is globally linear and has no breakpoint at all.  For the factorial
moments,

`S_1(lambda)=2k lambda`,

and, for every `rho>=2`,

`S_rho(lambda) = 2lambda sum_(d<=floor(k/rho)) (phi(d)/d)
                  sum_(m=rho..floor(k/d)) C(m-1,rho-1)/m
                  + 1_(rho=2) O_k(lambda)`.

Indeed, away from a primitive left-right overlap, `rho` simultaneous violations must be copies of
one primitive fraction.  If their largest multiplier is `m`, there are `C(m-1,rho-1)` choices and
their common one-sided length is `lambda/(md)`; the factor two is for the two sides of that primitive
center.  A primitive left-right overlap contributes exactly one extra pair and never a triple.  This
also proves the formulas for `m_r`: for `r>=3`, `{W>=r}` is precisely the disjoint union of the
radius-`lambda/(rd)` neighborhoods of primitive centers with `rd<=k`.

**Toothpick verdict.**  The A139250 analogy is structural imagery, not an exact self-similarity.
For one endpoint chain write `D_p={1/m:1<=m<=p}`.  Its exact doubling law is

`D_(2p) = D_p union (1/2)D_p union {1/m : p<m<=2p, m odd}`;

the last set has `p/2` genuinely new states when `p` is a power of two.  The first actual witness is
the gap `(0/1,1/4)`, whose switch `1/3` is not inherited from `D_2` and `(1/2)D_2`.  Globally, the
literal Farey-end accretion has additions `b(1)=1`, `b(n)=2phi(n)` for `n>=2`; it already violates
the toothpick recurrence at `n=5`, where `b(5)=8` but `2b(1)+b(2)=4`.  The totals
`2 sum_(d<=n)phi(d)-1` happen to equal A139250 at `n=1,2,4,8` (`1,3,11,43`), explaining the small
dyadic mirage, but at `n=16` they are `159`, not `171`.

**THM-826 cross-correction.**  Its claimed complete breakpoint range must omit `1/(2k-2)` when `k`
is even (as the later interval-totient census also noticed): denominator sum `2k-2` would require
`(k-2,k)`, `(k-1,k-1)`, or `(k,k-2)`, none a coprime distinct Farey-neighbor pair for even `k`.
The profile formula itself is unaffected.

Exact replay: `thm841_no_dyadic_breakpoint_codex_S14.py/.out` checks the structural inequalities
through `k=64` and every displayed `m_r,S_rho` identity against an independent exact circle sweep
through `k=20`.

## Evidence log

- [x] support lemma refereed (all violations lie on end divisor chains): k ≤ 10, dense λ grid, exact
- [x] m_r (r ≤ 4) formula vs direct measure: exact match on the grid (script tables)
- [x] S_ρ (ρ ≤ 3) closed-form vs direct integral: exact on the grid
- [x] full breakpoint collapse and all-`r`, all-`rho` formulas (codex-S14 exact replay)
