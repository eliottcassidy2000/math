# The Burnside core-kernel is a GCD/Euler-φ quadratic form — and the metric I was missing is the whole enumerator family

**Source:** kind-pasteur-2026-06-15-S7. Dispatch: make progress on exact recurrences for
tournament structure; survey which metagraph sequences have been explored; ask "what metrics
am I missing?"; let it get complicated and reframe. Anchored on the user's compression of
A000568: split each odd partition `λ = μ ∪ 1^r` (core `μ` = odd parts `≥3`, tail `1^r`), the
Burnside weight factors, and `a(n)=Σ_{m,t} B[m,t] 2^{C(n−m,2)+(n−m)t}/(n−m)!` with the
`n`-independent **core kernel** `B[m,t]` (834 active states at `n=100` vs 444793 odd
partitions — a 533× collapse).

## The reframe: the residual difficulty is a theta-function exponent

The 1-tail peels off cleanly; what remains hard is the cross-term `Σ_{i<j} gcd(λ_i,λ_j)` in
the edge-orbit count `e(μ)`, which refuses to factor over `(m,t)`. Smith's 1875 GCD identity
`gcd(p,q)=Σ_{d|p,d|q}φ(d)` dissolves it. With `M_d` = #parts of `μ` divisible by `d`:

> **`e(μ) = C(t,2) + ½ Σ_{d odd ≥3} φ(d) M_d²`** (verified, 1113 cores, 0 failures).

The Burnside exponent is a **positive-definite quadratic form on the divisor-multiplicity
lattice** `(M_3,M_5,M_7,…)` weighted by Euler's `φ` — i.e. `2^{e(μ)}` is a **theta-function
weight** on that lattice. The repo has met `φ` as the GCD-matrix eigenstructure before
(Smith determinants `∏φ(d)`); here it is the *exponent* that controls the growth of A000568
once the score/fixed-point layer is removed. The exact add-a-part recurrence
`Δe = (p−1)/2 + Σ_{d|p} φ(d) M'_d` makes the obstruction precise and honest: the increment
needs the **full divisor profile** `{M_d}`, so no `(m,t)`-only recurrence can close — the
difficulty is *exactly* the GCD cross-term, now named.

## The 1-tail is the score layer; the cores are the OCF non-spectral family

Two structural rhymes make this more than bookkeeping:

1. **The tail = the trivial/ranking layer.** The peeled `1^r` are the *fixed points* of the
   automorphism. This is the same split THM-511 found on the arc cube: the **level-1 / score
   (converse-odd)** "ranking" layer vs the higher cyclic content. The leading core `m=t=0`
   gives the rigid asymptotic `a(n) ~ 2^{C(n,2)}/n!`; the cores are the automorphism
   corrections.
2. **The cores ARE the OCF non-spectral carriers.** The survey flagged a standing mystery:
   *why does `A000009(n)−3` (partitions into odd parts ≥3) count the OCF non-spectral defect
   dimension (THM-505)?* Because that is **the same family** as the Burnside cores `μ`. Odd
   parts ≥3 index both (a) the automorphism core in the iso-class count and (b) the
   disjoint-odd-cycle packing in the OCF. Both are "the irreducible part after removing the
   trivial layer" (fixed points / score). The `1`-cycles of an automorphism and the
   length-`1`... there is no length-1 odd cycle in a packing — the floor `≥3` is the same in
   both stories because a 1-cycle is a fixed point (no edge-orbit freedom) and a "1-cycle" is
   no cycle (no OCF content). The two `3`-floors coincide for the same reason.

## What metric was I missing: the enumerator family, not the single value

The project tracks many metagraph sequences as *separate* objects (A000568 for `G_n`,
A002854 for the even-graph `E_n`, A000088 for graphs, the SC spine count, the level-edge
counts `0,0,1,15,136`, …). The missing metric is the **2-variable enumerator** that has each
as a single evaluation, and the **`(m,t)`-graded refinement** of all of them at once:

- `T_n(x) = Σ_{odd λ} x^{e}/z` — orientations; `x=2 → A000568`, `x=1·n! → A000246`.
- `G_n(x) = Σ_{all λ} x^{e_graph}/z` — graphs; `x=2 → A000088`.
- `SC(n)` — self-complementary; the **base-4** Burnside with a `2^{#parts}` fixed-vertex tax
  (survey) — *literally the same 1-tail peeling in base 4*.

and the integer polynomial **`P_n(x) = n!·T_n(x) = Σ_{odd-cycle σ} x^{#edge-orbits}`** — the
edge-orbit enumerator of odd-cycle permutations (`P_7(x)=720x³+504x⁵+280x⁷+70x¹¹+x²¹`, rows
of its coefficient triangle summing to A000246). This triangle was **not found in OEIS** this
session — a candidate new sequence (HYP-2539). The point: the `(m,t)` core-kernel compression
and the φ-reframe are **not special to counting tournaments** — they are a general tool for
*any* `S_n`-Pólya metagraph count whose weight factors over fixed points. The repo's separate
metagraph sequences are evaluations of one cycle-index polynomial, all compressible the same
way, all with the same `φ(d)M_d²` exponent (with the per-part rule set by the structure:
`(p−1)/2` for orientations, `⌊p/2⌋` for graphs, base-4 for self-complementary).

## The complicated reframe, in one line

> Peel the 1-tail (the score/fixed-point layer) off the Burnside sum; the rest is a
> theta-function `Σ 2^{½Σφ(d)M_d²}` over the divisor-multiplicity lattice of "odd parts ≥3"
> cores — the very family that also carries the OCF's non-spectral defect — and that one
> kernel, re-based, generates the whole metagraph-sequence family (tournaments, graphs,
> self-complementary, edge-colored) rather than a single count.

## Status / honesty

VERIFIED: the φ-reframe and recurrence (0 failures to mass 40/35), the `(m,t)` compression
(== A000568, `n≤16`, exact), the enumerator-family anchors (`n≤8`). This does **not** solve
A000568 — it isolates and *names* the difficulty (the GCD-quadratic theta sum). HYP-2538
(a clean generating function / θ-product for the core kernel), HYP-2539 (the `P_n(x)` triangle
is new; the family is genuinely unified). Cross-links: THM-514, THM-505 (A000009−3 = the same
core family), THM-511 (the score/level-1 = the 1-tail), [[triangle_foundation]] (Mode-A vertex
insertion = adding a fixed point = a `1`-part), the Γ-species quotient view (arXiv:1204.1402),
tournaments↔even-graphs equinumerous (arXiv:2204.01947).
