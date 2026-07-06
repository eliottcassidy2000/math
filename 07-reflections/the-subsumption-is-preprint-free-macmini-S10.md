# The (A)-subsumption is preprint-free: an elementary L/(2N) rate

**mac-mini-2026-07-06-S10 (HYP-4342).**  Strengthens the S6b subsumption
(HYP-4302 — the coupled 2-torus emptiness (A) is a corollary of the finite 1-D
census) by removing its one soft dependency: the Jain–Kravitz lift-limit
preprint.  The convergence is elementary.  Verification:
`lrc_subsumption_rate_macmini_S10.out`.

## The old chain and its weak link

S6b: a proper 2-torus U = ⟨r, ℓ⟩ is the lift-limit of the 1-D families
v⁽ᴺ⁾ = r + N·ℓ; (i) M(v⁽ᴺ⁾) ≤ M(U) [sub-torus, trivial]; (ii) M(v⁽ᴺ⁾) → M(U)
[J-K, a preprint]; (iii) opus's residue bridge [formal]; (iv) finite census
[the open crux].  Step (ii) was the only non-elementary, non-corpus input.

## The elementary lemma (no J-K, no equidistribution theorem)

Write F(t,s) = minᵢ ‖rᵢt + ℓᵢs‖ on T².  Then, **exactly**,

    M(v⁽ᴺ⁾) = max_τ minᵢ ‖(rᵢ + Nℓᵢ)τ‖ = max_τ F(τ, Nτ) = max_{C_N} F,

where C_N = {(τ, Nτ mod 1) : τ} is the rational-slope-N curve in T².  Two facts:

- **F is L-Lipschitz** (sup metric) with **L = maxᵢ(|rᵢ| + |ℓᵢ|)** — a *fixed*
  constant, independent of N.  (Each ‖rᵢt + ℓᵢs‖ is L_i-Lipschitz with
  L_i = |rᵢ|+|ℓᵢ|; the min of Lipschitz functions is Lipschitz with the max
  constant.)
- **C_N is 1/(2N)-dense** in T².  For any (t₀,s₀), set m = round(Nt₀ − s₀),
  τ = (s₀ + m)/N.  Then |τ − t₀| = |m − (Nt₀ − s₀)|/N ≤ 1/(2N) and Nτ ≡ s₀
  (mod 1), so dist∞((τ,Nτ), (t₀,s₀)) ≤ 1/(2N).

Combining (Lipschitz F on a 1/(2N)-net):

> **M(U) − L/(2N) ≤ M(v⁽ᴺ⁾) ≤ M(U).**

M(v⁽ᴺ⁾) → M(U) from below with an **explicit rate L/(2N)** — pure real
analysis, no lift-limit theorem.  (Verified: the gap M(U) − M(v⁽ᴺ⁾) sits well
under L/(2N) for L=12 lift families; and M(v⁽ᴺ⁾) is itself a sharper lower
bound on M(U) than a coarse 2-D grid.)

## The preprint-free subsumption

Suppose a proper 2-torus U has M(U) ∈ (1/13, 2/25].  By the lemma
M(v⁽ᴺ⁾) → M(U) with M(v⁽ᴺ⁾) < M(U) (or = only in the limit), so for all large
N, M(v⁽ᴺ⁾) ∈ (1/13, 2/25) — **strictly in the open window** (even when
M(U) = 2/25 exactly, the approaching values are strictly below).  These are
in-window 1-D families v⁽ᴺ⁾ at unbounded height N.  opus-S98's residue bridge
makes "M(v) ≥ 2/25" a residue-class property (a q ≤ 50 witness clears), decided
for all heights by the finite Q50/template census.  If that census is
gap-empty, no v⁽ᴺ⁾ is in the window — contradiction.  Therefore
**M(U) ∉ (1/13, 2/25]**, endpoint included.

> **(A) [no coupled 2-torus in the window] is a corollary of (the finite 1-D
> census is gap-empty), with NO preprint** — only the trivial sub-torus
> containment, this elementary L/(2N) rate, and opus's formal residue bridge.

## What this closes

- The subsumption's ledger is now: **trivial + elementary + formal + (the one
  open census)**.  The J-K lift-limit citation (which S4's whole (A)/(C)
  reduction rested on, and which was flagged as outside the ≤13 policy) is
  **removed** from the (A) leg entirely.
- Combined with S9 (CircleClearFloor splits three ways; its residual is the
  phase-coupling of tight-locus frequencies): the direct covering lane and the
  subsumption lane are BOTH clean, and both terminate at the same single finite
  object — the residue-family census.  The 2-torus / distinct-frequency covering
  work is confirmed *sufficient-not-necessary*.

## Lean feasibility (the containment + rate as a formal object)

The lemma is real analysis on T² (Lipschitz + net), Mathlib-standard:
`LipschitzWith`, a finite ε-net, `Real.dist`.  The subtlety is that the corpus's
`Lonely` predicate is a *pointwise* statement, not a `max_τ` real; formalizing
M(v⁽ᴺ⁾) = max_{C_N} F needs the max as an object (compactness of the circle,
`IsCompact.exists_forall_ge`).  The trivial containment M(v⁽ᴺ⁾) ≤ M(U) is a
one-liner once both maxima are objects (a max over a subset ≤ max over the set).

**DONE (the abstract core): `LRCTorusRate.lean` GREEN, kernel-pure** (standard
trio, first compile): `lipschitz_ge_at_near` (an L-Lipschitz f at a point within
ε of the maximizer loses ≤ L·ε), `two_point_rate` (|f x − f y| ≤ L·ε), and
`exists_net_ge` (the max transfers to any ε-dense set with a Lipschitz loss —
instantiated with the torus maximizer and the curve C_N, this IS
M(v⁽ᴺ⁾) ≥ M(U) − L/(2N)).  What remains for a full corpus object is only the
LRC-specific glue: F = minᵢ‖·‖ is L-Lipschitz, and C_N is 1/(2N)-dense — both
Mathlib-routine — plus wiring M(v⁽ᴺ⁾) = sup_{C_N} F to the corpus's M.  The
mathematically load-bearing inequality is now formal.
