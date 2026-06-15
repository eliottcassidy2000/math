---
id: THM-511
title: The converse-parity theorem — on the full arc cube, T↦T^op is the global sign flip, so converse-EVEN invariants (H, c_k, the OCF/conflict graph) live on EVEN Fourier levels and converse-ODD invariants (scores) live on ODD levels; the tiling model's fixed base path breaks this symmetry and is exactly what injects the level-1 "ranking signal" (the FKN dictator direction)
status: PROVED (elementary character theory of T^op = global flip; the level-parity = converse-eigenvalue identity). VERIFIED exhaustively on the full arc cube n=4,5,6 (odd-level mass of c3,H = 0 to machine precision; scores supported on level 1 only) and on the tiling cube n=4,5,6,7 (scores level ≤1; c3 level ≤2; symmetry-broken odd levels of H present).
source: kind-pasteur-2026-06-15-S6
depends_on:
  - THM-002   # OCF: H = I(Ω,2); the conflict graph Ω is the level-2 even content
  - THM-468   # determinant lens det(I+S)=∏(1+μ²); S=A-Aᵀ is converse-ODD, the skew spectrum converse-EVEN
related:
  - THM-507   # walk counts are spectral (the converse-shift A↦A-J=-(Aᵀ+I) is the same T↦T^op seam)
  - THM-510   # the n=4 Boolean square: complement = leg-swap = converse; SC classes = converse-fixed
  - HYP-2532  # H's even-level weights = the OCF α_k strata
  - HYP-2533  # quantitative-Arrow / FKN-defect monotonicity
  - reflection: the-converse-parity-of-the-arc-cube-is-the-fkn-dictator-hierarchy-kps
---

# THM-511 — the converse-parity theorem (the FKN dictator hierarchy of a tournament)

**One line.** Put each tournament on the **full arc cube** `{±1}^{C(n,2)}`: coordinate
`x_e = +1` iff arc `e=(i,j)` (`i>j`) is in the reference orientation `i→j`. Then the global
sign flip `x ↦ −x` is **exactly the converse involution `T ↦ T^op`** (reverse every arc).
Hence for any invariant `f`:

> **`f(T^op) = +f(T)` (converse-EVEN) ⟹ `f̂(S)=0` for all odd `|S|` ⟹ `f` lives on EVEN Fourier levels.**
> **`f(T^op) = −f(T)` (converse-ODD)  ⟹ `f̂(S)=0` for all even `|S|` ⟹ `f` lives on ODD Fourier levels.**

The Fourier **level-parity equals the converse-eigenvalue** of the invariant.

## The dictionary (each line verified n=4,5,6 on the full cube)

| invariant | `f(T^op)` vs `f(T)` | converse | Fourier levels |
|---|---|---|---|
| Hamiltonian-path count `H` | `H(T^op)=H(T)` | EVEN | `L0,L2,L4,…` only |
| 3-cycle count `c_3`, all `c_k`, OCF cycle data | invariant under reversal | EVEN | even only |
| conflict graph `Ω`, `α_k` | invariant | EVEN | even only |
| score `s_v` | `s_v(T^op)=(n−1)−s_v(T)` | ODD (after centring) | `L1` (affine) |
| skew `S=A−Aᵀ` | `S(T^op)=−S(T)` | ODD | odd |
| skew spectrum `{μ_j}`, `det(I+S)=∏(1+μ_j²)` | invariant | EVEN | even |

VERIFIED (`04-computation/fkn_converse_parity_kps.py`): on the full arc cube, odd-level
mass of `c_3` and `H` is `0` to machine precision at `n=4,5,6`; `score[0]` is supported on
level `1` only; `H` is supported on `L0,L2,L4`.

## Two clean by-products

1. **`Var(c_3) = 3·C(n,3)/16`** over uniform labelled tournaments, sitting **entirely at
   level 2** (verified `n=4,5,6`): each triple's cyclicity has Fourier mass only on `∅` and
   its three arc-pairs (`f = ¼ − ¼Σ_{pairs} x_a x_b`), and distinct triples' cyclicities are
   **uncorrelated** (the shared single arc carries no covariance because cyclicity is even in
   each triple's three arcs). A single triangle's 3-cycle indicator has **no degree-3 term** —
   its two 3-cycle orientations are full reversals of one another, so the indicator is even
   under reversing the triangle's three arcs.

2. **The tiling model's base path is the symmetry-breaker.** Fixing the base Hamiltonian
   path `n→…→1` (the `n−1` base arcs) breaks the global flip, so on the **tiling** cube
   `{±1}^{C(n−1,2)}` the scores become a genuine **level-1 affine signal** (verified: every
   coordinate's score is affine, `n=4..7`) and `H` acquires odd-level content
   (`L1,L3,…`). The chosen ranking IS the chosen converse-symmetry-breaking reference.

## Why this is the FKN content

Friedgut–Kalai–Naor: a Boolean function with `≥1−ε` of its Fourier mass on **level 1** is
`O(ε)`-close to a **dictator** `±x_i`. FKN's historical home is **quantitative Arrow**: a
pairwise-majority social-welfare function that rarely produces a Condorcet cycle is close
to a dictatorship. A tournament *is* the pairwise-majority outcome; "rarely cyclic" =
"small level-`≥2` (converse-even cycle) mass" = "level-1 / score-determined" = near-transitive
= "near-dictatorial." So:

> **level-1 (odd, score) mass = the near-transitive / dictatorial / "ranking" content;
> level-`≥2` (even, cycle) mass = the genuine Condorcet-cyclic content.** The **FKN-defect**
> `(mass ≥ L2)/(total variance)` of `H` on the tiling cube grows `0.25 → 0.51 → 0.68` for
> `n=4,5,6` — score-determination of `H` fails progressively, quantitatively.

## Scope / honesty

PROVED: the level-parity = converse-eigenvalue statement (elementary; `χ_S(−x)=(−1)^{|S|}χ_S(x)`)
and the `Var(c_3)` value. VERIFIED computationally `n=4,5,6` (full cube) / `n≤7` (tiling).
CONJECTURED (HYP-2532): `H`'s even-level weights are exactly the OCF `α_k` strata
(`L0↔1, L2↔` conflict-graph edges, `L4↔α_2`), and (HYP-2533) the FKN-defect of `H` is
monotone in `n` with a quantitative-Arrow lower bound on level-`≥2` mass.
