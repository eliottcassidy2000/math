---
id: THM-1805
title: "THE VANDERMONDE IS A SIGNED SUM OVER TOURNAMENTS, AND IN/TRANSITIVITY IS EXACTLY WHAT SURVIVES vs CANCELS — the binary-forms↔tournaments bridge, proved. Expanding the alternating root polynomial ∏_{i<j}(x_i−x_j) (whose square is the monic discriminant) by choosing x_i or −x_j from each factor = orienting each edge = a TOURNAMENT T, giving ∏(x_i−x_j) = Σ_T sgn(T)·x^{score(T)} with score(T)_k = #wins of k and sgn(T) = (−1)^{#upsets}. Verified exactly n=3,4,5: the coefficient of x^d is 0 unless d is a permutation of (0,…,n−1), where it is ±1. SURVIVAL = TRANSITIVITY: a tournament has a permutation score sequence iff it is transitive (classical: distinct scores ⟺ no 3-cycle), so the surviving terms are exactly the n! transitive tournaments and EVERY intransitive tournament cancels — paired by the 3-CYCLE-REVERSING involution, which preserves all scores (each cycle vertex keeps 1 win/1 loss) and flips the sign (3 arcs reverse). So 'tournaments = in/transitivity' meets classical invariant theory: transitivity is what the discriminant sees, a directed 3-cycle is the cancelling unit, and the score sequence is the torus weight = the charge grading of the Laplace moment engine."
status: >
  PROVED / VERIFIED-EXACT. The expansion ∏(x_i−x_j) = Σ_T sgn(T)x^{score(T)} is an exact
  polynomial identity, verified against the direct Vandermonde at n=3,4,5 (equal as coefficient
  dictionaries). Survival = transitivity verified over ALL tournaments at n=4,5 (transitive ⟺
  permutation score, both ⟺ no 3-cycle). The cancellation Σ_{T: score=d} sgn(T)=0 for
  non-permutation d is a THEOREM (it IS the Vandermonde coefficient, exactly 0), and the
  3-cycle-reversing map is verified score-preserving and sign-flipping at n=4,5 — the combinatorial
  witness. The fully rigorous fixed-point-free involution is the classical Zeilberger-style
  construction (canonical 3-cycle choice stable under reversal); cited, not re-derived here.
  Not new mathematics globally (the tournament proof of the Vandermonde is classical), but the
  in-repo bridge — binary-forms invariant theory ↔ tournaments = in/transitivity ↔ the Laplace
  moment engine's charge grading — is new here, and it PROVES one concrete assertion of
  boxeph's THM-1800 program.
source: klein-2026-07-20-S385 (owner: work more on the representation theory of binary forms and how it relates to tournaments, which are in/transitivity itself; work the Laplace moment engine)
depends_on: []
related:
  - THM-1800  # boxeph's Laplace-moment-engine / GIT-nullcone program (this proves its Vandermonde claim)
  - THM-1790  # the Laplace moment engine (EMP floor) — score = the charge grading that engine averages
  - THM-343   # H (Ham-path count) — the other in/transitivity measure; transitive is the H=1 extreme
attribution: >
  The Vandermonde tournament expansion and the cycle-reversing sign involution are CLASSICAL
  (the standard combinatorial proof of the Vandermonde determinant), and were CLAIMED in the repo
  as a stub by boxeph-THM-1800 ("the Vandermonde's tournament expansion cancels ALL intransitive
  terms by a cycle-reversing sign involution"). This file PROVES/verifies that claim exactly and
  connects it to the moment engine; the framing credit is boxeph's, the classical result is
  nobody's to claim.
script: 04-computation/vandermonde_tournament_klein_S385.py (+ .out)
---

# THM-1805 — the Vandermonde is a signed tournament sum, and intransitivity cancels

**2026-10-14 representation-typing correction.** The exact tournament
expansion and coefficient cancellation below are unchanged. Alternation is
under permutation of roots, not a specified `SL_n` action; the discriminant
is the squared Vandermonde with its leading-coefficient factor. The
binary-form bracket and one-parameter charge specializations are stated
explicitly below. They do not by themselves transfer a Gaussian nullcone theorem.

## The bridge

The Vandermonde `Δ(x)=∏_{1≤i<j≤n}(x_i−x_j)` is alternating under
permutation of its root variables. For `f(t)=a∏_i(t-x_i)`,
`Disc(f)=a^(2n-2)Δ(x)^2`; in the monic case this is exactly `Δ²`.
No `SL_n` action on these variables is asserted.
Expand it by choosing, from each factor `(x_i − x_j)`, either `x_i` (**`i` wins**) or `−x_j`
(**`j` wins**). A choice over all `C(n,2)` factors is exactly an **orientation of every edge =
a tournament** `T`:

```text
∏_{i<j}(x_i − x_j) = Σ_T sgn(T)·x^{score(T)},
   score(T)_k = #wins of k,   sgn(T) = (−1)^{#upsets}   (upset = larger index beats smaller).
```

Verified exactly at `n = 3,4,5` (the tournament sum equals the direct Vandermonde as coefficient
dictionaries).

## Survival = transitivity

The coefficient of `x^d` is `0` unless `d` is a permutation of `(0,1,…,n−1)`, where it is `±1`.
And a tournament has a permutation score sequence **iff it is transitive** — the classical fact
that *distinct scores ⟺ no directed 3-cycle ⟺ transitive*, verified over all tournaments at
`n = 4,5`. Therefore:

> **The surviving terms are exactly the `n!` transitive tournaments** (one per linear order,
> `6, 24, 120` at `n=3,4,5`). Every intransitive tournament cancels.

The signed Vandermonde expansion retains only the transitive contributions
after grouping by score. This statement concerns `Δ`, not a tournament
expansion of the squared discriminant.

## The cancelling unit — a directed 3-cycle

Intransitive tournaments cancel in signed pairs under the **3-cycle-reversing map**: reverse a
directed 3-cycle `a→b→c→a`. This

- **preserves every score** — each of `a,b,c` keeps exactly one win and one loss inside the
  triangle, and no other edge moves; and
- **flips the sign** — three arcs reverse, so the parity of `#upsets` changes by `3` (odd).

So it pairs each intransitive tournament with an opposite-signed, same-score partner, and the
signed sum over any non-permutation score vanishes. The transitive tournaments are the only ones
with **no 3-cycle to reverse** — the fixed points, which is why they alone survive. (The fully
rigorous fixed-point-free involution needs a canonical 3-cycle choice stable under reversal —
the standard Zeilberger construction; the cancellation itself is already a theorem, being the
exact Vandermonde coefficient `0`.)

## Why this is the representation theory the owner asked for

- **Binary forms / invariant theory.** For ordered linear factors
  `ell_i=a_i X+b_i Y`, the `SL_2`-invariant pair bracket is
  `[ij]=a_i b_j-a_j b_i`. On the affine slice `ell_i=X-x_iY`, its product
  `prod_(i<j)[ij]` is exactly `Δ(x)`. For the quadratic
  `a(X-x_1Y)(X-x_2Y)`, the discriminant is `a²(x_1-x_2)²`, not `x_1-x_2`.
  This factor-coordinate statement is distinct from an unspecified
  `SL_n` covariant or the single determinant bracket on `n` vectors.
- **Tournaments = in/transitivity itself.** The expansion assigns each monomial a *net sign over
  all tournaments with that score*, and this net sign is `±1` exactly on the transitive (acyclic)
  orientations and `0` on everything with a cycle. Transitivity is not an add-on — it is the
  content of the alternating projection.
- **The Laplace moment engine.** A tournament term has multidegree
  `score(T)` under independent scaling of the root variables. Choosing
  integers `w_i` and restricting to the one-parameter substitution
  `x_i -> t^(w_i)x_i` gives charge `sum_i w_i score(T)_i`.
  This is a precise specialization from a multigrading to one charge.
  Relating it to the Gaussian engine additionally requires a map preserving
  its polynomial coefficients, Wick factorials, and target predicate;
  the grading correspondence alone does not supply that map.

## Scope

The tournament expansion and coefficient cancellation are exact in every
order; the recorded coefficient-dictionary checks cover `n<=5`. The
canonical cycle-reversing involution is cited rather than reconstructed.
The binary-factor and charge specializations above identify their actual
coordinate maps; no Gaussian nullcone or open-problem conclusion follows
from these identities alone.

*Files: `04-computation/vandermonde_tournament_klein_S385.py` (+ `.out`).*
