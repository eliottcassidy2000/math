---
id: THM-1805
title: "THE VANDERMONDE IS A SIGNED SUM OVER TOURNAMENTS, AND IN/TRANSITIVITY IS EXACTLY WHAT SURVIVES vs CANCELS — the binary-forms↔tournaments bridge, proved. Expanding the fundamental alternating SL_n covariant ∏_{i<j}(x_i−x_j) (the √-discriminant) by choosing x_i or −x_j from each factor = orienting each edge = a TOURNAMENT T, giving ∏(x_i−x_j) = Σ_T sgn(T)·x^{score(T)} with score(T)_k = #wins of k and sgn(T) = (−1)^{#upsets}. Verified exactly n=3,4,5: the coefficient of x^d is 0 unless d is a permutation of (0,…,n−1), where it is ±1. SURVIVAL = TRANSITIVITY: a tournament has a permutation score sequence iff it is transitive (classical: distinct scores ⟺ no 3-cycle), so the surviving terms are exactly the n! transitive tournaments and EVERY intransitive tournament cancels — paired by the 3-CYCLE-REVERSING involution, which preserves all scores (each cycle vertex keeps 1 win/1 loss) and flips the sign (3 arcs reverse). So 'tournaments = in/transitivity' meets classical invariant theory: transitivity is what the discriminant sees, a directed 3-cycle is the cancelling unit, and the score sequence is the torus weight = the charge grading of the Laplace moment engine."
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

## The bridge

The Vandermonde `∏_{1≤i<j≤n}(x_i − x_j)` is the fundamental **alternating** covariant of `SL_n`
— the square root of the discriminant, the lowest-degree object of classical invariant theory.
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

The discriminant "sees" only transitivity; intransitivity is invisible to it.

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

- **Binary forms / invariant theory.** The Vandermonde is `∏(x_i−x_j)`, the alternating `SL_n`
  covariant; in the `n=2` case it is literally the discriminant of a binary quadratic, and the
  general one is the bracket `[1 2 … n]` of the first fundamental theorem. Tournaments are its
  monomial expansion.
- **Tournaments = in/transitivity itself.** The expansion assigns each monomial a *net sign over
  all tournaments with that score*, and this net sign is `±1` exactly on the transitive (acyclic)
  orientations and `0` on everything with a cycle. Transitivity is not an add-on — it is the
  content of the alternating projection.
- **The Laplace moment engine.** The `score(T)` grading is the **torus weight** — the same
  `U(1)` charge grading that `Λ_s(u) = Σ_q g_q(s)u^q` carries in the moment engine (THM-1790's
  frame): `x_i ↔ ` weight `+`, the score `= ` net weight. "One-sided in charge" (the GMC(2)
  nullcone) is the weight-extreme (transitive-like) locus, and the sign-twisted survival here is
  the finite-dimensional shadow of the moment engine's charge projection `CT_u`. This is the
  concrete `n`-vertex realization of boxeph's THM-1800 reading "one-sidedness = Hilbert–Mumford
  instability for the hyperbolic torus."

## Scope

An exact identity plus a classical mechanism, proved/verified `n ≤ 5`, framed inside the repo's
own objects. It settles boxeph's THM-1800 Vandermonde claim and pins the score↔weight↔charge
dictionary linking tournaments, binary-form invariant theory, and the Laplace moment engine. No
open problem is resolved; the bridge is now load-bearing rather than a stub.

*Files: `04-computation/vandermonde_tournament_klein_S385.py` (+ `.out`).*
