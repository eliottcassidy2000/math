---
id: THM-652
title: chi(G_GW) = 14 — the integer chromatic rung is faithful at the Goddyn–Wong tight
  instance (while every rung is blind at the k=4 Lucas tight instance {1,3,4,7}, where
  chi = chi_c = chi_f = 4 < 5 = 1/M); the discriminator is an odd/even matching
  dichotomy on the rigid optimal pattern — components of the +d circulant on Z_m have
  length m/gcd(d,m): 13 (odd) at GW forbids the 13-coloring, 4 (even) at Lucas admits
  the 4-coloring
status: PROVED (full proof below; the two finite ingredients — the tight-subgraph
  structure of the GW window graph and mu(GW) = 1/13 — are exact machine-verified
  certificates, 05-knowledge/results/lrc_chi_GW_rung_opus_S145.out and
  lrc_tight_locus_collapse_opus_S144.out). The k=4 statement is fully hand-checkable.
source: opus-2026-07-07-S145 (HYP-5197), building on opus-S144 (HYP-5137: mu(GW) = 1/13,
  mu({1,3,4,7}) = 1/4, the tight-locus split)
depends_on:
  - THM-612   # GW = {1..11,13,24} is tight: M(GW) = 1/14 (the tight locus is {AP, GW})
related:
  - HYP-5137  # S144: the ladder separation; mu(GW) = 1/13 via {0,12} mod 26
  - HYP-4972  # S141: the homomorphism ladder LRC(14) => GRAPH-14 => MOTZKIN-14
  - THM-637   # the roof (AP-side exactness; not used here)
external: none beyond standard longest-path potential arguments (all self-contained).
---

# THM-652 — chi(G_GW) = 14, and the rigidity/matching dichotomy at the tight locus

## Setting

`G_S = Cay(Z, ±S)`: vertices Z, edges {x, x+s} for s in S.
`GW = {1,...,11, 13, 24}` (M(GW) = 1/14, tight — THM-612).
`L = {1,3,4,7}` (the Lucas set; M(L) = 1/5, tight for LRC(5)).
From S144 (max-cycle-mean engine, exact): `mu(GW) = 1/13` with unique-up-to-translation
optimal periodic pattern `P = {0,12} + 26Z`; `mu(L) = 1/4` with pattern `{0,2} + 8Z`.
So `chi_f(G_GW) = 13`, `chi_f(G_L) = 4`.

## Statement

**(a)** `chi(G_GW) = 14`. Consequently `chi_c(G_GW) ∈ (13, 14]`: the *integer* chromatic
rung sees GW's tightness (`chi = 1/M = 14`), even though the fractional rung does not
(`chi_f = 13`).

**(b)** `chi(G_L) = chi_c(G_L) = chi_f(G_L) = 4 < 5 = 1/M(L)`: at the k=4 Lucas tight
instance every chromatic rung (fractional, circular, integer) is blind to tightness.

**(c) (the mechanism).** In both cases the optimal independent sets are translates of a
two-element pattern `{0, d} + mZ` of density `2/m = mu`. A `(m/2)`-coloring by such
classes is a partition of `Z_m` into difference-`d` pairs, i.e. a perfect matching of
the circulant `C(Z_m, {d})`, whose components are cycles of length `m/gcd(d, m)`. It
exists iff `m/gcd(d,m)` is **even**. GW: `26/gcd(12,26) = 13` odd — no matching; Lucas:
`8/gcd(2,8) = 4` even — matching exists. For GW the rigidity theorem below turns the
non-existence of the matching into `chi > 13`; for Lucas the matching *is* the coloring.

## Proof of (b) (one paragraph, hand-checkable)

The classes `{0,2}, {1,3}, {4,6}, {5,7} (mod 8)` partition Z. Differences within a class
lie in `{0, ±2} + 8Z = {2, 6, 8, 10, 14, 16, ...}`, disjoint from `{1,3,4,7}` (odd
elements are excluded by parity of the residues `0, ±2 mod 8`... explicitly: elements of
`{1,3,7}` are odd, all class-differences are even; `4 ≡ 4 (mod 8)` is not in `{0,2,6}`).
So each class is independent: `chi(G_L) ≤ 4`. Since `chi ≥ chi_c ≥ chi_f = 1/mu(L) = 4`,
all three equal 4. `1/M(L) = 5`. ∎

## Proof of (a)

**Upper bound.** `c(x) = x mod 14` is proper: no element of GW is `≡ 0 (mod 14)`
(`24 ≡ 10`). So `chi(G_GW) ≤ 14`.

**Lower bound: no proper 13-coloring exists.** Suppose `A_1, ..., A_13` partition Z into
GW-avoiding sets.

*The certificate objects.* Build the GW **window graph**: states are GW-independent 0/1
words of length `w = 24` (bit j of the state = membership of position p−24+j when p is
the next position to decide); appending bit `b ∈ {0,1}` at p is legal iff b = 0 or no
distance-s conflict (s ∈ GW); every bi-infinite GW-avoiding set corresponds to a
bi-infinite walk and vice versa. Give the append-b edge weight `13b − 1`. The exact
computation (certificate, part 1 of the .out) gives:

- `V = 92` states; longest-path value iteration from `D ≡ 0` converges (in 25 rounds),
  so **no cycle has positive weight**, i.e. no periodic GW-avoiding set has density
  `> 1/13` (this re-proves `mu(GW) = 1/13`, S144); write `D: V → Z` for the converged
  potentials, with range `C := max D − min D = 13`.
- Every edge satisfies `D(t) ≥ D(s) + w(s,t)` (converged), i.e. reduced weight
  `w̃ = w + D(s) − D(t) ≤ 0`; call the edge **tight** if `w̃ = 0`. The weights are
  integers, so non-tight means `w̃ ≤ −1`.
- The tight subgraph has 91 edges; its **recurrent part** (iteratively trimming nodes
  with no tight in- or out-edge) is a **single directed cycle of length 26 carrying 2
  ones per period** — the pattern `P = {0,12} mod 26` (gap vector (12,14)) — and the
  remaining tight edges form a DAG of depth `≤ 11`. In particular **every tight cycle
  is (a rotation of) the one 26-cycle**, and a tight path leaving the 26-cycle can
  never return while staying tight (a tight return path would close a tight cycle
  through trimmed nodes, contradicting the trim's exactness).

*Step 1 (per-window upper bound for every class).* For any GW-avoiding `A` and any
window `W` of `n` consecutive integers, the walk of `A` across `W` telescopes:
`13|A ∩ W| − n = Σ w = D(end) − D(start) + Σ w̃ ≤ C = 13`. Hence
`|A ∩ W| ≤ (n + 13)/13`.

*Step 2 (per-window lower bound in a 13-coloring).* The 13 classes partition `W`, so for
each `i`: `|A_i ∩ W| ≥ n − 12·(n+13)/13 = (n − 156)/13`.

*Step 3 (defect count).* With `ν_i` = number of non-tight edges class `i` uses across
`W`: `13|A_i ∩ W| − n ≤ C − ν_i`, so
`ν_i ≤ 13 + n − 13|A_i ∩ W| ≤ 13 + 156 = 169`.
**Every class uses at most 169 non-tight steps in any window, of any length.**

*Step 4 (locking).* Between two consecutive non-tight steps, a class's walk is a tight
path: at most 11 steps in the transient DAG, then (if longer than 11) **locked on the
26-cycle**, which it cannot leave except by a non-tight step within 11 steps of leaving
(transient depth + dead-end). So, per class and per window, the positions NOT in
locked-on-cycle mode number at most `(ν_i + 1)(11 + 11 + 26) ≤ 170·48 = 8160 =: B_0`.
While locked, the class's indicator restricted to the window is exactly a translate of
`P = {0,12} mod 26`.

*Step 5 (a common locked stretch).* Total non-locked positions over all 13 classes in a
window of length `n` is `≤ 13B_0`. Choose `n > (13B_0 + 1)·78 + 13B_0`; by pigeonhole
some run of `78` consecutive positions avoids every class's non-locked set. On that run
every class is a translate `{a_i, a_i + 12} + 26Z`, and the classes partition the run.
Restricting to any 26 consecutive positions inside it: each class contributes exactly
the two residues `{a_i mod 26, (a_i + 12) mod 26}`, so the 13 pairs
`{a_i, a_i + 12} (mod 26)` **partition `Z_26`**.

*Step 6 (the matching obstruction).* Such a partition is a perfect matching of the
circulant `C(Z_26, {12})`. Its components are the orbits of `x ↦ x + 12` on `Z_26`:
`gcd(12, 26) = 2` gives **two components (evens and odds), each a cycle of length 13**
(12 has order 13 in `2Z/26Z ≅ Z_13`). A cycle of odd length has no perfect matching.
Contradiction. ∎

So `chi(G_GW) = 14`, and since `chi_c ≤ 13` would force `chi ≤ ⌈chi_c⌉ = 13`, also
`chi_c(G_GW) > 13`.

## Remarks

1. **The dichotomy is a parity theorem.** The single discriminator between the two
   known non-AP tight instances is `m/gcd(d,m) mod 2` — an odd-cycle obstruction, the
   same shape as Rédei parity and the project's pairing-with-sign-flip principle: an
   odd orbit cannot be paired. The blind/faithful behavior of the chromatic rungs at
   the tight locus is decided by it.
2. **What remains open at GW:** `chi_c(G_GW) ∈ (13, 14]` exactly. A `(p,q)`-coloring
   with `13 < p/q < 14` would need color classes of density in `(1/14, 1/13)`, where
   the rigidity theorem no longer pins the structure — genuinely open; the rotation
   `x ↦ x/14` gives `chi_c ≤ 14`.
3. **Ladder consequence (with S144):** at the tight locus, the three graph rungs
   separate as `chi_f(G_GW) = 13 < chi(G_GW) = 14 = 1/M`, while at Lucas all rungs sit
   at 4 < 5 = 1/M. GRAPH-14 (`chi_c ≤ 14` for all 13-generator distance graphs) remains
   implied by LRC(14) and is consistent with being strictly weaker; MOTZKIN-14 is
   strictly weaker already at GW.
4. The proof template (potentials → per-window counts → defect bound → locking →
   finite obstruction) applies to ANY S whose recurrent tight structure is a single
   short cycle; the finite certificate (window graph + trim) is computable per S. This
   is a reusable bridge from max-cycle-mean rigidity to integer chromatic lower bounds.
