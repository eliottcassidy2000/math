---
id: THM-554
title: The three tiling recurrences combine into the score partition function; tiling-model mean 3-cycle count is (C(n,3)+(n-2))/4
status: PROVED (combination + mean formula); engine verified vs brute n<=7, reaches exact c3-distribution at n=10
source: kind-pasteur-2026-06-20-S20
depends_on:
  - THM-553   # codex two-clock (beta,tau) tile address
  - THM-550   # codex half-tiling even/odd parity recurrence
  - THM-549   # mac-mini complement-quotient / OCF complement-invariance
  - THM-280   # grid reflection = complement
related:
  - HYP-2690  # codex address DP program (this theorem is its concrete score engine)
  - THM-442   # full-staircase third-difference recursion
  - THM-027   # BIBD c3-extremes
---

# THM-554 - Tiling Score Partition Function

This realizes codex's HYP-2690 address-DP program with an explicit engine, and answers
"combine the three recurrences (full beta-strip; half-tiling even; half-tiling odd) at the
same time to understand any tile and compute structure-vs-n efficiently."

## The combination (the three recurrences = one partition function)

Canonical tiling model: base path `n -> n-1 -> ... -> 1`; free tile `(a,b)`, `1<=b<a<=n`,
`a-b>=2`; bit 0 => arc `a->b`, bit 1 => `b->a`. Codex's two-clock address (THM-553):
`beta(a,b)=a` (birth), `tau(a,b)=a+b-1` (mirror crossing), gap `=a-b-1=2beta-tau-2`.

Since each tile adds **exactly +1** to `score(a)` (bit 0) or `score(b)` (bit 1), the joint
score generating function over all `2^{C(n-1,2)}` tilings FACTORS:

```text
Z_n(x_1,...,x_n) = ( prod_{v=2}^n x_v ) * prod_{(a,b) tile} (x_a + x_b).      (Z)
```

The three recurrences are the three ways `Z_n` grows / folds:

1. **BETA-clock = full recurrence (THM-442/553) = the incremental update.** Adding vertex
   `n+1` multiplies by the birth strip `beta=n+1`:
   ```text
   Z_{n+1} = Z_n * x_{n+1} * prod_{b=1}^{n-1} (x_{n+1} + x_b).
   ```
   (This is the score-weighted refinement of the s95 endpoint-insertion cube recursion
   `Q_{m_{n+1}} = Q_{m_n} x Q_{n-1}`.)
2. **TAU-clock = half recurrence (THM-550) = the complement reflection.** `R:(a,b)->(n+1-b,n+1-a)`
   is tournament complement up to relabel (THM-280). Every **complement-invariant,
   score-determined** invariant is `R`-symmetric, so the half-tiling `tau<=n` representatives
   carry its whole distribution — the **address quotient**. Parity of `tau` is the even/odd split.
3. The **even/odd half-tiling recurrences** are this same `Z_n` seen through `parity(tau)`
   on the fixed-line crossing layer `tau=n` (size `floor((n-1)/2)`).

Consequently every score-determined complement-even invariant is a linear functional of
`Z_n`, computed by the beta-step **without enumerating `2^{C(n-1,2)}` tilings**:
- score-sequence census = the coefficients of `Z_n`;
- 3-cycle count `c3 = C(n,3) - sum_v C(s_v,2)` (the leading OCF/`alpha_1` term);
- hence the H-max / regular census = the central coefficient.

## The mean 3-cycle count (closed form, PROVED)

```text
E_tiling[c3] = ( C(n,3) + (n-2) ) / 4.
```

Verified exact (Fraction) `n=3..9`: `1/2, 3/2, 13/4, 6, 10, 31/2, 91/4`. The uniform-tournament
mean is `C(n,3)/4`; the tiling model's fixed base path adds exactly `(n-2)/4`.

**Proof.** By linearity, `E[c3] = sum over triples {i<j<k} of P(triple is a 3-cycle)`.
A triple is uniform (`P=1/4`) UNLESS the base path fixes arcs inside it. Exactly the `n-2`
**consecutive triples** `{v,v+1,v+2}` contain two fixed base arcs `v+2->v+1->v` (a transitive
2-path), so the triple is a 3-cycle iff the remaining tile gives `v->v+2`: `P=1/2`. Every
other triple has at most one fixed (consecutive) edge, which does not constrain cyclicity, so
`P=1/4`. Hence `E[c3] = (C(n,3)-(n-2))*(1/4) + (n-2)*(1/2) = C(n,3)/4 + (n-2)/4`. QED.

So the fixed Hamiltonian path biases the tiling ensemble toward MORE 3-cycles, by exactly one
quarter per consecutive transitive 2-path.

## Verification / engine

`04-computation/tile_address_score_gf_engine_kps.py` (output
`05-knowledge/results/tile_address_score_gf_engine_kps.out`):
- `Z`-engine c3-distribution == brute over all tilings, `n<=7` (exact);
- exact c3-distribution by the beta-clock to **n=10** (68.7e9 tilings, ~95 s, 29.5M `Z`-states,
  `c3_max=40 = (n^3-4n)/24` the regular max); `c3_max` matches the max-3-cycle formula all n;
- regular (H-max/Paley) score-sequence census `1, 3, 91, 29157` for `n=3,5,7,9`.

## Scope

`Z_n` is exact for score-determined invariants (scores, `c3=alpha_1`, the leading OCF term).
It does NOT give full `H`/OCF (those need the non-score cycle-pair data `alpha_2,...`; THM-442
already records `H` is not cell-affine). Its use is HYP-2690's: address-first, compute the
score/`c3` layer incrementally and complement-halved, then attach cycle packets.
