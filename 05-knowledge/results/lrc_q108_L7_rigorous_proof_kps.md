# A rigorous proof of the L7 torus-line discrepancy bound `D_{p,q} <= 14/p`

**kind-pasteur-2026-06-21.** This packages the elementary, fully rigorous proof of the analytic
heart of the L7 closure (the sole open lemma of the kps-S23 LRC(14)-S3 ledger). It strengthens the
THM-562 skeleton and the HYP-2730 closure document. Every constant here is verified exactly
(`04-computation/lrc_q108_L7_discrepancy_proof_kps.py`, and the `Var=2/7`, `int=L/7`, equally-spaced
checks). The only external input is **Koksma's inequality**, a classical theorem; everything else is
elementary.

## Setup and statement

Fix coprime integers `p > q >= 1` with `1 < p/q <= 43/20`. Let `I_t = [t/7, (t+1)/7)` for `t in Z/7`,
and for `y in R` write `sec(y) = floor(7*{y}) in Z/7` (`{y}` the fractional part). Define the
**cell-occupancy** of the `(q,p)` torus geodesic:
```
   mu(i,j) = Leb{ v in [0,1) : sec(qv) = i  and  sec(pv) = j },     i,j in Z/7,
```
and the **L1 cell-discrepancy**
```
   D_{p,q} = sum_{i,j in Z/7} | mu(i,j) - 1/49 |.
```

> **Theorem (L7 tail).** `D_{p,q} <= 14/p`.

Since the L7 resonance correction satisfies `|R(p/q)| <= max_{i,j} g_B(i,j) * D_{p,q} <= D_{p,q}`
(with `0 <= g_B <= 1`; HYP-2730), this gives `|R(p/q)| <= 14/p`, the rigorous tail bound.

## Two elementary lemmas

**Lemma 1 (uniform 1D marginals).** For every `i`, `sum_j mu(i,j) = 1/7`; for every `j`,
`sum_i mu(i,j) = 1/7`.

*Proof.* `sum_j mu(i,j) = Leb{ v : sec(qv) = i }`. As `v` runs over `[0,1)`, `qv` runs over `[0,q)`,
so `{qv}` sweeps the circle exactly `q` times; the preimage of `I_i` (length `1/7`) has measure
`q * (1/7) / q = 1/7`. Symmetrically for the columns. ∎

**Lemma 2 (equally-spaced sub-arc starts).** Fix `i`. The set `{ v : sec(qv) = i }` is the disjoint
union of `q` intervals `J_m = [v_m, v_m + 1/(7q))`, `v_m = (i + 7m)/(7q)`, `m = 0,...,q-1`. On `J_m`
the point `pv` ranges over `[p v_m, p v_m + L)`, `L := p/(7q)`. Put `a_m := {p v_m}`. Then
```
   { a_m : m = 0,...,q-1 }  =  { pi/(7q) + r/q  :  r = 0,...,q-1 }   (mod 1),
```
i.e. the `a_m` are exactly `q` equally spaced points of gap `1/q` (offset `pi/(7q)`).

*Proof.* `a_m = { p(i+7m)/(7q) } = { pi/(7q) + pm/q }`. Because `gcd(p,q) = 1`, `m -> pm mod q` is a
bijection of `{0,...,q-1}`, so `{ pm/q mod 1 } = { 0, 1/q, ..., (q-1)/q }`. ∎

## Reduction of `mu(i,j)` to a one-dimensional average

Fix `i, j`. On `J_m`, substitute `u = pv` (`du = p\,dv`); then
```
   mu_m(i,j) := Leb{ v in J_m : sec(pv) = j } = (1/p) * h_j(a_m),
   h_j(a) := Leb{ s in [0, L) : sec(a + s) = j } = overlap of [a, a+L) with the I_j-strips.
```
Summing over the `q` sub-arcs,
```
   mu(i,j) = sum_{m=0}^{q-1} (1/p) h_j(a_m) = (q/p) * ( (1/q) sum_{m=0}^{q-1} h_j(a_m) ).   (*)
```

`h_j` is the periodic **trapezoid**: as `a` slides, `[a, a+L)` enters, fully covers (up to `min(L,1/7)`),
and exits the strip `I_j`. Two exact facts (verified):
```
   integral_0^1 h_j(a) da = L * (1/7) = p/(49q),         (product-of-lengths)
   Var(h_j) = 2 * min(L, 1/7) = 2/7,                     (since L = p/(7q) > 1/7 because p/q > 1).
```

## Koksma step and the bound

By Koksma's inequality, for a periodic function of total variation `V` sampled at `N` points of star
discrepancy `D*`, `| (1/N) sum h(x_m) - integral h | <= V * D*`. The `q` equally-spaced points of
Lemma 2 have `D* <= 1/q`. With `V = Var(h_j) = 2/7`:
```
   | (1/q) sum_m h_j(a_m) - integral h_j |  <=  (2/7) * (1/q) = 2/(7q).
```
Insert into (*) with `integral h_j = p/(49q)`:
```
   mu(i,j) = (q/p) * ( p/(49q) + err ),  |err| <= 2/(7q)
           = 1/49 + (q/p) err,           => | mu(i,j) - 1/49 | = (q/p)|err| <= 2/(7p).
```
Summing the `49` cells:
```
   D_{p,q} = sum_{i,j} | mu(i,j) - 1/49 |  <=  49 * 2/(7p) = 14/p.    ∎
```

## Remarks (sharper constants, apex zeros, and the integer grid)

- **Verified constant.** `sup_{p/q in (1,43/20]} D_{p,q} * p = 20/7` (true), `sup D*q = 12/7` (codex
  HYP-2736); both far below the proven `14`. The proof above is intentionally lossy (`Var * D*` and the
  49-cell triangle bound) but elementary and analysis-light.
- **Apex zeros (HYP-2733).** `D_{p,q} = 0` iff `7 | pq`: the integer-grid form `D_{p,q} =
  (1/(7pq)) sum_{i,j} | c_{ij} - pq/7 |` (with `c_{ij}` the cell-counts on the common `7pq`-grid,
  `sum_j c_{ij} = sum_i c_{ij} = pq` by Lemma 1) is zero iff every `c_{ij} = pq/7`, possible iff
  `7 | pq`. So the tail bound is `R = 0` on the apex-aligned subset, which the small-denominator
  resonance atlas avoids entirely.
- **Tail threshold.** `|R(p/q)| <= 14/p < cap_k - P2(B)` for `p > p* := 14/(cap_k - P2(B))`
  (`p* <= 68` for `k = 8..12`, binding `k = 10`); the finite window `p <= p*` is the exact atlas
  (exhaustively `0` violations, k=8..12).

## How this closes L7

`L7 = ` [finite atlas `p <= p*`, exact `p0_inf < cap`] `+` [tail `p > p*`: this Theorem, elementary]
`+` [finite-`f1`: `p0(E) = p0_inf + O(1/f1)` via THM-546, PROVED] `+` [worst base = dense even AP,
verified over 400+ bounded bases; `r=2` dominates `r>=3` via the pairwise reduction HYP-2734]. The
analytic content -- the joint `r>=2` Erdos-Turan-Koksma constant the ledger called the one unproven
input -- is the elementary `D <= 14/p` above. Remaining for a full LRC(14) proof: the L1-L6 chain, the
S1/S2 cases, an end-to-end audit, and (for machine certification) Lean formalization of the finite
atlas + a mathlib Koksma. -> THM-562, HYP-2730/2733/2734/2736, OPEN-Q-108.
