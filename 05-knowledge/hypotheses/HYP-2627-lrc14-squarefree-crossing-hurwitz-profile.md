---
id: HYP-2627
title: LRC(14) squarefree crossing-Hurwitz profile - K14 is the first Hill four-block product with radical 210
status: OPEN
source: codex-2026-06-19-S20
depends_on:
  - HYP-2625
  - HYP-2626
  - HYP-2624
  - THM-538
related:
  - HYP-2598
  - HYP-2617
  - HYP-2619
  - THM-539
  - OPEN-Q-108
---

# HYP-2627 - LRC(14) Squarefree Crossing-Hurwitz Profile

## Claim

The useful relation between the Markov-Hurwitz equation

```text
w^2 + x^2 + y^2 + z^2 = wxyz
```

and the complete-graph crossing formula is not that the LRC(14) tuple solves
the equation.  It does not.

The useful relation is this:

```text
For a four-block tuple q=(w,x,y,z),
Hill crossing product = wxyz = 4*cr_block(q).
Markov-Hurwitz is the energy-critical surface wxyz = w^2+x^2+y^2+z^2.
```

Thus Markov-Hurwitz supplies the mutation/recurrence archetype, while Hill's
complete-graph crossing formula supplies the four-block squarefree product
profile.  For LRC(14), the live carrier is the latter:

```text
q_14 = (5,6,6,7)
P_14 = 5*6*6*7 = 1260 = 2^2 * 3^2 * 5 * 7
rad(P_14) = 210
cr(K_14) = P_14 / 4 = 315.
```

This packages the HYP-2625/HYP-2626 modular hierarchy in one tuple:

```text
6,6  -> mod-6 universal skeleton
5    -> mod-30 address layer
7    -> prime-7 coimage seam
rad -> mod 210 squarefree address.
```

## Computation

Script:

- `04-computation/lrc14_squarefree_hurwitz_crossing_codex_s20.py`
- output: `05-knowledge/results/lrc14_squarefree_hurwitz_crossing_codex_s20.out`

The script computes Hill four-block products

```text
P_n = floor(n/2) floor((n-1)/2) floor((n-2)/2) floor((n-3)/2)
```

for `5 <= n <= 22`, their radicals, and their Markov-Hurwitz pressure

```text
(w^2+x^2+y^2+z^2)/(wxyz).
```

It also generates the positive Markov-Hurwitz mutation orbit from `(2,2,2,2)`
through maximum coordinate `10^8`.

## Exact Findings

The first Hill product whose radical is divisible by `210` is exactly `n=14`.
The first Hill product equal to `1260` is also exactly `n=14`.

```text
n=13: q=(5,5,6,6), P=900,  rad(P)=30
n=14: q=(5,6,6,7), P=1260, rad(P)=210
n=15: q=(6,6,7,7), P=1764, rad(P)=42
```

So the `210` address appears precisely at the first open LRC runner count.  The
new prime `7` is not a late decoration; it is forced by the complete-graph
four-block crossing profile at `K_14`.

The K14 Hill tuple is not Markov-Hurwitz:

```text
(5^2+6^2+6^2+7^2)/(5*6*6*7) = 73/630.
```

So it is far product-dominant rather than energy-critical.

The Markov-Hurwitz mutation tree has a real hidden recurrence.  On the normalized
branch `(1,1,u,v)` after dividing original solutions by `2`, the sequence starts

```text
1, 3, 11, 41, 153, 571, 2131, 7953, 29681, 110771
```

and satisfies

```text
u_{j+1} = 4 u_j - u_{j-1}.
```

But that mutation tree is not the mod-210 carrier: among generated original
solutions with max coordinate `<=10^8`, `0` have a coordinate divisible by `5`,
while `10` have a coordinate divisible by `7`.

## Interpretation

This separates three roles that were easy to conflate:

1. **Markov-Hurwitz recurrence.** Good archetype for mutation and integer
   recurrences, but not the LRC14 squarefree carrier.
2. **Hill crossing four-block product.** The useful squarefree carrier.  At
   `n=14`, its product has radical `210` and product `1260`.
3. **LRC14 coimage seam.** HYP-2626 says the fixed support-six residual uses the
   prime-7 unit seam; HYP-2627 explains why the `7` belongs in the same
   four-block packet as the mod-6 and mod-30 clues.

This also explains why `1260` kept recurring in the LRC14 workspace.  It is not
only the known low-`L` denominator from the exact measure extremizer.  It is also
the Hill four-block product for `K_14`.

The match with the known THM-523 local champion is exact.  For the perturbation
`12 -> 36`, the half-gap computation telescopes as

```text
15/36 - 2/5 - 1/70 - 1/504
  = 5/12 - 2/5 - 1/(5*14) - 1/(36*14)
  = 1/60 - 1/70 - 1/504
  = 1/(5*6*14) - 1/(6*6*14)
  = 1/(5*6*6*14)
  = 1/(2*7*6*6*5).
```

The `tau <-> 1-tau` symmetry doubles this to `1/(7*6*6*5)=1/1260`.
So the raw Hill `K_14` denominator is not merely numerically equal to the
low-measure denominator; it is the denominator in the exact two-speed clash.

HYP-2628 and HYP-2629 add a totient-copy reason for keeping the raw product.
The user's copy rule `sum_{d|n} c(d)=n` is `c=phi`, so `P_14` can be decomposed
into exact-period/Euler-copy mass by squarefree masks.  Raw `P_14=1260` has
full `{2,3,5,7}` copy mass `576`, but `cr(K_14)=315` has full copy mass `0`:
dividing by four removes the entire dyadic gate.  Thus the raw product is not
only the right denominator; it is also the first quotient level that still
carries the complete mod-210 copy profile.

## Tournament Analysis

The computation uses proof quotients as vertices, not runners or raw crossings.

Hamiltonian path:

```text
squarefree_210_crossing_profile
> lrc14_prime7_coimage_seam
> mod30_address
> markov_hurwitz_mutation
> raw_hill_blocks
> raw_crossing_number
> raw_runner_vertices
```

Pairwise observable: preservation of the LRC14 predicate, squarefree transfer
value, recurrence legitimacy, and crossing/Hurwitz explanatory value.

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3_cycles = 0
```

## Status

Open / exploratory.  HYP-2627 does not prove LRC(14).  It identifies a new
four-block squarefree profile that explains why mod `210` and denominator
`1260` are natural at `n=14`:

```text
complete-graph crossing product for K14
-> radical 210
-> mod-210 LRC14 address
-> prime-7 coimage seam and repeated-residue tail.
```

Next test: explain why the `12 -> 36` perturbation selects the adjacent
four-block row `(5,6,6,7)`, and ask whether the repeated-residue packets in
HYP-2626 can be parameterized directly by the squarefree profile of that tuple
rather than by raw support tuples.
