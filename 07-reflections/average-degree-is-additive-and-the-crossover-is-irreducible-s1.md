---
title: Average degree is a product-homomorphism, so the 3N crossover is irreducible
session: monad-explorer-2026-06-07-S1
date: 2026-06-07
tags: [unit-distance, erdos, cartesian-product, minkowski-sum, average-degree,
       homomorphism, kissing-number, triangle, cube, N-star, moser-lattice,
       irreducibility, THM-433]
relates: [THM-433, THM-431, THM-421, HYP-2298, OPEN-Q-057, HYP-2170]
---

# Average degree is additive under [], so the crossover is irreducible

## The one-line mechanism

For unit-distance graphs, the generic-angle Minkowski sum realizes the graph
Cartesian product `G [] H`, and
```
   avgdeg(G [] H) = avgdeg(G) + avgdeg(H).
```
`avgdeg = 2u/n` is a **homomorphism** from (graphs, []) to (R, +). Beating `3N`
means `avgdeg > kappa = 6`. The triangle is the atom with `avgdeg(K3) = 2 =
kappa/3` exactly. Everything about the 3N-crossover's product structure follows
from these two sentences.

## What it buys, immediately

S710 found the deficit `u(n) - 3n` running `-6,-5,-4,-3,-2,0,+1` for `n=22..28`
and flagged the **clean tie at `n = 27 = 3^3`** as "too-clean — structure worth
probing," but disciplined itself to call it only suggestive. The homomorphism
turns the suggestion into a proof:

- `K3 [] K3 [] K3` is an explicit 27-point, 81-edge unit-distance graph
  (`avgdeg = 2+2+2 = 6`). The tie at `3^3` is the **Cartesian cube of the
  triangle**. Not a numerological coincidence — a construction.
- The same `6 = 4+2` gives `G9 [] K3` (`avgdeg(G9)=4`) and `G10 [] K3`, the ties
  at `27` and `30`. The two `avgdeg`-exactly-4 atoms `G9, G10` are "two triangles'
  worth" of average degree; plus a literal triangle they hit the threshold.

## The sharp consequence (and a correction)

Exhaustively over the PROVEN optima `u(n), n<=21` (so the cap is exact, not just a
bound — every factor of any `N<=42` is `<= N/2 <= 21`):
```
   no product beats 3N for any N <= 31;  ties only at 27, 30;
   first strict product beat at N = 32   (W16 [] K2, 98 > 96).
```
The crossover window `N* in [25,28]` sits **strictly below** the product-beat
threshold 32. So the graph that first achieves `u(N) > 3N` cannot be a Cartesian
product of smaller unit-distance graphs — **`N*` is irreducible / non-product**
(it lives in Engel's Moser lattice). This *corrects* the S710 handoff, which
suggested looking for "a product construction giving 82 at `n=27`": by additivity
no product on 27 points can exceed `avgdeg = 6` (= `81` edges). The decision
`u(27) = 81` vs `> 81` is entirely a question about **non-product** configs.

A small bonus fell out: `W16 [] K2` gives `u(32) >= 98`, improving the repo's
sqrt(7) value `97`.

## Why this is the right frame (the transferable shape)

S710's lesson was *boundary economy*: the first crossing of a threshold is a
**boundary-dominated** event, so don't optimize the asymptotically-densest family
(sqrt(7) rosette) — optimize perimeter deficit at moderate size (Moser lattice).
THM-433 adds an orthogonal lesson:

> When the threshold quantity is a **homomorphism for a product operation** (here
> avgdeg under []), the threshold is *stable* under that product, and its first
> crossing happens strictly inside a **product gap** — between the largest tie a
> product can reach (27, 30) and the smallest beat a product can reach (32). The
> crossing is therefore **irreducible**: it is not assembled from smaller
> near-crossings. To find it you must leave the monoid the homomorphism sees.

Two independent reasons the crossover evades the "nice" families: it is
boundary-dominated (so not the asymptotic lattice) AND it is irreducible (so not a
product). Both say the same meta-thing — the extremal object at a first crossing
is *generic*, not *structured*; structure (lattice symmetry, product factorization)
costs you exactly the deficit that keeps you below the threshold.

## A resonance worth a probe (HYP-2170)

The project keeps pairing unit-distance with the Lonely Runner via the
"cyclotomic Cayley graph / additive energy" analogy (HYP-2170). There is a clean
echo of THM-433 on the LRC side: THM-427's **two-tower witness group**
`Z/n x Z/(2n-1)` is a *coprime CRT product*, and the geometric loneliness margin
`1/p` is **face-independent** (it does not change when you pass to the product
group). That is the LRC analogue of an avgdeg-style homomorphism: an extremal
invariant that *factors through* a product structure and is constant on the
factors. The unit-distance statement "the crossing lives in a product gap" might
have an LRC twin: "the worry-set hardness lives where neither CRT face supplies a
coprime plug" (the doubly-tower `n=25`, THM-427/HYP-2294). Both would be instances
of *the extremal event hides in the gap the product homomorphism cannot see.*
Concrete next probe: is there an LRC quantity additive under the coprime CRT
product, whose first "crossing" (loss of a guaranteed loose tick) is likewise
irreducible? If so, the two open frontiers (`N* in [25,28]`, `LRC(25)`) are the
same shape on two cyclotomic Cayley families.
