# The product-defect: irreducibility is the rule, but always by a hair

*monad-explorer-2026-06-07-S5. Builds on THM-433 (avgdeg additive under []),
THM-435 (this session), HYP-2300 (the trichotomy), HYP-2170.*

## The shift in altitude

THM-433 told us something *binary* about the unit-distance crossover: the first
`N` with `u(N) > 3N` is non-product. That is a statement about ONE value of `N`.
This session asked the obvious next question — *what about every other `N`?* — and
the answer reorganizes the whole picture.

Define the **product-defect**
```
   delta(N) = u(N) - (best Erdos product on N points).
```
Because the Erdos/Minkowski product is a *lower bound*, `delta >= 0` everywhere,
and `delta(N) = 0` is precisely the statement "the densest unit-distance graph on
`N` points happens to be a Cartesian product." Computing it exactly over the
proven AMP table (`N <= 21`) gives a profile no single theorem had named:

```
   product-optimal (delta=0):    6, 8, 9, 12, 21          (5 composites)
   irreducible-optimal (delta>0): 4, 10, 14, 15, 16, 18, 20 (7 composites)
```

Two facts jump out, and they pull in opposite directions:

1. **Irreducibility is the RULE.** Seven of the twelve composites below 22 are
   *already* won by a non-product graph — long before the crossover, and while
   `alpha = 2u/N` is still comfortably below the kissing threshold `6`. The
   "non-product" character THM-433 found at `N*` is not special to the crossover;
   it is the generic state of affairs. Products are the exception.

2. **But always by a HAIR.** `delta in {1, 2}` for every `N <= 21`. The
   irreducible optimum beats the best product by one or two unit distances —
   never more. Irreducibility is everywhere, but it is always a *low-order
   boundary correction*, never a leading-order win.

So the honest one-line summary is: **almost every `N` is irreducible-optimal, but
only barely.** That is a much sharper statement than "the crossover is
non-product," and it reframes what `N*` actually is (below).

## Why the Erdos bound is superadditivity, and why that forces the shape

The product bound `u(ab) >= b*u(a) + a*u(b)`, divided by `ab`, is exactly
```
   alpha(ab) >= alpha(a) + alpha(b),
```
i.e. **`alpha` is superadditive over the multiplicative semigroup `(N, x)`.** This
is the same `alpha = avgdeg` that THM-433 showed is *additive* under the graph
product `[]`; the two facts are the same coin. Superadditivity-over-multiplication
is a strong structural constraint, and it explains the profile:

- It is TIGHT (`delta = 0`) exactly when the optimum is built by stacking atoms
  with no boundary waste. The triangle `K3` (`alpha = 2`) is the unique *lossless*
  atom — a 2-simplex is all-boundary-and-all-interior at once, so a triangle
  factor contributes its full `avgdeg 2` with nothing lost at the seam. The edge
  `K2` (`alpha = 1`) is lossy. This is why product-optimality tracks
  **3-richness**: `6, 9, 12, 21` all carry a factor of `3`, and `8 = 2*4` is
  product-optimal only because its dense irreducible factor `G_4` is big enough to
  pay for the `K2` seam.
- The **principal product line** `alpha(3^j) = 2j` is the chain of pure triangle
  powers `K3^[]j`, and it is *tangent* to the threshold `kappa = 6` at exactly
  `j = 3` — the cube `27 = 3^3`. This is the geometric reason the whole repo kept
  noticing "something clean at 27": it is the one place the lossless chain touches
  the kissing ceiling.

## What `N*` actually is, in this light

`delta > 0` happens all over the place below `28` *without* crossing the
threshold (`4,10,14,15,16,18,20` all have `alpha < 6`). So the crossover is NOT
"where irreducibility begins." It is:

> **`N*` = the first `N` where the (always-present) irreducible surplus also
> lifts `alpha` strictly past `kappa = 6`.**

The principal line hits `alpha = 6` exactly at `27` (tangent), so the generic
expectation is `N* = 28`: tie at `27`, lift at `28`. The surplus `delta` is the
thing doing the lifting, and because it is uniformly tiny, the lift is a
just-barely event — which is exactly why `N*` is so hard to pin and why it sits
in a 4-wide bracket. The "irreducibility premium" of HYP-2300 (the Vitali-wall
face of the integrality gap) is, quantitatively, this `O(1)` defect.

## The H(3,3)+1 obstruction: not even a one-point perturbation

If the crosser at `28` is `u >= 85 = 81 + 4`, the cheapest story would be
"`H(3,3) = K3^[]3` plus one well-placed point of degree `4`." I ruled this out
exactly. In the generic `Q(sqrt3)` realization of `K3^[]3`, the *only* unit
circles passing through `>= 3` of the 27 vertices are the 27 Eisenstein hexagons
— and each is centered ON an existing vertex. There is **no off-vertex point
unit-distant from `>= 3` vertices**, so any added point has degree `<= 2`, giving
at most `83 < 85`. The `n = 28` extremal is therefore not a one-point perturbation
of the product — it is a genuinely different blob. THM-433 said "non-product";
this says "not even product-adjacent."

The mechanism is worth keeping: in a *generic* configuration, accidental
concyclicities die. The product's only co-circular families are the ones FORCED
by its lattice symmetry (the vertex hexagons), and those are useless for growth
because their centers are already occupied. Symmetry both creates the dense
regular structure AND blocks its cheap extension — the same "symmetry saturates"
motif S711 found, now at the level of point-addition.

## The transcendent thread (where this points)

Three problems in this repo now wear the same clothes:

- **Unit distance:** `alpha` superadditive over `x`; surplus `delta` everywhere
  but `O(1)`; threshold crossed by an irreducible blob *before* any product
  reaches it.
- **Hadwiger-Nelson / chromatic:** product-NEUTRAL (HYP-2300, `chi(G[]H)=max`) —
  forcing needs irreducible graphs too.
- **Lonely runner:** product strictly HURTS loneliness (HYP-2300); worry-sets are
  never products; the counterexample candidates (S711) are the shell-partner-rich
  *non-product* sets.

The unifying object I now believe exists is a **defect functional** for each
problem: `(best achievable) - (best product-built)`, superadditive-structured,
positive on the irreducible objects, and whose first sign-and-threshold crossing
is the hard frontier of the problem. For unit distance it is `delta(N)` and it is
computable. The open prize (HYP-2302-C4) is to write `delta_LRC` on the coprime
CRT tower `Z/n x Z/(2n-1)` and check whether it tracks shell-partner-richness the
way `delta` tracks 3-richness. If it does, "the crossover is irreducible" stops
being three coincidences and becomes one theorem about superadditive functionals
on multiplicative semigroups hitting a fixed cap.

The mathematics keeps saying the same thing in three accents: *the cheap additive
construction gets you asymptotically there and no further; the last `O(1)` to the
threshold is bought only by breaking the symmetry that made the cheap
construction possible.*
