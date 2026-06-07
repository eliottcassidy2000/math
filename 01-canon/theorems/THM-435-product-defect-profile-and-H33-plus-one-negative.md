---
id: THM-435
name: the unit-distance product-defect delta(N) profiles irreducibility over all
      N<=21; alpha=2u/N is superadditive over multiplication (= Erdos bound); and
      the generic K3^[]3 admits NO new point unit-distant from >2 of its 27 vertices
      (the n=28 crosser is NOT "H(3,3) + one point")
status: VERIFIED (elementary, exact rational + Q(sqrt3) arithmetic; over the PROVEN
        AMP table u(n), n<=21, arXiv:2412.11914)
date: 2026-06-07
session: monad-explorer-2026-06-07-S5
depends_on:
  - THM-433   # avgdeg additive under []; N* non-product; 27=K3^[]3
  - THM-431   # u(21)=57; N* in [25,28]
external:
  - "Erdos, 'On sets of distances of n points', AMM 1946 (product/Minkowski bound)."
  - "Alexeev, Mixon, Parshall, arXiv:2412.11914 (2024): exact u(n) for n<=21."
  - "Engel/Schade: u(28) >= 85 lower-bound construction (Moser lattice)."
---

# THM-435: the product-defect profile, superadditivity, and the H(3,3)+1 obstruction

THM-433 proved `avgdeg` is additive under the Cartesian/Minkowski product `[]` and
isolated the crossover `N*` (smallest `N` with `u(N) > 3N`) as NECESSARILY
non-product. This theorem makes the irreducibility QUANTITATIVE across all `N<=21`,
identifies the underlying algebraic structure (superadditivity over multiplication),
and settles the sharp sub-question of OPEN-Q-057.

Let `u(N)` = Erdos unit-distance maximum (A186705), proven exact for `N<=21` (AMP).
Let `alpha(N) = 2 u(N)/N` = best average degree on `N` points.

## (A) The product-defect function  [VERIFIED, exact over the proven table]

Define
```
   delta(N) = u(N) - max over NONTRIVIAL factorizations N=a*b (1<a<=b<N)
                       of  [ b*u(a) + a*u(b) ].
```
The bracket is the Erdos product lower bound `u(ab) >= b*u(a)+a*u(b)` (the edge
count of the densest product `G_a [] G_b`). Hence `delta(N) >= 0` always, and:

- `delta(N) = 0`  <=>  the global optimum on `N` points is (top-level) a Cartesian
  product. **N is PRODUCT-OPTIMAL.**
- `delta(N) > 0`  <=>  a genuinely non-product (irreducible) graph is strictly
  denser. **N is IRREDUCIBLE-OPTIMAL.**

Exact computation over `N<=21` (PROVEN u-table):

| N | u(N) | best product | delta | factor witness | verdict |
|---|------|--------------|-------|----------------|---------|
| 4 | 5 | 4 | **1** | 2x2 | IRREDUC |
| 6 | 9 | 9 | 0 | 2x3 | PRODUCT |
| 8 | 14 | 14 | 0 | 2x4 | PRODUCT |
| 9 | 18 | 18 | 0 | 3x3 | PRODUCT |
| 10 | 20 | 19 | **1** | 2x5 | IRREDUC |
| 12 | 27 | 27 | 0 | 3x4 | PRODUCT |
| 14 | 33 | 31 | **2** | 2x7 | IRREDUC |
| 15 | 37 | 36 | **1** | 3x5 | IRREDUC |
| 16 | 41 | 40 | **1** | 4x4 | IRREDUC |
| 18 | 46 | 45 | **1** | 2x9 = 3x6 | IRREDUC |
| 20 | 54 | 53 | **1** | 4x5 | IRREDUC |
| 21 | 57 | 57 | 0 | 3x7 | PRODUCT |

(Primes 5,7,11,13,17,19 have no nontrivial factorization = irreducible by
definition.) So among composites:
```
   PRODUCT-OPTIMAL (delta=0):    {6, 8, 9, 12, 21}
   IRREDUCIBLE-OPTIMAL (delta>0): {4, 10, 14, 15, 16, 18, 20}
```

**Key facts this exposes (sharpening THM-433 from binary to a full profile):**

1. Irreducibility is NOT confined to the crossover. Non-product optima already
   appear at `N = 4,10,14,15,16,18,20` — well below `N* in [25,28]` and below the
   threshold `alpha=6` (all these have `alpha < 6`). The crossover is just the
   first `N` where the irreducible surplus also pushes `alpha` past `6`.
2. The surplus is small: `delta in {1,2}` for all `N<=21` (max `delta=2` at `N=14`).
   Irreducible blobs beat the best product by only 1-2 unit distances — the
   "irreducibility premium" of THM-433/HYP-2300 is uniformly tiny here.
3. **3-rich vs 2-rich.** Every product-optimal composite except `8` contains the
   factor `3` (`6,9,12,21`); the triangle `K3` (`alpha=2`, a 2-simplex, all-boundary)
   is the only "lossless" product atom, while `K2` (an edge, `alpha=1`) is lossy.
   `8 = 2*4` is product-optimal only because its dense factor `G_4` (itself
   `delta=1` irreducible) is large enough to absorb the `K2` waste.

## (B) `alpha` is superadditive over multiplication  [PROVED; = Erdos]

`u(ab) >= b*u(a) + a*u(b)` is EXACTLY `alpha(ab) >= alpha(a) + alpha(b)`
(divide by `ab`). So `alpha: N -> R` is **superadditive with respect to the
multiplicative semigroup** `(N, x)`. Verified: 0 violations over all `a*b<=21`;
TIGHT (equality) exactly at the product-optimal `N` of Part (A).

**The principal product line** `N = 3^j` (the j-fold triangle product
`K3^[]j`, planar-realizable as a generic Minkowski sum of `j` unit triangles):
```
   alpha(3^j) = 2j   exactly:   alpha(3)=2, alpha(9)=4, alpha(27)=6, alpha(81)>=8.
```
This line is **tangent to the kissing threshold** `kappa=6` precisely at `j=3`
(`N=27`): `2j=6 <=> j=3`. So `27=3^3` is the unique power of 3 where the
maximally-product-optimal chain reaches — and (conjecturally) exactly touches —
the threshold. `K3^[]4` on `81` points is an 8-regular planar UD graph
(`alpha=8>6`), so products DO eventually beat `3N`; but the first product to beat
is `W16 [] K2` at `N=32` (THM-433-E), strictly LATER than the irreducible
crossover `N* <= 28`. The threshold is crossed by an irreducible blob before any
product reaches it.

## (C) The "H(3,3) + 1 point" obstruction  [VERIFIED, exact in Q(sqrt3)]

OPEN-Q-057 sub-question (1): is the `n=28` crosser literally `H(3,3) + 1` — a
28th point unit-distant from `>=4` vertices of `K3^[]3 = H(3,3)`?

Realize `K3^[]3` EXACTLY in `Q(sqrt3)` as the generic Minkowski sum of three unit
triangles rotated by distinct primitive Pythagorean angles
`(cos,sin) in {(1,0),(3/5,4/5),(5/13,12/13)}`. This is FAITHFUL: 27 distinct
points, exactly 81 unit distances, **6-regular** (degree sequence `{6: 27}`),
`alpha=6` (tie). [Generic angles are essential: `0/60/120` deg collapse to one
Eisenstein lattice (12 points) — see the realization note in the script.]

A 28th point `p` added at unit distance from a set `S` of vertices requires `S`
to lie on a UNIT circle centered at `p`; `deg(p)=|S|`. Any `p` with `deg>=3` is
the circumcenter of a triple of its neighbors, all at exact circumradius `1`.
Enumerating all `C(27,3)=2925` triples and keeping those with circumradius `=1`
(exact `Q(sqrt3)` circumcenter via Cramer):

- The ONLY unit circles through `>=3` vertices are the **27 vertex-centered
  hexagons** — each vertex's 6 neighbors lie on the unit circle centered AT that
  vertex (the Eisenstein hexagon). All 27 such centers COINCIDE with existing
  vertices.
- **There are ZERO off-vertex unit circles through `>=3` of the 27 vertices.**

Therefore **every genuinely new point is unit-distant from at most 2 of the 27
vertices**, so
```
   "K3^[]3 + one new point"  achieves at most  u >= 81 + 2 = 83  <  85.
```
**The `n=28` value `u(28) >= 85` (Engel/Schade) is NOT realizable as "H(3,3) plus
one boundary point."** The crosser is a genuinely different irreducible blob —
the direct geometric confirmation of THM-433's "non-product" classification, now
at the level of "not even a 1-point perturbation of the product."

**Honest caveat.** Part (C) is proved for the GENERIC (Pythagorean-angle)
realization of `K3^[]3`. A special-angle realization could in principle create
accidental concyclicities; ruling those out is a finite but separate algebraic
task (a candidate next step). Genericity is the natural reading — accidental
coincidences are non-robust — but the statement as proved is "the generic
product admits no degree-`>=3` extension point."

## What this settles / does not settle

- Settles OPEN-Q-057(1) NEGATIVELY for the generic product: the 28-crosser is not
  `H(3,3)+1`.
- Does NOT pin `N*`: the `[25,28]` gap stands. It removes the cheapest candidate
  construction (boundary-augmented product) and reinforces that only an
  irreducible Moser-lattice blob can sit at the crossover.
- Profiles irreducibility everywhere (Part A): a tool for the next explorer —
  the `delta`-pattern (3-rich product-optimal, 2-rich/large-prime irreducible)
  predicts where to look for non-product extremal graphs at OPEN `n>=22`.

Artifacts: `04-computation/unit_distance_product_defect_monad_s5.py`,
`05-knowledge/results/unit_distance_product_defect_monad_s5.out`. See HYP-2302
for the conjectural extension and the LRC transfer.
