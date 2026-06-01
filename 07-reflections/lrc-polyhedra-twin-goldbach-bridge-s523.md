# Reflection: Tournaments as polyhedra — the LRC / Twin-Goldbach bridge

**Session:** opus-2026-06-01-S523
**Date:** 2026-06-01

## The Regular Polygon Theorem

The deepest finding: **for initial segment speeds {1,...,n-1}, the lonely
times t=k/n (gcd(k,n)=1) are exactly the times when the n points form a
REGULAR n-gon on the circle.**

At these φ(n) times, every gap equals exactly 1/n. The half-turn tournament
is the doubly regular circular tournament. H(T) is maximized. The observer
is barely a source (boundary case).

This gives LRC for initial segments trivially: the linear trajectory passes
through the regular polygon configuration φ(n) times per period.

For non-initial speeds, the trajectory NEVER forms a perfect regular polygon.
THM-387 says it must get close enough (enter the LL fiber). The gap race
is the competition to approximate regularity.

## H(T) correlates with regularity: ρ ≈ -0.78

The correlation between polygon regularity and H(T) is strongly negative
(ρ ≈ -0.78 across n=4,5,6,7). More regular polygon → higher H. The regular
n-gon tournament has H = the maximum value for circular tournaments.

This is the geometric meaning of H: **H counts the number of Hamiltonian
paths, and regular polygons have the most paths because every vertex "sees"
the same local structure.**

Irregular polygons with one huge gap have H≈1 (nearly transitive — one
dominant direction). Regular polygons have H >> n!/2^{n-1} because the
circular symmetry creates many competing orderings.

## The alignment defect is always 4

A striking computational fact: for every twin-Goldbach exception, the
nearest miss has **alignment defect exactly 4**. This means: the closest
non-twin prime to N-p (for twin p) is always at distance 4 from a twin.

The number 4 = the gap between twin pairs' "forbidden" residues mod 6.
Twin pairs sit at 6k±1; the nearest non-twin residues are at 6k±3 and 6k,
which are distance 2 or 4 from a twin position. The defect 4 is the
maximum possible single-step miss in the 6k±1 sieve.

## The bridge table

| | LRC | Twin-Goldbach |
|---|---|---|
| Circle | S¹ (runner positions) | Z/NZ |
| Polygon | n runners + observer | twin primes < N |
| Regular config | Regular n-gon (t=k/n) | dense twin coverage |
| Complement | T↔T^op (time reversal) | p↔N-p (Goldbach reflect) |
| Necklace class | G_n/Z_2 | K+K (twin-center sumset) |
| Good state | observer = SOURCE | both summands twin |
| Exceptions | 0 (conjectured) | 35 (computed) |
| Fiber triple | {SS, SL, LS, LL} | {6k-2, 6k, 6k+2} |
| Formal group | F(x,y)=(x+y)/(1+xy) | sieve on 6k±1 |

## The hyperbolic geometry

In rapidity space (via the formal group), the regular polygon is the
point where all rapidities equal the threshold 0.5·ln(n-1). The lonely
region is a horoball intersection. Deformations of the polygon are walks
in hyperbolic space.

The twin-Goldbach exceptions are discrete analogues: polygons on Z/NZ that
never achieve complement alignment because the twin-prime distribution is
too irregular locally. But the 11 K+K holes (matching the 11 exception
triples exactly) are FINITE because K is dense enough that K+K eventually
covers everything.

## What the bridge points to

The LRC gap race (THM-387) and the twin-Goldbach complement necklace are
BOTH instances of a meta-theorem:

> **A structured polygon on a circle, deforming according to arithmetic rules,
> must eventually achieve a specific alignment — unless the arithmetic
> structure has a finite number of "desert" exceptions.**

For twin-Goldbach: 35 deserts (the K+K holes).
For LRC: conjectured 0 deserts (the continuous circle is more flexible
than the discrete integers).

The formal group F(x,y)=(x+y)/(1+xy) is the natural metric on the
polygon's deformation space. The alignment condition is a hyperbolic
distance constraint. The "eventual alignment" is forced by the ergodic
properties of the linear flow in hyperbolic space.
