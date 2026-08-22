# Ternary ancestry and harmonic subsets remember where a level sits

**Status: elementary all-base proof candidate plus VERIFIED-EXACT companion,
pending independent audit.**  This extends the mechanism of proved THM-3510
from binary shortlex to every base `b>=2`, with the ternary Berggren tree as
the motivating case.  It supplies no cross-domain ancestry or LRC/JC result.

## 1. Every branch language is a harmonic subset, but the address matters

Give the empty word index one and order words first by length and then
lexicographically.  In a `b`-letter tree, the length-`n` words occupy

```text
I_n=[1+(b^n-1)/(b-1), 1+(b^(n+1)-1)/(b-1)-1].         (1)
```

Thus every branch language is literally a subset of the positive integers,
and therefore selects a subseries of the harmonic series.  This observation
is lossless only while the complete shortlex address is retained.  Replacing
the language by its level counts forgets where its selected integers sit
inside (1), and harmonic weights are not constant across a level.

## 2. Equal counts, different logarithmic densities in every base

For each `n>=1`, select either the first or the last block of `b^(n-1)`
indices in (1).  Both choices occupy exactly one `b`-th of the level.  Put

```text
c=1/(b-1).                                             (2)
```

After division by `b^n`, the level begins at

```text
c_n=c+(b-2)/((b-1)b^n).                               (3)
```

The first block is a mesh-`b^-n` Riemann sum on `[0,1/b]`, while the last
block is the same mesh on `[1-1/b,1]`, for the function

```text
f_c(x)=1/(c+x).                                        (4)
```

The shift `c_n-c` and the monotone Riemann-sum errors are both `O(b^-n)`;
their sum over all levels converges.  Hence the two one-level masses tend to

```text
mu_L=integral_0^(1/b) dx/(c+x)
    =log((2b-1)/b),

mu_R=integral_(1-1/b)^1 dx/(c+x)
    =log(b^2/(b^2-b+1)).                               (5)
```

They are strictly different for every `b>=2`, because cross multiplication
gives the exact certificate

```text
(2b-1)(b^2-b+1)-b^3=(b-1)^3>0.                        (6)
```

At the end of level `D`, the cutoff has logarithm `D log b+O(1)`.  Therefore
the two fixed regular languages have logarithmic densities

```text
delta_L=log((2b-1)/b)/log b,
delta_R=log(b^2/(b^2-b+1))/log b.                     (7)
```

The binary specialization is exactly THM-3510:

```text
b=2: (mu_L,mu_R)=(log(3/2),log(4/3)).                  (8)
```

For the ternary tree the new pair is

```text
b=3: (mu_L,mu_R)=(log(5/3),log(9/7)),
     (delta_L,delta_R)=(log(5/3)/log3,log(9/7)/log3).  (9)
```

So equal counts `3^(n-1)` at every ternary depth still do not determine the
harmonic/logarithmic density.

## 3. One ternary language with no logarithmic density

Let stage `k` contain

```text
r_k=b^(2^k)                                            (10)
```

consecutive levels.  Choose the first block throughout even stages and the
last block throughout odd stages.  The count on every level remains
`b^(n-1)`.  But

```text
(r_0+...+r_(k-1))/r_k ->0,                            (11)
```

so the newest stage dominates the complete-stage harmonic normalization.
The normalized sums tend to `delta_L` along even stages and `delta_R` along
odd stages.  By (6), the limits differ and the language has no logarithmic
density.

This is not pathology caused by irregular counts.  It is pure loss of the
within-level address.

## 4. The Berggren and Fibonacci consequences

The two Berggren parity trees of THM-3509 are ternary branch trees.  Shortlex
therefore sends any selected family of primitive triples to a harmonic
subset of `N`.  The induced weight is a weight on branch *addresses*, not on
Pythagorean size, hypotenuse, or Euclidean height; those are different
enumerations and need separate comparison theorems.

A bounded-width branch family has finite harmonic mass.  Indeed, the first
index at depth `n` is at least `b^n/(b-1)`, so at most `M` selected words on
that level contribute at most

```text
M(b-1)/b^n.                                            (12)
```

The sum of (12) converges.  Consequently every one-word-per-level ray has
logarithmic density zero in branch shortlex.  In particular, the selected
Fibonacci/Cassini ray of THM-3509 is harmonically thin even though it is
arithmetically distinguished.  Its golden recurrence controls *which* word
is chosen, while the exponential ternary level size controls the harmonic
mass.

Positive logarithmic density requires a macroscopically thick family—on
average, a positive fraction of the `b^n` addresses—not merely an infinite
or recurrent ray.

## 5. Where the tournaments of sizes four and six actually sit

Ternary depth two is

```text
F_3^2,
```

with nine words.  Lexicographic comparison gives the transitive tournament
`T9` and its 36 arcs.  It is not a tournament of size four or six.

There is, however, an honest Boolean face obtained by keeping the two extreme
digits:

```text
{0,2}^2={00,02,20,22}.                                 (13)
```

Lex order induces a `T4` on these four word vertices.  Its six comparisons
are the edges of `K4`.  After identifying `0,2` with `0,1`, the three nonzero
linear forms on `F_2^2` partition those six edges into the three perfect
matchings.  This is exactly the same finite representation grammar appearing
in THM-3510 and the tetrahedral atlas.

The six edges in (13) remain edge objects.  A separate tournament on those
six objects requires 15 pairwise comparisons and an added orientation.  The
ternary tree supplies neither automatically.  Likewise, the LRC chamber
`K4` and THM-3509 recurrence `K4` are typed carriers with different maps;
their shared size is not an ancestry transplant.

## 6. Reproduction and scope

Run

```text
python -B 04-computation/shortlex_bary_equal_count_harmonic_boundary_probe_20260816.py
python -B -O 04-computation/shortlex_bary_equal_count_harmonic_boundary_probe_20260816.py
```

The companion checks exact seams, counts, and rational harmonic masses in
bases `2..8`, the polynomial gap identity, superdominant-stage ratios, the
ternary `T9`/Boolean-`T4` face, its three matchings, and bounded-width ray
bounds.  Normal and optimized outputs agree.  Its semantic digest is

```text
8889cd4d5d3da976a099824e6242fef1a9b9f80fe83e8ae1a4c6fba5dbdf517f.
```

The claimed all-base theorem rests on the elementary Riemann-sum proof
(2)--(7), not finite extrapolation.  Independent audit is still required
before canon promotion.  No physical current, Pythagorean-height density,
cross-domain ancestry, automatic tournament, LRC, Keller, or Jacobian
consequence follows.
