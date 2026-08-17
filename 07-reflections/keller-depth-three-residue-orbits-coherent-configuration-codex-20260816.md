# Depth-three residue saturation: the tournament is a colored rooted-tree configuration

**SYNTHESIS around PROVED THM-3542; speculative transfers are labelled.**
LRC(14), depth-four residue saturation, full decomposition-group equality,
and fixed-map Keller factorization uniqueness remain open.

## Inheritance pass

- Closest proved mechanism: THM-3540 splits the depth-two predecessor cubic
  as the birational ancestry point plus one irreducible quadratic.
- Canonical hostile: THM-3539 permits the residue image to be inertia alone;
  full global wreath monodromy does not determine residue orbits.
- Corrected near miss: MISTAKE-421 forbids using etaleness as if it forced a
  chosen coordinate derivative to be nonzero.  THM-3542 instead checks the
  literal eliminant degree, denominator units, and squarefreeness.
- Least-used sidecar: an injective symmetric observable on unordered pairs.
  Here it is `x_i+x_j`, not a guessed tournament orientation.

The live concept board was:

1. newest-prime residue action;
2. marked ternary ancestry;
3. point/pair resultants;
4. tournament/XOR quotient loss;
5. Keller grade factorization versus fibre-orbit factorization; and
6. LRC same-stalk spectral closure.

THM-3542 materially changes the first three.  It supplies only an analogy and
a warning to the last three.

## What the `1+2+6` factorization really says

The nine depth-two words split relative to the marked ancestry word `00` as

```text
{00}                     size 1,
{01,02}                  size 2,
{10,11,12,20,21,22}      size 6.
```

The exact specialized eliminant factors with these degrees.  Moreover, its
degree-three marked-child divisor is exactly the `1+2` part, so the sextic is
not an accidental numerical remainder: it is the two off-ray root blocks,
each carrying one ternary child fibre.

This is the lawful version of the earlier branch-tree intuition.  Fibre
degrees multiply under composition, while marked-tree orbit degrees grow by
the shell rule

```text
1, 2, 2*3, 2*3^2, ... .
```

These numbers are orbit sizes, not numerical atoms in THM-3541's Keller grade
monoid.  In particular, the factor degree `2` here is not a Keller grade and
`6` here does not assert a degree-six map decomposition.

## Why a tournament is the wrong exact object—and what replaces it

There is an intrinsic binary observable on predecessor blocks: the rooted
isomorphism type of the minimal subtree spanned by the marked leaf and an
unordered pair.  It has no canonical orientation.  Forcing one would add a
gauge that neither the resultant nor the valuation uses.

The exact object is instead a **colored complete graph with colored
vertices**, or equivalently the symmetric two-orbit part of a marked rooted-
tree coherent configuration:

- vertex colors are the marked point and its ancestry shells;
- edge colors are the marked-pair, unequal-shell, and equal-shell/LCA types;
- the residue group may refine these colors;
- THM-3542 proves that it does not refine them at depth three; and
- a valuation packet is one color class together with its multiplicity.

At `n=3` there are three vertex colors and six edge colors, hence nine
packets.  This gives a precise replacement for “everything is a tournament
of size four”: the four-way/XOR pictures are small quotients of a richer
two-orbit configuration, and orientation is lost unless an extra sidecar
supplies it.  The V4 bidirected-square lesson from the Boolean atlas is the
same warning in another chart: a symmetric pair relation can have missing or
two-way edges and need not admit an invariant tournament orientation.

The degree-twelve pair factor sharpens the warning.  It is irreducible over
`Q`, hence its Galois action is transitive, but its recorded modular cycle
pattern retains a `6+6` block split.  Transitivity of the orbit is not the
existence of a twelve-cycle—just as strong connectivity, a Hamiltonian path,
and a tournament certificate are different predicates.

## The quadratic packet law has a boundary recurrence

THM-3539 gives `n` vertex colors and `n(n-1)` edge colors at depth `n`, for a
total of `n^2`.  Passing from `n` to `n+1` adds

```text
(n+1)^2-n^2=2n+1                                      (1)
```

new relation types.  This suggests a sharper induction target than “compute
the whole residue group”: preserve the inherited `n^2` colors, then prove
transitivity only on the `2n+1` boundary colors created by the new shell.

This is **a program, not a theorem about the actual residue image**.  Its
cheapest hostile is a subgroup that is saturated on all inherited colors but
splits one new equal-shell color.  Any induction must transport a sidecar
that detects exactly that split.

## A scalable depth-four experiment

For `n=4`, the marked stabilizer has point-orbit sizes

```text
1,2,6,18
```

and unordered-pair orbit sizes

```text
1,2,6,6,9,12,18,18,36,54,81,108,                    (2)
```

which sum to `351`.  A successful good specialization with four point factors
and twelve injective pair-observable factors would close depth four by the
same refinement argument.

The direct symmetric resultant has degree `27^2=729` and then requires a
square extraction.  A better exact construction is the additive compound.
If `C` is the companion matrix of the degree-27 point eliminant, the operator

```text
C^[2](u wedge v)=Cu wedge v + u wedge Cv              (3)
```

on `wedge^2 Q^27` has eigenvalues `r_i+r_j`, `i<j`.  Its degree-351
characteristic polynomial is the pair-sum resolvent directly.  If sums
collide, replace `(3)` by the symmetric generic observable

```text
C tensor I + I tensor C + kappa C tensor C            (4)
```

restricted to the exterior square; its eigenvalues are
`r_i+r_j+kappa r_i r_j`.  Testing a few small rational `kappa` values gives a
cheap separating-observable gate without labeling the roots.

The other bottleneck is coefficient height.  Search the small rational
normalization grid for `(tau,lambda)` whose first three iterates avoid all old
image primes while minimizing bit length.  The simple `lambda=0` axis is a
hostile: it falls back onto `L` after two steps, so its attractive small
coordinates are not a good depth-four fibre.

## Harmonic-subset lens: useful only after typing the universe

Every subset of the natural numbers can indeed be embedded as a subseries of
the harmonic series, but that lens distinguishes rather than identifies the
objects here.  The off-ray shell-size support is

```text
{2*3^k:k>=0},              sum 1/(2*3^k)=3/4,          (5)
```

so it is harmonically thin even though those shells exhaust every unmarked
vertex at every finite depth.  THM-3541's numerical atom support, by contrast,
is harmonic-divergent because it contains the odd primes and their doubles.
Thus “exhaustive in a fibre tree” and “large as a subset of N” are orthogonal
predicates.  Equation `(5)` is a useful sanity check against conflating the
ternary orbit tree with the grade monoid.

## LRC boundary

The coherent-configuration move may inspire LRC bookkeeping: retain vertex
types, pair types, and their sidecars instead of forcing a tournament.  It
does not transport THM-3542 to LRC.  The LRC one-leg Fubini computation still
has separate `y`, `v`, and `x-Q` variables; no same-stalk ancestry map or
physical current has been supplied.  The source, target, and preserved
predicate are therefore absent, and LRC(14) remains open.

The honest next cross-domain question is narrower: can the same-stalk LRC
response array be organized as a colored two-orbit configuration whose
colors retain endpoint ancestry and phase?  A positive answer would compress
bookkeeping; only an additional current/noncancellation theorem could turn it
into spectral closure.
