# The r=6 sharp horn contains a mini runner problem—but only on the equal-step slice

A self-similar observation from work toward the r=6 uniform bound (THM-1132,
renumbered after THM-1123 had already been claimed).  The sharp horn and the
one-variable `G(σ)` bands are proved. Because the residual safe set is closed while
each danger tooth is open, a component of length **at least** `1/(7k)` already
contains a `k`-safe point. THM-1144 closes the former worst core/step-two ray, and
THM-1134 closes the entire step-two family over all cores and scales while also
supplying a general multiplier-chart cone plus a separated-ratio gate. Arbitrary
five-killer shapes outside those gates remain open.

For killers `{b,b+d,...,b+4d}`, their phases at time `t` are

`bt, bt+dt, ..., bt+4dt`.

After writing `phi=bt` and `sigma=dt`, the offsets are the arithmetic progression
`{0,sigma,2sigma,3sigma,4sigma}`. In the `phi` coordinate the danger centres carry the opposite
sign, which reflection removes without changing gap lengths. Translation by `phi` must avoid
five danger arcs of width `1/7`. This is a mini view-obstruction problem: for a *fixed* sigma,
how large a translation window do those five centres leave? It is not a reduction of arbitrary
killers, whose offsets need not be equally spaced.

The clean object is not the union of arcs but the cyclic gap vector of the five centres. If
`H(sigma)` is its largest coordinate, then

`G(sigma)=H(sigma)-1/7`.

Reflection reduces to `u=min(sigma,1-sigma) in [0,1/2]`, and sorting the centres gives

```
H(u)=max(u,1-4u)             on [0,1/4],
H(u)=u                       on [1/4,1/3],
H(u)=max(3u-1,1-2u)         on [1/3,1/2].
```

Thus the strict set `G>1/7` is

```
[0,5/28) union (2/7,5/14) union (3/7,4/7)
  union (9/14,5/7) union (23/28,1],
```

of exact circle measure `9/14`. The equality points are the eight finite endpoints. The minima
`G=2/35` occur at all four nonzero fifths, `sigma=1/5,2/5,3/5,4/5`; the old band computation
missed threshold crossings inside its alleged combinatorial cells.

Three perspectives now fit together:

1. **Topology matters at the horn boundary.** A closed component cannot fill an open tooth of
   equal length. Equality is safe, although an eventual perturbative landing argument may want
   a strict margin.

2. **The useful “vertices” are gaps, not runners.** A tournament on the five centres records
   only an orientation or cyclic order. It forgets the metric gap sizes, and `G` depends on the
   largest size. Tournament fingerprints can classify order changes, but they do not preserve
   the certification predicate without metric edge labels.

3. **The missing bridge is two-dimensional.** Along an actual time interval, `phi=bt` and
   `sigma=dt` move simultaneously. The slogan `L*b approximately G(dt)` freezes sigma and is
   not an equivalence. A proof needs a landing rectangle or a quantitative drift estimate, and
   then still only handles equal-step offsets. The arbitrary-offset tail is a higher-dimensional
   phase configuration problem.

The self-similarity remains useful: the parent runner problem produces a five-centre circle
problem with the same `2*7` window arithmetic. But the sharper lesson is structural. Quotienting
to a one-variable AP orbit preserves the metric gap predicate only on that orbit; it destroys the
degrees of freedom that make the full five-killer problem hard.

## Multiplier charts resolve the fixed-chart arity mirage

The mini-runner recursion sharpened after allowing all core-safe charts
`t=u/13`, not only `u=1`. Every at-most-five-point residue pattern has some
nonzero multiplier with a six-unit cyclic gap; the exact proof has ten affine
orbits. This produces a fixed Kakeya rectangle and the cone
`B>=17 max(A,80)`. For the step-two AP pattern, a stronger 792-core rectangle
atlas plus exact finite complement closes every legal scale. The faithful
vertex is therefore not a killer at one frozen time, but an affine residue
orbit together with the selectable chart. A single-chart tournament forgets
precisely the maximization that supplies the wide gap.
