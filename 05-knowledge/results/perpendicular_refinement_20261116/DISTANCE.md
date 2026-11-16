# The weighted fork tree: sharp hyperbolic correction and exact integral decoding

The recovered live binary skeleton from the [earlier Collatz order geometry packet](https://github.com/eliottcassidy2000/math/blob/9dec4d1efcbb1389b41b1b4116154b39950ed050/05-knowledge/results/collatz_order_geometry_20261116/README.md)
has four unordered outgoing length pairs: `(2,2),(2,4),(1,4),(3,4)`. The
weights measure the number of shortcut inverse steps on a compressed
edge. This proof applies to skeleton **fork vertices** (not necessarily to
each original intermediate or sterile vertex under the same placement).

## The weighted thirds embedding

Give every vertex an interval, starting with `[0,1]` at the root of the
skeleton. A length-`e` edge to the left/right child selects the
leftmost/rightmost fraction `3^(-e)` of the interval. Place each vertex
at its interval midpoint with height equal to the interval length in
the hyperbolic upper half-plane. An edge of length `e` then has

```
cosh(hyperbolic edge length) = (5/8)(3^e + 3^(-e)) - 1/4.
```

Coordinates are inexpensive: along the path, append `e` ternary
zeros or `e` ternary twos and then add one half of the terminal
interval length to the resulting left endpoint.

For two vertices use coordinates relative to their least common ancestor's
interval, translated and scaled to `[0,1]`. Let `a,b` be the two
weighted depths from that ancestor, `A=3^(-a), B=3^(-b)`, and let
`u` be the absolute difference of normalized midpoints. Their weighted
tree distance is `D=a+b`. In curvature minus one,

```
Q = u^2 + A^2 + B^2
cosh d_H = Q/(2AB)
d_H - D ln 3 = ln Q + ln((1 + sqrt(1 - 4(AB)^2/Q^2))/2).
```

Only the suffix addresses after the common ancestor matter. On
incomparable pairs, keeping the first `m` weighted ternary digits
of each suffix leaves uncertainty at most `3^(-m)` in `u` (a
shorter suffix can be read exactly). The depth-dependent factors
`A,B` require only the total weights. Writing `g=53/81`, the
correction has derivative in `u` at most `2 sqrt(2)/g` on this
domain: `u>=g`, `D>=4`, and `4(AB)^2/Q^2<=1/2`.
Thus the horizontal correction converges exponentially from the two
short addresses.

## Sharp bounds and no-float distance recovery

**Incomparable endpoints.** The sibling gap in the common ancestor
chart is `1-3^(-e_L)-3^(-e_R)`, at least `g=53/81` by the four
allowed fork types. Also `u<=1-(A+B)/2`, and `A,B<=1/3`.
Hence `Q<=1`: expand the square and use
`(5/4)(A^2+B^2)+AB/2 <=(3/2)(A^2+B^2)<A+B`.
This implies `d_H<=D ln 3`. For the other direction,
`A^2+B^2 >= (AB)^2/u^2` (since `B<=u`), so comparison of coshes gives
`d_H>=D ln 3+2ln u>=D ln 3+2ln g`. The log argument in this
comparison is positive because `D>=4`.
Both bounds are sharp as limits on the full skeleton: use a (1,4)
split and follow inner extremes for the lower bound or follow outer
extremes for the upper bound. Such a split occurs at 512 on the
doubling ray; its odd child 341 is already live-binary and its
next live fork by doubling is four steps away.

**Ancestor endpoints.** Now `A=1`, `u<=(1-B)/2`. The vertical
distance gives `d_H>=D ln 3`. At `D ln3+ln(5/4)` the cosh numerator
on denominator `2B` is `5/4+(4/5)B^2`, which exceeds or equals
`1+B^2+(1-B)^2/4`. This proves the upper bound; an extremal
descendant ray makes the correction tend to `ln(5/4)`.

**Exact recovery.** The global signed error interval is
`[2ln(53/81), ln(5/4)]`. Its width is
`ln(32805/11236)<ln3`, since `32805<33708`. All weighted
tree distances are integers; hence

```
D = floor((d_H + 2 ln(81/53))/ln 3).
```

No ancestor, edge class or sign information is required by this
decoder once the hyperbolic distance is known. It can be implemented
entirely with rational arithmetic. Set `gamma=(53/81)^2`. Compute
`C=cosh(d_H)` exactly from the rational coordinates; for each
integer `k>=1` test

```
C >= (gamma*3^k + 1/(gamma*3^k))/2.
```

The largest passing `k`, or zero if none passes, is `D`.
These comparisons are monotone because `3 gamma>1`; exponential
bracketing then binary search suffices.

The **sign** identifies comparability between distinct endpoints. In
an ancestor pair, the descendant lies in a first outer subinterval,
hence its midpoint differs from the ancestor midpoint; the horizontal
component makes its distance strictly exceed the vertical distance
\(D\ln3\). In an incomparable pair, \(Q<1\) makes the inequality
strict in the opposite direction. Thus once \(D\) is decoded, a further
rational comparison with \((3^D+3^{-D})/2\) decides whether the pair is
comparable. If root distances are supplied too, they orient a comparable
pair and give the least common ancestor's depth by
\((D(r,v)+D(r,w)-D(v,w))/2\).

There is a uniform noise margin too. Let μ be the midpoint of the
signed error interval. Rounding `(d_H−μ)/ln 3` works even if the
measured hyperbolic distance differs from the exact distance by less
than \(	frac12\ln(33708/32805)\).

In this extreme-subinterval scheme, base three is the smallest
integer base that admits a single fixed shifted-floor correction
for the weighted tree. At base two the minimum gap is 7/16 at the
(1,4) split. The same limiting configurations give signed errors
approaching \(2\ln(7/16)\) and \(\ln(5/4)\); their separation
\(\ln(320/49)\) is greater than \(\ln 2\), so no offset can put
them both into one base-two logarithmic bin.

The sharp lower bound is tied to the exact length grammar, not
merely the alphabet `{1,2,3,4}`: an allowed `(1,1)` split would
give a much smaller gap. If instead one wants all original vertices
before compression, use outer **fifths** per ordinary edge. The
same argument with `g=3/5` gives the error interval
`[2ln(3/5), ln(5/4)]`, of width `ln(125/36)<ln5`. The
unweighted distance is then
`floor((d_H+2ln(5/3))/ln5)` on this different embedding.
