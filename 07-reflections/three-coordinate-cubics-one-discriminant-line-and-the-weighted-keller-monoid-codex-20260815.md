# Three coordinate cubics, one discriminant line, and the weighted Keller monoid

**Research synthesis / provenance, not a new theorem.**  The fixed-map facts
are in THM-2473, THM-2546, and THM-2582; the all-degree construction and grade
classification are in THM-3438.  This note separates those proved scopes.

## The three cubics are three coordinate views of one cubic extension

For the fixed sporadic degree-three Keller map, the generic fibre can be
described by cubics in the `x`, `r=b-y`, and `z` coordinates.  THM-2546 gives

```text
Disc_x=-4 S_x^2 L,
Disc_r=-4 S_r^2 L,
Disc_z=-4 S_z^2 L.                                    (1)
```

This is not three independent discriminant conditions.  The three source
coordinates lie in the same rank-three function-field extension, and their
global coordinate equations are cubic.  The exact formulas above establish
the common square class without requiring a new irreducibility claim for each
coordinate equation.  On any locus where `alpha,beta` are primitive elements,
the change between the power bases

```text
(1,alpha,alpha^2)  and  (1,beta,beta^2)
```

has determinant `I_(alpha,beta)`, and the trace-form discriminants obey

```text
Disc(beta)=I_(alpha,beta)^2 Disc(alpha).                (2)
```

The three square factors in `(1)` have the expected coordinate-projection /
index grammar.  Equation `(2)` explains that grammar where the chosen
coordinates generate; equation `(1)` is the proved input.  The invariant
content is the common square class

```text
[Disc]=[-L].                                           (3)
```

Even the visible `4` is a square factor; it does not supply a fourth sheet or
a degree-four map.  What is synchronized is one branch/Jelonek divisor and
one cubic `S_3` extension, seen through three coordinate projections.

## Why this does not classify the infinite family

The current repository classification has several levels which must not be
collapsed:

| object | degree / monodromy | discriminant carrier | what is classified |
|---|---|---|---|
| sporadic `F` | `3`, geometric `S_3` | three cubics with square class `[-L]` | one fixed map's fibre and boundary anatomy |
| `F o F` | `9`, imprimitive `S_3 wr S_3` | odd square-free divisor changes to `H` | one composite level; cubic law does not persist |
| weighted atom `F_n` | every `n>=3`, geometric `S_n` | `Disc_w(R(w)-Pw+cQ)` | one explicit composition atom in every grade |
| composites | degrees `ab`, `a,b>=3` | norm and cross-block discriminants | existence of decomposable maps in product grades |
| stable extensions | same degree in every dimension `m>=3` | same function extension plus identity variables | degree-value spectrum in higher dimensions |

THM-3438 proves the exact numerical spectrum

```text
KDeg(m)=N_{>=1} minus {2},             m>=3,            (4)
```

and the existence spectra

```text
AtomDeg(m)={n:n>=3},
DecompDeg(m)={ab:a,b>=3}.                              (5)
```

This is already an infinite all-dimension family and an exact classification
of degree values.  It is not a classification of maps up to tame conjugacy,
stable equivalence, or composition.  Mixed grades contain both the primitive
`S_n` weighted atom and a composition.  Distinguishing those maps needs the
inverse polynomial, monodromy block lattice, Jelonek tower, and oriented
boundary-effectivity data.

The degree-nine composite is the decisive hostile to extrapolating `(1)`:
THM-2582 shows that norm and cross-block terms change the odd discriminant
square class from `L` to `H`.  Thus “all three cubics have discriminant
`-4(square)^2L`” is a rigid coordinate theorem for the sporadic core, not the
generator relation of the degree-graded monoid.

## Tournament and XOR reading

Take the three coordinate projections as vertices and ask whether a pair has
the same discriminant divisor.  Every pair answers yes.  The resulting graph
is complete and unoriented; turning it into a tournament adds a gauge but no
information.  It forgets the individual index squares and, more seriously,
the signed valuations deciding whether a rational coordinate view is
polynomial at the boundary.

The faithful object is a small effectivity complex:

```text
generic cubic field / common branch divisor
        + three coordinate-projection / index squares
        + boundary valuation vectors.                 (6)
```

This exactly parallels the mask lessons on the LRC side.  XOR or ternary
intersection detects symmetric overlap.  It does not recover orientation,
endpoint scale, or current.  Here the missing sidecar is boundary
effectivity; there it is owner/mode/current data.

## The two appearances of discriminant `-4` stay separate

The Berggren U-spine quadratic

```text
c_t=2t^2+2t+1
```

has polynomial discriminant `-4`, encoding a Gaussian/parabolic quadratic
carrier.  The sporadic cubic has discriminant `-4 S^2 L`, whose nonsquare
content is the variable divisor `L` and whose monodromy is `S_3`.  Their
common printed factor is not a common cover, field, or cohomology class.

The cheapest lawful comparison is representation-theoretic: both calculations
separate a square index factor from a nonsquare carrier.  The carriers differ
(`-1` for the quadratic arithmetic, `-L` for the cubic function field), so no
transplant follows.

## Exact test executed provisionally

The reserved
[THM-3494](../01-canon/theorems/THM-3494-weighted-lift-primitive-coordinate-discriminant-atlas.md)
proof package now executes the proposed test, pending independent audit.  For
the canonical weighted seeds in degrees `3,4,5`, it computes and factors

```text
Disc_w(R_n(w)-Pw+Q),                                   (7)
```

and computes the actual `x`- and `y`-reconstructing coordinate eliminants in
each degree.  Their discriminant ratios are exact nonzero squares over
`Q(P,Q,C)`.  The degree-three row also checks `z`, giving three independent
cubic views.  The all-degree explanation is the trace-form change-of-basis
identity plus THM-3438's maximal point stabilizer: `S_(n-1)<S_n` has no
intermediate subgroup, so every nonconstant `x/y` view is primitive.

The result separates:

1. the reduced branch divisor;
2. the primitive-element index square;
3. the Jelonek/nonproper divisor;
4. the boundary valuation cone; and
5. the monodromy block lattice.

Items 1 and 2 are now exact on the generic chart.  Items 3 and 4 remain open
globally across `C=0`; item 5 is inherited from THM-3438.  A further exact
consequence is that four primitive views give six `K_4` edge indices forming
one multiplicative coboundary.  Tournament orientation is therefore gauge,
and XOR/square-class reduction loses the divisor-valued index sidecar.  This
turns the three-cubic observation into an all-degree primitive-element atlas
without pretending it classifies arbitrary Keller maps.  No LRC or physical
tournament consequence follows merely from the factor `-4`.
