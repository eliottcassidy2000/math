# Small rational data for JC(2): orbit height, coefficient varieties, and elliptic slices

**Status.**  The THM-1300 coefficient/gauge calculations and the finite Gauss
chart censuses below are **VERIFIED-EXACT / FINITE-EXACT**.  The planar degree,
factor, equivariance, and finite-group exclusions are **PROVED** or **CITED**
with their scopes stated below.  The proposed elliptic-slice search is a
**HEURISTIC SEARCH PROGRAM**.  It supplies neither a planar counterexample nor
a proof of `JC(2)`.

The exact companion is
[`jc2_small_rational_height_and_elliptic_slice_prior.py`](../04-computation/jc2_small_rational_height_and_elliptic_slice_prior.py),
with frozen
[`output`](../05-knowledge/results/jc2_small_rational_height_and_elliptic_slice_prior.out).
It replays identically under ordinary and optimized Python.

```text
script sha256: bd8c349d6fa6697b73a1af4f1962831b1349d6122c4c55264d50f41cd16d2f9c
output sha256: bfe629a238db327bc593e7abc94e32ba01b306f73ee0da7847566fa888ff2c9a
hash basis: raw LF bytes
```

## 1. The correction that organizes the search

The useful premise is not that a counterexample must have small expanded
coefficients.  It is this more careful statement:

> A high-degree exceptional map can have a small description in the correct
> generative language, even when expansion, normalization, or a rational-point
> parametrization makes its displayed coefficients enormous.

Three different heights must therefore be kept separate.

1. For a reduced rational `r=a/b`, put
   `H_Q(r)=max(|a|,b)`.  The **displayed height** of a polynomial tuple is the
   maximum `H_Q` of its expanded nonzero coefficients.
2. The **generator height** is the height of the parameters in the structural
   grammar that produced the tuple: factors, shear potentials, exponents,
   conductor corrections, or a point on a coefficient variety.
3. After fixing the origin and identity linear jet, the **diagonal-orbit
   height** is the minimum displayed height under rational diagonal
   conjugacy.

These are genuinely different.  The three-variable witness has displayed
height `12` in its natural identity-jet gauge and diagonal-orbit height `9`.
Forcing a particular rational collision to `0,(1,1,1)` raises its displayed
height to `945` and its support from `16` to `40` terms.  Thus collision
normalization is algebraically lawful but arithmetically hostile.

There is also a naming correction about the user-supplied elliptic object.
The link is to **curve number 273**, not to an elliptic curve of rank 273.  The
current primary database field is a **rank lower bound of 30**.  Nothing here
asserts exact rank `30`, much less rank `273`.

## 2. Inheritance pass and concept board

The closest proved mechanism is
[THM-3587](../01-canon/theorems/THM-3587-consecutive-keller-fibre-factor-toric-and-coefficient-span-gates.md):
write the derivative table as an `SL_2` plaquette, pay the determinant by a
Gauss word, and leave the two curls as the real integrability equations.  The
canonical hostile is the determinant-one matrix

```text
[xy+1  xy+2]
[xy    xy+1],
```

whose two curl defects are both `x-y`.

The corrected near miss is “normalize everything to `0` and `1`, then search
small coefficients.”  The exact height-`945` collision gauge shows why that
can discard the small representative of an orbit.  The least-used sidecar is
the character lattice of residual diagonal conjugacy: Smith reduction of
coefficient weights says exactly which coefficients can simultaneously be
set to `1` over `Q` and which square/power classes remain.

The active board is:

| object | representation | preserved predicate | missing coordinate |
|---|---|---|---|
| THM-1300 | seven-parameter factored support | Keller determinant and collision | none in its nondegenerate chart; all parameters reduce to gauge |
| planar candidate | non-elementary Cohn core with small elementary decorations | determinant and retained `SL_2/E_2` class | both curls, weighted transport holonomy, collision, nonproperness |
| coefficient locus | affine scheme over `Q` | all coefficient equations imposed in its ideal | global map semantics and excluded components |
| curve slice | normalization of a one-dimensional component | rational solutions of the sliced equations | points on other components and outside the slice |
| height quotient | diagonal character lattice | rational gauge orbit | nonlinear affine equivalences |
| symmetry hostile | exact torus/even finite-group action | equivariance | the asymmetric defect required by a counterexample |

## 3. Why the three-variable rational search was actually small

Put the THM-1300 map in tangent-to-identity order

```text
G=(F_3/2,F_2,F_1).
```

It has component term counts `(3,6,7)`, degrees `(4,6,7)`, determinant one,
displayed height `12`, and maximum denominator `2`.  Conjugate by
`D=diag(a,b,c)`.  Equivariance makes the result depend only on

```text
s=ab,                 t=a^2 c.
```

The complete coefficient list, including the three linear ones, is

```text
G_1:  1, -3s/2, -t/2,
G_2:  1, 3t/s, 6t, 3st, 12s, 9s^2,
G_3:  1, 3s, 3s^2, s^3, 4s^2/t, 7s^3/t, 3s^4/t.       (1)
```

This gives an elementary global height minimum.  If every coefficient had
height at most eight, `H_Q(s^3)<=8` would force the reduced numerator and
denominator of `s` to be at most two.  The coefficient `12s` then forces
`|s|=1/2`, but `9s^2=9/4` has height nine.  Conversely, exact bounded
enumeration gives all height-nine gauges:

```text
s=+/-1/2,       t=+/-1/4, +/-1/2, +/-3/4.              (2)
```

Thus the diagonal-orbit height is exactly `9`.  The secondary representative
`s=t=1/2` has maximum denominator `8` and the particularly small collision

```text
(0,0,-1/2), (1,-3,13), (-1,3,13)  ->  (0,0,-1/2).      (3)
```

The more important compression appears before expansion.  Consider the
seven-parameter torus-compatible grammar

```text
u=1+p xy,                 w=q+rxy,
G_1=x+C x^2y+D x^3z,
G_2=y+A x u^2z+Bxy^2w,
G_3=u^3z+y^2uw.                                           (4)
```

The equation `det DG=1` has twelve coefficient equations.  Exact lexicographic
Groebner elimination produces, among its ideal certificates,

```text
p^4(A-B),          p^4(Aq-12p),          p^3(Ar-9p^2),
-2Aq+2Bq+2C+3p,
-2A^2q+2ABq+Ap+6D.                                      (5)
```

On the nondegenerate branch `p!=0`, these force `A!=0` and

```text
B=A,       q=12p/A,       r=9p^2/A,
C=-3p/2,   D=-Ap/6.                                      (6)
```

Substitution verifies sufficiency.  Moreover `(6)` is exactly `(1)` under

```text
s=p,                   t=Ap/3.                           (7)
```

So the seven apparent parameters leave only two, and those two are entirely
diagonal gauge.  This is the clean version of “a small exhaustive search
could have found it”: the support, equivariance, and factor sharing do almost
all of the work.  A raw thirteen-coefficient search would not.

## 4. Curve #273: what is verified and what it teaches

The current ICARM database entry for
[curve #273](https://elliptic-rank.icarm.cloud/curve/273), exposed exactly by
its [public JSON API](https://elliptic-rank.icarm.cloud/curve/273.json), is

```text
y^2+xy=x^3+A x+B,
A=-201769035260418549083594900060734240952308696994802735114305555,
B=1151107939141058565733479426024323225135665982951300586808823640527729578307228357301072889377.
                                                                    (8)
```

As of the retrieved record (`updated_at=2026-08-20 11:24:48`), the database
lists `rank_lower_bound=30` and thirty witness points.  The database's
[API description](https://elliptic-rank.icarm.cloud/api) says submitted points
are checked exactly and independence is certified by exact 2-descent
quadratic characters.  That is a **CITED CURRENT DATABASE CERTIFICATE**, not
an independence calculation reproduced here.

The companion independently verifies only the following finite arithmetic:

- all `30/30` stored rational points satisfy `(8)` exactly;
- the frozen minimal discriminant agrees with the Weierstrass formula;
- `19` points are integral and `11` are nonintegral;
- every denominator has the standard form
  `den(x)=d^2, den(y)=d^3`, with
  `d in {1,2,3,5,13,23,37,97,1327,19963}`;
- the `a_4,a_6` values have `63,94` decimal digits, while the raw affine
  heights of the witness points have `46..60` digits.

The page commentary mentions an exact-rank conclusion conditional on
additional conjectures.  It is **not imported** here.  We assert neither an
unconditional upper bound nor exact rank.

The connection to JC(2) is methodological and must be typed:

```text
source: a rational/elliptic curve component of a coefficient scheme
map: rational reconstruction of the polynomial-map coefficients
preserved: the coefficient identities used to define that scheme
lost: unsliced components, excluded denominators, global collision/properness
sidecar: exact Keller/factor/span/symmetry/collision filters
test: substitute each reconstructed point and run the full exact passport.
```

The lesson of `(8)` is not “large rank guarantees a counterexample,” nor even
“the generators are small.”  The current entry supplies no small-generator
claim.  The lesson is that group structure and arithmetic geometry can make
rational points accessible even when the displayed equation and point
coordinates are far outside any naive coefficient box.

## 5. What can legitimately be normalized in a rational JC(2) search

Let `F=(P,Q) in Q[x,y]^2` have nonzero constant Jacobian.  Rational affine
changes always allow

```text
F(0)=0,                 DF(0)=I,                 Jac(F)=1. (9)
```

Constants become zero and the four linear coefficients become `0/1`.
After `(9)`, diagonal conjugacy acts on a coefficient of `x^i y^j` by

```text
P coefficient: c -> c a^(i-1)b^j,
Q coefficient: c -> c a^i b^(j-1).                     (10)
```

Thus two pivot coefficients can be normalized only when their two character
vectors permit it over `Q`.  A unimodular `2x2` character matrix permits
rational normalization to chosen signs.  Determinant of absolute value
greater than one leaves a genuine perfect-power or squareclass obstruction;
rank one leaves an invariant ratio.  Smith normal form, rather than informal
rescaling, is the correct quotient.

If a rational collision pair is already known, one can translate one point
to zero, use the derivative there to restore `(9)`, and use a residual
diagonal conjugacy to put the other point at `(1,1)` or on a coordinate axis.
But Section 3's height-`945` hostile shows that this must be a late semantic
normalization, not the coefficient-height gauge used to enumerate candidates.

Leading-form normalization is also conditional.  In a reduced pencil basis,
the common leading base may be made primitive and one of its coefficients may
be set to one.  A rational source `PGL_2` change can put two or three
**rational** roots at `0,infinity,1`; it cannot do this to nonrational Galois
orbits.  Split and nonsplit leading bases must be searched separately.

## 6. Why a generic low-degree or dense small-height search is the wrong universe

[THM-3550](../01-canon/theorems/THM-3550-prime-degree-exclusion-and-pencil-height-eight-floor.md)
proves internally that reduced pencil degrees are distinct composites with
lower degree at least six and height at least eight.  Its cited comparison
records the stronger imported gates: Nagata forces at least three prime
factors with multiplicity in every pencil degree, and the cited sub-125
classification gives reduced height at least `108`, with only `(72,108)`
below `125`.  These are degree bounds, not coefficient-height bounds.

Even the first cited degree cell is computationally absurd when expanded
densely.  A bivariate polynomial of degree at most `d` has
`(d+1)(d+2)/2` slots.  In a tangent-identity reduced `(72,108)` pair, fixing
the constants and linear jet leaves

```text
[(73*74)/2-3]+[(109*110)/2-3]=8690                     (11)
```

rational slots.  The rational height box

```text
R_B={a/b: gcd(a,b)=1, b>0, max(|a|,b)<=B}
```

has exact size

```text
|R_B|=4 sum_(n<=B) phi(n)-1.                            (12)
```

For `B=1,2,9` the dense universes in `(11)` have respectively `4147`,
`7344`, and `17774` decimal digits.  Small rational height does not rescue a
dense search.

By contrast, at `B=9`, fixed structural templates with `k=2,3,4,5`
parameters have exact sizes

```text
12,321;  1,367,631;  151,807,041;  16,850,558,151.      (13)
```

Three parameters are directly enumerable.  Four become realistic after
modular and symbolic pruning.  Five require that equations eliminate at
least one parameter before enumeration.

## 7. The rational/elliptic coefficient-variety program

The search should be an arithmetic-geometry search over **generators**, in
this order.

### 7.1 Quotient the full elementary class, then search a non-elementary core

The three-element Gauss chart was the first full-span generator:

```text
M=E_+(V)E_-(U)E_+(W)
 =[1+UV, V+W(1+UV); U,1+UW].                            (14)
```

It pays `det M=1` and the balanced matching-factor passport automatically.
Only the two curl equations remain:

```text
U_y=(UW)_x,
(1+UV)_y=[V+W(1+UV)]_x.                                (15)
```

It is no longer a live counterexample generator.  The
[four-Gauss search reflection](jc2-small-rational-four-gauss-search-boxeph-20260821.md)
supplies an all-degree characteristic-zero closure of this three-word chart.
This is **SESSION-PROVED + INDEPENDENTLY HOSTILE-AUDITED / PENDING CANON
PROMOTION**, so it is not used here as a canon dependency, but its elementary
mechanism is explicit.  If both rows are closed, either `U=0` and the map is
triangular, or there are univariate polynomials `A,Phi,Psi` such that

```text
W=A'(y),       s=x+A(y),       U=Phi'(s),
q=y+Phi(s),    V=Psi'(q),
F=(s+Psi(q),q).                                           (16)
```

The inverse is the sequential subtraction

```text
s=P-Psi(Q),       y=Q-Phi(s),       x=s-A(y).             (17)
```

Thus the known balanced span-four control is not merely one tame point: it
belongs to the all-degree sequential-shear solution class.

The companion exhausts all `3^9=19,683` affine `U,V,W` coefficient tuples in
`{-1,0,1}`.  Exactly `375` satisfy both curls, by two independent exact
implementations.  Every survivor has component degree at most two and
coefficient span at most one.  This is a useful hostile: the first meaningful
search must give at least one of `U,V,W` nonlinear and must not merely enlarge
the affine coefficient box.

The known three-shear map is a positive control: it solves `(15)`, has
coefficient span four and balanced `(2,2)` matching fibres, but has an
explicit polynomial inverse.  Any pipeline that discards it before the final
tameness test has over-pruned.

Wright's classical weak Jacobian theorem, independently replayed in the same
reflection, proves substantially more: every planar Keller map whose
normalized Jacobian lies in `E_2(k[x,y])` is tame, regardless of elementary
word length.  The reduced-word proof uses leading monomials, source-shear
cancellation, and an exact Bruhat shortening of internal constants.  Thus
neither a larger rational box nor a more complicated elementary support can
rescue

```text
M_5=E_+(R)E_-(Z)E_+(V)E_-(U)E_+(W).                    (18)
```

The [finite-exact companion](../04-computation/jc2_small_rational_four_gauss_search.py)
is retained as a replay/control: its `20,400` four-word first-curl systems
leave four nonzero one-dimensional families, all in the proved tame stratum,
and its `1,000` signed-
monomial five-word first corrections all fail.

The six-factor staged system

```text
M_6=E_-(S)E_+(R)E_-(Z)E_+(V)E_-(U)E_+(W).              (19)
```

For fixed `U,V,W,Z`, its first correction is linear in all coefficients of
`R`; only after that closes the top row is the second equation linear in `S`.
was the first finite-depth attempt with two nonzero staged defects.  The independent
[six-word companion](../04-computation/jc2_small_rational_six_gauss_search.py)
tests the first exact box: `U,V,W,Z` independently in
`+/-{x,y,x^2,xy,y^2}`, and the entire total-degree-at-most-six `R` space.  All
`10,000` systems have modular ranks `(28,29)` with maximal coefficient rank;
the added two-scale `W=y^2+/-x`, degree-eight correction box eliminates all
`2,000` systems with ranks `(45,46)`.  These are exact bounded hostiles inside
Wright's all-length tame class, not counterexample frontiers.

The first live structural generator is instead Cohn's certified
non-elementary matrix

```text
C=[1+xy x^2;-y^2 1-xy].                                (19a)
```

Its first repair equation forces the factorial tail
`a_i=(-1)^i/(i+2)!` and has no polynomial solution.  This is precisely where
small rational arithmetic remains useful: decorate `(19a)` on both sides,
build the monomial transport graph, and retain only cycles whose exact
multiplier product satisfies the holonomy equation.  Support parity without
edge gains is insufficient.  The full router is
[`jc2-cohn-parity-cycle-repair-codex-20260821.md`](jc2-cohn-parity-cycle-repair-codex-20260821.md).

Other typed generators worth slicing are:

- the rank-three trace/Hessian chart from THM-3587, with its stronger
  five-factor floor;
- reduced `(72,108)` leading bases
  `R_72=aH^(72/h), S_108=bH^(108/h)` for
  `h in {2,3,4,6,9,12,18,36}`, followed by sparse repair bands;
- a near-involution chart with one explicit asymmetric defect.  Exact
  reflection symmetry is not a live counterexample lane.

### 7.2 Construct and quotient the coefficient scheme

For fixed supports, form the coefficient ideal over `Q` from the curls,
remaining determinant equations, prescribed leading base, and any chosen
collision equations.  Then:

1. eliminate coefficients that appear linearly;
2. saturate every denominator and nonzero pivot used in the chart;
3. quotient the diagonal character lattice using Smith normal form;
4. primary-decompose when affordable, retaining the exclusions for every
   component;
5. compute dimensions before enumerating rational values.

Zero-dimensional components are solved exactly.  For positive-dimensional
components, systematic low-height affine planes give rational curves.  This
is a slice census, never an exhaustion of the parent component.

### 7.3 Route curves by genus

- A genus-zero component with a rational point is parametrized and searched in
  parameter height, not expanded coefficient height.
- A genus-one component with a rational point is converted to a Weierstrass
  model.  Descent, certified independent points, local information, and the
  Mordell--Weil group law then generate rational candidates far outside a
  rectangular coefficient box.
- Higher-genus components need their own arithmetic tool; a failed bounded
  point search is not an emptiness proof.

An elliptic search should enumerate bounded integer combinations of
**certified** independent points and map them back through the reconstruction
map.  This is an exact finite experiment.  It is not exhaustive on `E(Q)`
unless a full Mordell--Weil basis and a rank upper bound have also been proved.
Curve #273 illustrates why the lower-bound/exact-rank distinction must remain
visible.

Modular point counts, root numbers, and Nagao-style scores can prioritize
curve slices, but they are heuristics.  Every retained rational point must be
reconstructed and checked over `Q`.

### 7.4 Run the exact passport in the cheapest order

For each reconstructed pair:

1. verify both curls and `Jac(P,Q)=1` identically;
2. verify the intended degree and leading-form/Newton architecture;
3. compute the nonconstant coefficient-matrix span;
4. test the four cross ideals of `P_x,P_y,Q_x,Q_y`;
5. factor `T=P_yQ_x` and `T+1=P_xQ_y`, retaining multiplicities and distinct
   factors separately;
6. apply the scalar- and toric-carrier exclusions of THM-3587;
7. compute the exact symmetry stabilizer and reject inherited equivariant
   no-go lanes;
8. only then solve for a genuine collision, saturating away the diagonal;
9. test nonproperness and the independent every-line/fibre gates.

A rational collision is a decisive nonautomorphism certificate.  Absence of a
rational collision is not evidence of injectivity over `C`.

## 8. The symmetry boundary: fixed-line reflections and even commuting actions

[THM-1345](../01-canon/theorems/THM-1345-jc2-equivariant-category-poisson-reframing-dc1-shadow.md)
proves that planar Keller maps equivariant for a `C*` action are automorphisms.
Thus the exact mixed-weight symmetry that compresses THM-1300 to `(4)` cannot
be the planar generator.

Fixed-line reflections admit a direct reduction to the injective-line gate.
Suppose affine reflections `sigma,tau` have fixed affine lines `L,M` and a
planar Keller map satisfies

```text
F sigma = tau F.                                        (19)
```

For `p in L`, equation `(19)` gives `F(p) in M`, hence `F(L) subset M`.
If `v` is a nonzero tangent vector to `L`, then the derivative of the
restriction is `DF_p(v)`.  The Keller condition makes `DF_p` invertible, so
this vector is never zero.  In affine coordinates on `L` and `M`, the
restriction is therefore a one-variable complex polynomial whose derivative
has no zero.  Its derivative must be a nonzero constant, so `F|_L` is
injective.  The cited
[Gwozdziewicz injective-line theorem](https://arxiv.org/abs/alg-geom/9305008)
then makes `F` an automorphism.  This proves the reflection-intertwining
exclusion even when the source and target reflections are different; it does
not rely on classifying all plane involutions.

The fixed-point-only central involution belongs to a different cited gate.
Miyanishi's
[*Equivariant Jacobian Conjecture in dimension two*](https://arxiv.org/abs/2110.06709)
proves that an etale endomorphism of the complex affine plane commuting with
an effective finite group of even order is an automorphism.  Thus `-I`
equivariance, and more generally an effective even-order **commuting** action,
is fatal.  This statement does not cover arbitrary distinct source/target
central intertwiners, ineffective decorations, or symmetries that do not
commute with the map.  A rational search should therefore either:

- impose a computable certificate of trivial/effectively odd symmetry; or
- begin from a symmetric hostile and add an explicit symmetry-breaking layer
  before treating the object as a candidate.

This is where the two-dimensional program must differ from the
three-dimensional one.  Reusing small constants is reasonable; reusing its
exact torus or involution is not.

## 9. Search schedule

The recommended first exhaustive campaign is:

1. reject the entire elementary class by Wright's criterion, recording a
   factorization when available as a tame certificate;
2. enumerate elementary--Cohn--elementary double-coset supports with at most
   three residual rational parameters after curl elimination and Smith
   quotient, through generator heights `B=1,2,4,9,12`;
3. reject exposed factorial paths and cycles with wrong multiplier holonomy;
   attack four-parameter reciprocal-cycle cells with modular rejection before
   rational reconstruction;
4. slice every positive-dimensional survivor to genus-zero/genus-one curves;
5. run bounded Mordell--Weil combinations on certified elliptic slices;
6. only after a candidate passes span, factor, carrier, and symmetry gates,
   spend on collision and nonproperness elimination.

The first target is not “find the counterexample.”  It is a much sharper
finite object:

```text
a full-span, non-equivariant, row-integrable matrix in a certified nonzero
`SL_2/E_2` class, together with a rational or elliptic coefficient component
and one point passing every cheap gate.                              (20)
```

That outcome would justify raising heights or supports.  Failure should be
recorded per support component and slice, never as evidence for `JC(2)`.

## 10. Reproduction and scope ledger

Run

```bash
python3 04-computation/jc2_small_rational_height_and_elliptic_slice_prior.py
python3 -O 04-computation/jc2_small_rational_height_and_elliptic_slice_prior.py
```

The companion checks the determinant, identity jet, complete height-eight and
height-nine diagonal gauges, the twelve determinant equations and selected
Groebner ideal certificates, the exact two-parameter orbit, its collision,
rational-box sizes, the complete affine Gauss `B=1` census by two paths, a
balanced full-span tame control, the determinant-only curl hostile, all thirty
ICARM point incidences, denominator shapes, the minimal discriminant, and the
height-`945` normalization hostile.

It does **not** verify the ICARM independence certificate, prove a rank upper
bound, classify rational points on curve #273, exhaust any elliptic slice,
prove that a JC(2) coefficient variety has a genus-one component, or settle
`JC(2)`.
