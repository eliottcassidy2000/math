# The carrier grammar is now an ordinal task generator, not a conjectural closure

**Status (2026-08-24): SYNTHESIS OF PROVED + VERIFIED-EXACT RESULTS; JC(2)
REMAINS OPEN.**  The proof sources are
[THM-3920](../01-canon/theorems/THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction.md),
[THM-3929](../01-canon/theorems/THM-3929-root-regular-one-place-linear-color-cubic-is-monogenic.md),
[THM-3933](../01-canon/theorems/THM-3933-centered-degree-three-root-map-pole-partition-octic-nonentry.md),
[THM-3936](../01-canon/theorems/THM-3936-centered-degree-three-infinite-root-value-nonentry.md),
[THM-3938](../01-canon/theorems/THM-3938-centered-degree-four-root-map-nonentry.md),
and
[THM-3941](../01-canon/theorems/THM-3941-all-degree-centered-cubic-pole-carrier-routing.md).
This reflection records the research move and the next generated obligations;
it is not an additional theorem.

## 1. Inheritance pass and concept board

The anchor was the next centered one-place linear-color stratum after the
degree-three theorems.  The closest proved mechanism was THM-3920: a genuine
ramification component of a finite-flat cubic completion cannot coalesce two
normalization addresses at one point of an affine-plane boundary curve,
because that boundary curve is unibranch.  The canonical hostile example was
the root-regular exit, which is not a ramification fold at all and is instead
closed by the scalar-unit/monogenic gate of THM-3929.  The corrected near miss
was any argument that read multiplicity in the discriminant of an auxiliary
binary-cubic order before proving that the order was maximal.  The least-used
sidecar was the completed-local trace character at each pole.

The live board was:

1. completed-local trace and the `C2`/`C3` inertia characters;
2. the finite Riemann--Hurwitz budget of a cubic polynomial projection;
3. exact color divisibility after the pole carrier is fixed;
4. exponent-one discriminant factors versus square index factors;
5. normalization-address multiplicity and boundary unibranchness;
6. natural-number task addresses rather than oversized arithmetic labels;
7. elliptic-resolvent character rank as an orthogonal wildcard.

## 2. Independent convergence at degree four

The degree-four calculation was performed independently in several lanes.
All lanes obtained the same five collision-free non-root-regular pole rows,
the same single color family in each row, the same two scalar seams, and the
same two- or three-address obstruction for every non-scalar survivor.  While
that work was being prepared, an incoming session promoted the identical
result as THM-3938.  This was useful signal rather than duplicated canon: the
independent reconstructions became hostile audits of the incoming theorem.

One audit found a real proof obligation.  At the scalar seams, an index-form
value cannot be fed directly to THM-3801 unless the displayed binary-cubic
order is the maximal normalization.  The repaired proof computes

```text
Disc(O_E)=-2J_E,                 Disc(O_F)=-16J_F,
```

with each `J` irreducible.  Exponent one forbids a nontrivial codimension-one
index at the sole discriminant prime; away from it the order is etale.
Finite freeness supplies `S2`, and the height-one check supplies `R1`, so the
order is normal.  Only after that repair do the scalar values close by
THM-3801.  This is the reusable lesson: **separate order discriminant,
maximal-order discriminant, and geometric branch before using any one of
them as the other.**

## 3. The all-degree map

THM-3941 abstracts the degree-three and degree-four ledgers without claiming
all-degree color closure.  Its source and target are:

```text
trace-zero rational root map t(u)
    | forget lower Laurent jets, retain projection/support/orders
    v
(infinity order m, finite inertia carrier, finite pole orders)
    | completed-local trace + RH + shared-address routing
    v
shared A-address, C3, selected C2, C2 x C2, or the root-regular exit.
```

The preserved predicates are total degree, pole support, pole orders, local
trace character, and whether two poles share an `A`-fibre.  The map destroys
lower principal parts, polynomial-color coefficients, coefficient-ideal
data, maximal-order indices, and the actual target-address fibres.  Those
destroyed coordinates are exactly the sidecars required by the next stage.

For a collision-free non-root-regular degree `N` row the possibilities are:

- one pole of order nonzero modulo three at the `C3` point;
- one odd pole at a selected `C2` point; or
- two odd poles, one at each `C2` point;

together with infinity order zero or positive nonzero modulo three, and with
all orders summing to `N`.  Each such support/order tuple is nonempty at the
trace layer.  The `C3` block is a monomial.  At a `C2` point, the recurrence
for `A=u^3+u^2`, `v=1/u`,

```text
A v^3-v-1=0,
S_1=0, S_2=2/A, S_j=(S_(j-2)+S_(j-3))/A,
```

gives a triangular even-power correction to every odd leading pole.  Its
determinant is `2^n`, so the trace-zero block exists uniquely.  Adding blocks
over the two critical values and a trace-corrected polynomial infinity block
realizes every carrier.  This proves trace realizability, not color survival.
If two finite poles share an `A`-fibre, the bare equations force only their
common target address `a(A0)=C=0`; source non-unibranchness additionally needs
the genuine exact-double, reduced tame `(2,1)`, maximal-normal Keller
hypotheses. MISTAKE-470 and its infinite counterfamily guard this quantifier.

## 4. Natural-number addresses and the period-twelve law

Let `c_N` be the number of non-root-regular carriers.  The exact generating
function is

```text
M(R3+R2+P2),
M=1+(x+x^2)/(1-x^3),
R3=(x+x^2)/(1-x^3),
R2=x/(1-x^2),
P2=x^2/((1-x^2)(1-x^4)).
```

Writing `N=12j+r`, `0<=r<12`, gives the exact identity

```text
c_(12j+r)=6j^2+b_r j+d_r,
b=(16,10,14,16,16,14,22,16,20,22,22,20),
d=( 0, 2, 4, 5, 5, 7,10, 9,12,14,15,16).
```

This is the useful version of the user's ordinal principle.  A large odd
square, a triangular value, or a complicated carrier tuple may certify set
membership, but once the carrier has been selected its identity is simply
its rank.  At fixed `N`, order first by infinity order, then
`C3<C2<C2xC2`, then by the smaller odd part.  The carriers are exactly the
naturals `0,...,c_N-1`.  The rank is a procedural address, not a geometric
invariant, and overlaps between other encodings do no harm.

The companion proves the rational-function identity itself, not merely a
finite fit, and checks the deterministic ordinal bijection through degree
thirty.  An independent implementation agreed through degree 360.

## 5. The seven degree-five tasks

Degree five has the following zero-based work queue:

```text
0: m=0, C3(5)          1: m=0, C2(5)
2: m=1, C3(4)          3: m=1, C2xC2(1,3)
4: m=2, C2(3)
5: m=4, C3(1)          6: m=4, C2(1).
```

The separate `m=5` polynomial row is the root-regular THM-3929 exit.  For
each of tasks `0,...,6`, the generated obligation is now uniform:

1. solve the lower-principal-part trace equations;
2. divide the exact color and retain only polynomial families;
3. test the coefficient ideal for the unit ideal;
4. distinguish maximal-order ramification from square index factors; and
5. compute the normalization-address fibres of every genuine component.

This is a finite seven-row next frontier, but it is not yet seven solved
rows.  The cheapest decisive experiment is exact symbolic color division in
each row, with the inherited degree-three and degree-four survivors as
positive controls and a forbidden shared-fibre pole packet as hostile
control.

## 6. What the tournament inheritance does and does not say

There is no intrinsic tournament on these carriers.  The vertices would be
support/order tasks, but there is no canonical binary observable orienting
every pair, ties would be structural, and forgetting the trace/color
sidecars would destroy the target predicate.  Forcing an orientation would
therefore be cosmetic.

The useful inherited tournament idea is instead **common-target confluence**:
distinct normalization addresses that a finite-flat map sends to the same
target address must be reconciled on the source boundary.  In the genuine
ramification case, rank three coalesces them at one source point, and
unibranchness supplies the obstruction.  Thus the remembered tournament
lesson is to preserve the pairwise observable and its sidecar; here it leads
to an address-fibre incidence structure, not to a tournament.

## 7. Orthogonal signal and stopping boundary

THM-3937, THM-3939, and THM-3940 show that the tested two-boundary elliptic
resolvents have only the Cardano cubic character, even when the free
Mordell--Weil rank is raised.  Incoming results then tested the obvious
response rather than leaving it as advice.

[THM-3942](../01-canon/theorems/THM-3942-affine-linear-double-torus-factor-split-one-place-obstruction.md)
really obtains two independent smooth-locus Cardano characters, but its two
factor partitions have two or three infinity places.  Thus character rank is
not intrinsically scarce; it is expensive in the one-place coordinate.
[THM-3944](../01-canon/theorems/THM-3944-repeated-factor-double-torus-one-place-square-conductor-collapse.md)
internally splits a repeated factor and really obtains a smooth one-place
parabola, but the quadratic order is nonnormal, one radicand is a cube, and
the remaining character does not extend across the conductor to normalized
`A2`.
[THM-3947](../01-canon/theorems/THM-3947-scalar-weighted-repeated-square-split-trichotomy.md)
shows that arbitrary reciprocal scalar weighting merely produces three
parabolas (or a doubled endpoint component), never one irreducible full
branch.  Finally,
[THM-3948](../01-canon/theorems/THM-3948-classified-weight-eight-nine-sextics-have-no-polynomial-normalization.md)
closes affine-line normalization throughout the classified irreducible
weight-eight/nine sextics.

The wildcard is therefore a five-coordinate conjunction rather than a vague
request for "more characters": a **normal** quadratic surface, a **reduced
one-place component**, control of the **full discriminant**, at least **two
extendable** cubic characters, and compatibility with the polynomial/Keller
side. The affine-linear split, repeated-square split, scalar-weighted split,
and classified high-weight sextics each lose a different coordinate. The next
honest generator should kill the first lost coordinate cheaply before
attempting a full atlas.

A final incoming pair now proves a large part of that handoff.
[THM-3946](../01-canon/theorems/THM-3946-affine-internal-factor-split-two-end-conductor-collision-dichotomy.md)
closes the complete affine one-factor grammar: its irreducible rows have two
or three infinity places, while the remaining rows are reducible.  Moreover,
[THM-3949](../01-canon/theorems/THM-3949-coprime-one-variable-internal-factor-splits-are-reducible-or-multi-ended.md)
is now **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED**.  Its
Newton-polygon and residual-cubic analysis proves that coprime nonassociate
one-variable internal factors make
the full discriminant reducible or at least two-ended on the standard line at
infinity.  The wildcard generator can therefore delete the entire coprime
one-variable row and concentrate on genuinely bivariate factors,
gcd/multiplicity overlap, or a different affine boundary line.

[THM-3950](../01-canon/theorems/THM-3950-a1-internal-split-denominator-debt-and-equianharmonic-shadow.md)
then supplies the sharp survivor. Its explicit degree-one-ratio example has
a reduced `A1` graph component, a normal quadratic surface and two independent
extendable Cardano classes, but its full discriminant also has an irreducible
genus-one companion. Universally, the factor-ratio map has degree three and
branch values `{0,1,-omega,infinity}`; its `S3` closure is the fixed `j=0`
elliptic cover. This connects the THM-3947 scalar wall to a residual-genus
invoice. It does not provide a same-field cubic atlas or a Keller map.
THM-3951, not THM-3950, is the RESERVED empty namespace stub.

The session therefore ends at a sharp boundary.  **PROVED:** degree-four
centered nonentry and all-degree trace-carrier routing.  **VERIFIED-EXACT:**
the 2,884-gate carrier companion, its period-twelve identity, and its ordinal
map, plus THM-3950's 51-gate equianharmonic packet. **OPEN:** color survival
for the seven degree-five tasks, the THM-3950 same-field cubic atlas and full
boundary/source incidence, arbitrary
root gauges and coefficient planes, higher-degree non-centered strata, and
the planar Jacobian conjecture itself.
