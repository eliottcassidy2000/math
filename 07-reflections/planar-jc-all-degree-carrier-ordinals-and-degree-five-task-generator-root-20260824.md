# The carrier grammar is now an ordinal task generator, not a conjectural closure

**Status (2026-08-24): SYNTHESIS OF PROVED + VERIFIED-EXACT RESULTS; JC(2)
REMAINS OPEN.**  The proof sources are
[THM-3920](../01-canon/theorems/THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction.md),
[THM-3929](../01-canon/theorems/THM-3929-root-regular-one-place-linear-color-cubic-is-monogenic.md),
[THM-3933](../01-canon/theorems/THM-3933-centered-degree-three-root-map-pole-partition-octic-nonentry.md),
[THM-3936](../01-canon/theorems/THM-3936-centered-degree-three-infinite-root-value-nonentry.md),
[THM-3938](../01-canon/theorems/THM-3938-centered-degree-four-root-map-nonentry.md),
[THM-3941](../01-canon/theorems/THM-3941-all-degree-centered-cubic-pole-carrier-routing.md),
[THM-3959](../01-canon/theorems/THM-3959-centered-degree-five-root-map-color-closure.md),
and
[THM-3961](../01-canon/theorems/THM-3961-arbitrary-q-hidden-repetition-normality-and-conductor-debt.md).
This reflection records the research move and its degree-five resolution; it
is not an additional theorem.

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

## 5. The seven degree-five tasks are now closed

Degree five has the following zero-based work queue:

```text
0: m=0, C3(5)          1: m=0, C2(5)
2: m=1, C3(4)          3: m=1, C2xC2(1,3)
4: m=2, C2(3)
5: m=4, C3(1)          6: m=4, C2(1).
```

The separate `m=5` polynomial row is the root-regular THM-3929 exit.  For
each task `0,...,6`, the generated obligation was uniform:

1. solve the lower-principal-part trace equations;
2. divide the exact color and retain only polynomial families;
3. test the coefficient ideal for the unit ideal;
4. distinguish maximal-order ramification from square index factors; and
5. compute the normalization-address fibres of every genuine component.

[THM-3959](../01-canon/theorems/THM-3959-centered-degree-five-root-map-color-closure.md)
now executes this queue completely. Exact color division makes L/M/O empty.
The five surviving families are J1/J2/K/N/P. J1 has one scalar endpoint in an
actual normal maximal order; every nonscalar survivor has a genuine reduced
ramification arm with two or three normalization addresses over one target
point. Rank-three finite flatness coalesces them at one source point, so the
arm is non-unibranch and THM-3920 excludes it. The 105-gate companion audits
all resultants, seams, saturations, address fibres and the scalar endpoint.

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

An incoming chain now proves a large part of that handoff.
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
invoice.

[THM-3951](../01-canon/theorems/THM-3951-affine-plane-boundary-incidence-forest-and-equianharmonic-survivor-nonentry.md)
shows more generally that every affine-plane boundary prime is rational,
whereas the residual prime has positive genus for arbitrary nonzero common
debt. Thus every natural nonconstant-ratio same-field cubic is excluded. In a
clean row, two distinct smooth graph/residual meetings also create a cycle in
the common-resolution boundary tree. Then
[THM-3952](../01-canon/theorems/THM-3952-minimal-mobius-internal-split-carriers-are-four-critical-colors-and-nonentry.md)
classifies unit-debt degree-one ratios. Only four critical infinity colors
give polynomial `A1` carriers, and every one retains the forbidden repeated
incidence. Corrected/audited THM-3953 closes all three-distinct polynomial
roots, including constant ratios; THM-3954 proves the local common-debt
`A_(3m-1)`/non-unibranch refinement. THM-3956 closes every split hidden cubic,
THM-3958 closes the exactly-one-root case, independently audited THM-3960
closes the natural globally monogenic one-parameter family, and THM-3959
closes all seven centered degree-five carrier rows. THM-3961 then extends the
monogenic normality gate to irreducible arbitrary `q(P,t)`: adjusted hidden
squarefreeness is equivalent to normality, and every normal row is excluded
by the global different. Its repeated-factor and `P^2` nonnormal debts remain;
THM-3962 closes all coefficient-constant cylinders, including both debts,
while THM-3963 closes the moving scalar family `q=c(t)P^2`. General moving
repeated factors and `P^2q2(P,t)` remain. THM-3964 reserves a graph-double-root
subfamily; THM-3965 reserves a separate unit-ideal deformation lane.

The session therefore ends at a sharp boundary.  **PROVED:** degree-four
centered nonentry, all-degree trace-carrier routing, the arbitrary-common-debt
natural-cubic boundary obstruction, the four-color Mobius classification, and
the centered degree-five and arbitrary-`q` normal monogenic closures.
**VERIFIED-EXACT:** the 2,884-gate carrier companion and the
`51/62/44/75/92/62/161/105/42/51/46/45` gates in
THM-3950/51/52/53/54/56/58/59/60/61/62/63. **OPEN:** degree at least six,
arbitrary root gauges and coefficient planes, the remaining moving conductor
debts, nonmonogenic orders, normalization-parameter descent, higher-degree
non-centered strata, and the planar Jacobian conjecture itself.
