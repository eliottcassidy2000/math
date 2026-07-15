---
id: THM-807
title: Connected-return linear erosion selector and actual signed-complement grid-ladder limits
status: PROVED (uniform central/connected-return selector) + VERIFIED (exact signed-complement ladder limitations)
source: codex-2026-07-14-S10 post-THM-803 selector analysis
depends_on:
  - THM-772
  - THM-774
  - THM-789
  - THM-797
  - THM-803
related: [THM-817, HYP-6820]
verification:
  - 04-computation/lrc13_connected_return_linear_selector_codex_S10.py
  - 05-knowledge/results/lrc13_connected_return_linear_selector_codex_S10.out
---

# THM-807 — connected-return linear selector and grid-ladder limits

Use the two-sheet notation

```text
A=2U union {x,y},          |U|=10,
B=max(U),                  L=sum_(u in U)u,
gamma=1/11-1/13=2/143,
E_U={t:phi_U(t)>=1/11},
R_U={d:max_(u in U)||ud||<gamma},
a=(x+y)/2,                 b=|x-y|/2,
Q_(a,b)(t)=||at||+||bt||,
H_(x,y)={t:Q_(a,b)(t)>=11/13}.
```

THM-803 decides the full erosion predicate by an exact selector of size at
most `200B^2+22B`.  The first result below isolates a large branch on which
the quadratic term disappears.

## 1. The mandatory central return component

Put

```text
delta=gamma/B=2/(143B),       J_B=[-delta,delta] in R/Z,
C_U=E_U+J_B.
```

### Theorem 1 — a uniform linear necessary selector

One always has

```text
J_B subset closure(R_U),       C_U subset E_U+closure(R_U).       (1)
```

Decompose `C_U` into closed circular components.  If it is not the whole
circle, let

```text
Lambda(U;x,y)=
  {all component endpoints of C_U}
  union C_U intersect {k/(2a):0<=k<2a}
  union C_U intersect {k/(2b):0<=k<2b}.                           (2)
```

If `C_U` is the whole circle, omit the endpoint term.  Then

```text
C_U subset H_(x,y)
iff Q_(a,b)(t)>=11/13 for every t in Lambda(U;x,y),                (3)
```

and one may choose

```text
g=gcd(a,b)=gcd(x,y),
|Lambda(U;x,y)| <= 2L+2max(x,y)-2g.                               (4)
```

More precisely, if `C_U` is not the whole circle and has `s` components,
then the deduplicated selector has

```text
|Lambda|=2s+|(C_U intersect Cusp(a,b)) minus boundary(C_U)|,
Cusp(a,b)={k/(2a)} union {k/(2b)}.
```

If `C_U` is the whole circle, `|Lambda|=|Cusp(a,b)|`.

In a hypothetical tight two-sheet packet, THM-772 therefore gives the
**linear necessary selector**

```text
|Lambda(U;x,y)| <= 42B-2g <= 42B-2.                               (5)
```

Every non-cusp endpoint in (2) retains an exact owner label.  It comes from

```text
(11k+/-1)/(11u) +/- 2/(143B),       u in U,                        (6)
```

modulo one.  Thus (2) is an owner-labelled component/thickness selector, not
a sample mesh.

### Proof

For the representative `|d|<=delta` and every `u<=B`,

```text
||ud||=u|d|<=B delta=gamma.
```

The endpoints are approached from the interior toward zero, so the whole
closed interval belongs to `closure(R_U)`.  This proves (1).

The set `E_U` has at most `L` circular components: its complement is the
union of the `L` open `1/11`-danger teeth contributed by the ten speeds.
Adding one connected circular arc sends every component of `E_U` to one
connected arc (or the whole circle), and mergers can only reduce the number
of components.  Hence `C_U` has at most `L` components and at most `2L`
endpoints.

The two cusp grids in (2) have `2a` and `2b` points, with intersection of
size `gcd(2a,2b)=2g`.  Their union therefore has exactly
`2a+2b-2g=2max(x,y)-2g` points.  Between consecutive cusps `Q` is affine.
Its minimum on each component of `C_U` is consequently attained at an
endpoint or cusp, proving (3)--(4); endpoint/cusp coincidences and component
mergers can only improve the bound.  Since `x,y` are odd,
`gcd(a,b)=gcd(x,y)`.  THM-772 gives `L<=10B` and
`max(x,y)<=11B`, proving (5).  Boundaries of
`E_U` have the first form in (6), and thickening by `J_B` supplies the final
shift.  A boundary left after merging constituent arcs is still a boundary
of one constituent, proving the owner assertion. ∎

## 2. Connected returns collapse the full selector to linear size

### Theorem 2 — exact connected-return selector

The connected component of `closure(R_U)` containing zero is exactly `J_B`.
Consequently, if `closure(R_U)` is connected, then

```text
closure(R_U)=J_B                                                (7)
```

and the full erosion predicate is equivalent to the linear selector:

```text
E_U subset H_(x,y) minus R_U
iff E_U+closure(R_U) subset H_(x,y)
iff Q_(a,b)(t)>=11/13 for every t in Lambda(U;x,y).              (8)
```

Thus THM-803's general quadratic selector improves from
`200B^2+22B` to `42B-2gcd(x,y)` on every connected-return core.

### Proof

The interval `J_B` lies in the zero component by (1).  Immediately beyond
either endpoint, while `|d|<1/(2B)`, the speed `B` has nearest integer zero
and satisfies `||Bd||>gamma`.  Hence no larger interval about zero lies in
the closed return set; its zero component is exactly `J_B`.  If the full
closed return set is connected, there are no other components, proving (7).
THM-803's closure equivalence and Theorem 1 now give (8). ∎

The unconditional content is deliberately one-sided: (3) is always a
linear exact test for the mandatory central thickening, so its failure proves
the full erosion failure.  Passing (3) does not settle a core with noncentral
return components.  In the connected-return branch it is predicate-exact.

## 3. Actual signed-complement grids do not extend to a short fixed ladder

The next rows satisfy the **actual** residual arithmetic, not merely an
abstract residue relaxation.  In both, take

```text
(x,y)=(13,5),              (a,b)=(9,4).                         (9)
```

Each core is primitive, contains a multiple of every `m=2,...,12`, contains
no multiple of 13, has

```text
{u mod 13:u in U}=(Z/13Z)^* minus {5,8},
S(U)=C minus {5},          P(U)=C,                               (10)
```

and passes THM-772's `(A*)` tax and THM-797's sharpened metric pins.  Thus
these examples retain the exact signed complement and the parity bridge.

For a positive multiplier `d`, define the complete unit-grid profile

```text
P_d(U)={(p,phi_U(p/(13d)),Q_(9,4)(p/(13d))):
        1<=p<=13d/2, gcd(p,13d)=1, phi_U(p/(13d))>=1/11}.         (11)
```

### Theorem 3A — every multiplier through seven can be silent

Let

```text
U_all=(45,48,50,54,55,62,85,95,105,116).                        (12)
```

Then

```text
P_1(U_all)={(5,2/13,12/13)},
P_d(U_all)=empty                    for d=2,...,7.                (13)
```

Thus every `13d` unit grid through `d=7` is silent or trapped.  The very next
multiplier detects escape:

```text
P_8(U_all)=
 {(31,5/52,53/104),(45,11/104,3/8)}.                             (14)
```

Both displayed `Q` values are below `11/13`.

This core also passes the scalar taxes exactly:

```text
B=116,       M(U_all)=45/161,
rho=106/60697,
1/(13*5)+2rho=5729/303485 <= 36/845.                             (15)
```

Its return set is connected,

```text
closure(R_Uall)=[-1/8294,1/8294].                                (16)
```

The exact linear selector has 106 components and 212 points.  Its global
minimum is

```text
Q=709/28710 at t=709/373230,                                    (17)
```

so the connected-return selector sees a large escape which the initial grid
ladder omits.

### Theorem 3B — the even anti-grid ladder can be silent through eighteen

Let

```text
U_even=(6,9,20,24,30,36,42,54,66,90).                            (18)
```

For every even `d<=18`, the complete profile is empty except

```text
P_10(U_even)={(51,7/65,9/10)}.                                  (19)
```

The one surviving class is trapped because `9/10>11/13`.  Hence not only the
universal `d=2,4,6` anti-shells, but every full even multiplier grid through
18, can simultaneously fail to select escape under (10) and all the scalar
taxes.  Multiplier 20 finally gives

```text
P_20(U_even)={(73,7/65,31/52)},                                 (20)
```

which escapes.  An odd grid detects the row much earlier:

```text
p/(13d)=11/39,       phi_U=2/13,       Q_(9,4)=23/39.             (21)
```

Thus (18) is a boundary for the **even anti-grid ladder**, not evidence that
odd or component-generated denominators are dispensable.

The exact scalar and return data are

```text
B=90,        M(U_even)=2/13,
rho=1/1170,
1/(13*5)+2rho=2/117 <=36/845,
closure(R_Ueven)=[-1/6435,1/6435].                               (22)
```

The linear selector has 48 components and 96 points, with

```text
min Q=1159/5445 at t=1159/70785.                                (23)
```

### Verification and scope

The exact replay enumerates every unit numerator in (13)--(14) and
(19)--(21), constructs every closed component of `E_U` and `closure(R_U)`,
forms their circular sum, and evaluates every endpoint and folded cusp with
`Fraction` arithmetic.  It independently enumerates all pair crossings and
self-cusps to certify the two values of `M(U)`.  The canonical digest is

```text
39ad9ab6c2e7d932103d798b357b1f29beaa960bb4d4dd5b44a9cec895d04a7b.
```

The rows are loose method boundaries, not tight twelve-speed
counterexamples.  Their force is exact: no theorem using only the listed
necessary hypotheses can assert escape on one of `d=1,...,7`, or on one of
the even multipliers through 18.

## 4. Tournament Analysis and challenged assumptions

For (12), take multipliers `1,...,8` as vertices.  The pairwise observable is
the signed grid-detection margin; switch the gauge to the number of unit
classes killed by core danger, with increasing multiplier as the fixed tie
Hamiltonian path.  The two orders are

```text
(2,3,4,5,6,7,1,8),       (1,2,3,4,6,8,5,7),                     (24)
```

with nine edge flips.  Both tournaments have score histogram `0,...,7`, no
directed triangles, eight singleton SCCs, and one Hamiltonian path.

For (18), use even multipliers `2,...,20`.  The orders are

```text
(2,4,6,8,12,14,16,18,10,20),
(2,4,6,10,8,12,14,18,20,16),                                  (25)
```

with seven edge flips.  Both have score histogram `0,...,9`, no directed
triangles, ten singleton SCCs, and one Hamiltonian path.

These tournaments are telemetry only.  They preserve a ranking of finite
grid obligations but destroy the sign and location of non-grid component
escape, endpoint owners, and return incidence.  Runner vertices lose the
denominator obligations; denominator vertices lose component geometry;
folded residue vertices lose signed occupancy and parity twist; cusp vertices
alone lose which deep component reaches them.  Fixed circle sections and gaps
likewise forget the active owner unless decorated.

The assumption challenged here is that a sufficiently long *fixed* grid
ladder can replace THM-803's component carrier.  In the connected-return
branch, the predicate-preserving quotient is instead the owner-labelled
component/cusp incidence object in (2): it remembers exactly whether
`E_U+closure(R_U)` lies in the folded diamond.  Outside that branch one must
also retain the noncentral return-component incidence, and the general
quadratic selector remains the honest carrier.

**Forward sharpening.** THM-817 classifies those noncentral components as
signed maximum-speed cells and gives the exact adaptive carrier
`2c_E N_R+2W-2g`.  It also constructs a primitive signed-complement family
with `N_R=Theta(B)`, so connectedness cannot be inferred from the present
arithmetic/scalar gates.  This supersedes “quadratic” as a claim about
necessary complexity: satellites provide the second factor, but no quadratic
lower bound for a minimal selector is asserted.
