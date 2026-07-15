---
id: THM-803
title: The q=13 anti-grid ladder and the exact all-component erosion selector
status: PROVED (uniform half-divisor/parity-support and all-component selector theorems) + VERIFIED (complete 2-4-6 ladder, THM-797 correction, and sharp global-component limitation)
source: codex-2026-07-14-S10 q13-all-components analysis
depends_on:
  - THM-772
  - THM-774
  - THM-789
  - THM-797
related:
  - HYP-6820
verification:
  - 04-computation/lrc13_antigrid_all_component_selector_codex_S10.py
  - 05-knowledge/results/lrc13_antigrid_all_component_selector_codex_S10.out
---

# THM-803 — q=13 anti-grids and the all-component erosion selector

Use the two-sheet notation

```text
A=2U union {x,y},       |U|=10,       x,y odd,
E_U={t:phi_U(t)>=1/11},
R_U={d:max_(u in U)||ud||<2/143},
H_(x,y)={t:||at||+||bt||>=11/13},
a=(x+y)/2,              b=|x-y|/2.
```

THM-789's necessary erosion containment is

```text
E_U subset H_(x,y) minus R_U,                              (0)
```

where `H minus R={t:t+R subset H}`.  This theorem adds the
exception-divisor grids halfway between the integer grids of THM-797 and then
gives an exact finite selector for all of (0).

## 1. The half-divisor anti-grid

Let `q` be an odd divisor of `x`.  For a unit numerator modulo `2q`, taken
modulo sign, put

```text
D^(1/2)_q(U)={p mod +/- :
  gcd(p,2q)=1 and |[up]_(2q)|>=ceil(2q/11) for every u in U}.   (1)
```

### Theorem 1 — half-grid escape

Every class in `D^(1/2)_q(U)` gives an explicit failure of (0):

```text
t=p/(2q) lies in E_U minus H_(x,y).                         (2)
```

Consequently hypothetical two-sheet tightness forces

```text
D^(1/2)_q(U)=empty                                         (3)
```

for every odd divisor `q` of either exception.

### Proof

Write `x=qh`.  Since `x,q,h` are odd and `p` is a unit modulo `2q`, `hp` is
odd.  Therefore

```text
||xt||=||hp/2||=1/2>2/13.                                  (4)
```

Condition (1) is exactly `||ut||>=1/11` for every core speed, so `t in E_U`.
THM-774's folded diamond implies the individual exception eligibility
`||xt||<=2/13`; (4) therefore puts `t` outside `H_(x,y)`.  Finally
`0 in R_U`, hence `H minus R_U subset H`, proving (2)--(3).  Interchanging the
two exceptions proves the last assertion. ∎

## 2. The parity-twisted q=13 support gate

Now impose THM-772's conclusions: no member of `U` is divisible by 13, and a
deep two-sheet branch has a 13-divisible exception.  Let

```text
C=(Z/13Z)^*/{+/-1}={1,2,3,4,5,6}.
```

Besides THM-797's raw support

```text
S(U)={+/-u mod 13:u in U},                                 (5)
```

define the **parity-twisted support**

```text
P(U)={pi(u):u in U},
pi(u)=+/-u       mod 13,     u odd,
pi(u)=+/-(u/2)   mod 13,     u even.                        (6)
```

Here and below every residue in (5)--(6) is folded into `C`.

### Theorem 2 — exact half-grid support law

Identifying a unit numerator modulo 26 with its folded class in `C`, one has

```text
D^(1/2)_13(U)=C minus P(U)^(-1).                            (7)
```

Thus hypothetical containment forces the second full-support condition

```text
P(U)=C.                                                     (8)
```

For an even speed, the passage from its raw class to its twisted class is the
six-cycle

```text
sigma(c)=+/-7c mod 13=(1 6 3 5 4 2).                       (9)
```

Combining (8) with THM-797 gives

```text
S(U)=C or S(U)=C minus {+/-y}      if exactly x is 13-divisible,
S(U)=C                              if both exceptions are 13-divisible,
P(U)=C                              in either case.          (10)
```

In the aligned five-class case, if `c=+/-y`, there must be an even core speed
in raw class `+/-2c`.  This is the mandatory **parity bridge** into the missing
class.

### Proof

At a unit numerator `p mod 26`, an odd core speed has an odd nonzero balanced
residue.  Since `ceil(26/11)=3`, it fails to be `1/11`-deep exactly when

```text
up=+/-1 mod 26,
```

or equivalently when the folded class of `p` is the inverse of `pi(u)`.
For an even speed write `u=2v`.  Then

```text
||up/26||=||vp/13||,
```

and it fails the deep inequality exactly when `vp=+/-1 mod 13`, again the
inverse of `pi(u)`.  Taking the union of the ten forbidden numerator classes
proves (7), and Theorem 1 proves (8).

Because `2^(-1)=7 mod 13`, the even twist is (9).  The ordinary alternatives
in (10) are THM-797.  If raw class `c` is absent but twisted class `c` occurs,
its source cannot be odd.  It is therefore an even `u` with
`+/-(u/2)=c`, hence `+/-u=+/-2c`, proving the parity bridge. ∎

## 3. The complete universal even anti-shell ladder

Write the 13-divisible exception as `x=13X`, with `X` odd.  For even `d`, use
unit times

```text
t=p/(13d),       gcd(p,13d)=1.                             (11)
```

### Theorem 3 — exactly d=2,4,6 are universal

For each

```text
d in {2,4,6},                                                (12)
```

every unit time (11) makes `x` ineligible:

```text
||xt||=||Xp/d||>2/13.                                      (13)
```

Consequently (0) forces the ten core danger incidences to cover every signed
unit class modulo `13d`, or equivalently

```text
{p mod +/- : gcd(p,13d)=1 and p/(13d) in E_U}=empty.       (14)
```

The list (12) is complete among universal even denominators.  For every even
`d>=8`, the choice `X=p=1` has

```text
||Xp/d||=1/d<=1/8<2/13,                                   (15)
```

so exception ineligibility is no longer automatic.

### Proof

For `d=2`, the odd integer `Xp` has distance `1/2` from the nearest multiple
of two.  For `d=4`, its distance is `1/4`.  Modulo six an odd integer is
congruent to `+/-1` or `3`, giving distance `1/6` or `1/2`.  Each is strictly
larger than `2/13`; Theorem 1's argument proves (14).  Formula (15) proves
completeness. ∎

The `d=2` row is Theorem 2.  The `d=4,6` rows retain information destroyed by
reducing the time modulo 13, so they are genuine new gates rather than
restatements of raw folded support.

## 4. The exact all-component erosion selector

Let `closure(R_U)` denote the topological closure of the strict return set and
put

```text
K_U=E_U+closure(R_U) subset R/Z.                            (16)
```

Both terms are finite unions of closed circular intervals, where a degenerate
interval is allowed.  Decompose `K_U` into its closed circular components.  If
`K_U` is not the whole circle, define

```text
Sigma(U;x,y)=
  {all component endpoints of K_U}
  union K_U intersect {k/(2a):0<=k<2a}
  union K_U intersect {k/(2b):0<=k<2b}.                    (17)
```

If `K_U` is the whole circle, omit the endpoint term in (17).

### Theorem 4 — finite selector equivalence

The full erosion predicate is decided exactly by (17):

```text
E_U subset H_(x,y) minus R_U
iff K_U subset H_(x,y)
iff Q_(a,b)(t)>=11/13 for every t in Sigma(U;x,y),          (18)
```

where `Q_(a,b)(t)=||at||+||bt||`.

Let

```text
L=sum_(u in U)u,       W=max(x,y).
```

Then `E_U` and `closure(R_U)` each have at most `L` components, `K_U` is a
union of at most `L^2` circular arcs, and one may choose

```text
|Sigma(U;x,y)|<=2L^2+2W.                                  (19)
```

Under THM-772, with `B=max(U)`, this becomes the scale-normal quadratic bound

```text
|Sigma(U;x,y)|<=200B^2+22B.                               (20)
```

Every selector endpoint carries exact owner data.  An endpoint of `E_U` has
the form

```text
(11k+/-1)/(11u),                                          (21)
```

an endpoint of a return component has the form

```text
(143k+/-2)/(143v),                                        (22)
```

and every non-cusp endpoint of `K_U` is a sum of endpoints of the forms
(21)--(22).  Thus (17) is an owner-labelled component/return-incidence
selector, not a sample grid.

### Proof

The folded diamond `H_(x,y)` is closed.  For every `t`, therefore,

```text
t+R_U subset H_(x,y)
iff t+closure(R_U) subset H_(x,y).                         (23)
```

Taking the union over all `t in E_U` proves the first equivalence in (18).

The breakpoints of `||at||` are exactly `k/(2a)`, and those of `||bt||` are
exactly `k/(2b)`.  Between consecutive breakpoints their sum is affine.  On
each component of `K_U`, its minimum is therefore attained at a component
endpoint or at one of the breakpoints in (17).  This proves the second
equivalence.

The complement of `E_U` is the union, over `u in U`, of `u` open danger arcs;
the complement of `R_U` is the union of `u` closed nonreturn arcs.  Hence each
set has at most `L` components.  A sum of two connected circular arcs is a
connected arc or the whole circle, so the pairwise sums give at most `L^2`
arcs.  They contribute at most `2L^2` endpoints.  There are at most
`2a+2b=2W` folded cusps, proving (19).  THM-772 gives `L<=10B` and `W<=11B`,
which proves (20).  Finally every boundary of an interval intersection is a
boundary of an input tooth; this gives (21)--(22), and endpoints of a finite
union of pairwise sum arcs come from endpoints of a constituent sum. ∎

## 5. Correction to THM-797's sharp rows

THM-797's folded-support example is

```text
U_0=(1,2,3,4,7,9,10,11,12,16),       (x,y)=(13,5).         (24)
```

It is sharp for the **prime `q=13` support gate**, but not for the mandatory
anti-grid ladder.  At the quarter-grid point

```text
t=11/52,
phi_(U_0)(t)=5/52>1/11,
||13t||=1/4>2/13,
Q_(9,4)(t)=1/4<11/13.                                    (25)
```

The signed-wall-sharp row added to THM-797 concurrently,

```text
U_1=(1,2,4,6,7,9,10,11,12,16),       (x,y)=(13,5),
```

has exactly the same witness: `phi_(U_1)(11/52)=5/52` and
`Q_(9,4)(11/52)=1/4`.  Both examples remain sharp for their stated `q=13`
gates and retain the warning against selecting only global maximizers or
exception-divisor grids.  Neither survives the combined anti-grid ladder.

## 6. A sharp limitation after signed walls and anti-grids

Even the conjunction of THM-797's exact signed-wall theorem and the anti-grid
ladder does not replace the all-component selector.  Take

```text
U_*=(2,4,6,7,9,10,11,12,14,16),       (x,y)=(13,5).         (26)
```

This core is primitive, divisor-complete through 12, and has no 13-multiple.
Exact breakpoint evaluation gives

```text
M(U_*)=2/13,
argmax phi_(U_*)={5/13,8/13},
Q_(9,4)=12/13 at both maximizers.                           (27)
```

It satisfies THM-797's full signed-wall conclusion, not merely folded
alignment:

```text
{u mod 13:u in U_*}=(Z/13Z)^* minus {5,8},
S(U_*)=C minus {5}.                                        (28)
```

Its parity-twisted support is full; the required bridge into class `5` is
supplied by the even raw class `3`.  Every odd exception-divisor grid and every
universal anti-grid is silent or trapped.  Here `D_m` denotes the signed unit
numerators `p mod m` for which `p/m in E_(U_*)`:

```text
D_5(U_*)=empty,
D_13(U_*)={+/-5},        Q_(9,4)(5/13)=12/13,
D_26(U_*)=D_52(U_*)=D_78(U_*)=empty.                       (29)
```

The sharpened signed-wall metric pin and THM-772 determinant tax all pass:

```text
B=16,       rho=1/208,       x<=2B-1,       y<=B-1,
13B+2xy=338 <= 576=2B(x+y),
1/(13*5)+2rho=1/40 <= 36/845
                  =2/(13*13)+2/(13*5).                    (30)
```

Nevertheless the complete closed deep set has ten components:

```text
[1/22,5/88],               {7/22},
[23/66,27/77],             [67/176,43/110],
{9/22},                    {13/22},
[67/110,109/176],          [50/77,43/66],
{15/22},                   [83/88,21/22].                  (31)
```

The singleton `7/22` has endpoint owners `6` and `16`, and

```text
phi_(U_*)(7/22)=1/11,
Q_(9,4)(7/22)=9/22<11/13.                                 (32)
```

Thus the branch escapes only after retaining a nonmaximal closed component.
In this row

```text
closure(R_U*)=[-1/1144,1/1144],                            (33)
```

and the selector (17) has 24 exact entries.  Its global folded minimum is

```text
Q=463/1144 at t=365/1144.                                 (34)
```

The associated twelve-speed packet is loose but lies in the sporadic
max-peel inequality:

```text
2U_* union {13,5}=(4,5,8,12,13,14,18,20,22,24,28,32),
M=1/9,
M(after deleting 32)=1/8>1/12.                            (35)
```

Equations (26)--(35) are a method boundary, not a tight counterexample.  They
show that exception-divisor grids, THM-797's exact signed residue complement,
all three universal anti-grids, both support conditions, both global
maximizers, and all sharpened metric taxes do not select the required escape.
The owner-labelled all-component carrier in Theorem 4 is genuinely stronger.

## 7. Tournament Analysis and challenged assumptions

For (26), take the ten represented components of `K_U` as vertices.  The
pairwise observable is the signed erosion margin

```text
m(C)=11/13-min_(t in C)Q_(9,4)(t),                         (36)
```

computed by Theorem 4's selector.  The switch/gauge replaces it by component
width; left endpoint is the fixed tie Hamiltonian path.  The two orders are

```text
(3,6,4,5,0,9,2,7,1,8),
(1,4,5,8,2,7,3,6,0,9),                                   (37)
```

with 26 edge flips.  Both tournaments are transitive: score histogram
`(0,1,...,9)`, no directed triangles, ten singleton SCCs, and one Hamiltonian
path.  The first tournament becomes predicate-preserving only when every
vertex retains its exact interval, selector points, signed margin, endpoint
owners, and return incidence.  Its bare isomorphism class destroys the sign
of (36) and hence the containment verdict.

There is a second six-vertex telemetry tournament on folded residue
obligations.  Orient by raw multiplicity and switch to parity-twisted
multiplicity.  The tie paths are

```text
(5,1,2,3,4,6),       (1,3,4,2,5,6),                       (38)
```

with six edge flips; both are again transitive with score histogram
`(0,1,...,5)`, singleton SCCs, no cycles, and one Hamiltonian path.  The bare
tournament does not record which missing raw class is supplied by which even
source under the six-cycle (9).

The challenged assumption is that divisor grids, global maximizers, runners,
or bare component ranks can be the vertices of the final proof object.  The
predicate-preserving carrier is instead an owner-labelled component/return
incidence object whose finite obligations are exactly (17).

## Exact replay

The verifier exhausts the half-divisor identity through `q=101` (`71,673`
exact phase tests), all `2^12=4,096` signed support subsets modulo 26, and the
universal even anti-shell question through `d=40`.  It independently evaluates
(24)--(35), constructs every interval of `E_U`, `closure(R_U)`, and their
circle sum, checks all selector cusps and endpoints, and reproduces both
tournament fingerprints.  Its canonical digest is

```text
5e5851835ba5e753fd5776fe1c18c0639d778b037821aaa8b3e86f5358eb3495.
```

No floating-point or sampled-circle verdict enters the proof.
