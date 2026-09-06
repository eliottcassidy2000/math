---
id: THM-4411
title: "First-order collision-transgression equivalence and seminormal tradeoff"
status: >
  PROVED + VERIFIED-EXACT + TWO-PATH AUDITED. For nonzero branch tangents
  into a smooth surface, a normal deformation preserves a multiple-point
  collision to first order, up to branch reparametrization, iff every tangent
  relation kills its two-form transgression. In the THM-4381 exceptional
  quartic under Q(x)->Q(x)+s*h(x), this is exactly
  5h(-1)-18h(0)+13h(1)=0. Collision-preserving directions cannot fill the
  THM-4404 seminormal cokernel. The lawful direction h=x fills all 89
  quotient dimensions, but precisely because it splits the retained triple
  to first order. This is not the full moving-graph complex and proves no
  chart/seam entry, Keller pair, JC(2), or DC(2) claim.
source: root + variable_normal_transgression / JC2 and arXiv continuation session, 2026-09-03--05
audit: >
  PASS. The exact symbolic portion verifies the dual-number compiler for
  h=1,x,x^3-x, the surface relation, tangent and normal rows, wedge formula,
  unique incidence relation, collision solutions, and the coordinate-free
  equivalence. The lower bounds use both independently audited THM-4404
  implementations: carrier-target pairs at (p,alpha)=(421,126) and T-module
  Kahler generators at (443,112). They independently give ranks 88,89,88 for
  h=1,x,x^3-x. Finite-field minors are used only for sound characteristic-
  zero lower bounds; exact collision identities and dim E_x=89 give the
  upper bounds. Normal, optimized, and fixed-hash-seed replays byte-match the
  frozen output and perform 813 dynamic checks.
depends_on:
  - THM-4381-exceptional-quartic-seminormalization-and-conductor-fibre-classification
  - THM-4404-exceptional-quartic-descended-two-form-seminormal-cokernel
related:
  - THM-4408-rank-two-poisson-hamiltonian-primitive-nondescent-compensator-firewall
  - THM-4067-seminormal-period-kernel-and-figure-eight-completeness-obstruction
  - THM-4397-rank-two-poisson-counterexample-symplectic-gauge-equivalence
script: 04-computation/jc2_exceptional_quartic_variable_normal_transgression_scout_s616.py
output: 05-knowledge/results/jc2_exceptional_quartic_variable_normal_transgression_scout_s616.out
script_sha256: 23e46d84579706a018ba71382d3157c079488bf21d6d29ea42875c96e5de403e
output_sha256: dae47be81d6a4f29a32898bfbb572fe71c42f8516ffb7113efda6106487040c6
hash_basis: raw LF bytes
---

# THM-4411 -- collision persistence is exactly vanishing transgression period

**PROVED + VERIFIED-EXACT + TWO-PATH AUDITED, AT FIRST ORDER AND WITHIN THE
EXCEPTIONAL-QUARTIC VARIABLE FIXED-`x` NORMAL FAMILY SPECIFIED BELOW. THIS IS
NOT AN IDENTIFICATION WITH THE FULL MOVING-GRAPH COMPLEX, AN ADMISSIBLE
SOURCE-NORMAL `JC(2)` DIRECTION, CHART OR SEAM ENTRY, A KELLER PAIR, `JC(2)`,
OR `DC(2)`.**

**PROVED FORMAL-LOCAL EXTENSION — 2026-09-05.** The
[labelled-triple compensation theorem](../../05-knowledge/results/synthesis_20260905_transgression.md)
integrates the local coefficient constraint: one additional `x` coefficient
makes the endpoint/common-target Jacobian invertible with determinant `-288`,
giving unique formal compensation at every order. The compatible tangent
`h=4x^2-9x` nevertheless needs a nonzero second-order correction; its straight
pencil does not preserve the triple modulo `s^3`. This does not supply
polynomial termination, a descended primitive, or the full graph complex.

## 1. A coordinate-free collision--transgression lemma

Let `k` be a field, let `V` be a two-dimensional `k`-vector space, and fix a
nonzero area form

```text
omega in wedge^2 V^*.
```

For finitely many labelled branches `i in I`, let

```text
t_i in V\{0}        be the ordinary tangent,
n_i in V            be the normal velocity.            (1)
```

Say the collision persists to first order, allowing reparametrization of
each branch, if there are scalars `v_i` and one common target velocity `w`
such that

```text
n_i+v_i t_i=w                 for every i.              (2)
```

Define the branch transgression and target-motion map

```text
tau_i=omega(t_i,n_i),
L:V -> k^I,          L(w)=(omega(t_i,w))_i.             (3)
```

Then

```text
boxed: (2) holds iff tau belongs to image(L).           (4)
```

Indeed, `(2)` gives `tau_i=omega(t_i,w)`. Conversely, if
`tau_i=omega(t_i,w)`, then

```text
omega(t_i,n_i-w)=0.
```

Because `V` is two-dimensional, `omega` is nondegenerate and `t_i` is
nonzero, the kernel of `z |-> omega(t_i,z)` is exactly `k t_i`. Thus
`n_i-w` is a scalar multiple of `t_i`, giving `(2)`.

The left kernel of `L` is equally concrete. A row `ell=(ell_i)` annihilates
`image(L)` iff

```text
sum_i ell_i t_i=0.                                     (5)
```

Consequently `(4)` is equivalent to the period criterion

```text
boxed: sum_i ell_i tau_i=0 for every tangent relation
       sum_i ell_i t_i=0.                              (6)
```

This statement is invariant under changing surface coordinates or rescaling
the area form. It is purely first-order: it neither integrates the motion nor
asserts that it comes from a global polynomial graph.

## 2. Variable normal directions in the exceptional compiler

Use the THM-4404 exceptional-quartic compiler over its characteristic-zero
field `K`:

```text
N=K[x],
D=1+x^2 Q(x),
B=(D-1)(D+2)^2,
C=xD(D+2),
E=Q(D+3),
S=K[B,C,E] subset N.                                   (7)
```

For an arbitrary polynomial `h in N`, make the lawful dual-number
deformation

```text
Q(x) |-> Q(x)+s h(x),                  s^2=0.           (8)
```

The identity `C^2E=B(B+4)` remains true over `K[s]/(s^2)`. Differentiating
`(7)` gives

```text
gamma_(1,h)=h gamma_(1,1),
gamma_(1,1)(B)=3x^2D(D+2),
gamma_(1,1)(C)=2x^3(D+1),
gamma_(1,1)(E)=2(D+1).                                 (9)
```

At the retained normalization points `x=(-1,0,1)`, which all map to
`(B,C,E)=(0,0,-3)`, the ordinary rows and unit-normal rows are

```text
B'              =( 0,0,0),       gamma_(1,1)(B)=( 0,0, 0),
C'              =( 3,3,3),       gamma_(1,1)(C)=( 2,0,-2),
E'              =(-9,4,9),       gamma_(1,1)(E)=(-2,4,-2). (10)
```

Writing `h_i=h(i)`, contraction of `dC wedge dE` with the ordered
`(d/dx,d/ds)` tangent pair is

```text
W_h=(C' gamma_(1,h)(E)-E' gamma_(1,h)(C))_i
   =12(h_-1,h_0,h_1).                                  (11)
```

The unique tangent relation in `(10)` is

```text
ell=(5,-18,13),
5t_-1-18t_0+13t_1=0.                                  (12)
```

Thus the abstract criterion `(6)` becomes

```text
boxed: the retained triple persists to first order
       iff 5h(-1)-18h(0)+13h(1)=0.                    (13)
```

For completeness, let `(c,e)` be the common first variation of `(C,E)` and
let `v_i` move the three normalization parameters. The literal equations are

```text
3v_i+delta C_i=c,
E'_i v_i+delta E_i=e,
```

or, after eliminating `v_i`,

```text
e-E'_i c/3=4h_i.                                       (14)
```

Their compatibility is exactly `(13)`. When it holds, one solution is

```text
c=2(h_-1-h_1)/3,
e=2(h_-1+h_1),
v_i=(c-delta C_i)/3.                                   (15)
```

## 3. The seminormal period is the collision obstruction

Recall THM-4404's vector-space quotient

```text
E_x=N/dS,                    dim_K E_x=89,              (16)
```

and its normalized retained functional

```text
Lambda(v)=5v(-1)/18-v(0)+13v(1)/18.                   (17)
```

It is well defined on `E_x`; its cokernel line is generated by the derivative
of THM-4381's seminormal class `r`. For arbitrary target functions `P,Q`, the
normal in `(9)` gives

```text
tau_h(dP wedge dQ)=h tau_1(dP wedge dQ).               (18)
```

Equation `(18)` is an equality of raw densities in `N` before passage to
`E_x`; it does not assert that multiplication by `h` is a well-defined
endomorphism of the additive quotient `N/dS`.

At the retained triple, the chain rule evaluates both target gradients at
the same point. Therefore every descended target two-form has retained
transgression equal to a scalar multiple of `W_h`. Equations `(11)--(13)`
give the cutoff-free implication

```text
collision persists  =>  image(tau_h) subset ker Lambda. (19)
```

Since `ker Lambda` has dimension 88, no collision-preserving direction in
the entire polynomial family `(8)` can pay the missing seminormal line.
This is the structural content: the rank loss in THM-4404 was not an accident
of choosing the constant normal `h=1`.

The converse at the retained first-jet level is also exact. A nonzero
seminormal period is precisely the obstruction to solving the common-motion
equations `(14)`. It may enlarge the two-form image only by leaving the
collision-preserving locus.

## 4. Three sharp controls and exact ranks

The theorem companion tests three directions:

```text
h=1:       values (1,1,1),       ell(h)=0,
h=x:       values (-1,0,1),      ell(h)=8,
h=x^3-x:   values (0,0,0),       ell(h)=0.             (20)
```

For `h=1`, equation `(15)` gives

```text
(c,e)=(0,4),             (v_-1,v_0,v_1)=(-2/3,0,2/3), (21)
```

so the triple moves together; THM-4404 proves `rank tau_1=88`. The direction
`h=x^3-x` is a hostile control: it is nonzero globally but invisible on the
retained fibre. Equation `(19)` gives rank at most 88, and both positive-minor
calculations below give rank at least 88.

For `h=x`, the period in `(20)` is nonzero, so no endpoint velocities and
common target first jet solve `(14)`. Nevertheless its descended two-form
image has full rank:

```text
boxed: rank tau_1=88,
       rank tau_x=89,
       rank tau_(x^3-x)=88.                            (22)
```

The two independent lower-bound paths are:

```text
(p,alpha)=(421,126): carrier/target-pair columns,
(p,alpha)=(443,112): T-module multiples of three Kahler wedges. (23)
```

At each good reduction they obtain the rank sequence `(88,89,88)`. A
nonzero reduced minor proves the corresponding characteristic-zero minor is
nonzero. The exact bound `(19)` supplies both rank-88 upper bounds, while
`dim E_x=89` supplies the rank-89 upper bound. Thus `(22)` is exact; no rank
failure in finite characteristic is used as evidence.

The `h=x` direction is therefore the decisive hostile to a tempting but
false inference:

```text
"allowing variable normals repairs the seminormal miss while retaining the
triple"                                                    -- false. (24)
```

It repairs the miss **by destroying the triple at first order**.

## 5. Relation to Long's rank-two Poisson construction

Long's construction, audited in THM-4397/4408, pays a residual two-form by a
non-descended Hamiltonian primitive and an independent conjugate coordinate.
The lemma above explains why importing only the wedge correction into a
planar collision problem fails: in a surface, the same wedge period is the
complete first-order obstruction to moving all branches together.

The extra coordinate in a symplectic suspension can carry branch-sensitive
data without asking that it descend through the planar core. A planar map has
no such free compensator. This is an exact mechanism-level relation, not a
reduction from the Poisson construction to `JC(2)` and not a proof that every
possible planar deformation is captured by `(8)`.

## 6. Scope and reproduction

The family `(8)` fixes the normalization coordinate `x` and varies only the
compiler polynomial `Q`. It is lawful for the surface relation, but an
arbitrary `h` need not correspond to an admissible coefficient direction in
the full source-normal Jacobian chart. The theorem is first-order and says
nothing about integration, higher-order multiple-point persistence, the full
moving-endpoint graph complex of THM-4067, or chart/seam entry.

Replay from the repository root:

```text
python3 -B 04-computation/jc2_exceptional_quartic_variable_normal_transgression_scout_s616.py
python3 -B -O 04-computation/jc2_exceptional_quartic_variable_normal_transgression_scout_s616.py
PYTHONHASHSEED=2718 python3 -B 04-computation/jc2_exceptional_quartic_variable_normal_transgression_scout_s616.py
```

All streams byte-match the frozen LF output. The script performs `813`
dynamic exact and good-reduction checks; its hashes are in the frontmatter.
