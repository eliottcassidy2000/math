---
id: THM-3590
title: "A13 exponent-three three-by-three weight Darboux nonentry"
status: >
  PROVISIONAL + VERIFIED-EXACT / PENDING INDEPENDENT AUDIT.  Using the proved
  per-coordinate two-central-arm weight invoice of THM-3589, no polynomial
  Darboux pair on the A13
  exponent-three target has exactly three nonconstant weight pieces in each
  coordinate.  The complement-edge cases h=0,1,2,3 are complete.  This file
  remains outside the proved dependency graph until the present proof is
  independently audited and explicitly promoted.
source: root / delegated A13 three-by-three Darboux hostile, 2026-08-21
audit: >
  The arm-unit classification, complement matching, all three two-edge
  collision trees, unequal-gap six-cycle, and equal-gap congruence boundary
  have exact controls.  Normal and optimized runs are byte-identical to the
  stored 92,630-gate output; independent mathematical audit remains pending.
depends_on:
  - THM-3581-critical-value-multiarm-keller-compiler-and-A13-carrier
  - THM-3589-danielewski-central-arm-every-line-and-kummer-trace-darboux-gates
related:
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
  - THM-3579-equal-step-three-by-three-danielewski-darboux-nonentry
  - THM-3584-all-exponent-equal-step-three-by-three-danielewski-darboux-nonentry
external:
  - "Gwozdziewicz, Injectivity on one line, arXiv:alg-geom/9305008, Theorem 1.1."
script: 04-computation/jc2_a13_exponent_three_three_by_three_nonentry_thm3590.py
output: 05-knowledge/results/jc2_a13_exponent_three_three_by_three_nonentry_thm3590.out
script_sha256: 5d2864467984bcc9bd7929174f5c0a3b746b9238b14c528f9dfa4ac678751762
output_sha256: 1a7cae06671ec33499bd564735796791d6e7cd38429b6dde4735c8183e5a8702
hash_basis: LF-normalized bytes
---

# THM-3590 -- A13 exponent-three three-by-three weight Darboux nonentry

**PROVISIONAL + VERIFIED-EXACT / PENDING INDEPENDENT AUDIT.**  The explicit
two-central-arm invoice `(9)` is the proved THM-3589 input.  The new content
here is the complete complement-edge and singleton-graph proof.

This closes the complete `3 x 3` weight cell for the degree-thirteen carrier.
It constructs neither a Darboux pair nor a planar Jacobian counterexample.
The same arm/jet invoice closes every cell with a two-piece coordinate.  The
`3 x 4` and larger support cells remain open.

## 1. Target, grading, and complement graph

All rings are over `C`.  Fix `kappa!=0`, put

```text
Sigma(b)=b(b^2+kappa^2),
A=C[b,c,e]/(c^3 e-Sigma(b)),                            (1)
```

and use the symplectic bracket and grading

```text
{b,c}=c^3,       {c,e}=-Sigma'(b),       {b,e}=-3c^2e,
wt(b)=0,         wt(c)=1,                wt(e)=-3.     (2)
```

On `c!=0`, a homogeneous piece of weight `u` is uniquely `c^u f(b)`.
It is regular in `A` exactly when

```text
Sigma^ceil(-u/3) divides f                    if u<0.  (3)
```

For two homogeneous pieces define

```text
W_(u,v)(f,g)=v f'g-u f g'.                            (4)
```

Then

```text
{c^u f,c^v g}=c^(u+v+2) W_(u,v)(f,g).                (5)
```

Suppose, after subtracting constants, that

```text
P=sum_(i=0)^2 c^r_i f_i(b),      r_0<r_1<r_2,
Q=sum_(j=0)^2 c^s_j g_j(b),      s_0<s_1<s_2,         (6)
```

with all six displayed pieces nonzero.  A weight-zero coefficient retained
in `(6)` is nonconstant.  Assume for contradiction that `{P,Q}=1`.

The scalar coefficient is the sum over edges

```text
r_i+s_j=-2.                                           (7)
```

Because both supports are strictly ordered, these edges form an order-
reversing matching.  Let `h` be their number, so `0<=h<=3`.

## 2. Two-central-arm invoice from THM-3589

The proved THM-3589 per-coordinate statement says that, for each
`R in {P,Q}`, on the two central arms

```text
L_c={b=e=0},                  L_e={b=c=0},             (8)
```

one has

```text
deg(R|L_c)>=2,                 deg(R|L_e)>=2.          (9)
```

Here is the exact weight content of `(9)`.  On `L_c`, only nonnegative
pieces survive, and `c^u f(b)` restricts to `c^u f(0)`.  Thus every output
has a supported weight

```text
p>=2                  with the arm coefficient nonzero. (10)
```

For `u=-s<0`, write `m=ceil(s/3)` and `f=Sigma^m h`.  Then

```text
c^(-s)f=c^(3m-s)e^m h.                               (11)
```

This survives on `L_e` exactly when `s=3m`.  Its degree there is `m`.
Consequently every output also has a supported weight

```text
-3m,                    m>=2,                         (12)
```

whose coefficient after division by `Sigma^m` is nonzero at `b=0`.

There is a third, local requirement at the central intersection

```text
p_0=(b,c,e)=(0,0,0).                                  (12a)
```

Since `Sigma'(0)=kappa^2!=0`, the tangent coordinates at `p_0` are `c,e`.
At that point the Darboux equation reads

```text
-kappa^2(F_cG_e-F_eG_c)=1.                            (12b)
```

Thus the two linear jets are nonzero and independent.  A linear `c` term
has weight `1`, while a linear `e` term has weight `-3`.  Combined with
the distinct nonlinear arm pieces `(10),(12)`, every Darboux coordinate
has at least three weight pieces: a low exact `-3m` with `m>=2`, at least
one middle linear weight in `{-3,1}`, and a high weight at least `2`.
In a `3 x 3` cell each coordinate has exactly one of the two middle weights,
and invertibility of `(12b)` forces the two choices to be opposite.  After
reflection they are therefore canonically

```text
r_1=-3,                         s_1=1.                (12c)
```

This also gives the target-wide corollary that no Darboux pair on the A13
target can have only two weight pieces in either coordinate.

The proof mechanism of `(9)` is worth recalling.  In the THM-3581 compiler,
the source lines
`q=0` and `x=0` map isomorphically to `L_c` and `L_e`.  Composing a Darboux
map with the compiler gives a noninjective planar Keller map.  Gwozdziewicz's
one-line theorem says its restriction to neither source line can be
injective.  The Darboux map is etale, so both restricted polynomial curves
are immersed.  A constant restricted coordinate forces the other one to be
linear, while a linear coordinate is already injective; both are impossible.
This gives degree at least two for both coordinates on both arms.

The map/loss ledger is therefore

```text
source:       the two affine source lines q=0 and x=0
target:       the central c-arm and e-arm
map:          the THM-3581 compiler, isomorphic on each line
preserved:    polynomial restriction and immersion
lost:         inverse sheets, witnessed by the central d-fold collision
sidecar:      exact positive and negative weight surviving on each arm
cheap test:   degree 0 or 1 on either arm.                         (13)
```

## 3. The simple-arm unit channel

Let `beta` be any simple root of `Sigma` and use `z=b-beta`.  Consider one
complementary channel `u+v=-2`.

If `(u,v)=(-1,-1)`, both coefficients lie in `(z)`, so their Wronskian lies
in `(z)`.  If the weights are `(-2,0)` or `(0,-2)`, formula `(4)` again has
the negative coefficient as an undifferentiated factor and lies in `(z)`.

Every other orientation is, up to reflection,

```text
(u,v)=(-N,N-2),                 N>=3.                 (14)
```

Let the two coefficient orders at `beta` be

```text
a>=ceil(N/3),                   d>=0.                 (15)
```

The first term of their Wronskian has order and nonzero multiplier

```text
a+d-1,                     (N-2)a+Nd.                (16)
```

It is a unit exactly when

```text
N=3,                       a=1, d=0.                 (17)
```

Thus a scalar sum equal to one must contain a channel with weights

```text
(-3,1)                  or                  (1,-3).   (18)
```

The jet-and-arm invoice `(10)--(12c)` already makes `-3` and `1` the middle
weights in their respective three-element supports.  The local calculation
`(14)--(17)` independently verifies that this middle channel is the unique
kind that can supply the scalar arm unit.  Reflecting `(P,Q)` if needed,
write

```text
r=(-3m,-3,A),             s=(-3k,1,B),
m,k>=2,                   A,B>=2.                    (19)
```

The middle edge `(1,1)` is complementary.  The only other possible scalar
edges are `(0,2)` and `(2,0)`.

## 4. Zero or one complementary edge

If `h=0`, the scalar coefficient of the bracket is zero.

If `h=1`, the single scalar channel itself must have Wronskian one.  Section
3 reduces it to weights `(-3,1)`.  Its negative coefficient `f` is divisible
by `Sigma`, and

```text
W_(-3,1)(f,g)=f'g+3fg'.                               (20)
```

Its leading coefficient is

```text
(deg f+3deg g)lc(f)lc(g),                             (21)
```

which is nonzero, while `deg f>=deg Sigma=3`.  Hence `(20)` has positive
degree and cannot equal one.  Therefore

```text
h notin {0,1}.                                        (22)
```

## 5. Exactly two complementary edges

By reflection, it suffices to take `(0,2)` as the second scalar edge.  Then

```text
B=3m-2,                         A!=3k-2,              (23)
```

and the seven nonscalar edge sums are

```text
L_00=-3m-3k,       L_10=-3k-3,       L_01=-3m+1,
L_20=A-3k,         L_12=3m-5,        L_21=A+1,
L_22=A+3m-2.                                        (24)
```

The first and last are unique.  Among the other five, the only possible
collisions are

```text
C_1: L_20=L_12,       A=3(m+k)-5,
C_2: L_20=L_01,       A=3(k-m)+1,
C_3: L_12=L_21,       A=3m-6.                        (25)
```

At most one can occur.  Indeed `C_1+C_2` would force `m=1`, `C_1+C_3`
would force `3k=-1`, and `C_2+C_3` would give an integer divisible by three
equal to `-7`.

Every singleton nonscalar layer has zero Wronskian.  Since all weights in
`(19),(23)` are nonzero, a zero edge `(i,j)` gives

```text
f_i'/(r_i f_i)=g_j'/(s_j g_j).                        (26)
```

If there is no collision, the singleton edges contain the spanning tree

```text
{00,01,10,12,22}.                                     (27)
```

In the three collision cases, omit both colliding edges.  The remaining
singleton edges contain, respectively, the spanning trees

```text
C_1: {00,01,10,21,22},
C_2: {00,10,12,21,22},
C_3: {00,01,10,20,22}.                                (28)
```

Thus the bipartite graph on the three `P` pieces and three `Q` pieces is
connected in every case.  Equations `(26)` force one common normalized
logarithmic derivative for all six coefficients.  Every Wronskian, including
both scalar channels, is then zero.  This contradicts `{P,Q}=1`.

The case where `(2,0)` is the second edge is the reflected argument.  Hence

```text
h!=2.                                                 (29)
```

## 6. Three complementary edges

Full complementation gives

```text
B=3m-2,                    A=3k-2,                    (30)
```

and

```text
r=(-3m,-3,3k-2),
s=(-3k,1,3m-2).                                      (31)
```

Let the two gaps of `r` be

```text
p=3m-3,                         q=3k+1.               (32)
```

### 6.1 Unequal gaps

If `p!=q`, the six off-complement sums, relative to the scalar sum `-2`, are

```text
-(p+q),          -p,          -q,
 p,               q,           p+q.                  (33)
```

They are distinct, so all six off-matching Wronskians vanish separately.
The graph `K_(3,3)` minus the complementary matching is a connected six-
cycle.  Equation `(26)` again gives a common normalized logarithmic
derivative, making all three complementary Wronskians zero.  Their sum cannot
be one.

### 6.2 Equal gaps

The equal-step boundary cannot occur under the two-arm invoice.  Equations
`(32)` would give

```text
3m-3=3k+1,                    3(m-k)=4,               (34)
```

which is impossible.  Equivalently, beginning with the `-3` middle channel,
an arithmetic progression would be

```text
{-3m,-3,3m-6},
```

whose full complement has unique negative weight `4-3m`, congruent to one
modulo three rather than an exact `-3k` required by `(12)`.

This proves `h!=3`.  Together with Sections 4--5 it proves the theorem.
**QED.**

## 7. Hostiles, symmetry, and failure boundary

1. **One-sided exact-negative near miss.**  If `(12)` is imposed on only the
   output containing `-3`, the equal-step supports

   ```text
   {-9,-3,3},                    {-5,1,7}              (35)
   ```

   survive the support arithmetic and have the local unit channel `(-3,1)`.
   They are not claimed to satisfy the coefficient equations.

2. **No exact-negative invoice.**  If only the positive-arm lower bound is
   retained, the smaller formal near miss is

   ```text
   {-8,-3,2},                    {-4,1,6}.             (36)
   ```

3. **All two-edge collisions occur.**  Exact representatives for `C_1,C_2,
   C_3` are

   ```text
   (m,k,A)=(2,2,7),              (2,3,4),              (3,2,3). (37)
   ```

   Thus the spanning-tree split is load-bearing; a generic-sums argument
   alone is incomplete.

4. **Reflection only.**  The proof uses the bracket-preserving reflection
   `(P,Q)->(Q,-P)` to exchange orientations.  It asserts no geometric `S_3`
   symmetry of the surface or its arms.

5. **Zero weights.**  Once `(9)` and `(18)` are imposed, the six weights have
   the strict low/middle/high form `(19)` and none is zero.  Without the arm
   invoice, zero-weight pieces require a separate scalar-removal audit.

## 8. Exact scope and remaining counterexample frontier

```text
PROVISIONAL / PENDING INDEPENDENT AUDIT:
  no 3 x 3 weight support can carry a polynomial Darboux pair on
  c^3 e=b(b^2+kappa^2), kappa!=0;
  each Darboux coordinate needs at least three weight pieces, so every
  2 x m and m x 2 support cell is empty;

NOT CLAIMED:
  a target-wide Darboux nonentry theorem;
  emptiness of 3 x 4 or larger cells;
  existence or nonexistence of a polynomial Darboux pair on the target;
  a planar Jacobian counterexample or proof of JC(2).                 (38)
```

The mechanism is a three-scale pincer: infinity-sheet loss forces nonlinear
curves on both central arms; the arms force exact weights on opposite sides
of each support; and additive singleton layers propagate zero Wronskians
across a connected bipartite graph.  Dropping the exact negative arm sidecar
restores the formal near misses `(35)--(36)`.

## 9. Exact verification contract

Reproduce with

```bash
python3 04-computation/jc2_a13_exponent_three_three_by_three_nonentry_thm3590.py
python3 -O 04-computation/jc2_a13_exponent_three_three_by_three_nonentry_thm3590.py
```

The assertion-independent exact companion checks:

- the exponent-three bracket and negative-weight normal form;
- the central tangent-jet determinant and three-piece-per-output invoice;
- the complete simple-arm unit classification through an exact hostile atlas;
- the `h=0,1,2,3` matching split and both orientations;
- every possible two-edge collision and all four singleton spanning trees;
- the unequal-gap six-cycle and equal-gap congruence obstruction;
- the formal near misses `(35)--(36)` and collision representatives `(37)`;
- a bounded exhaustive parameter atlas `2<=m,k<=12`, `2<=A<=80`.

The finite atlas is a hostile control.  The theorem is the all-parameter
valuation, congruence, and connected-graph proof above, not an extrapolation
from the bounded range.
