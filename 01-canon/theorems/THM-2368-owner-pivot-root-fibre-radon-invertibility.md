---
id: THM-2368
title: "Owner-pivot root-fibre Radon invertibility"
status: >
  PROVED + VERIFIED-EXACT; CANDIDATE UNDER INDEPENDENT AUDIT. On every
  canonical first-depth-one owner row, choose the two owner-pivot grafts
  among the ordinary thirteen-units. On almost every thirteen-root fibre,
  the two graft words and the residual guard/three-unit branch word are
  nonempty proper Boolean words on C_13. The residual word has support
  between 3 and 10. Every nonempty proper rational Boolean word on C_13
  is a unit of Q[C_13], and the exact two-axis cyclic-correlation transform
  factors into three nonzero cyclotomic factors. Thus all 169 target
  Fourier modes of the root-fibre Radon table are nonzero, with integer
  directional, transverse, and mixed Dirichlet gaps at least 2/169.
  This pointwise theorem does not force THM-2365 bare H-drift: an exact
  C_13^2 nonnegative hostile has all 169 pointwise modes but integrates to
  H(r,s,t)=(12/169)1_(r!=t), hence zero drift. The missing coordinate is
  the rooted cyclotomic phase/event word through blocker multiplication
  and y-integration. No scalar row is excluded and LRC(14) remains open.
source: codex-2026-07-25-owner-pivot-root-fibre-radon
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
related:
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy
script: 04-computation/lrc14_owner_pivot_root_fibre_radon_thm2368.py
output: 05-knowledge/results/lrc14_owner_pivot_root_fibre_radon_thm2368.out
script_sha256: 8cdefa02bf1c11d9afb967195b1edc54ec3762db326748c728c98f7054dbdcb5
output_sha256: 479d52cb3884eb548a8c41c0ed752e3b87fc7bc318233c694773032f7db4afc6
hash_basis: working-tree bytes (LF)
---

# THM-2368 -- every canonical root-fibre target mode survives

**PROVED + VERIFIED-EXACT; CANDIDATE UNDER INDEPENDENT AUDIT.**

The owner-pivot target action has two balanced dipole axes. At the
thirteen-root scale, each axis is a cyclic translate of a one- or two-hole
Boolean word. The remaining four root-varying factors leave at least three
and at most ten branches. This makes the root-fibre transform much more
rigid than an arbitrary `13 x 13` table:

```text
every one of its 169 target characters is nonzero.                 (1)
```

The statement is pointwise in the base coordinate. THM-2365 instead
multiplies by moving blocker masks and then integrates that coordinate.
The distinction is load-bearing: Section 7 gives an exact nonnegative
hostile in which (1) holds at every point but the integrated tensor is
exactly circulant.

## 1. Canonical owner and root notation

Put

```text
p=13,
d(x)=1_(||x||<1/14),
g=1-d,
gamma(x)=1_(||x||>=1/7).                            (2)
```

Endpoints are assigned arbitrarily and consistently; all statements below
are unchanged off a finite null set.

Use any canonical first-depth-one scalar row. Its six nonblocker labels are

```text
U={H,q_1,...,q_5},                                  (3)
```

where `H,q_i` are units modulo thirteen, `gamma(Hx)` is the guard-safe
factor, and `g(q_i x)` are the five ordinary safe factors. Let

```text
B={j,a,c}
```

be the three blocker labels, all divisible by thirteen. The selected
exclusive owner has indicator

```text
1_(E_j)(x)
 =gamma(Hx) product_(i=1)^5 g(q_i x)
  d(w_j x)g(w_a x)g(w_c x).                        (4)
```

THM-2309 permits an omitted unit `u_0 in U` and two distinct graft rows.
Choose the graft labels

```text
k_a,k_c in {q_1,...,q_5} minus {u_0}.               (5)
```

There are always at least four available ordinary labels if `u_0` is
ordinary, and all five if `u_0=H`. THM-2350 identifies the quotient-dual
axes, in omitted-unit gauge, as

```text
eta=e_a-e_(k_a),
ell=e_c-e_(k_c).                                    (6)
```

For compatibility with THM-2365, use its minus-shift convention

```text
E_(s,t)(x)
 =product_i I_i(w_i x-(s eta_i+t ell_i)/p).         (7)
```

Thus the target blocker factors move by `-s,-t`, while their graft safe
factors move by `+s,+t`.

Write

```text
w_j=13A_j,       w_a=13A_a,       w_c=13A_c.        (8)
```

For `y in T` and `h in F_13`, the thirteen inverse branches are

```text
x=(y+h)/13.
```

The blocker factors in (8) are independent of `h`. Define the three
root words

```text
A_y(r)=g((k_a y+r)/13),
C_y(r)=g((k_c y+r)/13),                             (9)

P_y(h)
 =gamma(H(y+h)/13)
  product_(
    q in {q_1,...,q_5} minus {k_a,k_c}
  ) g(q(y+h)/13).
```

The exact two-axis root-fibre Radon table is

```text
Q_(s,t)(y)
 =1/13 sum_(h in F_13)
   A_y(k_a h+s) C_y(k_c h+t) P_y(h).               (10)
```

All arguments in (9)--(10) are read modulo thirteen where appropriate.
Changing the orientation gauge in (6) only replaces `s` or `t` by its
negative and does not alter any conclusion.

## 2. The canonical capacity is strictly between empty and full

A translated thirteen-grid meets an open circular interval of length
`1/7` in either one or two points. It meets one of length `2/7` in either
three or four points. Indeed the two counts are respectively

```text
floor(13/7) or ceil(13/7),
floor(26/7) or ceil(26/7).                          (11)
```

Multiplication by any unit modulo thirteen only permutes the grid.
Consequently, for almost every `y`,

```text
|supp A_y|,|supp C_y| in {11,12}.                  (12)
```

The complement of `P_y` is covered by one guard-danger word, of size at
most four, and three ordinary danger words, each of size at most two.
Therefore

```text
|supp P_y|>=13-(4+3*2)=3.                           (13)
```

On the other hand, `P_y` is contained in the guard-safe word, whose
complement has at least three points. Hence

```text
3<=|supp P_y|<=10.                                  (14)
```

Both inequalities are sharp from the count data alone: four guard holes
and three disjoint pairs leave three roots, while three guard holes
containing the three singleton unit holes leave ten. The theorem uses only
(14), not realizability of either equality pattern by one scalar row.

## 3. Prime-cyclotomic group-ring units

Let `C_13=<z>` and use convolution in `Q[C_13]`. If

```text
f=sum_(r=0)^12 f_r z^r,             f_r in {0,1},   (15)
```

is nonempty and proper, then `f` is a unit of `Q[C_13]`.

To prove this, use

```text
Q[C_13] isomorphic Q x Q(zeta_13),                  (16)
```

where the two components are evaluation at `1` and at a primitive
thirteenth root `zeta`. The first evaluation is `|supp f|`, which is
nonzero. If `f(zeta)=0`, the degree-at-most-twelve polynomial in (15) is
divisible by

```text
Phi_13(X)=1+X+...+X^12.
```

Its thirteen Boolean coefficients would then all be equal, making `f`
empty or full. This is a contradiction. Every conjugate `f(zeta^lambda)`,
`lambda!=0`, is nonzero as well. Thus both components in (16) are
invertible.

This proves, in particular, that `A_y,C_y,P_y` are units for almost every
`y`. It also gives the exact four canonical one/two-hole kernels. Up to a
cyclic rotation and a group automorphism, they are

```text
delta_0,
delta_0+delta_1,
1-delta_0,
1-delta_0-delta_1.                                  (17)
```

Their absolute convolution determinants are respectively

```text
1, 2, 12, 11.                                       (18)
```

Writing `S=sum_(r=0)^12 delta_r` and

```text
V=1/2(delta_0-delta_1+delta_2-...+delta_12),
```

explicit inverses are

```text
delta_0,
V,
-delta_0+S/12,
-V+S/22.                                           (19)
```

For example,

```text
(delta_0+delta_1)V
 =1/2(delta_0+delta_13)=delta_0.
```

Thus no canonical one- or two-root endpoint mask loses information under
**cyclic convolution**. Section 7 explains why this does not make
pointwise multiplication followed by integration invertible.

## 4. Exact Radon factorization

Use normalized transforms

```text
A_y^(lambda)
 =1/13 sum_r A_y(r)zeta^(-lambda r),

C_y^(mu)
 =1/13 sum_r C_y(r)zeta^(-mu r),

P_y^+(nu)
 =1/13 sum_h P_y(h)zeta^(nu h),                    (20)

Q_y^(lambda,mu)
 =1/13^2 sum_(s,t)
    Q_(s,t)(y)zeta^(-lambda s-mu t).
```

Substitute (10), then set

```text
r=k_a h+s,          q=k_c h+t.
```

The two finite sums separate exactly:

```text
Q_y^(lambda,mu)
 =A_y^(lambda) C_y^(mu)
  P_y^+(lambda k_a+mu k_c).                        (21)
```

Every factor on the right is nonzero by Sections 2--3, including its zero
character because each word has positive mass. Therefore, for almost every
`y`,

```text
Q_y^(lambda,mu)!=0
       for every (lambda,mu) in F_13^2.             (22)
```

This is all `169` target modes, not merely the thirteen modes on one affine
line.

Equation (21) is also an exact invertibility statement. If `A_y,C_y` are
known, `Q_y` determines every coefficient of `P_y`: for a requested `nu`,
take for instance

```text
lambda=nu k_a^(-1),        mu=0,
```

and divide (21) by its two nonzero graft factors. Fourier inversion then
recovers the complete rooted branch word. Thus the Radon table loses no
root-fibre information once the two labelled graft words are retained.

The one-axis specialization

```text
Q_s=1/13 sum_h A(kh+s)P(h)
```

has the analogous factorization

```text
Q^(lambda)=A^(lambda)P^+(lambda k),                 (23)
```

so all thirteen one-axis characters survive.

## 5. Uniform integer Dirichlet gaps

Put

```text
N_(s,t)(y)=13Q_(s,t)(y) in {0,1,...,13}.            (24)
```

By (22), none of the three arrays

```text
Delta_s N,
Delta_t N,
Delta_s Delta_t N                                  (25)
```

is identically zero. Each is integer-valued and has total sum zero on the
torus. A nonzero integer array of total sum zero has at least two nonzero
entries, each of squared magnitude at least one. Hence, pointwise in every
regular root fibre,

```text
sum_(s,t)|Q_(s+1,t)-Q_(s,t)|^2>=2/169,

sum_(s,t)|Q_(s,t+1)-Q_(s,t)|^2>=2/169,

sum_(s,t)|
 Q_(s+1,t+1)-Q_(s+1,t)-Q_(s,t+1)+Q_(s,t)
|^2>=2/169.                                        (26)
```

The one-axis table in (23) similarly obeys

```text
sum_s|Q_(s+1)-Q_s|^2>=2/169.                        (27)
```

The constant is the sharp lattice gap for a denominator-thirteen cyclic
array. No positive lower bound for an **integrated cyclotomic mode** follows
from (26): phases from different `y`-chambers may cancel exactly.

## 6. Exact map to THM-2365's bare tensor

THM-2365 defines

```text
Delta_r(x)=d(w_c x-r/13),

H_E(r,s,t)=integral_T E_(s,t)(x)Delta_r(x) dx.      (28)
```

Under `x=(y+h)/13`, all blocker factors are independent of `h`, and (28)
becomes exactly

```text
H_E(r,s,t)
 =integral_T M_(r,s,t)(y)Q_(s,t)(y) dy,             (29)

M_(r,s,t)(y)
 =d(A_j y)
  g(A_a y-s/13)
  g(A_c y-t/13)
  d(A_c y-r/13).
```

The signs in (29) are those of THM-2365. Its graft factors have the
opposite signs and are already inside (10).

Equation (22) says that `Q_y` has full target spectrum before the
multiplier `M` is applied. But multiplication by the two shifted blocker
masks is convolution in target Fourier space, and integration in `y`
then sums cyclotomic numbers with changing phases. Neither operation
preserves coordinatewise nonvanishing.

## 7. A minimal exact circulant hostile

The logical boundary is realized by a two-coordinate finite probability
space, not only by a signed formal array. Let

```text
Omega=F_13 x F_13
```

with normalized counting measure, and on `F_13` put

```text
d_0=delta_0,             g_0=1-delta_0.
```

For `(u,v) in Omega`, define the nonnegative pointwise root table and
blocker multiplier

```text
q_(s,t)(u,v)=g_0(u+s)g_0(v+t),

m_(r,s,t)(u,v)
 =g_0(u+s)g_0(v+t)d_0(v+r).                        (30)
```

For every fixed `(u,v)`, the two-dimensional Fourier transform of `q` is
the product of two transforms of translates of `g_0`. Since

```text
g_0^(0)=12/13,
g_0^(lambda)=-zeta^(lambda rho)/13,  lambda!=0,
```

all `169` pointwise modes are nonzero.

Nevertheless the integrated nonnegative tensor

```text
H(r,s,t)=E_(u,v)[q_(s,t)(u,v)m_(r,s,t)(u,v)]
```

is

```text
H(r,s,t)
 =12/169 * 1_(r!=t).                                (31)
```

Indeed the `u`-average contributes `12/13`; the `v`-average of
`g_0(v+t)d_0(v+r)` is zero on `r=t` and `1/13` otherwise. Thus

```text
H(t,s,t)=0,
H(r,s,t)=G(r-t),
D_H=0.                                              (32)
```

This hostile has exactly THM-2365's diagonal zero and circulant residual
form. It proves:

```text
full pointwise 169-mode root spectrum
  does NOT imply positive bare H-drift.             (33)
```

The same failure persists inside the exact Radon class of Section 4.
Let `R(s,t)` be any one of its nonnegative full-mode tables and average
all of its target translates:

```text
q_(s,t)(u,v)=R(s+u,t+v).
```

For any fixed Boolean word `b`, use

```text
m_(r,s,t)(u,v)
 =b(s+u)g_0(t+v)d_0(r+v).
```

After the change of variables `U=s+u,V=t+v`, the integral is

```text
1/13^2 sum_(U,V)
 R(U,V)b(U)g_0(V)d_0(V+r-t),                       (33a)
```

which depends only on `r-t` and vanishes when `r=t`. Choosing the
lower-capacity support-three Radon control from the companion makes all
twelve off-diagonal values positive. Hence even translates of an **exact
canonical-capacity Radon table** can form a zero-drift orbit packet.

The independent `u,v` product model is not asserted to be realized by one
canonical scalar row. Its role is sharper: it proves that Sections 1--5
alone cannot exclude the residual. Any canonical exclusion must use a
coordinate destroyed by the passage from the rooted fibres to their
integral.

This is the same alignment obstruction isolated abstractly in THM-2218:
an aligned packet and a packet running through all cyclic translates have
the same nonzero local spectra, while the latter averages every
nonconstant mode to zero. A tournament shadow on root branches or target
cells cannot repair this. There is no intrinsic binary relation, ties are
unavoidable, and orientation forgets the amplitudes, cyclotomic phases,
root gauges, and chamber weights which distinguish alignment from a full
rotating packet.

## 8. The exact phase/event sidecar

The hostile identifies a finite missing coordinate rather than an
unstructured analytic mystery. Every factor in (9), (10), and (29) is a
rational interval indicator. There is therefore a finite rational
partition

```text
0=y_0<y_1<...<y_L=1                                (34)
```

such that all rooted Boolean words and all blocker-shift words are constant
on each open chamber `(y_i,y_(i+1))`.

At a generic endpoint event, one named factor toggles one rooted residue
`rho` with sign `epsilon in {+1,-1}`. Its normalized cyclotomic transform
changes by

```text
epsilon zeta^(-lambda rho)/13.                      (35)
```

For a residual product word, the toggle survives precisely when every
other labelled factor is safe at that root. Simultaneous endpoint events
are retained as a labelled finite set of such toggles; they must not be
silently ordered.

Combining (21) and (35), a graft-`A` event changes a target mode by

```text
epsilon/13 zeta^(-lambda rho)
 C_y^(mu)P_y^+(lambda k_a+mu k_c),                 (36)
```

with the analogous formulas for `C` and `P`. The target blocker events in
`M` are stored in the same rooted gauge. Consequently the complete exact
finite sidecar for (29) is

```text
for each chamber:
  rational length y_(i+1)-y_i;
  labelled rooted words A_i,C_i,P_i;
  the four blocker words in M;
  their cyclotomic phase vectors;

for each boundary:
  simultaneous labelled (factor,residue,sign) toggles.             (37)
```

Summing chamber length times the resulting constant table reconstructs all
`13^3` rational masses `H_E(r,s,t)` exactly. Equivalently, for each target
character it produces a finite weighted walk in `Q(zeta_13)`; positive
drift is exactly the assertion that some nonzero-target walk has nonzero
endpoint.

THM-2201's labelled Hasse jets are a faithful alternative encoding of each
rooted Boolean snapshot and of the toggles in (37). They do not eliminate
the need for the rooted guard block, event order, winding/phase gauge, and
rational chamber lengths. THM-2218's integral group-algebra carrier is
closer to the present cyclic correlations and likewise shows why relative
phase must be retained. No theorem presently aligns the walk into a common
open half-plane or prevents it from making the full rotation in (30).

A useful sufficient certificate is now explicit. If, for some nonzero
target character, all nonzero chamber contributions can be multiplied by
one unit complex number so that their real parts are nonnegative and one is
strictly positive, then their sum is nonzero and THM-2365 gives positive
bare drift. This **phase-cone certificate** is not proved uniformly.

The cheapest decisive canonical computation is to construct (37) row by
row, quotient by simultaneous translations and factor-label symmetries,
and test whether every zero-drift word is forced into the rotating hostile
class. This is stronger and more diagnostic than recomputing only the
`13^3` terminal masses: a failure records the first cancelling phase event
and the missing anchor needed for a repair.

## 9. Consequence and boundary

The theorem proves a complete local invertibility statement on every
canonical owner-pivot root fibre. It rules out:

- a missing local target character;
- an empty or full residual branch word;
- a convolutional kernel in any canonical one/two-hole mask; and
- a locally flat target table.

It does **not** prove:

- `D_(H_E)>0` in THM-2365;
- a uniform positive lower bound after `y`-integration;
- positive drift for a delayed terminal word when the bare owner is
  circulant;
- a landing on the previously marked triangle;
- an all-coordinate `91`-unit address; or
- exclusion of any of the `165` scalar rows.

The exact remaining implication is therefore

```text
canonical arithmetic of the phase/event word (37)
  -> exclusion, or forced word-breaking, of the rotating circulant class.
                                                                  (38)
```

LRC(14) remains open.

## 10. Exact companion

The dependency-free companion verifies:

- all `8190` nonempty proper Boolean words and all `98,280` nonzero
  cyclotomic character tests;
- the four group-ring inverses and determinants `1,2,12,11`;
- the exact one/two ordinary and three/four guard root counts;
- the residual support interval `[3,10]`;
- the Radon factorization and all `169` modes on lower- and upper-capacity
  controls over an independent finite-field chart;
- the three denominator-thirteen Dirichlet gaps; and
- all `28,561` pointwise-mode and `2,197` orbit checks in the exact
  circulant hostile, plus a support-three exact Radon orbit hostile with all
  twelve off-diagonal values positive.

Run

```bash
python3 04-computation/lrc14_owner_pivot_root_fibre_radon_thm2368.py
python3 -O 04-computation/lrc14_owner_pivot_root_fibre_radon_thm2368.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_owner_pivot_root_fibre_radon_thm2368.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python. QED.
