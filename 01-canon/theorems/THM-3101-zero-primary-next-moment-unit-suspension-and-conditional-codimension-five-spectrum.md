---
id: THM-3101
title: "Scheme-zero-primary next-moment-unit suspension and conditional codimension-five spectrum"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-REAUDITED AFTER SECOND
  REPAIR.  A fixed zero-resultant physical child in the finite
  complete-intersection norm branch is healed by every
  sufficiently remote arbitrary-gap cluster when the vanishing moment is
  scheme-theoretically zero on its full zero-primary algebra and the next
  moment is a unit there.  Reduced/etale zero-primary children automatically
  satisfy the first condition.  The exact repair scale is the unique
  degree-m all-high scale; a generalized-alternant row replacement makes the
  repair form a unit on every composition face, and conjugate Artin norms
  make its coefficient positive.  For one appended slot, failure of the unit
  condition gives an exact persistent-null converse.  Under the repaired
  sidecar for every minimal bad prefix, the bad-support count improves
  conditionally from codimension four to codimension five.  The merely
  nilpotent nonreduced branch remains OPEN and requires weighted Rees/jet
  control.
source: zero-child-higher-jet-2026-08-02
audit: >
  The first hostile audit rejected the merely nilpotent version via the
  exact family C[e]/(e^2-epsilon), where Norm(e+tau)=tau^2-epsilon and an
  o(tau) lower deformation can dominate and reverse the proposed leading
  term.  The repaired theorem requires h_m A_0=0.  An independent re-audit
  accepted the continued Fitting summand, uniform o(tau) operator estimate,
  upper finite-free norm, row-replacement unit, conjugate-pair sign, one-slot
  converse, and codimension-five count.  A later audit caught that the
  perturbed lower quotient is finite free over the normal-coordinate ring,
  not a finite algebra over R.  The repaired proof now uses the lifted
  lower-ideal identity for h_m, compact upper-root reduction, and the fact
  that pure upper deformation cannot alter that identity.  This second
  repair was independently accepted after rederiving finite freeness over
  R[v], bounded upper reduction, the localized spectral projector, and the
  resulting uniform o(tau) operator estimate.  Normal, optimized, and stored
  exact companions agree; both LF hashes and documentation checks pass.
depends_on:
  - THM-3019-sfc-integral-handle-non-real-locus-and-extended-census
  - THM-3073-upper-triangular-resultant-norm-and-torsion-character-death-barcode
  - THM-3085-multi-normal-fixed-gap-cluster-and-unconditional-all-width-tail
  - THM-3093-arbitrary-gap-remote-cluster-monge-flag-compactification
  - THM-3097-translated-support-monge-compactification-and-cofinite-bad-set-induction
related:
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
  - THM-2960-local-smith-jet-fitting-barcode-and-negative-depth-chamber-atlas
  - THM-2985-multiparameter-normal-map-and-arc-factor-separation
  - THM-3099-finite-gamma-sum-quasianalyticity-and-remote-resultant-zero-dichotomy
script: 04-computation/gmc_zero_primary_next_moment_unit_suspension_thm3101.py
output: 05-knowledge/results/gmc_zero_primary_next_moment_unit_suspension_thm3101.out
script_sha256: 12a589ed47ca4e3860a0d4ba7144979c70ec5f2b789ec744acdb386369bd62dc
output_sha256: 353220f07f8da33466ec56fffaeabe9c70fbbed3832d6eaa476fd9665f5c5c97
hash_basis: LF-normalized bytes
---

# THM-3101 -- scheme-zero-primary next-moment-unit suspension

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-REAUDITED AFTER SECOND
REPAIR.**

THM-3093 suspends every fixed child whose first-window resultant is nonzero.
The leading carrier vanishes when the child resultant is zero, but that does
not make every remote extension singular.  In the finite-norm branch the
correct replacement is not one scalar derivative.  It is the full
zero-primary child algebra, together with the requirement that the next
moment be a unit on every one of its local factors.

That sidecar exposes one more physical Newton layer.  Its unique lower
response is a unit on the upper Monge algebra.  When the old vanishing moment
is zero as an element of the zero-primary algebra, every lower deformation is
`o(tau)` there and the unit response supplies the first determinant power.
This gives a uniform arbitrary-gap suspension theorem above the reduced, and
slightly more general scheme-zero, branch and a sharp persistent-null
converse for one appended slot.

## 1. Child algebra and the pointwise next-jet hypothesis

Fix `m>=3` and a real physical `m`-slot factorial support.  Eliminate its
linear first moment by the determinant-one convention of THM-3097.  In the
remaining homogeneous coordinates

```text
y=(y_1,...,y_(m-1))
```

write its restricted moment forms as

```text
H_r(y),                         r>=2,                    (1)
```

and put

```text
S_m=Res_y(H_2,...,H_m).                                 (2)
```

This theorem concerns a zero child:

```text
S_m=0.                                                   (3)
```

Assume that

```text
X=V(H_2,...,H_(m-1)) subset P^(m-2)                    (4)
```

is zero-dimensional.  Since the polynomial ring is Cohen--Macaulay, the
forms in `(4)` are a regular sequence and

```text
length_C X=2*3*...*(m-1)=(m-1)!.                       (5)
```

Choose a real linear form `ell` missing `X`, dehomogenize by `ell=1`, and
let

```text
A=R[y]/(H_2,...,H_(m-1),ell-1).                         (6)
```

Write `h_r` for the class of `H_r` and `M_r` for multiplication by `h_r`.
The Fitting decomposition for `M_m` is

```text
A=A_0 x A_x,
h_m nilpotent on A_0,             h_m a unit on A_x.    (7)
```

For example, one may take

```text
A_0=ker(M_m^[(m-1)!]).                                  (8)
```

with its canonical Fitting idempotent.  Equation `(3)` says `A_0!=0`.

The first load-bearing hypothesis is the **scheme-zero-primary condition**

```text
h_m A_0=0.                                             (9a)
```

This is stronger than the automatic nilpotence in `(7)`.  It holds whenever
`A_0` is reduced, in particular whenever the zero-primary part of `X` is
finite etale over `R`.  It can also hold on a nonreduced factor annihilated by
`h_m`.

The second load-bearing hypothesis is

```text
h_(m+1)|_(A_0) is a unit.                              (9b)
```

Equivalently, under `(9a)`,

```text
det(M_(m+1)|A_0)!=0
 iff V(H_2,...,H_m,H_(m+1))=empty.                     (10)
```

Thus `(9b)` is a pointwise next-moment detector on the entire zero-primary
scheme.  A nonzero trace, one detected component, or one selected root is
strictly weaker.  Both `(9a)` and `(9b)` are required for the suspension
claim below.

The real algebra in `(6)` has no real residue field.  Indeed, for every
nonzero real mean-zero coefficient vector `y`, THM-3019's integral handle
gives

```text
H_2(y)=integral_0^infinity P_y(s)^2 exp(-s) ds>0.       (11)
```

Hence `X(R)=empty`, and the complex local factors of `A` occur in disjoint
conjugate pairs, including their nilpotent structures.

## 2. Uniform arbitrary-gap suspension above the zero child

Fix `q>=1`, and put

```text
p=m+1,                  k=m+q,
r_i=p+i-1               (1<=i<=q),
D_z=product_(i=1)^q r_i=k!/m!.                          (12)
```

For every sufficiently large remote offset `C`, choose an arbitrary distinct
gap vector

```text
0=h_1<h_2<...<h_q                                      (13)
```

with no upper bound and no prescribed dependence on `C`, and append the
physical support offsets

```text
C+h_1,...,C+h_q.                                       (14)
```

Let `R_k(C;h)` be the canonical normalized physical first-window resultant.
Under `(4)`, `(9a)`, and `(9b)`, there is a threshold

```text
C_0=C_0(child,m,q)<infinity                            (15)
```

independent of the entire gap vector, such that

```text
R_k(C;h)>0                         for every C>=C_0.    (16)
```

There is also a quantitative normalized form.  With the positive carriers

```text
U_r(C)=L(f_(n+C)^r)/L(f_n^r)
      =(rn+1)_(rC)/(n+1)_C^r,                          (17)

tau_C=U_m(C)/U_p(C)^(m/p)>0,                           (18)
```

put

```text
Rhat_k(C;h)=R_k(C;h)/
 product_(i=1)^q U_(r_i)(C+h_i)^[k!/r_i].              (19)
```

If

```text
d_0=dim_R A_0,                  nu=d_0 D_z,             (20)
```

then there are positive child-, `m`-, and `q`-dependent constants `c,C'`
such that, uniformly in `(13)`,

```text
c tau_C^nu <=Rhat_k(C;h)<=C' tau_C^nu                 (21)
```

for every sufficiently large `C`.  Gauss multiplication gives the exact
repair scale

```text
tau_C
 ~kappa_(child,m) (m/(m+1))^(mC) C^[1/(2(m+1))],
kappa_(child,m)>0.                                     (22)
```

The bounds `(21)` are non-effective because the composition compactification
and resultant continuity are used by contradiction.

## 3. The unique repair layer and the equality-cell correction

Apply the fixed-child coordinates and arbitrary-gap Monge flag of THM-3093.
Before the internal column potentials, THM-3085 assigns to a physical cell
with `j` actual remote factors and normal degree `ell>=j` the exponential
base

```text
j^j/p^ell                              in a row r<=m,
j^j p^(r-ell)/r^r                      in a row r>=p.   (23)
```

The base-one cells are the fixed child forms, `H_p`, and the full all-high
upper blocks.  Convexity of `j log(j/p)` shows that the largest remaining
base is

```text
rho_m=(m/p)^m<1.                                      (24)
```

There are exactly two equality cells:

```text
lower:  (r,j,ell)=(m,m,m),
upper:  (r,j,ell)=(p,m,m).                            (25)
```

The first is the genuine repair of `H_m`.  The second was easy to miss: it
deforms the first upper equation at the same scale and must not be placed in
an `o(tau_C)` remainder.  Every other nonmodel cell has a base
`sigma_m<rho_m`, uniformly after the Jensen/flag bounds.

Consequently, after passing to any composition-face subsequence, the system
has the coefficientwise form

```text
H_r+o(tau_C),                              2<=r<m,
H_m+tau_C P_m(v)+o(tau_C),

H_p+B_1(v)+O(tau_C),
B_2(v)+o(1),...,B_q(v)+o(1).              (26)
```

Here `P_m` is the limiting degree-`m` all-high lower response.  The
`O(tau_C)` term in the first upper equation is retained in the upper algebra
below.  Formula `(22)` follows from

```text
U_r(C)~kappa_(r,n) r^(rC) C^[-(r-1)/2]
```

and the exact exponent identity

```text
-(m-1)/2 + m(p-1)/(2p)=1/(2p).                         (27)
```

## 4. The lower module and the compact upper quotient

Write the first `m-2` equations in `(26)` as

```text
H_r(y)+delta_(r,C)(y,v),                 2<=r<m,       (28a)
```

where every coefficient of `delta_(r,C)` is `o(tau_C)`, uniformly on the
composition compactification.  Quotienting by `(28a)` gives

```text
mathcal L_C
 =R[v][y]/(H_2+delta_(2,C),...,H_(m-1)+delta_(m-1,C)). (28b)
```

After finitely many fixed affine/Macaulay charts, `mathcal L_C` is a finite
free `R[v]`-module of rank `(m-1)!`.  It is **not** a finite `R`-algebra, and
no such claim is used below.  At the limiting face,

```text
mathcal L_0=A tensor_R R[v].                           (28c)
```

Let `e_0 in A` be the Fitting idempotent of `(7)` and choose polynomial lifts
of `e_0` and the lower ideal.  The scheme-zero hypothesis `(9a)` gives an
exact identity

```text
e_0 h_m=sum_(r=2)^(m-1) a_r(y) H_r(y)                 (28d)
```

in the chosen affine coordinate ring.  In `mathcal L_C`, the same lift says

```text
e_0 h_m=-sum_r a_r(y) delta_(r,C)(y,v).                (28e)
```

The analogous lifted identity for `e_0^2-e_0` makes multiplication by `e_0`
an idempotent up to `o(tau_C)`.  No spectral correction is made yet inside
the infinite `R[v]`-module `mathcal L_C`.

The upper equations have homogeneous normal resultant uniformly bounded
away from zero by THM-3093.  Their roots lie in finitely many bounded normal
charts on the composition compactification, and quotienting `mathcal L_C`
by them gives a finite-dimensional flat algebra `mathcal B_C` of the fixed
total rank.  Reduction in the corresponding fixed monomial bases is
uniformly bounded.  On `mathcal B_C`, let `T_C=M_(e_0)`.  It converges
uniformly on each compact chart to the exact limiting idempotent.  A fixed
small spectral contour therefore corrects `T_C` to a projector `Pi_C` of
rank `nu=d_0D_z`, with `Pi_C-T_C=o(tau_C)`.  Since multiplication by `h_m`
is uniformly bounded, `(28e)` gives

```text
Pi_C M_(h_m)
 =T_C M_(h_m)+(Pi_C-T_C)M_(h_m)=o(tau_C).              (28f)
```

Thus multiplication by `h_m` on the continued zero-primary upper summand
is `o(tau_C)`, uniformly in every arbitrary gap vector.  This is the only
finite-dimensional spectral projection used in the proof.

At a composition-face limit the upper algebra is

```text
mathcal B=
A[v_1,...,v_q]/(h_p+B_1,B_2,...,B_q).                  (28)
```

and only this limiting algebra is asserted to be finite free over `A`, of
rank `D_z`.  The `O(tau_C)` upper equality cell changes the finite quotient,
but cannot create a pure-upper contribution to `(28f)`: if every
`delta_(r,C)` is set to zero, identity `(28d)` already vanishes in
`A_0 tensor R[v]` before imposing any upper relation.  Mixed lower/upper
terms remain `o(tau_C)` on the compact charts.  The repair multiplication
matrix converges to the unit proved below.

## 5. The repair form is a unit on every Monge face

It remains to prove that `P_m` is a unit on

```text
mathcal B_0=mathcal B tensor_A A_0.                    (29)
```

First consider a finite first composition block `I=[1,b]`.  Its line matrix
for `B_1,...,B_b` has consecutive row degrees `p,...,p+b-1`.  Remove the
first row, as in `(28)`, and insert the degree-`m` row belonging to `P_m`.
The new ordered row set is

```text
m,p+1,...,p+b-1.                                       (30)
```

Against the strictly increasing fixed gaps in the block, its determinant is
a positive generalized alternant: an ordinary Vandermonde times a positive
Schur polynomial.  The ordered monomial of coefficient one supplies the same
gap-independent floor as THM-3093.  Every later diagonal block is unchanged.
THM-3073 removes arbitrary forward entries, so

```text
Res(P_m,B_2,...,B_q)!=0                                (31)
```

on every composition face, uniformly up to positive row and column factors.

Scheme-theoretically, if a geometric point of `(29)` killed `P_m`, then
`(31)` together with `B_2=...=B_q=0` would force `v=0`.  The first relation
in `(28)` would then say `h_p=0`, contradicting `(9b)`.  Thus `P_m` is a unit
on `(29)`.  This argument includes `q=1`, when `P_m` is a nonzero scalar
multiple of `v_1^m` and there are no forms `B_2,...,B_q`.

## 6. Scheme-zero-to-unit determinant and positive sign

Use the fixed monomial trivializations from Section 4 on the continued
zero-primary upper summand of rank `nu`.  Equations `(26)` and `(28f)` and
the unit conclusion of Section 5 give

```text
M_(h_m)+tau_C M_(P_m)+o(tau_C)
 =tau_C (M_(P_m)+o(1)),

det=tau_C^nu (det(M_(P_m))+o(1)).                     (32)
```

On `mathcal B tensor_A A_x`, multiplication by `h_m` is already a unit and
has a nonzero limiting determinant.  The smaller lower cells are already
included in `(28a)`, while the repaired upper equality cell in `(25)` only
changes the convergent repair matrix in `(32)`.

Both determinants have positive real sign.  Any real residue of `mathcal B`
would restrict to a real residue of `A`, contrary to `(11)`.  Hence every
complex local factor of `mathcal B` is paired with its conjugate, and the
determinant of multiplication by any real unit is

```text
product_({P,Pbar})
 |det_C(M_u on mathcal B_P)|^2>0.                      (33)
```

The standard physical resultant orientation agrees with this norm
orientation; form permutations contribute only powers of the even total
degree `k!`.  The positive carriers in `(19)` restore the physical scale.

Now take any sequence with `C->infinity` and arbitrary gaps.  THM-3093's
compactness gives a composition-face subsequence.  Equations `(31)--(33)`
give a positive normalized limit after division by `tau_C^nu`, and the
Schur floors prevent that limit from approaching zero.  This proves the two
uniform inequalities `(21)` and hence `(16)`.

## 7. Sharp one-slot converse

For one appended slot, failure of `(9b)` gives exact persistence rather than
merely a missing proof.  If `h_(m+1)` is not a unit on `A_0`, some geometric
child point satisfies

```text
H_2=...=H_m=H_(m+1)=0.                                (34)
```

For every appended exponent `C`, set its new coefficient equal to zero.
The enlarged physical forms restrict exactly to the child forms, so `(34)`
is a common projective zero of the entire `(m+1)`-slot system.  Therefore

```text
R_(m+1)(C)=0                     for every admissible C. (35)
```

Combining `(16)` and `(35)`, the scheme-zero-primary branch has an exact
remote one-slot dichotomy:

```text
h_m A_0=0 and next moment a unit     => eventually positive;
next moment not a unit               => identically zero.  (36)
```

The pointwise word "unit" is sharp.  In a product of two conjugate local
pairs, the element `(1,0)` has nonzero trace but is a zero divisor and leaves
the second pair undetected.  Its multiplication matrix has

```text
trace=2,                    determinant=0.              (37)
```

Thus one scalar higher moment or one nonzero component cannot replace
`(9b)`.

## 8. Conditional codimension-five bad-support spectrum

Use THM-3097's notation `B_t` for bad `t`-slot supports in `[0,X]` and
`C_m` for the finite bank of minimal bad prefixes of width `m`.  Assume,
conditionally, that every `c in C_m`, `m>=4`, satisfies the finite-norm,
scheme-zero-primary, and next-moment-unit hypotheses `(4)`, `(9a)`, and
`(9b)`.

Fix a target width `t`.  If `c` has width `m<t`, then `(16)`, with
`q=t-m`, says that a bad extension must choose its first new support point
from a finite `c,t`-dependent set; the threshold is independent of every
later gap.  Only `t-m-1` tail points remain free.  Prefixes with `m=t`
contribute only their one fixed support.  Therefore

```text
#B_t(X)=O_t(X^max(t-5,0)).                              (38)
```

In particular, for `t>=5`,

```text
#B_t(X)=O_t(X^(t-5)),
density(B_t in the t-slot chamber)=O_t(X^-5).           (39)
```

The already proved finiteness of `B_4` is unchanged.  More generally, if
SFC is known through width `b` and every later minimal bad prefix has the
same next-jet sidecar, then

```text
#B_t(X)=O_t(X^max(t-b-2,0)).                            (40)
```

This is a conditional one-coordinate improvement over THM-3097, not an
enumeration or an emptiness theorem.

## 9. Boundaries and the next higher-jet object

Every hypothesis above is visible in a sharp hostile.

1. If a child witness kills all moments through a target width, setting all
   new coefficients to zero saturates the entire old prefix cylinder.
2. If `(9b)` fails, `(35)` already gives a full one-slot cylinder.
3. A nonzero trace or aggregate can miss a zero-primary component, as in
   `(37)` and the two-component hostile.
4. Nilpotence cannot replace `(9a)`.  In the flat family

   ```text
   A_epsilon=C[e]/(e^2-epsilon),
   h_m=e,                         P_m=1,
   ```

   the multiplication norm is

   ```text
   Norm(e+tau)=tau^2-epsilon.                            (41)
   ```

   Set `tau=s^2` and `epsilon=s^3`.  Then
   `epsilon=o(tau)`, but `(41)` is negative for `0<s<1` and dominates the
   predicted `tau^2` term.  At `s=0`, the limiting algebra is
   `C[e]/(e^2)` and `e` is nilpotent but nonzero.  This is the first failed
   implication in the unrepaired candidate: an `o(tau)` deformation of the
   lower algebra need not be `o(tau^nu)` in its norm.
5. The physical Newton bases already expose the same boundary.  The
   `H_(m-1)` all-high lower base and repair base are

   ```text
   alpha_m=((m-1)/(m+1))^(m-1),
   rho_m=(m/(m+1))^m,
   alpha_m/rho_m^2
    =(1-1/m^2)^m (m+1)/(m-1)>1.                        (42)
   ```

   Thus even a length-two Smith/Jordan block requires a genuinely mixed
   weighted face; the scalar inequality `sigma<rho_m` is insufficient.
   THM-2960's local Smith-jet Fitting barcode is exactly the relevant
   sidecar: it must measure the arc valuation of every characteristic
   coefficient, not only the nilpotent special fibre.
6. If `(4)` has excess dimension, the finite Artin norm and flat rank `D_z`
   used above are unavailable.  A local normal cone or excess-intersection
   sidecar is then required.
7. If `H_(m+1)` is a zero divisor but a later moment detects the residual
   scheme, the correct source is the nested zero-primary flag for
   `h_m,h_(m+1),...`, not one scalar "first nonzero moment."  Its physical
   target is an iterated Rees/Monge normal algebra.  That deeper theorem is
   not proved here.

There is an exact sufficient sidecar for the open nilpotent branch.  On the
continued rank-`nu` zero-primary upper summand, factor out the convergent unit
repair matrix and write the full pre-repair operator as `Q_C`.  If

```text
det(z I+Q_C)=z^nu+a_1(C)z^(nu-1)+...+a_nu(C)            (43)
```

obeys the weighted characteristic-ideal estimates

```text
a_j(C)=o(tau_C^j)                 for every 1<=j<=nu,   (44)
```

then the proof of Section 6 goes through verbatim.  Conditions `(43)--(44)`
are the missing Rees/jet datum; they are not proved for physical zero
children here.

The theorem fixes the child and `q`; it does not handle the merely nilpotent
nonreduced branch, growing width, effective thresholds, enumerate any `C_m`,
prove `SFC(4)`, prove SFC at any new width, or extend arbitrary-radial GMC to
a zero support without both sidecars `(9a)--(9b)`.  It proves no LRC, NC2,
JC, or DC statement.

## 10. Exact evidence

The dependency-free companion verifies:

1. the exact `tau_C` power `1/[2(m+1)]` and all lower/upper exponential
   cells through `m=10`, including precisely the two equality cells `(25)`;
2. positive degree-`m` row-replacement generalized alternants on every
   composition face through rank five;
3. exact scheme-zero-to-unit determinant factorizations, the
   `e^2-epsilon` nilpotent deformation hostile, the physical
   `alpha_m>rho_m^2` boundary, and a two-component aggregate hostile;
4. conjugate-pair norm positivity and an indefinite-real-root sign hostile;
5. the `D_z`, zero-primary length, and higher-moment exponent ledgers for
   one through four normal directions;
6. exact zero-coefficient persistent-root controls; and
7. the binomial finite-difference degree drop giving `(38)--(40)`.

Run

```text
python 04-computation/gmc_zero_primary_next_moment_unit_suspension_thm3101.py
python -O 04-computation/gmc_zero_primary_next_moment_unit_suspension_thm3101.py
```

Both modes must equal the stored transcript after LF normalization.

**QED.**
