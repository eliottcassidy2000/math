---
id: THM-3082
title: "Admissible suspension-word simultaneous chambers and scale-tree holotopy"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.  Every
  finite rooted suspension word built from the one-normal operation of
  THM-3069 and the fixed-gap two-normal operation of THM-3063 thickens from
  an ordered limit to an explicit nonempty simultaneous scale chamber.  The
  proof first takes the exact lower-error-ideal quotient and only then pays a
  multislow factorial-entropy invoice.  It gives the closed carrier and a
  strict exponential margin in every fixed width, but not a width-uniform
  angle, maximal fan, arbitrary support, or arbitrary-radial GMC theorem.
source: root-gmc-all-width-multiscale-2026-08-01
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-3063-terminal-suspension-transverse-resultant-and-five-slot-tail-holotopy
  - THM-3069-one-normal-remote-terminal-suspension-and-physical-tropical-flag
  - THM-3073-upper-triangular-resultant-norm-and-torsion-character-death-barcode
  - THM-3075-two-scale-entropy-newton-chambers-and-rooted-scale-tree-holotopy
related:
  - THM-2985-multiparameter-normal-map-and-arc-factor-separation
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3054-affine-moving-lower-tropical-beta-gamma-tail-holotopy
script: 04-computation/gmc_admissible_suspension_word_chambers_thm3082.py
output: 05-knowledge/results/gmc_admissible_suspension_word_chambers_thm3082.out
script_sha256: afc93fa647422010ae4cd9d76abe71efcbbaa503e56f31e65b36ade760f3465d
output_sha256: b8dcd63825011dab488b5272ec7146ee98e1ecfa42119139e8088b4d719d6a13
hash_basis: LF-normalized bytes
---

# THM-3082 -- admissible suspension-word simultaneous chambers

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.**

THM-3069 constructs physical support towers by taking one remote limit after
another.  THM-3063 supplies a second graft that adds a fixed-gap terminal
pair.  An ordered tower does not by itself give one simultaneous ray: lower
coefficients may grow exponentially as soon as their scale ratio is positive.

Nevertheless every finite word in these two grafts has a nonempty
simultaneous chamber.  The order of operations is load-bearing.  First take
the exact associated-graded quotient by the lower-error ideal; only then use
entropy to control its complement.

## 1. Admissible suspension words

Fix

```text
n>=1,             0<=a<b<c,                              (1)
```

and let `S_3=R_3(a,b,c)>0` be the standard-orientation three-slot resultant.
Its positivity is THM-2824.

An **admissible suspension word** is a sequence

```text
3=k_0<k_1<...<k_s=K,       k_i-k_(i-1) in {1,2}.          (2)
```

Choose scales

```text
0<sigma_1<...<sigma_s,
L_i(T)=sigma_i T+O(1) in Z,                               (3)
```

as `T->infinity`.  At a step ending in width `k=k_i`:

* `k_i-k_(i-1)=1` is an `O_k` node: append `L_i(T)`;
* `k_i-k_(i-1)=2` is a `P_k` node: append
  `L_i(T),L_i(T)+h_i`, where `h_i>=1` is fixed.

Thus fixed pair gaps have slope zero.  Put

```text
sigma_0=0,
delta_i=max_(j<i) sigma_j/sigma_i=sigma_(i-1)/sigma_i.   (4)
```

Every resulting support is strictly increasing for large `T`.

For `k>=4`, define the full factorial bill and entropy tax

```text
B_k=(k-1)k!,
J_k(delta)=B_k delta[log(k/delta)+1],
J_k(0)=0.                                                  (5)
```

The node gaps proved by THM-3069 and THM-3063 are

```text
gamma_O(k)=log[k^k/(k-1)^(k-1)],
gamma_P(k)=(k-2)log[(k-1)/(k-2)].                         (6)
```

Require at every node

```text
J_(k_i)(delta_i)<gamma_O(k_i)       for O_(k_i),
J_(k_i)(delta_i)<gamma_P(k_i)       for P_(k_i).          (7)
```

Each condition defines a nonempty interval `0<=delta_i<epsilon_(type,k)`.
Indeed `J_k` is continuous, vanishes at zero, and has derivative
`B_k log(k/delta)>0` on `(0,1)`.  Its endpoint has the closed form

```text
epsilon_(type,k)=
 -gamma_(type,k)/
 [B_k W_(-1)(-gamma_(type,k)/(e B_k k))],                 (8)
```

where `W_(-1)` is the negative real Lambert branch.  Formula `(8)` is only a
description of the aperture; the proof uses `(7)` directly.

## 2. The multislow entropy lemma

For nonnegative slopes `t=(t_1,...,t_r)`, write

```text
Phi(t)=(sum t_j)log(sum t_j)-sum t_j log t_j,              (9)
```

with `0 log 0=0`.  This is the exponential rate of a normalized `r`-factor
factorial atom.

Suppose some coordinates are fast, hence equal to one, some are fixed zero,
and every remaining slow coordinate lies in `[0,delta]`, with
`0<delta<1/e` and `r<=k`.  If `t^0` is obtained by setting all slow
coordinates to zero, then

```text
0<=Phi(t)-Phi(t^0)
  <=r delta[log(k/delta)+1].                              (10)
```

To prove this, add the slow coordinates one at a time.  If the current sum of
the other coordinates is `A`, the derivative of

```text
(A+x)log(A+x)-A log A-x log x
```

is `log((A+x)/x)`, which is nonnegative and at most `log(k/x)`.  Integration
from zero to `x<=delta`, followed by monotonicity of
`x[log(k/x)+1]`, proves `(10)`.

At width `k`, the resultant of forms of degrees `2,...,k` has coefficient
multidegrees

```text
mu_r=k!/r.                                                 (11)
```

Consequently one fully expanded resultant monomial contains exactly

```text
sum_(r=2)^k r mu_r=(k-1)k!=B_k                            (12)
```

factorial factors, with multiplicity.  Summing `(10)` proves that turning on
all child scales can improve the exponential weight of any outer resultant
atom by at most `J_k(delta)` per unit parent scale.  This bound permits
arbitrarily many distinct child slopes; replacing them all by their maximum
is lawful by `(10)`.

## 3. Group before taking the asymptotic

At every outer node, pivot the determinant-one normal coordinates on the
original fixed offset `c`, never on a moving child.  Permuting `c` into the
distinguished position changes coordinate orientation only; the physical
resultant covariance exponent `k!` is even, so no sign survives.

Write each lower equation as

```text
F_r=H_r+E_r,                                               (13)
```

where `H_r` is its pure child restriction and `E_r` contains a new normal
variable.  Let `I_i` be the ideal generated by all coefficients of these
lower errors.  Specializing `I_i=0` is an exact algebraic operation.

For an `O_k` node the lower degrees are `2,...,k-1` and the degree-`k`
equation is arbitrary.  THM-3073 gives the exact quotient

```text
R_k mod I_i = R_(k-1)^k times the one-normal norm.          (14)
```

After the physical scaling of THM-3069, its grouped carrier is

```text
C_(O_k)=R_(k-1)^k U_(k,L_i)^[(k-1)!].                     (15)
```

For a `P_k` node the pure lower degrees are `2,...,k-2` and the upper
degrees are `k-1,k`.  THM-3073 instead gives child exponent `k(k-1)` and
normal exponent `(k-2)!`.  Combining that exact quotient with THM-3063's
physical two-normal scaling gives

```text
C_(P_k)=R_(k-2)^[k(k-1)] E_(k-1,k;L_i,h_i)^[(k-2)!]
         U_(k-1,L_i)^[k!/(k-1)] U_(k,L_i)^[(k-1)!].       (16)
```

Here `E_(k-1,k;L,h)` is the exact intrinsic all-high binary normal
resultant of THM-3063.  It satisfies

```text
E_(k-1,k;L,h)
 =(k^h-(k-1)^h)^[k(k-1)](1+O(L^-1))>0                    (17)
```

for all sufficiently large `L`.  Fixed-layer normal leakage is strictly
smaller than the gap in `(6)` and remains in the complementary bank.  The
exact `E`, including its `O(L^-1)` symbol correction, stays in the carrier;
it is not part of that exponentially small leakage.

Equations `(14)--(16)` are grouped identities.  They make no claim that the
arbitrary upper coefficients converge.  Every resultant term outside the
displayed carrier lies in `I_i` or in the smaller fixed-layer leakage bank.
With the child frozen, THM-3069 bounds it by `exp(-gamma_O(k)L_i)` at an
`O_k` node, and THM-3063 bounds it by `exp(-gamma_P(k)L_i)` at a `P_k`
node.  Turning on every child scale costs at most `(5)`.  Thus its relative
cost is

```text
O(poly(T) exp(-sigma_i[gamma_type(k_i)-J_(k_i)(delta_i)]T)).  (18)
```

There is no hidden factor two in `(18)`: each error atom starts at boundary
rate at most `-gamma_type`, and `(10)--(12)` raise its absolute rate by at
most `J_k`; the exact grouped carrier has nonnegative child-scale rate by
induction because its `U` factors grow, its `E` factors are subexponential,
and `S_3` is fixed.  Hence the relative rate is at most `-gamma_type+J_k`.

The strict inequality `(7)` is therefore the exact sufficient chamber gate.

## 4. Simultaneous chamber and closed carrier

Define

```text
eta=min_i sigma_i[gamma_type(k_i)-J_(k_i)(delta_i)]>0.     (19)
```

Iterate `(15)` and `(16)` from the inside out.  If a child already equals its
carrier times `1+O(poly(T)e^(-eta_child T))`, raising it to the fixed node
degree preserves that relative rate.  Equation `(18)` adds the new rate, so
induction takes their minimum.  The exact exponent recursion telescopes to

```text
R_K(T)=S_3^[K!/6]
       product_(r=4)^K U_(r,L_(node(r)))^[K!/r]
       product_(P_k nodes)
          E_(k-1,k;L_k,h_k)^[K!/(k(k-1))]
       [1+O(poly(T)e^(-eta T))].                           (20)
```

The node containing degree `r` determines `L_(node(r))`; both degrees in a
`P_k` node use the same leading scale.  The exponents follow either from the
induction or directly from resultant multidegrees:

```text
S_3: K!/6,             U_r: K!/r,
E at P_k: K!/[k(k-1)].                                  (21)
```

All factors in `(20)` are positive for large `T`.  Hence the physical
first-window resultant is positive throughout every chamber `(7)` for all
sufficiently large `T`.

There are `F_(K-2)` admissible words from width three to width `K`, where
`F_1=F_2=1`.  They give relative-open fixed-gap strata.  They are not proved
disjoint and are not claimed to be the maximal entropy-Newton fan.

Taking only `O` nodes proves the promised all-width corollary: every finite
lexicographic one-normal tower of THM-3069 thickens to an honest simultaneous
open scale cone.  If positivity of the intermediate single-point prefix
inside every `P_k` node is also desired, impose in addition

```text
J_(k-1)(delta_i)<gamma_O(k-1).                            (22)
```

These finitely many extra inequalities still define a nonempty cone, and
THM-3069 plus the same argument makes every individual first-window prefix
positive.

## 5. A seven-slot chamber

Take the word

```text
O_4, P_6, O_7,                                             (23)
```

with scales `D,C,R`.  The two nontrivial gates are

```text
D/C < epsilon_(P,6)
    =0.00001808125357979315927305509...,
C/R < epsilon_(O,7)
    =0.00000636703607580763243071763... .                 (24)
```

For example, `D/C=1/100000` and `C/R=1/200000` are rigorously safe.  Formula
`(20)` becomes

```text
R_7(T)=S_3^840 U_(4,D)^1260 E_(56,C)^168
       U_(5,C)^1008 U_(6,C)^840 U_(7,R)^720
       [1+O(poly(T)e^(-eta T))]>0.                        (25)
```

This is an honest three-scale seven-slot physical chamber, not an ordered
sequence of limits.

## 6. Two sharp failure boundaries

Grouping before estimation is essential.  In the binary homogeneous system

```text
H=y^d,                    G=z^e+M_T y^e,                  (26)
```

the upper coefficient `M_T` may diverge arbitrarily while

```text
Res(H,G)=1.                                                   (27)
```

It disappears only after the triangular quotient.  This is the universal
algebraic version of THM-3075's physical pure-upper divergence hostiles.

The strict margin is equally essential.  For

```text
F=y+epsilon_T z,          G=z+c_T y,                      (28)
```

the lower-error quotient has resultant one, but

```text
Res(F,G)=1-epsilon_T c_T.                                  (29)
```

Choosing `epsilon_T=e^(-gamma T)` and `c_T=e^(gamma T)` cancels the carrier
identically at the zero-margin boundary.  Beyond `(7)`, an associated-graded
face alone is insufficient; one must compute the grouped signed entropy bank
or the iterated normal cone of THM-2985.

## 7. Scale-tree holotopy, phase, and scope

Normalize the root scale and use the edge ratios
`q_i=sigma_(i-1)/sigma_i` as coordinates.  The strict inequalities `(7)`
form a product, hence star-shaped, rate chamber; contracting an edge
`q_i->0` recovers the corresponding ordered bracketing.  On every compact
strict subcone the estimates are uniform.  This is a rate-space rooted-tree
holotopy, not a claim of a literal continuous path through integer supports.
It preserves the grouped carrier, resultant multidegrees,
complete-intersection lengths, fixed-gap normal determinants, and positivity.

The word itself is an order-one sidecar.  Every word has the same tropical
`S_3` and `U_r` exponent ledger, so the exponential rate forgets its
parenthesization.  Its `P_k` positions and gaps are remembered by the exact
`E_(k-1,k)` factors in `(20)`; each has asymptotic gap-symbol exponent `K!`.

For the physical factorial family the `U` and `E` factors are positive.  If
one allows their phases to vary abstractly, the same exponent ledger is a
weighted norm and THM-3077's augmentation/carry sidecar is still necessary.
Scalar positivity alone must not be read as phase reconstruction.

The apertures in `(8)` shrink roughly like
`gamma/[B log(Bk/gamma)]`, hence superfactorially with width.  There is no
width-uniform angle.  This theorem also does not cross a chamber wall, prove
arbitrary equal-scale or arbitrary supports, a maximal fan, a uniform remote
threshold, an all-order Stieltjes tail, arbitrary-radial GMC(2), NC2,
LRC(14), JC(2), or DC(2).

## 8. Exact evidence

The exact companion verifies:

1. `B_k=(k-1)k!` for every width `4<=k<=14`;
2. all `376` admissible words through width fourteen, including all `144`
   words at width fourteen, against the closed exponent ledger `(21)`;
3. rigorous rational logarithm certificates for the `P_6` and `O_7`
   thresholds and every rational safe ray in `(24)`;
4. the complete `O_4,P_6,O_7` carrier `(25)`;
5. forty-eight exact Sylvester-resultant divergence hostiles `(26)--(27)`;
   and
6. four exact zero-margin cancellation hostiles `(28)--(29)`.

Run

```text
python 04-computation/gmc_admissible_suspension_word_chambers_thm3082.py
python -O 04-computation/gmc_admissible_suspension_word_chambers_thm3082.py
```

Both modes must equal the stored transcript after LF normalization.

**QED, pending independent audit and status promotion.**
