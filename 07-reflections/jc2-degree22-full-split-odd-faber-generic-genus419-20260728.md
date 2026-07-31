# The full split odd-Faber deformation has generic normalization genus 419

**Status:** PROMOTED AS PROVED THM-2719 + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.  This is a generic theorem about the
full chosen-sheet coefficient family.  It does not exclude exceptional
coefficient strata, close the split degree-twenty-two branch, or prove
`JC(2)` or `DC(2)`.

Companion:

```text
04-computation/jc2_degree22_full_split_odd_generic_genus419_scout.py
05-knowledge/results/jc2_degree22_full_split_odd_generic_genus419_scout.out
script LF SHA-256 59852a54589b4b3f8b68e6aba9a4f3d7d022032244ccabbe0875791a99b5bea9
output LF SHA-256 18e8346b3cfa8858a456318ba836de433f921bacc54704697400957986ef7611
```

Reproduce with ordinary and optimized Python.  Both paths use explicit
`require` exceptions and have byte-identical output.

## 1. Inheritance and the missing family

THM-2704 proves normalization genus `89` on the nonzero-first-flux
**even-Faber** section, where all eleven odd seeds vanish.  THM-2230 and
THM-2202 show that the full target-translated degree-twenty-two Faber gauge
on one split sheet is

```text
Q=E_22+sum_(j in J) a_j E_j,

J={1,2,3,5,6,7,9,10,11,13,14,15,17,19,21}.          (1)
```

The `E_18` coefficient has been killed by the legal translation `P -> P+c`;
the degrees divisible by four are target shears.  The even members of `(1)`
are THM-2411's `B,C,D,E`, while the eleven odd coefficients are the omitted
split bank.  On the two split sheets an odd coefficient is anti-invariant;
after choosing one sheet it is an ordinary complex scalar `a_j`.

The canonical hostile to an elimination claim is the new exact degree-six
kernel

```text
128 E_5+32 E_3+E_1.                                  (2)
```

Boundary, first flux, second flux, the fourth Hamiltonian row modulo its
free constant, and its first relative differential do not kill `(2)` on the
first quadratic augmentation grade.  Thus the odd bank should be deformed,
not silently deleted.

## 2. Exact weighted homogeneity of every odd direction

Write the depressed quartic as

```text
P=w^4+2d w^2+q w+(d^2-s).                            (3)
```

Assign the intrinsic quartic weights

```text
wt(d,q,s)=(2,3,4).                                   (4)
```

For every Faber degree `j`, if `Phi_j,Psi_j` are THM-2129's first two
Laurent observables, the recurrence gives the identities

```text
wt(Phi_j)=j+1,                  wt(Psi_j)=j+2,        (5)

Phi_j(d,-q,s)=(-1)^(j+1)Phi_j(d,q,s),
Psi_j(d,-q,s)=(-1)^j Psi_j(d,q,s).                   (6)
```

These statements are termwise polynomial identities for all sixteen degrees
in `(1)`, including the top degree.  The companion derives every row both
from the Laurent recurrence and from an independent finite multinomial sum;
they are not inferred from a finite specialization.  Formula `(6)` also
types the split constants correctly:
`a_j Phi_j` is anti-invariant and `a_j Psi_j` is invariant for both odd and
even `j`.  Hence on a chosen sheet the two constant-flux equations are

```text
sum a_j Phi_j=lambda,              sum a_j Psi_j=W,  (7)
```

with `a_22=1`, `a_18=0` and the sum otherwise over `J`.

Adjoin `h` of weight one.  The exact homogenizations are

```text
F_23=Phi_22+sum_(j in J) a_j h^(22-j)Phi_j-lambda h^23,

G_24=Psi_22+sum_(j in J) a_j h^(22-j)Psi_j-W h^24.  (8)
```

Thus the full odd bank is one family of `(23,24)` complete intersections in

```text
P(1,2,3,4)_[h,d,q,s].                                (9)
```

Equivalently, before projectivizing, the coefficient weights are

```text
wt(a_j)=22-j,             wt(lambda)=23, wt(W)=24.   (10)
```

This is the promised full-bank extension of THM-2704's prime-23 grading.
The number `23` remains the first-flux weight; odd coefficients do not alter
it.

## 3. Arithmetic genus 425

On the regular-sequence locus, the homogeneous coordinate ring has Hilbert
series

```text
(1-t^23)(1-t^24)
-----------------.                                   (11)
(1-t)(1-t^2)(1-t^3)(1-t^4)
```

Its `a`-invariant is

```text
23+24-(1+2+3+4)=37.                                  (12)
```

The exact degree-`37` coefficient is

```text
511-47-39=425.                                       (13)
```

Graded duality therefore gives

```text
p_a=425.                                             (14)
```

As an independent control, for every multiple of twelve from `60` through
`132`, the companion checks the Hilbert polynomial

```text
H(n)=23n+1-425.                                      (15)
```

## 4. The all-even fibre is an integral specialization

Take

```text
a_j=0 for odd j,
B=C=D=E=W=1,                    lambda=1/7496192.     (16)
```

Then THM-2704's parameter `eta` is one.  On `s!=0`, the exact change

```text
y=11s,       u=dq^2,       Z=q^4,
t=1/y,       v=u/y^2,      zeta=Z/y^3                (17)
```

identifies `(7)` birationally with THM-2704's all-one prime-23 curve.  The
inverse retains the chosen sheet because THM-2704 reconstructs

```text
q=-7496192 lambda t^5/f1.                            (18)
```

That curve is geometrically integral and has normalization genus `89`.
The companion additionally checks two possible lost-component boundaries:

```text
gcd(F_23|_(h=1,s=0),G_24|_(h=1,s=0))=1,
gcd(F_23|_(h=1),G_24|_(h=1))=1,                     (19)
```

while the top forms at `h=0` are coprime.  Hence no component is trapped on
`s=0` or at projective infinity.  The direct weighted fibre `(16)` is itself
geometrically integral.  It lies in the regular-sequence/flat locus and has

```text
delta_total=425-89=336.                              (20)
```

This supplies more than a nonempty fibre: geometric integrality is open in
the resulting proper flat family.

## 5. The unique forced infinity point

At `h=0`, all lower columns disappear.  Up to nonzero constants the two top
forms are

```text
Phi_22=q(-84d^2q^4s+3dq^6+280dq^2s^3-21q^4s^2-84s^5),

Psi_22=-224d^3q^6+3360d^2q^4s^2-336dq^6s
       -3360dq^2s^4+3q^8+560q^4s^3+224s^6.          (21)
```

Exact Groebner elimination of `(21)` together with the three Jacobian
minors gives:

```text
q=1 patch: empty,
s=1 patch: empty,
d=1 patch: support q=s=0.                             (22)
```

Thus every member has exactly one possible singular point at infinity,

```text
P_infty=[h:d:q:s]=[0:1:0:0].                        (23)
```

Every other common top zero is transverse, independently of the lower
coefficients.

The point `(23)` lies on the index-two quotient chart.  On its index-one
cover, the involution is

```text
(h,q,s)->(-h,-q,s).                                  (24)
```

The `E_21` first-flux coefficient at the origin is the nonzero number

```text
Phi_21(1,0,0)=88179/131072.                          (25)
```

Therefore `a_21!=0` lets `F_23=0` solve analytically for `h`.  The top first
flux begins in ordinary order six, so substituting this solution changes no
order-six term of the second equation.  That term is exactly

```text
-231/256 (q-s)(q+s)(q^2-4qs+s^2)(q^2+4qs+s^2).      (26)
```

It consists of six distinct lines on the cover.  Passing to the coarse
invariant `Q=q^2`, `(26)` becomes

```text
-231/256 (Q-s^2)(Q^2-14Qs^2+s^4).                   (27)
```

Hence the coarse curve has three smooth branches

```text
Q=c s^2+O(s^3),        c=1, 7+4sqrt(3), 7-4sqrt(3). (28)
```

Every pair has intersection multiplicity two.  Thus the forced singularity
has the exact generic delta

```text
delta(P_infty)=3*2=6.                                (29)
```

This is the sharp role of the highest odd seed: `E_21` reduces the special
even-section delta budget from `336` to the unavoidable quadratic-deck
quotient delta `6`.

## 6. Generic genus 419

Let `A` be the affine parameter space of the fifteen lower Faber
coefficients together with `lambda,W`.  On `h=1`, the universal incidence
defined by `(7)` is the graph obtained by solving for `lambda,W`; it is
smooth.  The top-only two-flux map has the exact rank-two minor

```text
det partial(Phi_22,Psi_22)/partial(d,q)|_(d,q,s)=(0,1,0)
 =9801/524288 !=0.                                   (30)
```

Thus the incidence projection is dominant.  Generic smoothness gives a
nonempty open parameter locus whose affine fibres are smooth.  Intersect it
with:

1. the open regular-sequence/flat locus containing `(16)`;
2. the open geometrically integral locus supplied by `(16)`; and
3. the open condition `a_21!=0`.

The intersection is nonempty because the parameter base is irreducible.
Equations `(22)` and `(29)` show that on this intersection the only
singularity is `(23)`, with delta six.  Combining `(14)` and `(29)` gives

```text
g(normalization)=425-6=419.                          (31)
```

In particular, the full eleven-dimensional odd bank is not merely compatible
with positive genus: on a nonempty Zariski-open locus it raises the exact
generic genus from the even-section value `89` to `419`.
Consequently no nonconstant rational Keller trajectory can lie in that open
coefficient locus.

## 7. Boundary and cheapest next test

This result is generic, not uniform.  A hypothetical Keller trajectory may
land in the proper discriminant/genus-drop locus.  The all-even point itself
is the sharp warning: its delta jumps from `6` to `336`.  Neither openness nor
lower semicontinuity can exclude such a specialization.

The exact connection ledger is

```text
source:       full split chosen-sheet odd Faber coefficient bank;
target:       weighted (23,24) flux complete intersection;
map:          Laurent observables plus h^(22-j) homogenization;
preserved:    all eleven odd directions, split parity, lambda, W;
destroyed:    third flux, Keller one-form, and exceptional-locus data;
new invariant:
              arithmetic genus 425 and forced mu_2 delta 6;
hostile:      a_21=0 even section, delta 336, genus 89;
needed sidecar:
              prove an actual Keller coefficient point cannot lie in the
              genus-drop locus, or spend the third flux there.
```

The cheapest decisive next computation is the discriminant order transverse
to `a_21=0`: compute the first nonzero `a_21`-term in the `336`-unit delta
conductor of the THM-2704 point, then intersect that initial conductor with
the retained third flux.  A unit initial term would show the even slice is a
high-codimension collision rather than a stable exceptional component; a
vanishing term would identify the next odd coefficient that controls the
genus-drop stratum.
