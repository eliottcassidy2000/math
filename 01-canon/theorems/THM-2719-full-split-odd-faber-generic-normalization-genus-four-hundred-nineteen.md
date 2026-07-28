---
id: THM-2719
title: "Full split odd-Faber generic normalization genus four hundred nineteen"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The full chosen-
  sheet split degree-twenty-two Faber coefficient family is a weighted
  (23,24) complete intersection in P(1,2,3,4), of arithmetic genus 425.  On
  a nonempty Zariski-open locus with a21 nonzero, its affine part is smooth,
  it is geometrically integral, and its sole infinity singularity has three
  smooth coarse branches of total delta six.  Its normalization therefore
  has genus 419.  The all-even THM-2704 point has genus 89 and delta 336, so
  exceptional genus-drop strata, the third flux, JC(2), and DC(2) remain
  open.
source: odd-faber-deformation-scout/root-2026-07-28
audit: thm2694-full-lift-fibre-scout-2026-07-28
depends_on:
  - THM-2202-uniform-all-degree-quartic-pole-closure
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2704-split-even-prime23-generic-genus-eighty-nine
related:
  - THM-2713-split-prime23-component-divisor-budget-and-perfect-power-normal-form
script: 04-computation/jc2_degree22_full_split_odd_generic_genus419_scout.py
output: 05-knowledge/results/jc2_degree22_full_split_odd_generic_genus419_scout.out
script_sha256: 59852a54589b4b3f8b68e6aba9a4f3d7d022032244ccabbe0875791a99b5bea9
output_sha256: 18e8346b3cfa8858a456318ba836de433f921bacc54704697400957986ef7611
hash_basis: LF-normalized bytes
---

# THM-2719 -- the full split odd-Faber family has generic genus 419

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2704 proves genus `89` only after setting all eleven odd Faber seeds to
zero.  Keeping those seeds changes neither the prime weight `23` nor the
complete-intersection mechanism.  It instead unfolds almost all of the
special fibre's delta invariant: generically only the quadratic-deck
singularity at infinity remains.

## 1. Statement

Write the depressed quartic as

```text
P=w^4+2d w^2+q w+(d^2-s),             wt(d,q,s)=(2,3,4). (1)
```

In the target-translated chosen split sheet put

```text
J={1,2,3,5,6,7,9,10,11,13,14,15,17,19,21},
Q=E_22+sum_(j in J) a_j E_j,           a_18=0.            (2)
```

Let `Phi_j,Psi_j` be the first two Laurent observables of `E_j`, and let
`C_a` be the weighted projective closure of the two constant-flux equations.
Then there is a nonempty Zariski-open subset of the full coefficient space
on which

```text
C_a is geometrically integral,
p_a(C_a)=425,
delta(C_a)=6,
g(normalization(C_a))=419.                               (3)
```

On this open set the affine curve is smooth.  Its only singular point is

```text
P_infty=[h:d:q:s]=[0:1:0:0],                            (4)
```

and the three branches there are smooth.  Consequently no nonconstant
rational Keller trajectory lies in this generic coefficient locus.

This is a generic statement.  It neither excludes the exceptional
genus-drop locus nor closes the split degree-twenty-two branch.

## 2. The full weighted family

The Laurent recurrence and a direct multinomial expansion independently give

```text
wt(Phi_j)=j+1,                    wt(Psi_j)=j+2,
Phi_j(d,-q,s)=(-1)^(j+1) Phi_j(d,q,s),
Psi_j(d,-q,s)=(-1)^j Psi_j(d,q,s).                     (5)
```

Adjoin `h` of weight one.  With `lambda,W` the two flux constants, the exact
homogenizations are

```text
F_23=Phi_22+sum_(j in J) a_j h^(22-j)Phi_j-lambda h^23,
G_24=Psi_22+sum_(j in J) a_j h^(22-j)Psi_j-W h^24.      (6)
```

Thus every member is a `(23,24)` intersection in

```text
P(1,2,3,4)_[h,d,q,s].                                  (7)
```

Equivalently `wt(a_j)=22-j`, `wt(lambda)=23`, and `wt(W)=24`.
Formula `(5)` also checks the chosen-sheet deck typing: the odd coefficient
directions are retained rather than quotient-killed.

The top forms `Phi_22,Psi_22` are coprime.  Hence no specialization of `(6)`
has a common hypersurface factor, and the equations form a regular sequence.
The family is proper and flat with Hilbert series

```text
          (1-z^23)(1-z^24)
------------------------------------.                   (8)
(1-z)(1-z^2)(1-z^3)(1-z^4)
```

Its `a`-invariant is `23+24-(1+2+3+4)=37`.  Exact coefficient extraction
gives

```text
dim (coordinate ring)_37=511-47-39=425.                (9)
```

Equivalently, for every sufficiently large multiple of twelve,

```text
H(n)=23n+1-425.                                        (10)
```

Therefore every member of the regular-sequence family has arithmetic genus
`425`.

## 3. A geometrically integral fibre

Set the odd coefficients to zero and take

```text
a_2=a_6=a_10=a_14=W=1,              lambda=1/7496192. (11)
```

On `s!=0`, the exact change

```text
y=11s,   u=dq^2,   Z=q^4,
t=1/y,   v=u/y^2,  zeta=Z/y^3                         (12)
```

is birational to THM-2704's all-one prime-23 curve; its chosen sheet is
recovered because `lambda!=0`.  That normalization has genus `89` and the
curve is geometrically integral.  Exact gcd checks show that the direct
weighted model `(6)` has no extra component trapped on `s=0`, in the affine
chart, or at `h=0`.  Cohen--Macaulay purity then makes `(11)` a geometrically
integral member of the present family.

Geometric integrality is open in a proper flat family.  Thus it holds on a
nonempty open subset of the full coefficient space.  The special member has

```text
delta_total=425-89=336,                                (13)
```

which is the sharp warning that the later generic statement is not uniform.

## 4. The only possible infinity singularity

At `h=0`, lower columns disappear.  Up to nonzero constants the two forms
are

```text
Phi_22=q(-84d^2q^4s+3dq^6+280dq^2s^3-21q^4s^2-84s^5),

Psi_22=-224d^3q^6+3360d^2q^4s^2-336dq^6s
       -3360dq^2s^4+3q^8+560q^4s^3+224s^6.            (14)
```

Exact Groebner elimination of `(14)` and the three Jacobian minors is empty
on the `q=1` and `s=1` patches.  On `d=1` its support is `q=s=0`.
Consequently `(4)` is the only possible singular point at infinity; every
other top intersection is transverse for every lower coefficient choice.

The `d=1` chart is the index-one cover of an index-two quotient, with

```text
(h,q,s) -> (-h,-q,s).                                  (15)
```

The highest odd direction is transverse:

```text
Phi_21(1,0,0)=88179/131072 !=0.                        (16)
```

If `a_21!=0`, the equation `F_23=0` solves analytically for `h`.  The solution
starts in ordinary order six.  Every lower contribution to `G_24` after this
substitution has order at least seven, so the order-six face is exactly

```text
-231/256 (q-s)(q+s)
           (q^2-4qs+s^2)(q^2+4qs+s^2).                (17)
```

It has six distinct cover lines.  In the coarse invariant `Q=q^2`, `(17)` is

```text
-231/256 (Q-s^2)(Q^2-14Qs^2+s^4).                     (18)
```

Weighted Hensel lifting gives three smooth branches

```text
Q=c s^2+O(s^3),               c=1, 7+4sqrt(3), 7-4sqrt(3). (19)
```

For distinct constants in `(19)`, the branch difference has order two.
Thus every pair has intersection multiplicity two and

```text
delta(P_infty)=3*2=6.                                  (20)
```

This local conclusion holds throughout `a_21!=0`, independently of the
other lower coefficients.

## 5. Generic affine smoothness and genus

On `h=1`, the universal incidence defined by `(6)` is a smooth graph: solve
the two equations for `lambda,W`.  Its projection to coefficient space is
dominant because the top two-flux map has the exact rank-two minor

```text
det partial(Phi_22,Psi_22)/partial(d,q) at (d,q,s)=(0,1,0)
   =9801/524288 !=0.                                   (21)
```

Generic smoothness therefore supplies a nonempty open with smooth affine
fibres.  Intersect that open with the nonempty geometric-integrality locus
from Section 3 and with `a_21!=0`.  The coefficient space is irreducible, so
the intersection remains nonempty.  Sections 4 and 5 then leave exactly the
single singularity `(4)` with delta six.  Hence

```text
g(normalization)=p_a-delta=425-6=419.                  (22)
```

## 6. Scope and next target

The all-even point shows why generic positive genus is not a uniform
exclusion: its delta jumps from `6` to `336`.  A hypothetical Keller point
may still lie in the proper discriminant/genus-drop locus.  The theorem also
does not retain the third flux or Keller one-form, treat the split `y=0`
boundary, or prove `JC(2)` or `DC(2)`.

The next exact target is therefore transverse rather than another genus
count: compute the first nonzero odd-parameter term in the conductor of the
THM-2704 delta-`336` fibre and intersect that initial conductor with the
third flux.  This distinguishes a high-codimension collision from a stable
exceptional Keller component.

## 7. Exact companion

Run

```text
python 04-computation/jc2_degree22_full_split_odd_generic_genus419_scout.py
python -O 04-computation/jc2_degree22_full_split_odd_generic_genus419_scout.py
```

and compare both transcripts with

```text
05-knowledge/results/jc2_degree22_full_split_odd_generic_genus419_scout.out.
```

The companion reconstructs all `32` `Phi/Psi` rows independently by Laurent
recurrence and finite multinomial sums; verifies every weight and deck parity,
the Hilbert series and genus, top singular support, transverse `E_21` row,
six-line/three-branch quotient face, the all-even boundary gcds, and the
dominance minor.  Normal and optimized executions byte-match the stored
`18`-line output, use no Python `assert`, and have the declared LF hashes.

QED.
