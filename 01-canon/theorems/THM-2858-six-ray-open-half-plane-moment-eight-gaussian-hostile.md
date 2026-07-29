---
id: THM-2858
title: "Six-ray open-half-plane moment-eight Gaussian hostile"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT ACCEPTED,
  FINAL IMMUTABLE REPLAY PENDING.  Six adjacent-difference coefficients
  lying strictly in one open half-plane annihilate factorial moments one
  through four.  Equivalently the witness lies in the complex span of two
  independent strictly positive adjacent cones.  Its degree-thirteen
  two-charge Gaussian lift has moments one through eight zero, while
  moment ten is nonzero.  This refutes a finite moment-four sector
  detector but is not a GMC2 counterexample.
source: root/audit-2809-six-ray-sector-hostile-2026-07-28
depends_on: []
related:
  - THM-2841-all-order-adjacent-difference-factorial-tensor-positivity
  - THM-2846-arbitrary-positive-cone-moment-three-transverse-boundary
  - THM-2848-whitened-moving-plane-multipole-and-pearson-boundary
script: 04-computation/gmc_six_ray_half_plane_moment8_hostile_thm2858.py
output: 05-knowledge/results/gmc_six_ray_half_plane_moment8_hostile_thm2858.out
script_sha256: 47b2b08362c835f97e4af24b3df16e957819515d6522c4ed46f843971cbc5e97
output_sha256: e739c7316721c14672b9122bb5593ffce2bb20208cbda7eed0fa60c3c71f7e58
hash_basis: LF-normalized bytes
---

# THM-2858 -- six-ray open-half-plane moment-eight Gaussian hostile

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT ACCEPTED,
FINAL IMMUTABLE REPLAY PENDING.**

This is the next finite-moment obstruction beyond THM-2846.  It also
refutes the tempting strengthening suggested by the four-slot searches:
an open-half-plane coefficient sector, equivalently a binary positive-cone
plane, need not be detected by factorial moment four.

## 1. Statement

Put

```text
L(s^n)=n!,                 f_n=s^n/n!,
d_n=f_(n+1)-f_n.                                      (1)
```

There are coefficients `c_1,...,c_6` such that

```text
2000 Re(c_j)>|Im(c_j)|,          1<=j<=6,             (2)
```

and, for

```text
H=sum_(j=1)^6 c_j d_j,                                (3)
```

one has

```text
L(H)=L(H^2)=L(H^3)=L(H^4)=0,
L(H^5)!=0.                                             (4)
```

The construction is exact: the coefficients are the unique zero in a
displayed rational box, certified by a rational preconditioned-Newton
contraction.

## 2. The exact box

Normalize `c_1=1` and freeze the rational coordinates

```text
c_2=
 0.0006557567088377382
-1.2141776252864656 i,

c_4=
 0.009935469352174497
+1.2542668065909375 i.                                 (5)
```

Let the six real variables be the real and imaginary parts of
`c_3,c_5,c_6`.  Their rational center is

```text
c_3=
 0.000063621434994120405029506867116333343987061889208839
+0.040009898798925760784579140591476738268555395921997 i,

c_5=
 0.0045193541665913666254482593704175376727873260493792
-0.65481952791873310288768173269952383004730694141286 i,

c_6=
 0.01313872825314743624425793558960064291574202146893
+0.098948300256372361843344915268029197245827621460675 i.
                                                               (6)
```

Take the coordinatewise closed box `X` of radius

```text
r=10^(-45)                                               (7)
```

about this center.

For

```text
Phi=
 (Re L(H^2),Im L(H^2),
  Re L(H^3),Im L(H^3),
  Re L(H^4),Im L(H^4)),                                 (8)
```

the center Jacobian `J_0` is invertible.  Let `C=J_0^(-1)` and form the
preconditioned Newton map

```text
T(v)=v-C Phi(v).                                         (9)
```

Every entry is rational.  Exact natural-interval evaluation of the full
Jacobian on `X` proves

```text
T(X) subset int(X),
sup_(J in J(X)) ||I-CJ||_infinity<1.                    (10)
```

The minimum coordinate inclusion margin exceeds

```text
(99999/100000)r.                                        (11)
```

Thus `T` is a strict contraction of the complete box `X` into itself.
Its unique fixed point satisfies `C Phi=0`, and `C` is invertible, so it
is the unique zero of `(8)` in `X`.  This proves the last three
vanishings in `(4)`; the first is automatic from `L(d_j)=0`.

The same exact box calculation gives

```text
Re L(H^5)>647000000000,
Im L(H^5)<-3290000000000,                              (12)
```

so the fifth-moment exit is separated far from zero.

## 3. Two strictly positive cones

Inequality `(2)` is checked on the whole box, not just at its center.
Define

```text
a_j=(Re(c_j)-Im(c_j)/2000)/2,
b_j=(Re(c_j)+Im(c_j)/2000)/2.                          (13)
```

Then every `a_j,b_j` is strictly positive and

```text
c_j=a_j(1-2000i)+b_j(1+2000i).                         (14)
```

Consequently

```text
U=sum_j a_j d_j,             V=sum_j b_j d_j,
H=(1-2000i)U+(1+2000i)V.                               (15)
```

Both `U` and `V` are strictly positive adjacent-difference cones.  They
are independent: `c_1=1` gives `a_1=b_1`, whereas the nonzero imaginary
part of `c_2` gives `a_2!=b_2`.

The two rays in `(15)` have angular separation

```text
2 arctan(2000)<pi.                                     (16)
```

Thus the coefficient vector lies strictly in an open half-plane and in
the positive cone generated by two nonopposite rays.

## 4. Gaussian lift

Every `d_j` in `(3)` is divisible by `s`; write `H=sh`.  For a standard
centered complex Gaussian `Z`, put `W=conj(Z)` and

```text
P(Z,W)=W+Z h(ZW).                                      (17)
```

More explicitly,

```text
P=
 W-Z
 +(c_1-c_2) Z^2W/2!
 +(c_2-c_3) Z^3W^2/3!
 +(c_3-c_4) Z^4W^3/4!
 +(c_4-c_5) Z^5W^4/5!
 +(c_5-c_6) Z^6W^5/6!
 +c_6 Z^7W^6/7!.                                      (18)
```

Since `Re(c_6)>0`, this is a nonzero polynomial of exact total degree
thirteen.  Charge balance gives, for every `k>=1`,

```text
E[P^(2k)]=binom(2k,k)L(H^k),
E[P^(2k-1)]=0.                                         (19)
```

Equations `(4),(12),(19)` therefore prove

```text
E[P^m]=0,                 1<=m<=8,
E[P^10]=252L(H^5)!=0.                                  (20)
```

This is a finite-moment hostile, not a counterexample to the Gaussian
Moment Conjecture.

## 5. Mechanism and boundary

THM-2841 proves that every atomic tensor

```text
L(d_(i_1)...d_(i_k)),                 k>=2,
```

is strictly positive.  THM-2858 shows that even an open-half-plane
complex combination of those positive entries can cancel simultaneously
at orders two, three, and four.  The failure is not an atomic sign
failure; it is a transverse six-real-equation phase cancellation very
near, but strictly away from, the opposite-ray boundary.

THM-2846 supplied the corresponding moment-three/Gaussian-six hostile on
three adjacent rays and exited at moment four.  The present witness needs
six rays, survives one further factorial rung, and exits at moment five.
No minimality claim is made: the theorem does not prove that five rays
are impossible, nor does it address arbitrary shifted windows.

The Pearson conic of THM-2848 is compatible with this result.  It measures
strict real transverse departure around a bad plane; it does not exclude
the complex null points themselves.  The new hostile rules out using an
open-half-plane sector as the missing scalar-moment sidecar.

## 6. Exact evidence and audit boundary

The exact companion:

1. reconstructs every adjacent tensor through the direct `2^k` factorial
   expansion;
2. constructs the moment forms by symmetric multiindex expansion;
3. inverts the center Jacobian over `Q`;
4. interval-evaluates all thirty-six Jacobian entries and proves the
   strict inclusion and contraction bounds `(10)--(11)`;
5. checks the full-box sector inequalities `(2)`; and
6. proves the separated fifth-moment bounds `(12)`.

Every truth-bearing check uses an explicit exception and survives
optimized execution.  Normal, optimized, and stored transcripts agree.

An independent hostile audit rebuilt `H` directly as a polynomial,
computed `L(H^k)` from its factorial coefficients rather than importing
the companion's adjacent-tensor engine, reran the rational enclosure, and
checked the contraction inference, the two-cone decomposition, the
degree-thirteen lift, and every Gaussian charge factor.  It accepted the
mathematics.  Final immutable replay of the packaged paths and hashes is
still pending, so this file remains a proved candidate rather than a
proved dependency.
