---
id: THM-2906
title: "Atomic tensor positivity and local TP3 do not orient mixed endpoint holonomy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the positive family
  U=d_1+x d_3, V=d_2+y d_3, all adjacent tensors are positive and both
  local response curves have strict TP3, but the mixed endpoint-holonomy
  determinant is positive at (x,y)=(1,1) and negative at (1,2).
  Both cells lie outside the simultaneous cubic-null ideal.  Thus the
  cubic selector or an equivalent mixed-ideal sidecar is load-bearing;
  neither atomic positivity nor separate response curvature determines
  the endpoint sign.
source: root/mixed-endpoint-holonomy-no-go-2026-07-29
audit: >
  Two independent hostile audits reconstructed the tensors directly from
  the polynomials d_j, verified both quotient-chart derivatives, witness
  cells, cubic invariants, TP3 values, polynomial census and hashes, and
  accepted the cubic-ideal scope without qualification.
depends_on:
  - THM-2853-gamma-adjacent-tensor-cycle-weighted-positivity
  - THM-2873-two-ray-factorial-response-tp3-curvature
related:
  - THM-2872-four-slot-shared-multipole-quartic-norm-and-response-secant-reduction
  - THM-2879-all-shift-cubic-null-endpoint-holonomy-exit
script: 04-computation/gmc_atomic_tp3_mixed_holonomy_no_go_thm2906.py
output: 05-knowledge/results/gmc_atomic_tp3_mixed_holonomy_no_go_thm2906.out
script_sha256: 10a0e72b084ff54f6eb81b65690bd3cd04b1df7af92f713973e1008dbcf60732
output_sha256: 9e84dd2ac5b86fdc48365b005dbb9e8ae87bce3216bf940b0b1ecef8cd2f8d04
hash_basis: LF-normalized bytes
---

# THM-2906 -- atomic TP3 does not orient mixed endpoint holonomy

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Put

```text
L(s^m)=m!,                 f_j=s^j/j!,
d_j=f_(j+1)-f_j.                                      (1)
```

For `x,y>0`, define

```text
U_x=d_1+x d_3,                  V_y=d_2+y d_3.         (2)
```

THM-2853 gives positivity of every tensor containing at least two
nonzero positive adjacent-cone elements.  THM-2873 gives strict local
TP3 for both response curves in `(2)`: `U_x` is a gap-two cone and `V_y`
is a gap-one cone.

Nevertheless the mixed endpoint-holonomy determinant has both signs
inside `(2)`.  This proves that those two positivity statements, even
together, do not orient the mixed endpoint comparison of THM-2872.

## 1. The mixed endpoint statistic

Write

```text
g_0=L(U^2),              g_1=L(UV),             g_2=L(V^2),
A_i=L(U^(4-i)V^i),                              0<=i<=4. (3)
```

The left and right determinations of the linear coefficient in a possible
quartic quotient are

```text
r_L=(2A_1g_0-A_0g_1)/g_0^2,
r_R=(2A_3g_2-A_4g_1)/g_2^2.                          (4)
```

Their denominator-cleared difference is

```text
J=g_0^2g_2^2(r_L-r_R).                               (5)
```

For the symbolic family `(2)`, exact expansion gives

```text
deg_x J=deg_y J=5,
number of monomials=30,
positive coefficients=13,
negative coefficients=17.                            (6)
```

Four low coefficients already expose the alternating boundary:

```text
[1]J=273888,              [x]J=-5839200,
[y]J=5559264,             [xy]J=-82033920.            (7)
```

Thus coefficientwise positivity of the tensors entering `(3)` is lost
in the endpoint subtraction itself.

## 2. Two exact positive cells with opposite holonomy

At `(x,y)=(1,1)`,

```text
(g_0,g_1,g_2)=(30,37,46),
g_0g_2-g_1^2=11,

r_L=73097746/75,          r_R=515413080/529,
J=610878432>0.                                        (8)
```

At `(x,y)=(1,2)`,

```text
(g_0,g_1,g_2)=(30,61,126),
g_0g_2-g_1^2=59,

r_L=142883488/75,         r_R=120168140/63,
J=-33115086144<0.                                     (9)
```

Both quadratic Gram determinants are positive, as they must be for
independent real polynomials.

For a positive cone `W`, let

```text
Delta_1(W)=det [
 L(d_1W)   L(d_1W^2)   L(d_1W^3)
 L(d_2W)   L(d_2W^2)   L(d_2W^3)
 L(d_3W)   L(d_3W^2)   L(d_3W^3)
].                                                    (10)
```

The exact local curvatures in the two cells are

```text
Delta_1(U_1)=36467293512,

Delta_1(V_1)=86548240800,
Delta_1(V_2)=2255452579800.                           (11)
```

They are all strictly positive.  Equations `(8)--(11)` are therefore
literal witnesses to

```text
atomic tensor positivity + separate local TP3
does not determine sign(J).                           (12)
```

This is stronger than observing that the proof of THM-2873 does not
mention `J`: the desired implication is false on the same positive
two-ray support architecture.

## 3. Why the statistic is genuinely cross-chart

Put

```text
q(z)=g_0+2g_1z+g_2z^2,
f(z)=A_0+4A_1z+6A_2z^2+4A_3z^3+A_4z^4.              (13)
```

In the left chart let `Phi=f/q`.  In the reciprocal right chart let

```text
q^vee(w)=w^2q(1/w),             f^vee(w)=w^4f(1/w),
Psi=f^vee/q^vee.                                      (14)
```

Direct differentiation gives

```text
Phi'(0)=2r_L,                  Psi'(0)=2r_R.           (15)
```

Consequently `J` compares derivatives of two separately normalized
quotient charts.  Local TP3 orients consecutive secants *within* each
fixed-`W` response curve; it supplies no comparison between the two
normalizations in `(15)`.

## 4. The load-bearing cubic ideal

With

```text
t_0=L(U^3),     t_1=L(U^2V),     t_2=L(UV^2),     t_3=L(V^3),
```

the two quadratic-divisibility remainders are

```text
I_1=3t_1g_0g_2-t_3g_0^2-2t_0g_1g_2,
I_2=3t_2g_0g_2-2t_3g_1g_0-t_0g_2^2.                 (16)
```

Their exact values are

```text
(I_1,I_2)_(1,1)=(141000,169320),
(I_1,I_2)_(1,2)=(957960,2436120).                    (17)
```

Thus neither hostile cell is cubic-null.  The sign reversal does not
contradict THM-2879.  Rather, it proves that THM-2879's use of the exact
cubic selector is essential: on the unique positive common zero of
`I_1,I_2`, its coefficientwise-negative pseudo-remainder gives `J<0`,
whereas the ambient positive chamber contains both signs.

The repaired next general target is therefore mixed and ideal-relative.
One sufficient form would be a division-free identity, on an explicitly
typed positive chamber,

```text
-D A^m J=M_1I_1+M_2I_2+C,                            (18)
```

where `D,A` are positive and `C` has a positive certificate.  Equivalently
one may retain an exact selector and prove that the pseudo-remainder of
`J` modulo the cubic branch eliminant has one sign.  Identity `(18)` is
an **OPEN design target**, not a consequence of this theorem.

## 5. Scope and exact verification

This theorem is a stopping boundary.  It proves no universal TP3
statement, no sign claim on an arbitrary cubic-null branch, no new
four-slot closure, and no SFC, GMC, or Gaussian-nullcone theorem.  It
does prove that extending atomic positivity or within-curve TP3 alone is
the wrong logical target for the remaining general mixed four-slot line.

The exact companion:

1. reconstructs adjacent tensors by literal `2^k` factorial expansion;
2. checks `31` positive atomic cells of orders two through four;
3. constructs the full endpoint polynomial and its sign census;
4. verifies both quotient-derivative identities `(15)`;
5. checks every rational value in `(8)--(11)` and `(17)`; and
6. retains explicit `require` gates under optimized Python.

Run

```text
python3 04-computation/gmc_atomic_tp3_mixed_holonomy_no_go_thm2906.py
python3 -O 04-computation/gmc_atomic_tp3_mixed_holonomy_no_go_thm2906.py
```

Both modes byte-match

```text
05-knowledge/results/gmc_atomic_tp3_mixed_holonomy_no_go_thm2906.out.
```

Two independent hostile audits reconstructed the moments directly from the
polynomials `d_j`, rather than the companion's tensor routine.  They matched
both witness cells, derivative identities, cubic invariants, TP3 values,
the full endpoint-polynomial census and digest, and the declared LF hashes.
Both accepted the logical scope exactly as stated.

**QED.**
