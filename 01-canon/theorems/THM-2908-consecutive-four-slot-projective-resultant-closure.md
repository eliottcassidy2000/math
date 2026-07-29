---
id: THM-2908
title: "Consecutive four-slot projective resultant closure"
status: >
  PROVED CANDIDATE; EXACT FACTOR REPLAY AND FINAL INDEPENDENT AUDIT
  PENDING.  For every n>=0, the first four factorial moments have no
  common nonzero projective zero on the consecutive support
  {n,n+1,n+2,n+3}.  A sign-blind degree-2804 direct resultant closes the
  entire finite moving-plane chart, including selector-zero fibres, and a
  separate degree-49 resultant closes the projective line at infinity.
  Thus consecutive first-window SFC(4) holds at every depth.
source: root/consecutive-projective-resultant-closure-2026-07-29
depends_on:
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
  - THM-2872-four-slot-shared-multipole-quartic-norm-and-response-secant-reduction
  - THM-2879-all-shift-cubic-null-endpoint-holonomy-exit
related:
  - THM-2836-sfc3-arbitrary-support-shifted-window-census
  - THM-2848-whitened-moving-plane-multipole-and-pearson-boundary
  - THM-2866-positive-factorial-difference-semiring-and-cubic-pascal-response-ladder
  - THM-2890-discrete-newton-closure-of-strict-consecutive-gmc-wedges
  - THM-2891-zero-normal-coordinate-consecutive-cone-boundary-cubic-emptiness
  - THM-2906-atomic-tp3-does-not-orient-mixed-endpoint-holonomy
script: 04-computation/gmc_consecutive_four_slot_projective_closure_thm2908.py
output: 05-knowledge/results/gmc_consecutive_four_slot_projective_closure_thm2908.out
---

# THM-2908 -- consecutive four-slot projective resultant closure

**PROVED CANDIDATE; EXACT FACTOR REPLAY AND FINAL INDEPENDENT AUDIT
PENDING.**

Put

```text
L(s^m)=m!,                 f_j=s^j/j!,
d_j=f_(j+1)-f_j.                                      (1)
```

For every integer `n>=0` and every nonzero complex coefficient vector
`z=(z_0,z_1,z_2,z_3)`, define

```text
H_z=sum_(j=0)^3 z_j f_(n+j).                         (2)
```

Then

```text
L(H_z), L(H_z^2), L(H_z^3), L(H_z^4)
```

do not vanish simultaneously.  Equivalently, the Strong Factorial
Conjecture holds on every consecutive four-slot support

```text
{n,n+1,n+2,n+3}
```

in its first window.

The proof is projective and sign-blind.  It excludes every real
moving plane in the mean-zero three-space, not only a positive,
cone-cutting, or cone-avoiding sector.

## 1. From a common point to a real moving plane

The real mean-zero space of `(2)` is

```text
W_n=span_R{d_n,d_(n+1),d_(n+2)}.                     (3)
```

On `W_n` put

```text
Q(H)=L(H^2),       C(H)=L(H^3),       K(H)=L(H^4).   (4)
```

The factorial integral makes `Q` positive definite over `R`.
THM-2843 proves the exact moving-plane equivalence:

```text
Q=C=K has a nonzero complex projective point
```

if and only if there is a real two-plane `E subset W_n` such that

```text
Q|_E divides C|_E,              Q|_E divides K|_E.    (5)
```

Indeed, a common point `u+iv` has independent real and imaginary parts,
and `E=span_R(u,v)` is the required plane; conversely either complex
root of `Q|_E` is common to all three forms.

For context, THM-2848 gives the equivalent Maxwell/multipole description,
and THM-2866 removes its global vanishing-quartic-harmonic branch.  Those
results are not needed for the direct reduction here: THM-2843 alone
shows that it is enough to prove that no real plane satisfies `(5)`.

## 2. The complete projective atlas

Write

```text
e_0=d_n,              e_1=d_(n+1),          e_2=d_(n+2). (6)
```

A real plane in `(3)` has a projective normal `[a:b:c] in RP^2`.
When `c!=0`, it has the unique finite high-chart form

```text
E_(X,Y)=span_R{U,V},
U=e_0+X e_2,                    V=e_1+Y e_2,          (7)

X=-a/c,                         Y=-b/c.
```

The line `c=0` has a second affine chart.  When `b!=0`, write

```text
E_t=span_R{e_0-t e_1,e_2},       t=a/b.              (8)
```

Its remaining point `b=0` is

```text
E_infinity=span_R{e_1,e_2}.                            (9)
```

Equations `(7)--(9)` are the whole real projective plane of moving
planes.  In particular, no sign assumption on `X,Y,t` is available or
needed.

## 3. Cubic divisibility and the finite selector

For a basis `U,V`, write

```text
g_0=L(U^2),       g_1=L(UV),        g_2=L(V^2),
t_0=L(U^3),       t_1=L(U^2V),
t_2=L(UV^2),      t_3=L(V^3).                       (10)
```

As in THM-2879, `Q|_E` divides `C|_E` exactly when

```text
I_1=3t_1g_0g_2-t_3g_0^2-2t_0g_1g_2=0,
I_2=3t_2g_0g_2-2t_3g_1g_0-t_0g_2^2=0.              (11)
```

Substitute the finite chart `(7)`, clear the positive common tensor
denominators

```text
72(2n+1)^2(2n+3)^2(3n+1)(3n+2)(3n+4)(3n+5),
36(2n+1)(2n+3)^2(3n+1)(3n+2)(3n+4)(3n+5),           (12)
```

and call the resulting polynomials `P_1(n,X,Y)` and `P_2(n,X,Y)`.
Both denominators are positive for `n>=0`.  Their subresultant degrees
in `X` are

```text
4,3,2,1,0.                                           (13)
```

The terminal resultant is, up to a nonzero rational scalar,

```text
Res_X(P_1,P_2)=C_+(n) G(n,Y)^2 P(n,Y),               (14)
```

where `C_+(n)` has no zero for `n>=0`,

```text
G(n,Y)
 =(4n+6)Y^2+(4n+6)Y+n+2
 =(4n+6)(Y+1/2)^2+1/2>0,                             (15)
```

and `P` has bidegree `(21,15)` and `352` terms.  Thus every real
cubic-divisible finite-chart plane satisfies

```text
P(n,Y)=0.                                             (16)
```

Here the full exact cofactor, apart from a nonzero rational scalar, is

```text
C_+(n)
 =(n+1)^2(n+2)^5(2n+1)^4(3n+1)^3(3n+2)^2 C_10(n),
```

and all eleven coefficients of `C_10` are strictly positive.  This
includes the depth `n=0`; it is not inherited from a positive-sector
argument.

The primitive linear subresultant is

```text
A(n,Y)X-N(n,Y),                                      (17)
```

with profiles

```text
A: bidegree (22,10), 251 terms;
N: bidegree (22,10), 253 terms.                       (18)
```

The removed common content is, up to a nonzero rational scalar,

```text
864(n+1)(n+2)^3(2n+1)^4(3n+1)^3(3n+2)^2 G(n,Y).     (19)
```

It is nonzero for real `Y` and `n>=0`.  Therefore every actual real
common zero of `(11)` satisfies `(17)`, including any fibre on which
`A=0`.

## 4. Homogeneous endpoint clearing

Let `J(n,X,Y)` be the denominator-cleared endpoint-holonomy statistic
of THM-2872/2879.  It has degree five in `X`.  Write

```text
J=sum_(j=0)^5 e_j(n,Y)X^j
```

and define the homogeneous selector clearing

```text
F(n,Y)
 =sum_(j=0)^5 e_j(n,Y)N(n,Y)^j A(n,Y)^(5-j).         (20)
```

Its exact profile is

```text
deg_n F=121,             deg_Y F=55,
terms(F)=6820.                                           (21)
```

The normalized endpoint statistic has the positive denominator

```text
1024(n+3)(2n+1)^2(2n+3)^4
    (4n+1)(4n+3)(4n+5)(4n+7).                         (22)
```

Quartic divisibility in `(5)` requires endpoint agreement, hence

```text
J(n,X,Y)=0.                                            (23)
```

If `A!=0`, equations `(17),(20),(23)` give

```text
F=A^5 J(n,N/A,Y)=0.                                    (24)
```

If `A=0`, the actual linear subresultant `(17)` forces `N=0`, and every
term of the homogeneous expression `(20)` again vanishes.  Thus the
selector-zero fibre is not divided away:

```text
every finite shared cubic--quartic plane gives P=F=0. (25)
```

This is the reason for using `(20)` rather than a rational substitution.

## 5. The direct affine resultant

Exact elimination over `Z[n]` gives

```text
R(n)=Res_Y(P,F),

deg R=2804,                  terms(R)=2805.            (26)
```

Its complete nonlinear degree/multiplicity profile is

```text
203^1, 65^10, 19^3, 10^50, 5^161,                    (27)
```

and its linear part is

```text
(n+3)^15 (2n+1)^24 (3n+1)^37
(n+1)^48 (3n+2)^105 (n+2)^360.                       (28)
```

The degree `5`, `10`, and `19` factors have positive value at `n=0`
and strictly positive Gregory--Newton differences at base `1`.
For any such degree-`d` factor `f`, the exact identity

```text
f(n)=sum_(k=0)^d Delta^k f(1) binom(n-1,k),           n>=1,
```

therefore makes `f(n)>0`; the separately checked `f(0)>0` includes
depth zero.  Consequently these factors are positive at every integer
`n>=0`.  The degree
`65` and `203` factors are sign-indefinite in that basis, but the exact
companion proves respectively that they have no residue root modulo

```text
43 and 83.                                             (29)
```

An integer root would reduce to a residue root, so neither mixed factor
has an integer zero.

Every factor in `(26)` is therefore nonzero for integer `n>=0`:

```text
R(n)!=0.                                               (30)
```

Equations `(25),(30)` exclude every finite-chart plane satisfying
`(5)`.  Notice that this conclusion is stronger than any sectorwise
sign result: the resultant is used only for nonvanishing.

## 6. The line at infinity

For the finite points `(8)`, direct substitution into `(11)` gives
polynomials of profiles

```text
I_1: deg_t=4, deg_n=8, 45 terms;
I_2: deg_t=3, deg_n=7, 32 terms.                      (31)
```

Their exact resultant is

```text
Res_t(I_1,I_2)
 =-62208
   (n+1)^3(n+2)^10(2n+1)^4(2n+3)^2
   (3n+1)^3(3n+2)^3 P_5(n)P_19(n),                  (32)
```

where

```text
P_5(n)
 =16n^5+344n^4+865n^3+731n^2+247n+29,               (33)
```

and every coefficient of the degree-`19` factor `P_19` is strictly
positive.  Hence `(29)` is nonzero for every `n>=0`, and no finite `t`
plane is even cubic-divisible.

At the missing point `(9)`, the two cleared cubic invariants are

```text
-(n+2)^2(28n^2+87n+66),
-2(n+2)(4n+5).                                       (34)
```

Both are nonzero for `n>=0`.  This closes `t=infinity` and therefore
the entire projective boundary.

## 7. Consequence and scope

Sections 5--6 exclude `(5)` on every real plane in `(3)`.  The
THM-2843 moving-plane equivalence now proves the claim following `(2)`.
As an explanatory translation into THM-2848's language, THM-2866
removes the global
`F^o=0` branch and this theorem removes every shared
cubic--quartic multipole line.

This is an infinite exact SFC(4) family, but its scope is sharp:

```text
support:  four consecutive normalized factorial monomials;
window:   the first four moments only;
depth:    every integer n>=0.                          (35)
```

No arbitrary four-slot support, shifted moment window, general SFC(4),
full Strong Factorial Conjecture, GMC(2), Gaussian nullcone, or planar
Jacobian-conjecture conclusion is asserted.  THM-2906 remains the sharp
warning that ambient tensor positivity or separate TP3 cannot replace
the cubic selector used here.

## 8. Exact verification

The exact companion hash-pins THM-2879, reconstructs `(12)--(22)`,
computes and refactors `(26)` over `Z[n]`, verifies every
Gregory--Newton and root-free-prime certificate, reconstructs
`(31)--(34)`, and uses explicit `require` gates under optimized Python.
It requires the exactly pinned `python-flint==0.9.0` for the degree-2804
resultant.  Install that isolated computation dependency with

```text
python -m pip install -r 04-computation/requirements-gmc-projective-resultant.txt
```

Then run

```text
python 04-computation/gmc_consecutive_four_slot_projective_closure_thm2908.py
python -O 04-computation/gmc_consecutive_four_slot_projective_closure_thm2908.py
```

Immutable hashes will be inserted after normal, optimized and stored
replay.

**QED (candidate pending exact replay and final independent audit).**
