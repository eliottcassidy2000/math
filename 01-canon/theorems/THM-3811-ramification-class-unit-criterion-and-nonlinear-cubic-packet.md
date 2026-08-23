---
id: THM-3811
title: "Ramification class unit criterion and nonlinear cubic packet"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT, with one
  explicitly OPEN decisive class gate.  Units on the maximal etale open of
  a normal finite completion are controlled exactly by integer relations
  among the deleted ramification classes.  The nonlinear Delone--Faddeev
  packet (a,b,c,d)=(A,C,7A,C^2-3A) is a normal nonmonogenic S3 cubic with a
  squarefree irreducible rational branch: the four Veronese rays glue into
  one ramification curve, with an explicit simple companion.  Whether that
  ramification class has infinite order, and whether its complement is
  affine with only constant units, remain OPEN.  No Jacobian counterexample
  is claimed.
source: jc_quartic_c3_construct / nonlinear binary-index design lane, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The divisor exact sequence, norm/companion
  torsion criterion, Delone--Faddeev signs, irreducibility and normality,
  squarefree discriminant, rational ramification normalization, four-branch
  vertex gluing, simple companion, and the bounded affine norm-Pell no-go
  were checked independently inside the proof.  The exact companion has a
  direct coefficient Groebner gate for the affine exponent-one norm cell.
  Independent hostile audit remains due; the infinite-order class and
  affineness gates are deliberately not promoted.
depends_on: []
related:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3808-homogeneous-linear-binary-cubic-veronese-unit-trap
script: 04-computation/jc2_nonlinear_binary_cubic_ramification_class_thm3811.py
output: 05-knowledge/results/jc2_nonlinear_binary_cubic_ramification_class_thm3811.out
script_sha256: 58a7f24e757ec8f5cbe9d39fbcc508662989585eb056e06f87a29e7528ac740b
output_sha256: 766083ef1042888ee4d86dfe5f15dff351037e1084f85e5760c9e980c91c700f
semantic_sha256: 583470aeadabb6be8a80609612756c969f44de29722a5115f746a3e2b098c963
hash_basis: raw LF bytes
---

# THM-3811 -- constant units are a ramification-class problem

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT; OPEN class
gate.**  Let `k` be algebraically closed of characteristic zero.

## 1. The exact deleted-divisor unit criterion

Let `Xbar=Spec S` be a normal integral affine surface and let
`E_1,...,E_r` be distinct prime divisors.  Put

```text
U=Xbar minus union_i E_i,                 B=Gamma(U,O_U).           (1)
```

Removing an additional codimension-two subset does not change `B`, because
`Xbar` is normal.  There is an exact sequence

```text
1 -> S* -> B* --div_E--> direct_sum_i Z E_i
  -> Cl(S) -> Cl(U) -> 0.                                        (2)
```

The middle map sends a unit on `U` to its Weil divisor on `Xbar`, which is
supported on the deleted primes.  Therefore, if `S*=k*`,

```text
B*=k*
iff the classes [E_1],...,[E_r] are Z-linearly independent in Cl(S). (3)
```

This is the exact obstruction behind THM-3808.  Its four ramification rays
lie in `Cl(S)=Z/3`, so relations are unavoidable.  A minimal escape glues
the four local rays into one global prime `E`; then `(3)` becomes

```text
B*=k* iff [E] has infinite order,             provided S*=k*.      (4)
```

Affineness of `U` is a separate requirement.  Sequence `(2)` computes its
global units whether or not `U` is affine; a counterexample factor surface
needs both `(3)` and affineness.

There is a useful norm form of the same gate.  Suppose `S` is finite cubic
over `R=k[A,C]`, one irreducible squarefree branch is `Gamma=V(delta)`, and
its generic pullback is

```text
pi^*Gamma=2E+F,                                                     (5)
```

with `E` ramified and `F` the simple companion.  Then `[E]` has finite
order if and only if there are `n>=1` and `u in S minus {0}` such that

```text
div_Xbar(u)=nE.                                                     (6)
```

Equivalently, `(6)` is the simultaneous norm/companion packet

```text
Norm_(K/k(A,C))(u)=kappa delta^n,        kappa in k*,
v_F(u)=0.                                                       (7)
```

Indeed a principal relation `n[E]=0` has an effective divisor `nE`, so its
principalizer has no poles and belongs to the normal ring `S`.  Pushing its
divisor to `A2` gives `(7)`.  Conversely, the norm equation together with
the companion valuation says that every height-one zero of `u` is on `E`,
and the total valuation is `n`.  Thus `(7)`, not a discriminant computation
alone, is the cheapest decisive constant-unit test.

## 2. The first nonlinear packet

Let `R=k[A,C]` and give the free module

```text
S=R*1 direct_sum R*omega direct_sum R*theta                         (8)
```

the Delone--Faddeev multiplication for

```text
(a,b,c,d)=(A,C,7A,C^2-3A):                                        (9)

omega^2=-7A^2+C omega-A theta,
omega theta=3A^2-AC^2,
theta^2=3AC-C^3+(C^2-3A)omega-7A theta.                           (10)
```

Its binary index form is

```text
I(x,y)=-(A x^3+C x^2y+7Axy^2+(C^2-3A)y^3).                       (11)
```

Thus every index value lies in `(A,C)` and no global generator exists.  The
only common zero of the four coefficients is the origin, where

```text
S/(A,C)S=k[omega,theta]/(omega,theta)^2.                            (12)
```

The trace discriminant is

```text
Delta = A(C+5A)(4C+19A)(3C-17A)
      + C^2(162A^3+126A^2C-4C^3)-27A^2C^4.                        (13)
```

Exact differentiation gives

```text
gcd(Delta,Delta_A,Delta_C)=1,                                     (14)
```

so the affine branch is reduced.  Unlike THM-3808, `(13)` is not a union of
radial lines.

The characteristic polynomial of `omega` is

```text
f(T)=T^3-CT^2+7A^2T+3A^3-A^2C^2.                                 (15)
```

It is irreducible over `k(A,C)`.  Because it is monic, a generic
factorization would specialize at `A=1`.  A root `v(C) in k(C)` of

```text
v^3-Cv^2+7v+3-C^2=0                                               (16)
```

cannot have a finite pole: its cubic term would have strictly largest pole
order.  Hence `v` is a polynomial.  Degree comparison leaves only
`v=C+beta`; substitution gives

```text
(beta-1)C^2+(2beta^2+7)C+(beta^3+7beta+3),                         (17)
```

which cannot vanish because `beta=1` makes the linear coefficient `9`.
Thus `(15)` is irreducible.  The free ring `(8)` is a domain, and the same
squarefree-discriminant DVR index argument as in THM-3808 gives `R1`; it is
free over `R`, hence `S2`.  Therefore `S` is normal.  Its nonsquare reduced
discriminant and simple inertia give generic Galois group `S3`.

## 3. Four rays glue into one rational ramification curve

On `A!=0`, put

```text
q=omega/A,                 G(q)=q^3+7q+3.                          (18)
```

Equation `(15)` becomes the rational chart

```text
A G(q)=C(C+q^2).                                                   (19)
```

At a double root, `f=f_T=0`.  Solving these two equations gives

```text
R(q)=(q-3)(q+1)(q+2)=q^3-7q-6,

A(q)=-2q^2R(q)/(3q^2+7)^2,
C(q)=-qR(q)/(3q^2+7),
omega(q)=A(q)q,
theta(q)=A(q)(q^2-7)/2.                                          (20)
```

Direct substitution verifies `(10)`, `(15)`, `f_T=0`, and
`Delta(A(q),C(q))=0`.  Conversely, at every point of `V(Delta)` with
`A!=0` and a nontriple cubic, the unique double root recovers `q`, and
`(20)` follows.  No branch component lies in `A=0`, because

```text
Delta(0,C)=-4C^5.                                                  (21)
```

Hence the closure of `(20)` is the whole branch curve.  It is irreducible
and rational; squarefreeness in `(14)` then makes `Delta` its irreducible
equation.

Most importantly,

```text
q=-2,-1,0,3                                                     (22)
```

are four distinct normalization points mapping to

```text
(A,C,omega,theta)=(0,0,0,0).                                     (23)
```

The degree-four initial form of `(13)` is exactly the four-line
discriminant of THM-3808.  Thus the nonlinear term `C^2` performs the
minimal geometric repair: it glues the four radial ramification rays into
one global curve while retaining them as four local branches at the
square-zero vertex.

The simple companion is equally explicit.  On the branch, the roots of
`(15)` are

```text
T=Aq  (double),                 T_comp=C-2Aq  (simple).             (24)
```

The corresponding `theta` is recovered from `(10)`.  Therefore the generic
pullback is exactly `(5)`, with one rational ramification curve and one
residue-degree-one companion.  The packet passes the nonmonogenic,
squarefree, `S3`, and companion gates of THM-3801 without the Veronese
arrangement symmetry.

## 4. First exact norm-Pell obstruction

Write an arbitrary element uniquely as

```text
u=p+q_1 omega+q_2 theta.                                           (25)
```

For the universal Delone--Faddeev law, its norm is

```text
a^2d q_1^3+2abd q_1^2q_2-ac^2 q_1^2q_2-2acd q_1q_2^2
+acp q_1^2-ad^2q_2^3+3adp q_1q_2+b^2d q_1q_2^2
-bcp q_1q_2+bdp q_2^2+bp^2q_1-cp^2q_2+p^3.                       (26)
```

The exact companion substitutes `(9)`, lets each of `p,q_1,q_2` be an
arbitrary affine polynomial in `(A,C)`, and equates `(26)` to
`kappa Delta` with `kappa!=0`.  The resulting 38 coefficient equations in
ten coefficient variables, plus the inverse for `kappa`, have reduced
Groebner basis

```text
[1].                                                               (27)
```

Thus no exponent-one principalizer in `(7)` has affine coordinate profiles.
This is a decomposability obstruction, not merely a linear span test.

There is also an all-exponent degree floor.  If
`deg p,deg q_1,deg q_2<=D`, then `(26)` has total degree at most

```text
3D+5,                                                              (28)
```

whereas `Delta^n` has degree `6n`.  Therefore an order-`n` torsion witness
must have

```text
D>=2n-1;                                                           (29)
```

for `n>=2`, while `(27)` strengthens the `n=1` floor to `D>=2`.

The OPEN endpoint is now exact:

```text
OPEN-3811:
  prove [E] has infinite order in Cl(S), prove S*=k*, and prove
  U=Xbar minus E is affine; or exhibit the first norm/companion witness
  (7) and explain the resulting nonconstant unit.                  (30)
```

Until all three positive clauses in `(30)` are proved, this is not an
etale constant-unit factor surface and not a planar Jacobian counterexample.
The proved gain is the first nonlinear normal `S3` completion that glues the
four fatal Veronese rays while preserving the exact companion grammar, plus
the decisive class/norm test that remains.  **QED for `(2)--(29)`;
`OPEN-3811` remains open.**
