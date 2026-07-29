---
id: THM-2946
title: "Full Macaulay maximal-minor gcd and chart-free resultant"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  Over the universal coefficient ring, the primitive
  polynomial gcd of all 36-by-36 maximal minors of the full degree-seven
  Macaulay map for ternary forms of degrees (2,3,4) is exactly the
  irreducible resultant, with exponent one.  The radical statement follows
  from complete-intersection regularity and socle degree six; the exponent
  is fixed by THM-2942's explicit minor
  q200^6*c300*K*Res.  This is a universal chart-free statement.  After
  specialization to a one-parameter factorial family, new common factors
  may appear and must still be audited.
source: codex-gmc-uniform-width-extension-2026-07-29
depends_on:
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
related:
  - THM-2927-general-width-flagged-macaulay-leading-coefficient
  - THM-2945-nonnegative-complete-intersection-norm-and-repeated-divisor-gate
script: 04-computation/gmc_full_macaulay_maximal_minor_gcd_thm2946.py
output: 05-knowledge/results/gmc_full_macaulay_maximal_minor_gcd_thm2946.out
script_sha256: da3486a0e982b26dc7b11b04ffe7370b24ad0ec13f12cad089815901eff9bbc2
output_sha256: 29c09f7c801eb17c17fb25a8c43ceaad2ac0dc1beda4fc84d8d0918620a7836f
hash_basis: LF-normalized bytes
---

# THM-2946 -- full Macaulay maximal-minor gcd

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

Let

```text
S=Z[x0,x1,x2],
A=Z[all coefficients of Q,C,F],                      (1)
```

where `Q,C,F` are universal ternary forms of degrees `2,3,4`.
Consider the full degree-seven multiplication map

```text
Phi_7:
 (S_5 tensor A) direct-sum (S_4 tensor A)
                    direct-sum (S_3 tensor A)
       -> S_7 tensor A,
(u,v,w) |-> uQ+vC+wF.                                (2)
```

Its matrix has size

```text
36 by (21+15+10)=36 by 46.                            (3)
```

Let `I_36(Phi_7)` be the ideal of all its maximal minors.  Then the
primitive polynomial gcd of the nonzero members of this ideal is

```text
gcd_primitive{all 36-by-36 minors}
   =Res(Q,C,F),                                       (4)
```

up to the fixed sign convention for the ternary resultant.

Equation `(4)` is a gcd statement, not the stronger and generally false
claim

```text
I_36(Phi_7)=(Res).                                    (5)
```

Individual minors retain chart/flag cofactors.  Taking their universal
gcd removes those cofactors and leaves the chart-free common-zero
divisor.

## 1. Rank loss is exactly the resultant hypersurface

Work first over an algebraically closed field `k` of characteristic
zero.  If

```text
Res(Q,C,F)=0,                                         (6)
```

the three forms have a common projective zero `P`.  Evaluation at `P`
is a nonzero functional on `S_7`, and it annihilates every element in
the image of `(2)`.  Hence

```text
rank Phi_7<36.                                        (7)
```

Conversely suppose the resultant is nonzero.  Then `(Q,C,F)` has no
common projective zero.  Its homogeneous ideal has height three in the
three-variable Cohen--Macaulay ring `k[x0,x1,x2]`; therefore `Q,C,F`
form a regular sequence.  The complete-intersection Hilbert series is

```text
(1-t^2)(1-t^3)(1-t^4)/(1-t)^3
 =1+3t+5t^2+6t^3+5t^4+3t^5+t^6.                    (8)
```

In particular the quotient has socle degree

```text
(2-1)+(3-1)+(4-1)=6                                  (9)
```

and its degree-seven piece vanishes.  Thus `(2)` is surjective and

```text
rank Phi_7=36.                                        (10)
```

Equations `(7)--(10)` prove the set-theoretic identity

```text
V(I_36(Phi_7))=V(Res).                               (11)
```

The exact dimension invoice is consistent with the Koszul complex:

```text
dim domain=46,
dim target=36,
dim ker Phi_7
 =dim S_2+dim S_1+dim S_0
 =6+3+1=10.                                          (12)
```

There is no degree-seven contribution from the triple Koszul syzygy,
because `7-(2+3+4)<0`.

## 2. The only possible common irreducible divisor

The universal resultant is primitive and irreducible.  Every maximal
minor vanishes on its hypersurface by Section 1, so the resultant
divides every maximal minor over the universal UFD `A`.

Now let `H` be any nonconstant irreducible polynomial dividing all
maximal minors.  Then

```text
V(H) subset V(I_36)=V(Res).                           (13)
```

Both sides of `(13)` contain irreducible hypersurfaces.  The
Nullstellensatz gives `Res in radical(H)=(H)`, hence `H` is associated
to `Res`.  Therefore for some integer `e>=1`,

```text
gcd_primitive{maximal minors}=Res^e.                  (14)
```

This is the Fitting-support step: the zeroth Fitting ideal of
`coker(Phi_7)` may have higher-codimension structure, but it has only
the resultant as a divisorial component.

## 3. The exponent is exactly one

THM-2942 computes one literal maximal minor in the same full matrix:

```text
Delta_J0
 =q200^6*c300*K(Q,C)*Res(Q,C,F),                      (15)
```

where

```text
K
 =c120*q200^2-c210*q110*q200
  -c300*q020*q200+c300*q110^2.                       (16)
```

The cofactor in `(15)` is nonzero and independent of every coefficient
of `F`.  The resultant has degree six in the coefficients of `F`, so it
cannot divide that cofactor.  Consequently

```text
valuation_Res(Delta_J0)=1.                            (17)
```

Equations `(14)` and `(17)` force `e=1`, proving `(4)`.

This argument is stronger than a generic-rank heuristic: it proves the
universal polynomial multiplicity exactly and retains the selected flag
factor rather than silently setting it to one.

## 4. One simple common-root control

The exponent-one mechanism has a small independent specialization.
Put

```text
Q=x^2-z^2,
C=y^3-z^3,
F_t=z^3(x+2y-3z)+t z^4.                              (18)
```

The complete intersection `Q=C=0` consists of the six reduced points

```text
(x,y,z)=(+-1,omega,1),             omega^3=1.         (19)
```

At `t=0`, only `(1,1,1)` lies on `F_0=0`, and the vanishing is
transverse.  Poisson multiplication over the six points gives, up to
one nonzero convention constant (written `doteq`),

```text
Res(Q,C,F_t)
 doteq [(t-2)^3+8][(t-4)^3+8]
 =t(t^2-6t+12)(t^3-12t^2+48t-56).                   (20)
```

Hence

```text
ord_(t=0) Res=1,
derivative at zero of the displayed norm=-672!=0.    (21)
```

The full matrix in `(2)` has rank `35` at `t=0` and rank `36` at
`t=1`.  Thus the specialization exhibits literally the generic
corank-one transverse crossing used abstractly above.

## 5. Specialization boundary

Equation `(4)` lives in the universal coefficient ring.  Polynomial gcd
does not commute with specialization.  The elementary pair

```text
gcd(u,v)=1 in Z[u,v],
but after v=u: gcd(u,u)=u                              (22)
```

is the sharp model.

Therefore, after inserting a one-parameter factorial family

```text
Q=Q_n,             C=C_n,             F=F_n,          (23)
```

the specialized maximal minors may acquire an additional common factor.
The theorem proves that such a factor is not a second universal
Macaulay/Pluecker divisor.  It does **not** prove that the specialized
resultant is nonzero, squarefree, or root-free on `n>=0`.

In particular `(4)` does not by itself close arbitrary-width SFC(4), the
remaining GMC four-slot branch, the Gaussian Moment Conjecture, or any
Jacobian-conjecture chart.  Those applications still require a
specialized resultant or repeated-divisor audit such as the gates exposed
in THM-2945.

## 6. Exact verification

The exact companion:

1. verifies all graded dimensions in `(3)`, `(8)` and `(12)`;
2. constructs the full `36`-by-`46` Macaulay matrix for `(18)` and checks
   ranks `35` and `36` at `t=0,1`;
3. derives both cubic norm factors in `(20)` by exact univariate
   resultants and checks `(21)`;
4. checks that all five nonselected complete-intersection points have
   nonzero `F_0` product; and
5. verifies the specialization boundary `(22)`.

All truth-bearing gates use explicit exceptions, so optimized execution
retains the full audit.

Run

```text
python 04-computation/gmc_full_macaulay_maximal_minor_gcd_thm2946.py
python -O 04-computation/gmc_full_macaulay_maximal_minor_gcd_thm2946.py
```

Both modes must LF-normalize to the stored transcript and the hashes in
the frontmatter.

Promotion requires an independent audit of the Fitting-support direction,
the exponent-one use of `(15)`, and the universal-versus-specialized
scope.
