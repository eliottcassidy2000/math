---
id: THM-2947
title: "Conjugate-pair corank parity and one-minor resultant gate"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  For a real ternary (2,3) complete intersection with
  no real geometric support, the kernel of quartic multiplication on
  its length-six Artin algebra has even real dimension, including in
  the nonreduced case.  Hence the full degree-seven Macaulay map has
  rank in {30,32,34,36}; one nonzero 35-minor forces rank 36 and
  nonvanishing of the genuine (2,3,4) resultant.  Equivalently, the
  sum of squares of all 35-minors is a positive exact resultant gate
  on this real locus.  On factorial first-window forms, any one of the
  216 quartic-row cofactors in the inherited chart is a sufficient
  degree-at-most 55M-35 certificate.  No fixed cofactor is proved
  uniformly nonzero and no new support width is closed.
source: codex-gmc-uniform-width-extension-2026-07-29
audit: Pending independent hostile audit.
depends_on:
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
related:
  - THM-2945-nonnegative-complete-intersection-norm-and-repeated-divisor-gate
  - THM-2946-full-macaulay-maximal-minor-gcd-and-chart-free-resultant
script: 04-computation/gmc_conjugate_pair_corank_parity_thm2947.py
output: 05-knowledge/results/gmc_conjugate_pair_corank_parity_thm2947.out
script_sha256: d1bd09ff20925183f5488fcd8850469867f1dfad2bdb808504fc896708605744
output_sha256: 7369dce44b44d45dc435e8fc26ad119d96179689efde07780c66ad0dbb9c7eea
hash_basis: LF-normalized bytes
---

# THM-2947 -- conjugate-pair corank parity and one-minor resultant gate

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. The length-six algebra and its real parity

Let

```text
S=R[x0,x1,x2],                                      (1)
```

and let `Q,C,F` be homogeneous forms of degrees `2,3,4`.  Assume

```text
Q,C are a regular sequence,
V_R(Q,C)=empty.                                     (2)
```

Choose a real linear form `L` missing the finite complex scheme
`V(Q,C)`, and put

```text
A=R[x0,x1,x2]/(Q,C,L-1),             dim_R A=6,
f=F/L^4 in A.                                        (3)
```

Then

```text
dim_R ker(m_f:A->A) is even.                         (4)
```

Indeed, after complexification the Artin algebra splits as the product
of its local factors.  Condition `(2)` says that complex conjugation
has no fixed maximal ideal, so the local factors occur in distinct
pairs

```text
(A_P,A_Pbar).                                        (5)
```

Because `f` is real, the two multiplication matrices in each pair are
complex conjugates and their kernels have the same complex dimension.
Moreover,

```text
(ker_R m_f) tensor_R C = ker_C(m_f tensor C).         (6)
```

Each pair therefore contributes twice one complex kernel dimension to
the real kernel.  This proof does not assume that the six geometric
points are reduced.

The four possible rank values are consequently

```text
rank_R(m_f) in {6,4,2,0}.                             (7)
```

## 2. The degree-seven Macaulay rank gap

Consider the full degree-seven map

```text
Phi_7:S_5 (+) S_4 (+) S_3 -> S_7,
(A_5,B_4,D_3) |-> A_5 Q+B_4 C+D_3 F.                 (8)
```

The Hilbert series of the `(2,3)` complete intersection is

```text
(1-t^2)(1-t^3)/(1-t)^3.                              (9)
```

Its degree-three and degree-seven coefficients are both `6`.  Hence

```text
dim (Q,C)_7=36-6=30.                                 (10)
```

Modulo `(Q,C)_7`, the last summand in `(8)` is precisely

```text
mu_F:(S/(Q,C))_3 -> (S/(Q,C))_7.                     (11)
```

Multiplication by the nonzerodivisor `L` identifies both six-dimensional
graded pieces with `(3)`, and `(11)` becomes `m_f`.  Therefore

```text
rank(Phi_7)=30+rank(m_f)
           in {36,34,32,30}.                         (12)
```

In particular, rank `35` cannot occur.

Since an element of a finite Artin algebra is a unit exactly when it is
absent from every maximal ideal,

```text
Res(Q,C,F)!=0
 iff m_f is invertible
 iff rank(Phi_7)=36
 iff rank(Phi_7)>=35.                                (13)
```

This proves the pointwise one-minor gate:

```text
one nonzero 35-by-35 minor of Phi_7
       ==> Res(Q,C,F)!=0.                            (14)
```

Conversely, if the resultant vanishes, `(4)` forces corank at least
two, so every `35`-minor vanishes.

## 3. A canonical positive fifth-compound certificate

Fix monomial bases in `(8)` and define

```text
E_35(Phi_7)
 =sum_(|I|=|J|=35) det(Phi_7[I,J])^2.                (15)
```

There are exactly

```text
C(36,35) C(46,35)=480268195056                       (16)
```

summands.  Equations `(12)--(13)` give the exact real-locus identity

```text
E_35(Phi_7)>0 iff Res(Q,C,F)!=0.                     (17)
```

After choosing bases of the two quotient spaces in `(11)`, the smaller
equivalent object is the squared norm of the fifth compound,

```text
E_5(mu_F)=sum_(|I|=|J|=5) det(mu_F[I,J])^2.          (18)
```

Thus the determinant problem has a lower-order positive carrier:
singularity cannot hide at corank one.  Formula `(18)` is basis
dependent as a number but its strict positivity is not.

## 4. Factorial first-window application

For a translated four-slot support

```text
(n,n+a,n+b,n+M),             0<a<b<M,                (19)
```

THM-2843 identifies the first-window forms with a positive-definite
ternary quadratic `Q`, a cubic `C`, and a quartic `F`.  Positive
definiteness gives `V_R(Q)=empty`.  The restriction-to-a-three-slot-face
argument from THM-2824, already used in THM-2843, gives

```text
Q does not divide C,                                 (20)
```

so `(2)` holds at every physical depth.

Take the inherited `36`-row chart with

```text
20 Q-rows, 10 C-rows, 6 F-rows.                     (21)
```

Deleting one of the six `F` rows and any one of the `36` columns gives
`216` distinguished `35`-minors.  Any one which is nonzero closes that
depth by `(14)`.  Under the THM-2925 denominator clearing, its degree is
bounded by

```text
20(M-1)+10(2M-1)+5(3M-1)=55M-35.                    (22)
```

This is lower than the `58M-36` invoice for a full selected
determinant.  Testing these `216` cofactors is therefore the cheapest
fixed-chart experiment suggested by the theorem.

The implication is only sufficient for this selected bank.  If all
`216` vanish because that chart degenerates, another `35`-minor of the
full `36`-by-`46` map may still survive.

## 5. Sharp boundaries

The absence of real support is load-bearing.  In the reduced real
algebra

```text
A=R^6,              m_f=diag(0,1,1,1,1,1),          (23)
```

the determinant vanishes but the rank is `5`, and the fifth-compound
energy equals `1`.  Thus `(14)` is false when a single real local
factor can die.

Nor does `(14)` select a universal cofactor.  An invertible matrix may
have many zero `5`-cofactors; the identity matrix already has zero
off-diagonal cofactors.  A factorial closure still needs one of:

```text
* a fixed cofactor proved nonzero;
* a finite cofactor atlas with no common physical zero; or
* a direct positive lower bound for E_5 or E_35.      (24)
```

No such uniform statement, no new width closure, and no
coefficientwise positivity claim is made here.

## 6. Exact companion

The companion

```text
python 04-computation/gmc_conjugate_pair_corank_parity_thm2947.py
python -O 04-computation/gmc_conjugate_pair_corank_parity_thm2947.py
```

checks with integer arithmetic:

1. the degree-three and degree-seven Hilbert dimensions `6,6`;
2. the rank ladders `(6,4,2,0)` and `(36,34,32,30)`;
3. fifth-compound positivity for a full conjugate-block model;
4. fifth-compound vanishing after one conjugate pair dies;
5. the sharp real-summand hostile `(23)`; and
6. the exact full-bank size in `(16)`.

Both executions reproduce the stored transcript byte-for-byte.
