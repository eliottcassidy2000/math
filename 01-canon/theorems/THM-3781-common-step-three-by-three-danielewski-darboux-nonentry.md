---
id: THM-3781
title: "Complete scalar-centred three-by-three arithmetic-progression Danielewski Darboux nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY DERIVED; PENDING HOSTILE AUDIT OF
  THE FULL GENERALIZATION.  On every exponent-two squarefree
  Danielewski surface with at least two arms, no Darboux pair can have three
  nonzero homogeneous weight pieces in each output, arithmetic-progression
  supports with arbitrary positive steps, and scalar-centred convolution.
  Unequal steps leave a forbidden lonely mixed-sign bracket.  For equal
  steps, endpoint gcds are unrestricted: algebraic roots of the two endpoint
  owners make both adjacent equations integrable, and infinity contradicts
  the negative middle membership invoice.  The independently audited
  step-three theorem is a strict special case.  This closes the complete
  scalar-centred three-by-three AP census, not planar JC(2).
source: root / planar-Jacobian Danielewski Darboux session, 2026-08-23
audit: >
  The d=3 specialization was independently hostile-audited by
  jc_zero_debt_lift (2026-08-23, bfad661005), including support typing,
  endpoint UFD laws, adjacent integrations, arm valuations, and terminal
  degree exits.  The same auditor independently rederived both the all-d
  radical equations/infinity case split and the unequal-step lonely-weight
  reduction.  Exact-text and frozen-hash audit of the strengthened statement
  remains pending.  The companion checks 2,470 unequal cells through step
  20 (including 100 doubling hostiles) and every equal-step cell through
  d=30 (including 135 nonprimitive endpoints), with 47,466 active gates
  under normal and optimized Python.
depends_on:
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
related:
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
script: 04-computation/jc2_common_step_three_by_three_danielewski_thm3781.py
output: 05-knowledge/results/jc2_common_step_three_by_three_danielewski_thm3781.out
script_sha256: 066fab6dee42f56bb32eb3916e2f05c130059750cc27c30ed8de9c5a831495cb
output_sha256: 44b00c69be0493b4f4e49f2eb6c0cc7cd10cfad7c43c8bfca14c248cc11dc425
semantic_sha256: 0f8e172b3c6ebc8679753322e1c5b5ee717cc029bcaff517c6f87b2f1ca96073
hash_basis: raw LF bytes
---

# THM-3781 -- every scalar-centred three-by-three AP cell is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY DERIVED; PENDING HOSTILE AUDIT OF
THE FULL GENERALIZATION.**  A lonely-weight argument closes unequal steps;
a common endpoint-owner calculation closes equal steps.  The essential move
in the latter case is to retain endpoint gcds by working temporarily in a
finite algebraic extension of `C(b)`.

Work over `C`.  Let `Sigma in C[b]` be squarefree with at least two distinct
roots, and put

```text
A_Sigma=C[b,c,e]/(c^2e-Sigma(b)).                       (1)
```

Use the Poisson bracket and weights

```text
{b,c}=c^2,          {c,e}=-Sigma'(b),
{b,e}=-2ce,
wt(b)=0,            wt(c)=1,             wt(e)=-2.     (2)
```

The weight-`r` part and homogeneous bracket are

```text
(A_Sigma)_r={c^r f(b):
  Sigma^ceil(-r/2) divides f when r<0},                 (3)

{c^r f,c^s g}=c^(r+s+1)(s f'g-r f g').                (4)
```

Subtract constants from two proposed outputs.  Suppose that each remaining
output has exactly three nonzero homogeneous pieces, its weight support is
an arithmetic progression with positive step `r` or `s`, and the middle
weights `p,q` satisfy

```text
p+q+1=0.                                               (5)
```

Thus scalar output weight zero is the central coefficient of the support
convolution.  The theorem says that no such `P,Q` satisfy

```text
{P,Q}=1.                                               (6)
```

## 1. The complete support family and the unequal-step exit

The extreme output weights occur once, so their brackets must vanish.  A
negative and a positive homogeneous piece cannot have zero bracket: at a
root of `Sigma`, if their coefficient orders are `R>0,S>=0`, then the first
surviving coefficient in `(4)` is proportional to

```text
sR-rS>0                 when r<0<s.                    (7)
```

If one extreme has weight zero and the other does not, `(4)` instead forces
the zero-weight coefficient to be constant, which was removed.  Hence both
lower endpoints are negative and both upper endpoints are positive.

After exchanging the outputs and negating one if needed to restore bracket
`+1`, write

```text
p=-a,       q=a-1,       1<=a<=min(s,r-1).            (8)
```

Every remaining support cell is therefore

```text
supp(P)=(-a-r, -a, -a+r),
supp(Q)=(a-1-s, a-1, a-1+s).                          (9)
```

The nine output weights, before collisions, are

```text
ir+js,                  i,j in {-1,0,1}.              (10)
```

Suppose first that `r!=s`.  Output weight `-r` comes only from the lower
piece of `P` and middle piece of `Q`, except when

```text
s=r  or  s=2r.                                        (11)
```

If `a>1`, that unique pair has negative and positive weights, so its bracket
cannot vanish by `(7)`.  If `a=1`, the second piece has weight zero; its
coefficient is nonconstant because constants were removed, and `(4)` again
makes the bracket nonzero.

Likewise, output weight `+s` comes only from the middle piece of `P` and
upper piece of `Q`, except when

```text
r=s  or  r=2s.                                        (12)
```

This pair always has negative and positive weights, hence cannot bracket to
zero.  An unequal pair of steps could evade the first lonely weight only by
`s=2r` and the second only by `r=2s`; it cannot evade both.  Thus every
unequal-step cell is empty before its endpoint equations are used.

It remains to take

```text
r=s=d>=2,       1<=a<=d-1.                            (13)
```

Then `(9)` becomes

```text
supp(P)=(-(d+a), -a, d-a),
supp(Q)=(-(d-a+1), a-1, d+a-1),                       (14)
```

with multiplicities `1,2,3,2,1` at weights
`-2d,-d,0,d,2d`.  There are `d-1` cells up to exchange.  No endpoint
coprimality has been imposed.

## 2. Endpoint owners, including nonprimitive exponents

Abbreviate

```text
m=d+a,       n=d-a+1,       R=d-a,       S=d+a-1,
u=gcd(m,n),  v=gcd(R,S).                               (15)
```

Write the proposed pair as

```text
P=c^-m f+c^-a alpha+c^R F,
Q=c^-n g+c^(a-1) beta+c^S H.                          (16)
```

The endpoint equations at weights `-2d` and `2d` are

```text
-n f'g+m fg'=0,                 S F'H-R FH'=0.         (17)
```

Unique factorization in `C[b]` gives nonzero constants `A,B,L,M` and
polynomials `k,ell` such that

```text
f=A k^(m/u),       g=B k^(n/u),
F=L ell^(R/v),     H=M ell^(S/v).                     (18)
```

The negative-weight membership invoices in `(3)` force every root of
`Sigma` to divide `k`; in particular

```text
Sigma|k,                  deg(k)>0.                    (19)
```

Choose a finite algebraic function-field extension `E/C(b)` containing
elements `K,N` with

```text
K^u=k,                    N^v=ell.                     (20)
```

Then `(18)` becomes

```text
f=A K^m,       g=B K^n,       F=L N^R,       H=M N^S. (21)
```

The derivation `d/db` extends uniquely to `E`.  Its constant field is `C`:
it is algebraic over the constant field of `C(b)`, and `C` is algebraically
closed.  This is the only fact needed from the radical extension.

## 3. The two adjacent equations integrate uniformly

At output weight `-d`, formula `(4)` and `(21)` give, after removing the
nonzero factor `K^(d+1)`,

```text
A m (K^(a-1) beta)' - B n (alpha K^-a)'=0.             (22)
```

At output weight `d`, removing `N^(d-1)` gives

```text
M S (alpha N^a)' - L R (beta N^(-(a-1)))'=0.           (23)
```

Since the constants of `E` are `C`, these integrate to

```text
alpha/K^a=lambda K^(a-1) beta+C,
beta/N^(a-1)=mu alpha N^a+D,                           (24)
```

where

```text
lambda=A m/(B n),       mu=M S/(L R)                  (25)
```

are nonzero and `C,D in C`.  Put `Z=KN` and `nu=lambda mu`.  Eliminating
`beta` produces the universal obstruction

```text
alpha(1-nu Z^(2a-1))
  =K^a(lambda D Z^(a-1)+C).                           (26)
```

This equation contains no `d`; the common step has disappeared.

## 4. Infinity contradicts the negative middle invoice

Choose any place `w` of `E` above the infinite place of `C(b)`, and divide
its pole order by the ramification index.  Denote the resulting normalized
pole degrees of `K,N` by

```text
kappa=deg(k)/u>0,             eta=deg(ell)/v>=0.       (27)
```

Thus `Z` has pole degree `kappa+eta>0`, and the left factor
`1-nu Z^(2a-1)` in `(26)` has pole degree
`(2a-1)(kappa+eta)`.

If `D!=0` and `a>1`, the `D Z^(a-1)` term strictly dominates `C` at this
place.  Comparing pole orders in `(26)` gives

```text
deg(alpha)+(2a-1)(kappa+eta)
  =a kappa+(a-1)(kappa+eta),

deg(alpha)=-a eta<=0.                                 (28)
```

If the surviving parenthesis on the right of `(26)` is a nonzero constant
-- either `D=0,C!=0`, or `a=1,lambda D+C!=0` -- the same comparison gives

```text
deg(alpha)=-(a-1)kappa-(2a-1)eta<=0.                  (29)
```

If that parenthesis is zero, `(26)` forces `alpha=0`, because `Z` has a
pole and hence `1-nu Z^(2a-1)` is not zero.

Every case is impossible.  Indeed, `c^-a alpha` was declared a nonzero
weight piece, and `(3)` requires

```text
Sigma^ceil(a/2) divides alpha.                         (30)
```

Consequently `deg(alpha)>0`, contradicting `(28)`, `(29)`, or
`alpha=0`.  The scalar weight equation is never reached.

## 5. Boundary and surviving route

For `d=3`, `(14)` gives precisely the two previously audited cells

```text
(-4,-1,2) x (-3,0,3),
(-5,-2,1) x (-2,1,4).                                 (31)
```

The new radical-owner argument also covers every nonprimitive endpoint
pair, beginning with `(d,a)=(4,2)`, as well as all larger common steps.

This is the full scalar-centred three-by-three arithmetic-progression
census.  A noncentral occurrence of scalar weight zero, genuinely gapped
supports, and supports with at least four pieces remain open.  It constructs
no Darboux pair and no planar Jacobian counterexample.  **QED, with the full
exact-text hostile audit pending.**
