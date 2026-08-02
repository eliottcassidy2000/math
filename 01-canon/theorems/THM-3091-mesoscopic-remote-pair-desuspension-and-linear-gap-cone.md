---
id: THM-3091
title: "Mesoscopic remote-pair desuspension and linear-gap cone"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.  Above any
  fixed nonzero physical child resultant, a remote pair at C,C+H is positive
  uniformly throughout a nonempty linear wedge 1<=H<=kappa_m C.  A second
  exact binary normalization makes the far equation tend triangularly to a
  pure power and converts the pair normal factor into the far one-normal
  carrier.  For H->infinity the resulting ledger is exactly the iterated
  one-normal ledger.  The binary normal face remains stable at every fixed
  H/C; the present restriction comes from the outer child-to-pair
  contraction.  This is fixed-width and not an arbitrary two-scale theorem.
source: root-gmc-mesoscopic-pair-2026-08-01
depends_on:
  - THM-3063-terminal-suspension-transverse-resultant-and-five-slot-tail-holotopy
  - THM-3073-upper-triangular-resultant-norm-and-torsion-character-death-barcode
  - THM-3085-multi-normal-fixed-gap-cluster-and-unconditional-all-width-tail
  - THM-3089-logarithmic-moving-gap-cluster-cone-and-condition-number-boundary
related:
  - THM-3069-one-normal-remote-terminal-suspension-and-physical-tropical-flag
  - THM-3082-admissible-suspension-word-simultaneous-chambers-and-scale-tree-holotopy
  - THM-3086-arbitrary-cluster-composition-chambers-and-alternant-clutch-holotopy
script: 04-computation/gmc_mesoscopic_remote_pair_thm3091.py
output: 05-knowledge/results/gmc_mesoscopic_remote_pair_thm3091.out
script_sha256: b0f6be980638f3c33bad26097972bebfd5751a54a4e29bfdf82fe715cb6cf188
output_sha256: 4854ee02b27351adeda8e0e87d2a82d9dd8df1251661c6163d088bf7069df340
hash_basis: LF-normalized bytes
---

# THM-3091 -- mesoscopic remote-pair desuspension

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.**

THM-3063 fixes the gap inside one remote pair.  THM-3089 moves a general
cluster only through square-root diameter because its all-high form must
remain close to a powered linear system.  A pair has a second triangular
quotient unavailable in higher normal rank.  It absorbs every mixed term of
the first upper equation and pushes the physical cone all the way to a
positive linear angle.

## 1. Linear-wedge theorem

Fix a physical lower support of width `m>=3` with normalized first-window
resultant

```text
S_m!=0.                                                     (1)
```

Put

```text
p=m+1,                   k=m+2.                           (2)
```

Append the two offsets

```text
C,                       C+H,                             (3)
```

where `C,H` are positive integers.  There is a constant
`kappa_m>0`, depending only on the fixed width, such that for all sufficiently
large `C`, **uniformly** for

```text
1<=H<=kappa_m C,                                           (4)
```

the physical width-`k` first-window resultant satisfies

```text
R_k(C,H)>0.                                                (5)
```

The threshold in `C` can depend on the fixed child support and on the chosen
`kappa_m`, but not on `H` inside `(4)`.  In particular every arbitrary gap
sequence

```text
H(C)=o(C)                                                  (6)
```

is covered.

Let

```text
U_r(C)=L(f_(n+C)^r)/L(f_n^r),
V_r(N,H)=L(f_(N+H)^r)/L(f_N^r),       N=n+C.             (7)
```

Along every sequence with `H->infinity` in `(4)`, the full carrier simplifies
to

```text
R_k(C,H)=S_m^[pk]
 U_p(C)^[k!/p] U_k(C+H)^[(k-1)!]
 [1+O(poly(C,H)e^(-c_m H)+poly(C,H)e^(-c'_m C))].        (8)
```

Thus the pair cluster has **desuspended** into the exact carrier ledger of
two iterated one-normal points.  Formula `(8)` is not asserted when `H`
stays bounded; that finite bank is positive by THM-3085 and retains its
exact fixed-gap normal factor.

When `m=3`, THM-2824 supplies `S_3>0` for every arbitrary physical
three-slot base.  Hence every such base followed by `(C,C+H(C))` has positive
five-slot resultant for every sublinear gap sequence, and indeed throughout
one fixed small linear wedge.

## 2. Exact binary all-high block

Work first in the direct variables `u,v` corresponding to `f_N` and
`f_(N+H)`.  The fixed-pivot binary normal basis differs by a determinant-one
change, so its resultant is the same.  For `r=p,k`, write the exact normalized
all-high form as

```text
A_r(u,v)=sum_(j=0)^r binom(r,j) Q_(r,j) u^(r-j)v^j,

Q_(r,j)=(rN+1)_(jH)/(N+1)_H^j.                          (9)
```

Here `(x)_s` is rising factorial.  The endpoint coefficients are

```text
Q_(r,0)=1,              Q_(r,r)=V_r(N,H).               (10)
```

Let

```text
E_(C,H)=Res_bin(A_p,A_k)                                (11)
```

be the exact pair normal factor in THM-3085's outer carrier.

## 3. Secondary normalization and exact covariance

Put

```text
lambda=V_p^(-1/p),
s=V_p^(k/p)/V_k,

B_p(u,v)=A_p(u,lambda v),
B_k(u,v)=s A_k(u,lambda v).                             (12)
```

The `v^p` coefficient of `B_p` and the `v^k` coefficient of `B_k` are both
exactly one; the `u^p` coefficient of `B_p` is also exactly one.

The binary resultant has variable-covariance degree `pk`, and its degree in
the coefficients of the degree-`k` equation is `p`.  Therefore

```text
Res(B_p,B_k)
 =lambda^(pk)s^p E_(C,H)
 =E_(C,H)/V_k(N,H)^p.                                  (13)
```

The powers of `V_p` cancel exactly.  This identity is not an asymptotic
normalization.

## 4. Coefficient rates and the triangular quotient

The same rising-factorial calculation as THM-3089 gives, uniformly for
`0<=j<=r`,

```text
log Q_(r,j)=jH log r+O_r(H^2/N).                        (14)
```

Consequently every coefficient of `B_p` has logarithmic size

```text
O_m(H^2/N),                                              (15)
```

while the coefficient of `u^(k-j)v^j` in `B_k` satisfies

```text
log |coeff_(k,j)(B_k)|
 =-(k-j)H log(k/p)+O_m(H^2/N),        0<=j<k.            (16)
```

Choose `kappa_m` small enough.  For `H<=kappa_m N`, the error in `(16)` is
at most a fixed small fraction of its negative main term.  Coefficients of
`B_p` can grow, but only like `exp(O_m(kappa_m H))`.

Now take the exact quotient in which every non-top coefficient of `B_k`
vanishes.  It is triangular:

```text
Res(B_p,v^k)=B_p(1,0)^k=1,                              (17)
```

independently of **all** other coefficients of `B_p`.  Every resultant term
outside `(17)` contains at least one non-top coefficient from `B_k` and only
a fixed number of coefficients from `B_p`.  Shrinking `kappa_m` once more,
the finite resultant bank and `(15)--(16)` give, uniformly when
`H>=H_0` and `(4)` holds,

```text
E_(C,H)/V_k(N,H)^p
 =1+O(poly(C,H)e^(-c_m H))>0.                           (18)
```

This is grouping before estimation: coefficientwise convergence of `B_k`
alone would not justify `(18)` if the arbitrary `B_p` coefficients were
allowed to grow without the triangular quotient.

## 5. Outer contraction and bounded-gap patch

THM-3085's outer contraction at the common base `C` has gap

```text
rho_m=(m/(m+1))^m<1.                                    (19)
```

THM-3089's pre-cancellation moving-gap invoice bounds every complementary
layer by

```text
poly(C,H)exp(a_m H)rho_m^C.                              (20)
```

Shrink `kappa_m` so that `a_m kappa_m<-(1/2)log rho_m`.
Then `(20)=O(poly(C,H)e^(-c'_m C))`.  The exact separated quotient is

```text
S_m^[pk] E_(C,H)^m!,                                    (21)
```

and the physical covariance contributes

```text
U_p(C)^[k!/p] U_k(C)^[(k-1)!].                          (22)
```

Since `pk` and `m!` are even, all carrier factors in `(21)--(22)` are
positive.  Equations `(18)--(22)` prove `(5)` uniformly for `H>=H_0`.

The remaining gaps `1<=H<H_0` form a finite bank.  THM-3085 applies to each
fixed value, and the maximum of their finitely many thresholds is uniform.
This completes the entire wedge `(4)`.

Finally,

```text
U_k(C)V_k(n+C,H)=U_k(C+H),
p m!=(k-1)!.                                             (23)
```

Raising `(18)` to `m!` and inserting `(23)` into `(21)--(22)` proves the
desuspended carrier `(8)`.

## 6. The binary normal face survives fixed ratios

The smallness of `kappa_m` comes from the **outer** contraction, not from the
binary normal block.  To see the distinction, suppose `H/N->delta>0`.  The
per-`N` rate of `Q_(r,j)` is

```text
q_(r,j)(delta)
 =(r+j delta)log(r+j delta)-r log r
  -j(1+delta)log(1+delta).                               (24)
```

After `(12)`, the mixed degree-`p` coefficient rate is

```text
S_j(delta)=q_(p,j)(delta)-j delta log p.                (24a)
```

It is strictly convex in `j`, with `S_0=S_p=0`; hence `S_j<0` for
`0<j<p`.  The degree-`k` coefficient rate is

```text
T_j(delta)=q_(k,j)(delta)
 +(k-j)delta log p-k delta log k.                        (25)
```

It obeys

```text
T_0=-k delta log(k/p)<0,
T_k=0,
d^2 S_j/dj^2=delta^2/(p+j delta)>0,
d^2 T_j/dj^2=delta^2/(k+j delta)>0.                      (26)
```

Strict convexity places every `T_j`, `0<j<k`, below the chord joining the
endpoints, hence `T_j<0`.  Thus at every fixed positive ratio the whole
secondary binary system tends coefficientwise to

```text
(u^p+v^p, v^k),                                          (27)
```

whose resultant is one.  Extending the physical theorem to a larger linear
cone requires a sharper signed analysis of the outer child-to-pair entropy
face, not another binary normal argument.

## 7. Holotopy, boundaries, and exact evidence

At fixed gap, the pair is one two-normal node with an exact alternant
sidecar.  As `H->infinity` inside `(4)`, `(23)` changes the bracketing to two
one-normal nodes without changing the final factorial ledger.  This is a
physical associativity/desuspension holotopy: the pair alternant is converted
into the far `U_k` carrier rather than discarded.

The theorem fixes the child support and width.  It proves only the existence
of a nonzero linear angle, not its optimal value.  It does not cover every
fixed `H/C`, arbitrary two-scale supports, a general `q>=3` linear cluster,
wall crossing of the outer entropy face, a width-uniform cone, arbitrary
supports, arbitrary-radial GMC(2), NC2, LRC(14), JC(2), or DC(2).

The exact companion verifies:

1. `112` exact secondary coefficient controls in rational `p`-th powers,
   including the degree-`p` growth and every degree-`k` decay rate;
2. eighteen exact covariance ledgers giving `(13)`;
3. fifty-six independent Sylvester determinants for the arbitrary
   triangular quotient `(17)`;
4. eighty exact carrier identities `(23)` and every factorial exponent;
5. `208` rigorous rational-log degree-`p` mixed cells and `272` degree-`k`
   non-top cells, all strictly negative, together with both convexities in
   `(26)`.

Run

```text
python 04-computation/gmc_mesoscopic_remote_pair_thm3091.py
python -O 04-computation/gmc_mesoscopic_remote_pair_thm3091.py
```

Both modes must equal the stored transcript after LF normalization.

**QED, pending independent audit and status promotion.**
