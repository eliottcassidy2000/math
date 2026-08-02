---
id: THM-3091
title: "Arbitrary-gap remote-pair desuspension and exact Jensen contraction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Above any
  fixed nonzero physical child resultant, a remote pair at C,C+H is positive
  uniformly for every integer H>=1 once C is large, with the threshold in C
  independent of H.  A second exact binary normalization contracts every
  physical outer error atom by Jensen's inequality, makes the far equation
  triangular, and converts the pair normal factor into the far one-normal
  carrier.  For H->infinity the resulting ledger is exactly the iterated
  one-normal ledger.  This is fixed-width and not an arbitrary two-scale
  support theorem.
source: root-gmc-arbitrary-gap-pair-2026-08-01
audit: >
  Two independent audits rederived the exact binary coefficients, response
  bound, Jensen inequalities, whole-system covariance, raw physical atom
  typing before signed cancellation, uniform finite-H lower-bound patch,
  separated resultant, carrier merger, parity, and scope.  Normal and
  optimized replay matched the 693-byte stored transcript and declared LF
  hashes; documentation passed.
depends_on:
  - THM-3063-terminal-suspension-transverse-resultant-and-five-slot-tail-holotopy
  - THM-3073-upper-triangular-resultant-norm-and-torsion-character-death-barcode
  - THM-3085-multi-normal-fixed-gap-cluster-and-unconditional-all-width-tail
  - THM-3089-logarithmic-moving-gap-cluster-cone-and-condition-number-boundary
related:
  - THM-3069-one-normal-remote-terminal-suspension-and-physical-tropical-flag
  - THM-3082-admissible-suspension-word-simultaneous-chambers-and-scale-tree-holotopy
  - THM-3086-arbitrary-cluster-composition-chambers-and-alternant-clutch-holotopy
script: 04-computation/gmc_arbitrary_gap_remote_pair_thm3091.py
output: 05-knowledge/results/gmc_arbitrary_gap_remote_pair_thm3091.out
script_sha256: 4a5abb0053f4d087044a2398c0fa3b8b43f197e1ba5f3677814ae4ebce167f2e
output_sha256: 32a747aad1ff2a6692c6a6c9713bdf4d9596732ca4e205fcf9040e88747184f6
hash_basis: LF-normalized bytes
---

# THM-3091 -- arbitrary-gap remote-pair desuspension

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3063 fixes the gap inside one remote pair.  THM-3089 moves a general
cluster only through square-root diameter because its all-high form must
remain close to a powered linear system.  A pair has a second triangular
quotient unavailable in higher normal rank.  More importantly, the same
secondary scaling acts on the **whole physical system**: exact Jensen bounds
show that it contracts, rather than enlarges, every far-gap multiplier in an
outer error atom.  The apparent linear boundary disappears completely.

## 1. Arbitrary-gap theorem

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

where `C,H` are positive integers.  For all sufficiently large `C`,
**uniformly for every**

```text
H>=1,                                                       (4)
```

the physical width-`k` first-window resultant satisfies

```text
R_k(C,H)>0.                                                (5)
```

The threshold in `C` can depend on the fixed child support and width, but is
independent of `H`.  Thus every arbitrary gap sequence

```text
H=H(C)>=1                                                  (6)
```

is covered.

Let

```text
U_r(C)=L(f_(n+C)^r)/L(f_n^r),
V_r(N,H)=L(f_(N+H)^r)/L(f_N^r),       N=n+C.             (7)
```

Put

```text
c_p=p(p+1)/(p(p+1)+1)<1,
rho_m=(m/(m+1))^m<1.                                      (7a)
```

Along every sequence with `H->infinity`, at an arbitrary rate relative to
`C`, the full carrier simplifies to

```text
R_k(C,H)=S_m^[pk]
 U_p(C)^[k!/p] U_k(C+H)^[(k-1)!]
 [1+O_(child,m)(c_p^H+poly(C)rho_m^C)].                  (8)
```

Thus the pair cluster has **desuspended** into the exact carrier ledger of
two iterated one-normal points.  Formula `(8)` is not asserted when `H`
stays bounded; Section 4 instead retains its positive exact normalized
binary factor.

When `m=3`, THM-2824 supplies `S_3>0` for every arbitrary physical
three-slot base.  Hence every such base followed by `(C,C+H)` has positive
five-slot resultant for **every** positive gap `H`, uniformly once `C` passes
one base-dependent threshold.

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

## 3. Secondary normalization and whole-system covariance

Put

```text
lambda=V_p^(-1/p),
s=V_p^(k/p)/V_k,

B_p(u,v)=A_p(u,lambda v),
B_k(u,v)=s A_k(u,lambda v).                             (12)
```

The `v^p` coefficient of `B_p` and the `v^k` coefficient of `B_k` are both
exactly one; the `u^p` coefficient of `B_p` is also exactly one.  Apply this
change not merely to the binary face but to the **whole outer-transformed
physical system**: scale its far direct coordinate `v` by `lambda` and its
degree-`k` equation by `s`.

The binary resultant has variable-covariance degree `pk`, and its degree in
the coefficients of the degree-`k` equation is `p`.  Hence, exactly,

```text
Etilde_(C,H):=Res(B_p,B_k)
 =lambda^(pk)s^p E_(C,H)
 =E_(C,H)/V_k(N,H)^p.                                  (13)
```

On the full system, the second covariance factor is

```text
lambda^(k!)s^[(k-1)!]=V_k(N,H)^[-(k-1)!].               (13a)
```

Combining `(13a)` with THM-3085's first outer scaling gives

```text
Rhat_k(C,H)=R_k(C,H)/[
 U_p(C)^[k!/p] U_k(C)^[(k-1)!] V_k(N,H)^[(k-1)!]].      (13b)
```

All powers of `V_p` cancel.  These are exact covariance identities, not
asymptotic normalizations.

## 4. Exact Jensen contraction and the triangular quotient

Let `D_H=(N+1)_H`.  Blocking the numerator of `V_r` into `H` packets gives

```text
V_r^(1/r)=product_(t=0)^(H-1) G_r(N+t)/(N+t+1),
G_r(x)=[product_(a=1)^r (rx+a)]^(1/r).                  (14)
```

Since `k=p+1` and `x>=p`,

```text
G_p(x)<=p(x+1),
G_k(x)>=kx+1,
G_p(x)/G_k(x)
 <=p(x+1)/[(p+1)x+1]
 <=p(p+1)/[p(p+1)+1]=c_p.                              (15)
```

The last rational function decreases in `x`.  Therefore the response ratio
obeys the uniform exponential bound

```text
R_(N,H):=V_p^(1/p)/V_k^(1/k)<=c_p^H<1                  (16)
```

for every `N>=p` and every `H>=1`.

There is a second exact inequality.  The function

```text
F_r(t)=log (rN+1)_t
```

has increasing discrete increments and `F_r(0)=0`.  Convexity on
`0<=j<=r` gives

```text
F_r(jH)<=(j/r)F_r(rH),
Q_(r,j)^r<=V_r^j.                                       (17)
```

Consequently every coefficient of `B_p` is bounded by its binomial
coefficient, while for `0<=j<k`,

```text
|[u^(k-j)v^j]B_k|
 <=binom(k,j)V_k^(j/k)V_p^((k-j)/p)/V_k
 =binom(k,j)R_(N,H)^(k-j)
 <=binom(k,j)c_p^[H(k-j)].                              (18)
```

In the quotient where the non-top coefficients of `B_k` vanish,

```text
Res(B_p,v^k)=B_p(1,0)^k=1.                              (19)
```

Every other monomial of the fixed universal resultant contains a non-top
coefficient of `B_k`, whereas all coefficients of `B_p` remain uniformly
bounded.  Hence

```text
Etilde_(C,H)=1+O_m(c_p^H)                               (20)
```

for every `N>=p`.  This is grouping before estimation; no growing
coefficient of the first equation has been hidden in the error.

For a uniform lower bound, choose `H_0` so that the error in `(20)` is at
most `1/2` for `H>=H_0`.  For each of the finitely many `1<=H<H_0`, as
`N->infinity`,

```text
(B_p,B_k)->((u+v)^p,((p/k)^H u+v)^k),
Etilde_(C,H)->[1-(p/k)^H]^(pk)>0.                       (21)
```

Taking one maximum threshold and half the finite minimum in `(21)` yields
constants `C_0` and `eta_m>0` such that

```text
Etilde_(C,H)>=eta_m                                     (22)
```

for every `C>=C_0` and every integer `H>=1`.

## 5. Exact outer atom invoice and positivity

It remains to justify applying `(12)` to the whole physical system.  Work in
THM-3085's fixed-pivot/direct coordinates before any signed inclusion sum.
An expansion atom in a row `r` has `j` actual remote factors, `x` of them at
the far point, total normal degree `ell`, and far-coordinate degree `t`.
The coordinate construction gives

```text
ell>=j,                    t>=x.                         (23)
```

Relative to moving those `x` factors back from `C+H` to `C`, its exact extra
factor is

```text
W=(jN+D+1)_(xH)/(N+1)_H^x,                              (24)
```

where, if `A` is the total offset of the lower physical factors,

```text
D=(r-j)n+A
```

is fixed by the child atom.  Once `C` exceeds one child-dependent threshold,
`jN+D<=pN` in every row `r<=p`, and `jN+D<=kN` in the degree-`k` row.  Thus

```text
W<=Q_(p,x)                    (r<=p),
W<=Q_(k,x)                    (r=k).                     (25)
```

Every factor in the blocked numerator for `V_p` is strictly larger than its
matching denominator factor, so `V_p>1` and `lambda<=1`.  Equation `(17)`
therefore implies, atom by atom,

```text
W lambda^t
 <=Q_(p,x)V_p^(-x/p)<=1                 (r<=p),

W lambda^t s
 <=Q_(k,x)V_p^((k-x)/p)/V_k
 <=R_(N,H)^(k-x)<=1                     (r=k).           (26)
```

This estimate is unsigned and precedes every inclusion cancellation.  Hence
the second normalization cannot enlarge any nonsurviving layer.  We use the
atomwise `(j,ell)` proof inside THM-3085, not a repeated-gap theorem at a
formal `H=0`: its declared lower block and complete all-high forms are the
base-one banks, and every other atom pays the strict outer gap.  That gap
survives with no `exp(O(H))` invoice:

```text
rho_m=(m/(m+1))^m<1,

(whole transformed system)
 =(H_2,...,H_m,H_p+B_p,B_k)
   +O_(child,m)(poly(C)rho_m^C),                         (27)
```

uniformly for **all** `H>=1`.  THM-3073 gives the exact reference resultant

```text
S_m^[pk] Etilde_(C,H)^m!.                                (28)
```

Its modulus is uniformly bounded below by `(22)`; all reference coefficients
are uniformly bounded by `(17)--(18)`.  Fixed-degree resultant continuity
therefore turns `(27)` into a uniform relative error.  Undoing `(13b)` and
using

```text
U_k(C)V_k(n+C,H)=U_k(C+H),
p m!=(k-1)!                                              (29)
```

gives, for every `H>=1`,

```text
R_k(C,H)=S_m^[pk] Etilde_(C,H)^m!
 U_p(C)^[k!/p] U_k(C+H)^[(k-1)!]
 [1+O_(child,m)(poly(C)rho_m^C)]>0.                     (30)
```

Both `pk` and `m!` are even and every carrier is positive.  Equations
`(20)` and `(30)` also prove the desuspended asymptotic `(8)`.

## 6. The macroscopic binary face as a consistency check

The exact Jensen proof already covers every relative scale.  The earlier
macroscopic calculation remains useful because it identifies the limiting
normal face directly.  Suppose `H/N->delta>0`.  The per-`N` rate of
`Q_(r,j)` is

```text
q_(r,j)(delta)
 =(r+j delta)log(r+j delta)-r log r
  -j(1+delta)log(1+delta).                               (31)
```

After `(12)`, the mixed degree-`p` coefficient rate is

```text
S_j(delta)=q_(p,j)(delta)-j delta log p.                (32)
```

It is strictly convex in `j`, with `S_0=S_p=0`; hence `S_j<0` for
`0<j<p`.  The degree-`k` coefficient rate is

```text
T_j(delta)=q_(k,j)(delta)
 +(k-j)delta log p-k delta log k.                        (33)
```

It obeys

```text
T_0=-k delta log(k/p)<0,
T_k=0,
d^2 S_j/dj^2=delta^2/(p+j delta)>0,
d^2 T_j/dj^2=delta^2/(k+j delta)>0.                      (34)
```

Strict convexity places every `T_j`, `0<j<k`, below the chord joining the
endpoints, hence `T_j<0`.  Thus at every fixed positive ratio the whole
secondary binary system tends coefficientwise to

```text
(u^p+v^p, v^k),                                          (35)
```

whose resultant is one.  The uniform proof in Sections 4--5 shows why this
face remains lawful even when `H/N` has no limit or tends to infinity.

## 7. Holotopy, boundaries, and exact evidence

At fixed gap, the pair is one two-normal node with an exact alternant
sidecar.  As `H->infinity`, `(29)` changes the bracketing to two one-normal
nodes without changing the final factorial ledger.  This is a
physical associativity/desuspension holotopy: the pair alternant is converted
into the far `U_k` carrier rather than discarded.

The theorem fixes the child support and width, takes `C->infinity`, and
requires the two remote offsets to remain distinct.  It does not cover a
moving child, growing width, repeated pair offsets, `S_m=0`, arbitrary
two-scale supports, or a general `q>=3` cluster; the latter is isolated in
reserved THM-3093.  It does not prove arbitrary supports, arbitrary-radial
GMC(2), NC2, LRC(14), JC(2), or DC(2).

The exact companion verifies:

1. exact secondary coefficient controls in rational `p`-th powers,
   including every Jensen inequality and the response bound `(16)`;
2. `15,705` exact whole-system outer-atom controls for `(25)--(26)`,
   including raw physical starts, maximal admissible lower displacement, and
   `t>x` normal degrees;
3. eighteen exact covariance ledgers giving `(13)`;
4. fifty-six independent Sylvester determinants for the arbitrary
   triangular quotient `(19)`;
5. eighty exact carrier identities `(29)` and every factorial exponent;
6. `208` rigorous rational-log degree-`p` mixed cells and `272` degree-`k`
   non-top cells, all strictly negative, together with both convexities in
   `(34)`.

Run

```text
python 04-computation/gmc_arbitrary_gap_remote_pair_thm3091.py
python -O 04-computation/gmc_arbitrary_gap_remote_pair_thm3091.py
```

Both modes must equal the stored transcript after LF normalization.

**QED.**
