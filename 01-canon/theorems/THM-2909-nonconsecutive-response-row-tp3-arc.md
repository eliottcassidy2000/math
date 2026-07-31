---
id: THM-2909
title: "Nonconsecutive response-row TP3 arc"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For every nonzero nonnegative two-ray cone with support gap one or two,
  every ordered triple of factorial-response rows at least one has positive
  three-column determinant.  The mechanism is an exact weighted-chord
  decomposition: THM-2873's consecutive slope gaps make every long left
  chord flatter than every later right chord.  The row cutoff and cone
  positivity are sharp: d_4 has negative determinant on rows 0,1,2, while
  d_0-(1/5)d_1 has positive first-column values on rows 1,2,3 but determinant
  -1728/15625.
source: root/tp3-arbitrary-row-audit-2026-07-29
audit: >
  An independent hostile audit rederived the weighted-chord determinant,
  both sharp failure matrices, every dependency direction and scope
  boundary, and replayed normal, optimized and stored output with the
  declared LF hashes.
depends_on:
  - THM-2830-disjoint-positive-adjacent-cone-factorial-moment-three-detection
  - THM-2866-positive-factorial-difference-semiring-and-cubic-pascal-response-ladder
  - THM-2873-two-ray-factorial-response-tp3-curvature
related:
  - THM-2906-atomic-tp3-does-not-orient-mixed-endpoint-holonomy
script: 04-computation/gmc_nonconsecutive_response_row_tp3_thm2909.py
output: 05-knowledge/results/gmc_nonconsecutive_response_row_tp3_thm2909.out
script_sha256: 70ae887f60cc9e1428d158e668e89342a409671190bb2e6cb77f7b5f88387a14
output_sha256: 3b086ca71d33f96b71c00da3c57cd708293e791814c84a134bdcee979eaee50d
hash_basis: LF-normalized bytes
---

# THM-2909 -- nonconsecutive response-row TP3 arc

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Put

```text
L(s^k)=k!,                    f_j=s^j/j!,
d_j=f_(j+1)-f_j.                                      (1)
```

For a fixed polynomial `W`, define its first three response columns

```text
A_r=L(d_rW),             B_r=L(d_rW^2),
C_r=L(d_rW^3).                                      (2)
```

The theorem applies to

```text
W=x d_p+y d_(p+h),

p>=0,             h in {1,2},
x,y>=0,           (x,y)!=(0,0).                       (3)
```

For every three integers

```text
1<=r_0<r_1<r_2,                                    (4)
```

one has

```text
D_(r_0,r_1,r_2)(W)
 :=det [
   A_(r_0)  B_(r_0)  C_(r_0)
   A_(r_1)  B_(r_1)  C_(r_1)
   A_(r_2)  B_(r_2)  C_(r_2)
 ] >0.                                                (5)
```

Thus THM-2873's consecutive three-row curvature extends to every
nonconsecutive response-row triple.  This is a theorem about the row
arc for the same gap-one and gap-two cones.  It is not an extension to
arbitrary support gaps.

## 1. The weighted-chord lemma

The mechanism is an abstract discrete-convexity fact.  Suppose

```text
A_r>0,                R_r=B_r/A_r,
S_r=C_r/A_r,          delta_r=R_(r+1)-R_r>0.          (6)
```

Define the consecutive response slope

```text
lambda_r=(S_(r+1)-S_r)/delta_r.                       (7)
```

For `i<j`, the chord slope is the positive weighted average

```text
Lambda_(i,j)
 =(S_j-S_i)/(R_j-R_i)
 =[
    sum_(t=i)^(j-1) delta_t lambda_t
  ]/
  [
    sum_(t=i)^(j-1) delta_t
  ].                                                  (8)
```

If

```text
lambda_(r+1)>lambda_r                 for every r>=1, (9)
```

then, for `1<=r_0<r_1<r_2`, every slope occurring in the
left average `Lambda_(r_0,r_1)` is strictly smaller than every slope
occurring in the right average `Lambda_(r_1,r_2)`.  Hence

```text
Lambda_(r_0,r_1)<Lambda_(r_1,r_2).                   (10)
```

There is an exact denominator-free version.  Normalizing the three rows
of `(5)` by their positive first entries gives

```text
D_(r_0,r_1,r_2)/(A_(r_0)A_(r_1)A_(r_2))

 =sum_(u=r_0)^(r_1-1) sum_(v=r_1)^(r_2-1)
    delta_u delta_v (lambda_v-lambda_u).              (11)
```

Indeed, the normalized determinant is

```text
(R_(r_1)-R_(r_0))(S_(r_2)-S_(r_1))
 -(S_(r_1)-S_(r_0))(R_(r_2)-R_(r_1)),                (12)
```

and expanding the four differences gives `(11)`.  Every summand in
`(11)` is strictly positive under `(6)` and `(9)`.  This proves the
lemma without a limiting argument or a hidden continuity assumption.

## 2. Application to the factorial response

For every cone `(3)`, THM-2830 gives

```text
A_r>0,                 R_(r+1)>R_r.                  (13)
```

THM-2866 similarly gives

```text
S_(r+1)>S_r,                                      (14)
```

so the slopes `(7)` are well-defined and positive.  THM-2873 proves, for
every `r>=1`,

```text
det [
 A_r      B_r      C_r
 A_(r+1)  B_(r+1)  C_(r+1)
 A_(r+2)  B_(r+2)  C_(r+2)
] >0.                                                 (15)
```

Dividing `(15)` by `A_rA_(r+1)A_(r+2)` and using `(6)--(7)` gives

```text
delta_r delta_(r+1)(lambda_(r+1)-lambda_r)>0.         (16)
```

The first two factors are positive by `(13)`, so `(9)` follows.
The weighted-chord lemma now proves `(5)`.

Equivalently, the normalized response points

```text
(R_r,S_r),                  r=1,2,3,...               (17)
```

form a strictly rising, strictly counterclockwise discrete arc: any
chord across an earlier block of rows has smaller slope than any chord
across a later disjoint block.  The proof uses only consecutive
curvature and positive horizontal increments; no new six-label
coefficient certificate is required.

## 3. Sharp failure boundaries

### 3.1. Row zero cannot be admitted uniformly

Take `W=d_4` and rows `0,1,2`.  Direct factorial evaluation gives

```text
[
 1    812    3854466
 5   5040   33243210
15  22260  201170970
].                                                     (18)
```

Its determinant and first-column-normalized determinant are

```text
-339570000,                    -4527600,               (19)
```

respectively.  Thus the cutoff `r_0>=1` is sharp even for a single
positive adjacent ray.

### 3.2. Cone positivity cannot be replaced by entry positivity

Take the signed polynomial

```text
W=d_0-(1/5)d_1.                                       (20)
```

On rows `1,2,3`, its exact response matrix is

```text
[
 3/5   6/5    351/125
 2/5  28/25    66/25
 1/5   4/5     33/25
].                                                     (21)
```

In particular, its first column is

```text
(A_1,A_2,A_3)=(3/5,2/5,1/5)>0,                       (22)
```

and in fact every displayed entry is positive.  Nevertheless

```text
D_(1,2,3)(W)=-1728/15625<0.                           (23)
```

The failure is therefore not a zero denominator or a negative moment
entry.  It occurs at the first logical input: signed adjacent
coefficients need not satisfy the positive response ladders
`(13)--(16)`.

## 4. Exact companion and finite gap scout

The exact companion reconstructs every tensor from

```text
L(prod_i d_(a_i))
 =1/(prod_i(a_i+1)!)
  sum_e (-1)^e e_e(a_i+1)
       (sum_i(a_i+1)-e)!,                             (24)
```

then:

1. verifies the weighted-chord identity `(11)` exactly;
2. reconstructs both matrices `(18)` and `(21)`;
3. checks both determinants and the normalized row-zero curvature;
4. checks `17,640` theorem-scope cells with
   `0<=p<=6`, `h in {1,2}`, integer `0<=x,y<=3` not both zero,
   and every row triple in `{1,...,9}`; and
5. separately checks `52,920` cells with `3<=h<=8` in the same finite
   box.

All `70,560` scanned determinants are positive; the minimum raw value is
`12`, attained at `p=0,h=1,(x,y)=(1,0)` on rows `1,2,3`.  The
`h=3,...,8` scan is **FINITE-EXACT evidence only**.  It neither proves
an arbitrary-gap theorem nor enters the proof of `(5)`.

Reproduce with

```text
python3 04-computation/gmc_nonconsecutive_response_row_tp3_thm2909.py
python3 -O 04-computation/gmc_nonconsecutive_response_row_tp3_thm2909.py
```

Both executions byte-match

```text
05-knowledge/results/gmc_nonconsecutive_response_row_tp3_thm2909.out.
```

The LF-normalized hashes are recorded in the frontmatter.

## 5. Exact scope

This theorem upgrades the **response-row choice**, not the support
architecture.  It proves `(5)` only for the two-ray gaps `h=1,2`
already covered by THM-2873.  It does not prove arbitrary-gap TP3,
universal TP3 for an arbitrary nonnegative cone, or any signed-cone
version.  Here “TP3 arc” refers to the ordered maximal three-column
minors `(5)`; it does not silently assert every `2 x 2` minor required
by a full matrix-total-positivity convention.

It is also a one-curve statement.  THM-2906 gives two positive response
curves with strict local TP3 but opposite signs for the mixed
endpoint-holonomy statistic.  Therefore `(5)` supplies no comparison
between different `W`, no cubic-ideal selector, and no new closure of
the general shared cubic--quartic multipole line, SFC, GMC, or a Gaussian
nullcone.

An independent hostile audit rederived `(11)` directly, checked that
THM-2830 and THM-2873 are used in the required directions, reconstructed
both hostile matrices without importing the companion's tensor routine,
and replayed normal, optimized and stored output with the declared LF
hashes.

**QED.**
