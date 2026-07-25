---
id: THM-2205
title: "Scalar depth-(1,1,3) exact lift-capacity exclusion"
status: >
  PROVED + VERIFIED-EXACT. In the scalar five-unit/three-blocker branch, the
  actual blocker valuation profile (1,1,3) is empty. On the primitive
  13^4 guard-safe annulus, the two depth-one blocker masks are pullbacks
  from the 1,014 unit sign classes modulo 13^3. An exact cyclic root-fibre
  matrix records, for every one of the 13,182 unit sign classes modulo
  13^4 and every quotient phase, its zero, one, or two guard-surviving
  endpoints. Exact transpose/intersection accumulation audits all 514,605
  unordered shallow pairs. The unique minimum conditional top-five margin
  is 1,286. The previous depth-three audit is recovered as a positive
  control with margin 86. Together with THM-2204 this closes the repeated
  shallower-depth profiles (1,1,3) and (2,2,3), but the mixed profile
  (1,2,3) and every profile with deepest valuation at least four remain
  open; this is not a proof of LRC(14).
source: codex-2026-07-24-scalar-depth113-exact-lift
depends_on:
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
related:
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2204-scalar-depth-223-thirteen-lift-capacity-law
script: 04-computation/lrc14_scalar_five_plus_three_depth_three_lift_exact_codex_20260724.py
output: 05-knowledge/results/lrc14_scalar_five_plus_three_depth_three_lift_exact_codex_20260724.out
script_sha256: d885ecf2557fe8832b4f2f198be202dc7a40fc26c00a350818880fba28726ba6
output_sha256: 2308347d6816677b7fa6c459c2b8fa60ed78356172da7e4f20b356215ed119f8
hash_basis: working-tree bytes (LF)
---

# THM-2205 -- scalar depth-`(1,1,3)` exclusion

Put

```text
D_a={t in R/Z:||at||<1/14},
C_H={t in R/Z:||Ht||>1/7}.                           (1)
```

In the scalar `5+3` branch of THM-2192 and THM-2198, suppose

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(j=1)^3 D_(c_j)             (2)
```

almost everywhere.  The coefficients `H,q_1,...,q_5` are positive
thirteen-units and the three actual blockers are positive multiples of
thirteen.  After relabelling, THM-2192 gives

```text
1<=lambda_1<=lambda_2<lambda_3,
lambda_j=nu_13(c_j).                                 (3)
```

This theorem excludes

```text
(lambda_1,lambda_2,lambda_3)=(1,1,3).                (4)
```

The proof is a finite exact lift argument, but its finite object is not a
flat mask table.  It is the labelled root-fibre incidence matrix singled
out abstractly by THM-2201 and by THM-2204's guard-hole correlation vector.

## 1. The primitive depth-four layer

Assume (4) and put

```text
N=13^4=28561,                    Q=13^3=2197.         (5)
```

Multiplication of primitive numerators by `H modulo N` normalizes the guard
to one and replaces every terminal coefficient by its product with
`H^(-1) modulo N`.  This is a bijection and preserves all three
thirteen-valuations.

The primitive guard-safe universe is

```text
U_N={z mod N:13 does not divide z and 7||z||_N>N},
|U_N|=18830,                                         (6)
```

where `||z||_N=min(z mod N,-z mod N)`.

The depth-three blocker is safe throughout `U_N`.  Indeed, if its normalized
coefficient is `13^3w`, then at a primitive numerator

```text
13^3w*z/13^4=wz/13
```

is a nonzero thirteenth root.  Its norm is at least
`1/13>1/14`.

Each depth-one blocker has the form `13u`.  Its mask depends only on the
unit part `u modulo Q`, up to sign, so there are

```text
phi(Q)/2=1014                                       (7)
```

possible shallow masks.  Distinct actual blockers may reduce to the same
mask.  Thus the exact pair universe includes repetitions and has

```text
C(1014+1,2)=514605                                  (8)
```

unordered pairs.  The five unit masks range through

```text
phi(N)/2=13182                                      (9)
```

sign classes modulo `N`.

## 2. Cyclic exponent coordinates

The integer `2` is a primitive root modulo every power of thirteen.  At
the two depths in (5), define

```text
G_4=(Z/NZ)^*/{+/-1}=<2>,       |G_4|=13182,
G_3=(Z/QZ)^*/{+/-1}=<2>,       |G_3|=1014.           (10)
```

We identify these cyclic groups with exponent classes

```text
G_4=Z/13182Z,                  G_3=Z/1014Z.          (11)
```

Reduction modulo `Q` is exponent reduction modulo `1014`.  Consequently
the fibre over `r in G_3` is

```text
F_r={r+j*1014:j in F_13} subset G_4.                (12)
```

All guard and terminal predicates are invariant under sign.  Counting in
the quotient therefore counts exactly one point from every `{z,-z}` pair;
all final cardinalities are doubled.

Define the guard and unit-mask predicates

```text
C(e)=1_[7||2^e||_N>N],
M(q,e)=1_[14||2^(q+e)||_N<N].                       (13)
```

The root-sheet matrix is

```text
W(q,r)=sum_(e in F_r) C(e)M(q,e),
q in G_4, r in G_3.                                 (14)
```

THM-2198's exact integer-window law gives

```text
W(q,r) in {0,1,2}.                                  (15)
```

It is the number of endpoints of the labelled unit mask `q` which survive
the guard on the rooted fibre over `r`.  The fibre height

```text
h(r)=sum_(e in F_r) C(e)                            (16)
```

belongs to `{9,10}`.  On the sign quotient the exact histogram is

```text
725 fibres of height 9,
289 fibres of height 10,
725*9+289*10=9415.                                  (17)
```

Doubling (17) recovers (6).

For a shallow exponent `a in G_3`, put

```text
A_a(r)=1_[14||2^(a+r)||_Q<Q].                       (18)
```

The mask of the actual coefficient `13*2^a` is constant on `F_r` and is
active precisely when `A_a(r)=1`.

## 3. The exact transpose/intersection accumulation

For a shallow pair `(a,b)`, let

```text
P_(a,b)={r in G_3:A_a(r)=A_b(r)=0}.                 (19)
```

The sign-quotient residual size and the conditional capacity of a unit
class `q` are

```text
R(a,b)=sum_(r in P_(a,b)) h(r),
C_q(a,b)=sum_(r in P_(a,b)) W(q,r).                 (20)
```

These are not approximations to the torsion masks: equations (12)--(18)
give an exact disintegration of those masks.

The computation evaluates (20) by a reusable integer dynamic program.
Put

```text
U_q=sum_r W(q,r),
L_(a,q)=sum_r A_a(r)W(q,r),
I_(a,b)={r:A_a(r)=A_b(r)=1}.                         (21)
```

Then exact inclusion--exclusion gives

```text
C_q(a,b)
 =U_q-L_(a,q)-L_(b,q)+sum_(r in I_(a,b)) W(q,r).    (22)
```

There is an identical formula for `R(a,b)` with `W(q,r)` replaced by
`h(r)`.

Equation (22) is the structural reason the large audit is feasible and
faithful.  The matrix `L` is the integer product of the shallow-incidence
matrix with `W^T`.  For each pair, only the intersection rows
`I_(a,b)` must be restored.  Thus the operation is a
**transpose/intersection DP**:

```text
root-sheet matrix W
 -> all one-shallow losses L
 -> restore the exact rows deleted twice
 -> retain the complete 13182-vector of unit capacities. (23)
```

Every accumulator is a nonnegative integer.  The companion uses only
`uint8`, `uint16`, and `int32`, after checking the corresponding maximum
possible sums; it uses no floating point, FFT, or optimization-sensitive
assertion.

This is precisely the sidecar absent from a scalar lift average.
THM-2204 proves that the sum of the thirteen coefficient-lift capacities
has a simple recurrence, but also exhibits a family whose capacities range
from zero to `2460`.  Equation (22) keeps the labelled guard-hole
correlations which that family sum destroys.  THM-2201's full labelled
Hasse jets can reconstruct the same phasewise incidence; here the raw
integer augmentation is sufficient because each entry in (15) is at most
two.

## 4. Exact capacity deficit

For each shallow pair, order the `13182` integers in (20) decreasingly:

```text
C_(1)(a,b)>=...>=C_(13182)(a,b).                    (24)
```

The exhaustive result is

```text
R(a,b)-sum_(i=1)^5 C_(i)(a,b)>=643                  (25)
```

on the sign quotient, for all `514605` pairs.  Equivalently, on the full
primitive annulus,

```text
|R_full(a,b)|-sum_(i=1)^5 C_(i,full)(a,b)>=1286.    (26)
```

The minimum is unique in the fixed exponent convention.  Its exponent and
least-positive sign representatives are

```text
(a,b)=(30,599),                 labels=(183,799).    (27)
```

At this row,

```text
|R_full|=13946,

((capacity,unit label))_(i=1)^5
 =((2700,2380),(2652,5193),(2506,10386),
   (2434,1190),(2368,7139)).                         (28)
```

Their sum is `12660`, and

```text
13946-12660=1286.                                   (29)
```

To freeze both the carrier and the scan, the exact digests are

```text
root-sheet matrix W:
  6c3a3cf58a37cb1e728864487e7fca97c5772d39a2fbace2b0df56df1fa75d97,

one-shallow loss matrix L:
  dd37b559a3b811df8e17bc3cb1238b27ddabb36b50fb0e227c716d5baf40cde0,

all-pair summary table:
  7c05788531e385a83c641c49e119a13e1a301fe1c17dd8bcf394e45819e5b186.
                                                               (30)
```

The companion independently reconstructs (27)--(29) from the definition
on all `18830` full torsion residues.  Before the new scan, it also reruns
the `N=13^3` universe and recovers THM-2198's unique hostile row

```text
shallow labels=(14,46),
residual=1046,
top five=(204,200,190,186,180),
margin=86.                                           (31)
```

Thus the old theorem is an inherited positive control, not merely a shared
code path.

## 5. Exclusion of `(1,1,3)`

Reduce the five actual unit coefficients modulo `N` and sign.  Repeated
classes do not enlarge their union, so discard repetitions.  Inside the
residual belonging to the actual shallow pair, the union of the remaining
at most five masks has cardinality at most the sum of their individual
capacities.  This is at most the sum of the five largest entries in (24),
which is strictly below the residual size by (26).

Hence some primitive residue is safe from all five unit masks and both
depth-one blockers.  It is safe from the depth-three blocker by Section 1
and lies in the guard-danger set by construction.  A power of thirteen is
coprime to both seven and fourteen, so no guard or terminal equality occurs
at this residue.  It thickens to an uncovered open interval, contradicting
the almost-everywhere cover (2).  This proves that (4) is empty.

## 6. Recursive signal and scope

The exact minimizers expose a nonaccidental lift:

```text
depth-three control:
  worst shallow pair (14,46);
  first two extremal unit labels (183,799);

depth-four theorem:
  worst shallow pair (183,799).                      (32)
```

Thus the next hostile shallow pair is obtained by retaining the two largest
unit-capacity labels from the preceding hostile row.  Its new extremal labels
begin `(2380,5193,...)`.  This is a precise recursive signal, but it is not
yet an all-depth theorem: THM-2204 shows that scalar lift-family sums do not
control the lifted top-five order statistic, and the size of the labelled
correlation carrier grows with depth.

Together, THM-2198, THM-2204, and the present theorem exclude

```text
(1,1,2),       (2,2,3),       (1,1,3).               (33)
```

The mixed profile `(1,2,3)` is not covered by either repeated-depth audit.
All profiles with `lambda_3>=4` also remain open.  The next exact object is
the finite lift transition on the labelled correlation/Hasse-jet carrier,
not its scalar family sum.  This theorem does not prove LRC(14).  QED.
