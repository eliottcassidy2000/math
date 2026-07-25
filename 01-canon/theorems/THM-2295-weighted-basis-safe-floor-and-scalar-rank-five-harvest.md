---
id: THM-2295
title: "Weighted basis-safe floor and scalar rank-five harvest"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. On a nine-coordinate
  relation torus with one guard-safe coordinate of mass 5/7 and eight
  ordinary safe coordinates of mass 6/7, every saturated rank-r relation
  lattice with nonzero coordinate characters has safe Haar mass at least
  ((5-r)/7)(6/7)^(8-r) for 0<=r<=4. A Machin-certified squared-Fejer
  comparison at degree 196 then proves that every live scalar
  five-unit/three-blocker cover has relation rank at least five by
  coefficient height 196. The fixed original-row section has rank at least
  five by height 392. No scalar profile is excluded, no exact Fourier atom
  or owner ancestry is selected, and LRC(14) remains open.
source: codex-2026-07-25-weighted-basis-rank-five
depends_on:
  - THM-2185-rank-two-safe-cube-floor-and-height-500-rank-three-harvest
  - THM-2190-basis-safe-floor-and-height-500-rank-six-harvest
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2283-mixed-rank-two-safe-torus-floor-and-scalar-rank-three-harvest
related:
  - THM-2282-thirteen-adic-saturation-and-unit-anchored-minor
  - THM-2284-thirteen-adic-anchored-rank-three-plucker-lift
script: 04-computation/lrc14_weighted_basis_rank_five_thm2295.py
output: 05-knowledge/results/lrc14_weighted_basis_rank_five_thm2295.out
script_sha256: 63cb31b4c5102b198f9140ad571a72ffc270bbc770741903e57d58e86e3660dc
output_sha256: b6e1490db811e5bb926137cb2f406f49a23d2ef24a89966928240422e7e954b4
hash_basis: working-tree bytes (LF)
---

# THM-2295 -- weighted basis-safe floor and scalar rank five

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2283 proves a positive safe-mass floor for every saturated rank-two
mixed scalar torus by classifying all sparse relation patterns. An older,
less-used mechanism from THM-2190 gives a stronger rank ladder here:
choose coordinate characters as a rational basis, retain their safe product
exactly, and union-bound only the remaining relation-rank many coordinates.

The guard has a different safe mass from the other eight coordinates, but
the two possible locations of the guard give the same formula:

```text
alpha_r=((5-r)/7)(6/7)^(8-r),          0<=r<=4.       (1)
```

At rank four this floor barely clears `1/13`. The exact Jackson error at
degree `196` is smaller still, forcing a fifth bounded relation.

## 1. The mixed scalar torus

Let

```text
w=(w_0,w_1,...,w_8) in Z^9,
w_i!=0 for every i,                                  (2)

Lambda(w)={m in Z^9:m.w=0}.                          (3)
```

Coordinate zero is the guard. Write its safe and danger sets as

```text
J_0, I_0,
measure(J_0)=5/7,
measure(I_0)=2/7.                                    (4)
```

The other eight coordinates have safe and danger sets

```text
J_i, I_i,
measure(J_i)=6/7,
measure(I_i)=1/7,                    1<=i<=8.         (5)
```

These are the translated circular intervals in the scalar
five-unit/three-blocker reduction. Put

```text
Safe_mix=J_0 x J_1 x ... x J_8.                      (6)
```

Let `L` be a saturated rank-`r` sublattice of `Lambda(w)` and define

```text
K_L={x in (R/Z)^9:l.x=0 mod 1 for every l in L}.     (7)
```

Saturation makes `K_L` connected and

```text
ann(K_L)=L.                                          (8)
```

Every coordinate character is nonzero on `K_L`. Indeed,

```text
e_i in L subset Lambda(w)
 -> w_i=0,                                           (9)
```

contrary to (2).

## 2. Weighted coordinate-basis floor

Put

```text
d=9-r,
X=Z^9/L.                                             (10)
```

The coordinate classes generate the rank-`d` free abelian group `X`.
Choose `d` of them as a basis of `X tensor Q`. Exactly as in THM-2190,
their finite-index span gives a finite cover

```text
rho:(R/Z)^d -> K_L                                   (11)
```

that pushes Haar probability to Haar probability. On this cover, the
chosen coordinate characters are

```text
D y_1,...,D y_d                                      (12)
```

for one nonzero integer `D`, while each extra coordinate pulls back to a
nonzero character

```text
a_i.y,                a_i in Z^d minus {0}.           (13)
```

Let `B` be the event that every basis coordinate is safe. If an extra
coordinate `i` has danger mass `p_i`, choose any `j` for which
`(a_i)_j!=0`. Integrate first in `y_j`. Multiplication by the nonzero
integer `(a_i)_j` pushes Haar to Haar, so the danger fibre has total mass
`p_i`; intersecting it with the safe condition on `D y_j` can only
decrease it. Therefore

```text
measure({a_i.y in I_i} intersection B)
 <=p_i product_(h!=j) q_h,                           (14)
```

where `q_h` are the safe masses of the basis coordinates.

There are two cases.

### 2.1 The guard is in the basis

The basis-safe mass is

```text
measure(B)=(5/7)(6/7)^(d-1).                         (15)
```

All `r` extra coordinates are ordinary. In (14), using the guard basis
coordinate costs

```text
(1/7)(6/7)^(d-1),                                    (16)
```

and using an ordinary coordinate costs less. The union bound gives

```text
measure_(K_L)(Safe_mix)
 >=(5/7)(6/7)^(d-1)-r(1/7)(6/7)^(d-1)
 =((5-r)/7)(6/7)^(8-r).                              (17)
```

### 2.2 The guard is an extra coordinate

Now all basis coordinates are ordinary, so

```text
measure(B)=(6/7)^d.                                  (18)
```

The extra guard costs at most

```text
(2/7)(6/7)^(d-1),                                    (19)
```

and the other `r-1` ordinary extras cost at most

```text
((r-1)/7)(6/7)^(d-1).                                (20)
```

Again,

```text
measure_(K_L)(Safe_mix)
 >=(6/7)^d-((r+1)/7)(6/7)^(d-1)
 =((5-r)/7)(6/7)^(8-r).                              (21)
```

This proves (1). The exact floors are

```text
r=0: 8398080/40353607,
r=1: 1119744/5764801,
r=2: 139968/823543,
r=3: 15552/117649,
r=4: 1296/16807.                                    (22)
```

They strictly decrease. In particular,

```text
alpha_r>=alpha_4=1296/16807,

alpha_4-1/13=41/218491>0.                            (23)
```

The proof uses no positivity or distinctness of the scalar values beyond
their nonvanishing in (2).

## 3. The degree-196 Jackson ledger

Use the normalized squared-Fejer approximants from THM-2185 and THM-2283.
For a bandwidth `N`, each safe interval has an approximant `q_(i,N)` with

```text
0<=q_(i,N)<=1,
degree(q_(i,N))<=2N-2,
||q_(i,N)-1_(J_i)||_1<=eta_N.                        (24)
```

THM-2283's Machin-series certificate gives the strict rational upper cap
on `pi`

```text
pi<104348/33215.                                     (25)
```

Substitution into the exact odd-mode formula for (24), followed by exact
rational summation, gives

```text
N=98:
  eta_bar_98>43/10000>1/234,

N=99:
  eta_99<eta_bar_99<213/50000<1/234.                 (26)
```

The first displayed cap does not prove the true error at `N=98` is too
large. It says only that this exact rational-upper certificate does not
clear the rank-four ledger there.

At `N=99`,

```text
H=2N-2=196.                                          (27)
```

Two nine-factor product telescopes cost at most `18 eta_99`. The convenient
strict margin is

```text
alpha_4-18(213/50000)
 =180981/420175000
 >0.                                                 (28)
```

For comparison, the adjacent displayed certificate fails:

```text
alpha_4-18(43/10000)
 =-24309/84035000
 <0.                                                 (29)
```

The companion independently reconstructs every closed Jackson coefficient
at `N=98,99` by convolution, scans all `N=2,...,99`, and verifies that `99`
is the first passing bandwidth for this exact rational cap and ledger. For
the exact cap, rather than the convenient estimate in (28), it also checks

```text
481/1000000
 <alpha_4-18 eta_bar_99
 <482/1000000,
```

so the exact margin is greater than `1/2500`. This is certificate
optimality, not an optimality theorem over all kernels or exact values of
`pi`.

## 4. Scalar rank five

For `H>=0`, put

```text
W_H^*(w)
 =span_Q(Lambda(w) intersection [-H,H]^9).           (30)
```

Assume the scalar mixed safe event is null:

```text
measure{t in R/Z:
  w_0t in J_0 and w_it in J_i for 1<=i<=8}=0.        (31)
```

Suppose for contradiction that

```text
r=dim_Q W_196^*(w)<=4.                               (32)
```

Take the full bounded-relation lattice

```text
L=W_196^*(w) intersection Z^9.                       (33)
```

This choice is load-bearing. It is saturated, has rank `r`, and lies in
`Lambda(w)`.

Let

```text
Q_99(x)=product_(i=0)^8 q_(i,99)(x_i).               (34)
```

Every Fourier frequency in (34) belongs to `[-196,196]^9`. Its line average
survives exactly when

```text
m.w=0,
```

which is exactly when

```text
m in Lambda(w) intersection [-196,196]^9
 subset L.                                           (35)
```

Conversely, every frequency in `L` annihilates both the line and `K_L`.
Since `ann(K_L)=L`, the complete finite survivor sets agree:

```text
integral_(R/Z) Q_99(wt)dt
 =integral_(K_L) Q_99(x)dx.                          (36)
```

Every coordinate marginal is Haar on both averaging groups. On the line
this follows from `w_i!=0`; on `K_L` it follows from (9). Equations (24)
and (31), with a nine-factor telescope, give

```text
integral_(R/Z) Q_99(wt)dt<=9 eta_99.                 (37)
```

The weighted basis floor and the torus telescope give

```text
integral_(K_L) Q_99(x)dx
 >=alpha_r-9 eta_99
 >=alpha_4-9 eta_99.                                 (38)
```

Equations (36)--(38) would imply

```text
alpha_4<=18 eta_99,                                  (39)
```

contradicting (28). Therefore

```text
dim_Q W_196^*(w)>=5.                                 (40)
```

This abstract statement applies to every nonzero nine-coordinate integer
row with the scalar mixed null-cover geometry. In particular it applies to
all `165` live scalar first-depth-one profiles.

## 5. Fixed-section original-row lift

THM-2203's fixed scalar section lifts a scalar relation by

```text
(x_H,x_rest) |->(2x_H,x_rest).                       (41)
```

The map is integral and injective, so it preserves five-dimensional
relation rank and doubles coefficient height by at most two. Equation (40)
therefore gives

```text
fixed-section original relation rank>=5
by coefficient height 392.                           (42)
```

## 6. Comparison and exact stopping boundary

This theorem strictly strengthens the relation-rank conclusion of THM-2283:

```text
THM-2283: scalar rank>=3 by height 3540;
THM-2295: scalar rank>=5 by height 196.               (43)
```

THM-2283 remains independently valuable for its sharp rank-two mixed-torus
floor, sparse-line classification, delicate mean-one core, and exact
Jackson machinery used in Section 3.

The elementary weighted basis floor stops exactly after rank four:

```text
alpha_5=0.                                           (44)
```

At relation rank five, the six basis coordinates no longer have enough
safe budget to pay all five extras by this union bound. A strict
facet-intersection theorem, analogous to the rank-six boundary of
THM-2190 and its uniform repair in THM-2193, is needed to force rank six.

High relation rank alone does not select support, coefficient signs, a
thirteen-adic pivot, an exact Fourier lift, or owner ancestry. The
alternating-matroid/large-cogirth models remain hostile controls: many
independent relations can all have large support.

The connection and loss ledger is

```text
source:
  THM-2190's coordinate-basis finite cover, the scalar mixed interval
  masses, and THM-2283's exact Jackson cap;

target:
  a uniform mixed safe floor through relation rank four and five
  independent scalar relations at height 196;

map:
  retain a product-safe coordinate basis, union-bound only the r
  extra characters, and compare the full bounded survivor lattice
  on the line and its relation torus;

preserved:
  all nine scalar labels, guard/ordinary interval type, the complete
  degree-196 Fourier survivor set, rational rank, and fixed-section lift;

destroyed:
  the individual relation basis, support, coefficient signs, modular
  pivot labels, exact frequency, root sheet, and owner ancestry;

cheapest hostile probes:
  guard inside versus outside the coordinate basis, rank four at the
  41/218491 coarse margin, the N=98 certificate failure, and a
  high-cogirth rank-five code;

needed sidecar:
  a positive rank-five mixed facet floor or an owner-labelled operation
  that exploits five relations without scalarizing their support.       (45)
```

No scalar profile is excluded. LRC(14) remains open.

## 7. Exact reproduction

Run

```bash
python3 04-computation/lrc14_weighted_basis_rank_five_thm2295.py
python3 -O 04-computation/lrc14_weighted_basis_rank_five_thm2295.py
```

Both executions must match

```text
05-knowledge/results/lrc14_weighted_basis_rank_five_thm2295.out
```

byte-for-byte after LF normalization. The companion uses only explicit
integer and `Fraction` arithmetic; optimized execution retains every
validity check. QED.
