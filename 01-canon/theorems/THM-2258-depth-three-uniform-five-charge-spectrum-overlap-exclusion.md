---
id: THM-2258
title: "Depth-three uniform five-charge spectrum/overlap exclusion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the scalar
  five-unit/three-blocker branch, THM-2243's uniform Bellman bound on the
  sole depth-three survivor (3,4,5) forces the five individual
  guard/danger charges to have sum at least 5092517/8911032. THM-2080's
  exact odd-guard spectrum then forces all five distinct reduced ratios
  into a 13-row bank, of whose 1287 five-subsets only 15 meet that sum.
  Exact literal union sweeps give residual at least 317/1155 in every row,
  exceeding the Bellman bound 8907541/62377224 by at least
  150561431/1143582440. Equivalently the forced near-extremal marginals
  incur overlap rebate at least 29/220. This independently excludes
  (3,4,5), reproducing the scalar-ledger reduction 166 to 165 already
  obtained by THM-2257; it is not an additional decrement and does not
  prove LRC(14).
source: codex-2026-07-25-uniform-five-charge-spectrum
depends_on:
  - THM-2080-unequal-comb-overlap-removes-depth-five
  - THM-2243-composite-union-transfer-dual-depth-three-profile-exclusion
related:
  - THM-2091-centered-triple-energy-obstruction
  - THM-2137-deep-scalar-tail-boundary-complexity
  - THM-2250-depth-three-pair-incidence-partition-reduction
  - THM-2257-depth-three-common-core-169-image-sieve-exclusion
script: 04-computation/lrc14_depth_three_uniform_five_charge_overlap_thm2258.py
output: 05-knowledge/results/lrc14_depth_three_uniform_five_charge_overlap_thm2258.out
script_sha256: 521ac39d8e2636c6636ab70ea7f5f600cbf828d65819b7c3cd1144511b2df548
output_sha256: 204837e97fc3275c2c2aa8cf92fbe535c2ebf219fd897973078141bf033e7702
hash_basis: working-tree bytes (LF)
---

# THM-2258 -- near-maximal five-charge marginals force catastrophic overlap

On `T=R/Z`, put

```text
D_q={x:||qx||<1/14},
C_H={x:||Hx||>1/7},
E_H=T minus C_H={x:||Hx||<=1/7} a.e.                 (1)
```

Endpoints are null throughout, so the displayed weak/strict convention for
`E_H` does not change a measure. Assume for contradiction that a surviving
scalar five-unit/three-blocker row has valuation profile `(3,4,5)`. Thus
`H,q_1,...,q_5` are positive thirteen-units, `H` is odd, and the five `q_i`
are distinct.

Define the literal five-unit residual

```text
A=C_H minus union_(i=1)^5 D_(q_i),       p=measure(A). (2)
```

Under this survivor hypothesis, the argument below proves

```text
p>=317/1155
  =8907541/62377224+150561431/1143582440
  >8907541/62377224.                                  (3)
```

THM-2243 proves that a surviving `(3,4,5)` row would instead satisfy

```text
p<=B:=8907541/62377224.                               (4)
```

Equations (3) and (4) are contradictory. Hence the scalar depth-three
profile `(3,4,5)` is empty.

This proof is uniform in all three normalized blocker cores. It does not use
THM-2250's subsequent equality reduction or THM-2257's common-core image
sieve.

## 1. Turn the Bellman ceiling into a marginal-charge floor

For each unit danger define its charge inside the guard complement,

```text
A_i=C_H intersect D_(q_i),             c_i=measure(A_i). (5)
```

The positive part of THM-2243's signed residual is exactly `1_A`: inside
`C_H` it equals one precisely when no unit danger is active, and outside
`C_H` it vanishes. Thus its `p` is literally (2), not merely a score
majorizing (2).

The elementary union bound gives

```text
p=5/7-measure(union_i A_i)
 >=5/7-sum_i c_i.                                     (6)
```

Combining (4) and (6), any survivor must pay

```text
sum_i c_i
 >=S:=5/7-B
    =5092517/8911032.                                 (7)
```

This implication is only necessary: it discards every higher intersection.
The discarded overlap will be the decisive coordinate.

## 2. The exact reduced-ratio charge spectrum

For each `i`, reduce the pair `(H,q_i)`:

```text
g_i=gcd(H,q_i),       a_i=H/g_i,       b_i=q_i/g_i.  (8)
```

Then

```text
gcd(a_i,b_i)=1,       a_i is odd,       13 does not divide a_i b_i. (9)
```

The reduced pairs are distinct because `b_i/a_i=q_i/H` and the `q_i` are
distinct.

THM-2080 gives, with

```text
x=(a mod 14)/14,                 y=(b mod 7)/7,
F(x,y)=min(x,y)+(x+y-1)_+-2xy,                       (10)
```

the exact charge

```text
c(a,b)
 =1/7-measure(E_a intersect D_b)
 =5/49-(2/(ab))F(x,y).                               (11)
```

Since `F>=-1/8`,

```text
c(a,b)<=5/49+1/(4ab).                                (12)
```

The four largest capacities at distinct admissible reduced ratios are

```text
capacity       (a,b)

5/42           (1,6)
9/77           (11,1)
4/35           (1,5)
4/35           (3,5).                                (13)
```

Here is the finite-completeness argument behind (13), rather than a
hardcoded appeal to the displayed rows. A capacity at least `4/35` obeys,
by (12),

```text
ab<=1/[4(4/35-5/49)]=245/12<21.                     (14)
```

The companion evaluates (11) on every coprime pair with `a` odd,
`13 does not divide ab`, and `ab<=20`; there are 35. Its first four rows
are exactly (13). Thus the sum of the four largest possible distinct-ratio
charges is

```text
T_4=5/42+9/77+4/35+4/35=1073/2310.                  (15)
```

Order the five actual charges decreasingly. Equations (7) and (15) force
even the smallest one to be at least

```text
kappa=S-T_4
 =122343163/1143582440.                              (16)
```

In particular, all five charges are at least `kappa`. Since

```text
1/[4(kappa-5/49)]
 =2001269270/39557541
 =50.591346...<51,                                  (17)
```

equation (12) forces

```text
a_i b_i<=50                                  for every i. (18)
```

This is the spectrum collapse: an a priori infinite five-ratio problem has
become a small exact bank.

## 3. The 13 ratios and 15 marginally viable five-sets

Evaluating (11) on the complete universe

```text
a odd, gcd(a,b)=1, 13 does not divide ab, ab<=50     (19)
```

and retaining `c(a,b)>=kappa` gives exactly

```text
capacity       reduced ratios (a,b)

5/42           (1,6)
9/77           (11,1)
4/35           (1,5), (3,5)
1/9            (9,1), (9,2)
17/154         (11,2)
19/175         (25,1)
3/28           (1,4), (1,12), (1,20), (3,4), (5,4). (20)
```

There are 102 pairs in (19), 13 survivors in (20), and

```text
binomial(13,5)=1287                                  (21)
```

five-subsets. Exact addition of the five capacities leaves only 15 subsets
whose sum is at least `S`.

The following table is the complete normalized atlas. `Delta` is
`sum_i c_i-S`. The column `Omega` is the overlap rebate

```text
Omega=sum_i c_i-measure(union_i A_i)
      =integral_(C_H) (number of active D_(q_i)-1)_+. (22)
```

```text
 H      normalized q_i                    Delta                 p          Omega

 33     3,55,132,165,198                  183527/1143582440     131/462    31/220
 33     3,55,165,198,396                  183527/1143582440     733/2310   269/1540
 33     3,55,165,198,660                  183527/1143582440     661/2310   221/1540
 33     3,44,55,165,198                   183527/1143582440     317/1155   29/220
165     15,132,275,825,990                183527/1143582440     1121/3850  3431/23100
 99     9,11,165,495,594                  42493973/10292241960  6103/20790 46/297
 99     9,22,165,495,594                  42493973/10292241960  5977/20790 221/1485
 33     3,6,55,165,198                    3896457/1143582440    237/770    389/2310
825     33,75,1375,4125,4950              9086081/5717912200   2407/8250  207/1375
 99     9,11,22,495,594                   9820189/10292241960   3041/10395 3133/20790
 99     9,11,18,495,594                   342047/1470320280     3092/10395 46/297
 99     9,18,22,495,594                   342047/1470320280     2047/6930  353/2310
 99     9,11,22,165,594                   9820189/10292241960   3203/10395 3457/20790
 99     9,11,18,165,594                   342047/1470320280     6451/20790 317/1890
 99     9,18,22,165,594                   342047/1470320280     2113/6930  25/154
                                                               (23)
```

The smallest literal residual and the smallest overlap rebate occur in the
fourth row:

```text
p_min=317/1155,                    Omega_min=29/220. (24)
```

This is the structural point hidden by (6). Marginal near-extremality does
not make the five danger pieces efficient. The arithmetic ratios that make
all five individual charges large force a duplicate-cover debt of at least
`29/220`. Indeed, from (22),

```text
p=5/7-sum_i c_i+Omega.                              (25)
```

For the worst row the marginal surplus is only
`183527/1143582440`, but the overlap rebate is `29/220`; their difference
is exactly

```text
29/220-183527/1143582440
 =150561431/1143582440
 =317/1155-B.                                        (26)
```

Thus even the row most favorable to a small literal residual misses the
Bellman ceiling by the stated positive margin.

## 4. Why the normalized atlas is uniform in the original scale

Fix one of the 15 reduced-ratio sets and put

```text
L=lcm(a_1,...,a_5).                                  (27)
```

Every `a_i` divides the original `H`, so `d=H/L` is a positive integer.
Equation (8) gives

```text
H=dL,                 q_i=d(Lb_i/a_i).               (28)
```

The normalized atlas therefore uses

```text
H_0=L,                q_(i,0)=Lb_i/a_i.              (29)
```

Multiplication by `d` is a Haar-preserving surjective endomorphism of the
circle, and

```text
C_H=[d]^(-1)C_(H_0),       D_(q_i)=[d]^(-1)D_(q_(i,0)). (30)
```

Hence the literal residual in (2) has exactly the normalized measure in
(23). No coprimality of `d` with the normalized coefficients is needed.
This proves that the 15-row computation covers every original fixed-`H`
configuration, not just primitive examples.

## 5. Exact verification and independent path

The companion uses rational arithmetic only.

1. **Formula path.** It evaluates (10)--(12), proves the global top-four
   reduction (14), checks all 102 rows of (19), constructs the 13-row bank,
   and tests all 1287 five-subsets against (7).
2. **Literal interval path.** For a comb of speed `s`, its only relevant
   endpoints are

   ```text
   (14j+-2)/(14s) for C_s,
   (14j+-1)/(14s) for D_s.                           (31)
   ```

   The script forms the common rational endpoint partition and evaluates a
   midpoint of every atom. It independently checks (11) on all 102
   small-product rows, then computes the literal five-comb residual for all
   15 normalized rows.

Normal and `python3 -O` runs are byte-identical. An independent common-mesh
implementation, not importing the companion, reproduced all 13 candidates,
15 qualifying subsets, and the worst count

```text
7608/27720=317/1155.                                 (32)
```

Reproduction:

```bash
python3 04-computation/lrc14_depth_three_uniform_five_charge_overlap_thm2258.py
python3 -O 04-computation/lrc14_depth_three_uniform_five_charge_overlap_thm2258.py
```

## 6. Consequence and nonconsequence

Combining (3) with THM-2243's necessary bound (4) proves that `(3,4,5)` is
empty. Starting from THM-2239's pre-closure scalar ledger of 166, this is an
independent proof of the reduction to 165 first-depth-one rows.

THM-2257 already made that decrement in the current canonical proof graph by
a different common-core image argument. The present theorem is therefore an
independent proof and a mechanism refinement, not a second subtraction.
The 165 first-depth-one scalar rows and the routed non-scalar frontiers
remain; LRC(14) is not proved.

The reusable pattern is:

```text
near-extremal scalar upper bound
 -> rank the individual edge/charge spectrum
 -> force every edge into a finite high-charge atlas
 -> restore the literal union geometry
 -> charge the unavoidable duplicate-cover debt.                    (33)
```

It applies when a marginal relaxation is nearly sharp but the simultaneous
arithmetic realization of its best rows is highly resonant. QED.
