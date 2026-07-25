---
id: THM-2226
title: "Three-checkpoint Bellman sieve and eight-profile geometric residue"
status: >
  PROVED + VERIFIED-EXACT. In THM-2222's scalar five-unit/three-blocker
  branch, lambda_1>=4 forces three consecutive 169-checkpoints. After
  primitive normalization, the three-clause Bellman transfer relaxation
  depends on only 16 capped 169-adic offset schedules. Exactly one schedule,
  (0,1,2), fails the 961/6930 target; all other 15 pass, with worst passing
  bound 7444136/62748517. Consequently 217 of the 225 profiles with
  lambda_1 in {4,5} are empty, leaving exactly eight explicitly listed
  profiles. Together with the completed depth-four ledger THM-2219 and the
  high-first-depth closure THM-2224, this reduces the current scalar
  valuation ledger from 675 to 458. The geometric chain genuinely defeats
  a uniform three-checkpoint inequality, and LRC(14) remains open.
source: codex-2026-07-24-three-checkpoint-clause-bellman
depends_on:
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2224-transfer-owner-word-temporal-union-bound
related:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2219-scalar-depth-four-sparse-tail-exclusion
script: 04-computation/lrc14_three_checkpoint_clause_bellman_thm2226.py
output: 05-knowledge/results/lrc14_three_checkpoint_clause_bellman_thm2226.out
script_sha256: 6afb6f32905753f45710cbe6cdbc096f3e2d891751b9c0894f34864028277784
output_sha256: 7b891802b76f40cb450daa25fbfaf12ef8c6823d2b8f1c548b02549d1863267a
hash_basis: working-tree bytes (LF)
---

# THM-2226 -- the three-checkpoint Bellman sieve

THM-2224 closes every scalar valuation profile with `lambda_1>=6` by four
checkpoints. At first depths four and five, three checkpoints already close
all but one `169`-adic offset schedule. The surviving schedule is not an
artifact: THM-2222's geometric chain exceeds the target there.

## 1. Three checkpoints from the parity tower

Retain THM-2222's notation

```text
A_+={R>0},
delta_5=measure(A_+) floor =961/6930,
U_k=union_(j=1)^3 D_(c_j/13^k).
```

Its even parity inclusions, assembled in equation (18) there, give

```text
A_+ subset U_k                         for every even k<=lambda_1. (1)
```

If `lambda_1>=4`, retain the three checkpoints `k=0,2,4` and put

```text
d_j=c_j/13^4.                                             (2)
```

Then

```text
U_4=union_j D_(d_j),
U_2=union_j D_(169d_j),
U_0=union_j D_(169^2d_j).                                (3)
```

Consequently

```text
A_+ subset K_3(d),

K_3(d)
 =intersection_(r=0)^2 union_(j=1)^3D_(169^r d_j),       (4)

measure K_3(d)>=961/6930.                                (5)
```

The common-gcd pullback is Haar preserving. Divide the triple `d` by its
gcd and reorder it without changing (4). For the resulting primitive triple,

```text
v_j=nu_13(d_j)=lambda_j-lambda_1,
min_j v_j=0.                                             (6)
```

Here and below `d` denotes the normalized triple.

## 2. Three-clause Bellman carrier

Put `A=169` and write

```text
d_j=A^(h_j)b_j,
h_j=floor(v_j/2),
nu_13(b_j) in {0,1}.                                    (7)
```

Exactly as in THM-2224, for almost every terminal point the number of
`A`-preimages in `D_(b_j)` is at most

```text
25,                  if 13 does not divide b_j,
26,                  if 13 divides b_j but 169 does not. (8)
```

Use the uniform cap `q=26/169`. Formula (4) is the three-clause event

```text
C_r=OR_j Z_(j,r+h_j),                  0<=r<=2,
Z_(j,t)(x)=1_(D_(b_j))(A^t x).                         (9)
```

At time `t`, clause `r` has conditional marginal cap

```text
q_(t,r)=26#{j:h_j=t-r}/169.                           (10)
```

Let `Omega={0,1,2}`. For a function `w:2^Omega -> R`, define the one-slice
Bellman relaxation

```text
(E_t w)(S)
 =max sum_Z p_Z w(S union Z),                         (11)
```

over all joint laws with

```text
p_Z>=0,
sum_Z p_Z=1,
sum_(Z:r in Z)p_Z<=q_(t,r)       for active r.        (12)
```

Starting from `w_(-1)(S)=1_(S=Omega)`, iterate

```text
w_t=E_t w_(t-1).                                     (13)
```

The backward-preimage-chain induction in THM-2224 is independent of the
number of clauses. Thus `w_t(S)` is a uniform a.e. conditional upper bound
for chains through times `0,...,t`, given that times greater than `t` have
already satisfied the clause set `S`. If `H=max_j h_j+2`, it gives

```text
measure K_3(d)<=w_H(empty).                           (14)
```

No independence or extra intersection constraint is assumed: the LP admits
every joint root-induced clause law satisfying the marginal caps.

## 3. Sixteen schedules and the exact boundary

Sort `0=h_1<=h_2<=h_3`. Each component contributes three consecutive slices

```text
h_j,h_j+1,h_j+2.                                     (15)
```

An offset gap larger than three inserts only empty identity slices. Therefore
the Bellman value depends only on

```text
(min(h_2-h_1,3),min(h_3-h_2,3))
 in {0,1,2,3}^2.                                     (16)
```

There are 16 schedules.

Run

```bash
python3 04-computation/lrc14_three_checkpoint_clause_bellman_thm2226.py
python3 -O 04-computation/lrc14_three_checkpoint_clause_bellman_thm2226.py
```

The exact referee independently enumerates every primal basic feasible
solution and every dual vertex for every clause state and schedule, requiring
exact equality. Ordinary and optimized outputs are byte-identical. Exactly
one schedule fails:

```text
offsets (0,1,2):
Bellman bound =5934/28561
               =961/6930+13675499/197927730.          (17)
```

Every other schedule passes. The worst passing row is

```text
offsets (0,2,4):
Bellman bound =7444136/62748517
               =961/6930-8713462357/434847222810.     (18)
```

The stored output records all 16 exact rows, three gap-compression controls,
and the full valuation census.

The failure in (17) reflects a genuine boundary. For the geometric triple

```text
d=(1,169,169^2),
```

THM-2222 computes

```text
measure K_3(d)=916159/4826809>961/6930.                (19)
```

Thus no uniform three-checkpoint estimate at the target can close the
`(0,1,2)` schedule. No coefficient triple is claimed to realize the larger
relaxed value in (17).

## 4. Exact valuation consequence

THM-2224 already forces `lambda_1<=5`. Apply the three-checkpoint sieve to
`lambda_1 in {4,5}`. By (6)--(7), the unique failed offset schedule is
equivalent to

```text
floor((lambda_2-lambda_1)/2)=1,
floor((lambda_3-lambda_1)/2)=2,                       (20)
```

or

```text
lambda_2-lambda_1 in {2,3},
lambda_3-lambda_1 in {4,5}.                           (21)
```

There are

```text
sum_(lambda_2=4)^18(19-lambda_2)=120
```

profiles with first depth four and

```text
sum_(lambda_2=5)^18(19-lambda_2)=105
```

with first depth five. The sieve closes `116/120` and `101/105`,
respectively. The only eight not closed here are

```text
(4,6,8), (4,6,9), (4,7,8), (4,7,9),
(5,7,9), (5,7,10), (5,8,9), (5,8,10).                (22)
```

Thus this theorem empties `217` additional profiles. These all have first
depth four or five and hence are disjoint from THM-2219's four profiles with
deepest depth four. THM-2219 and the two earlier depth-four theorems leave
`1130` profiles; THM-2224 removes `455`; therefore the present theorem leaves

```text
1130-455-217=458.                                    (23)
```

The eight rows in (22), the other low-first-depth profiles, owner/current
data, and LRC(14) remain open. QED.
