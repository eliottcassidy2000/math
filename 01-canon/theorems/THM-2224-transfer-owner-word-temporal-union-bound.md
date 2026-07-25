---
id: THM-2224
title: "Clause-state Bellman transfer bound for four danger-comb checkpoints"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED. For every positive primitive
  coefficient triple, the four consecutive 169-checkpoint intersection of
  its three danger-comb union has measure at most 3272/28561, strictly below
  961/6930. The proof strips deterministic 169-adic shifts, retains the four
  labelled clauses, and applies a uniform backward-chain Bellman relaxation
  using the 26/169 root-sheet cap. An exact 25-schedule rational LP referee
  independently enumerates primal and dual vertices; ordinary and optimized
  outputs are byte-identical. Consequently THM-2222's scalar branch has
  lambda_1<=5, removing all 455 profiles with lambda_1>=6. Together with
  the complete depth-four exclusion THM-2219, this leaves 675 profiles
  before the later THM-2226 sieve. The relaxed equality value is not claimed
  realizable, and LRC(14) remains open.
source: codex-2026-07-24-four-checkpoint-clause-bellman
depends_on:
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
related:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2210-nested-binomial-minorant-and-adaptive-moment-lp-hierarchy
  - THM-2219-scalar-depth-four-sparse-tail-exclusion
  - THM-2226-three-checkpoint-bellman-sieve-and-eight-profile-residue
script: 04-computation/lrc14_four_checkpoint_clause_bellman_thm2224.py
output: 05-knowledge/results/lrc14_four_checkpoint_clause_bellman_thm2224.out
script_sha256: 579010b917bc87f447886480a5fdb314ae4a086f8fa901c67c37edadceb3ea55
output_sha256: dc47521eccc3ba313a85f67c2fb5fef8984406e13731ac48865a6e8af8c761eb
hash_basis: working-tree bytes (LF)
---

# THM-2224 -- the clause-state Bellman transfer bound

The open extremal inequality isolated by THM-2222 holds uniformly, with no
height restriction:

```text
S_4(B) <= 3272/28561
       < 961/6930                         for every B.             (1)
```

The exact gap is

```text
961/6930 - 3272/28561 = 4772161/197927730 > 0.                    (2)
```

The mechanism is a labelled-clause transfer supermartingale. Powers of `169`
remain as time offsets, while the remaining arithmetic is relaxed to a
uniform conditional root-sheet cap.

## 1. Strip the deterministic `169` shifts

Put `A=169`, and let

```text
D_a={x in R/Z : ||a x||<1/14}.
```

The common-gcd pullback is Haar preserving, so normalize a positive triple
`d=(d_1,d_2,d_3)` by its gcd. Then `min_j nu_13(d_j)=0`. Write

```text
d_j=A^(h_j)b_j,
h_j=floor(nu_13(d_j)/2),
nu_13(b_j) in {0,1}.                                  (3)
```

For checkpoints `r=0,1,2,3`,

```text
1_(D_(A^r d_j))(x)=1_(D_(b_j))(A^(r+h_j)x).           (4)
```

Thus the four-checkpoint event is the four-clause Boolean formula

```text
C_r = OR_(j=1)^3 Z_(j,r+h_j),                         (5)
Z_(j,t)(x)=1_(D_(b_j))(A^t x),
K_4(d)={C_0=C_1=C_2=C_3=1}.                           (6)
```

This representation retains the clause labels. Expanding the clauses into
`3^4` owner words loses their overlap and gives a weaker bound.

## 2. Uniform root-sheet cap

For the circle endomorphism `T(x)=Ax`, condition on a terminal point `y`.
Its `A` preimages are `(y+k)/A`, `0<=k<A`.

If `13` does not divide `b`, the two successive root-count identities give,
for almost every `y`,

```text
sum_(Tx=y) 1_(D_b)(x)=24+1_(D_b)(y)<=25.              (7)
```

If `b=13u` with `13` not dividing `u`, grouping the `A` roots by their
residue modulo `13` gives, again for almost every `y`,

```text
sum_(Tx=y) 1_(D_b)(x)
 =13 sum_(13z=y)1_(D_u)(z)
 =13(2-1_(D_u)(y))
 <=26.                                                (8)
```

Only finitely many terminal endpoints are exceptional at each finite stage;
their images and preimages remain null. Therefore every stripped current
danger atom occupies at most

```text
q=26/169                                               (9)
```

of almost every root sheet, uniformly in the terminal point. We deliberately
use `26` also in the unit case; the lost unit saving is not needed.

At time `t`, clause `r` has

```text
m_(t,r)=#{j:h_j=t-r}                                  (10)
```

current atoms. Their union therefore has conditional marginal at most

```text
q_(t,r)=26m_(t,r)/169.                                (11)
```

Since there are only three components, this number is always below one.
No independence between clauses is asserted.

## 3. The exact clause-state Bellman relaxation

Let `Omega={0,1,2,3}` be the clause set. For a function
`v:2^Omega -> R`, and the active clauses at time `t`, define

```text
(E_t v)(S)
 = max sum_(Z subset Omega) p_Z v(S union Z),         (12)
```

where the maximum is over

```text
p_Z>=0,
sum_Z p_Z=1,
sum_(Z:r in Z)p_Z <= q_(t,r)    for every active r,  (13)
```

and inactive clauses never occur in `Z`. This LP permits every coupling
consistent with the one-clause marginal caps. In particular, the actual
root-induced distribution is feasible for almost every terminal point;
allowing nonintegral multiples of `1/169` only enlarges the feasible set.

Set

```text
v_(-1)(S)=1_(S=Omega),
v_t=E_t v_(t-1).                                     (14)
```

Here `v_t(S)` is a uniform conditional upper bound for a backward preimage
chain through times `0,...,t`, given that variables at times strictly greater
than `t` have already satisfied exactly the clauses recorded by `S`. This is
a conditional chain induction, not a pointwise comparison of indicators.

Fix a nonexceptional terminal `x_(t+1)`. Average uniformly over all backward
chains

```text
x_0 -> x_1 -> ... -> x_t -> x_(t+1),
T x_s=x_(s+1).                                       (15)
```

For a chosen root `x_t`, let `Z_t(x_t)` be the set of clauses satisfied by
current atoms at time `t`. By induction, the conditional success probability
over the earlier roots is at most

```text
v_(t-1)(S union Z_t(x_t)).                            (16)
```

Across the `169` choices of `x_t`, the law of `Z_t` obeys every marginal
constraint (13), by (11). Averaging (16) and relaxing to all distributions
in (13) proves the uniform a.e. bound `v_t(S)`.

Let

```text
H=max_j h_j+3.                                       (17)
```

After time `H` there are no variables, so the future-satisfied set is empty.
Averaging the uniform-in-terminal bound over the terminal point gives

```text
measure K_4(d) <= v_H(empty).                         (18)
```

This proof needs only marginals because the LP includes arbitrary joint laws.
Shared root atoms, coincident cores, and correlations can only select one of
the already admitted joint distributions.

## 4. Only 25 offset schedules

The cap (9) is the same for all three components, so sort

```text
0=h_1<=h_2<=h_3.                                     (19)
```

The first equality follows from primitive normalization. Component `j`
contributes the four consecutive time slices

```text
h_j,h_j+1,h_j+2,h_j+3,                               (20)
```

with clause labels `0,1,2,3`.

Write the two start-time gaps as

```text
g_1=h_2-h_1,
g_2=h_3-h_2.                                         (21)
```

If either gap exceeds `4`, the extra slices between the two corresponding
four-slice blocks are empty. For an empty slice `E_t` is the identity.
Deleting those empty slices changes neither the ordered list of nonempty
clause-cap packets nor (18). Hence the Bellman bound depends only on

```text
(min(g_1,4),min(g_2,4)) in {0,1,2,3,4}^2.            (22)
```

There are exactly 25 schedules.

## 5. Exact LP audit

For a time slice with `m<=3` active clauses, (12)--(13) has at most
`2^m<=8` probability variables. Add one slack variable to every marginal
constraint. The resulting standard-form system has rank `m+1<=4` and at
most `2^m+m<=11` columns. Every optimum occurs at a basic feasible solution,
so enumerating at most

```text
binom(11,4)=330                                       (23)
```

bases solves each LP exactly.

Run

```bash
python3 04-computation/lrc14_four_checkpoint_clause_bellman_thm2224.py
python3 -O 04-computation/lrc14_four_checkpoint_clause_bellman_thm2224.py
```

The companion uses rational arithmetic, enumerates every primal basis,
independently enumerates every dual vertex, and requires exact primal-dual
equality for every state at every time in all 25 schedules. Ordinary and
optimized outputs are byte-identical. It also checks three hostile
gap-compression controls and the geometric-chain control from THM-2222.

The maximum relaxed value is unique:

```text
(h_1,h_2,h_3)=(0,1,2),
v_H(empty)=3272/28561.                               (24)
```

The next bounds are

```text
(0,2,4): 73891216/815730721,
(0,2,3),(0,1,3): 377276/4826809,
(0,3,4),(0,1,4): 54096584/815730721.                (25)
```

All are smaller than (24). The stored output prints all 25 exact rows.
Combining (18), (22), (24), and (2) proves (1).

No coefficient triple is claimed to realize equality in the relaxed Bellman
bound. At least one stripped core is in fact a `13`-unit and has the sharper
cap `25/169`, but the strict relaxed estimate already suffices.

## 6. Scalar consequence and scope

THM-2222 proves that any scalar survivor with `lambda_1>=6` would supply a
primitive triple with

```text
measure K_4(d)>=961/6930.
```

This contradicts (1). Therefore

```text
lambda_1<=5                                           (26)
```

in every scalar survivor. The `455` profiles with `lambda_1>=6` are empty.
After the complete depth-four exclusions THM-2213/2215/2219, the current
scalar valuation ledger falls from `1,130` to `675`, all satisfying

```text
lambda_1<=5<=lambda_3.
```

The theorem is uniform in `B`; it does not use THM-2222's enormous finite
height once the four-checkpoint reduction has been reached. It does not by
itself settle the remaining `675` scalar profiles after THM-2219, and the
later THM-2226 sieve still does not settle LRC(14). QED.
