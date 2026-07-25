---
id: THM-2233
title: "Guard-danger hidden-state Bellman profile exclusion"
status: >
  PROVED + VERIFIED-EXACT. In the scalar five-unit/three-blocker branch,
  retain the guard bit and the three unit-core danger bits at every
  multiplication-by-thirteen time. Their exact backward root marginals are
  (10-G_next)/13 and (2-X_next)/13. An arbitrary-coupling hidden-state
  Bellman relaxation, with the guard at time zero and every available even
  divided-blocker union as labelled clauses, closes 452 shallow profiles.
  Relative to THM-2229's 240-profile ledger it newly excludes exactly the
  fifteen rows (2,2,c), 5<=c<=19, and the geometric row (5,7,9). The exact
  combined scalar ledger is therefore 224: 165 first-depth-one rows, 29
  first-depth-two rows, 29 first-depth-three rows, and (4,6,8). Every one of
  184,749 distinct rational Bellman LPs has matching exact primal and dual
  certificates. LRC(14) remains open.
source: codex-2026-07-25-guard-danger-hidden-state
depends_on:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2219-scalar-depth-four-sparse-tail-exclusion
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2224-transfer-owner-word-temporal-union-bound
  - THM-2229-unit-time-positive-set-bellman-profile-exclusion
related:
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
  - THM-2227-sharp-parity-three-checkpoint-bellman-profile-exclusion
  - THM-2232-same-core-signed-eigen-markov-dual-exclusion
script: 04-computation/lrc14_guard_danger_hidden_state_bellman_thm2233.py
output: 05-knowledge/results/lrc14_guard_danger_hidden_state_bellman_thm2233.out
script_sha256: ac67cb17681feb09a4a3c4ae6bc4190ec3302ae8f2850a00e7814580b8e2e556
output_sha256: 35475450148a913b32aff11aba2244c32392c7951e5b2504a41fd116b6440e1f
hash_basis: working-tree bytes (LF)
---

# THM-2233 -- keep the guard's parent bit

THM-2229 retains the positive set as a clause but replaces it on each
thirteen-root fibre by the uniform cap `10/13`. The cap forgets whether the
terminal point is itself guard-safe. Keeping that one parent bit, together
with the three corresponding danger parent bits, closes sixteen more rows.

## 1. The guard/checkpoint event

Use the scalar `5+3` notation of THM-2198 and THM-2222. Thus

```text
C_H={x:||Hx||>1/7},
D_a={x:||ax||<1/14},

c_j=13^(lambda_j)u_j,       13 does not divide u_j,
1<=lambda_1<=lambda_2<lambda_3<=19,                 (1)
```

where `H,q_1,...,q_5,u_1,u_2,u_3` are thirteen-units. The five unit
terminal coefficients define

```text
A_+=C_H setminus union_(i=1)^5 D_(q_i).             (2)
```

Here `setminus` in (2) denotes set difference. THM-2198 gives

```text
measure(A_+)>=delta_5=961/6930.                     (3)
```

For every even `k` with `0<=k<=lambda_1`, THM-2222 gives

```text
A_+ subset U_k,
U_k=union_(j=1)^3D_(c_j/13^k).                      (4)
```

Consequently

```text
A_+ subset E(lambda),

E(lambda)
 =C_H intersection
   intersection_(0<=k<=lambda_1, k even) U_k.       (5)
```

Replacing `A_+` by `C_H` in (5) is a safe enlargement. Unlike the cap used
in THM-2229, however, the construction below retains the exact guard parent
bit along every backward root step.

## 2. Exact four-bit root law

Put `Sx=13x mod 1` and define

```text
G_t(x)=1_(C_H)(S^t x),
X_(j,t)(x)=1_(D_(u_j))(S^t x),       1<=j<=3.       (6)
```

The divided-blocker clauses in (4) become literal unit-time clauses:

```text
U_k={x: OR_(j=1)^3 X_(j,lambda_j-k)(x)=1}.          (7)
```

For almost every terminal point `y`, exact root counting gives

```text
sum_(Sx=y)1_(C_H)(x)=10-1_(C_H)(y),
sum_(Sx=y)1_(D_(u_j))(x)=2-1_(D_(u_j))(y).          (8)
```

The first identity is THM-2222's guard-root law; the second is its
unit-danger root law. Therefore, conditional on the future bit vector

```text
xi=(xi_0,xi_1,xi_2,xi_3)
   =(G_(t+1),X_(1,t+1),X_(2,t+1),X_(3,t+1)),        (9)
```

the actual current vector

```text
eta=(G_t,X_(1,t),X_(2,t),X_(3,t))                   (10)
```

has exact one-coordinate marginals

```text
P(eta_0=1|xi)=(10-xi_0)/13,
P(eta_j=1|xi)=(2-xi_j)/13,          1<=j<=3.       (11)
```

The four coordinates share the same root digit and need not be independent.
We therefore allow every joint law on `{0,1}^4` with the marginals (11).
The actual root-induced law is feasible.

Under Haar measure the terminal vector has arbitrary dependence but exact
marginals

```text
P(G=1)=5/7,
P(X_j=1)=1/7,                       1<=j<=3.        (12)
```

The guard value `5/7` also follows by taking expectations in the first
identity of (8).

## 3. Clause-and-parent-bit Bellman recurrence

Let

```text
K={0,2,...,2 floor(lambda_1/2)},
Omega={star} union K.                                (13)
```

Clause `star` is hit by the guard bit at time zero. Clause `k in K` is hit
by blocker `j` at time `lambda_j-k`. For a current bit vector `eta`, let
`Z_t(eta)` be the set of clauses hit at time `t` by these incidences.

For `B subset Omega` and `xi in {0,1}^4`, start with

```text
V_(-1)(B,xi)=1_(B=Omega).                            (14)
```

Recursively, for `0<=t<=lambda_3`, set

```text
V_t(B,xi)
 =max_p sum_(eta in {0,1}^4)
       p_eta V_(t-1)(B union Z_t(eta),eta),          (15)
```

where

```text
p_eta>=0,       sum_eta p_eta=1,

sum_(eta:eta_0=1)p_eta=(10-xi_0)/13,
sum_(eta:eta_j=1)p_eta=(2-xi_j)/13,  1<=j<=3.       (16)
```

Finally maximize

```text
B(lambda)
 =max_pi sum_xi pi_xi V_(lambda_3)(emptyset,xi),     (17)
```

over all terminal laws `pi` satisfying (12).

The backward-root induction proves

```text
measure(E(lambda))<=B(lambda).                       (18)
```

Indeed, for a fixed generic terminal point, the empirical distribution on
its thirteen roots satisfies (16) exactly by (8). Conditional on a chosen
root, the induction hypothesis bounds the earlier chain by the corresponding
term in (15). Averaging and enlarging to every coupling allowed by (16)
proves the pointwise step. At the final time, the actual Haar joint law is
feasible in (17). The finitely many interval endpoints and all of their
finitely many images and preimages are null.

This proof uses neither independence nor a hidden equality between different
coefficients. Coincident root sheets, every cross-core correlation, and
arbitrary terminal dependence are admitted by the relaxation.

Combining (3), (5), and (18), a scalar profile is impossible whenever

```text
B(lambda)<961/6930.                                 (19)
```

## 4. Exact rational LP certificate

For each future vector `xi`, the feasible laws in (16) form a rational
polytope with sixteen variables and five equality constraints. Its column
matrix is

```text
A_eta=(1,eta_0,eta_1,eta_2,eta_3)^T.                (20)
```

There are exactly `3008` invertible five-column bases among the
`binom(16,5)` candidates. The companion enumerates them once and, for every
one of the sixteen conditional right sides and the terminal right side,
retains every distinct nonnegative basic distribution. Their exact vertex
counts are

```text
(94,73,73,41,73,41,41,30,
 34,71,71,36,71,36,36,30; 34 terminal).             (21)
```

Every Bellman maximum is first obtained as an exact dot-product maximum over
these vertices. It is then independently certified from an optimal basis
`I`: the script solves

```text
A_I^T y=c_I                                         (22)
```

and checks, with rational arithmetic,

```text
A_eta^T y>=c_eta            for all sixteen eta,
b^T y=the primal value.                             (23)
```

Thus (23) is a feasible dual certificate with the same value. Across the
complete census there are `184749` distinct objective/right-side pairs, and
every one receives such a certificate. Repeated Bellman calls reuse only an
already certified exact pair.

## 5. The geometric boundary and the full shallow census

For the two profiles left at first depth at least four by THM-2229, the exact
bounds are

```text
B(4,6,8)
 =800688435/5710115047
 =961/6930+8764327769/5653013896530,                 (24)

B(5,7,9)
 =10149510390/74231495611
 =961/6930-142908611353/73489180654890.              (25)
```

Therefore (19) excludes `(5,7,9)` but not `(4,6,8)`.

The companion evaluates all `685` legal profiles with `lambda_1<=5`. The
strictly passing counts are

```text
lambda_1=1:   0/171,
lambda_1=2: 121/153,
lambda_1=3: 107/136,
lambda_1=4: 119/120,
lambda_1=5: 105/105.                                (26)
```

The exact closure classification is as follows.

At first depth two, the passing rows are exactly

```text
lambda_2=2 and lambda_3>=5,

or lambda_2=3 and lambda_3>=6,

or lambda_2>=5 and lambda_3!=lambda_2+2.             (27)
```

At first depth three they are exactly

```text
lambda_2=3 and (lambda_3=4 or lambda_3>=6),

or lambda_2=4 and lambda_3>=7,

or lambda_2>=6 and lambda_3!=lambda_2+2.             (28)
```

At first depth four every row except `(4,6,8)` passes, and every first-depth
five row passes. No first-depth-one row passes this test. These set
classifications, not merely their counts, are checked against frozen
SHA-256 digests.

The raw passing set has `452` rows. Its intersection with THM-2229's raw
passing set is all `436` rows of that theorem, and it contains all six
THM-2227 rows. Its overlap with the ten shallow exact exclusions is the
single row `(3,3,4)`. Relative to the exact `240`-profile combined ledger
before this theorem, the novel set is exactly

```text
{(2,2,c):5<=c<=19} union {(5,7,9)},                  (29)
```

of size sixteen. Hence the new combined scalar ledger has `224` rows:

```text
lambda_1=1:
  all 165 rows with lambda_3>=5;

lambda_1=2:
  (2,3,5);
  (2,4,c), 5<=c<=19;
  (2,b,b+2), 5<=b<=17;                              (30)

lambda_1=3:
  (3,3,5), (3,4,5), (3,4,6);
  (3,5,c), 6<=c<=19;
  (3,b,b+2), 6<=b<=17;                              (31)

lambda_1>=4:
  (4,6,8).                                          (32)
```

The counts in (30) and (31) are twenty-nine each. Together with the 165
first-depth-one rows and (32), they total

```text
165+29+29+1=224.                                    (33)
```

## 6. Connection and loss ledger

The successful carrier is

```text
source:
  the positive residual A_+, the exact guard-root law, and the even
  divided-blocker tower;

map:
  enlarge A_+ to its guard C_H, express every divided blocker at its exact
  unit time, and retain the guard plus three blocker parent bits;

preserved:
  every even checkpoint label, exact integer times, the exact 9-versus-10
  guard count, exact 1-versus-2 danger counts, and terminal marginals;

destroyed:
  avoidance of the five unit terminal dangers inside A_+, the actual common
  root digit, coefficient phases, root-sheet incidence, owners, and all
  joint correlations beyond one-coordinate conditional marginals;

cheapest hostile control:
  (4,6,8), whose exact relaxed excess is displayed in (24);

needed continuation:
  common-root phase/incidence, the chord/deficiency carrier of THM-2192 and
  THM-2197, or another sidecar that constrains the arbitrary coupling.     (34)
```

The improvement over THM-2229 is not a smaller static guard cap. It is the
dependence of the current guard capacity on the future guard bit. Conversely,
(24) proves a precise stopping obstruction: guard plus blocker parent bits
and their exact one-coordinate Markov laws do not by themselves close the
last all-even geometric row.

## 7. Reproduction

Run

```bash
python3 04-computation/lrc14_guard_danger_hidden_state_bellman_thm2233.py
python3 -O 04-computation/lrc14_guard_danger_hidden_state_bellman_thm2233.py
```

Both modes produce the stored transcript byte for byte. All load-bearing
checks use explicit exceptions rather than Python assertions. The script
freezes the complete universe, exact overlap and remaining sets, three set
digests, both geometric fractions and margins, every polytope census, and
exact primal/dual equality for every distinct LP. This theorem reduces only
the scalar valuation ledger. Owner/current branches and LRC(14) remain open.
QED.
