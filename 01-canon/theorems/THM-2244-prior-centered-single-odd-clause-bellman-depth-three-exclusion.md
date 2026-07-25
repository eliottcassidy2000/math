---
id: THM-2244
title: "Prior-centered single-odd-clause Bellman exclusion of depth-three scalar profiles"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED. In the scalar five-unit plus
  three-blocker branch, retain every even checkpoint and the single odd
  checkpoint k=1. For each normalized blocker core, center that odd atom at
  its earliest strictly prior time of the opposite atom parity. The resulting
  affine score is response-null and pointwise nonpositive on the negative
  carrier. Its exact three-bit arbitrary-coupling Bellman recursion
  independently closes 57 low-first profiles in the 223-row reference
  census: all 29 first-depth-two profiles and 28 of the 29 first-depth-three
  profiles. These rows are now all also closed by the stronger nonnegative
  charge in THM-2239, so this theorem supplies an alternate mechanism and
  hostile audit but no additional current ledger decrement. Both routes
  leave the same 166-row ledger: the 165 first-depth-one rows and (3,4,5).
  The closest new pass is (3,5,6), with exact bound
  1511656180038/12545122758259<961/6930. LRC(14) remains open.
source: codex-2026-07-25-prior-centered-single-odd-clause
depends_on:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2233-guard-danger-hidden-state-bellman-profile-exclusion
  - THM-2239-unrestricted-multicore-signed-dual-profile-exclusion
related:
  - THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy
  - THM-2232-same-core-signed-eigen-markov-dual-exclusion
  - THM-2234-first-depth-one-private-two-owner-mass-and-two-step-expansion
script: 04-computation/lrc14_prior_centered_single_odd_clause_bellman_thm2244.py
output: 05-knowledge/results/lrc14_prior_centered_single_odd_clause_bellman_thm2244.out
script_sha256: 55c009410f347820abff74d6741edc310bb45a250b01b0ec5034efcfef3cac14
output_sha256: 0cc0d068e903c3b0c9623e10475bfaeb9208b4aec2a80d9a6c765cd50c528e26
hash_basis: working-tree bytes (LF)
---

# THM-2244 -- the first odd clause closes all but one depth-three row

An earlier, weaker form of THM-2239 left the intermediate scalar valuation
ledger

```text
194 = 165 first-depth-one profiles + 29 first-depth-three profiles. (1)
```

A single odd checkpoint, centered separately and strictly earlier on each
normalized blocker core, excludes 28 of the 29 depth-three rows. The sole
survivor is `(3,4,5)`. The strengthened current THM-2239 subsequently
subsumes all 57 exclusions proved here.

The point of the construction is not the number of retained clauses. It is
the freedom to choose a different prior response base for each labelled core.
That coordinate was absent from the earlier common-time score and remains a
useful independent control on THM-2239's stronger shifted charge.

## 1. Scalar residual and the intermediate ledger

On `T=R/Z`, put

```text
D_b={x:||bx||<1/14},       C_H={x:||Hx||>1/7}.       (2)
```

Use the scalar `5+3` branch of THM-2198 and write

```text
R=1_(C_H)-sum_(i=1)^5 1_(D_(q_i)),

c_j=13^(lambda_j)u_j,       13 does not divide u_j,
1<=lambda_1<=lambda_2<lambda_3<=19.                 (3)
```

The guard, the five `q_i`, and the three `u_j` are thirteen-units. No
relation among the `u_j` is assumed. As in THM-2239,

```text
p:=integral R_+=integral R_- >=delta_5:=961/6930,
0<=R_+<=1.                                          (4)
```

Let `Sx=13x mod 1` and define

```text
X_(j,t)(x)=1_(D_(u_j))(S^t x),
B_k=sum_(j=1)^3X_(j,lambda_j-k).                    (5)
```

For

```text
E={0,2,...,2 floor(lambda_1/2)},

P=intersection_(k in E){B_k>=1},
N_1={B_1>=1},                                       (6)
```

THM-2222 gives, almost everywhere,

```text
support(R_+) subset P,
support(R_-) subset N_1.                            (7)
```

The guard condition in the positive support has been discarded in `P`.
This is a safe enlargement, and it is why the Bellman state below needs only
the three labelled danger bits rather than THM-2239's guard bit as well.

THM-2239's 29 surviving depth-three profiles are

```text
(3,3,5), (3,4,5), (3,4,6);
(3,5,c), 6<=c<=19;
(3,b,b+2), 6<=b<=17.                                (8)
```

## 2. Separate prior centering

Assume first that `lambda_1` is two or three. For each core put

```text
t_j=lambda_j-1,
s_j=lambda_j mod 2 in {0,1},
d_j=t_j-s_j.                                        (9)
```

Because `lambda_j>=lambda_1>=2`, every `d_j` is a positive odd integer.
In particular, `s_j<t_j`. Set

```text
rho=-1/13,
c_j=rho^(d_j)<0,
Y_j=X_(j,t_j)-c_jX_(j,s_j),

q=1-sum_(j=1)^3Y_j.                                 (10)
```

The unnormalized transfer identity `L^tR=(-1)^tR` from THM-2222 and transfer
duality imply, for every labelled core and every `t>=s`,

```text
integral R X_(j,t)
 =rho^(t-s) integral R X_(j,s).                     (11)
```

Applying (11) with `(t,s)=(t_j,s_j)` gives

```text
integral RY_j=0,
integral Rq=0.                                      (12)
```

This cancellation is core-by-core. Coincident blockers, arbitrary arithmetic
relations among the `u_j`, and arbitrary cross-core root dependence do not
change it.

## 3. The negative carrier has the right sign

Since every `c_j` is negative,

```text
Y_j=X_(j,t_j)+|c_j|X_(j,s_j)>=X_(j,t_j).            (13)
```

On `N_1`, at least one of the three raw atoms `X_(j,t_j)` is one. Hence

```text
q<=1-sum_jX_(j,t_j)<=0                 on N_1.       (14)
```

Equations (4), (7), (12), and (14) now give

```text
p
 =integral [R_+(1-q)+R_-q]
 <=integral 1_P(1-q)_+.                             (15)
```

No independence assumption, guard cap, or same-core hypothesis enters (15).
The positive part cannot in general be removed before expanding the score,
although for this one-clause choice

```text
1-q=sum_jX_(j,t_j)-sum_jc_jX_(j,s_j)>=0.            (16)
```

Thus the remaining task is an exact upper bound for the expectation in (15).

## 4. Three-bit arbitrary-coupling Bellman recursion

Write

```text
Z_t=(X_(1,t),X_(2,t),X_(3,t)) in {0,1}^3.           (17)
```

For a fixed future state `xi`, exact thirteen-root counting gives the current
one-coordinate marginals

```text
P(Z_current,j=1 | Z_future=xi)
 =(2-xi_j)/13.                                      (18)
```

Let `C(xi)` be the polytope of all laws on `{0,1}^3` with the three marginals
in (18). The actual common-root law is one feasible point of `C(xi)`.
Maximizing over the full polytope therefore discards all cross-core joint
information but remains an upper bound.

The automaton summary is

```text
(M,r,b),                                             (19)
```

where

```text
M = the set of even clauses in E already hit,
r = the number of the three raw atoms X_(j,t_j) hit,
b_j = the retained value of X_(j,s_j).              (20)
```

When every scheduled time has been processed, its payoff is exactly

```text
1_(M=E) (r-sum_jc_jb_j).                            (21)
```

For each time `t`, state `z`, and summary `sigma`, let
`Phi_t(sigma,z)` update all incidences in (20) scheduled at `t`. For a vector
of continuation values `v_z`, one Bellman step is

```text
max_(pi in C(xi)) sum_(z in {0,1}^3) pi_z v_z.       (22)
```

The recursion eliminates times `0,1,...,lambda_3-1` one root layer at a time.
At time `lambda_3` it updates (19) with the terminal state and maximizes over
every terminal joint law whose three marginals are `1/7`.

### Why the recursion is nonanticipative

Fix an exact parent point before choosing one of its thirteen roots. Its three
conditional marginals are already fixed by (18). The future automaton summary
is a function of the parent point and its forward orbit; it does not reveal
which root will be chosen. The actual distribution of the current three bits
is therefore feasible in (22). Conditional on a chosen root, the previous
Bellman value bounds the remaining earlier path. Iterating this pointwise
argument and then integrating proves that the terminal Bellman value bounds
the right side of (15).

The optimizer may select a different arbitrary coupling for every future
state and every automaton summary. This enlarges the actual arithmetic
process; it does not impose a Markov or independence model on the three
cores.

## 5. Exact primal-dual evaluation

For fixed marginals `m`, the primal polytope in (22) has eight nonnegative
variables and four equality constraints. Its columns are

```text
(1,z_1,z_2,z_3)^T,       z in {0,1}^3.              (23)
```

The companion enumerates all `58` invertible four-column bases. For every
distinct objective and marginal vector it first maximizes over all exact
nonnegative basic laws. It then solves the dual equations on a maximizing
primal basis and checks

```text
alpha+sum_j beta_j z_j >=v_z       for all z,

alpha+sum_jm_jbeta_j = primal maximum.              (24)
```

Thus every local maximum has an independently checked exact dual certificate.
Across the full package there are

```text
39,017 distinct exact primal-dual LP checks.         (25)
```

The package evaluates the 58 first-depth-two/three profiles in the 223-row
low-first reference census. Exactly 57 pass: all 29 first-depth-two profiles
and the following 28 first-depth-three profiles:

```text
(3,3,5), (3,4,6);
(3,5,c), 6<=c<=19;
(3,b,b+2), 6<=b<=17.                                (26)
```

The exact digest of every profile, base triple, correction triple, Bellman
bound, and target margin is

```text
de457ace002dc402475993cf0050c03a028f93eaac38171b2cc79d5614f5e511. (27)
```

The closest newly excluded row is `(3,5,6)`:

```text
B(3,5,6)
 =1511656180038/12545122758259,

delta_5-B(3,5,6)
 =32246645775991/1774238790096630
 >0.                                                (28)
```

Combining (4), (15), and (28) with the exact profile census proves that every
row in (26) is empty. The 29 depth-two passes are an independent overlap
control for THM-2239 and are not counted again. Relative to THM-2239's current
ledger at the time this route was developed, the then-new decrement was

```text
194-28=166.                                         (29)
```

The exact scalar residue is therefore

```text
all 165 profiles with lambda_1=1 and lambda_3>=5;
(3,4,5).                                            (30)
```

## 6. Exact failure boundary

### The exceptional depth-three row

For `(lambda_1,lambda_2,lambda_3)=(3,4,5)`, checkpoint one has atom times
`(2,3,4)`. All strictly prior bases giving odd relative exponents are

```text
(s_1,s_2,s_3) in
{(1,0,1),(1,2,1),(1,0,3),(1,2,3)}.                 (31)
```

Their exact Bellman bounds, in increasing order, are

```text
(1,0,1): 1393650030/5710115047,
(1,2,1): 108834462/439239619,
(1,0,3): 1456253766/5710115047,
(1,2,3): 672486/2599051.                            (32)
```

Every value in (32) exceeds `delta_5`. The best one misses by

```text
1393650030/5710115047-delta_5
 =85113758117/807573413790>0.                        (33)
```

Thus merely moving the three prior bases cannot close `(3,4,5)` inside this
arbitrary-coupling relaxation.

### Why first depth one is absent

At first depth one the only odd checkpoint is `k=1`, and the first labelled
atom occurs at

```text
t_1=lambda_1-k=0.                                  (34)
```

There is no nonnegative time `s_1<t_1`, so no strictly prior response center
exists. This is a structural boundary, not a missing loop iteration.

If one instead uses the least later center `s_1=1`, its coefficient is
`rho^(-1)=-13`. The positive carrier for even checkpoint zero contains
`X_(1,1)`, and the resulting pointwise majorant is at least
`13X_(1,1)` there. Its integral contribution alone is

```text
13 measure(X_(1,1)=1)=13/7>delta_5.                 (35)
```

Thus the naïve later-center repair is far outside the scalar target.
Depth one needs a different signed observable or the owner/private-mass
information of THM-2234.

## 7. Connection and loss ledger

```text
source:
  the signed scalar residual and the first odd transfer checkpoint;

map:
  center each labelled odd atom at its own earliest strictly prior
  same-parity base and retain the three base bits in the Bellman automaton;

preserved:
  residual sign, every even clause, the selected odd clause, labelled core
  identity, exact integer times, exact root marginals, and response nullity;

destroyed:
  the guard bit, all arithmetic phase, the actual common root digit,
  cross-core joint laws, owner/carry data, and odd checkpoints beyond k=1;

hostile boundary:
  all four prior centers fail at (3,4,5), while depth one has no prior base;

needed continuation:
  a second response-null direction, realized-root incidence, guard/owner
  coupling, or the private two-owner carrier of THM-2234.                 (36)
```

## 8. Reproduction

Run

```bash
python3 04-computation/lrc14_prior_centered_single_odd_clause_bellman_thm2244.py --exact-ledger
python3 -O 04-computation/lrc14_prior_centered_single_odd_clause_bellman_thm2244.py --exact-ledger
```

Both modes reproduce the stored transcript byte for byte. The companion
freezes:

```text
the 223-row low-first reference census and the intermediate 194-row ledger;
all 57 raw passes and the 28-row then-new intersection;
the complete profile/bound/margin digest;
the closest raw and closest novel passes;
all four prior-center values for (3,4,5);
the weakened guard-free three-bit (4,6,8) geometric control;
39,017 exact primal-dual LP certificates;
the exact 166-row consequence.                       (37)
```

This theorem supplies an independent scalar exclusion certificate. It adds
no decrement beyond current THM-2239; owner/current branches and LRC(14)
remain open. QED.
