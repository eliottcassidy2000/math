---
id: THM-2239
title: "Unrestricted multicore signed-dual profile exclusion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED + HOSTILE-AUDITED. In the
  scalar five-unit/three-blocker branch, separately center the three blocker
  danger processes along the signed transfer eigenline and add the least
  pointwise shift making every centered atom nonnegative and at least one
  when active. Selected odd-checkpoint clauses then annihilate the entire
  negative-support term in the signed dual. An exact arbitrary-coupling
  three-bit Bellman bound closes all 29 remaining first-depth-two profiles,
  28 of the 29 remaining first-depth-three profiles, and the unrestricted
  profile (4,6,8), with no equality or overlap hypothesis on the normalized
  cores. Relative to THM-2233, exactly 58 profiles are newly excluded and
  the exact scalar ledger falls from 224 to 166: the 165 first-depth-one
  profiles with deepest depth at least five, and (3,4,5). Every one of
  39,295 distinct rational coupling LPs has matching exact primal and dual
  certificates. LRC(14) remains open.
source: klein-codex-2026-07-25-unrestricted-multicore-signed-dual
depends_on:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2233-guard-danger-hidden-state-bellman-profile-exclusion
related:
  - THM-2080-unequal-comb-overlap-removes-depth-five
  - THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy
  - THM-2227-sharp-parity-three-checkpoint-bellman-profile-exclusion
  - THM-2229-unit-time-positive-set-bellman-profile-exclusion
  - THM-2232-same-core-signed-eigen-markov-dual-exclusion
script: 04-computation/lrc14_unrestricted_multicore_signed_dual_thm2239.py
output: 05-knowledge/results/lrc14_unrestricted_multicore_signed_dual_thm2239.out
script_sha256: 37ee166558bf6b90c4ce3be2e79c485eb205710a7e5b5fbe66285d6864f7d15f
output_sha256: b7d7b3d353945215be716b0b29671195585069e3415bffd5a8374635f8cffdf6
hash_basis: working-tree bytes (LF)
---

# THM-2239 -- a nonnegative multicore signed dual

The same-core calculation in THM-2232 shows that unsigned checkpoint mass is
not the right object: the positive and negative parts of the *same* residual
must be transported with opposite parity. The obstruction there appeared to
require a shared normalized blocker core. It does not. The missing move is to
center each core separately and then make the centered atoms pointwise
nonnegative.

## 1. Scalar residual and checkpoint carriers

On `T=R/Z`, put

```text
D_a={x:||ax||<1/14},       C_H={x:||Hx||>1/7}.       (1)
```

Work in THM-2198's scalar five-unit/three-blocker branch. Thus
`H,q_1,...,q_5,u_1,u_2,u_3` are thirteen-units,

```text
c_j=13^(lambda_j)u_j,

1<=lambda_1<=lambda_2<lambda_3<=19,                 (2)
```

and the signed five-unit residual is

```text
R=1_(C_H)-sum_(i=1)^5 1_(D_(q_i)).                  (3)
```

Write `R_+=max(R,0)`, `R_-=max(-R,0)`, and

```text
p=integral R_+=integral R_- >= delta_5:=961/6930.   (4)
```

The equality follows from `integral R=5/7-5/7=0`; the lower bound is
THM-2198. Pointwise,

```text
0<=R_+<=1,       0<=R_-<=5.                         (5)
```

Let `Sx=13x mod 1` and define three distinct danger processes, without any
relation among their cores:

```text
X_(j,t)(x)=1_(D_(u_j))(S^t x),       j=1,2,3.        (6)
```

For every checkpoint `0<=k<=lambda_1`, put

```text
U_k={OR_(j=1)^3 X_(j,lambda_j-k)=1}.                (7)
```

THM-2222's transfer-parity inclusions give

```text
{R_+>0} subset U_k       when k is even,
{R_->0} subset U_k       when k is odd.             (8)
```

Consequently, if

```text
P_lambda
 =intersection_(0<=k<=lambda_1, k even) U_k,        (9)
```

and `O` is any nonempty set of available odd checkpoints, then

```text
{R_+>0} subset P_lambda,
{R_->0} subset N_O:=intersection_(k in O)U_k.       (10)
```

Only the selected odd clauses in `O` will be needed below.

## 2. The nonnegative centered charge

Let `L` be the unnormalized transfer operator for `S`. THM-2222 proves

```text
L^tR=(-1)^tR.                                       (11)
```

Transfer duality applied separately to every unit core gives, with
`rho=-1/13`,

```text
integral R X_(j,t)
 =rho^t integral R X_(j,0).                         (12)
```

The raw centered atom

```text
Y_(j,t)=X_(j,t)-rho^tX_(j,0)                        (13)
```

has zero signed correlation with `R`, but for even `t` it can be negative.
Add exactly the parity-dependent constant needed to repair this:

```text
Z_(j,t)
 =Y_(j,t)+max(rho^t,0)
 =X_(j,t)-rho^tX_(j,0)+max(rho^t,0).                (14)
```

Since `integral R=0`, equations (12)--(14) imply

```text
integral R Z_(j,t)=0.                               (15)
```

More importantly, an exhaustive two-bit check gives the pointwise facts

```text
Z_(j,t)>=0,
X_(j,t)=1  implies  Z_(j,t)>=1.                     (16)
```

Indeed, for even `t`, `Z=X_t+rho^t(1-X_0)`; for odd `t`,
`Z=X_t+|rho|^tX_0`.

For `m=|O|`, define

```text
q_O
 =m-sum_(k in O)sum_(j=1)^3 Z_(j,lambda_j-k).       (17)
```

Equation (15) gives `integral Rq_O=0`. On `N_O`, each selected clause
contains at least one active danger atom. By (16), each clause contributes
at least one to the double sum in (17), with multiplicity retained.
Therefore

```text
q_O<=0 on N_O.                                      (18)
```

Now use `integral Rq_O=0` and (4):

```text
p
 =integral [R_+(1-q_O)+R_-q_O]
 <=integral 1_(P_lambda)(1-q_O)_+.                  (19)
```

The negative-support term is nonpositive by (10) and (18); the positive
term is bounded using (5). This is the key unrestricted-core dual
inequality. It neither assumes independence nor compares two distinct
cores.

## 3. Exact arbitrary-coupling Bellman majorant

For a future bit vector `xi in {0,1}^3`, root counting gives the exact
one-coordinate current marginals

```text
P(X_(j,t)=1 | X_(j,t+1)=xi_j)
 =(2-xi_j)/13.                                      (20)
```

The three current bits use the same root digit, so their joint law is not
independent and can depend on more than `xi`. Enlarge it pointwise to the
full coupling polytope

```text
K_xi={
 pi_eta>=0:
 sum_eta pi_eta=1,
 sum_(eta:eta_j=1)pi_eta=(2-xi_j)/13, j=1,2,3
}.                                                  (21)
```

At a terminal time, enlarge the actual joint law in the same way to every
law on `{0,1}^3` with all three marginals equal to `1/7`.

The companion runs a backward Bellman recurrence over:

```text
time t;
the labelled mask of even clauses in P_lambda already hit;
the total multiplicity of selected odd-charge atoms already hit;
the future bit vector xi.                           (22)
```

At every transition it maximizes over `K_xi`. At the terminal state it
returns exactly

```text
1_(P_lambda)(1-q_O)_+.                              (23)
```

Call the resulting rational number `B(lambda;O)`. The actual three-core
root process is feasible in every enlargement in (21), fibre by fibre.
Backward root induction and then the terminal Haar average therefore prove

```text
integral 1_(P_lambda)(1-q_O)_+
 <=B(lambda;O).                                     (24)
```

No Markov property for the *joint* three-core process is asserted. The
Bellman process is an adversarial relaxation retaining only the exact
single-core conditional marginals. Combining (4), (19), and (24), any
profile satisfying

```text
B(lambda;O)<961/6930                                (25)
```

is impossible.

## 4. Exact rational LP certificates

The coupling polytope in (21) has eight variables and four equality
constraints. Its columns are

```text
A_eta=(1,eta_1,eta_2,eta_3)^T.                      (26)
```

Exactly `58` of the `binom(8,4)` four-column bases are invertible. The
companion enumerates all nonnegative basic distributions, deduplicates
their vertices, and obtains the following vertex counts for the eight
future vectors in lexicographic bit order:

```text
(6,9,9,6,9,6,6,6).                                 (27)
```

The terminal marginal polytope has six vertices. Every Bellman maximum is
first computed over all exact rational vertices. For an optimal basis `I`,
the companion then solves

```text
A_I^T y=c_I                                         (28)
```

and verifies exactly

```text
A_eta^T y>=c_eta       for every eta,
b^T y=the primal optimum.                           (29)
```

Thus every maximum has an independent feasible dual certificate of the
same value. The full computation makes `45,617` LP calls and encounters
`39,295` distinct objective/right-side pairs; every distinct pair passes
the exact primal-dual check.

## 5. The complete newly closed census

THM-2233 leaves exactly 224 scalar profiles. Its 29 first-depth-two rows
are

```text
S_2
 ={(2,3,5)}
  union {(2,4,c):5<=c<=19}
  union {(2,b,b+2):5<=b<=17},                       (30)
```

and its 29 first-depth-three rows are

```text
S_3
 ={(3,3,5),(3,4,5),(3,4,6)}
  union {(3,5,c):6<=c<=19}
  union {(3,b,b+2):6<=b<=17}.                       (31)
```

For every row in `S_2 union S_3`, choose the single odd checkpoint

```text
O={1},

q=1-sum_(j=1)^3 Z_(j,lambda_j-1).                   (32)
```

The exact Bellman computation gives

```text
B(lambda;{1})<delta_5      for every lambda in S_2, (33)

B(lambda;{1})<delta_5
  for every lambda in S_3\{(3,4,5)}.                (34)
```

The largest bound in (33) occurs at `(2,4,5)`:

```text
B((2,4,5);{1})
 =8945166533/74231495611,

delta_5-B((2,4,5);{1})
 =1335209029783/73489180654890>0.                   (35)
```

The largest passing bound in (34) occurs at `(3,5,6)`:

```text
B((3,5,6);{1})
 =1471733046268/12545122758259,

delta_5-B((3,5,6);{1})
 =265250422864237/12419671530676410>0.              (36)
```

For the high row `(4,6,8)`, choose both available odd checkpoints:

```text
O={1,3},

q=2
  -Z_(1,3)-Z_(1,1)
  -Z_(2,5)-Z_(2,3)
  -Z_(3,7)-Z_(3,5).                                (37)
```

All six times in (37) are odd. The exact robust bound is

```text
B((4,6,8);{1,3})
 =17322925655936326/358301251098635299,

delta_5-B((4,6,8);{1,3})
 =32039946787164254737/354718238587648946010>0.     (38)
```

Equations (25), (33), (34), and (38) exclude all 29 rows in `S_2`,
28 rows in `S_3`, and `(4,6,8)`: exactly 58 profiles.

The six odd-time root corrections in (37), grouped by core, are

```text
(170/2197, 170/371293, 170/62748517).               (39)
```

Their sum is strictly below one. Thus the zero-charge terminal atom has
zero positive-part cost, including the worst choice of the three time-zero
bits. This is frozen as a boundary check in the companion.

The resulting exact scalar ledger is

```text
165 first-depth-one profiles
  (1,b,c),  5<=c<=19, 1<=b<c;

the single profile (3,4,5).                         (40)
```

It has `165+1=166` rows. This is a scalar-branch reduction, not a proof of
LRC(14).

## 6. Hostile and equality controls

The same carrier deliberately fails on the nearest surviving shallow row:

```text
B((3,4,5);{1})
 =17878637620/74231495611
 =delta_5
  +1072703906621/10498454379270.                    (41)
```

Equation (41) is a failure of this relaxed certificate, not a realizable
counterexample.

At first depth one, the only first-blocker charge available from checkpoint
one has time zero. Formula (14) then gives

```text
Z_(1,0)=1                                           (42)
```

identically. It certifies the odd clause without recording whether its
first-blocker literal fired, so the present nonnegative-charge carrier loses
its discriminating coordinate at exactly the remaining large family.

A tempting alternative is to scale the even-time raw centered atom by its
minimum value on the active set. At time two this leaves the inactive atom

```text
-rho^2/(1-rho^2)=-1/168<0.                          (43)
```

It therefore invalidates the clausewise implication used in (18).
The additive shift in (14), not this scaling, is the hostile-safe repair.

Finally, impose the common-core specialization `u_1=u_2=u_3` only as an
independent control for `(4,6,8)`. Direct exact Markov enumeration gives

```text
unsigned even-checkpoint mass =916159/4826809,

signed capacity
 =2303649491556761/51185893014090757.                (44)
```

The signed value in (44) is THM-2232's certificate and is no larger than
the unrestricted robust bound in (38), confirming the direction of the
arbitrary-coupling relaxation.

## 7. Connection and loss ledger

The proved connection is

```text
source:
  the signed five-unit residual, all even checkpoint supports, selected
  odd checkpoint supports, and the separate transfer eigen-correlations
  for each normalized blocker core;

target:
  a finite three-bit clause-and-charge Bellman capacity;

map:
  center each danger atom along rho=-1/13, add its least nonnegative
  parity shift, sum one unit of charge per selected odd clause, and use
  the resulting q<=0 certificate to discard the negative-support term;

preserved:
  every even clause label, every selected odd-clause multiplicity, exact
  profile times, the three individual time-zero correction bits, exact
  1-versus-2 root counts, and exact stationary marginals;

destroyed:
  the shared root digit, cross-core phase and incidence, guard and five-unit
  danger bits, owners/current, unselected odd clauses, and every joint
  restriction beyond the three one-coordinate root marginals;

cheapest hostile probes:
  (3,4,5), the time-zero identity Z_(1,0)=1, the invalid -1/168 scaled
  atom, and the exact common-core specialization;

needed sidecar:
  a first-depth-one literal/owner coordinate that survives time-zero
  centering, or cross-core incidence strong enough to sharpen (3,4,5).
                                                               (45)
```

The structural gain is that common-core alignment was never the essential
part of THM-2232. The essential object is a separately centered,
pointwise-positive charge whose clause sum has a fixed sign on the negative
residual. The precise remaining obstruction is now visible: at time zero
that charge becomes a constant and forgets the literal it was meant to
certify.

## 8. Reproduction

Run

```bash
python3 04-computation/lrc14_unrestricted_multicore_signed_dual_thm2239.py
python3 -O 04-computation/lrc14_unrestricted_multicore_signed_dual_thm2239.py
```

Both modes produce the stored transcript byte for byte. All load-bearing
checks use explicit exceptions rather than Python assertions. The companion
freezes the post-THM-2233 universe, every low-row and high-row bound, the
closed and remaining sets, SHA-256 digests of the full bound and profile
ledgers, the hostile values (41)--(44), every coupling-polytope census, and
exact primal-dual equality for every distinct LP. QED.
