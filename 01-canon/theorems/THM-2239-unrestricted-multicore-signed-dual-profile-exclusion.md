---
id: THM-2239
title: "Unrestricted multicore signed-dual profile exclusion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In the scalar five-unit
  plus three-blocker branch, a separately centered odd-checkpoint score and
  the four-bit arbitrary-coupling Bellman law close the sole remaining
  high-first profile (4,6,8), with no equality, common-core, or overlap
  assumption on the normalized blockers. The same score closes all 29
  profiles of first depth two left by THM-2233. The exact scalar ledger drops
  from 224 to 194: 165 profiles of first depth one and 29 of first depth
  three remain. The worst of the 30 new exclusions is (2,3,5), with capacity
  7265183507/74231495611 < 961/6930. LRC(14) remains open.
source: klein-2026-07-25-unrestricted-multicore-signed-dual
depends_on:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2233-guard-danger-hidden-state-bellman-profile-exclusion
related:
  - THM-2080-unequal-comb-overlap-removes-depth-five
  - THM-2232-same-core-signed-eigen-markov-dual-exclusion
script: 04-computation/lrc14_unrestricted_multicore_signed_dual_thm2239.py
output: 05-knowledge/results/lrc14_unrestricted_multicore_signed_dual_thm2239.out
script_sha256: e88978746022bdddb417ba3588337030bcb135456f465c10a07e405cd2808782
output_sha256: 7ca34e6f6215c7c71030e1bb8de24ce4f4dfbc77a50ef111fca21cbc35080774
hash_basis: working-tree bytes (LF)
---

# THM-2239 -- unrestricted multicore signed-dual exclusion

THM-2233 leaves `224` scalar valuation profiles, including one profile whose
first depth is at least four:

```text
(4,6,8).                                               (1)
```

The unsigned guard-hidden Bellman capacity of (1) is still above the scalar
residual floor. The missing coordinate is the sign of the original
five-unit residual. Keeping that sign through a separately centered score
does more than close (1): it removes every first-depth-two profile still in
the ledger.

## 1. Signed residual and checkpoint carriers

On `T=R/Z`, put

```text
D_b={x:||bx||<1/14},       C_H={x:||Hx||>1/7}.         (2)
```

Work in the scalar `5+3` branch of THM-2198. Thus the guard `H` and the five
unit coefficients `q_i` are thirteen-units, while

```text
c_j=13^(lambda_j)u_j,       13 does not divide u_j,
1<=lambda_1<=lambda_2<lambda_3<=19.                   (3)
```

No relation among `u_1,u_2,u_3` is assumed. Define

```text
R=1_(C_H)-sum_(i=1)^5 1_(D_(q_i)),
R_+=max(R,0),                 R_-=max(-R,0).           (4)
```

Since `R<=1`, its positive part is the indicator of

```text
A_+=C_H\union_(i=1)^5D_(q_i).                         (5)
```

THM-2198 and the zero integral of `R` give

```text
p:=integral R_+=integral R_- >=delta_5:=961/6930,
0<=R_+<=1,                       0<=R_-<=5.            (6)
```

Let `S(x)=13x mod 1` and write

```text
X_(j,t)(x)=1_(D_(u_j))(S^t x),
B_k(x)=sum_(j=1)^3 X_(j,lambda_j-k)(x)                (7)
```

for every `0<=k<=lambda_1`. Put

```text
E={0,2,...,2 floor(lambda_1/2)},
O={1,3,...,2 ceil(lambda_1/2)-1},                     (8)

P=C_H intersection intersection_(k in E){B_k>=1},
N=intersection_(k in O){B_k>=1}.                      (9)
```

THM-2222's transfer-parity tower says, almost everywhere,

```text
support(R_+) subset P,
support(R_-) subset N.                               (10)
```

The guard in `P` is not an extra hypothesis: positive residual already
means guard membership.

## 2. The multicore centered score

Put

```text
rho=-1/13,
Y_(j,t)=X_(j,t)-rho^t X_(j,0).                       (11)
```

The unnormalized transfer operator satisfies `L^tR=(-1)^tR`. Transfer
duality therefore gives, separately for every core `u_j`,

```text
integral R X_(j,t)
 =rho^t integral R X_(j,0),
integral R Y_(j,t)=0.                                (12)
```

Let `o=|O|` and define the canonical odd-checkpoint score

```text
q=o-sum_(k in O)sum_(j=1)^3
       Y_(j,lambda_j-k).                             (13)
```

Equations (6), (12), and `integral R=0` imply

```text
integral Rq=0.                                       (14)
```

This is the decisive difference from the common-core dual of THM-2232:
every normalized blocker has its own centered eigenline. Equation (14)
therefore survives arbitrary cross-core relations, coincidences, and
root-digit couplings.

## 3. A sign-safe pointwise majorant

From (6) and (14),

```text
p=integral [R_+(1-q)+R_-q].                          (15)
```

Using (10) pointwise, without assuming that `q` lies in `[0,1]`,

```text
p<=integral [
     1_P(1-q)_+ + 5 1_N q_+
   ].                                                (16)
```

Indeed, if `1-q<0`, then `R_+(1-q)<=0`; otherwise use `R_+<=1_P`.
Likewise, if `q<0`, then `R_-q<=0`; otherwise use `R_-<=5 1_N`.
Thus no absolute-value or termwise sign estimate enters (16).

For the high-first profile (1), `O={1,3}` and all six times
`lambda_j-k` are odd. Hence

```text
Y_(j,t)=X_(j,t)+13^(-t)X_(j,0)>=X_(j,t).             (17)
```

On `N`, each of `B_1,B_3` is at least one, so

```text
q=2-sum_(k in {1,3})sum_jY_(j,lambda_j-k)
 <=2-B_1-B_3<=0.                                    (18)
```

The negative term in (16) is therefore identically zero for (1). This is
an algebraic sign consequence, not a numerical feature of the optimizer.

## 4. Exact four-bit score Bellman recursion

Retain the simultaneous hidden state

```text
Z_t=(G_t,X_(1,t),X_(2,t),X_(3,t)),
G_t=1_(C_H)(S^t x).                                  (19)
```

For a fixed future state

```text
xi=(g,x_1,x_2,x_3) in {0,1}^4,                       (20)
```

the exact thirteen-root identities give the current-coordinate marginals

```text
P(G_current=1 | xi)=(10-g)/13,
P(X_(j,current)=1 | xi)=(2-x_j)/13.                  (21)
```

The three danger bits and guard bit share one root digit and need not be
independent. Let `C(xi)` be the polytope of **all** joint laws on
`{0,1}^4` having the four marginals (21). Every actual common-root law is
feasible, so maximizing over `C(xi)` is a safe relaxation for arbitrary
unit cores.

The Bellman state retains:

```text
a set A of already hit guard/even/odd clauses;
an exact rational score s accumulated from (13);
the future four-bit state xi.                        (22)
```

At unit time `t`, a bit vector `z` hits the guard clause when `t=0` and its
guard bit is one, and hits checkpoint clause `k` when

```text
t=lambda_j-k and z_j=1                               (23)
```

for some `j`. Let `H_t(z)` be this clause mask. The score increment is

```text
d_t(z)
 =-sum_({(k,j):k in O, lambda_j-k=t}) z_j
  +1_(t=0) sum_(k in O)sum_j
       rho^(lambda_j-k) z_j.                         (24)
```

The two terms in (24) cancel automatically when `lambda_j-k=0`.
Starting with

```text
V_(-1)(A,s,xi)
 =1_(all P clauses in A)(1-o-s)_+
  +5 1_(all N clauses in A)(o+s)_+,                  (25)
```

define

```text
V_t(A,s,xi)
 =max_(pi in C(xi)) sum_z pi_z
    V_(t-1)(A union H_t(z),s+d_t(z),z).              (26)
```

After time `lambda_3`, maximize over every terminal joint bit law with
stationary marginals

```text
(5/7,1/7,1/7,1/7).                                  (27)
```

Backward root induction, exactly as in THM-2233, shows that this terminal
value bounds the right side of (16). The accumulator in (22) preserves the
signed affine observable that an ordinary clause-only Bellman state loses.

Every coupling polytope has sixteen nonnegative atom variables and five
independent equalities. The companion enumerates all exact basic feasible
laws and, for every distinct maximization in (26), constructs an exact
dual-feasible basic certificate with the same value.

## 5. Exact profile consequence

THM-2233's exact `224`-profile residue is

```text
all 165 profiles with lambda_1=1 and lambda_3>=5;

(2,3,5);
(2,4,c), 5<=c<=19;
(2,b,b+2), 5<=b<=17;

(3,3,5), (3,4,5), (3,4,6);
(3,5,c), 6<=c<=19;
(3,b,b+2), 6<=b<=17;

(4,6,8).                                             (28)
```

The exact recurrence (25)--(27) closes the following `30` profiles:

```text
(2,3,5);
(2,4,c), 5<=c<=19;
(2,b,b+2), 5<=b<=17;
(4,6,8).                                             (29)
```

The worst capacity among (29) is

```text
B(2,3,5)
 =7265183507/74231495611
 =0.097871980716543...,

delta_5-B(2,3,5)
 =2998392225523/73489180654890>0.                    (30)
```

For the sole high-first row,

```text
B(4,6,8)
 =13063069772240011/358301251098635299
 =0.036458342615845...,

delta_5-B(4,6,8)
 =36257204112023606587/354718238587648946010>0.      (31)
```

Combining (6), (16), and the strict bounds (30)--(31) proves that every
profile in (29) is empty. The current scalar ledger becomes

```text
224-30=194,                                          (32)
```

consisting exactly of the `165` first-depth-one rows and the following `29`
first-depth-three rows:

```text
(3,3,5), (3,4,5), (3,4,6);
(3,5,c), 6<=c<=19;
(3,b,b+2), 6<=b<=17.                                (33)
```

In particular, every surviving scalar profile now has first depth either
one or three. This theorem closes the unrestricted high-first branch, but
does not close the scalar lane or LRC(14).

## 6. Failure boundary of the canonical score

The all-ones score (13) is not a universal closure. As an exact hostile
control, its capacity at the surviving profile `(3,17,19)` is

```text
229985997674091853361506127343787104565593
--------------------------------------------------
1150805888298461593768523506367492533368331

 =0.199847776251945... > delta_5.                    (34)
```

The exact excess is

```text
69696929318090707454812338981897758581799057
---------------------------------------------------
1139297829415476977830838271303817608034647690 >0.  (35)
```

Thus the next depth-three step must change the affine weights, add another
signed observable, or retain more realized root information. Merely rerunning
the canonical odd-checkpoint score cannot be presented as a complete scalar
closure.

THM-2080's guard-danger pair cap is compatible with the four-bit state, but
is not used in (30)--(31). The proved mechanism is stronger and cleaner:
the exact `-1/13` transfer phase kills the negative carrier while the robust
hidden-state recursion prices the remaining positive carrier.

## 7. Connection and loss ledger

```text
source:
  the original signed five-unit residual and every transfer-parity
  checkpoint;

map:
  center each blocker core separately on the -1/13 eigenline, sum the
  odd-checkpoint atoms into q, and attach q as a Bellman score;

preserved:
  residual sign, guard membership, all even and odd checkpoint clauses,
  exact one-coordinate root laws, and arbitrary cross-core coupling;

destroyed:
  the actual joint root pattern inside each marginal polytope, owner labels,
  arithmetic relations among distinct cores, and all score directions
  orthogonal to the chosen affine q;

needed sidecar for (33):
  optimized checkpoint weights, a complementary even-checkpoint score,
  a second transfer eigen-observable, or a realizability restriction on
  the relaxed root couplings.                                (36)
```

## 8. Exact and independent audit

Run

```bash
python3 04-computation/lrc14_unrestricted_multicore_signed_dual_thm2239.py
python3 -O 04-computation/lrc14_unrestricted_multicore_signed_dual_thm2239.py
```

The exact companion checks:

```text
the 224-profile input residue and all set digests;
all 30 claimed strict exclusions and their exact worst row;
the pointwise-zero negative term for (4,6,8);
the exact 194-profile output residue;
one exact surviving hostile control;
61,748 score-Bellman states;
54,483 local LP calls;
44,631 distinct exact primal/dual certificates;
ordinary/optimized/stored transcript identity.       (37)
```

An independent implementation rebuilt every local coupling LP directly
with SciPy/HiGHS rather than using the exact vertex enumerator. For (1) it
returned

```text
0.036458342615845143,
```

within `1.4e-17` of (31), and independently audited the schedule, signs,
support inclusions, terminal law, and arbitrary-core quantifiers. QED.
