---
id: THM-2243
title: "Composite-union transfer dual for depth-three scalar profiles"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED. In the scalar five-unit plus
  three-blocker branch, the four depth-three checkpoint unions are shifts
  of one composite Boolean danger event. Centering that whole union across
  the two-step transfer eigenline gives a nonlinear signed score which
  vanishes pointwise on the negative carrier. The exact four-bit
  arbitrary-coupling Bellman recursion independently excludes 28 of the 29
  depth-three profiles in THM-2233's 224-row residue. The sole hostile row
  is (3,4,5), with bound 8907541/62377224 > 961/6930. These 28 rows overlap
  the independently proved THM-2239 and THM-2244 exclusions and are not
  counted again; the scalar ledger remains 166. LRC(14) remains open.
source: codex-2026-07-25-composite-union-transfer-dual
depends_on:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2233-guard-danger-hidden-state-bellman-profile-exclusion
related:
  - THM-2237-truncated-boolean-moment-interval-and-parity-top-atom-majorants
  - THM-2239-unrestricted-multicore-signed-dual-profile-exclusion
  - THM-2244-prior-centered-single-odd-clause-bellman-depth-three-exclusion
script: 04-computation/lrc14_composite_union_transfer_dual_thm2243.py
output: 05-knowledge/results/lrc14_composite_union_transfer_dual_thm2243.out
script_sha256: 756b85f6dae672ecaf8ce6f011f5d8c61053a191f90adb14b195857369c6e40c
output_sha256: b7cc692350a2e330429ad6510aaf40b28c5d475d48d31749c01b9c9e1265c0cb
hash_basis: working-tree bytes (LF)
---

# THM-2243 -- composite-union transfer dual

THM-2233 leaves the exact scalar ledger

```text
224 = 165 first-depth-one + 29 first-depth-two
      + 29 first-depth-three + (4,6,8).                    (1)
```

The present theorem gives an independent treatment of the 29 depth-three
rows: take the union of the three blocker atoms *before* applying transfer
duality. This nonlinear order of operations closes all but `(3,4,5)`.
THM-2239 and THM-2244 now close the same 28 rows by different atomwise
scores, so the point here is mechanism separation, a sharper hostile bound,
and an exact cross-check rather than another ledger decrement.

## 1. The four checkpoint unions form one orbit

On `T=R/Z`, put

```text
D_b={x:||bx||<1/14},       C_H={x:||Hx||>1/7}.       (2)
```

Work in the scalar `5+3` branch of THM-2198. Thus

```text
R=1_(C_H)-sum_(i=1)^5 1_(D_(q_i)),

c_j=13^(lambda_j)u_j,       13 does not divide u_j,
3=lambda_1<=lambda_2<lambda_3<=19.                  (3)
```

The guard, five unit coefficients, and normalized blocker cores are
thirteen-units. No relation among `u_1,u_2,u_3` is assumed. As in THM-2239,

```text
p:=integral R_+=integral R_- >=delta_5:=961/6930,
0<=R_+<=1,                       0<=R_-<=5.           (4)
```

Let `Sx=13x mod 1`, and define

```text
X_(j,t)=1_(D_(u_j)) o S^t,
B_k=sum_(j=1)^3 X_(j,lambda_j-k),
U_k=1_{B_k>=1},                         0<=k<=3.     (5)
```

Now form the single composite event

```text
U
 =1_(union_(j=1)^3 D_(13^(lambda_j-3)u_j)).          (6)
```

For every checkpoint `k`,

```text
X_(j,lambda_j-k)
 =1_(D_(13^(lambda_j-3)u_j)) o S^(3-k),

U_k=U o S^(3-k).                                    (7)
```

Thus the four clauses are not four unrelated unions: they are four points
on one transfer orbit. This identity is exact for arbitrary normalized
cores, including coincidences and arithmetic relations among them.

The depth-three support inclusions from THM-2222 are

```text
support(R_+) subset
 P:=C_H intersection {U_0=1} intersection {U_2=1},

support(R_-) subset
 N:={U_1=1} intersection {U_3=1}.                   (8)
```

## 2. Center the nonlinear union, not its atoms

Put

```text
rho=-1/13,                    r=rho^2=1/169.         (9)
```

THM-2222 proves the unnormalized transfer eigenidentity

```text
L^tR=(-1)^tR.                                       (10)
```

Transfer duality therefore gives, for every bounded measurable `f`,

```text
integral R(f o S^t)=rho^t integral Rf.              (11)
```

This universal quantifier is load-bearing: apply (11) to the composite
Boolean function `U`, not separately to the three atoms. Equations (7) and
(11) yield

```text
integral R(U_1-rU_3)=0.                             (12)
```

Define

```text
q=1-(U_1-rU_3)/(1-r).                              (13)
```

Since `integral R=0`, equation (12) implies

```text
integral Rq=0.                                      (14)
```

On the negative carrier `N`, both `U_1` and `U_3` equal one, so

```text
q=0                         pointwise on N.          (15)
```

This is stronger than merely making the negative contribution small.

## 3. The sign-safe pointwise reduction

Equations (4) and (14) give

```text
p=integral [R_+(1-q)+R_-q].                         (16)
```

By (8) and (15), the second term is identically zero. The first term obeys
the sign-safe bound

```text
p<=integral 1_P(1-q)_+.                             (17)
```

The complete two-bit truth table of `(U_1,U_3)` gives the exact hinge

```text
(1-q)_+
 =[(U_1-rU_3)/(1-r)]_+
 =U_1(169-U_3)/168.                                 (18)
```

Consequently

```text
p<=E[
  1_(C_H) U_0 U_2 U_1(169-U_3)/168
].                                                  (19)
```

No absolute-value estimate, independence assumption, common-core
hypothesis, or termwise replacement of the union appears in (19).

## 4. What information the composite event retains

At a fixed time, write `X_j=X_(j,t)`. Boolean inclusion-exclusion gives

```text
U
 =X_1+X_2+X_3
  -X_1X_2-X_1X_3-X_2X_3
  +X_1X_2X_3.                                      (20)
```

THM-2239's score is degree one in the labelled atoms. Equation (20) retains
the pair intersections and the cubic top atom. This is exactly the kind of
coordinate that THM-2237 identifies as invisible to a packet of proper
Boolean marginals: parity fibres can agree below top degree and split only
at the top Walsh coordinate. The root-imported
`GMC2ParityPolarization.lean` calculation supplies the same explanatory
control—quadratic data agrees on its parity fibres while the cubic monomial
separates them. Neither result is a proof dependency here; the actual proof
is the all-functions transfer identity (11).

## 5. Exact four-bit Bellman evaluation

Retain the hidden state

```text
Z_t=(G_t,X_(1,t),X_(2,t),X_(3,t)),
G_t=1_(C_H) o S^t.                                  (21)
```

For a fixed future state `xi=(g,x_1,x_2,x_3)`, exact root counting gives
the current-coordinate marginals

```text
P(G_current=1 | xi)=(10-g)/13,
P(X_(j,current)=1 | xi)=(2-x_j)/13.                 (22)
```

Let `C(xi)` be the polytope of all laws on `{0,1}^4` with the four
marginals (22). Every actual common-root law is feasible. Maximizing over
the whole polytope forgets joint root incidence but is therefore a safe
upper bound for arbitrary cores.

The automaton stores which of the five clauses

```text
guard, U_0,U_1,U_2,U_3                               (23)
```

have been hit. At unit time `t`, the bit `X_(j,t)` hits clause `U_k`
exactly when

```text
t=lambda_j-k.                                       (24)
```

After all times have been processed, the terminal payoff is precisely

```text
1_(guard,U_0,U_2 all hit)
  1_(U_1 hit) [169-1_(U_3 hit)]/168.                (25)
```

Backward root induction, with (22) at each step and terminal marginals

```text
(5/7,1/7,1/7,1/7),                                  (26)
```

proves that the resulting Bellman value bounds (19). The optimizer may
choose a different arbitrary coupling for every future state and every
clause mask; this enlarges the actual arithmetic process and does not impose
a Markov or independence model on the cores.

Each local polytope has sixteen nonnegative variables and five independent
equalities. The companion enumerates its exact basic feasible laws and, for
every distinct objective, constructs an exact dual-feasible basic
certificate with the same value.

## 6. Exact profile consequence and overlap ledger

The 29 depth-three rows in THM-2233's residue are

```text
(3,3,5), (3,4,5), (3,4,6);
(3,5,c), 6<=c<=19;
(3,b,b+2), 6<=b<=17.                                (27)
```

The exact recurrence closes the following 28:

```text
(3,3,5), (3,4,6);
(3,5,c), 6<=c<=19;
(3,b,b+2), 6<=b<=17.                                (28)
```

The largest bound among (28) occurs at `(3,4,6)`:

```text
B(3,4,6)
 =18956451/270301304
 =0.070130815943085...,

delta_5-B(3,4,6)
 =1310115793/19114163640
 >0.                                                (29)
```

Applied by itself to THM-2233's residue, (29) leaves

```text
224-28=196                                          (30)
```

profiles: all 165 first-depth-one rows, all 29 first-depth-two rows,
`(3,4,5)`, and `(4,6,8)`.

THM-2239 and THM-2244 independently exclude exactly the same 28
depth-three profiles by different atomwise centering mechanisms. The set in
(28) is therefore an overlap control and must not be subtracted again. In
the current proof graph, THM-2239 already gives the canonical `166` ledger;
THM-2243 sharpens the mechanism and quantitative hostile boundary without
reducing that ledger further.

## 7. Exact failure boundary

For the sole survivor, the composite-union Bellman gives

```text
B(3,4,5)
 =8907541/62377224
 =0.142801176916754...,

B(3,4,5)-delta_5
 =42493973/10292241960
 >0.                                                (31)
```

Thus the score (13) and the arbitrary-coupling state (21) do not by
themselves close `(3,4,5)`. This hostile value is nevertheless much sharper
than THM-2244's best prior-centered atom bound on the same row. The remaining
gap is about `0.004129`; the next step must retain information discarded by
`C(xi)`, such as realized guard/owner incidence or within-union labelled
multiplicity, rather than counting (28) again.

The construction also has no first-depth-one analogue of (12): there is
only one odd checkpoint, so there is no same-parity pair of odd union
clauses separated by two transfer steps. First depth one remains routed to
THM-2234's private two-owner expansion and the carry sidecars.

## 8. Typed connection and loss ledger

```text
source:
  the four transfer-parity checkpoint unions and the signed scalar residual;

target:
  one nonlinear transfer-orbit score attached to the four-bit Bellman state;

map:
  identify U_k=U o S^(3-k), center U_1 against U_3 on the rho^2
  eigenline, and price the exact hinge on the positive carrier;

preserved:
  residual sign, guard membership, all four union clauses, pair/triple
  Boolean intersections, the cubic top atom, exact times, and exact
  one-coordinate root laws;

destroyed:
  which labelled blocker owns a hit, hit multiplicity inside U_k, the
  actual common root digit, arithmetic relations among cores, and joint
  incidence beyond marginal coupling polytopes;

cheapest decisive continuation test:
  add one realized guard/owner or within-union multiplicity coordinate to
  the (3,4,5) Bellman and demand a strict bound below 961/6930.           (32)
```

## 9. Exact audit and reproduction

Run

```bash
python3 04-computation/lrc14_composite_union_transfer_dual_thm2243.py
python3 -O 04-computation/lrc14_composite_union_transfer_dual_thm2243.py
```

The two modes reproduce the stored LF transcript byte for byte after the
platform's newline convention is normalized. The audit checks

```text
the exact 224-profile post-THM-2233 input and 29-row depth-three census;
the complete U_0,...,U_3 score truth table;
all 28 strict exclusions and the hostile (3,4,5) row;
the exact 196-row standalone consequence, the overlapping 166-row
canonical consequence, and all profile/bound digests;
195,584 Bellman states and 195,613 local LP calls;
65,437 distinct exact primal/dual certificates;
ordinary/optimized/stored LF-normalized transcript identity. (33)
```

This theorem concerns only the scalar valuation ledger. It does not close
the first-depth-one rows, the exceptional `(3,4,5)` row, any owner/current
branch, or LRC(14). QED.
