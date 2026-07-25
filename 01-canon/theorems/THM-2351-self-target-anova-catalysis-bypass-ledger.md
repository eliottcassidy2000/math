---
id: THM-2351
title: "Self-target allocation ANOVA catalysis--bypass ledger"
status: >
  CLAIMED + VERIFIED-EXACT, PROOF CANDIDATE UNDER INDEPENDENT AUDIT. For
  two nontrivial prime knots P,Q and X=P#Q, the complete two-token
  allocation energy of the canonical source packet `(U,X)` has a unique
  self-target optimum. Its zero-marginal pair coefficient is exactly minus
  one quarter of the sum of the two directional bypass terms. Together
  with the root defect it detects whether either summand catalyses the other,
  while corrected unary fields recover both directional catalytic terms
  separately. The two continuation-row Mobius coefficients recover the same
  catalysis--bypass split. Thus the full self-target ANOVA is an affine
  reparameterization of the existing directional pair ledger, not a
  stronger invariant; the optimum and the skew pair sector are blind. The
  Brittenham--Hermiller mirror pair lies sharply on the pure-bypass boundary.
  No new Gordian distance, positive catalyst, unknotting number, stable
  diagonal, or knot classification is asserted.
source: codex-2026-07-25-self-target-anova-ledger
depends_on:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2339-prime-owner-deletion-and-target-allocation-hypergraph
  - THM-2346-global-allocation-anova-normal-form
  - THM-2348-prime-type-rectangularity-and-target-token-conditioning
related:
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2294-anchored-plucker-tournament-and-kakeya-address-bank
  - THM-2336-prime-target-gordian-owner-diagram-and-bypass-split
script: 04-computation/self_target_anova_catalysis_bypass_thm2351.py
output: 05-knowledge/results/self_target_anova_catalysis_bypass_thm2351.out
script_sha256: 080751d05d71027a9ef847fd6bb9031dee35c64709a2c44909ffe9e74b6724db
output_sha256: 40dc845b7dbf8562dee38992d0d5eaf9cfee056bc0fc3ad51c5500f0ff9a56af
hash_basis: working-tree bytes (LF)
---

# THM-2351 -- the self-target ANOVA catalysis--bypass ledger

**CLAIMED + VERIFIED-EXACT; PROOF CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2346 turns a composite-target allocation energy into unique global
zero-marginal tensors. THM-2348 shows that its two-token cohabitation
coefficient is a continuation of the connected-sum defect. The most
canonical place to evaluate that structure is the source packet consisting
of the unknot and the target itself.

At that self-target packet, the allocation optimum is automatic and
therefore uninformative. The full ANOVA is not. Its symmetric pair field
measures total geodesic bypass, and its corrected unary imbalance measures
directional translation catalysis. This gives an exact catalyst test for the
two displayed summands, but the resulting data are linearly equivalent to
THM-2176's existing pair ledger. The theorem is therefore both a positive
connection and a stopping boundary.

## 1. The directional pair ledger

Let `P,Q` be nontrivial oriented prime knots and put

```text
X=P#Q.
```

For the main allocation statement assume `P,Q` are nonisomorphic, so the
two target-prime occurrences have distinct types. Section 6 records the
equal-prime quotient.

Write

```text
a=u(P),                 b=u(Q),                 L=u(X),

sigma=a+b-L.                                             (1)
```

Use THM-2176's directional translation and bypass terms

```text
c_P=C_Q(P)=a-d_G(X,Q),

c_Q=C_P(Q)=b-d_G(X,P),

b_P=B_Q(P)=sigma-c_P,

b_Q=B_P(Q)=sigma-c_Q.                                  (2)
```

The subscripts on `c_P,b_P` name the summand whose root cost is being
transported; the other prime is its context. Every quantity in (1)--(2) is
nonnegative, and

```text
0<=c_P,c_Q,b_P,b_Q<=sigma,

c_P+b_P=c_Q+b_Q=sigma.                                (3)
```

No exact value of any unknown knot distance is assumed.

## 2. The canonical self-target allocation table

Take the two-block source packet

```text
(U,X)
```

and the composite target `P#Q=X`. Give a target token colour `0` when it is
assigned to the source block `U`, and colour `1` when it is assigned to the
source block `X`.

THM-2339's normalized allocation energy subtracts the block-root baseline

```text
u(U)+u(X)=L.
```

Index rows by the colour of the `P` token and columns by the colour of the
`Q` token.

> **Self-target table.**
>
> ```text
> E=
>   [ L                 2a-c_P-L ]
>   [ 2b-c_Q-L              -L  ].                    (4)
> ```
>
> Before baseline subtraction, the four lift costs are
>
> ```text
>   [ 2L                2a-c_P ]
>   [ 2b-c_Q                 0  ].                    (5)
> ```

### Proof

If both tokens go to the `U` block, the endpoint is `(X,U)`, at product
distance

```text
d_G(U,X)+d_G(X,U)=2L.
```

If `P` goes to `U` and `Q` goes to `X`, the endpoint is `(P,Q)`, at cost

```text
u(P)+d_G(X,Q)=a+(a-c_P)=2a-c_P.
```

The other split has cost `2b-c_Q`. Assigning both tokens to `X` leaves the
source packet `(U,X)` unchanged and costs zero. Subtract `L` to obtain
(4). QED.

Because `P,Q` are nontrivial,

```text
L>0,

2a-c_P>=a>0,

2b-c_Q>=b>0.                                          (6)
```

Therefore

```text
Opt_(U,X)(P#Q)={(1,1)}.                               (7)
```

The tied optimum is the same for every possible directional split. In
particular, this composite self-owner is a sharp instance of THM-2348's
warning that one Cartesian optimum does not imply decoupling.

THM-2339's prime-owner deletion law does not apply to (7): the distinguished
block is the **composite** target `P#Q`, not a prime target. Fixing both
target tokens at that block is instead THM-2348's labelled conditioning
operation.

## 3. Exact ANOVA coordinates

Put

```text
h(0)=1,                    h(1)=-1.
```

Write THM-2346's unique uniform two-token ANOVA expansion as

```text
E(i,j)
 =m+r h(i)+s h(j)+t h(i)h(j).                        (8)
```

Here `m` is the mean, `r,s` are the `P,Q` unary coefficients, and `t` is
the symmetric Ising/Potts pair coefficient.

> **Catalysis--bypass ANOVA ledger.**
>
> ```text
> m=-t=(b_P+b_Q)/4,                                  (9)
>
> r+s=L,                                             (10)
>
> r-s=a-b+(c_Q-c_P)/2,                               (11)
>
> t=-(b_P+b_Q)/4.                                    (12)
> ```

### Proof

For a `2x2` table, the four coefficients are

```text
m=(E_00+E_01+E_10+E_11)/4,

r=(E_00+E_01-E_10-E_11)/4,

s=(E_00-E_01+E_10-E_11)/4,

t=(E_00-E_01-E_10+E_11)/4.                          (13)
```

Substitute (4). Since

```text
b_P+b_Q
 =2sigma-c_P-c_Q
 =2(a+b-L)-c_P-c_Q,                                 (14)
```

the mean is `(b_P+b_Q)/4` and the pair coefficient is its negative.
The two remaining sums simplify to (10)--(11). QED.

Thus the allocation pair tensor is

```text
H_{P,Q}=t h tensor h,                               (15)
```

with

```text
-sigma/2<=t<=0.                                     (16)
```

It is always symmetric and nonpositive. In statistical-mechanical
language it is a ferromagnetic cohabitation term. It is never the skew
observable of a tournament on the two source blocks.

The boundary cases are exact:

```text
t=-sigma/2
 iff c_P=c_Q=0
 iff neither P nor Q catalyses the other;

t=0
 iff b_P=b_Q=0
 iff both directional interactions are pure translation;

-sigma/2<t<0
 iff the pair has both some translation saving and some bypass
     after the two directions are combined.                         (17)
```

In particular,

```text
positive translation catalysis in at least one direction
 iff t>-sigma/2.                                     (18)
```

This test requires the root sidecar `sigma`. The pair tensor alone measures
total bypass, not total catalysis.

## 4. Recovering both directions and the tournament shadow

The root lengths and the full ANOVA recover the directional terms:

```text
c_P
 =sigma+2t-r+s+a-b,

c_Q
 =sigma+2t+r-s-a+b,                                 (19)

b_P
 =-2t+r-s-a+b,

b_Q
 =-2t-r+s+a-b.                                      (20)
```

Indeed, (9)--(12) give

```text
c_P+c_Q=2sigma+4t,

c_Q-c_P=2[r-s-(a-b)],                               (21)
```

and (19)--(20) follow by solving the two linear equations and then using
(2).

The directional catalytic imbalance is therefore

```text
c_P-c_Q
 =2[(a-b)-(r-s)].                                   (22)
```

Its sign is THM-2294's tie-aware catalytic tournament shadow on the knot
pair. Equation (22) shows where that antisymmetry lives in allocation
coordinates: it is a root-length-corrected **unary** contrast. The pair
tensor (15) has zero alternating sector.

Under `P<->Q`,

```text
r<->s,          c_P<->c_Q,          t->t.
```

Thus `(a-b)-(r-s)` is skew while `t` is symmetric. A tournament retaining
only the sign of (22) discards both the common positive catalytic saving
and the bypass magnitude.

The inequalities (3) become the exact ANOVA diamond

```text
|r-s-(a-b)|
 <=min(sigma+2t,-2t),                               (23)
```

together with (16). This is just the four directional nonnegativity
conditions in new coordinates; no additional knot inequality is claimed.

## 5. The two continuation rows

Use THM-2348's target-row cohabitation coefficient

```text
mu_K(P,Q)
 =d_G(K,P#Q)-d_G(K,P)-d_G(K,Q)+u(K).                (24)
```

At the two canonical source rows,

```text
mu_U(P,Q)=-sigma,                                   (25)

mu_X(P,Q)=-sigma+c_P+c_Q.                           (26)
```

Consequently

```text
mu_X-mu_U=c_P+c_Q=2sigma+4t,                        (27)

-mu_U-mu_X=b_P+b_Q=-4t,                             (28)

t=[mu_U+mu_X]/4.                                    (29)
```

Equation (27) is an exact symmetric catalysis detector:

```text
mu_X>mu_U
 iff c_P+c_Q>0
 iff one of P,Q catalyses the other.                (30)
```

Equation (28) gives the complementary total bypass. The interval

```text
mu_U<=mu_X<=-mu_U                                   (31)
```

is equivalent to (3). Its lower equality is pure bypass in both
directions, and its upper equality is zero bypass in both directions.

These identities hold abstractly in any commutative nonexpansive metric
monoid if (24) is read as the four-point target-row difference and
`X=P+Q`. The prime-knot hypotheses are needed only to identify the four
allocations in (4) with the complete THM-2339 target fibre.

## 6. Equal-prime quotient

If `P=Q`, use two labelled occurrences of the same oriented prime and
then quotient by their transposition as in THM-2339 and THM-2346. The two
split allocations in (4) have the same endpoint and the same cost. In the
labelled ANOVA,

```text
a=b,              c_P=c_Q,              r=s.         (32)
```

The corrected unary imbalance (22) vanishes, as it must for
indistinguishable tokens, while

```text
t=-B_P(P)/2                                       (33)
```

still records the common bypass. The equal-prime quotient therefore keeps
the symmetric Potts interaction and correctly kills the artificial
direction.

## 7. Brittenham--Hermiller calibration

Let

```text
P=T(2,7),                    Q=mirror(P),

X=P#Q.
```

THM-2176 proves

```text
a=b=3,

c_P=c_Q=0,

sigma=6-u(X)>=1,                                   (34)
```

without computing `u(X)`. The normalized self-target table and ANOVA are

```text
E=
  [ u(X)       sigma ]
  [ sigma      -u(X) ],                             (35)

m=sigma/2,

r=s=u(X)/2,

t=-sigma/2.                                        (36)
```

The raw lift-cost table is

```text
  [ 2u(X)      6 ]
  [ 6          0 ].                                (37)
```

Thus the first nonadditive connected-sum pair lies sharply on the lower
boundary in (16): its pair ANOVA is entirely bypass, its catalytic
tournament is tied, and the composite self-owner remains unique. No exact
value of `u(X)` is inferred.

## 8. Exact stopping theorem

Fix the two root lengths `a,b`. Equations (1), (9)--(12), and (19) give
mutually inverse affine maps

```text
(sigma,c_P,c_Q)
 <->
 (m,r,s,t) subject to m=-t and r+s=a+b-sigma.       (38)
```

Therefore:

> **Self-target stopping theorem.** The full self-target allocation table,
> its ANOVA normal form, and THM-2176's directional pair ledger contain
> exactly the same numerical information once `u(P),u(Q)` are retained.
> The ANOVA is a canonical realization of that ledger, not a strictly
> stronger pair invariant.

Three common shadows are strictly weaker.

- The optimum (7) is constant and loses the whole ledger.
- The pair tensor (15) retains only `b_P+b_Q`.
- The tournament sign from (22) retains only which directional catalytic
  term is larger and keeps genuine ties.

Accordingly, a new obstruction to catalysis by a third knot, a value of the
catalytic localization `u_cat`, or the stable connected-sum diagonal must
use another continuation row, target, scale, or external calibrator. It
cannot come from re-expressing this one self-target table again.

## 9. Exact companion and scope

The optimization-safe exact companion:

- reconstructs every table and ANOVA coefficient over a finite exact bank
  of algebraically admissible directional ledgers;
- verifies (9)--(31), all boundary iff statements, and the affine inverse;
- checks the unique self-target optimum throughout the bank;
- checks equal-prime symmetry and quotient compatibility;
- checks the Brittenham--Hermiller formulas for every unresolved integral
  value allowed by its known upper bound.

Reproduce with

```bash
python3 04-computation/self_target_anova_catalysis_bypass_thm2351.py
python3 -O 04-computation/self_target_anova_catalysis_bypass_thm2351.py
```

and compare both transcripts with

```text
05-knowledge/results/self_target_anova_catalysis_bypass_thm2351.out
```

The finite ledger bank is a control for the identities, not an assertion
that every algebraically admissible row is realized by knots. The theorem
computes no new Gordian distance, positive catalyst, unknotting number,
stable diagonal, or knot classification.
