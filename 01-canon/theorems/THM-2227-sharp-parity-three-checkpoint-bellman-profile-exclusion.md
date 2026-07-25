---
id: THM-2227
title: "Half-time hidden-parity Bellman and six-profile exclusion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED + ALTERNATE CARRIER. For the eight
  three-checkpoint profiles left by THM-2226, retain the three unit-core
  danger bits at every multiplication-by-13 half-time. The exact one-core
  backward marginal is (2-X_next)/13; optimizing over every joint current-bit
  coupling with these marginals gives a safe hidden-state Bellman relaxation.
  It closes all six profiles having an odd relative valuation and leaves only
  (4,6,8) and (5,7,9), both with relative valuation (0,2,4). Exact primal and
  dual vertex enumerations agree, ordinary and optimized outputs are
  byte-identical, and an independent floating-point implementation reproduces
  every bound. Relative to the THM-2219/2226 ledger alone this would reduce
  458 to 452; THM-2229 independently subsumes the same six exclusions and
  closes additional profiles, so this is not an additive current-ledger
  decrement. LRC(14) remains open.
source: codex-2026-07-24-half-time-hidden-parity-bellman
depends_on:
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2219-scalar-depth-four-sparse-tail-exclusion
  - THM-2226-three-checkpoint-bellman-sieve-and-eight-profile-residue
related:
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2224-transfer-owner-word-temporal-union-bound
  - THM-2229-unit-time-positive-set-bellman-profile-exclusion
script: 04-computation/lrc14_half_time_hidden_parity_bellman_thm2227.py
output: 05-knowledge/results/lrc14_half_time_hidden_parity_bellman_thm2227.out
script_sha256: 86a8eb3557adc5d1613ba7604cc78f8e3875d42ea85b913aff7d8cdb44957d6e
output_sha256: 268520822faeca8ed2babe43ab85abf0d842f1a4cc033c28619104ecc68f4ae8
hash_basis: working-tree bytes (LF)
---

# THM-2227 -- the half-time hidden-parity Bellman sieve

THM-2226 leaves eight scalar profiles because its uniform three-checkpoint
Bellman carrier strips each coefficient to a `169`-adic time offset and then
uses the common cap `26/169`. The missing bit is parity: a coefficient with
one residual factor of `13` is observed halfway between two `169` checkpoints.
Retaining that intermediate joint state closes six of the eight profiles.

## 1. Half-time representation

Use the primitive normalized three-checkpoint triple from THM-2226. Write

```text
d_j=13^(v_j)u_j,
13 does not divide u_j,
0=v_1<=v_2<v_3.                                      (1)
```

The three checkpoint clauses are

```text
C_r=OR_(j=1)^3 1_(D_(169^r d_j)),        0<=r<=2.   (2)
```

Define the unit-core bit at every multiplication-by-13 half-time:

```text
X_j(n,x)=1_(D_(u_j))(13^n x).                        (3)
```

Then clause `r` is exactly

```text
C_r=OR_j X_j(2r+v_j,x).                              (4)
```

Thus relative valuation parity is a literal time coordinate, not merely a
choice between the static caps `25/169` and `26/169`.

## 2. Exact one-step hidden state

Let `T(x)=13x`. For every `13`-unit `u`, root counting gives, for almost
every terminal point `y`,

```text
sum_(Tx=y)1_(D_u)(x)=2-1_(D_u)(y).                   (5)
```

Consequently, if the future bit is `xi_j=X_j(n+1)`, the current bit
`eta_j=X_j(n)` has exact conditional marginal

```text
P(eta_j=1 | xi_j)=(2-xi_j)/13
 =2/13,  if xi_j=0,
 =1/13,  if xi_j=1.                                  (6)
```

For a fixed terminal point, the three current bits share the same root digit
and need not be independent. The safe relaxation therefore permits every
joint law of

```text
eta=(eta_1,eta_2,eta_3) in {0,1}^3                  (7)
```

having the three marginals (6). The actual root-induced law is feasible.
The endpoint-exceptional terminal set and its finitely many images and
preimages are null.

## 3. Clause-and-bit Bellman recurrence

Let `Omega={0,1,2}` be the clause set. At half-time `n`, define the active
incidences

```text
A_n={(j,r): n=2r+v_j, 0<=r<=2}.                      (8)
```

For a current bit vector `eta`, let

```text
Z_n(eta)={r: some (j,r) in A_n has eta_j=1}.         (9)
```

The Bellman state retains both a future-satisfied clause set `S` and the
future joint bit vector `xi`. Start with

```text
V_(-1)(S,xi)=1_(S=Omega).                            (10)
```

Recursively define

```text
V_n(S,xi)
 =max sum_eta p_eta V_(n-1)(S union Z_n(eta),eta),  (11)
```

where

```text
p_eta>=0,
sum_eta p_eta=1,
sum_(eta:eta_j=1)p_eta=(2-xi_j)/13,   j=1,2,3.      (12)
```

This has a uniform conditional interpretation. Fix an almost-everywhere
terminal point at time `n+1`, hence its bit vector `xi`. Across its thirteen
roots, the actual current-vector law satisfies (12). Given a current vector
`eta`, the induction hypothesis bounds all earlier roots by the corresponding
term in (11). Averaging and relaxing over every feasible coupling proves the
bound pointwise in the terminal cell.

This half-time composition is strictly stronger than coupling an entire
`169`-root pair in one step: two successive applications of (11) force the
existence of one common intermediate joint bit law. The one-step pair
relaxation forgets that compatibility.

Let

```text
H=max_j(v_j+4).                                      (13)
```

After half-time `H` no clause atom remains. Under Haar measure, the terminal
joint bit vector has arbitrary dependence but exact marginals

```text
P(xi_j=1)=1/7.                                       (14)
```

Hence

```text
measure K_3(d)<=B(v),                                (15)

B(v)
 =max sum_xi pi_xi V_H(empty,xi),                   (16)
```

where `pi` ranges over all joint laws with the three marginals (14).

## 4. Exact finite LP

Every transition LP has eight variables `p_eta` and the four independent
equalities in (12): total mass and three bit marginals. Every primal optimum
therefore has a basic feasible representative obtained from at most

```text
binom(8,4)=70                                        (17)
```

bases. The dual has four free variables and eight subset constraints, so its
vertices are obtained from the same 70 active-constraint choices.

Run

```bash
python3 04-computation/lrc14_half_time_hidden_parity_bellman_thm2227.py
python3 -O 04-computation/lrc14_half_time_hidden_parity_bellman_thm2227.py
```

The companion uses exact rational arithmetic and independently enumerates
both sides at every Bellman state, requiring exact primal-dual equality. It
also computes the exact same-core power-tower control from the reversible
single-core transition matrix

```text
[[11,2],
 [12,1]]/13.                                         (18)
```

Ordinary and optimized outputs are byte-identical. An independent
floating-point implementation gives the same decimals.

The four relative-valuation rows are:

```text
v          hidden Bellman bound

(0,2,4)    1086371907/5710115047 =0.1902539437...  FAIL
(0,2,5)     541866697/5710115047 =0.0948959333...  PASS
(0,3,4)     464433205/5710115047 =0.0813351747...  PASS
(0,3,5)   1005912252/10604499373 =0.0948571183...  PASS.   (19)
```

Here PASS means strict comparison with

```text
delta_5=961/6930=0.1386724386... .                   (20)
```

The exact positive margins in the three passing rows are respectively

```text
247469192851/5653013896530,
324128349931/5653013896530,
3219951991093/73489180654890.                        (21)
```

For coincident unit cores, the exact power-tower checkpoint measures are

```text
(0,2,4): 916159/4826809,
(0,2,5): 895649983/10604499373,
(0,3,4): 57121111/815730721,
(0,3,5): 895649983/10604499373.                      (22)
```

The first exceeds the target and is THM-2226's genuine geometric
obstruction; the other three lie below it. No equality case for the robust
Bellman relaxation is claimed realizable.

## 5. Profile consequence

THM-2226 leaves exactly

```text
(4,6,8), (4,6,9), (4,7,8), (4,7,9),
(5,7,9), (5,7,10), (5,8,9), (5,8,10).               (23)
```

Subtracting the first coordinate sends these to the four rows in (19), each
appearing twice. Therefore this theorem closes

```text
(4,6,9), (4,7,8), (4,7,9),
(5,7,10), (5,8,9), (5,8,10),                        (24)
```

and leaves only

```text
(4,6,8), (5,7,9).                                    (25)
```

Before adding the stronger concurrent carrier THM-2229, reconciling the
complete depth-four closure THM-2219 with THM-2226 gives `458` profiles. The
six exclusions above would reduce that ledger to

```text
458-6=452.                                            (26)
```

THM-2229 independently subsumes all six exclusions in (24) and closes
additional low-first-depth profiles. Thus (26) records this theorem's
standalone consequence relative to the THM-2219/2226 base; it must not be
subtracted again from the post-THM-2229 current ledger. The value of this
alternate carrier is structural: it exhibits the common intermediate joint
bit law that the coarser `169`-step coupling loses.

The two geometric-parity profiles (25), the other low-first-depth rows,
owner/current data, and LRC(14) remain open. QED.
