---
id: THM-827
title: Every active n=8 black curvature cell is outward under source conditioning and reverse under target conditioning
status: FINITE-EXACT (atlas-free n=8 joint endpoint-q disintegration) + PROVED MARGINAL CONSEQUENCES
source: codex-2026-07-15-S13-postjoin
depends_on: [THM-643, THM-811, THM-814]
related: [THM-785, THM-801, THM-825, HYP-6880]
verification:
  - 04-computation/black_curvature_joint_q_flow_codex_S13_postjoin.cpp
  - 05-knowledge/results/black_curvature_joint_q_flow_codex_S13_postjoin.out
---

# THM-827 — black curvature flow reverses with the conditioned endpoint

Orient every non-tied size-eight mixed--pure-black complement line toward
increasing cyclic-triangle count `C3`.  Let

```text
O_ij = # {mixed source of q=i -> pure target of q=j},
R_ij = # {pure source of q=i  -> mixed target of q=j}.    (1)
```

Let `M_i` and `P_i` be the complete masses of mixed and pure-black endpoints
of curvature `i`.  The two matrices in (1) have the same 17-cell support.
Every active cell satisfies the strict odds sandwich

```text
M_i/P_i < O_ij/R_ij < P_j/M_j.                           (2)
```

Equivalently,

```text
O_ij/M_i > R_ij/P_i        (source-conditioned outward),
O_ij/P_j < R_ij/M_j        (target-conditioned reverse). (3)
```

Thus THM-814's source-normalized outward dominance is not merely a marginal
effect: it holds separately in all 17 joint curvature cells.  But the same
17 cells all reverse when normalized by the target fibre.  Curvature
disintegration therefore has a genuine endpoint-direction asymmetry even
though the literal line relation is undirected before `C3` orientation.

## 1. Exact joint matrices

With rows indexed by source `q=0,...,6` and columns by target `q=0,...,6`, the
atlas-free census gives

```text
O = [1160 2826 1310    0    0 0 0]
    [   0 1600 2874  692    0 0 0]
    [   0    0 1002  324    0 0 0]
    [   0  704  940  176    0 0 0]
    [   0    0 1382  756   20 0 0]
    [   0    0    0  168   12 0 0]
    [   0    0    0    0  172 0 0],                     (4)

R = [2044 4442 2134    0    0 0 0]
    [   0 2784 5226 1076    0 0 0]
    [   0    0  932  424    0 0 0]
    [   0  784 1448  264    0 0 0]
    [   0    0 1542  848   32 0 0]
    [   0    0    0  316   52 0 0]
    [   0    0    0    0   22 0 0].                     (5)
```

All entries are even because staircase reflection acts freely on black lines
and preserves both endpoint curvatures, endpoint category, and the `C3`
orientation.  The common support lies in `|i-j|<=2`.  This band is structural:
for the apex-zero orientation THM-811/814 gives

```text
q(kappa x)-q(x)=b_left+b_right in {0,1,2},               (6)
```

and reversing the `C3` orientation only changes the sign of the endpoint-q
difference.

The raw totals are

```text
sum O_ij=16,118,          sum R_ij=24,370.               (7)
```

So the unnormalized aggregate is reverse-heavy, exactly as THM-814 found.

## 2. Both endpoint normalizations

The endpoint denominators are

```text
i             0       1       2       3      4      5    6
M_i        8,282  16,240  19,160   7,844  2,734    180  176
P_i      404,710 702,608 559,720 274,780 81,682 14,156  784. (8)
```

Dividing each row of (4) by `M_i` and each row of (5) by `P_i` verifies the
left inequality in (3) in every one of their 17 common nonzero cells.
Dividing each column of (4) by `P_j` and each column of (5) by `M_j` verifies
the right inequality in (3), again in every cell.  All comparisons are made
by integer cross multiplication in the verifier.

The source marginals reproduce THM-814:

```text
source q      0     1     2     3     4     5     6
sum_j O_ij  5296  5166  1326  1820  2158   180   172
sum_j R_ij  8620  9086  1356  2496  2422   368    22,  (9)
```

and `(sum O)/M_i > (sum R)/P_i` for all seven rows.  This now follows directly
by summing the cellwise left inequalities, since the denominators are fixed
within a row.

The target marginals are

```text
target q      0     1      2     3    4  5  6
sum_i O_ij  1160  5130   7508  2116  204  0  0
sum_i R_ij  2044  8010  11282  2928  106  0  0.         (10)
```

Here `(sum O)/P_j < (sum R)/M_j` for every active target stratum
`j=0,...,4`; `j=5,6` are ties only because no oriented cross-category edge
enters them.  Again this follows by summing the cellwise right inequalities.

Notice that target `q=4` has raw `204>106` but is still target-normalized
reverse because `P_4/M_4` is large.  Neither the raw sign nor source
conditioning predicts a target incidence rate.

## 3. Exact method

The verifier streams all `2^20` complement lines and uses no size-eight node
atlas.  For each of the `2^21` endpoints it constructs the marked-path
tournament, first checks score-complement symmetry, and then uses exact
degree-bucket anti-isomorphism backtracking to decide self-converseness.
THM-643 turns that predicate into the mixed/pure-black endpoint category.
The replay independently recovers

```text
score-symmetric candidates   744,678
self-converse tilings         58,712
black lines                1,046,528
tied cross-category lines     12,584.                    (11)
```

For a non-tied cross-category line, the sign of
`C3(kappa x)-C3(x)` chooses its source, after which both endpoint `q` values
are entered in (4) or (5).  Assertions recover every THM-814 source row,
both raw totals, common support, reflection parity, the 17/17 source-strict
and target-reverse counts, and the curvature band.  Runtime is under one
second with about 5 MB maximum RSS on the recorded run.

## 4. Structural reading

The data distinguish three measures on the same relation:

```text
edge counting measure,
edge incidence / source endpoint mass,
edge incidence / target endpoint mass.                  (12)
```

The first makes reverse edges more numerous.  The second makes every joint
cell outward.  The third makes every joint cell reverse.  The apparent
contradiction disappears because (2) places the same raw edge ratio between
two very different endpoint-category odds.  The mixed fibres are tiny
relative to the pure fibres at both ends, but their sizes vary strongly with
`q`; changing which endpoint supplies the denominator changes the question.

This rules out a tempting monotonic interpretation of black-edge curvature:
there is no endpoint-free statement that “curvature flow points outward.”
The correct object is the marked transport cell

```text
(source category, source q, target category, target q, C3 direction), (13)
```

together with the endpoint measure used for normalization.

## 5. Tournament Analysis and preservation boundary

Tournament Analysis uses the seven curvature strata as diagnostic vertices.
Its pairwise observable is the larger conditional outward-minus-reverse
marginal.  Switching from source to target conditioning flips one of the 21
pair orientations; both gauges are transitive, with score histogram
`{0:1,...,6:1}`, no directed triangle, singleton SCCs, and one Hamiltonian
path.  Increasing `q` is the tie Hamiltonian path.  The small gauge flip does
not weaken the cellwise reversal: the tournament ranks the *sizes* of seven
already-signed signals.

The challenged assumption is that vertices should be unmarked tournament
classes or runners.  For this question the proof-bearing vertices are
oriented complement-line incidences, while q-strata are only their
disintegration bins.  The calculation preserves both literal endpoints,
self-converse category, endpoint curvatures, `C3` direction, multiplicity,
reflection pairing, and the choice of endpoint denominator.  It destroys
node identity, positional layer data, owner labels, metric wall position,
continued-fraction phase, the LRC loneliness predicate, and all future
continuation behaviour. ∎
