---
id: THM-1261
title: The bounded-variation needle certificate closes carrier 44 exactly
status: PROVED FINITE-EXACT for carrier c=44, every slow-gap phase, and every faster integer speed.  Twenty-two reflection-representative 200-bin probability densities have exact closed-comb load below 1/6 on a finite vertical ladder and an exact BV-tail cap below 1/6 thereafter.  Consequently no six faster combs cover a complete 44-gap, even with repeated speeds.  Together with THM-1182/1197/1255/1257/1259 this closes every carrier 1<=c<=44; carriers c>=45 remain open
source: codex-2026-07-19 carrier-44 exact census
depends_on: [THM-1182, THM-1197, THM-1255, THM-1257, THM-1259]
related: [THM-1176, THM-1198, THM-1232, THM-1250, THM-1254]
script: 04-computation/lrc14_slow_gap_bv_needle_carrier44_thm1261.py
independent_verifier: 04-computation/lrc14_slow_gap_bv_needle_carrier44_thm1261_verify.py
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCCarrier44BVNeedle.lean
certificate: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier44_thm1261_certificate.json
output: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier44_thm1261.out
independent_output: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier44_thm1261_verify.out
---

# THM-1261 -- carrier 44 has no six-comb slow-gap cover

Let

```text
Dbar_d={t in R: ||dt||<=1/14},
G_k=[(14k+1)/(14*44),(14k+13)/(14*44)].              (1)
```

> **Theorem.**  For every integer `k` there is an absolutely continuous
> probability measure `mu_k`, supported on `G_k`, such that
>
> ```text
> mu_k(Dbar_d)<1/6                                    (2)
> ```
>
> for every integer `d>44`.  Hence no six faster danger combs cover a
> complete `44`-gap.

The statement permits repeated speeds and therefore excludes all primitive
or distinct six-comb packets as special cases.

## 1. The 22 exact density rows

Each representative gap is divided into `200` equal bins and assigned a
stored nonnegative integer probability density.  For a row-specific
`B_k=44M_k`, exact arithmetic verifies

```text
max_(45<=d<=B_k) mu_k(Dbar_d)<1/6,                    (3)

1/7+3 TV_R(mu_k)/(49(B_k+1))<1/6.                    (4)
```

The cutoff census is

```text
M          12  20  28  44  60  76
rows        7   1   8   4   1   1.                  (5)
```

All phases close inside the predeclared bank
`{12,16,20,28,44,60,76,108}`; neither `16` nor `108` is selected.  The
largest finite and tail loads are respectively

```text
3624398201278/21749999997651
 =1/6-401198887/14499999998434,                       (6)

51367660709987/308249999966709
 =1/6-14678568929/616499999933418.                    (7)
```

The floating LP only discovers integer weights.  Ordinary and optimized
replays check all stored loads with `Fraction` arithmetic and are
byte-identical.

THM-1182's BV estimate

```text
|mu_k(Dbar_d)-1/7|<=3 TV_R(mu_k)/(49d)                (8)
```

combines with (4) to control every `d>B_k`; equation (3) controls the finite
range.  An independent standard-library verifier recomputes the interval
loads using

```text
F(y)=floor(y)/7+min(frac(y),1/7),
|[u,v] intersect Dbar_d|
 =[F(dv+1/14)-F(du+1/14)]/d.                         (9)
```

It reproduces (6)--(7) exactly.  The finite maximum is at
`(k,M,d)=(20,28,58)`, with `d=44*1+14`; the tail maximum is the same density
row, whose normalized zero-extension variation is

```text
2332892045600/249999999973.                          (10)
```

## 2. Reflection has no fixed phase

After reducing `k` modulo `44`, reflection `t -> 1-t` sends

```text
G_k -> G_(43-k),                                     (11)
```

with endpoint identities

```text
1-(14k+13)/(14*44)=(14(43-k)+1)/(14*44),
1-(14k+1 )/(14*44)=(14(43-k)+13)/(14*44).            (12)
```

Every integer danger comb is reflection invariant.  The equation
`k=43-k` has no integral solution, so the `44` phases split into `22`
two-element orbits, represented exactly by

```text
k=0,1,...,21.                                        (13)
```

The Lean consumer proves (12)--(13), including representative coverage of
every `0<=k<=43`.

If six combs covered `G_k`, the exact measure would give

```text
1<=mu_k(union_(i=1)^6 Dbar_(d_i))
 <=sum_(i=1)^6 mu_k(Dbar_(d_i))<1,                   (14)
```

a contradiction.

## 3. What the depth jump means

The first substantial cutoff-depth jump in this carrier sequence occurs at

```text
(k,M)=(6,76), (7,60).                                (15)
```

These rows are not the margin bottleneck: their finite gaps below `1/6` are
about `6.33e-4` and `3.86e-4`, whereas the global finite gap at `k=20` is only
about `2.77e-5`.  The large cutoffs therefore finance the variation estimate
rather than reveal a remote dangerous speed.  This distinction matters:
extending a finite multiplier bank is not itself structural progress, while a
rule choosing densities with controlled variation would be.

Write `d=44q+s` and parameterize a gap by `t=(k+u)/44`.  Then

```text
dt=(q+s/44)u+sk/44                  modulo 1.         (16)
```

Thus a row is a density-selection problem on one normalized interval, driven
by the coupled slope `q+s/44` and phase shift `sk/44`.  The new mainline
certificate-rung analysis warns that bounded certificates repeatedly strand
on arithmetically structured blockers.  In this BV setting, (15) is the
analogous warning sign: the next high-leverage target is a parametric
selection/transport lemma for (16), not an inference that a longer scan has
made the unbounded carrier branch small.

## 4. Tournament loss audit

We challenged runners, comb teeth, walls, residues, density bins, gap phases,
wall-crossing events, and proof obligations as vertices.  The faithful proof
objects are density-decorated phase-obligation stalks: they retain both the
finite loads and the variation budget.  Orienting the 22 stalks by their
exact bottleneck load, using `k` only as a deterministic tie gauge, produces
a transitive tournament with score histogram `0,...,21`, no directed
3-cycles, singleton SCCs, one Hamiltonian path, and ordering

```text
14,0,21,15,13,11,17,6,7,2,4,8,10,12,9,18,3,5,19,16,1,20.  (17)
```

This quotient records relative difficulty but destroys the density weights
and therefore cannot carry the proof.  Tournament data here is a loss audit,
not a disguised runner tournament.

## 5. Formalization, hashes, and frontier

The independent replay passes all `22` rows.  The Lean module checks the two
rational margins, reflection endpoints, fixed-point-freeness, representative
coverage, the common load cap, and the six-load contradiction.  JSON
completeness, interval integration, and the analytic BV estimate remain
explicit external providers.

```text
generator_sha256   = c37c3dae32e0d13e664010f3f3fbafa2067a6f0a1925f51e2f68a7eb20ffe5aa
independent_sha256 = 341c976c7e83acf76d8bbf37e1a42b364003e13defd5deb173f2680ade2dbf38
certificate_sha256 = 05c71960933a7002ea51fb2719227f2a281c67efda3e9c6596e792ca9ce71644
output_sha256      = d95186d610169d260e13a05eeae3d354a99bcef9674eaf8d54721c587138b669
ind_output_sha256  = 9b224c83a411b65cfd02116c36bd03b9906c076360c7cb4a0f82f2c8fba991a0
lean_sha256        = 4bd64c4b1ef19fb6857e03c0dc112dfb5ebcce0617661dcf7622c172eb9a4724
certificate_bytes  = 53567
```

Together with THM-1182/1197/1255/1257/1259, this closes all carriers through
`44`.  It does not prove LRC(14) or uniform six-comb noncoverage.  Carrier
`45` is the next exact rung, but the depth profile in (15) says the more
valuable objective is the parametric phase-orbit density-selection theorem
that could absorb THM-1254's whole bounded `c<=1171` branch.
