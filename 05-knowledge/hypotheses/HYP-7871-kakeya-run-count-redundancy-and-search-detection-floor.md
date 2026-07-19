---
id: HYP-7871
title: The D>=3P Kakeya run-count route is redundant, while its witness law and branch-2 detection-floor audit remain informative
status: MIXED/RESCOPED. The intended 2/21 continuum ceiling was already proved more strongly by THM-1203. THM-1251 then confirms that the eroded start-complex supplier is moot because THM-1214 already closes its clustered r=5 stratum. The gap-coordinate ordering, non-self-mirror law, and lattice-tube interpretation are proved; the periodic witness law on (1,D-1,D) is computationally verified rather than formally proved. The adversarial-search audit proves that several reported search floors lie above 1/14 and hence cannot support negative evidence. No uniform branch-2 or n=12 conclusion follows.
source: kind-pasteur-2026-07-19-S128c84; namespace repaired by codex-2026-07-19-S82
depends_on: [THM-1203, THM-1214, THM-1245, THM-1251, MISTAKE-183]
related: [HYP-7870, HYP-7875]
---

# HYP-7871 — run-count redundancy and the search detection floor

This detail file repairs the index-only namespace collision recorded at the
top of `05-knowledge/hypotheses/INDEX.md`.  The first-pushed `HYP-7870` is the
seven-wall positioned-transfer hypothesis.  The later kind-pasteur entry is
therefore `HYP-7871`.

Its main proposed route was

```text
D>=3P   equivalently   number of runs=2P<=2D/3,
```

for a three-gap Kakeya geodesic.  That route is no longer a live closing
target: THM-1203 already proves the desired continuum bound
`mu(BAD_4)<=2/21`, with equality exactly at four-term arithmetic progressions,
by additive-triangle deletion.

Three pieces remain useful:

1. On the family `(1,D-1,D)`, computation through `D=700` finds

   ```text
   P(14q+r)=3q+e(r),
   e(r)=0 on {0,1,2,4,5},
        =1 on {3,6,8,9},
        =2 on {7,10,11,12,13}.
   ```

   This makes `D>=3P` immediate on that family, with equality only at
   `D=3`; the periodic formula itself still needs a proof if reused.

2. In ordered gap coordinates, every run stays in one permutation chamber,
   equal rates are excluded, and `u->1-u` has no admissible fixed run.  Hence
   the number of runs is exactly even.  Unrolling the geodesic identifies the
   run count with a lattice-point count in a convex tube.

3. A search for the `q0>14` branch must first rediscover a known `1/14`
   extremizer from inside its own domain.  The audited stochastic searches
   fail that gate: their best observed floors (`4/43`, `6/61`, `2/19`, and
   `1/10`) all exceed `1/14`.  Therefore a report of no counterexamples from
   those instruments is vacuous.  Individual grid witnesses remain rigorous
   positive certificates; `1,200/1,200` sampled structural rows had one, but
   this is not a uniform proof.

THM-1251 retires the eroded-start supplier: it lives entirely on the clustered
`r=5` stratum already closed uniformly by THM-1214.  The surviving targets are
only a genuinely separate certificate-producing structural theorem for the
whole `q0>14`/branch-2 regime, if that regime remains live after dependency
audit, and the separate all-height non-AP `n=12` problem.  This hypothesis
makes no LRC(14) closure claim.
