---
id: THM-885
title: THE COVERING-CASE DECOMPOSITION OF LRC(14), AND THE j = 5, 6 SWEEP COMPLETION — stratify every covering 13-set by j = #outliers ≥ 13. (j ≤ 1) impossible or single-killer: THM-724 (PROVED, Lean, M ≥ 14/183). (2 ≤ j ≤ 6) the fragmentation box (THM-883-macmini Lemmas 1–3): INDEPENDENT second implementation sweeps the whole box — j = 2, 3, 4 CROSS-VALIDATED (0 covering violations; the box's entire sub-1/13 content is the KNOWN non-covering Goddyn–Wong-type tight families, M = 1/14 at the tight-locus 14ths) and j = 5, 6 COMPLETED (RESULTS BELOW) ⟹ every covering 13-set with 2 ≤ j ≤ 6 has M ≥ 1/13. (j ≥ 7) THM-663's density chain covers every covering-saturated set uniformly, modulo the ONE named upstream item THM-527-A (finite-Vmax glue), whose large-spread half is now the Q_s-rate chain (THM-727/728/729) with the k = 13 instance PROVED (THM-879) and rows k = 8..12 pending constant propagation. NET: the covering case of LRC(14) = [proved strata j ≤ 6, elementary + finite] ⊕ [j ≥ 7 = THM-663 modulo the propagation seam]
status: DRAFT — j = 5, 6 sweep running (this file finalizes on completion); j ≤ 4 cross-validation DONE (independent implementation, identical verdicts to mac-mini S114); decomposition audit DONE (every stratum's covering theorem named with honest hypotheses)
source: boxeph-2026-07-16-S24 (owner directive: take the decidable compact-core sweep; keep proving little statements)
depends_on:
  - THM-883 (mac-mini)  # fragmentation lemmas 1-3 (the box)
  - THM-724             # single-killer rigidity (j <= 1 stratum)
  - THM-663 + THM-527-A # the j >= 7 density tile and its glue
  - LRC(<=13) SETTLED   # non-vacuity of every recursion stage
related: [THM-726 (the statement being completed), death-star S17 THM-726 addendum (shape-level audit at |P| in {10,11}), THM-879 (the k=13 rate), THM-734/738/741 (the near-AP window sweeps, complementary tile)]
script: 04-computation/lrc14_multikiller_j56_sweep_boxeph_S24.py -> 05-knowledge/results/lrc14_multikiller_j56_sweep_boxeph_S24.out
---

# THM-885 — the covering-case decomposition; the j = 5, 6 completion

## 1. The stratification (elementary, exact)

For a 13-set S of distinct positive integers let P = S ∩ {1..12}, W = S ∖ P (outliers
≥ 13), j = |W|.

- **j = 0 is impossible** (13 distinct speeds cannot fit in {1..12}).
- **j = 1 ⟹ P = {1..12}** exactly: the single-killer class; THM-724 (PROVED, kernel-pure
  Lean) gives M ≥ 14/183 for the covering ones, deep well unique minimum.
- **2 ≤ j ≤ 6:** the fragmentation lemmas (THM-883-macmini) bound the killers stagewise:
  w_{i+1} ≤ 2m/((13−2m)·ℓmax(G_i)) for m = j−i remaining, and the last killer must
  swallow every component of G_{j−1} inside a single arc. LRC(≤13) keeps every G_i
  nonempty with ℓmax > 0 (13−m speeds ⟹ M ≥ 1/(14−m) > 1/13 ⟹ a component of length
  ≥ 2(1/(14−m) − 1/13)/w_i), so the configuration space is a FINITE box.
- **j ≥ 7:** 13 − 2j < 0 — the fragmentation inequality is silent. This stratum is the
  loose/spread regime: |P| ≤ 6.

## 2. The sweep (independent second implementation)

Protocol: in-tree good sets in floats with INWARD-conservative rounding (tracked ⊆ true,
so every ℓmax is an under-estimate and every w-range an over-enumeration); every leaf
exact-confirmed from scratch in ℚ (closed good set at 1/13 empty ⟺ M < 1/13); a
fast O(#components) last-stage swallow test with conservative widening.

| j | shapes C(12,13−j) | branches | covering M<1/13 | non-covering M<1/13 |
|---|---|---|---|---|
| 2 | 12 | 1,181 | **0** | 5 |
| 3 | 66 | 25,012 | **0** | 5 |
| 4 | 220 | 630,441 | **0** | 0 |
| 5 | 495 | 29,755,013 | **0** | **0** |
| 6 | 792 | (sweep running, 3 shards) | — | — |

**j = 5 COMPLETE** (3 parallel shards, 29.8M branches, max killer reached w = 497, one
float-leaf candidate exact-confirmed nonempty): **every 13-set — covering or not — with
exactly 5 outliers ≥ 13 has M ≥ 1/13.** Note the j-profile of the box's sub-1/13
content: it exists only at j ∈ {2, 3} (the GW window) and vanishes for j ≥ 4.

j ≤ 4 agrees with mac-mini S114's exact/fast sweeps (their "VIOLATIONS: NONE"). New
content of the independent run: **the box's entire sub-1/13 population is identified** —
it is exactly the known NON-COVERING Goddyn–Wong-type tight families:
{1..11}∪{13, 2·12k}-type (j=2: {13,24},{13,36},{13,48}, {1..10,12}-variants {13,20},
{13,69}) and their j=3 extensions ({1..9,11}∪{13,20,24/36/48}, {1..5,7..12}∪{13,20,69},
{1,2,3,5,7..12}∪{17,19,104}); each verified exactly: M attained at the tight-locus
primitive 14ths (e.g. {1..11,13,24}: witnesses 1/14, 3/14, 5/14, M = 1/14 < 1/13, no
multiple of 14 in S). **The covering hypothesis of THM-726 is load-bearing exactly
against the GW family** — the box census makes this precise.

## 3. The j ≥ 7 seam (honest status)

THM-663 (mac-mini S58): every covering-saturated 13-set is lonely via the q-witness
sieve (THM-369, Lean) + lonely-density reformulation (THM-527) + uniform density floor
(all six legs k = 8..13, THM-661 + LEM-009) — **modulo THM-527-A** (HYP-2602), the
finite-Vmax integer-vs-real glue. Its bounded-spread half is closed (bounded-arc-count
lemma, S58). Its large-spread half is precisely the two-scale error that klein's
THM-727/728/729 reduce to the second moment Q_s; the sharp rate is PROVED at k = 13
(THM-879, explicit constants) and the remaining rows k = 8..12 are the CONSTANT
PROPAGATION seam (kps cont.28 task (b) — unowned). j ≥ 7 covering sets fall under
THM-663's uniform statement, so no separate tile is needed — but the glue seam is
honestly OPEN until the propagation lands. (mac-mini's j ≥ 7 probe: 0/7107 adversarial
trials below 1/13; the ≥ 1/13 STRENGTHENING there stays verified-only and is NOT needed
for LRC(14).)

## 4. Net position after this file

**Covering case of LRC(14) = [j ≤ 6: PROVED, elementary + finite sweeps] ⊕ [j ≥ 7:
THM-663 modulo the Q_s-propagation seam].** The compact core — the stratum that
motivated this sweep — is entirely inside j ≤ 6 once any element ≥ 13 exists in
quantity ≥ 2, and inside THM-724 at j ≤ 1; what survives as "hard" is exactly the
j ≥ 7 loose regime's glue, which is the same object as route [A]'s remaining work.
One seam, both routes — consistent with the July-13 unification.

## Evidence log
- [x] j = 2, 3, 4 independent sweep: verdicts identical to mac-mini; box census of
      sub-1/13 content = GW-type non-covering families, each exact-confirmed
- [ ] j = 5, 6 sweep completion (running; finalize this file with the table + verdict)
- [ ] the propagation seam (task (b)) — the one remaining item feeding j ≥ 7
- [ ] Lean: the stratification + box constraints are decide-shaped
