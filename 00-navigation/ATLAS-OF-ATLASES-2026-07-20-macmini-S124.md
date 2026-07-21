# THE ATLAS OF ATLASES — the navigation meta-index, and a measured dormancy census

> **DATED SNAPSHOT.** Its counts and “latest” labels were already stale by
> 2026-07-21. Use [`START-HERE.md`](START-HERE.md) as the live router and retain
> this file for its dormancy method and historical census only.

**mac-mini-2026-07-20-S124.** Owner directive: find the problems this repo mentioned
and then left behind; be comprehensive; build an atlas, or an atlas of atlases.

Two things were needed, and only one of them was "a list of ideas". The navigation
layer has grown to **41 files** with heavy supersession and no map of itself — that is
Layer 0. And "left behind" had never been *measured* — that is Layer 1. The headline
result is a negative that corrects the premise, so it goes first.

---

## HEADLINE — the repo does not have an abandoned-ideas problem

Measured, not asserted (method in Layer 3):

| layer | population | genuinely dormant | rate |
|---|---|---|---|
| tangents (`TANGENTS.md`) | 1,380 markers / 572 parsed entries | **50** raw → **~22** after Layer 2 | **3.6% → ~1.6%** |
| hypotheses (`hypotheses/INDEX.md`) | 1,946–2,707 parsed entries* | **0 open** (46 index-only, all closed) | **0%** |

Every one of the 46 index-only hypotheses carries a **closed** status — finished work,
not dropped work. And the 50 dormant tangents are not scattered: **~28 of them are one
cluster** — which, on inspection (Layer 2), turned out to be **live Lean code with stale
tangent labels**, not abandoned work at all. Follow-up discipline in this corpus is far
higher than its size would predict. **There is no graveyard.** What the census found is a
bookkeeping lag in one log, and a navigation layer that needed a map of itself.

---

## LAYER 0 — the atlas of atlases (41 navigation files, classified)

### 0.1 LIVE — canonical, updated 2026-07-19/20
| file | size | role |
|---|---|---|
| `SESSION-LOG` | 9.8 MB | the chronological spine |
| `TANGENTS` | 1.4 MB | 1,380 tangent markers, the idea log |
| `INVESTIGATION-BACKLOG` | 991 KB | the master lead list |
| `OPEN-QUESTIONS` | 835 KB | the live frontier |
| **`PROBLEM-LEDGER`** | 16 KB | **canonical** problem/results ledger |
| `LRC14-PROOF-MAP` | 102 KB | the proof skeleton |
| `LRC14-FORMALIZATION-MANIFEST-2026-07-17` | 134 KB | Lean surface |
| `CONSTANTS-INDEX`, `LRC14-FINITE-CHECK-FEASIBILITY-LEDGER-2026-07-19`, `LRC14-NEAR-MISS-LEDGER-…-klein-S319` | | live ledgers |

### 0.2 SUPERSEDED — the problem-ledger pile-up (one owner prompt, four drafts)
A fleet-wide prompt on 2026-07-20 produced four overlapping files. `PROBLEM-LEDGER`
(death-star-S59u) **explicitly consolidates the other three with credit**:
- `PROBLEM-ATLAS-2026-07-20` (kind-pasteur-S128c104) → merged
- `PROBLEM-PORTFOLIO-2026-07-20` (mac-mini-S140) → merged
- `PROBLEM-LEDGER-2026-07-20-klein-S332` → merged

> **Read `PROBLEM-LEDGER` only.** The other three are drafts of it.

### 0.3 SUPERSEDED — the frontier/status series (10 dated snapshots)
`LRC14-STATUS-07-08` → `-07-09` → `FINISH-MAP-07-11` → `-07-13` → `FRONTIER-07-14` →
`FRONTIER-AND-AVENUES-07-14` → `TRIANGULATION-07-14` → `FRONTIER-07-15` →
`FRONTIER-ASSESSMENT-07-16` → **`UNIFIED-FRONTIER-SYNTHESIS-2026-07-18`** (latest).

> **Read the 07-18 synthesis.** The nine earlier snapshots are history, not state.

### 0.4 STALE — not touched in ≥ 3 weeks (verify before trusting)
| file | last touched | risk |
|---|---|---|
| `CONCEPT-MAP` (430 KB) | 2026-06-28 | the big concept graph, ~700 theorems out of date |
| `HISTORIAN-TANGENT-INDEX` (90 KB) | 2026-06-29 | self-declares coverage "as of 2026-05-30" |
| `LRC-TECHNIQUE-MULTIVERSE-INDEX`, `LRC-CARRIER-PULLBACK-INDEX` | 2026-06-26 | technique maps |
| `CONCURRENT-SESSIONS` | 2026-06-12 | protocol doc |

### 0.5 DOMAIN ATLASES (topic-scoped, current enough)
`LRC-TECHNIQUE-INDEX` (799 KB), `LRC-TOURNAMENT-TECHNIQUE-INDEX` (737 KB),
`LRC-LENS-MAP`, `METAGRAPH-ATLAS`, `METAGRAPH-PRESERVATION-AVENUES`,
`N9-DEFECT-CONTINUATION-AVENUES`, `LRC14-CONTINUED-FRACTION-FRONTIER`.

---

## LAYER 1 — the dormancy census (the new measurement)

**Definition used.** An idea is *dormant* if its identifier appears **nowhere outside
the log that raised it** — raised once, never picked up. Fresh items (< 3 days) are
excluded: they are new, not abandoned.

**Tangents.** 50 dormant of 1,380. Distribution by era shows the corpus got *more*
disciplined over time, not less:

```
T 800– 899   1 / 60   ( 2%)
T 900– 999  48 / 98   (49%)   <-- the anomaly (Layer 2)
T1000–1099   4 /125   ( 3%)
T1100–1499   0 /249   ( 0%)
T1500–1599   0 / 40   ( 0%)   (5 excluded as fresh)
```

**Hypotheses.** 0 dormant-and-open. The hypothesis discipline mandated in `CLAUDE.md`
("log every hypothesis, confirmed OR refuted") is demonstrably being followed.

---

## LAYER 2 — the one genuine abandoned cluster

> **The codex Lean p0 / witness-floor formalization sprint, 2026-06-21/22.**
> ~28 of the 50 dormant tangents, in one thread: S86g (12), S81 (4), S79 (3),
> S71 (3), S77 (2), S78 (2), S86g2 (2). Tags: `#lean` 28, `#formalization` 20,
> `#p0` 12, `#witness-floor` 11, `#goodset` 7, `#bonferroni` 4.

The mathematical content that was built and then never cited again — a Lean stack for
the Part-A / positive-p0 route:

- `LRCWitnessPartA.lean` — finite-Part-A bridge; `#arcs(GOOD(E))` period-bounded signal
- `LRCGoodSet.lean` — the `phaseGapSet → goodSet` quotient (`fract(u−v)` identity)
- `LRCDenseCovers.lean` — `exists_phase_arc`, `coverSet_compl_subset_denseSet_compl`
- `LRCWitnessFloorConcrete.lean` — the concrete witness floor; `witness_pos_from_strict_cover_bound`
- `LRCWitnessBonferroni.lean` — the corrected positive-p0 route after the S27 retraction
- `LRCCoverBound.lean` (kps S31/S31b) — the elementary `hp0cap` bound
- `LRCGapReach.lean` — gap-reach / near-integer geometry
- `LRCP0Concrete` — the concrete-p0 surface

### RESOLVED — it is not abandoned at all

I flagged two readings and then ran the check rather than leaving it open. Every one of
the eight modules **still exists**, and **seven of the eight are still imported by
`TournamentH7/TournamentH7.lean`** (the build root):

| module | file present | imported in root |
|---|---|---|
| `LRCWitnessPartA`, `LRCGoodSet`, `LRCDenseCovers`, `LRCWitnessFloorConcrete`, `LRCWitnessBonferroni`, `LRCCoverBound`, `LRCGapReach` | yes | **yes** |
| `LRCP0Concrete` | yes | not in root (the single possible exception) |

**So the code is live and the tangent entries are stale *labels* on it — not orphaned
scaffolding for a superseded plan.** The p0/witness-floor stack was absorbed into the
current Lean surface; only the tangent log failed to record that it had been.

**Consequence for the headline.** With this cluster reclassified, the repo's genuinely
abandoned *mathematical content* is close to nil. What the census actually detected is a
**bookkeeping lag in the tangent log**, not lost work. The one concrete follow-up is to
check `LRCP0Concrete` — the only module not reached from the root — and either wire it in
or retire it.

---

## LAYER 3 — honest limits of this census

- **Blind spot.** The method detects ideas whose *identifier* is never re-cited. An idea
  cited repeatedly but never *advanced* looks healthy to it. So this measures dropped
  **threads**, not stalled **problems** — for the latter, `PROBLEM-LEDGER` and
  `OPEN-QUESTIONS` are the right instruments, and they are current.
- **Parsing.** *The hypothesis entry count varies with the regex (1,946 vs 2,707 on two passes); the decisive figure — index-only-and-OPEN = **0** — is stable under both. 572 of 1,380 tangent markers matched the `**T####** [src]` entry form;
  the rest are cross-references. Dormancy rates are over the 572 parsed entries.
- **Fresh-exclusion.** Items after 2026-07-17 are excluded as new. Five of today's
  tangents (T1551–T1555) would otherwise have been miscounted as abandoned — a false
  positive worth naming, since it is the obvious way to over-report a graveyard.
- **Not covered.** `04-computation` (10,672 files) was scanned for citations but not
  audited for orphaned scripts; a script that no document references is a different and
  probably larger category.

---

*Cross-links: `PROBLEM-LEDGER` (canonical problem/results ledger — read that, not its
three drafts); `LRC14-UNIFIED-FRONTIER-SYNTHESIS-2026-07-18` (current frontier — not the
nine earlier snapshots); `TANGENTS.md`, `INVESTIGATION-BACKLOG.md`, `OPEN-QUESTIONS.md`
(live). Method: `04-computation/atlas_dormancy_census_macmini_S124.py`.*
