# Lane E findings: independent audit of the LRC(14) mixed-sector closure chain

Auditor: death-star delegate (lane E), 2026-07-30.
Scope: MSG-2937/2938/2939/2940, MISTAKE-329, THM-2928 mixed four-aligned/
three-drift closure, frozen artifacts in 04-computation and
05-knowledge/results.  Repo read-only; all scripts in this directory.

## 1. The two claimed missed rows: ARITHMETIC VERIFIED-EXACT

Script: `laneE_verify_rows.py` (this directory).  Fully independent
reconstruction: the danger set is built pointwise from the raw comb
definition (boolean array over Z/L), NOT via the repo's range-merge code;
a separate arc-merge reimplementation drives the global scan; the repo's
`support_size_bitset` is imported only as a third cross-check.

Row 1: `F=(1,4,5,7,9,11)`
- `L = lcm(14v) = 2^3*3^2*5*7^2*11 = 194040 = D`  (exact)
- `|S_D| = 55392`  (pointwise count; repo bitset agrees)
- Reflection `x -> L-1-x` maps the safe set to itself (verified pointwise);
  parity classes: even `27696`, odd `27696` — exactly `S/2` each.
- `d3 = 194040`: `ell = ceil(194040/7) = 27720` (exact division),
  `g = gcd(d3,2) = 2`, `H = D/lcm(d3,2) = 1`,
  `fibre_cap(D,d3,2) = H*ceil(ell/g) = 13860`.
- EXACT parity test: `27696 > 13860` — KILL.  (claim verified)
- COARSE test `|S_D| > 2*ambient_capacity`: `55392 > 55440` — FALSE, missed.
- Common-u does NOT fire (`55392 <= 97020 + 27720`), and the full
  phase-free Lorenz screen for shape `(2,2,194040)` over ALL 142 proper
  nontrivial divisor quotients `q|D` (exhaustive in `s`, stronger than the
  breakpoint set) finds NO violation.  So the row genuinely reaches the
  double-2 branch as a Lorenz survivor, and the miss is a real missed kill.

Row 2: `F=(1,5,7,8,9,11)`
- `L = 2^4*3^2*5*7^2*11 = 388080 = D`; `|S_D| = 109044` (all three methods).
- Parity classes `54522/54522`; reflection-invariant.
- `d3 = 388080`: `ell = 55440`, cap `= 27720`.
- EXACT: `54522 > 27720` — KILL.  COARSE: `109044 > 110880` — FALSE, missed.
- Common-u no (`109044 <= 194040 + 55440`); full Lorenz screen over all
  178 proper quotients: no violation.  Real missed kill.

Verdict: the ambient-half-cap vs exact-parity-fibre discrepancy is REAL,
and both claimed witnesses are exact in every digit (L, S, parity loads,
caps, both test outcomes).

## 2. "Misses exactly these 2 kills": PROVED at universe level

Global superset scan (same script): over all `3,003` bodies and `251,536`
body/divisor rows (arc arithmetic; row counts match the frozen referee's
`body_count=3003`, `divisor_rows=251536`), for every support-hard row
(`S/D <= 26/31`) and every legal shape `(2,2,d3)` with `lcm(2,2,d3)=D`
(candidates `d3 in {D} u {D/2 if odd}`), the coarse and exact double-2
tests DISAGREE on exactly 2 occurrences — precisely the two rows above,
both with `d3=D`, one occurrence each.  No disagreement in the other
direction (the old coarse TSV had no unsound kill; it was exactly 2
survivors too fat).  Since the actual referee stream is a subset of this
universe and both witnesses are proved to be in the stream (Section 1),
the coarse referee misses exactly these 2 kills:
`29,221 - 2 = 29,219` and `6,754 + 2 = 6,756`.  MSG-2937's counts are
consistent and now fully independently confirmed.

Soundness note on the cap itself: the enlarged denominator-d needle is an
`ell`-term unit-step progression of classes mod `d`; for even `d` a unit
step is odd, so parities alternate along the progression and any parity
class contains at most `ceil(ell/2)` of the classes, each of ambient
multiplicity `H = D/lcm(d,2)`.  So `H*ceil(ell/g)` is a correct (and
sharp) `q=2` fibre cap; the reflection argument makes each parity load of
`S_D` exactly `S/2` (structural: tooth centers `k*(L/v)` form a
negation-closed set and half-open teeth reflect onto themselves, and the
reflection descends to every even `D | L`).

## 3. Repo state: the repair has ALREADY LANDED; no stale 29,221 anywhere

- Commit `650533ec3` (2026-07-29 14:51:33 -0600, "lrc14: close mixed
  three-drift sector") — 34-38 minutes AFTER MSG-2937/2938 — introduced
  `lrc14_three_drift_mixed_lorenz_activity_thm2928.py` already containing
  the exact parity test (`third_parity_cap = height*ceil(ell/common)`,
  compared against `support_count//2`, with a parity-balance `require`).
  `git log --follow` shows no earlier buggy version was ever committed.
- Frozen `.out` asserts and records
  `activity_status_counts={'common-u': 383389, 'double-2': 6756,
  'survive': 29219}` — the corrected numbers.
- Grep of `01-canon` and `00-navigation` for `29,221/29221/29,364/29364`:
  the only `29,221` is in `01-canon/MISTAKES.md` (MISTAKE-329), correctly
  recorded as SUPERSEDED, with the same two witness rows and the corrected
  chain `544,571 -> 419,364 -> 29,219`.  THM-2928 canon states only the
  sound chain `544,571 -> 419,511 -> 29,364 -> 19 -> 0` and explicitly
  notes "even total cardinality is not used as a substitute for reflection
  balance".  `00-navigation/CURRENT-FRONTIER.md` cites THM-2928/2941
  ("literal zero/one/two/three-drift sectors are empty") without the raw
  counts.  Nothing stale to repair.

## 4. Counts and hash chain

Ledger reconciliation (frozen `.out` files, exact):
- Activity/fibre referee: `544,571` all-top; `125,060` single-fibre kills
  -> `419,511`; common-u `383,391`; `double_two_after_general = 6,756`,
  all killed by the exact parity test -> `final_occurrences = 29,364`
  (rows 2,974, shapes 4,993, bodies 2,878, divisors 110).
- Lorenz sidecar: `544,571 -> 419,364` (extra `147` higher-Lorenz kills),
  common-u `383,389`, double-2 `6,756`, survive `29,219`.
- Cross-consistency: `419,511 - 419,364 = 147`; common-u delta `2` plus
  survivor delta `145` sums to `147`.  Exact.
- Threshold transport: `inherited_activity_residual = 29,364`, kills
  `29,345`, leaves `19` on five bodies (2+4+2+1+10, matching THM-2928
  table (37th)); terminal local-incidence audit kills `19/19` at moduli
  `Z/99` (18) and `Z/98` (1), survivors `0`.
- IMPORTANT STRUCTURE: the closure chain inherits the ACTIVITY residual
  (29,364), not the Lorenz TSV; the Lorenz sidecar (29,219) is a
  non-load-bearing tightening.  This is exactly the architecture MSG-2938
  demanded.
- SHA-256 hash DAG: every `EXPECTED_*` constant embedded downstream
  matches the actual file hashes upstream —
  support script `778842c0...`/out `648327d3...` (in combined script);
  combined script `42dc1657...`/out `2e211620...` (in activity, Lorenz,
  threshold, local scripts and recorded in all four `.out`s);
  activity script `067424a0...`/out `cab8e7b4...` (in threshold script and
  `.out`); threshold script `13e524e7...`/out `435c34b2...` (in local
  script and `.out`); semantic hashes `88ac1b6c...` (all-top),
  `87fc8e09...` (q-survivor), `6a0b57e9...` (29,364 residual),
  `088b4e16...`/`ca19dda6...` (threshold), `3aa5d8b3...` (terminal-19),
  `107548c4...` (audit witness) all agree across producer and consumer.
  Verified with `sha256sum`; no dangling or mismatched link.
- Independent corroboration of a claimed S: `S=109044` for
  `(1,5,7,8,9,11)@388080` also appears verbatim in the threshold `.out`
  terminal witnesses (10 of the 19), computed by a different pipeline.
- Full recompute of `29,364`: rerun of the activity/fibre referee was
  attempted under a 10-minute cap (`activity_rerun.out` here); the script
  self-asserts every ledger (it hard-fails on any drift) — see rerun note
  at the bottom for whether it finished inside the cap.  The Lorenz,
  threshold, and terminal referees were not rerun (each rebuilds the
  251,536-row bitset universe; well over the cheap budget).  Their frozen
  outputs are consistent with each other, with canon, and with my
  independent 2-row + flip-set computations.

## 5. What is solid / what needs the root team / what death-star can do

SOLID (my independent verification):
- MISTAKE-329's two witnesses, every digit; the flip set is exactly 2.
- The corrected numbers 29,219/6,756 and the superseded status of 29,221.
- The canonical closure chain's hash DAG and count reconciliation.
- Canon and navigation contain no stale number.

DEFERS TO ROOT-TEAM AUDIT (not independently rerun here):
- Full replay of the four big referees (each rebuilds the 251,536-row
  universe; MSG-2938 says root is hostile-auditing and controls the
  canonical package — do not collide).
- Soundness review of the 8-state transport LP as an upper relaxation
  ((37tg)/(37tg+)) and the Z/99, Z/98 local set-cover mechanics; the
  frozen scripts carry positive/hostile controls, but I did not re-derive
  the LP duality certificates.

USEFUL NON-COLLIDING NEXT STEPS FOR DEATH-STAR:
- Nothing to repair in canon for this sector; do NOT rebuild the Lorenz
  TSV story.  THM-2958 is reserved by opus (MSG-2939) — send them this
  audit as supporting material rather than opening a competing stub.
- The transferable content for death-star's own lanes is the MISTAKE-329
  rule (symmetry-balanced load must be compared against the exact cap in
  the symmetry fibre, not an ambient cap divided by orbit size) — the
  same reflection/parity pattern appears in the GMC/LRC positivity
  manoeuvre residuals (arc-endpoint jump sums).
- MSG-2940 (rank>=35 next-prime signal false) is a separate lane (opus);
  no action from this audit.

## Artifacts in this directory
- `laneE_verify_rows.py` — independent two-row verification + global flip
  scan + repo cross-check (runtime ~3 min).  Output ends
  `LANE-E VERDICT: PASS`.
- `activity_rerun.out` — attempted 10-min-capped replay of the
  activity/fibre referee.

## Rerun note
The capped replay of `lrc14_three_drift_mixed_activity_fiber_thm2928.py`
COMPLETED in under 10 minutes and its output is BYTE-IDENTICAL to the
frozen `05-knowledge/results/lrc14_three_drift_mixed_activity_fiber_thm2928.out`
(39 lines, `final_occurrences=29364`, `double_two_parity_kills=6756`,
`final_semantic_sha256=6a0b57e9...`, `all_exact_controls=PASS`).  The
load-bearing 29,364 residual is therefore VERIFIED-EXACT by full replay,
not just by hash chain.
