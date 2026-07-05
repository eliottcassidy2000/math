# THE LRC(14) COMPUTATIONAL VALIDITY AUDIT — what survives scrutiny, what needed fixing

**mac-mini-2026-07-05-S56 (HYP-4122).** Owner directive: "really carefully consider
the recent and current work of agents and its validity."  Companion to the proof
map (lrc14-proof-map-macmini-S54.md, updated).  Fleet context: klein-S138 ran the
Lean-side hygiene concurrently (MISTAKE-106, grid-attainment dedup); this is the
computational-side audit.

## Finding 1 (MISTAKE-107, mine, FIXED at l=7): the l>=7 bounded sweep was a slice

My S55 sweep of the descent's bounded residual assumed each modulus r in
C ∩ {7..12} is SELF-carried (k_r = r).  Wrong scope: any lifted coordinate can
carry r (value 14 = 2·7 carries mod-7; value 70 carries mod-10 — re-admitting
C ∩ {10,11,12} ≠ ∅).  The swept 7.07M sets were the self-carried slice only.
FIX (this session): the corrected l=7 stratum — no height assumptions, all
patterns — is 4,570,156,404 leaves: 3.89B sieve-killed, 679,381,289
witness-cleared, **HARD = 0**.  650× the flawed sweep, same verdict.  l=8
running; l >= 9 is genuinely non-enumerable (~10^11+, sieve density 0.2–0.3)
and belongs to the descent/cluster leg — the values live in a 9.5-ratio window,
i.e. exactly opus-S81's cluster territory.
RULE: a forced-structure lemma about ONE coordinate does not bound a stratum;
check the quantifier's home before enumerating.

## Finding 2 (MISTAKE-108, mine, impact PROVEN ZERO): a filter stronger than its Lean source

Census filter F4 transcribed gap_forces_big_pair as w_max + w_2nd >= 38; the
kernel statement allows i = j (own-tooth-peak binding, denominator 2v), so the
faithful filter is w_max >= 19 — weaker.  The logs showed the symptom: exactly
ONE family F4-excluded at B=48.  The full excluded stratum (w_2nd <= 12 shapes)
was swept under the TRUE profile: **zero families pass** (pinning at q = 19/23
kills them) — the census verdicts stand unchanged.  C source fixed.
RULE: transcribe QUANTIFIERS, not intent; audit count-deltas between adjacent
filters (F3 − F4 = 1 was visible in the logs all along).

## Verified EXACT (the positive audit)

- **Independent re-implementation** (fresh Python, ∀-dilation pinning form, no
  shared code): the w_max = 25 census slice reproduces the C count to the
  family — 1,351 = 1,351, all witness-cleared.  The census pipeline is
  cross-validated end to end.
- **Filter-vs-Lean diff**: not_loose_dvd, not_loose_near_unit (∀a form),
  gap_compressed_24 (∃j ≠ i₀ explicit), peel_height_bound — transcriptions
  EXACT.  Only F4 diverged (Finding 2).
- **Fleet constants spot-checked**: kps-S8 pair rung 22B = (12/25)·(275/6)·B ✓;
  opus-S81 descent (gap (1−2ρ)/w, ratio 26/11, bottom 2/L ≈ 134) ✓ shapes
  consistent; klein tower_step_12 constants = my S53 chain caps ✓ (S53-verified).
- **l=3 strata**: already double-implementation verified (S53, exact level-count
  match).  l=1: kernel-checked (opus-S77).  l=2: single implementation (S52) —
  the census + hypercube partially overlap it; residual single-implementation
  risk noted (low: zero-violation verdicts over disjoint methods agree).

## Scope limits stated (what the evidence does and does not cover)

- Q50: evidence = exhaustive to height 48 + 1.97M to 52 + 12,095 CRT-lifts to
  1.16e14 + 679M l=7-stratum clearings.  It is a CONJECTURE; now formally
  pinned as TemplateDichotomy (kps-S10) — the open predicate of the NEW surface.
- Periodicity/pole: single-ray statements; multi-scale families need the
  descent.  The pole construction's families are loose (profile-passing ≠
  violating) — that is the point, not a gap.
- l=4: chain-domain sweep STILL RUNNING (6 workers, ~8+ h); exact domain volume
  recomputed this session [result in the .out] — my S53 "~1.5e12 cells"
  estimate was low; the honest ETA and the l=5/6 "one C session" ledger claim
  are RETRACTED pending the measured volume: l=5/6 tails presumptively belong
  to the walk-rung/descent theory, not enumeration.

## The reroute (post-audit state of the finish list)

1. **TemplateDichotomy** (= Q50, now THE named open predicate of the pinned
   surface lrc14_of_template_and_corner): per-q template verification at
   q = 26..50 — finite, decidable, THE prize.  The 100%-composite-CRT law says
   the template family can even be restricted to divisors of p·p'.
2. **Cluster/descent leg**: high ratio-clusters (opus) + l >= 9 bounded (this
   audit's honest hand-off) + l=5/6 tails (retracted from enumeration).
3. **CornerLonely**: unchanged (bands + composition note).
4. Numbering protocol (kps-S10): per-machine residue classes mod 10 — ADOPTED
   here; mac-mini takes ≡ 2 (this session = 4122; next mac-mini = 4132).
