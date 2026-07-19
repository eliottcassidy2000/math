---
id: THM-1235
title: SLACK 1 CANNOT BE BEATEN — IT IS MINIMAL FOR NON-EXTREMAL FAMILIES — AND THE TRUE n=14 GAP EDGE IS 3/41, NOT 2/27 — slack := 14D − s is an integer, and slack = 0 ⟺ M = D/(14D) = 1/14 ⟺ EXTREMAL, so among non-extremal families the minimum is slack = 1, attained by {1,…,11,13,36} at 3/41 (D=3, s=41). Slack is NOT dilation-invariant — D→kD, s→ks sends slack→k·slack, so 2·{1,…,11,13,36} has slack 2 and D≥2 "extremals" are just dilates (2·{1,…,13}: D=2, s=28, slack 0) — hence slack is meaningful only for PRIMITIVE families. The slack-1 ladder M = D/(14D−1) = 1/13, 2/27, 3/41, 4/55, … → 1/14 has only D=1 and D=3 realised in testing; and a 1552-family scan over four natural shapes with x ≤ 400 found NOTHING strictly inside (1/14, 3/41). So the attained gap edge is 3/41, refining THM-1230's (1/14, 2/27)
status: the slack = 0 ⟺ extremal equivalence is exact algebra; the dilation scaling is exact; the 3/41 witness is verified (THM-1230). The emptiness of (1/14, 3/41) is a SEARCH result over four shapes and 1552 families, not a theorem — and my previous session's framing of 2/27 as the gap edge was the Farey neighbour, not the attained one, which this corrects
source: opus-2026-07-19-S396 (owner: work the D=3 stratum and whether slack 1 can be beaten)
depends_on: [THM-1230 (the 3/41 witness and the (1/14,2/27) framing this refines), THM-1205 (the D-coordinate), THM-1050 (dilation), boxeph-S123 (the stratification)]
scripts: 04-computation/slack_ladder_opus_S396.py -> 05-knowledge/results/slack_ladder_opus_S396.out
---

# THM-1235 — slack 1 is the floor for non-extremals

> **⛔ AMENDMENT (kind-pasteur-2026-07-19-S128c86, exact, gate-verified): the D=2 slack-1
> rung IS realised — "only D=1 and D=3 realised" is corrected.**
> M({1,…,12, 26}) = **2/27 EXACTLY**, witness t = 2/27, active pair (1, 26), s = 27,
> D = 2, slack 14·2 − 27 = 1, primitive (contains 1). The family lies INSIDE this
> theorem's own scanned shape {1..12, x} at x = 26. One-line mechanism: at q = 13m+1,
> a = m, the far element 26 = 2·13 sits at distance exactly 2/q (since 13m ≡ −1), the
> base {1..12} sits in-band (residues m..12m), and nothing beats the m = 2 rung — the
> direct transfer of THM-633's ladder law to N = 13. Script:
> `04-computation/lrc14_ladder_realization_crossN_kps_S128c86.py` (+.out), gates include
> reproducing 1/14, 1/13, 2/25, 3/37, 3/41, 3/23, 5/33, 3/35. **Scope of the correction:**
> the INTERVAL claim (nothing strictly inside (1/14, 3/41)) is UNAFFECTED — 2/27 > 3/41
> lies outside the interval, which is presumably how the survey missed it (a rung search
> restricted to in-interval families is vacuous at D = 2, since 2/27 cannot lie there).
> Slack-1 realization now reads: **D=1 ✓ {1..12,14}, D=2 ✓ {1..12,26}, D=3 ✓
> {1..11,13,36}, D ≥ 4 open (4/55 canonical).** Three consecutive realised rungs sharpen
> the isolated-vs-accumulation question this file poses. See MISTAKE-187 and HYP-7890;
> the realization mechanism for WHICH rungs are attained is the 07-06 binder-gate theory
> (HYP-4506/4516/4572), the same object in an older vocabulary.

> **PROGRESS, NOT SETTLEMENT (opus-S397), see THM-1240.** The open question below is still open. What is now proved: M = D/s lies inside (1/14, 3/41) iff 41D/3 < s < 14D, and D = 1, 2, 3 admit NO integer s -- so the interval FORCES D >= 4 (one step beyond boxeph's D >= 3 for the wider Farey interval). Also 1/14 and 3/41 are Farey NEIGHBOURS, so 4/55 is the mediant and the unique least-denominator fraction inside -- the canonical target. ~12,400 families now tested with zero hits, including a residue-band construction built specifically for 4/55. But the candidate (D,s) list is INFINITE, so settling needs a BOUND ON D, which is the named missing ingredient.

## Slack 1 cannot be beaten

slack := 14D − s is an **integer**, and

> slack ≥ 0 ⟺ LRC(14) holds;  **slack = 0 ⟺ M = D/(14D) = 1/14 ⟺ extremal.**

So the only way to "beat" slack 1 is to be extremal. Among non-extremal
families **slack = 1 is minimal**, and {1,…,11,13,36} attains it at
M = 3/41, D = 3, s = 41 (41 ≤ 42). The question as posed therefore has a
clean negative answer, and 3/41 is optimal in its class.

## Slack is not dilation-invariant

Under V → kV, D → kD and s → ks, so **slack → k·slack**:

| family | M | D | s | slack |
|---|---|---|---|---|
| {1,…,13} | 1/14 | 1 | 14 | **0** |
| 2·{1,…,13} | 1/14 | 2 | 28 | **0** |
| {1,…,11,13,36} | 3/41 | 3 | 41 | **1** |
| 2·{1,…,11,13,36} | 3/41 | 6 | 82 | **2** |

So slack is meaningful only on **primitive** families, and the apparent
"D ≥ 2 extremals" are simply dilates — 2·{1,…,13} sits at D=2, s=28, still
slack 0. This is dilation invariance of M (THM-1050) expressed in the
(D,s) coordinates.

## The slack-1 ladder, and the corrected gap edge

Slack 1 means **M = D/(14D − 1)**: 1/13, 2/27, 3/41, 4/55, 5/69, … → 1/14.
Testing which rungs are realised:

| D | M | realised? |
|---|---|---|
| 1 | 1/13 | **yes** — {1,…,12,14} |
| 2 | 2/27 | not found |
| **3** | **3/41** | **yes** — {1,…,11,13,36} |
| 4…8 | 4/55 … 8/111 | not found |

And a targeted scan of 1552 families over four natural shapes —
{1,…,11,13,x}, {1,…,12,x}, {1,…,10,12,13,x}, {1,…,11,14,x} with x ≤ 400 —
found **nothing strictly inside (1/14, 3/41)**.

**This corrects my own framing in THM-1230.** There I called (1/14, 2/27) the
stability gap and reported it populated by 3/41. That is true but reads the
wrong interval as the gap: 2/27 is the *Farey neighbour* of 1/14, whereas the
*attained* edge is **3/41**. The honest statement is that the gap above the
floor is (1/14, 3/41), with 3/41 attained and its interior empty as far as
tested.

## The sharp open question

The unrealised slack-1 rungs D/(14D−1) for D ≥ 4 all lie inside (1/14, 3/41)
and accumulate at 1/14. So:

> **Is (1/14, 3/41) genuinely empty, or are the D ≥ 4 slack-1 rungs realisable
> by families outside the shapes tested?**

If empty, 3/41 is the true second value of the LRC(14) spectrum and the floor
is isolated. If some D ≥ 4 rung is realisable, 1/14 is an accumulation point
from above and no gap exists at all. That is a sharp, finite-flavoured target
and the natural successor to this session.
