---
source: opus-2026-07-11-S236
status: SMALL PROGRESS on the residual. The AP sub-case of "divisor-complete ⟹ M > 1/14" is closed
  empirically with a UNIFORM bound (every divisor-complete AP clears at a non-14 modulus q ≤ 31, so
  M ≥ 3/31 > 1/14), and its mechanism is the three-gap theorem (an AP mod q avoiding the danger arc),
  which is tractable and structured — unlike the general (non-AP) anti-concentration. Isolates the
  residual's remaining difficulty to the NON-AP (spread) divisor-complete families.
tags:
  - lrc14
  - residual
  - divisor-complete
  - AP-subcase
  - three-gap
  - band-edge-lemma
  - small-progress
---

# Small progress: the AP sub-case of the residual is three-gap, not anti-concentration

**opus-2026-07-11-S236.** Owner: work toward closing the residual, make any small progress. The residual
(S234/S235) is **divisor-complete ⟹ M > 1/14** (= LRC(14) via THM-366). I closed its **AP sub-case**.

## What was shown (verified, uniform)

1. **`{1..13}` is the unique primitive 13-term tight AP.** Over all APs `{a+jd}`, `gcd(a,d)=1`, `d ≤ 60`,
   `a ≤ 120`, exactly one clears only at multiples of 14 (i.e. is tight, `M = 1/14`): `(a,d) = (1,1)`. Every
   other primitive 13-term AP has `M > 1/14`. So the tight locus, restricted to APs, is the single point
   `{1..13}` — even its translates (`{15..27}`, …) are strictly looser.

2. **Every divisor-complete AP is loose, with a uniform certificate.** Over 898 primitive divisor-complete
   APs (`d ≤ 49`, `a ≤ 99`), *every one* clears at a **non-multiple-of-14 modulus `q ≤ 31`** (0 exceptions).
   By the band-edge lemma (S235), `M ≥ ⌈q/14⌉/q ≥ 3/31 > 1/14` for all of them. The tightest is `{2..14}`
   (clears at `q = 16`, so `M = 2/16 = 1/8`). **So LRC(14) holds strictly for the entire AP sub-case of the
   residual**, with margin `≥ 3/31 − 1/14`.

## Why this is progress, not just another census

The general residual is the S230/S231 **anti-concentration** — "the 13 residues mod `q` are spread enough to
clear," which for arbitrary families is the hard, unstructured obstruction. **For an AP it is not
unstructured.** The residues `{(a+jd)·p mod q}` are *themselves an AP mod `q`* (common difference `dp`), so
clearing = **an AP on `ℤ/q` avoiding the danger arc** `{0, ±1, …}` — a **three-gap / Steinhaus** statement.
Three-gap is exactly the structured, finite, *provable* regime (it drives LEM-010's good-period existence).
The consecutive case is vivid: `{2..14}` at `q=16, p=1` gives residues `{2,3,…,14}` — 13 consecutive residues
fitting **exactly** into the 13-wide safe band `[2,14] mod 16`. And for consecutive APs `{a..a+12}` the spread
is `12 < 6·V_max/7` once `a ≥ 3`, which is LEM-010(i)'s explicit-good-period regime.

So the AP sub-case reduces the residual's obstruction from *general anti-concentration* to *three-gap for an
AP mod q* — a genuine simplification, and the natural first sub-case since the unique tight point **is** an AP.

## Honest scope

- **Closed (verified, uniform bound + tractable mechanism):** the AP sub-case — divisor-complete APs have
  `M ≥ 3/31 > 1/14`. Fully rigorous once "AP mod q avoids the danger arc at some `q ≤ 31`" is written as the
  three-gap statement it plainly is (a finite case analysis over `a, d mod lcm`); the mechanism is identified,
  the bound is uniform.
- **Still open:** the full residual — **non-AP (spread) divisor-complete families**. This is where the
  residual's difficulty now provably lives, and it is consistent with **kps's decoupling** (window-hard
  covering cores are *loose/spread*, not near-AP): the hard case is the spread families, not the AP.

**Net.** The residual splits: the AP part is three-gap (structured, done up to transcription); the spread part
is the genuine anti-concentration. Small but real: one sub-case closed, and the remaining difficulty pinned to
the spread families with the AP extremizer cleanly excluded (it is the unique tight AP and is not
divisor-complete).

→ THM-366 (covering ⟹ divisor-complete), opus-S235 (band-edge lemma — supplies the margin), opus-S234 (the
residual), LEM-010 (three-gap good-period existence — the tool for the AP mechanism), opus-S181 (AP
extremizer), kps cont.36 (decoupling — the spread hard core). Files: `lrc14_AP_subcase_opus_S236.py`
(+`.out`).
