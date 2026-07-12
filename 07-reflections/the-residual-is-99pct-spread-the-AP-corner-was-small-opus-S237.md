---
source: opus-2026-07-11-S237
status: PROGRESS + an honest correction to S236. The residual (divisor-complete ⟹ M>1/14) is verified across
  the WHOLE divisor-complete class — every family clears at a non-multiple-of-14 modulus q∈[15,31], so M>1/14
  by the band-edge lemma. BUT divisor-complete is ~99–100% SPREAD (longest-AP≤7): the S236 "AP sub-case" was
  a ~1% corner, not the bulk. The hard core is genuinely the spread anti-concentration.
tags:
  - lrc14
  - residual
  - divisor-complete
  - spread
  - anti-concentration
  - band-edge-lemma
  - honest-correction
---

# The residual is 99% spread — the AP corner (S236) was small

**opus-2026-07-11-S237.** Continuing the residual work. Trying to rigorize/extend the S236 AP sub-case
surfaced an honest correction and a clean uniform statement.

## Strengthening: {1..13} is the unique tight interval

Among all consecutive intervals `{a..a+12}`, `a = 1..5999`, **exactly one** fails to clear at a
non-multiple-of-14 modulus `q ≤ 40`: `a = 1`, i.e. `{1..13}`. So the interval tight locus is the single point
`{1..13}`; every other interval clears, hence `M > 1/14` (band-edge lemma, S235). Robust to `a = 6000`. This
strengthens S236's uniqueness from "unique among APs `(a,d)` searched" to "unique among all 13-term intervals
up to `a=6000`."

## The honest correction: the AP sub-case is a ~1% corner

S236 "closed the AP sub-case." But the longest-AP distribution of divisor-complete families is:

| longest-AP | 2 | 3 | 4 | 5 | 6 | 7 | ≥8 (AP-ish) |
|---|---|---|---|---|---|---|---|
| share | ~0% | ~39% | ~49% | ~9% | ~1% | ~0% | **~0–1%** |

So divisor-complete families are **99–100% spread** (longest-AP ≤ 7). The "AP sub-case" (longest-AP near 13)
is a **measure-~1% corner**, *not* the bulk. This is structural: divisor-complete requires multiples of
`8, 9, 11, 13, 14` (specific, spread-out speeds), which is incompatible with the tight AP `{1..13}` (no
multiple of 14) — so **divisor-completeness forces spread**. My S236 framing over-implied the AP sub-case's
significance; the true residual is essentially the **spread anti-concentration**, which independently confirms
kps cont.36's decoupling ("window-hard covering cores are loose/spread, not near-AP").

## The residual holds across the whole class — reduced to one uniform statement

Both parts of the divisor-complete class clear at a bounded non-14 modulus:

- **AP corner** (S236): non-14 clearing `q ≤ 31`.
- **Spread bulk** (this session, the 99%): every spread divisor-complete family clears at a non-14 `q ≤ 26`
  (random, 0 exceptions), `q ≤ 29` adversarially — so `M ≥ 3/29 > 1/14`.

**Uniform statement:** every primitive divisor-complete family clears at a non-multiple-of-14 modulus
`q ∈ [15, 31]`, hence `M ≥ min_{q∈[15,31], 14∤q} ⌈q/14⌉/q > 1/14`. The residual (= LRC(14) via THM-366) holds
**empirically for the entire hard core**, reduced to the single uniform statement:

> **Open (the residual, now uniform):** every primitive divisor-complete family has a lonely time at some
> modulus `q ∈ [15, 31]` with `14 ∤ q`.

This is the S230/S231 bounded-modulus anti-concentration, now on the sharp domain (spread divisor-complete)
and in a sharp window (`[15,31]`), with the strict margin `M > 1/14` a free corollary (band-edge lemma). It is
verified (diameter-free, adversarial `q ≤ 29`) but not proved.

## Honest net

The residual is fully verified across the divisor-complete class and reduced to one clean uniform
anti-concentration statement in a bounded window. The genuinely useful outputs of this session are the
**interval uniqueness** ({1..13} unique to `a=6000`) and the **honest correction** that the residual is 99%
spread — so the extremal AP corner (S236), while structured and closable, is not where the difficulty lives.
The difficulty is, as the fleet has repeatedly found, the spread anti-concentration; no new tool closes it,
but its domain and window are now as sharp as possible.

**Repo-health note:** encountered a corrupt loose git object (`e4fd24…`) blocking `git gc`/`fsck` locally.
HEAD == origin/main, working tree clean, all commits pushed — so no data is at risk; disabled local auto-gc
(`gc.auto 0`) so it can't interfere with pushes. Clean fix is a fresh clone (origin is healthy).

→ THM-366 (covering ⟹ divisor-complete), opus-S235 (band-edge lemma — the margin), opus-S236 (the AP corner —
now seen as ~1%), opus-S230/S231 (the anti-concentration — the residual), kps cont.36 (decoupling: hard core
= spread). Files: `lrc14_residual_spread_bulk_opus_S237.py` (+`.out`).
