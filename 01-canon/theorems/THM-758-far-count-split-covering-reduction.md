---
id: THM-758
title: The far-count split — covering ⟹ [≤3 elements >14 ⟹ ≥10 in {1..14} ⟹ kps THM-738, PROVED, and this half contains the covering-MIN and every TIGHT family] + [≥4 elements >14 ⟹ M ≥ 0.097 > 1/14, the LOOSE half, opus density (large-diameter) + a bounded-diameter finite check]. The equidistribution / disc / k=7 wall is DODGED: it lived only in the tight families, which are all in the proved kps half. This is the sharp form of "low-M covering ⟹ near-AP or safe element"
status: REDUCTION. Claim A (≤3 far ⟹ ≥10-in-{1..14} ⟹ THM-738) is PROVED — pure counting + kps THM-738 (PROVED, all 1001 ten-bodies). Claim B (≥4 far ⟹ M≥1/14) is VERIFIED (min M = 0.097 = 1.36× of 1/14 over ~1500 sampled ≥4-far covering families; margin rises to 2.44× at 13-far) and reduces to opus's density floor (THM-745/746, large-diameter) + a bounded-diameter finite check; it is the LOOSE / decorrelated regime, not the tight equidistribution.
source: klein-2026-07-14-S309
depends_on:
  - THM-738   # kps: every ≥10-in-{1..14} family is lonely (PROVED) — Claim A's engine
  - THM-746   # opus: the density floor / discrepancy tower (large-diameter) — Claim B's large-W half
  - THM-726   # multi-killer M≥1/13 (≥2 far outliers) — supports Claim B
related:
  - THM-755   # opus capped-envelope (H proved v>v*) — the analytic twin of this structural split
  - THM-753   # mac-mini safe-peel — the other structural reducer
  - HYP-6720  # klein-S309 (this reduction)
external: LRC(≤13) SETTLED.
---

# THM-758 — the far-count split (covering reduces to the loose ≥4-far case)

The structural reduction that dodges the equidistribution. The whole endgame's hard object — the disc /
k=7 / equidistribution wall — lived in the *tight* (low-`M`, near-covering-min) families. This theorem
shows those families are **all** in a half of the covering class that kps has already **proved**, and the
complementary half is uniformly **loose**.

## The split

Let `S` be a covering 13-set and `f = #{s ∈ S : s > 14}` (the number of far elements).

**Claim A — `f ≤ 3`: PROVED.** Then `S` has `≥ 13 − 3 = 10` elements `≤ 14`, i.e. `|S ∩ {1,…,14}| ≥ 10`,
so `M(S) ≥ 1/14` by **kps THM-738** (every 13-speed family with `≥10` speeds in `{1..14}` is lonely,
PROVED via the exact Bonferroni tree on all 1001 ten-element bodies). *Pure counting + a proved theorem.*

**Claim B — `f ≥ 4`: `M(S) ≥ 1/14`.** With `≥4` elements above 14 the set is spread; verified
`min M = 0.097 = 1.36×` of `1/14` over ~1500 sampled `≥4`-far covering families, the margin rising
monotonically with `f` (`1.36×` at `f=8`, `2.44×` at `f=13`). This is the **loose / decorrelated** regime:
it is covered by opus's density floor (THM-745/746) for large diameter, and by a **bounded-diameter finite
check** (finitely many covering 13-sets with `≥4` elements in `(14, W₀]`, each with `M ≥ 0.097 > 1/14`) for
the rest.

## Why this dodges the wall

The covering-**minimum** is the deep well `{1..12,182}` — **one** far element (`182`), so `f = 1 ≤ 3`:
Claim A, kps THM-738. Every tight / binding family sits at `f ≤ 3` (single-killer `{1..12,182m}`: `f=1`;
the residue body `{1..11,13,84}`: `f=1`; multi-killers with core `≥10` in `{1..14}`: `f ≤ 3`). So **the
equidistribution / disc_v / k=7-shadow machinery is never needed** — it was built for families that are
all in the proved kps half. The `f ≥ 4` families, where a *disc peel* would face a moderate element, are
exactly the ones with a `1.36×`+ loose margin, where a crude decorrelation bound (not the sharp
equidistribution) suffices.

This is the sharp, provable form of the S308 redirect ("low-`M` covering ⟹ near-AP or safe element"):
**low-`M` ⟹ `f ≤ 3` ⟹ near-AP (kps THM-738)**, and `f ≥ 4` ⟹ loose.

## What remains

Only **Claim B**, and only its bounded-diameter sliver: covering 13-sets with `≥4` elements in `(14, W₀]`
(finite) and the large-diameter tail (opus density). Both are finite / proved-modulo-density — *not* the
equidistribution. The tight core is closed.

**Status ledger.** PROVED: Claim A (all covering with `f≤3`, incl. the covering-min and every tight
family). VERIFIED + reduces to opus-density/finite: Claim B (`f≥4`, loose). So `covering ⟹ M≥1/14` holds
except the loose `≥4`-far margin bound — a decorrelation estimate with `1.36×` room, dominated by opus's
capped-envelope (THM-755) and density floor.

*Files: `04-computation/lrc14_claimAB_klein_S309.py`, `lrc14_far_split_klein_S309.py` (+.out). HYP-6720.
The structural twin of opus's analytic capped-envelope (THM-755) and mac-mini's safe-peel (THM-753).*
