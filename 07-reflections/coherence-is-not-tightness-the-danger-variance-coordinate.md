# Coherence is not tightness: the danger-count variance, and why no pairwise invariant can characterize the AP

*klein-2026-07-18-S324. I set a filter in S323 — a transfer candidate must distinguish the AP from a
random covering set. This finds one that does, spectacularly, and then shows it still fails, for a reason
that generalizes to every pairwise-only invariant.*

## The candidate

`X(t) = #{i : ‖v_i t‖ < 1/14}` is the **danger count** — mac-mini-S89's object, playing the tournament
score/adjacency role in the shared partition function `Q(w) = E_t[w^X]` (`Q(2)` tournaments, `Q(0) = L`
loneliness). Its mean is `13·(2/14) = 13/7` for **every** 13-set, so the mean carries nothing and all
content sits in `Var_t(X)`. The tournament analog is score variance / eigenvalue flatness — what Paley
extremizes (THM-133).

`Var(X)` is a pure **pairwise** functional: expanding, `Var(X) = 13/7 + Σ_{c≠c'} μ(D_c ∩ D_{c'}) − (13/7)²`,
with the independent-events reference `78/49 = 1.5918`.

## It passes the filter, decisively

| family | `Var(X)` | vs independent |
|---|---|---|
| deep well `{1..12,182}` | 2.8342 | +78% |
| `2·{1..12} ∪ {13}` | 2.8328 | +78% |
| **AP `{1..13}`** | **2.8116** | **+77%** |
| GW `{1..11,13,24}` | 2.5977 | +63% |
| `{1..11,13,84}` | 2.5248 | +59% |
| **random primitive covering** | **1.55 – 1.92** | **≈ +5%** |

**The structure-vs-randomness answer:** random covering sets are *essentially independent* (1.66 against
the independent value 1.59, only +5%); the AP and every near-tight object are strongly positively
correlated (+60–78%). Coherence — when one runner is near an integer, many are — is exactly what the
near-tight families have and random ones lack. Every named low-`M` object sits at the **100th percentile**
of the random cloud.

## But it is a coordinate, not a threshold — and not a detector

Three checks, each refining the picture:

1. **No gap.** Interpolating AP → random by replacing `k` of the 13 speeds gives a *smooth* decay:
   `2.81, 2.65, 2.49, 2.32, 2.20, 2.08, 1.95, 1.87, …, 1.68`. My "coherence gap" hypothesis is **wrong**;
   the apparent gap in the census was a *sampling* gap — random 13-sets are essentially never near-AP, so
   they all pile up at the bottom of a continuous scale.
2. **Not predictive.** Within the random cloud, `corr(Var(X), M) = −0.065`; top-decile vs bottom-decile
   variance give mean `M` of 0.2119 vs 0.2201 (ratio 1.04). It separates the resonant cluster from the
   generic population but does not rank inside it.
3. **Not variational.** The AP does **not** maximize it: **585** single-replacement families beat
   `Var(AP) = 2.8116`, topping out at **3.0576**. The winners drop a speed coprime to 14 and add one of
   `14, 28, 42, …` — multiples of 14 sitting *exactly on the 1/14 grid*. They are **non-covering**, with
   `M ≈ 0.0769–0.0909`, comfortably above tightness. Pure resonance artefact.

## The structural conclusion — and it generalizes

`Var(X)` is a complete summary of the **pairwise** overlap data. It cannot characterize the AP. Therefore

> **no invariant built from pairwise overlaps alone can characterize tightness.**

The reason is visible in check 3: coherence and tightness come apart. A speed on the `1/14` grid maximizes
pairwise correlation while *helping* loneliness (it vacates the dangerous window), so pairwise statistics
reward exactly the wrong thing. Tightness needs the higher-order/covering structure — which speeds are
divisible by which `d` — that pairwise data discards.

This independently corroborates **opus-S358 / THM-1026** ("pairs insufficient at 13", the five-slot
ledger) from a completely different direction: opus reached it by counting slots, this by exhibiting 585
pairwise-superior non-tight families. Two routes, one wall.

## Status of the transfer programme

- S322: SC-keyed transfers are dead *a priori* (tournament `Z₂` has fixed points, LRC's is free).
- S323: the QR/NQR candidate is dead empirically (LRC extremals are QR-agnostic).
- S324 (this): the danger-variance candidate **passes** the AP-vs-random filter but fails as a detector,
  and its failure rules out the whole pairwise class.

Sharper filter for the next candidate: it must distinguish the AP from a random covering set **and** from
`{1,…,12,14}` — the resonance artefact that beats the AP on every pairwise statistic. Anything that cannot
separate those two is measuring coherence, not tightness.

*Files: `04-computation/lrc_dangervar_klein_S324.py` (+ 4 `.out`). Extends mac-mini-S89, THM-133,
THM-739/1012 (pairwise overlaps), opus THM-1026.*
