# The G-K pinning: a withdrawn theorem, a surviving phrase, and the gridmax window

**Instance:** opus-2026-07-19-S402 (owner: do the G-K pinning session to promote HYP-7930).
**Outcome:** SPLIT promotion. THM-1289 (floor isolated from above — citation-grade, published
source); MISTAKE-190 (my S401 headline propagated a mis-extracted chain); HYP-7930-UPDATE
(window-finiteness demoted to conditional on their Conjecture 1.5, and RELOCALIZED to a
repo-native object). Controls: `lrc_gk_pinning_controls_opus_S402.py` (Fan–Sun 4/4 exact;
gridmax boundary value).

## 1. The detective story (why this citation needed archaeology)

arXiv:2304.01462 has five lives: v1/v2 (2023) claimed **acc(S(n)) = S(n−1)** as the main
theorem — the exact sentence the repo's proof-map quoted in the MISTAKE-117 era. **v3
(2023-12) is a withdrawal notice.** v4 (2025-03) is the repaired version, now PUBLISHED
(Math. Proc. Camb. Phil. Soc., DOI 10.1017/S0305004125101497), which proves the weaker
**Theorem 1.4** and demotes the equality to **Conjecture 1.5**. Their footnote 2 names the
error: intersecting subtori with coordinate subspaces produces *disconnected subgroups*,
not subtori — exactly the star-index subtlety this pinning session had flagged in advance
as the make-or-break question. The /abs fetch that started HYP-7930 served the withdrawn
abstract; the repo has now quoted a withdrawn theorem twice, thirteen days apart, from two
different directions (MISTAKE-117's sup-misuse; my S401 acc-use of the equality). The
difference this time: the pinning process was already scheduled, and it caught the
demotion within one day (MISTAKE-190).

## 2. What is actually proven, and what survives of the corollaries

**Proven (published Theorem 1.4):** S_k(n) and S*_{k'}(n) have **only upper accumulation
points** (x is upper-acc if S ∩ (x, x+ε) ≠ ∅ ∀ε; lower analogous), and
acc(S_k(n)) ⊆ S_{k+1}(n) ⊆ S*_{k+1}(n) = S*₀(n−k−1). Mirror cross-check: their Theorem
1.3 (Kravitz) has the tilde (M-side) spectrum with lower accumulation ⊇ S̃(n−1) — the
rising s/(ns+1) ladders — and S(n) = 1/2 − S̃(n): from-below is the allowed direction,
from-above the forbidden one. Internal consistency across two theorems and the abstract:
the extraction is as pinned as HTML allows (residual: one human glance at the printed
phrase; flagged in THM-1289's status).

- **(C1) SURVIVES, UPGRADED — THM-1289:** in M-language, NO value is approached from
  above by 13-speed families; (1/14, 1/14+δ) is empty (δ ineffective, all heights).
  THM-1268's "no gap at all" horn is dead. Likewise (1/13, 1/13+δ′) at 12 speeds — the
  ineffective all-heights sibling of HYP-7920's effective height-258,276 micro-gap.
- **(C2)/(C3) DEMOTED:** full finiteness of (1/14, 3/41] (resp. the 12-speed gap as a
  finite list) needs Conjecture 1.5, since the proven chain's terminal set S*₀(n−2) is
  **grid-loneliness land** — D-values of finite proper subgroups of the (n−2)-torus —
  which is not floored at 1/13: gridmax((1,…,11); 14) = 1/14 exactly (verified), and
  interior window values are realizable.

## 3. The relocalization (what the correction buys)

Under the published theorem, any accumulation inside (1/14, 1/13) is from below only, and
must sit AT a point of **G*₀(11) := the grid-loneliness value set of finite proper
subgroups of the 11-torus**. So the window question localizes to a repo-native object:

> **NEW LEVER: characterize G*₀(11) ∩ (1/14, 1/13 − ε).** If it is finite or structurally
> thin (rung-quantized — these values ARE certificate-rung objects), then
> near-finiteness of the window follows from the PUBLISHED theorem alone, without
> Conjecture 1.5. The repo's gridmax machinery, certificate-rung profiler, and (D,s)
> quantization are exactly the tools; nobody outside this repo is better equipped.

Second lever: the paper states the sufficient condition for Conjecture 1.5 — **S₂(n) =
S(n−1)**, a rank-2-subtorus statement. The repo's dusty Route-2 Lean artifacts
(LRCRankRigidity, LRCTorusProjection, torus_loose_of_rank2) are machinery aimed at
precisely this object; the corpse of the MISTAKE-117 route is equipment for the
legitimate revival's missing half.

## 4. Bonus pins (the spectrum program IS the repo's rung program)

- **G-K Conjecture 1.1 (Kravitz's original) is already refuted in the literature** by
  Fan–Sun: ML(8, 4r+3, 4r+11, 4r+19) = (2r+7)/(8r+30) — k=2 rungs at n=4. Verified 4/4
  exactly this session (doubling as a positive control on the repo scanner). The live
  conjecture is **Amended 1.2: S̃(n) ∩ (0, 1/n) ⊆ {s/(ns+k) : 1 ≤ k ≤ n}** — which is
  the repo's (D, s)/slack program verbatim (their s = our D; their k = our slack;
  THM-1205/1210/1230/1235/1268/1269 are its local theory at n = 13).
- **death-star's D-graded gate tower is the k = D−1 diagonal of Conjecture 1.2's table**
  (4/127: k=3 at n=31; …; 12/116396303: k=11 at the 9699690-primorial) — systematic
  realizability data far beyond Fan–Sun's k=2, with eight out-of-sample confirmations.
  Their Theorem 1.3 conversely predicts the from-below pile at each (n−1)-value — the
  escape-family behavior the repo has observed at 1/13 all along.
- The note-to-authors package (backlog, owner decision) now spans: C7.1 refuted +
  repaired version (THM-1288), the tightness cage (HYP-7920), the I(13,p,1) bridge, the
  k=D−1 realizability tower, and the G*₀-window localization — a genuine two-way research
  exchange, all sourced.

## 5. Cross-links

THM-1289 (the promotion) · MISTAKE-190 · HYP-7930 + UPDATE · THM-1268 (horn-2 banner) ·
HYP-7920 (the effective sibling) · THM-1288 (C7.1) · MISTAKE-117 (the boundary, now fully
mapped: sup-misuse retracted, acc-use pinned to the published version) · death-star
S59/THM-1256/1276/1287 (the k=D−1 diagonal) · Fan–Sun (proof-map prior anchor) ·
`lrc_gk_pinning_controls_opus_S402.py` + frozen out.
