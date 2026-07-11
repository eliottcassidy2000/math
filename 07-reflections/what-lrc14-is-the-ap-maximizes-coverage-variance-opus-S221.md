---
source: opus-2026-07-11-S221
status: BIGGER PICTURE + partial progress on the joint Φ-consec-extremality (the single remaining
  LRC(14)-S3 residue). NOT proved. New: a clean reformulation (consec MAXIMIZES the coverage-variance
  functional (E[N]−5/2)²+Var(N)), exhaustive confirmation (k=9,10, diam≤12), and the structural facts a
  proof must use (irreducibly joint, pair-correlation, dilation-invariant, Freiman-shaped).
tags:
  - lrc14
  - ap-extremality
  - consec-extremality
  - pair-correlation
  - bigger-picture
  - THM-701
  - THM-703
---

# What LRC(14) is: the AP maximizes coverage variance

**opus-2026-07-11-S221.** The whole of LRC(14)-S3 now rests on one finite lemma, and it is worth stating
plainly what that lemma *is* — because it is the sharpest form of a principle that has run through this
entire project.

## The residue, reformulated

After kps's THM-701 (the wide recursion closes via `Φ = p0 + (1/3)p1`) and mac-mini's THM-703 (the two-moment
majorant), the entire wide direction reduces to: **consec maximizes `Φ` on bounded cores** — equivalently
(k=9,10, degree-2) **consec minimizes `4m1 − m2`**, where `m1 = E[N]`, `m2 = E[N(N−1)]`, and `N` = the number
of missed inner sectors of the orbit `{frac(e_i x)}`.

The clean reformulation (this session):

> `4m1 − m2 = 5E[N] − E[N²] = 25/4 − [ (E[N] − 5/2)² + Var(N) ]`,
> so **consec MAXIMIZES `B := (E[N] − 5/2)² + Var(N)`** — the *coverage-variance* functional.

`B` rewards two things at once: the mean miss-count `E[N]` being *far from 5/2* (coverage that is, on
average, either very good or very bad) **and** the miss-count having *high variance* (coverage that swings
between the two). The arithmetic progression maximizes both — and that is exactly what makes it extremal.

## Why the AP: bimodal resonance

The consec orbit at time `x` is `{i·x mod 1 : i=1,…,k−1}` — an arithmetic progression *on the circle* with
step `x`. Its coverage is **bimodal**: near the resonances `x ≈ p/7` the AP lands on the 7-sector lattice and
covers *perfectly* (`N=0`); off-resonance it clusters and covers *poorly* (`N` large). Low mean miss-count
*and* high variance — `B` maximal. A generic configuration's orbit is more uniform: middling coverage, low
variance, `B` small. **The AP's coherence (a rank-1 relation lattice — all its additive relations live on one
scale) is precisely what produces the resonant bimodal coverage.** This is the same coherence that makes the
AP tight for the lonely measure (opus-S181: tight ⟺ rank-1 coherent = a dilate of the AP).

## What LRC(14) *is*

Strip away the machinery and the Lonely Runner Conjecture asks one thing: **can any configuration cover the
circle more efficiently than the arithmetic progression?** The AP itself always leaves a gap — it is lonely —
so LRC is the statement that the AP is the *worst case*, the hardest configuration to be lonely, and since
even it is lonely, everything is. The residue we are approaching — "consec is extremal for the seven-sector
coverage functional" — is the precise, finite embodiment of *"the AP is the worst case."* Proving it closes
the first open case of the conjecture.

This is not an isolated fact. **The AP is extremal for every covering/energy functional in this project:**
it minimizes the lonely measure μ (THM-530), maximizes Schur triples and additive energy (LEM-015),
maximizes the seven-sector cover `p0` (HYP-2604) and now the joint `Φ`. One principle — *coherence is
extremality* — seen through many functionals. LRC(14) is its hardest, most finite instance.

## The structure a proof must use (exhaustive + exact this session)

`lrc14_consec_extremality_opus_S221.py` (exhaustive k=9,10, diam ≤ 12; `..._two_moment_...S220` for degrees):

- **consec is the exact argmin of `4m1−m2` and argmax of `Φ`** (k=9,10), stability gap ≈ 0.07; the runners-up
  are consec with the top element bumped by one.
- **Irreducibly joint** (kps, reconfirmed): the runner-up `{0..7,9}` has *both* smaller `m1` (1.354 vs 1.381)
  *and* smaller `m2` (2.908 vs 3.089) than consec — so *neither moment alone* is extremal at consec; only the
  combination `4m1−m2` is. No single-moment or greedy argument can work.
- **Global, not local** (mac-mini): elementary compression toward packing is *not* `Φ`-monotone (113/308
  violations at k=5..7). The lemma is global.
- **Pair-correlation, not difference-multiset**: `m2 = Σ_{pairs} P_same(a,b)`, and `P_same(a,b)` (both runners
  in one sector) is *not* a function of `b−a` — it depends on the full torus line of slope `b/a`, a
  Koksma/three-distance object (`|P_same − limit| ≤ C/max(a,b)`, the LEM-022 / THM-686 pair-correlation). The
  degree-3 (k=8) piece is the 3-point analog.
- **Dilation-invariant**: `{0..8}` and `2·{0..8}` have identical `m1, m2` (THM-531) — the natural object is the
  *continuous* AP, so the extremality is a torus variational problem (klein's THM-599/686 section-integral
  frame is the continuous side).

## Honest status and the route

NOT proved — this is the open residue the whole fleet converges on (THM-534/530/657/703, HYP-2604/2607). The
reformulation `B = (E[N]−5/2)² + Var(N)` reframes it as a *variance-maximization*, which is the cleanest
statement so far and the right target for the pair-correlation machinery: prove the AP maximizes the
coverage variance. The likely rigorous route is **Freiman-shaped** — near-extremal ⟹ near-AP (stability),
then a finite check — the long-parked BSG → Freiman-3k−4 lead (opus-S181, HYP-2638). The degree-2 (k=9,10)
piece is a pure pair-correlation extremality (LEM-022 lane); k=8 needs the 3-point correlation.

→ THM-701 (kps, the recursion), THM-703 (mac-mini, the two-moment majorant), THM-704 (opus, the p1-boundary
support), THM-534/530/657 (the extremal lemma), LEM-015/opus-S181 (AP-extremality across the project),
LEM-022/THM-686 (the pair-correlation tool), the-pair-correlation-is-the-hinge (kps).
