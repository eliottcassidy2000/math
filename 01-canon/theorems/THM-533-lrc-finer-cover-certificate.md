---
id: THM-533
title: The finer-cover certificate for LRC(14) — replacing the seven fixed sectors by L equally-spaced length-1/7 arcs (N⊆S_L⊆S_7, meas(S_L)↓meas(N)) buys ~5× slack and makes the crude relation-height certificate CLOSE at L=14 where it failed at L=7; the universal weight bound W(E)≤W(consec_k) is elementary (gap monotonicity); reduces the dangerous rows k=8..11 to ONE analytic bound corr_L(E)≤C_L·W(E)
status: MIXED. PROVED: N(E)⊆S_L(E)⊆S_7(E) and meas(S_L) decreasing to meas(N) (so meas(S_L)≤cap_k ⟹ LRC); the universal weight bound W_raw(E)≤W_raw(consec_k) (elementary, term-by-term via e_j−e_i≥j−i). VERIFIED (exact/thorough-bank): the finer-cover certificate M_L+C_L·W(consec_k)≤cap_k CLOSES at L=14 for the dangerous rows k=8,9,10 (and k=11), margin ~1.6–2×, whereas L=7 (S_7) FAILS (bound 0.408>cap 0.382); consec maximizes meas(S_L) and the primitive W. OPEN (the one analytic piece): the rigorous corr_L(E)≤C_L·W(E) (the L-arc sector-Fourier absolute bound, HYP-2601 flavor) with explicit C_L; an exact M_L; primitive-W extremality; k=12,13 (the codex HYP-2604 AP-frontier). LRC(14) NOT proved.
source: mac-mini-2026-06-18-S7
depends_on:
  - THM-532   # the seven-sector relation-height split (this fixes its honest gap)
  - HYP-2603  # codex's seven-sector net-cap reduction
  - HYP-2602  # the 1/7-spread bound (the crux)
  - HYP-2601  # the high-relation-height absolute certificate (the C_L analytic model)
related:
  - HYP-2604  # the AP-frontier envelope (k=12,13)
  - OPEN-Q-108
external: Lonely Runner Conjecture (first open case = 13 speeds). Beurling–Selberg extremal functions; Bonferroni inequalities.
---

# THM-533 — The finer-cover certificate

THM-532's seven-sector route reduced LRC(14)-S3 to `meas(S7(E)) ≤ cap_k` but left an honest
gap: the crude relation-height certificate `corr ≤ C·W` does **not** close (`C*·W(consec_8) =
0.384 > margin_8 = 0.357`) because `S7` — seven *fixed* sectors — is the **crudest** finite
cover of the 1/7-net, leaving only `0.054` slack at the extremal AP. This theorem fixes that
by refining the cover.

## A. The finer cover (PROVED reduction)

For `L ≥ 7` let `S_L(E) = {x : every arc A_j=[j/L, j/L+1/7) (j=0..L−1) is hit by some
frac(e_i x)}`. Since a 1/7-net hits *every* length-1/7 arc,
> `N(E) ⊆ S_L(E) ⊆ S_7(E)`, and `meas(S_L) ↓ meas(N)` as `L→∞`.

So `meas(S_L(E)) ≤ cap_k` still implies the crux `μ_{1/7}(E) ≥ thr_k` (HYP-2602), hence
LRC(14)-S3. The finer cover is a **tighter** upper bound on the net, with **more slack**
against `cap_k`. Exact (k=8, consecutive cluster):

| cover | meas(S_L(consec_8)) | slack `cap_8 − ·` |
|---|---|---|
| `S_7` | 0.327 | 0.054 |
| `S_14` | 0.196 | 0.185 |
| `S_21` | 0.157 | 0.225 |
| `S_42` | 0.107 | 0.275 |
| `meas(N)` | 0.060 | 0.322 |

`meas(S_L)` is scale-invariant; the iid main term `M_L → 0` as `L` grows (a dissociated
orbit threads `L` arcs ever more rarely).

## B. The certificate CLOSES at L=14 (VERIFIED, the advance)

Write `meas(S_L(E)) = M_L + corr_L(E)`. The correction is the offset relation-lattice tail,
bounded (support-3 dominant, exactly HYP-2601's calculation in the `L`-arc Fourier basis) by
`corr_L(E) ≤ C_L·W(E)`, `W(E) = Σ_{primitive triples} 1/H`. Hence
`meas(S_L(E)) ≤ M_L + C_L·W(E) ≤ M_L + C_L·W(consec_k)` (part C). The closure
`M_L + C_L·W(consec_k) ≤ cap_k` (true iid `M_L`, worst `C_L` over a thorough 300–500-shape
bank incl. all perforations):

| k | L=7 | L=14 |
|---|---|---|
| 8 | bound 0.408 > cap 0.381 — **FAILS** | bound 0.196 < cap 0.381 — **CLOSES** (`W0=19` vs `W(consec)=9.7`) |
| 9 | 0.524 > 0.494 — FAILS | 0.341 < 0.494 — CLOSES |
| 10 | 0.50 (cert=meas) | 0.389 < 0.604 — CLOSES (`W0=26.7` vs `16.9`) |
| 11 | — | CLOSES |

So **the finer cover converts the failed S_7 certificate into a closing one** with a
comfortable `~1.6–2×` margin (`W0/W(consec)`). No sharp measure-extremal and no finite
low-height residual are needed — only crude bounds.

## C. The universal weight bound (PROVED for the raw weight)

> **`W_raw(E) ≤ W_raw(consec_k)` for every set `E` of `k` distinct integers** — elementary.
> *Proof.* For sorted distinct integers `e_0<…<e_{k−1}`, `e_j−e_i ≥ j−i` (each step ≥1), so
> every triple's height `(e_j−e_i)(e_l−e_j)(e_l−e_i) ≥ (j−i)(l−j)(l−i)`, the consecutive
> triple's height. Summing `1/height` over the common index-triples gives the bound
> term-by-term. ∎

This is exactly the extremal piece the `S7` route needed ("consecutive is extremal"), now a
one-line monotonicity. The *primitive*-weight version `W(E)≤W(consec_k)` (the one the
singular series actually uses, with the `ζ`-factor from relation multiples) is VERIFIED (0
violations over 3000 random shapes per `k=8..11`) but not yet given the one-line proof (the
gcd of a triple complicates the term-by-term step); it is a clean combinatorial extremal,
strictly easier than the measure extremal.

## D. What remains (one analytic bound)

The dangerous rows `k=8..11` now reduce to the single rigorous statement
> `corr_L(E) = meas(S_L(E)) − M_L ≤ C_L·W(E)` with an explicit constant `C_L` (the `L`-arc
> sector-Fourier absolute bound — a support-3 sum + a geometric support-≥4 tail, the HYP-2601
> calculation re-run for the `L`-arc test function),

together with: an exact `M_L` (replace the Monte-Carlo value by the inclusion–exclusion
closed form), and the primitive-`W` extremality (C). Because the closure margin is now
comfortable (`~2×`), even a *loose* rigorous `C_L` suffices — the opposite of the `S7` knife-
edge. `k=12,13` stay on the codex HYP-2604 AP-frontier (large slack, separate handling).

## Net

Refining the test set from 7 fixed sectors to `L=14` equally-spaced length-1/7 arcs buys back
the slack that `S7` threw away, turning THM-532's failed crude certificate into a closing one
and reducing the extremal obligation from a measure inequality to the elementary weight
monotonicity `W_raw(E)≤W_raw(consec)`. The dangerous rows of LRC(14)-S3 now hinge on **one**
HYP-2601-style absolute bound, with margin to spare. **LRC(14) is NOT proved.**
