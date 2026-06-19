---
id: THM-532
title: The seven-sector relation-height split for LRC(14) — meas(S7(E))=M7(k)+corr(E) with a TINY main term M7(k)≪cap_k (certificate margin ≥0.357), corr the offset relation-lattice tail bounded by C·W(E) (short-relation weight, the THM-503 7-vanishing triple bound); consecutive maximizes BOTH meas(S7) and W; finite check meas(S7(consec_k))<cap_k all k=8..13; reduces LRC(14) to a high-relation-height sector certificate + a finite low-height (AP-rich) residual
status: MIXED. PROVED: meas(S7) scale-invariance; M7(k) closed form + M7(k)≪cap_k (exact); the finite check meas(S7(consec_k))<cap_k for all k=8..13 (exact rational). VERIFIED (exact/exhaustive): meas(S7(E))≤cap_k everywhere tested (codex HYP-2603); consecutive maximizes meas(S7) (k=8,9,10 exhaustive) and W(E) (k=8..11); corr(E)≤C·W(E) with C*≈0.0395; the near-binding residual at k=8 is ONLY consecutive. PROOF PROGRAM (not closed): the high-relation-height sector certificate (corr≤C·W, rigorous absolute bound to be written) + a finite low-height residual ({W>W0}, AP-rich, finite mod scale, all ≤cap_k). LRC(14) NOT proved.
source: mac-mini-2026-06-18-S6
depends_on:
  - HYP-2603   # codex's seven-sector net-cap reduction meas(S7)<=cap_k => HYP-2602
  - HYP-2602   # the 1/7-spread bound (mu_{1/7}(E)>=thr_k) — the crux S7 strengthens
  - THM-531    # the exact 1/7 union closure + AP-orbit invariance
  - THM-530    # the global-witness two-branch floor
related:
  - HYP-2601   # the high-relation-height certificate B(E)<(5/7)^k (the 5/7-window analog)
  - HYP-2600   # the subtorus-theta singular series (codex/mac-mini-S4)
  - OPEN-Q-108
external: Lonely Runner Conjecture (first open case = 13 speeds). Steinhaus three-gap; Weyl/Erdős–Turán discrepancy.
---

# THM-532 — The seven-sector relation-height split

Codex's **HYP-2603** reduces the LRC(14) crux to a clean combinatorial bound. Let
`S7(E) = {x∈[0,1) : every fixed sector [j/7,(j+1)/7), j=0..6, is hit by some frac(e x), e∈E}`
be the **seven-sector cover set**. Since any `1/7`-net `N(E)={maxgap≤1/7}` must hit every
sector, `N(E) ⊆ S7(E)`, so `meas(S7(E)) ≤ cap_k := min_{|P|=13−k}meas(G_P)` implies
`μ_{1/7}(E) = 1−meas(N(E)) ≥ 1−cap_k = thr_k` (HYP-2602) implies LRC(14)-S3. This theorem
makes the **relation-height split** of `meas(S7)` concrete.

## A. The main term is tiny — a huge certificate margin (PROVED)

By inclusion–exclusion over missed sectors (sector `0` is always hit, since `0∈E` pins a
point at `0`),
`meas(S7(E)) = Σ_{T⊆{1..6}}(−1)^{|T|}∫_0^1 ∏_{e∈E}1_{B_T}(ex)dx`, `B_T = [0,1)∖∪_{j∈T}S_j`.
The `n=0` (independent) term gives the **closed-form main term**
> `M7(k) := Σ_{T⊆{1..6}}(−1)^{|T|}(1−|T|/7)^{k−1} = Σ_{t=0}^{6}(−1)^tC(6,t)((7−t)/7)^{k−1}`.

Exact values: `M7(8)=20160/823543≈0.0245`, rising to `M7(13)≈0.297`. Crucially
**`M7(k) ≪ cap_k`**, giving the **certificate margin** `margin_k := cap_k − M7(k)`:

| k | M7(k) | cap_k | margin_k |
|---|---|---|---|
| 8 | 0.0245 | 0.3815 | **0.357** |
| 9 | 0.0577 | 0.4943 | 0.437 |
| 10 | 0.1049 | 0.6044 | 0.499 |
| 11 | 0.1631 | 0.7253 | 0.562 |
| 12 | 0.2285 | 0.8571 | 0.629 |
| 13 | 0.2973 | 1.0000 | 0.703 |

(Contrast HYP-2601's `(5/7)^k` budget — `(5/7)^{13}=0.013` — the sector margin is `~30×`
larger.) **`meas(S7)` is scale-invariant**: `meas(S7(dE))=meas(S7(E))` (substitute
`y=frac(dx)`, measure-preserving), so dilated clusters reduce to their primitive shape.

## B. The correction is the relation-lattice tail ≈ C·W(E) (the height split)

`corr(E) := meas(S7(E)) − M7(k)` is the sum of the `n≠0` terms over the **offset relation
lattice** `Λ°(E)={n: Σ_{e≠0}n_e e=0}`. Its dominant (support-3) part is bounded, exactly as
in HYP-2601, by a sum over primitive offset triples of the sector-Fourier products; the
sector coefficients carry the **THM-503 7-vanishing** (`â_T(7m)=0`), and decay like
`1/|n|`, so the per-triple contribution is `≲ 1/H_{triple}` (`H=|abc|` = relation height):
> `corr(E) ≤ C · W(E)`,  `W(E) := Σ_{primitive triples}1/H_{triple}`  (short-relation weight).

Verified (k=8): `corr` tracks `W` cleanly — `consec: corr 0.303, W 9.73`; `perforated:
0.099, 8.5`; `dissociated: 0.010, 1.4`; `generic: 0.003, 0.76`. The empirical constant is
`C* = max corr/W ≈ 0.0395`. **High relation-height (`W` small) ⟹ `corr≈0` ⟹
`meas(S7)≈M7(k)≪cap_k`, certified with the full margin.** This is the user's "high-height
sector discrepancy": a dissociated cluster's orbit rarely fills all seven sectors at once.

## C. Consecutive is the extremiser of BOTH meas(S7) and W (VERIFIED)

`meas(S7(consec_k)) < cap_k` for **every** `k=8..13` (exact; slack `0.054 → 0.324`), and
**consecutive maximises `meas(S7)`** (exhaustive `k=8,9,10`, `0` shapes exceed it) and the
short-relation weight `W` (`k=8..11`). The near-binding residual at `k=8` (`meas(S7) ≥
0.85·cap_8`) is the **single** shape `{0..7}` — the AP is the unique extremal, and
everything else falls away.

## D. The proof program and its honest gap

The split gives a two-piece program for `meas(S7(E)) ≤ cap_k` (= codex's
"finite low-height patterns + high-height tail"):
- **High-height tail (certificate):** `C·W(E) ≤ margin_k ⟹ meas(S7)≤cap_k`. Rigorous once the
  absolute bound `corr ≤ C·W` is written with the explicit sector-Fourier `C` (the support-3
  sum + a geometric tail for support ≥4). Fires for all but the most relation-rich shapes.
- **Low-height residual (finite):** `{E : W(E) > W0_k := margin_k/C}` is **AP-rich**; since a
  large `W` forces many low-height triples, it is **finite modulo scale+translation** (the
  relation pattern bounds the offset ratios). Each is checked `≤ cap_k` (consecutive is the
  max, e.g. `0.327 < 0.381` at `k=8`).

**Honest gap.** The *crude* product bound does NOT close on its own: `C*·W(consec_8) = 0.384 >
margin_8 = 0.357` (the worst ratio `C*` and the max weight `W(consec)` are attained at
different shapes, so their product overshoots the true `corr(consec)=0.303`). Hence the
finite low-height check is genuinely needed — but the residual is a *narrow* band near the
single AP (at `k=8`, only `W∈(W0,9.73]`), not an infinite family. Remaining to finish:
(i) write the rigorous `corr ≤ C·W` with explicit `C`; (ii) prove `{W>W0}` is finite mod
scale and enumerate it; (iii) the upstream finite-`Vmax` glue (THM-527-A/HYP-2589).

## Net

The seven-sector cover is the cleanest carrier of the residual found so far: a *closed-form
tiny main term*, a *relation-height-graded correction* with the explicit short-relation
weight `W(E)`, and `consecutive` the unique extremiser of every relevant functional. The
crux `μ_{1/7}(E)≥thr_k` (HYP-2602) is now a **high-relation-height sector certificate plus a
finite AP-rich residual** — codex's HYP-2603 program, concretised with the margin quantified
and the correction split. **LRC(14) is NOT proved.**
