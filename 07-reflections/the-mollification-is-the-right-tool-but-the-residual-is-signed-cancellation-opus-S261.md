---
source: opus-2026-07-11-S261
status: EXPERIMENT (ran the Beurling-Selberg mollification of G' against the coprime core). Findings: (1) a
  FINITE degree K~50 captures the full discrepancy eps_v (fast tail decay) -- the mollification is the correct
  tool; (2) the L2 bound |eps_v| <= sqrt(tail_v)/(sqrt6*|G'|) improves the naive Erdős–Turán ~17x (< 1 per arc
  for large core, vs ~14) and reduces to tail_v = high-freq mass of G'; (3) BUT it is still ~40x too weak vs
  the actual eps ~0.02, because magnitude bounds discard the SIGNED CANCELLATION in Sum_h b_h ghat(-hv) that
  makes eps small. Runner 1 (low-freq) is the isolated exception => S255. The residual is now precisely named:
  a cancellation/bilinear estimate for the coprime core against the non-core resonance lattice.
tags:
  - lrc14
  - covering-min
  - anti-concentration
  - beurling-selberg
  - mollification
  - signed-cancellation
  - bilinear
---

# The mollification is the right tool; the residual is signed cancellation

**opus-2026-07-11-S261.** Owner: run the Beurling–Selberg mollification of `G'` against the core (the S260
refined path). Done — via FFT on the good set — and it sharpens the crux to a precise cancellation estimate.

## What the mollification shows

For a covering family, `density(D_v in G') = 1/7 + ε_v`, `ε_v = (Σ_{h≠0} b_h ĝ(−hv))/|G'|`,
`b_h = sin(πh/7)/(πh)`. Running the mollification (degree-`K` truncation + tail):

1. **A finite degree suffices.** The truncated discrepancy `Σ_{|h|≤K} b_h ĝ(−hv)` converges to `ε_v` **fast** —
   essentially exact at `K = 50` (tail beyond negligible). So a degree ~50 Beurling–Selberg majorant captures
   the *entire* discrepancy. The mollification is the correct, finite tool.

2. **Two regimes.** **Large core runners** (`v ≥ 17`, coprime): `ε_v` **small** (0.01–0.09) — their frequencies
   `hv ≥ 17` land on the *high-frequency* (small) part of `ĝ`, so they equidistribute. **Runner 1** (`v = 1`):
   `ε_1` **large** (`0.57` at the deep well) — its frequencies `h` land on the *low-frequency* (large) part of
   `ĝ`. Runner 1 is the exception; when it is the only core runner (deep well), `coreCover = density(runner 1)
   < 1` is the **near-AP case, handled by S255**.

3. **The L² bound improves the naive ~17× — but ignores cancellation.** `|ε_v| ≤ √(tail_v)/(√6·|G'|)`,
   `tail_v = Σ_{|m|≥v}|ĝ(m)|²` (high-frequency mass of `G'`). For large core `v` this is **0.4–0.9** (< 1 per
   arc, versus the naive `N/(6v|G'|) ≈ 14`). *But* the actual `ε_v ≈ 0.02` — the L² bound is still **~40×**
   too weak, because it uses `|ĝ|` magnitudes and **discards the signed cancellation** in `Σ_h b_h ĝ(−hv)`.

## The residual, precisely named

The true smallness of `ε_v` is a **signed cancellation**: `Σ_h b_h ĝ(−hv)` is small because, for `v` coprime to
the non-core, the frequencies `−hv` are *generic* relative to the resonance lattice of `ĝ` (the integer
combinations of non-core speeds), so the **signed** sum cancels — not because the terms are individually tiny.
Every magnitude bound (naive `N/m`, L² Cauchy–Schwarz) throws this away and is `~40–700×` too weak. Capturing
it needs a **bilinear / cancellation estimate** on `Σ_h b_h ĝ(−hv)` exploiting `gcd(v, non-core) = 1` — the same
*kind* of object as the fleet's **LRCFourierCompletion** completion identity (`|C_w − b²/q|`, a signed
cancellation bound) and the resolved `t≥3` signed-cancellation thread (task #41). So the crux is now a specific
cancellation bound for the coprime core against the non-core resonance lattice.

## Net (honest)

Running the Beurling–Selberg mollification: it is **the right tool** (finite degree `K~50` captures the
discrepancy), it **improves the bound ~17×** (L² `< 1` per arc for large core), it **reduces to** `tail_v` =
high-frequency mass of `G'`, and it **isolates runner 1** (low-frequency ⟹ S255). But it does **not** close the
proof: the residual is the **signed cancellation** in `Σ_h b_h ĝ(−hv)`, a bilinear estimate for `v` coprime to
the non-core lattice — exactly the arithmetic that magnitude bounds cannot see. The target remains the
independent model `coreCover ≈ 1 − (6/7)^{core} < 1` (margin `(6/7)^{core}`); the cancellation is what must be
proven to reach it. The honest state of the covering-min after the S253–S261 arc: extremizer proved (S255);
loose stratum = the coprime core's **signed-cancellation** discrepancy against `G'`, wired into the
LRCFourierCompletion cancellation machinery, with runner 1 folded into S255.

→ opus-S260 (mollification path — run here), opus-S259 (equidistribution), opus-S255 (runner-1/near-AP),
LRCFourierCompletion / task #41 (signed cancellation of this type), s558o. Files:
`lrc14_beurling_selberg_mollification_opus_S261.py` (+`.out`).
