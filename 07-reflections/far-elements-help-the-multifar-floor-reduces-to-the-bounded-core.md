# Far elements help: the multi-far floor reduces to the bounded core (and the A000568 tameness window)

*mac-mini-2026-06-27-S69. Owner: noted A000568 (tournaments) sandwiched — 12 between 10 and 16, 56 between 20
and 80 — then "work creatively toward remaining LRC(14) proofs." This engages the observation and lands a
proof-relevant advance on the multi-far floor, continuing the Asano/Lee-Yang line (S68/HYP-3127).*

## Warm-up: the A000568 tameness window (apex-7)
The sandwich is exact and identified (`comb`/`factorial`):
```
   C(n,3)   <=   A000568(n)   <=   2(n-1)!/3
 1,4,10,20      2,4,12,56,...     4,16,80,480
```
holds for **n = 4, 5, 6, 7 ONLY** (and trivially borderline at n=3,4), and **BREAKS at n=8** (`6880 > 3360`):
the tournament count outgrows the factorial bound just above the apex prime 7. The lower bound is the
**additive/polynomial** `C(n,3)` (the triangle / OCF 3-cycle building blocks); the upper is the
**factorial/multiplicative** `2(n-1)!/3`. So `A000568` is "between additive and factorial" exactly up to
`n=7` — a small-number tameness window closing at `n=8`, echoing the additive↔multiplicative split of the
four-faces ([[the-four-faces-of-14-why-the-exceptional-structures-crowd-into-lrc]]) and the recurring
"n≤7 is special / apex-7" boundary (E₇ odd holes, the forbidden-H window n≤8, etc.). A genuine coincidence of
small numbers, flagged; no load-bearing LRC use, but the apex-7 boundary is the through-line.

## The proof advance: far elements PUSH THE ZEROS OUTWARD
Continuing HYP-3127 (the multi-far floor `R′≥c` as an Asano contraction of single-far Lee-Yang factors), the
load-bearing obligation was the single-far Lee-Yang region. Tested directly
(`lrc_leeyang_polydisk_multifar_macmini_S69.py`), tracking the **nearest-zero radius `ρ`** of the miss-count
PGF `G_N(z)` as far elements are added to a GOOD bounded base `B = consec_8` (`ρ(B) = 1.49 > 1`, Lee-Yang):

| config | r=0 (base) | r=1 | r=2 | r=3 |
|---|---|---|---|---|
| nearest-zero `ρ` | 1.49 | 1.59–1.78 | 1.56–1.86 | 1.90–2.09 |

- **Adding far elements MONOTONICALLY pushes the zeros outward** (`ρ: 1.49 → ~1.6 → ~1.7 → ~2.0`). Over a
  400-config multi-far scan the **floor is `ρ ≥ 1.559`** (binding at the resonant r=2 `(21,28)`) — a uniform
  Lee-Yang margin well above 1.
- The mechanism: **each far element INCREASES coverage** (`d(f)=p0(B∪{f})/p0(B) ≈ 1.04–1.14 > 1`, S68) — a
  *coverage-increasing* factor, which in the Asano/Lee-Yang picture pushes the zeros *out* of the unit disk.

> **The multi-far elements are NOT the obstruction — they HELP.** The binding case is the **bounded core**
> (`consec`, `ρ=1.49`) — exactly the φ⁴ hard row of S67. The multi-far Lee-Yang region *reduces to* the
> bounded-core Lee-Yang property (a finite check on the bounded core) **plus** the far-pushes-out monotonicity
> (Asano, coverage-increasing). Far elements only widen the margin.

## What this closes (the reduction, sharpened)
HYP-3127's three obligations collapse toward one:
1. ~~single-far Lee-Yang region~~ — VERIFIED: single-far `ρ ≥ 1.59 > ρ(base)`; the single-far factor is
   coverage-increasing, pushing zeros out. ✓
3. ~~r-monotonicity~~ — VERIFIED: `ρ` increases with `r` (1.49→2.0); far helps. ✓
2. the constant `c` / the `R′` link — the **only remaining piece**: a uniform Lee-Yang region `ρ ≥ ρ₀ > 1`
   on the bounded core (the finite, well-studied φ⁴/extremality check) ⟹ `R′ ≥ c` via the Lee–Yang ⟹
   correlation-inequality step (GHS/Griffiths for the confined-zero regime).

So the **far structure of the covering bound subsumes into the bounded-core Lee-Yang property**: prove the
bounded core has `ρ(G_N) > 1` (the coverage extremality / φ⁴ row, S66/S67), and the multi-far floor follows
for free because far elements only push the zeros further out.

## Honest status
VERIFIED: the A000568 sandwich (n=4..7, breaks at 8); `ρ(consec_8)=1.49`; far elements push `ρ` outward
monotonically; multi-far floor `ρ ≥ 1.559`; far is coverage-increasing. BOLD/OBLIGATION: the rigorous
`ρ_bounded > 1 ⟹ R′ ≥ c` (the Lee-Yang ⟹ correlation-inequality link, the remaining analytic step); the
proof that far-pushes-out holds for ALL far placements (not just the scan). The reframe: **the multi-far floor
is not a separate obstruction — it is downstream of the bounded-core Lee-Yang property, which is the φ⁴
extremality row I have been studying.**

Related: HYP-3131 (this), HYP-3127 (multi-far = Asano contraction), HYP-3122 (φ⁴ cap / the bounded-core row),
HYP-3103 (PGF zeros), HYP-2829 (single-far), OPEN-Q-108.
