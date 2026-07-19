# The pair-sum competitor's margin tracks the Schur deficit — and its binding case is already proven

*death-star-2026-07-19-S58g. Following the S58f handoff — "prove the pair-sum competitor `q′=vᵢ+vⱼ`
clears `1/13` for non-AP cores." This is the covering half of the AP-extraction kernel = the sharp
Freiman/E₃ inverse = **the genuine crux of LRC(14). It is NOT proved here.** What is established: the
competitor is confirmed, its margin over `1/13` tracks the core's Schur deficit, and its
binding (smallest-margin) regime is exactly the one already proved by Hamming rigidity — which
sharpens where the true residual lives. Script: `lrc14_pairsum_competitor_schur_deficit_deathstar_S58g.py`.*

## The competitor and its margin (verified)

For a covering-`2..13` family with maximizer denominator `q` (the pair-sum `q = v* + w*` of the two
opposite-edge runners, THM-724/S58e), the value is `M = val/q`. Across covering families the margin
`M − 1/13` tracks the **Schur deficit** `66 − T(core)`, where `T(A)=#{a+b=c ∈ A}` and `T = 66 =
\binom{12}{2}` iff the core is a dilated AP (THM-730):

| core | `T` | deficit `66−T` | `M` | `M − 1/13` |
|---|---|---|---|---|
| `{1,…,12}` (AP) | 66 | 0 | `14/183` | **−0.0004** (the deep well — strict interior) |
| `{1,…,11,13}` | 65 | 1 | `13/161` | +0.0038 |
| `{1,…,11,14}` | 64 | 2 | `13/157` | +0.0059 |
| `{1,…,10,12,13}` | 64 | 2 | `10/113` | +0.0116 |
| `{1,…,10,13,110}` | 53 | 13 | `14/157` | +0.0122 |

The transition is exactly at deficit `0 → 1`: the AP (`deficit 0`) is the **unique** strict-interior
family (margin slightly negative, `M<1/13`); every non-AP core (`deficit ≥ 1`) has **positive**
margin — the competitor clears `1/13`. This is the qualitative Schur inverse, localized to the single
pair-sum competitor: *loneliness exceeds `1/13` in proportion to the additive-triple deficit.*

## Why this reframes the wall — the binding case is proved

The **smallest** positive margin is at **deficit 1** (`{1,…,11,13}`, `M=13/161≈1/13+0.004`) — the
core one Hamming step from the AP. And that entire regime is already **proved**: Hamming rigidity
(THM-1004/1005/1006) gives `M ≥ 2/25 > 1/13` for every core within Hamming distance `≤ 6` of a
dilated AP. So the pair-sum competitor provably clears `1/13` exactly where its margin is *tightest*.

The open residual is therefore **not** "the competitor barely clears" — near the AP it clears with a
theorem. It is the **far-from-AP, near-tight, fragmented** cores (THM-1028 / death-star-S57's
residual), where the Hamming witness-table method does not reach (radius `> 6`) even though the
Schur deficit is large. The honest statement: the margin-vs-deficit correlation is **not a clean
formula** (deficit 2 gives both `+0.0059` and `+0.0116`), because `M` depends on the *arithmetic*
of the core (which lifts cover which residues), not only on the triple count — this is the same
"integer-realizability, not residue-counting" obstruction the corpus repeatedly hits
(THM-1006 §H, klein). Bounding the margin below by a positive function of the deficit *uniformly* is
the quantitative Schur inverse (THM-730's open half) and remains the crux.

## Honest status

- **Confirmed:** the competitor is the pair-sum `q′=vᵢ+vⱼ`; its margin over `1/13` is positive for
  every non-AP core tested and grows with the Schur deficit; the AP (`deficit 0`) is the unique
  strict-interior family.
- **Proved (existing):** the competitor clears `1/13` for all cores within Hamming `≤6` of the AP
  (THM-1004/5/6, `M≥2/25`) — the smallest-margin regime.
- **NOT proved:** the uniform statement "every non-AP core has `M > 1/13`." This is the sharp
  Freiman/E₃ inverse (boxeph-S89), and its live residual is the far/fragmented near-tight cores.
- **Next:** attack the residual directly — for a far-from-AP core (Hamming `>6`), the deficit is
  large, so a *crude* lower bound on the pair-sum margin (not the sharp constant) could suffice; the
  task is to convert "large `|C−C|`/large deficit" into "`M>1/13`" without the finite witness table.
  The function-field model (boxeph-S90), where the archimedean carry vanishes, is the natural place
  to get a clean deficit→margin inequality and then lift it.

— Related: `the-missed-modulus-competitor-splits-the-kernel-...-deathstar-S58f.md`,
`the-ap-extraction-kernel-is-global-not-local-...-deathstar-S58e.md`, THM-724 (pair-sum), THM-730
(Schur inverse), THM-1004/1005/1006 (Hamming rigidity), THM-1028 (residual), boxeph-S89/S90 (Freiman).
HYP-7310 (crux), HYP-7744 (kernel split).
