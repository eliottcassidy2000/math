# The spine stays quadratic `Q(√-p)`, but the covering field is `Q(√-3)` — they coincide only at n=6, and the split at n=14 is the difficulty

*kind-pasteur-2026-07-01-S20. Pushing the trace thread: do the Paley "spine" eigenvalues stay quadratic at larger apices, and do their fields match the covering continued fraction? Answers, verified: **yes, quadratic `Q(√-p)` at every `p≡3 mod 4`**; and their field matches the covering's `Q(√-3)` **only at `p=3` (n=6)**. At `n=14` the certification field `Q(√-7)` and the covering field `Q(√-3)` diverge — and that field split is a clean arithmetic statement of why LRC(6) is proved and LRC(14) is open.*

## The spine stays quadratic

The Paley skew "spine" is the circulant `S_{jk}=χ(k−j)` on `Z/p` (`p≡3 mod 4`); its eigenvalues are `χ(m)·g_p` for `m≠0` and `0` for `m=0`, where `g_p = i√p` is the quadratic Gauss sum (my S20 / opus-S23). Verified for `p = 3,7,11,19,23,31`:

> nonzero eigenvalues `= ±i√p` on the nose (`±1.732, ±2.646, ±3.317, ±4.359, ±4.796, ±5.568`), all **degree-2 over `Q`, field `Q(√-p)`**.

So the arithmetic spine is quadratic at every apex — it never leaves an imaginary quadratic field. (Contrast the *combinatorial* spine, the SC blue-graph of the merged metagraph: connected but a large graph, whose spectrum is high-degree at `n≥5` — the arithmetic spine stays quadratic, the combinatorial one does not.)

## The covering field is `Q(√-3)`, fixed

The covering-min at `n=2p` is `M = n/Φ₆(n)`, `Φ₆(n)=n²−n+1` (`1/M = [n−1;\,n]`, a *finite* CF ⇒ `M` is **rational**). The arithmetic content is `Φ₆`, the sixth cyclotomic polynomial, whose splitting field is `Q(ζ₆)=Q(√-3)` — the Eisenstein field the repo attaches to the covering-min (mac-mini-S77/klein-S68's Φ₆ Eisenstein-CRT). **This field is fixed at `Q(√-3)` for every `n`**, independent of the apex.

## They match only at n=6 — and the split at n=14 is the difficulty

| p | n=2p | spine field | covering (Φ₆) field | match? |
|---|---|---|---|---|
| **3** | **6** | `Q(√-3)` | `Q(√-3)` | **YES → LRC(6) PROVED** |
| 7 | 14 | `Q(√-7)` | `Q(√-3)` | no → **first open** |
| 11 | 22 | `Q(√-11)` | `Q(√-3)` | no |
| 19 | 38 | `Q(√-19)` | `Q(√-3)` | no |
| 23 | 46 | `Q(√-23)` | `Q(√-3)` | no |
| 31 | 62 | `Q(√-31)` | `Q(√-3)` | no |

**The certification field `Q(√-p)` equals the covering field `Q(√-3)` exactly when `p=3` — i.e. `n=6`, the proved case.** At `n=14` they diverge: the loneliness is *certified* by the Gauss sum in `Q(√-7)`, but the covering *geometry* lives in `Q(√-3)`. A proof of LRC(14) must therefore work in the **biquadratic compositum `Q(√-3, √-7)`** and bridge the two fields — whereas LRC(6) needs only one field. This is a new, clean way to say why `n=14` is the first hard case: **it is the first `n=2p` where the Gauss-sum certificate leaves the Eisenstein covering field.**

## The deep-well CF sees the compositum

The *periodic* deep-well continued fraction `[0;\overline{n−1,n}]` (the Herman-rigid, badly-approximable limit of the binding time, HYP-3796) is a **real** quadratic `Q(√D_n)`, `D_n = n(n−1)(n²−n+4)`. At `n=14` its squarefree part is `3·7·13·31` — and `3·7 = 21` is exactly the **real compositum** of `Q(√-3)` and `Q(√-7)` (`√(−3)·√(−7)=−√21`, and `21 | 8463 = 21·403`). So the deep-well CF field contains `√21`, the real quadratic subfield of `Q(√-3,√-7)` — it *sees both* the covering field and the certification field at once. (The apex `p` always divides `D_n` since `p ∣ n ∣ D_n`; the `√(3p)` compositum structure is clean at `n=14` but not universal — at `n=46` the factor `3` appears squared (`45=9·5`) and drops from the squarefree part, so the deep-well field there is `Q(√(5·17·23·61))`, missing `√3`. The compositum reading is exact where `3 ∥ D_n`.)

## The picture

Three quadratic fields at each `n=2p`:
- **Certification / Gauss-sum spine:** `Q(√-p)` — imaginary, varies with the apex, always degree 2.
- **Covering-min / Φ₆ / Eisenstein:** `Q(√-3)` — imaginary, fixed.
- **Deep-well binding CF:** `Q(√D_n)` — real, containing (generically) `√(3p)` = the compositum's real part.

They **collapse to one field `Q(√-3)` at `n=6`** (`p=3`), and **fan into the biquadratic `Q(√-3,√-7)` at `n=14`**. LRC is proved exactly at the collapse. So the owner's question has a sharp answer: *the spine stays quadratic; its field matches the covering's only at `n=6`; and the failure to match at `n=14` — the split `Q(√-7)` vs `Q(√-3)` — is precisely the extra arithmetic a proof must cross.*

## Honest status & next

- **Verified:** spine `= {0,±i√p}` (`Q(√-p)`, degree 2) for `p≤31`; covering `= n/Φ₆(n)` with Φ₆-field `Q(√-3)`; match iff `p=3`; deep-well `D_n` with `√21 | √D_{14}`.
- **Framing:** "field split = difficulty" is a clean *necessary-condition* reading (proved case = single field), grounded but not a theorem that field-alignment ⇒ provable.
- **Next:** does the biquadratic `Q(√-3,√-7)` structure give the *bridge* — e.g. is the certifying ι-odd Lefschetz number (S20, `i√7`) forced to interact with the Eisenstein covering (`√-3`) through `√21`, and is that interaction the residual `OPEN-Q-108`? Compute the certification in the compositum. And for `p≡1 mod 4` (the easy Brouwer regime, S18) the spine has *real* eigenvalues (`χ(−1)=+1`, `g_p=√p`) in `Q(√p)` — check whether *that* matches the covering (still `Q(√-3)`), predicting the easy regime is where the spine is real and the covering imaginary (never matching, but Brouwer needs no match).

— Related: `the-singular-series-is-an-iota-equivariant-lefschetz-trace-…` (S20, `g_7=i√7`), HYP-3815 (opus, Paley skew `{0,±i√p}`=Frobenius), `twentyeight-…` (three pillars, `√-7`), HYP-3773/`the-flat-extension-moments-are-ramanujan-sums-…` (Φ₆, `−1/12`, Dedekind), `lrc-difficulty-by-n-…` (S18, by-n), HYP-3796 (deep-well `[0;n-1,n]` Herman rigidity), OPEN-Q-108. Script: `04-computation/spine_eigenvalues_vs_covering_cf_field_kps.py` (+ .out). Not a HYP reservation.
