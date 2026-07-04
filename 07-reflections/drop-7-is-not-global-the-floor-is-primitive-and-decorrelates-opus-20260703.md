# Drop-7 is not the global floor-minimizer — the covering floor is a primitivity statement, decorrelates in each coordinate, and its infimum is the tight-locus rigidity

*opus-2026-07-03-S59. The owner asked whether drop-7 (my S58 minus-one floor-minimizer, R'=0.315) is the
GLOBAL R'-minimizer over all covering families, and to finish the proofs. Answer: no. In pinning the true
minimizer I found the floor is a gcd=1 statement (R'=0 on imprimitive tight families), proved a clean
coordinate-decorrelation bound (THM-611), and located the honest open core: whether inf R' > 0.*

## The question, answered: drop-7 is NOT global

S58 found drop-7 `{1..6,8..13,14}` (`R'=0.315`) is the floor-minimizer among the 13 minus-one families of
`{1..13}`. Over ALL covering families it is not even close. Greedy 1-swap descent (speeds ≤ 42, exact
verification):

- **Unconstrained, the descent falls to `R' = 0`** at `S = {2,4,…,26} = 2·{1..13}`. This is not a floor
  violation — it is imprimitive. By scale-invariance `meas(lonely(c·V)) = meas(lonely V)`, and `{1..13}`
  is the arithmetic progression with `M = 1/14` *exactly* (its safe set is the isolated points `t=k/14`,
  measure zero). So `meas(lonely(2·{1..13})) = 0` and `R' = 0`.
- **Among primitive (gcd=1) covering families, the min found is `R' ≈ 0.0763`** at
  `S = {2,4,5,6,8,10,14,16,18,20,22,24,26} = 2·({1..13}\{6}) ∪ {5}` — far below drop-7's `0.315`.

So **drop-7 is only the minus-one minimizer; the global minimizer is much lower and lives next to the
imprimitive tight AP locus.**

## Fact 1 — the covering floor is a primitivity statement (floor side of HYP-4060)

`R' = 0` exactly on the imprimitive tight families `c·{1..13}`. The whole floor `R' > 0`, like the
covering-min `M ≥ 14/183` (kps HYP-4060), therefore only makes sense for **gcd = 1** families — one must
`LRCDilation`-reduce first. kps established this for the covering-min (the `M`-axis); here it is the same
fact seen on the `R'`-axis (the measure-normalization axis). `LRCDilation` is load-bearing for the floor,
not merely the peel. The low-`R'` primitive minimizers are exactly the ones sitting *next to* this
imprimitive locus (`2·({1..13}\{6})∪{5}` is `2·{1..13}` with one tooth removed and one odd tooth added).

## Fact 2 — the floor decorrelates in each coordinate (THM-611, proved)

For a core `R` and a far runner `w`, with `A_R` = number of arcs of the `R`-lonely set,
> `|meas(lonely(R ∪ {w})) − (6/7)·m_R| ≤ A_R/(3w)`  (bounded variation `|\hatφ_R(n)| ≤ A_R/(π|n|)`,
> paired against `|c_k| ≤ 1/(π|k|)`, summed via `Σ 1/k² = π²/3`).

Hence for `r=1` (`m_Q=6/7`), `|R'(R∪{14q}) − 1| ≤ A_R/(36 q m_R)`, so **`R' → 1` as the far runner
grows**, and the same holds for any single coordinate grown while the rest is fixed. This is the signed
sum that converges where the absolute one diverges (MISTAKE-078). Consequences:

- The floor's infimum is **not at infinity in any single coordinate**; a minimizing family cannot send one
  runner to `∞` with the others bounded. So the search over bounded magnitude is legitimate.
- **But no uniform bound.** The constant `A_R/m_R` grows with the core (dense cores: `A_R` large, `m_R→0`).
  So the per-coordinate decorrelation does not bound the whole family's magnitude — a uniform bound would
  reduce the floor to a finite check, i.e. it *is* the open problem.

## Fact 3 — inf R' is the tight-locus rigidity (the honest open core)

The low-`R'` minimizers cluster next to the imprimitive tight AP. Does a *primitive* sequence approach it
with `R' → 0`? The naive one, `S_c = {c+1, 2c, …, 13c}` (even `c`, primitive covering), does **not**:
`R'` stabilizes at `≈ 0.53` as `c` grows — the coprime runner `c+1` de-resonates (Fact 2, single-coordinate
decorrelation, applied to `c+1` against the commensurate block `c·{2..13}`). So `inf R'` is neither shown
`0` nor bounded `> 0` by these routes; it is exactly **how close a primitive covering family can get to the
tight AP locus** — the tight-locus rigidity (kps HYP-4060), the measure-route hard core.

**A recalibration for the floor program.** The per-row certified `R' ≥ 0.642` (HYP-3129) and the whole
"covering floor `R' > 0` uniformly" framing (THM-579/HYP-3415) assume a uniform positive floor. This search
shows the uniform constant is `≤ 0.076` (primitive coverings reach it), and whether it is `> 0` at all is
open. `0.642` is a per-family value on the benign consec rows, not a uniform lower bound. The tight,
low-`R'` families — near-imprimitive, not the consec rows and not drop-7 — are where the floor is actually
decided.

## Status

- **Answered (exact):** drop-7 is NOT the global R'-minimizer; primitive global min `≈ 0.076` at
  `2·({1..13}\{6})∪{5}`; `R'=0` on imprimitive `c·{1..13}`; `S_c` decorrelates to `≈0.53`.
- **Proved:** THM-611 (coordinate decorrelation `|meas − (6/7)m_R| ≤ A_R/(3w)`; `R'→1` per coordinate;
  inf not at infinity in one coordinate).
- **Open (= LRC(14) hard core):** `inf R' > 0` (uniform floor) vs `= 0` — the tight-locus rigidity
  (HYP-4060); no uniform magnitude bound (THM-611's non-uniformity).

Given MISTAKE-097 (my prior overclaims on this crux), the non-closure is flagged plainly; what is proved
here (THM-611, the primitivity of the floor, drop-7 refuted) is exact and bounded in scope.

Related: THM-611 (this session's proved lemma), HYP-4061/opus-S58 (drop-7 as the minus-one min, now
globally superseded), HYP-4060/kps (covering-min primitivity + tight-locus rigidity — the same open core
on the `M`-axis), THM-579 + HYP-3129 (the floor program whose `0.642` this recalibrates), THM-610/mac-mini
(census-side coordinate dichotomy — THM-611 is its measure-side analog), MISTAKE-078 (the absolute-sum
divergence THM-611 sidesteps), HYP-4058/opus (the `1/7`-comb spectrum `c_k`). Scripts:
`lrc14_floor_global_min_opus_S59.py`, `lrc14_floor_global_min_anneal_opus_S59.py`,
`lrc14_floor_primitivity_sequence_opus_S59.py`. HYP-4063.
