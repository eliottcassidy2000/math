---
source: oracle-2026-06-03-S578o
status: progress (3-term vs 4-term in LRC; hardness = 3-term folds; 4-term energy is the translation-invariant shadow; how 3 encodes 4; corrected a BSG guess; implications for Lemma A)
tags:
  - lonely-runner
  - three-term
  - four-term
  - additive-energy
  - fold
  - summand-graph
  - deformation
  - lemma-A
---

# 3-Term Folds Carry the Hardness; 4-Term Energy Is the Translation-Invariant Shadow

Working the user's 3-vs-4 decomposition (Lemma A: circuit-free ⟹ `G ≥ δ`; Lemma B: a
3-term relation `v_c=v_a+v_b` is a fold). Going back and forth between the two, the **shifted
AP** is the particularity that settles where the hardness lives — and where it *hides*.

## The two layers, and their energy weights

In the resonance energy `E(v)=Σ_{Σ mᵢvᵢ=0} ∏|ĝ(mᵢ)|` (S550):
- a **3-term fold** is the relation `(1,1,−1)`, weight `(1−2/n)^{k−3}|ĝ(1)|³` — the heavy layer;
- a **4-term** is `(1,1,−1,−1)`, weight `(1−2/n)^{k−4}|ĝ(1)|⁴` — lighter by `(1−2/n)/|ĝ(1)| ≈ 6`.

So one fold is worth ~6 four-terms (n=14). Folds dominate the energy; 4-terms are the next,
softer layer.

## The smoking gun: the shifted AP (hardness = 3-term, *not* 4-term)

`lrc_three_four_term_energy_encoding_s578.py`, the table that says everything
(`{N..N+12}`, n=14, `δ=1/14`):

| N | set | #3-term | #4-term | G | G/δ |
|---|-----|---------|---------|---|-----|
| 1 | {1..13} | 36 | **125** | 1/14 | 1.00 (TIGHT) |
| 3 | {3..15} | 25 | **125** | 0.167 | 2.33 |
| 7 | {7..19} | 9 | **125** | 0.269 | 3.77 |
| 15 | {15..27} | **0** | **125** | 0.357 | 5.00 |
| 100 | {100..112} | **0** | **125** | 0.472 | 6.60 |

**Additive (4-term) energy is translation-invariant** (`a+b=c+d ⇔ (a+s)+(b+s)=(c+s)+(d+s)`),
so it stays pinned at **125** the whole way. **3-term folds are not** (`v_c=v_a+v_b` needs the
sum to *land on a vertex*, destroyed by shifting up), so they fall `36 → 0`. And **G rises
monotonically from `δ` (tight) to `~6δ` (very safe)**. The hardness tracks the 3-term folds
exactly while the 4-term energy is *blind* to it.

> **Hardness is carried by 3-term folds, not 4-term energy.** This is the de-risking of
> Lemma A made mechanistic: 4-term-rich configs are safe because 4-energy is the wrong
> invariant — the shifted AP has *maximal* 4-energy and is *very* safe.

## Information hidden in the deformation

The translation `S → S+s·𝟙` is a **deformation that hides the hardness**: it preserves the
entire 4-term/additive-energy structure (the AP's 125 quadruples are untouched) while
silently destroying every 3-term fold. **An observer looking only through the
energy/4-term/`L⁴` lens sees nothing change** as the AP goes from tight to safe. The hardness
information lives precisely in the part the 4-term lens cannot see — the *absolute-position*
3-term folds.

## How 3 encodes 4 (both directions)

- **4 = depth-2 fold (3 in `S∪(S+S)`).** A 4-term `v_a+v_b=v_c+v_d=C` is a summand-graph
  node `C` of in-degree ≥2; *add* `C` as a vertex and it splits into two 3-term folds. So
  **4-term(S) = 3-term folds of `S∪(S+S)`** — a 4-term is a fold whose apex is another
  pair-sum instead of a vertex. (AP and shifted-AP have the *same* 19 collision nodes / 74
  folds-into-collisions, though the shifted AP has zero folds *in S itself*.)
- **3 generates 4.** Two folds sharing a summand — `v_a+v_b=v_c`, `v_a+v_d=v_e` — force
  `v_c−v_b = v_e−v_d`, i.e. the 4-term `v_c+v_d=v_e+v_b`. Verified: **94 of the AP's 125
  four-terms** arise from shared-summand fold pairs (the rest from longer fold chains). So
  the fold structure *generates* the additive energy.
- **4 = the translation-invariant shadow of 3.** Shifting kills the folds but keeps all 125
  four-terms — the 4-term relations are exactly the part of the fold structure that depends
  only on *differences*, hence survives translation. 3-term = folds + absolute anchoring;
  4-term = the difference-only residue.

## Correction (a guess refuted by the particularity)

I first conjectured a **Balog–Szemerédi–Gowers ceiling**: 3-term-freeness should *cap* the
4-term energy (high energy ⟹ near-AP ⟹ folds). **False** — the shifted AP is 3-free with
*maximal* 4-energy. The flaw: BSG/Freiman gives "structurally an AP," but an AP can be
3-term-free (shift it up); the 3-term property is translation-sensitive, additive energy is
not. **3-freeness does not bound 4-energy.** Using the particularity (4-energy's translation
invariance) as noise refuted the clean-but-wrong idea.

## Implication for Lemma A (the open part)

Lemma A ("3-term-free ⟹ `G ≥ δ`") **cannot** be proved by bounding additive energy / an
`L⁴`-discrepancy bound — the shifted AP shows 4-energy can be maximal in a 3-free safe set.
The discrepancy input must come from the **absolute** 3-term-free structure, not energy. And
even "equidistribution" is subtle: the shifted AP's positions `{(N+j)t}` are themselves an
AP-in-position (poorly spread as a set) yet `G` is large — so safety is about the *observer
landing in a large gap*, which 3-term folds are what destroy (a fold pins a runner near the
observer at the critical time). **The right Lemma-A statement is gap-structural, not
equidistribution- or energy-based.**

## Verdict / next
- **Hardness = 3-term folds; 4-term energy is the translation-invariant shadow** (red herring
  for hardness). The translation deformation hides the hardness from the energy lens.
- **3 encodes 4** three ways: 4 = depth-2 fold in `S∪(S+S)`; folds *generate* 4-terms
  (shared-summand); 4 = difference-only residue of folds.
- Concrete next: (1) a **gap-structural** Lemma A — prove 3-term-free ⟹ the observer can sit
  in a `≥δ` gap (not via energy); (2) quantify the `G ↔ #3-term` monotonicity (the shifted-AP
  curve) — is `G` a decreasing function of a *weighted* fold count?; (3) Lemma B: which folds
  actually pin `G=δ` (only the AP's *low* folds `1+2=3,…` do — `{2..14}` has 30 folds but
  `G=1.75δ`), tying to the summand-graph node arithmetic (S560o) and the apex (S577o).

## Artifacts
```
04-computation/lrc_three_four_term_energy_encoding_s578.py
05-knowledge/results/lrc_three_four_term_energy_encoding_s578.out
```
Related: S550 (resonance energy / sufficient condition), S559o/S560o (summand graph: pair-sum
node = pinch denominator; AP = Freiman joint-extremum), S577o (the apex / mod-7 tie-wall),
HYP-2075 (pinch witnesses), `additive_comb_bridge.py` (degree-4 Walsh = additive quadruples).
