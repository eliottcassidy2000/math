---
id: HYP-2125
status: EXPLORATION + RESULTS (verified + core formalized) — the recursive fractal of translated
  APs: Riesz-product factorization, multiplicative additive energy, self-similar 3-term count, and
  a SCALE-INVARIANT loneliness gap (the digit gap is a fixed point). Sum-free preservation
  formalized. The lacunary scale-invariance proof is OPEN.
source: claudebox-2026-06-03-S587
related:
  - HYP-2120
  - HYP-2064
---

# HYP-2125: the recursive fractal of translated APs has a scale-invariant gap

Recursing S585's AP-translation flip (HYP-2120): use a sum-free **translated-AP digit set** as the
digits of a base-`b` construction with `b > 2·max(digit)` (no carries),

```
S_d = { Σ_{i<d} a_i · b^i : a_i ∈ digits }   — a self-similar / Cantor-like set, |digits|^d speeds.
```

No carries ⇒ the additive structure factorizes across scales, two lenses:
- **rep theory:** the character polynomial is a **Riesz product** `p_{S_d}(x) = Π_{i<d} D(x^{b^i})`,
  `D(y)=Σ_{a∈digits} y^a` — a lacunary spectrum.
- **Haskell:** `S_d = map (base-b eval) (replicateM d digits)` — a `d`-fold fold.

## Verified factorization laws

- **Additive energy is multiplicative:** `E(S_d) = E(digits)^d` — exact (19, 361, 6859, 130321 for
  digits `{4,5,6}` or `{1,2,3}`, `E(digits)=19`). Energy is a *tensor power* across scales.
- **3-term count is self-similar:** sum-free digits ⇒ `N₃(S_d) = 0` at **every** scale (fractal
  circuit-free). For digits `{1,2,3}`: `N₃(S_d) = (3^d+1)/2` (2,5,14,41,122; recurrence
  `N₃(d)=3N₃(d−1)−1`) — a clean self-similar growth law.
- **Fractal difference set:** `|S_d − S_d| = |digits−digits|^d`, so `S−S` is itself self-similar of
  fractal dimension `log|D−D| / log|digits| = log 5 / log 3 ≈ 1.465`.

## The headline: a scale-invariant loneliness gap (verified)

`G(S_d)` **stabilizes to the digit gap** `G(digits)`, a fixed point of the recursion:
- digits `{4,5,6}`: `G(S_d) = 0.400, 0.399, 0.3965` (`= G(digits)` to 0.9%);
- digits `{1,2,3}`: `G(S_d) = 0.250, 0.250, 0.2494`.

Since `δ = 1/(|S_d|+1) → 0` as the set grows exponentially, the **margin `G−δ → G(digits) > 0`** —
fractal sets are *increasingly* safe with a fixed gap. (`{4,5,6}` margin: +0.15 → +0.30 → +0.36;
`{1,2,3}`: only the `d=1` base `{1,2,3}` is tight, `d≥2` are safe.)

**Mechanism (Riesz/lacunary):** at the digit witness `t*` (where `G(digits)` is achieved), the
finer factors `D(e^{2πi b^i t*})` have rapidly varying (lacunary) arguments, so the fine-scale
runners land away from 0 too — the coarse witness survives with an `O(small)` correction. The
lacunarity makes the scales asymptotically independent (Lemma A at each scale), pinning
`G(S_d) ≈ G(digits)`.

## Two fractal families (the regime picture)

- **Upper-half digits** `{h+1,…,2h}` → recursively **sum-free** → Lemma A: safe at all scales,
  exponentially large circuit-free sets with a fixed positive gap.
- **Lower digits** `{1,…,k}` → recursively **3-term-rich** (`N₃ ~ 3^d/2`) → Lemma B structure, but
  still eventually safe because `δ` shrinks faster than `G` drops. The tight base `{1..k}` does not
  propagate its tightness up the fractal.

## Formalized (math-lean `8443da1`, sorry-free)

`Math/LonelyRunner/FractalSumFree.lean` — `combine_threeTermFree`: combining a sum-free digit set
`A` (with `2a < b`, no carry) with any high set `B` at base `b` preserves 3-term-freeness (the
lowest digit alone blocks every relation `c+c'=c''` via the mod-`b` defect `∈(−b,b)`). Iterating is
the recursive fractal; this is its inductive heart, machine-checked.

## Open / next

- **Prove the scale-invariant gap** `lim_d G(S_d) = G(digits)` (or `|G(S_d) − G(digits)|` small) via
  the Riesz product / lacunary independence — a clean target, and a way to manufacture arbitrarily
  large LRC-easy sets with a uniform gap.
- **Energy-multiplicativity** `E(S_d)=E(digits)^d` is a clean next Lean lemma (the no-carry quadruple
  bijection, like `addEnergy_map_addRight`).
- The frontier of the conjecture is the **boundary** between fractal/lacunary-safe and
  dense-3-term-closed-hard sets; the tight `{1..k}` is maximally non-fractal.

**Artifacts:** formal — `claude-monad/math-lean` `Math/LonelyRunner/FractalSumFree.lean`
(`8443da1`); informal — `04-computation/lrc_fractal_translated_ap_s587.py` (+`.out`),
`07-reflections/lrc-fractal-translated-ap-s587.md`. Builds on HYP-2120/S585 (the AP-translation
flip, fusion, energy symmetry), HYP-2064 (sieve / Lemma A).
