# HYP-2325 — Friendliness = "never having been lonely yet": the first lonely time τ, the survival curve, and the 1/n floor

**Session:** S647
**Status:** CONFIRMED (1/n floor + monotone "yet" formalized; survival behaviour verified + charted)
**Provenance forward:** math-lean `Math/LonelyRunner/Friendliness.lean` (sorry-free)
**Prompt:** make friendliness "the property of never having been lonely yet" and run a long session.

---

## 0. The definition and what it makes LRC into

Instantaneous friendliness (S647 first chart) was *how many runners you crowd now*. The user's sharper
definition makes it a **survival property**: you are **friendly** at time `t` iff you have **never been
lonely** on `[0, t]`. Since "has been lonely by `t`" only grows (the **"yet"**), friendliness is
**one-way / first-passage**: it holds on `[0, τ)` and then ends forever, where

```
  τ  :=  the FIRST lonely time  =  inf{ t > 0 : lonely t }
       =  the first GAP in the danger-arc covering  ⋃ᵢ {t : ‖vᵢ t‖ < δ}   (the Vitali view, HYP-2200)
```

So this friendliness recasts LRC as a **first-passage / covering** problem:
- `τ` = the right endpoint of the first covered block from `0` = where the danger arcs first leave a gap.
- **LRC ⟺ `τ < 1`** (a gap exists within the lap ⟺ the arcs don't cover `[0,1)` ⟺ everyone eventually
  gets a lonely moment). The survival curve reaching `0` is the conjecture.

---

## 1. Formalized (math-lean, sorry-free): `Math/LonelyRunner/Friendliness.lean`

- **`everLonely_monotone`** — `everLonely t := ∃ s ∈ [0,t], lonely s` is **monotone** in `t`; so
  friendliness (its antitone negation) is a genuine first-passage/survival property — the "yet".
- **`friendly_until_inv_n`** — the **1/n floor**: with a runner of speed `1` present, for every
  `t ∈ [0, 1/n)` the clock distance `dZ(1·t) = dZ t ≤ t < 1/n`, so that runner is inside your gap and
  you are not lonely. Hence **`τ ≥ 1/n`**: friendliness is *guaranteed* for the first `1/n` of the lap.
  (General floor: `τ ≥ 1/(n·v_min)`; with the unit speed it is `1/n`.) Uses the `dZ`/`dZ_le_dist_int`
  covering-depth machinery (HYP-2195).

---

## 2. The survival curve and the first lonely time (verified + charted)

`friendliness_survival_s647.py` computes `τ` per config and the **survival function**
`S(t) = P(never lonely yet by t)` over random `n`-runner configs (charted, PNG+SVG):

- Each curve is **flat at `1` until `t = 1/n`** (the formalized floor), then decays to `0`.
- **Median first-lonely time shrinks with `n`:** `τ_med ≈ 0.24, 0.13, 0.09` for `n = 5, 8, 12`
  (`≈ 1/n`-ish) — *more runners ⟹ loneliness comes sooner*. The survival drops faster and earlier for
  larger `n`.
- Every curve reaches `0` within the lap — the empirical face of LRC (every config gets a lonely moment).

**Two regimes of τ** (verified, panel A):
- **Generic config** (e.g. the `n=14` wall `{1..11,13,14}`): positive lonely measure `p₀ = 0.0122`, with
  a definite first lonely time `τ = 0.0824` (friendly for `8.2%` of the lap; `τ ≥ 1/14 = 0.071` ✓).
- **Tight extremal `{1,…,n−1}`**: lonely set is **measure zero** (lonely only at the isolated tight
  instant `t = 1/n`, where all `‖kt‖ = 1/n` with equality). So you are "never lonely yet" *almost
  everywhere the whole lap* — friendliness survives **maximally**, `τ` is a single point. The extremal
  LRC case is the friendliest one: it touches loneliness on a null set.

---

## 3. Ties & synthesis
- **First-passage = covering = Vitali** (HYP-2200/S618): `τ` is the first gap in the arc cover;
  friendliness-survives ⟺ no gap yet ⟺ the partial cover `[0,t]` is still full. LRC = a gap exists.
- **The 1/n floor is the slowest runner's arc**: near `0` every runner's danger arc contains `0`, and the
  *widest* is the slowest's, half-width `1/(n v_min)`; so the first block is `[0, 1/(n v_min))`. The unit
  speed makes this `1/n` — the gap itself. The floor is the conjecture's own constant.
- **The tight case is the friendliest** (measure-zero loneliness) — the extremal config for LRC
  maximizes friendliness-survival; the resonance/collapse family (`p₀ = 0`, S617) is exactly "friendly
  a.e. the whole lap."
- **Survival = `1 − everLonely`** is the monotone complement of the lonely-by-`t` process; the lonely
  measure `p₀` (the whole arc's master object) is the *instantaneous* bottom bar, while `τ` is its
  *first-passage* shadow.

## 4. New threads / handoffs
- **τ as a config invariant:** is `τ(V) ≥ 1/n` always (the LRC bound is a first-passage floor), and is
  `τ ≤` something clean (an upper bound = the conjecture)? Formalize `τ < 1` ⟺ LRC.
- **Median `τ ≈ c/n`:** pin the constant `c` and the survival-curve shape (does `S(t)` have a scaling
  limit `S(nt) → ?`).
- **The tight/extremal = measure-zero-lonely = maximal friendliness:** connect to the collapse/additive-
  chain family (S617) — friendliness-survival is maximized exactly on the resonance-rich configs.
- Formalize the covering statement `τ = inf (⋃ᵢ Aᵢ)ᶜ` and `LRC ⟺ τ < 1`.
