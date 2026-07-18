# The inverse theorem is a function-field Freiman 3k-4 — and the far element is the second dimension

**death-star-2026-07-18-S56.** Working the LRC(14) wall (boxeph THM-1017: `M<1/13` covering ⟹ the 12
non-max speeds are a dilated AP) through the Freiman 3k-4 lens the owner suggested — specifically the
**function-field / coset-progression** version. The payoff: a sharp additive reformulation, a precise
diagnosis of why the *ordinary* (ℤ) 3k-4 cannot apply, and the identification of the right tool. **The
wall is not closed** — but it is now a bounded-dimension statement, one step past boxeph's one-gap lemma.

---

## 1. The residue picture (from boxeph THM-1017 / difference-closure lemma)

At the maximizer `t = a/q`, `M = val/q`, the 13 residues `r_i = v_i·a mod q` lie in the band
`[val, q−val]` (length `< 12·val`, since `M<1/13 ⟺ q < 14·val`). boxeph's lemma: some pair is within
`val` (one aligned gap). The deep well `{1..12,182}` (`q=183, val=14`) makes it concrete:
```
residues = {14, 28, …, 168} ∪ {169} = 14·{1,…,12}  ∪  {169},   169 ≡ 1 (mod 14).
```

## 2. The additive signature of `M<1/13` (verified, sharp)

Reducing the residues **mod val** across the whole deep-well tower `{1..12,182m}` and dilates:

| family | `M` | `val` | # residue classes mod `val` | core `= val·{1..12}`? |
|---|---|---|---|---|
| deep well | `14/183` | 14 | **2** (sizes 12,1) | ✓ |
| ladder m=2,3,7 | `<1/13` | 28,42,98 | **2** | ✓ |
| `5·core+killer` | `14/183` | 14 | **2** | ✓ |
| near-miss `{1..11,13,84}` (`M>1/13`) | `0.0787` | 7 | **6** (scattered) | — |

**`M<1/13` ⟹ the residues occupy exactly 2 classes mod `val`, split 12+1**; `M>1/13` scatters them.
So the additive signature of dropping below `1/13` is: the residue set is a **coset progression of
dimension 2** — twelve residues in the coset `0 + val·ℤ`, one in a second coset.

## 3. The band forces the AP (rigorous)

The twelve residues `≡ 0 (mod val)` are `val·c_i` with `c_i` distinct and `val ≤ val·c_i ≤ q−val`, so
`c_i ∈ [1, q/val − 1) ⊂ [1,13)` because `q < 14·val`. Twelve distinct integers in `{1,…,12}` are all
of them: **the class-0 residues are exactly `val·{1,…,12}` — a dilated AP.** So the inverse theorem is
*equivalent* to the additive statement:

> **`M(V) < 1/13` ⟹ the residue set `a·V mod q` has additive dimension ≤ 2** (lies in two cosets of the
> subgroup `val·ℤ`, twelve in one).

This is strictly sharper than boxeph's difference-closure lemma (which delivers *one* aligned gap, i.e.
dimension `< 13`, not `≤ 2`). The whole gap between the two is the Freiman step.

## 4. Why ordinary Freiman 3k-4 fails — and why the function-field version is the tool

The temptation is to bound the doubling of the residue set `R` and quote 3k-4 (`|R−R| ≤ 3k−4 ⟹ AP`).
**It fails on the nose.** For the deep well, `k=13` and
```
|R−R| = 47  >  3k−4 = 35,
```
because the far residue `169` contributes a whole second progression `1 + val·ℤ` to the difference set.
The *raw doubling is too large* — the deep well is not an AP and never will be (it's genuinely
dimension 2). The ordinary ℤ 3k-4, keyed on `|A+A|`, cannot see the structure.

The **function-field / coset-progression** Freiman theorem is exactly the tool that can: it characterizes
sets by **additive dimension relative to a subgroup**, not by raw doubling. In `F_q^n` (or `F_q[t]`),
"small doubling" is replaced by "few cosets of a subspace," and 3k-4 becomes "dimension 1 unless a
genuine second direction is present." Here the subgroup is `val·ℤ` (or `val·ℤ / qℤ`), the first
direction is the AP core, and the far element is the **second dimension**. The residue set is a
Green–Ruzsa coset progression of dimension exactly 2, and the inverse theorem is the assertion that
`M<1/13` forbids dimension 3+. That is a clean function-field-3k-4-shaped statement, and it is the one
worth proving — not a doubling bound.

## 5. What this buys, and the residual

- **Reframes the wall** from "core is a dilated AP" to "residue set has additive dimension ≤ 2" — a
  bounded-dimension inverse statement, the natural home of function-field Freiman.
- **Diagnoses the obstruction to the naive route:** raw doubling is `> 3k−4` (the far element), so any
  attempt via ℤ 3k-4 is doomed; the dimension/coset framing is mandatory.
- **Sharp reduction:** dimension ≤ 2 + the band `q<14·val` ⟹ AP core (§3, rigorous). So the entire open
  content is: **`M<1/13` ⟹ additive dimension ≤ 2.**

**Residual (the wall, unchanged in hardness):** prove `M<1/13 ⟹` the residues lie in ≤ 2 cosets of
`val·ℤ`. boxeph's pigeonhole gives one aligned pair; the function-field 3k-4 target is to promote that
to the full 2-coset structure. The likely mechanism: `M<1/13` is a near-extremal (high additive-energy)
condition on `a·V`, and the covering constraint forces the second direction to be the `lcm(13,14)=182`
killer — so the second coset is *arithmetically pinned* (`≡ 1·val`-ish, the `169 = 13²` alignment),
which is the extra rigidity a coset-progression argument can exploit that a raw-doubling one cannot.

→ THM-1017, THM-1028, boxeph difference-closure lemma (169 reflection), HYP-7362; owner cue: Freiman
3k-4 for function fields.
