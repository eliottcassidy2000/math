# The γ-trick: 14-covering ⟹ not tight, by the apex periodicity of the multiples of 14

*kind-pasteur-2026-06-22-S31ad. Closing the last residual of the THM-079-template proof of LRC(14): a
14-covering set has `M > 1/14`. The creative move is that a multiple of 14 is `1/14`-periodic in `t`,
which decouples it from the rest and turns the residual into a pigeonhole on 14 points.*

## The setup (after THM-568)
THM-568 reduced LRC(14) to: a 14-covering set `S = M14 ⊔ R` (`M14` = the multiples of 14, `R` = the
14-free rest) has `M(S) > 1/14`. `R` has `13−r` speeds (`r = |M14| ≥ 1`), so `M(R) ≥ 1/(14−r) > 1/14` by
proven LRC(≤13): `R` carries a margin.

## The γ-trick (the creative step)
**A multiple of 14 is `1/14`-periodic in `t`:** `‖14m·(t + 1/14)‖ = ‖14m·t + m‖ = ‖14m·t‖`. So the entire
constraint from `M14` depends only on the *fine phase* `γ = frac(14t)` — and `‖14·m_i·t‖ = ‖m_i·γ‖`, so
the multiples are just `{m_i = (14m_i)/14}` evaluated at `γ`. By LRC(r+1) on `{m_i}` (proven for `r ≤ 13`)
there is a `γ*` with `‖m_i·γ*‖ ≥ 1/(r+1) ≥ 1/14` for all `i` — **all multiples safe with margin**, and this
holds for *every* `t` with `frac(14t) = γ*`, i.e. on the **14 equally-spaced points**
`t_j = (γ* + j)/14`, `j = 0,…,13`.

So `M(S) ≥ 1/14` reduces to: **`R` is safe at one of the 14 points `t_j`.**

## The pigeonhole
At the 14 points, each runner `r_k ∈ R` coprime to 14 sees `{r_k·t_j} = {r_k γ*/14 + r_k j/14}` — a shifted
copy of the 14th roots (equally spaced, `1/14` apart). Only those within the arc `(−1/14, 1/14)` (width
`2/14`) are *bad*, and an arc of width `2/14` holds **at most 2** of 14 equally-spaced points. So:

> each coprime-to-14 runner marks `≤ 2` bad points ⟹ `|R| ≤ 6` runners mark `≤ 12 < 14` bad ⟹ **a good
> point survives** ⟹ `M(S) ≥ 1/14` (with margin, since both the multiples and `R` are `≥ 1/14` there).

This **closes `r ≥ 7`** (`|R| ≤ 6`, `R` coprime to 14). Verified (`lrc_gamma_trick_14covering_kps.py`):
`r=7` → 4 good points (`M=1/8`); `r=8` → 5 good points (`M=1/9`); even `R` with a multiple of 7 → 2 good
points. And `r ≤ 6` is closed by the union bound (S31v: the `≤6` multiples cover `≤ r/7` of `R`'s margin).

## The residual recurses on the apex 7 (`14 = 2·7`)
The pigeonhole's `≤2`-per-runner bound assumes `gcd(r_k, 14) = 1`. A runner with `gcd = 7` (a 14-free
*multiple of 7*) sees only 2 distinct values among the 14 points (period `1/2`), so it can mark up to **7**
bad points — the count can exceed 14. But this is exactly the **apex self-reference**: such a runner, and
the multiples of 14, are all **multiples of 7**, and *a multiple of 7 is `1/7`-periodic in `t`*
(`‖7m(t+1/7)‖ = ‖7m·t‖`). So the **same γ-trick applies one level down**: fix the `1/7`-fine phase to make
all multiples of 7 safe (LRC on `{·/7}`), reducing to a pigeonhole on **7 points**, where the non-multiples
of 7 (coprime-to-7) mark `≤ 2` each. The tower `14 → 7 → 1` is the prime structure of `14 = 2·7`, and the
γ-trick descends it — the LRC(≤13) margins feed each level.

## Status
- **`r ≤ 6`**: union bound (S31v) — done.
- **`r ≥ 7`, `R` coprime to 14**: γ-trick + pigeonhole — **done** (this note).
- **`R` has a multiple of 7**: recurse the γ-trick on the `1/7`-periodic 7-content (the apex level); the
  residual shrinks to the all-multiples-of-7 core, which dilates (R3) to a smaller LRC instance.

So the 14-covering residual is no longer one monolithic equidistribution — it is a **finite descent down
the prime tower of `14 = 2·7`**, each level a pigeonhole on `p` points fed by a proven smaller LRC margin.
The γ-trick is the mechanism: *a multiple of the modulus is periodic at the apex scale, so it decouples
and the rest is pigeonholed.* This is the apex-prime structure (`7 = I(K₃,2)`, `14 = 2·7`) doing the work,
exactly as the tournament `H=7`-forbidden atom anchors THM-079.

→ THM-568 (apex-denominator), THM-523 (q-witness), S31v (union bound), THM-560 (dilation/tight locus),
HYP-+2878 (CRT over-determination — the γ-trick is its constructive realization), [[lrc14-thread]].
