# Why the covering-minimum is exactly 14/183 — the equioscillating ladder, and what 169 was marking

*klein-2026-07-18-S317. Written as a depth-first chase of every "169" in the repo; the thread closed into
a one-line explanation of a number the project has treated as a brute fact for months.*

## The number that kept appearing

`169 = 13²` shows up in the corpus in four unrelated-looking places:

1. **MISTAKE-104**: the `n=13` deep well `{1,…,11, 168}`, killer `168 = 13²−1`, witness `t = 14/169`,
   equioscillating (`v=1` and `v=168` both at distance exactly `14/169`); single-lift rigidity gap `1/169`.
2. **THM-1007**: the lacunary chain's tightest strict value `169/2352 = (1/12)(13/14)²`.
3. **THM-897**: the triple-beat discrepancy `|μ(W∩D₃) − (2/13)μ(W)| ≤ 22κ_W/(169·x₃)`.
4. **MISTAKES.md**: an "`a/169` grid stratum" left behind by a replaced route.

Two of these turn out to be the *same* structure seen twice, one is a genuine obstruction-measure, and one
is a coincidence of arithmetic. Separating them is the point.

## The equioscillating ladder

Both deep wells are instances of one family. Take the body `{1,…,b}`, a killer `k ≡ −1 (mod q)`, and the
witness `t = m/q`. Then `v=1` sits at distance `m/q` and `v=k` sits at `‖−m/q‖ = m/q` — an
**equioscillation**. The body stays safe iff `mj ≤ q−m` for `j ≤ b`, i.e. `q ≥ m(b+1)`. Taking the
smallest such `q` and `b = 12` (so the set has 13 speeds) gives `q = 13m+1`, `k = q−1 = 13m`:

```text
V_m = {1,…,12} ∪ {13m}          M(V_m) = m/(13m+1)      (exact, verified m = 1…16, 28, 42)
```

`M` is **increasing** in `m`, rising to `1/13` from below. And `V_m` is primitive for every `m`, while

```text
V_m is covering  ⟺  14 | 13m  ⟺  14 | m.
```

## The explanation

Read the ladder in order:

| `m` | killer | `M` | status |
|---|---|---|---|
| 1 | 13 | **1/14** | the tight AP `{1,…,13}` — the LRC(14) extremal. Non-covering. |
| … | … | increasing | all non-covering |
| 13 | **169 = 13²** | **13/170 = 0.0764706** | *below* the covering-min — but non-covering |
| **14** | **182** | **14/183 = 0.0765027** | **covering — THE COVERING-MINIMUM** (THM-724) |
| 28 | 364 | 28/365 | covering (tower) |

So the covering-minimum is not a coincidence and not an isolated well. It is **the first rung of a
monotone ladder at which the covering condition switches on**. The whole ladder below `m = 14` sits *lower*
— `m = 13` is lower by `1/(170·183)` — and is excluded for one reason only: `14 ∤ m`, so the set has no
multiple of 14 and LRC(≤13) dispatches it. **`14 | m` is precisely what lifts `m` from 13 to 14, and
`14/183` is what you get when it does.**

And that is what `169` was marking: it is the **last non-covering rung before the covering-minimum**, the
near-miss that would have been the answer if divisibility had not intervened. Its being `13²` is a numerical
accident (`13·13 = 169` and `14·12+1 = 169` collide), which is why it looked structural in the `n=13` deep
well `{1,…,11,168}` — there `q = 14·12+1 = 169` by the *same* ladder at `b = 11`, and it merely happens to
be a perfect square.

## Cross-check: this re-derives THM-1014 from the other side

THM-1014 proved, from the witness-predicate direction, that `V = d·{1,…,12} ∪ {k}` is sub-`1/13` exactly
when `13d | k`, that primitivity forces `d = 1`, covering forces `182 | k`, and hence every exception is
non-compact (`ρ = k/12 ≥ 182/12 > 13`). The ladder says the same thing from the parametric side: the
exceptions **are** `V_m` for `14 | m`, i.e. the tower `{1,…,12, 182m′}`, and their `ρ = 13m/12` is forced
past 13 by `m ≥ 14`. Two independent derivations, same tower.

## The other 169 — a real obstruction, and the lesson that closed a stratum

`169/2352 = (1/12)(13/14)²` is a different animal. The balance-ladder route (THM-724 → THM-1007) pays a
factor `13/14` **per killer**, so `j` killers cost `(13/14)^j`; at `j = 2` that is `169/196`. This `169` is
a *per-step multiplicative loss*, and it is exactly why the **clustered** stratum stayed open: clustering
means the running maximum barely shrinks between steps, so the product degrades toward `1/2` and the
telescoping dies.

Naming the obstruction that way said what to do. An additive certificate has no per-step cost. The
interval-survival tail (THM-1004) lets `r` killers enter only through `Σ 1/k_i` — additively and
symmetrically, with **spacing invisible** — and closes the clustered stratum for large killers
(THM-1015). The moral is worth keeping:

> When an obstruction is multiplicative *per step*, do not sharpen the step. Replace the ladder with a
> certificate in which the steps enter additively.

## What this suggests next

- The ladder is the `b = 12` case of `V_{m,b} = {1,…,b} ∪ {(b+1)m}` with `M = m/((b+1)m+1)`. The same
  reading should locate the covering-minimum at every `n`: it is the first rung where `(n) | m`. Worth
  checking whether the known small-`n` covering minima all sit at that rung.
- If the covering-minimum is always "the first divisible rung", then the *margin* between the covering-min
  and the true infimum over the ladder is `1/((13m+1)(13m+14))`-sized — an explicit, shrinking quantity.
  That is a cleaner handle on the covering-min's tightness than the census route.
- `169`'s two roles should not be conflated in future write-ups: the deep-well `169` is an accident of
  `14·12+1`, the balance-loss `169` is `13²` and is structural.

*Files: `04-computation/lrc14_equioscillating_family_klein_S317.py` (+ .out); THM-1014, THM-1015.*
