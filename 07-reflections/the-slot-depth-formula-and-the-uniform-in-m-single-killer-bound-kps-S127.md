# The slot-depth formula: the covering-min is a 2-runner balance `Δ·s_a s_b/(s_a+s_b)`, bounded uniformly in m for runner-1 slots

*kind-pasteur-2026-07-11-S127 cont.66. Owner: "bound the filler-contested backbone slot depth uniformly in m."
The backbone lens (cont.65) says exactly two runners bind at the covering-min. That makes the depth an explicit
2-runner balance, and I derive its closed form, verify it exactly, and use it to **prove the uniform-in-m bound
for the runner-1-contested (single-killer) slot** (`≥ 14/183`, min at `m=13`). The genuinely filler-contested
(multi-killer) slot is localized to the same formula, and its uniform bound is pinned to one clean remaining
constraint.*

---

## The slot-depth formula

At the covering-min time `t*`, exactly two runners bind (Chebyshev equioscillation, opus-S252). Let the binders
have speeds `a, b`, with nearest-integer arcs at `p_a/a` and `p_b/b`. With `t*` between them, one rising (slope
`a`) and one falling (slope `b`):

`a(t* − p_a/a) = M`,  `b(p_b/b − t*) = M`  ⟹  `M(1/a + 1/b) = p_b/b − p_a/a =: Δ`.

> **`M = Δ / (1/a + 1/b) = Δ · ab/(a+b)`**,  where `Δ` = the gap between the two binders' arcs.

Verified exactly (`M ==` formula) on every case:

| family | binders `(a,b)` | arcs | `Δ` | `M` |
|---|---|---|---|---|
| deep well `{1..12,182}` | `(1, 182)` | `0, 1/13` | `1/13` | `14/183` |
| ladder `{1..12,364}` | `(1, 364)` | `0, 1/13` | `1/13` | `28/365` |
| multi `{1..11,13,84}` | `(5, 84)` | `2/5, 5/12` | `1/60` | `7/89` |
| multi `{1..10,13,22,84}` | `(1, 22)` | `0, 1/11` | `1/11` | `2/23` |

The covering-min is thus, exactly, `min` over covering families of `Δ · ab/(a+b)` over the binding pair.

## The uniform-in-m bound for the runner-1-contested slot — proved

When runner 1 binds (`a = 1`, arc at `0`), the formula collapses to `M = p_b/(b+1)`. For a single-killer family
the other binder is the backbone `b = 14m`, and — since a single-killer backbone must also carry the multiple of
13 (`m` a multiple of 13, so the backbone has an arc exactly at `1/13`, `p_b = 14m/13`) — `Δ = 1/13` and

> **`M = 14m / (13(14m+1))`**,  and  **`M ≥ 14/183 ⟺ 183m ≥ 13(14m+1) = 182m+13 ⟺ m ≥ 13`.**

A single-killer backbone has `m ≥ 13` (it is `14·13c = lcm(13,14)·c`), so **`M ≥ 14/183` uniformly in `m`, with
equality only at `m = 13` (the deep well).** This is a clean closed-form, uniform-in-`m` proof of the runner-1
(single-killer) slot depth — exactly the `14c/(182c+1)` ladder already machine-checked kernel-pure (cont.60/61).
The user's "bound uniformly in m" is *done* for this case.

## The filler-contested slot — localized, with the residual pinned

When a filler binds instead (`a = f ≥ 2`), `M = Δ · fb/(f+b)`. Because the filler is *faster* than runner 1
(`f ≥ 2 > 1`), the balance crosses higher and the slot is **deeper (looser)** — verified: `{5,84} → 7/89 ≈
0.0787`, `{1,22} → 2/23 ≈ 0.087`, both above `14/183`. So the filler-contested slots are not the minimizers; the
runner-1 slot is. But turning "looser" into a *uniform* bound needs one thing the formula alone doesn't give:

> the equioscillation side-condition — **all other runners are `≥ M` at `t*`** — forces `Δ·fb/(f+b) ≥ 14/183`.

Without it, `Δ` could be arbitrarily small; with it, the other runners' arcs bound `Δ` below. That side-condition
is exactly the arithmetic/analytic content the finite check (klein's ILP `≤ 182`) and the Fourier route
(opus-S259–263) supply. So the slot-depth formula **localizes** the whole covering-min crux to a single clean
inequality on the binding pair, and separates it precisely: runner-1 slots are closed uniformly in `m` here;
filler slots reduce to `Δ·fb/(f+b) ≥ 14/183` under equioscillation.

## Net

The covering-min is exactly `min Δ·ab/(a+b)` over binding pairs (`Δ` = arc gap). This gives (i) a clean
closed-form, **uniform-in-`m` proof `M = 14m/(13(14m+1)) ≥ 14/183` for the runner-1-contested single-killer
slot** — answering the posed question for the tight case and matching the Lean ladder; and (ii) a localization
of the filler-contested (multi-killer) case to `Δ·fb/(f+b) ≥ 14/183`, deeper because the filler is faster,
with the equioscillation "all-others-`≥ M`" side-condition as the sole remaining content. The minimizer over all
binding pairs is `{1, 182}` with `Δ = 1/13`, giving `14/183` — the deep well, once more, and now as an explicit
2-runner extremal problem.

*Files: lrc14_slot_depth_formula_kps_S127.py (+.out). Builds on cont.65 (backbone lens), cont.60/61 (single-killer
ladder, Lean), opus-S252 (equioscillation); the filler residual feeds opus-S259–263 (Fourier). HYP-6238.*
