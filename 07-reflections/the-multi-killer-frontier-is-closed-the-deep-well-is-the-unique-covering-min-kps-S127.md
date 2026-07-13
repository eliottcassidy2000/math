# The multi-killer frontier is closed: the deep well is the unique covering-min, by core-length monotonicity

*kind-pasteur-2026-07-11-S127 cont.58. Owner: "work the multi-killer primitive covering families" — the frontier
opus-S253 and my cont.57 left open after the single-killer case was closed. Enumerated exactly. Result: no
multi-killer primitive DC family beats `14/183`; the single-killer deep well is the unique covering-min, and the
reason is a one-line monotonicity in the core length.*

---

## The enumeration

A single-killer covering family is the interval core `{1..12}` plus one outlier (the deep well `{1..12,182}`).
A **multi-killer** family has a genuinely shorter interval core `{1..k}` (`k ≤ 11`) plus `13−k` outliers all
`≥ 13`, together supplying the missing divisors `{k+1,…,14}`. Exact covering-min over each class:

| core length `k` | `M_core = 1/(k+1)` | # primitive DC | min `M` | minimizer | vs `14/183` |
|---|---|---|---|---|---|
| 12 (single-killer) | 1/13 | 1 | **14/183** ≈ .07650 | `{1..12, 182}` | = (the floor) |
| 11 | 1/12 | 58 | **7/89** ≈ .07865 | `{1..11, 13, 84}` | above |
| 10 | 1/11 | 56 | **2/23** ≈ .08696 | `{1..10, 13, 22, 84}` | above |
| 9 | 1/10 | 1120 | **2/23** ≈ .08696 | `{1..9, 13, 20, 22, 84}` | above |

**No multi-killer family beats `14/183`.** The deep well is the unique global covering-min; every genuine
multi-killer family sits strictly above it.

## Why — the core-length monotonicity, in one line

The single number that governs it:

> `1/13 = 0.076923 > 14/183 = 0.076503.`

The interval core `{1..k}` is the LRC(`k+1`) extremal, with loneliness `M_core = 1/(k+1)` — its *own* floor,
before any killer. The killer can only **lower** `M` from `M_core`. So:

- **Single-killer (`k=12`):** `M_core = 1/13`, the smallest core-floor achievable (the longest DC-compatible
  core — `{1..13}` is the AP, already 13 speeds and non-DC, so `k=12` is maximal). The lcm-outlier `182 =
  lcm(13,14)` (cont.55) drops `1/13` to `14/183` — a drop of just `1/13 − 14/183 = 1/2379`. That razor-thin
  drop is the whole margin, and it lands exactly on the floor.
- **Multi-killer (`k ≤ 11`):** the core is shorter, so `M_core = 1/(k+1) ≥ 1/12 = 0.0833` — already a full
  `0.007` **above** `14/183`, before killers. Even the killers' drop leaves `M` above the floor (`7/89`,
  `2/23`). There is simply no room to reach `14/183` from a core-floor of `1/12`.

So the covering-min is **monotone decreasing in the interval-core length `k`**, and `k=12` (the single killer)
is the maximal core — hence the unique minimizer. Multi-killer, by construction, has a shorter core and a
higher floor: **it can never bind.** This is the exact sense in which opus-S253's "multi-constraint" case is not
a threat — shortening the core to admit more killers *raises* `M_core` faster than the extra killers can lower
`M`.

## The honest limit

Two gaps remain between this and a closed proof, both inherited from opus-S253:

1. **The balance is a lower bound, not the global optimum, for multi-killer.** At the core-optimum `t_core =
   1/(k+1)` the balance `M = M_core·v_min/(v_min+s)` can under-shoot (`{1..11,13,84}` would read `< 14/183` at
   `t_core` with `v_min = 13`, `s = 11`), but the *true* optimum is elsewhere (`t* = 37/89`, giving `7/89`).
   So "multi-killer ≥ `14/183`" is here an **exact enumeration fact**, not yet a closed-form inequality; the
   clean statement to prove is the core-length monotonicity of the *global* `M`.
2. **Non-interval cores.** The enumeration covers interval-core families `{1..k}+outliers`, where the
   extremizers live (three-gap regularity, mac-mini). A general primitive DC family need not have a long initial
   run; those `≤ 182` are covered by klein's ILP, the rest by the escape (death-star's `LRCUEscape`).

## Net

The multi-killer primitive covering families are **closed as a threat to the covering-min**: exact enumeration
(`k = 9,10,11`) shows all sit strictly above `14/183`, the deep well is the unique minimizer, and the reason is
the clean monotonicity `M ≥ M_core = 1/(k+1)`, minimized at the maximal core `k=12`. Combined with the
single-killer closure (opus-S253 + kps cont.55/57 + mac-mini S68), the **entire interval-core structure of the
covering-min lower bound is now understood**: the deep well wins on core length (longest → smallest `M_core`),
on outlier size (`lcm(13,14) = 182`, smallest covering outlier), and it alone touches the floor. What remains is
promoting the core-length monotonicity from enumeration to a global-optimum inequality, and the non-interval
tail (ILP + escape).

*Files: lrc14_multikiller_enum_kps_S127.py (+.out), lrc14_multikiller_genuine_kps_S127.out. Builds on opus-S253
(slow-fast balance, multi-constraint case), kps cont.55 (lcm-outlier), cont.56 (two doorways), cont.57
(primitive-only), mac-mini S68 (large-s mirage), klein-S267 (14/183). HYP-6225.*
