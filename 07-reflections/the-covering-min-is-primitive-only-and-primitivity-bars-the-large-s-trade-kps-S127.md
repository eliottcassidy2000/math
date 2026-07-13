# The covering-min is primitive-only, and primitivity is exactly what bars opus's large-s trade

*kind-pasteur-2026-07-11-S127 cont.57. Owner: "keep pushing the odd-doorway packing argument toward a proof."
Pressure-testing opus-S253's slow-fast balance `M = M_core·v_f/(v_f+s)` against dilation invariance surfaced a
concrete family that is DC yet sits below `14/183` — which turns out not to be a refutation but the precise
witness for why the crux needs primitivity, and why opus's one isolated open obstacle (the "large-`s` trade")
never actually materializes.*

---

## The witnesses: imprimitive DC reaches the global floor 1/14

`M(cv) = M(v)` (dilation invariance, THM-531), but **divisor-completeness is not dilation-invariant**. So a
*dilated* family can be DC while its primitive core is not. Three exact witnesses (all genuinely DC — a multiple
of every `d ∈ {2..14}` — all with `M < 14/183`):

| family | gcd | DC? | primitive reduction | `M` | clears at base |
|---|---|---|---|---|---|
| `2·{1..13} = {2,4,…,26}` | 2 | **yes** | `{1..13}` (AP, non-DC) | **1/14** | 28 (even) |
| `2·GW = {2,…,22,26,48}` | 2 | **yes** | GW (non-DC) | **1/14** | 28 (even) |
| `2·{1..12,91} = {2,…,24,182}` | 2 | **yes** | `{1..12,91}` (non-DC) | **7/92** | 184 (even) |
| deep well `{1..12,182}` | 1 | yes | itself (DC) | 14/183 | 183 (odd) |

`2·{1..13}` is DC and has `M = 1/14` — the *global* LRC floor, a full `1/84` below the covering-min. So:

> **The crux is "*primitive* DC ⟹ `M ≥ 14/183`."** Drop "primitive" and it is flatly false: imprimitive DC
> families reach all the way down to `1/14`. The fleet's standing "*primitive* q-covering" qualifier
> (HYP-2566/2579/2589/…) is **load-bearing, not a normalization convenience.**

Why this is not a refutation of anything: an imprimitive DC family `cv` dilates to a primitive `v` that is
**non-DC** (the dilation manufactured the divisors), so `v` lands in the THM-366 bucket with `M ≥ 1/14` — LRC
holds. The logical order matters: for the covering-min you may **not** say "take a DC family, WLOG primitive"
(dilating breaks DC); you must say "take a *primitive* DC family." For LRC(14) itself (target `1/14`) the
reduction is innocent; for the `14/183` crux it is not.

## The parity/primitivity link (sharpening cont.56)

This snaps onto the two-doorways parity picture. Every witness clears at an **even** base (28, 28, 184) — because
each is (a multiple of) a **non-DC** primitive, which lives at the even doorway. Only the **primitive** DC family
sits at the **odd** doorway `Φ₆ = 183`. So the odd doorway *is* exactly the primitive-DC bucket:

> **odd doorway (base `Φ₆`, band-packing, `M ≥ 14/183`) = the primitive-DC bucket.** Imprimitivity is an
> **even-base escape hatch** that slides a DC family down to its non-DC primitive core, i.e. down to `1/14`.

## The push: primitivity dissolves the large-s trade

opus-S253 proved the interval-core single-killer case (`M = M_core·v_f/(v_f+s) ≥ n/Φ₆`, `s=1`) and isolated the
sole open obstacle for the general single-killer bound: a **large-`s` trade** — a core whose binding runner is
fast (large `s`), shrinking `v_f/(v_f+s)`, against a larger `v_f`. The balance shows that even a modest `s`
would break the bound: with `M_core = 1/13` and `v_f = 182`,

> `M ≥ 14/183` needs `182/(13(182+s)) ≥ 14/183`, i.e. **`s ≤ 1`.** An `s = 2` extremal core gives `M = 7/92 <
> 14/183` — exactly the witness above.

So the whole question is whether a primitive covering family can have an `M_core = 1/13` core with binding speed
`s ≥ 2`. **It cannot.** An `M_core = 1/13` core is an LRC(13) extremal, i.e. essentially `c·{1..12}` (up to the
GW variants), with binding speed `s = c` and core optimum `t_core = 1/(13c)`. A resonant killer needs
`13c ∣ v_f`, and covering the divisors `c·{1..12}` misses (for `c=2`: only `d=13`, since `2·7=14` is already
present) needs `v_f` a multiple of `13`; together `v_f` is a multiple of `lcm(13, 26) = 26`, hence **even** —
sharing the factor `2` with the core, so the whole family is **imprimitive**. Primitivity bars it. The same
factor-sharing obstruction kills every `c ≥ 2`.

> **Primitivity forces `s = 1` for the extremal-core single-killer family.** opus's proved `s=1` case is therefore
> not a lucky special assumption — it is the *only* case primitivity permits at the binding `M_core = 1/13`. The
> large-`s` trade opus worried about does not materialize for single-killer primitive covering families; the
> mechanism is the coupling of the core's gcd to the killer's parity through the lcm resonance.

## Net — and the honest limit

The odd-doorway packing argument gains two things. **Scope:** the theorem to prove is `primitive DC ⟹ M ≥
14/183`, and the primitivity is necessary, with `1/14` witnesses proving it. **Content:** for the extremal-core
single-killer case — the binding case opus reduced everything to — primitivity *forces* `s=1`, closing opus's
one isolated obstacle there. What remains is exactly what remains for opus-S253: **multi-killer** primitive
covering families (several resonant runners cleared simultaneously) and **non-extremal cores** (`M_core > 1/13`,
which have slack). Those are the genuine open frontier; the single-killer binding case is now controlled on both
the `v_f ≥ n(n−1)` side (the lcm-outlier bound, cont.55) and the `s = 1` side (primitivity, here).

*Files: lrc14_covering_min_needs_primitive_kps_S127.py (+.out). Sharpens opus-S253 (slow-fast balance; the
large-s trade), builds on kps cont.55 (lcm-outlier v_f ≥ n(n−1)) and cont.56 (two doorways / parity), rests on
THM-366 (non-DC ⟹ 1/14), THM-531 (dilation), klein-S267 (14/183 primitive covering-min). HYP-6222.*
