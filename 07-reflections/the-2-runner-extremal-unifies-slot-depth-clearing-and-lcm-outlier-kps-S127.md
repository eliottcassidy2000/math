# The 2-runner extremal problem unifies the slot-depth formula, klein's clearing modulus, and the lcm-outlier

*kind-pasteur-2026-07-11-S127 cont.67. Owner: "work the 2-runner extremal problem." Working the covering-min in
its 2-runner form `min Δ·ab/(a+b)` (cont.66), I found it is the same object as klein's clearing-modulus
certificate and my own lcm-outlier mechanism — three frames the fleet had developed separately are one
extremal problem. The covering-min is `min` over binding pairs, and each frame is a different reading of it.*

---

## The 2-runner extremal problem

The covering-min is (cont.66): over primitive covering families, minimize `M = Δ·ab/(a+b)`, where `{a,b}` is the
binding pair (speeds), `Δ` the gap between their arcs `p_a/a, p_b/b`, subject to equioscillation (all other
runners `≥ M` at `t*`). Verified minimizer: `{1, 182}`, `Δ = 1/13`, `M = 14/183` — the deep well.

## Reading 1 = Reading 2: slot-depth **is** klein's clearing modulus

When runner 1 binds (`a = 1`, arc at `0`), the slot formula gives `M = p_b/(b+1)`. That is **exactly a
band-clearing certificate at modulus `q = b+1`**, band-edge `μ = p_b`, `M = μ/q`. Verified:

| family | binding pair | `M` | `q = b+1` | `μ = p_b` | `M = μ/q` |
|---|---|---|---|---|---|
| deep well `{1..12,182}` | `(1,182)` | 14/183 | **183** | 14 | ✓ |
| ladder `{1..12,364}` | `(1,364)` | 28/365 | **365** | 28 | ✓ |
| multi `{1..10,13,22,84}` | `(1,22)` | 2/23 | **23** | 2 | ✓ |

So **the co-binder is the clearing modulus minus one, `b = q−1`.** My slot-depth formula (kps cont.66) and
klein's clearing-modulus framework (the `q ≤ 183` certificate) are the *same* 2-runner extremal — the co-binder
`b` and the modulus `q` are the same datum. The deep well's `b = 182` is klein's `q = 183` minus one.

## Reading 3: the co-binder **is** the lcm-outlier

Why is the minimizing co-binder `182`? In the 2-runner frame, at the extremal `t* = 14/183` (`q = 183`), the
covering requirement `d = 13` forces a multiple of 13 that is lonely (`≥ M`): `13u` with `‖13u·14/183‖ ≥
14/183`. Since `182u ≡ −u (mod 183)`, this holds iff `u ∈ [14, 169]`, so the **smallest lonely multiple of 13
is `13·14 = 182 = lcm(13,14)`** — which is also the multiple of 14, and the co-binder. This is exactly the
lcm-outlier mechanism (kps cont.55) — *the minimal single speed repairing both divisors the AP misses* — now
seen as the co-binder constraint of the 2-runner extremal. The three frames coincide on the single number `182`.

## The one problem, three readings — and where each closes it

> **Covering-min = the 2-runner extremal `min Δ·ab/(a+b)`.** Its co-binder `b` is simultaneously the slot's
> second binder (kps slot-depth), `q−1` for a clearing modulus `q` (klein), and the lcm-outlier `lcm(13,14) =
> 182` at the minimizer (kps). One extremal problem, minimizer `{1,182}`, value `14/183`.

The split is clean and the proved boundary sharp:

- **Co-binder a multiple of `lcm(13,14) = 182`** (it repairs 13 and 14 at once → single-killer): the extremal is
  `M = 14m/(13(14m+1)) ≥ 14/183` uniformly in `m`, **proved and machine-checked kernel-pure** (cont.60/61/66).
  Here all three readings agree in closed form.
- **Co-binder a multiple of 14 only, with a separate multiple of 13** (multi-killer): the two clearing
  runners differ, `t*` moves mid-interval, a filler binds, and the bound is klein's finite clearing certificate
  (`q ≤ 183`, ILP `≤ 182`) plus the Fourier route (opus-S259–263). This is the open residual, now identified
  across all three frames as *the same* case.

## Net

The 2-runner extremal problem is not a new object — it is the covering-min — but working it shows that the
fleet's three separate handles (my slot-depth balance, klein's clearing modulus, my lcm-outlier) are one and the
same, tied through `b = q−1 = lcm(13,14)` at the minimizer. This unifies the bookkeeping (a clearing certificate
*is* a binding-pair balance) and pins the open case identically in every frame: the co-binder splitting `13` and
`14` into separate runners. The proved half — co-binder `= lcm`-multiple — is closed in closed form and in Lean;
the open half is the finite clearing check klein certifies to `182`.

*Files: lrc14_two_runner_extremal_kps_S127.py (+.out). Unifies kps cont.66 (slot-depth), cont.55 (lcm-outlier),
klein's clearing-modulus certificates (q ≤ 183), and the single-killer Lean (cont.60/61). HYP-6242.*
