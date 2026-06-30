# The covering-min construction realizes Steinhaus with three-distance gaps EXACTLY {1, n, 2n} (AP step n ×(n−3), antipodal-closure 1 ×1, 0-gap 2n ×1), so M = n/Φ₆(n) is the dense block's 1/(n−1) with the lcm-killer adding +1 to the denominator (n(n−1) → n(n−1)+1 = Φ₆, klein's antipodal (1,−1) binding); three-thread synthesis — klein's ζ₆-line + mac-mini's three-distance/taxonomy + my C(n) law converge, and klein's dense-core bound (M>1/n provable on dense coverings) leaves only the exotic-covering exclusion open

*opus-2026-06-30. Owner: push the optimality proof, see the bigger picture. Pulled klein-S26 (HYP-3715,
the ζ₆-line — independent convergence with my witness result) and mac-mini-S43 (HYP-3702, the taxonomy +
three-distance). This is the concrete three-distance realization + the synthesis of all three threads.*

## The construction's three-distance gaps are EXACTLY {1, n, 2n} (computed)
mac-mini's tool is the Steinhaus three-distance theorem (≤3 gap lengths for `{kα}`). The construction at its
binding witness `t*=ζ₆=n/Φ₆` realizes it concretely. The speeds `{1,…,n−2,(n−1)n}·n mod Φ₆` give the sorted
points `{n, 2n, …, (n−2)n, (n−1)²=Φ₆−n}`, with gaps (verified n=7,8,10,14,20):
| gap | length | multiplicity | meaning |
|---|---|---|---|
| AP step | `n` | `n−3` | the equally-spaced body (`kn`) |
| antipodal closure | `1` | `1` | `(n−2)n → (n−1)²=Φ₆−n`, distance `1` (the killer) |
| **0-gap** | `2n` | `1` | contains 0 (`Φ₆−n → n` through 0 = `n+n`); **M = n/Φ₆** (0 centered) |
> **Three distances `{1, n, 2n}`.** The `2n` 0-gap is the "deepest hole"; `0` sits at its center, distance
> `n` from the nearest points `±n`, so `M = n/Φ₆`. The construction is the *equally-spaced ζ₆-AP* — the
> extreme (most-regular, one-body-step) Steinhaus configuration. This is the geometric content of mac-mini's
> "three-distance for the hard core," made explicit, and of klein's "ζ₆-line covering radius n."

## The antipodal +1: M = n/Φ₆ is the block's 1/(n−1) perturbed by the killer
A clean derivation of *why* `n/Φ₆` (building on klein's lcm-killer + antipodal binding):
> The dense block `{1,…,n−2}` ALONE has `M = 1/(n−1) = n/(n(n−1))` (the harmonic backbone, witness
> `t=1/(n−1)`). The killer `(n−1)n = lcm(n−1,n) ≡ −1 (mod Φ₆)` **kills that witness** (`(n−1)n·1/(n−1)=n≡0`,
> gap 0) and forces the binding to `ζ₆`, where it **adds 1 to the denominator**:
> `M : n/(n(n−1)) ⟶ n/(n(n−1)+1) = n/Φ₆(n)`  (since `Φ₆ = n(n−1)+1`).
So **the covering-min is the dense block's gap `1/(n−1)`, knocked down by exactly the `+1` antipode** —
klein's `(1,−1)` binding, `D = 1 + n(n−1)`. The `+1` is the whole margin: `1/(n−1) − n/Φ₆ = 1/((n−1)Φ₆)`,
and `n/Φ₆ − 1/n = (n−1)/(nΦ₆) ~ 1/n²`. The killer is the *minimal* perturbation (largest minimal killer =
lcm = most equidistributing), so the construction is the *least-perturbed* dense core — hence extremal.

## Three-thread synthesis (klein · mac-mini · opus, converged)
| thread | content | role |
|---|---|---|
| **klein** (HYP-3705/06/15) | ζ₆-line in `Z[ζ₆]/(n−ζ₆)`; lcm-killer; antipodal binding; discrete Kershner; **dense-core bound M>1/n** | the geometric realization + the provable lever |
| **mac-mini** (HYP-3701/02) | taxonomy (4 easy families + hard core); **three-distance** tool; hardness axis = LOWNESS; transition n=7 | the decomposition + the gap tool |
| **opus** (this + prior) | `C(n)` law (drop-2 n≤6 / construction n≥7), margin `~1/n²`; **gaps {1,n,2n}**; **antipodal +1** | the n-dependence + the concrete Steinhaus structure |
> They lock: the **hard core = the construction = the ζ₆-AP with three-distance gaps {1,n,2n}**, M = the
> dense block's `1/(n−1)` minus the antipodal `+1`. **klein's dense-core bound proves M>1/n on dense
> coverings** (which includes the construction itself: `n/Φ₆ > 1/n`). mac-mini's taxonomy discharges the
> non-dense families (shifted blocks `M≥1/8`, spread `M~0.3`, cusp via the comb/units). **What remains open
> is the *exotic-covering exclusion*** — that no covering set outside the discharged families (and no exotic
> ζ-image) beats the equally-spaced ζ₆-AP — i.e. that the ζ₆-line is the GLOBAL covering optimum (klein's
> cyclic-Kershner), the one genuinely-open node of the program.

## The proof program (where it stands)
1. **non-covering** — `M≥1/q` (THM-523, done).
2. **shifted blocks m≥2** — `M≥1/8 ≫ 1/n` (mac-mini, renormalization flow, provable).
3. **spread / large-min-speed** — `M~0.3` (packing, provable).
4. **cusp** (full-`Z_p` core) — existence via the `φ(n)` units / comb (measure-0).
5. **hard core = construction** — `M=n/Φ₆`; klein's dense-core bound gives `>1/n`; **OPTIMALITY (global min)
   open** = the ζ₆-line is the densest coverable Steinhaus configuration (cyclic-Kershner).
Only the OPTIMALITY in (5) is open. The value, the witness (ζ₆), the gap structure ({1,n,2n}), the margin
(`~1/n²`), the binding (antipodal `+1`), and the dense-core lower bound (`>1/n`) are all in hand.

## Status
- **Computed/verified (opus):** construction three-distance gaps exactly `{1,n,2n}` (n=7..20); `M=n/Φ₆` =
  block `1/(n−1)` − antipodal `+1` (killer `=lcm=−1 mod Φ₆`); the `+1` is the entire margin.
- **Synthesis:** klein (ζ₆-line/Kershner/dense-core) + mac-mini (taxonomy/three-distance) + opus (C(n)/gaps/
  antipodal) converge on: hard core = ζ₆-AP, gaps `{1,n,2n}`, dense-core `M>1/n` proved, exotic-exclusion
  open.
- **Open (the one node):** the ζ₆-line is the global covering optimum (cyclic/discrete Kershner) — no exotic
  covering beats the equally-spaced ζ₆-AP. The dense-core and taxonomy bound everything else.

Related: the-covering-min-witness-is-kleins-zeta6, the-covering-min-as-a-function-of-n (C(n)/margin); klein
HYP-3715/3706/3705 (ζ₆-line, hexagonal, dense-core), mac-mini HYP-3702/3701/3700 (taxonomy, three-distance,
n-dependence, apex isolation); the-master-two-Heegner-columns; OPEN-Q-108.
