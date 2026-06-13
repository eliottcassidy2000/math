---
source: opus-2026-06-02-S571 (remote-control)
status: REDUCTION (clean) + rigorous partial (large multiples) + verification n=4..14; the lift lemma collapses to C' "n|v ⟹ loose"
tags: [LRC, lift-lemma, C-prime, multiple-of-n, no-multiple-witness, THM-369, S564, dodge, LRC-induction, n14]
---

# The lift lemma collapses to C′: "a multiple of n forces positive measure"

**Prompt (user):** work on the lift lemma; look for simplifications and speedups.

The lift lemma (`measure-0 ⟹ unblocked small pair`, S569/S570) **splits cleanly by
one bit — whether the config contains a multiple of `n`** — and the whole of LRC(n)
collapses onto a single, sharper, *already-conjectured* statement.

## 1. The split (the simplification)

For a speed set `S` (`n-1` runners), test `t = 1/n`:

- **No multiple of `n` in `S`:** every `‖v_i/n‖ = min(v_i mod n, n-v_i mod n)/n ≥ 1/n`,
  so **`t=1/n` is a witness, `M(S) ≥ 1/n`.** (THM-369, the divisibility sieve.) Free.
- **A multiple of `n` in `S` (`v=nw`):** `t=1/n` (and every `j/n`) is *killed* by `v`
  (it sits at `0`). This is the only hard case.

So the lift lemma's entire content is the hard case, and that case is exactly:

> **C′ (already conjectured, S564 as a property of tight sets):**
> `n | v` for some `v ∈ S` `⟹ M(S) > 1/n` (positive measure / loose).

## 2. The reduction (the headline)

> **C′ ⟹ LRC(n).** Every config has no multiple of `n` (→ `1/n` witness, `M≥1/n`)
> **or** has one (→ `M>1/n` by C′). Either way `M≥1/n`. ∎

S564 noticed C′ as a *symptom* of tight sets ("tight ⟹ no multiple of `n`"). The new
point is the **converse direction is the whole conjecture**: proving "a multiple of
`n` forces loose" *is* proving LRC(n). All of LRC concentrates into one structural
statement about a single distinguished runner — the one that vanishes on the entire
`n`-clock.

## 3. Verification — C′ holds, robustly, every multiplier size (n=4..14)

`lrc_lift_lemma_measure_bound_s571.py`, `lrc_lift_Cprime_residual_s571.py`:

| n | `v=n` exactly (hardest, w=1) | small mult (w≤3) | large mult (w≥4) |
|---|---|---|---|
| 6 | loose 790/790, minμ .025 | 794/794, .013 | 800/800, .036 |
| 8 | 800/800, .020 | 800/800, .022 | 800/800, .039 |
| 10 | 800/800, .022 | 800/800, .026 | 800/800, .023 |
| 12 | 800/800, .035 | 800/800, .034 | 800/800, — |
| 14 | 800/800, .033 | 800/800, .023 | 800/800, .024 |

**0 tight-with-multiple anywhere** (also 0 over the STEP-1 exhaustive small boxes,
n=4..14). The safe measure is always strictly positive — even for `v=n` (`w=1`), the
runner that kills every clock point.

## 4. Why the easy proof fails, and the rigorous partial that works

**The crude global bound fails.** `μ(safe S) ≥ μ(safe S') − μ(D_v) = μ(safe S') − 2/n`
(removing the multiple `v`, `S'` the other `n-2` runners). But the min safe-measure
of `(n-2)`-runner configs at level `1/n` is **far below `2/n`** (e.g. `n=6`,
`S'={1,3,4,5}`: `μ'=0.05 ≪ 1/3`). So subtracting the whole `2/n` is too lossy — the
real mechanism is that `D_v`'s **thin, evenly-spaced arcs cannot cover** `safe(S')`.

**Rigorous partial — large multiples are dodgeable (using the *proven* LRC(n−1)):**

> **Lemma.** If `v = nw > (n-1)·max{other speeds}`, then `M(S) > 1/n`.
>
> *Proof.* `S' = S\{v}` has `n-2` runners; by LRC(n−1) (proven for `n-1 ≤ 13`)
> there is `t₀` with `min_i ‖v_i t₀‖ ≥ 1/(n-1)`. On the interval `I` of half-width
> `δ = 1/(n(n-1)V')` (`V' = max other speed`) around `t₀`, every `v_i` stays `> 1/n`
> (it falls at rate `≤ v_i`). The danger arcs of `v=nw` have radius `ρ = 1/(n²w)`.
> When `δ > ρ` — i.e. `nw > (n-1)V'`, exactly the hypothesis — `I` is wider than one
> arc, so it contains a `v`-safe sub-interval on which all of `S` exceed `1/n`:
> positive measure, `M(S) > 1/n`. ∎*

This **cleanly uses the literature's proven `n=13`** to push on `n=14`'s
multiple-of-14 configs: any multiple `14w > 13·(max other speed)` is handled.

## 5. The residual (sharp)

What remains of C′ (hence of LRC(n)) is the **small multiples**: `n | v` with
`v ≤ (n-1)·max(others)`, down to `v=n` itself. There the arc radius `ρ` is no longer
below the dodge window `δ`, so the one-arc dodge is not guaranteed — yet the configs
are still loose (verified). The missing step is **equidistribution**: the evenly-
spaced arcs of `nw` (period `1/(nw)`, total danger `2/n`) meet `safe(S')` in
`≈ (2/n)·μ(safe S') < μ(safe S')`, so they cannot cover it. Proving the arcs can't
*align to cover* `safe(S')` is the honest open core — a three-distance / Weyl-
equidistribution statement about a single arithmetic progression against a fixed
union of intervals.

## 6. Speedups delivered

- **Verification speedup:** the lift lemma no longer needs `M(S)` or measure on all
  configs. Split: no-multiple-of-`n` configs are **free** (`1/n` witness); only the
  **multiple-of-`n` slice** needs checking — and there only "loose," via a fast
  positive-measure certificate (one safe interval), not exact `M`.
- **Measure speedup:** `safe_measure` via **endpoint enumeration** (the
  `(kn±1)/(nv)` breakpoints) is exact and grid-free — `O(Σv_i · m)` per config.
- **Scope speedup:** the residual is now *small multiples only* (`v ≤ (n-1)·max`),
  a thin band, with large multiples proven dodgeable from LRC(n−1).

## 7. Honest status

- **Reduction `C′ ⟹ LRC(n)`: rigorous** (THM-369 + the split).
- **C′ verified** n=4..14, all multiplier sizes, 0 exceptions (incl. the `v=n` worst case).
- **Large-multiple dodge: rigorous** (Lemma §4, on proven LRC(n−1)).
- **Small-multiple residual: open** — the arc-equidistribution non-covering. This is
  the genuine remaining heart, now *isolated and named*.

This is the simplification asked for: LRC(n) ≡ "a multiple of `n` can't be tight,"
the no-multiple half is trivial, large multiples are dodgeable from the proven
predecessor, and only the small-multiple arc-cover remains.

**Artifacts:** `04-computation/lrc_lift_lemma_measure_bound_s571.py` (+`.out`),
`lrc_lift_Cprime_residual_s571.py` (+`.out`). Builds on S569/HYP-2095 (lift lemma),
S570/HYP-2097 (64-class container), S564 (C′ for tight sets), THM-369 (divisibility
sieve), proven LRC(n−1). New: **HYP-2102**.
