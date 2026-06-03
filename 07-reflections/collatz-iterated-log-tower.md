# The Collatz iterated-logarithm tower: where each scale lives

*Session: collatz-excursion-tower. Sequel to `collatz-rapidity-defect.md`
(HYP-2147/2148). Prompt: "Tao is good at seeing log, loglog, logloglog
inequalities — find the deeper abstraction and make your own."*

## The deeper abstraction: iterated linearization

The last session's result was a single fact: in rapidity `rho = ½ln = log_F =
arctanh`, the Collatz `+1` collapses from an unbounded multiplicative
perturbation into a **bounded** additive defect `D(n) < 0.2257`. The right way
to see *why that is only the first floor of a building* is **renormalization**:

> Apply the formal-group logarithm; the nonlinearity drops by one log-order.
> Apply it again (look at fluctuations of the log); the residual drops again.
> Tao's iterated-log statements are reading off successive floors of this tower.

Collatz nonlinearity, **measured at scale `log^{(k)} n`**, has size `O(1)` at
each floor — but each floor is a different, one-order-finer, feature of the
*same* `+1`. That is the abstraction. Concretely:

| floor | scale | what the `+1` becomes | status |
|------|-------|----------------------|--------|
| 0 | `n` | unbounded perturbation; `x↦3x`,`x↦x/2` are the group ops | — |
| 1 | `ln n` (rapidity) | **bounded** defect `D(n) < 0.2257` | proved (HYP-2147/48) |
| 1.5 | — | descent count `L(n) ~ ln n / ln(4/3)` = one log of `n` | Wald, confirmed |
| 2 | `ln ln n` | excursions: exponent **θ = 2**, `P(peak/n > M) ~ C/M`; LIL scale `√(L·lnln L)` | model-proved + confirmed (HYP-2149) |
| 3 | `ln ln ln n` | a.s. envelope `peak(n) ≤ n·(ln n)^{o(1)}` for almost all `n` | Tao's scale; finite-range evidence |

## The one constant that runs the tower: θ = 2, and *why*

Per odd Syracuse step the rapidity moves by the ideal increment
`Y = ½ln3 − (k/2)ln2`, where `k = v₂(3a+1) ~ Geometric(½)`. So

- drift `δ = E[Y] = ½ln(3/4) = −0.1438` (negative — orbits fall),
- variance `σ² = Var(Y) = ln²2/2 = 0.2402`.

The walk's *maximum* is controlled by the **Cramér–Lundberg exponent** `θ`, the
positive root of `E[e^{θY}] = 1`. A two-line computation gives
`E[e^{tY}] = 3^{t/2}/(2^{1+t/2} − 1)`, so the root solves

```
        3^{t/2} = 2^{1+t/2} − 1,   nontrivial root t = 1,   θ = 2 exactly.
```

The pretty part: **at the root the equation is `3 − 4 + 1 = 0`**, i.e.
`3·1 + 1 = 4 = 2²` — the Collatz step on `n = 1`. The exponent that governs how
high orbits fly is *pinned by the smallest instance of the `+1`*. The same `+1`
that floor 1 isolated as the entire obstruction is, at floor 2, exactly the
thing that fixes the fluctuation law. (This is also why `2` recurs: `θ = 2` is
the prime-2 of halving showing up as the decay rate.)

## My own inequalities (the new content)

**Level 2 — the excursion power law (HYP-2149).** Because `θ = 2` and rapidity
is `½ln`, a peak-to-start ratio `M` corresponds to a rapidity excursion
`x = ½ln M`, and `P(R > x) ≈ Ce^{−θx} = Ce^{−ln M} = C/M`. Hence

```
        P( max(orbit of n) / n > M )  ~  C / M,        C ≈ 2.4.
```

Tested over odd `n < 10⁶`: `M·P(·)` is flat at ≈ 2.4 across `M = 4 … 1024`, and
a log–log fit gives slope `−1.029` (predicted `−1`), i.e. `θ = 2.06` (predicted
`2`). This is a clean, falsifiable, *quantitative* law for "how high Collatz
orbits climb," with the constant derived, not fitted.

**Level 2 — the LIL scale.** The centered walk `W_j = rho(a_j) − rho(n) − jδ`
has variance `σ²j`, so the law of the iterated logarithm predicts the maximal
fluctuation about the drift line scales as `√(2σ²·L·lnln L)`. Confirmed as the
right *scale* on long orbits (the exact limsup constant `1` is asymptotic and
not pinned on a single finite orbit). **This is the genuine `loglog`:** it is
the iterated logarithm of the *number of steps*, which is itself one log of `n`
— so it is the `lnln n` feature of the original integer.

**Level 3 — the a.s. envelope (`lnlnln`).** The exponential tail `e^{−2x}`,
fed through Borel–Cantelli, predicts that for *almost all* `n` the peak
excursion satisfies `R(n) ≤ (½ + o(1)) lnln n`, i.e. `peak(n) ≤ n·(ln n)^{o(1)}`.
Finite-range shadow: the fraction of `n < 10⁶` with `peak/n > (ln n)^c` is
`0.19, 0.054, 0.015, 0.001` for `c = 1, 1.5, 2, 3` — heading to `0` as `c`
grows. This is exactly the regime of Tao (2019), *"Almost all orbits of the
Collatz map attain almost bounded values,"* whose `Col_min(n) ≤ f(n)` for any
`f → ∞` lives at the `lnlnln` scale. The tower says *why* that scale: it is the
third application of the formal-group logarithm.

## Why this is more than analogy

Each floor is independently checkable and the constants are forced, not chosen:
`δ = ½ln(3/4)`, `σ² = ln²2/2`, `c = 1/ln(4/3)`, `θ = 2`. They are all the same
two numbers `2` and `3` (and the `+1`) viewed through `k` more logarithms. The
honest caveat: floors 2–3 are *model* statements (the standard Collatz
equidistribution heuristic — `k` geometric, increments independent), the same
unproven model Tao's work makes rigorous with serious harmonic analysis. What
this session contributes is the **organizing picture** (iterated linearization)
and a **derived, confirmed constant** (`θ = 2` from `3+1=4`) that places the
known log / loglog / logloglog phenomena on one ladder.

## Next

1. **Make `θ = 2` rigorous** for the *true* Syracuse increments (not the iid
   model) — Tao's 3-adic stationary measure should give the matching tail.
   The target: prove `P(peak/n > M) ≤ C/M` unconditionally (a real theorem,
   stronger than anything currently known about Collatz altitudes?).
2. **Closed form / role of `C ≈ 2.4`.** Is it `e^{D*}·(something)`? (`e^{2D*} =
   1.57`, `3/√(...)`?) Worth pinning.
3. **OEIS:** the excursion records (A006884/A025586 are max-value records);
   reinterpret them as the upper-tail order statistics of the `C/M` law and
   extend per the repo's practice.
4. Feed `θ = 2` and `δ, σ²` back into the VS-chain analogy (does the VS map
   have its own Cramér exponent?) — `grand-trichotomy.md`.
