# The metagraph H concentrates with CV² ~ 2/n (proved, the dispersion = a succession average), but the LRC covering ensemble does NOT — its dispersion plateaus. The contrast PINPOINTS why the rehearsal reaches the floor and not the cap: the metagraph has no tail; the LRC's tight tail is the whole difficulty

*opus-2026-06-29. Owner: return to the Var(H) succession closed form and push the G(n)/n!→1
concentration rate that feeds the LRC second-moment floor. The rate is exactly 2/n; transferring it to
the LRC reveals a sharp DISanalogy that locates the hard part.*

## The concentration rate (proved)
From `Var(H) = (n!/4^{n-1})(G(n)−n!)` (previous reflection):
> **`CV²(H) = Var(H)/E[H]² = G(n)/n! − 1 ~ 2/n`**, verified `n·CV² → 2` (1.0, 1.33, …, 1.985 at n=50,
> monotone to 2). Leading term `2(n−2)/(n(n−1))`; the full series sums to constant **2**.
- Exact dispersion as a **succession average:** `E[H²]/E[H]² = G(n)/n! = E_σ[\,2^{\#\text{asc-succ}(σ)}·
  \mathbf 1[\text{no desc-succ}(σ)]\,]` over a uniform random permutation `σ`. (`G(n)=1,2,8,32,158,928,…`,
  Riordan succession theory.)
- So **`H` concentrates fully**: the Hamiltonian-path count is `n!/2^{n-1}` with relative spread
  `√(2/n) → 0`. The "2" is the single-succession correction (weight `2^1`).

## The LRC does NOT concentrate the same way (computed)
The witness count `D(q,S)` over the **covering ensemble** has dispersion `CV²_S(D/q)` that **plateaus**,
it does not vanish:
| q | mean `D/q` | `CV²_S` | min `D/q` (tail) |
|---|---|---|---|
| 101 | 0.071 | 0.190 | 0.000 |
| 503 | 0.086 | 0.083 | 0.004 |
| 1009 | 0.079 | 0.083 | 0.006 |
| 2003 | 0.080 | 0.083 | 0.005 |
| 4001 | 0.081 | **0.082** | 0.005 |
> **`CV²_S → ~0.083`, a positive plateau** = the variance of the **singular series `L(S)` across covering
> sets**. Unlike the metagraph (`CV² → 0`, no tail), the LRC ensemble keeps a **persistent tight TAIL**
> (the near-tight sets like `{1..11,13,84}`, `min D/q ≈ 0.005 ≪ mean 0.08`).

## The diagnosis (why the rehearsal reaches the floor but not the cap)
This contrast is the cleanest statement yet of where the metagraph rehearsal stops:
> **The metagraph fully concentrates (`CV²→0`), so it models only the BULK — the generic floor (the
> `SL(2)`/Han-Lee "almost every covering set is lonely"). The LRC's difficulty is the persistent TAIL
> (`CV²_S→0.083`), the near-tight covering sets, which the metagraph — having NO tail — cannot model.
> The plateau value IS the floor-to-cap gap: the spread of `L(S)` that the bulk-concentration argument
> averages away but the tight cap must control.**
So `CV²(H)~2/n` is a genuine, exact concentration result, and it *feeds the floor* exactly as hoped —
but only the generic floor. The tight cap (`SL(4)`, the persistent tail) is provably beyond a
concentration/second-moment argument that treats the ensemble as homogeneous, because **the ensemble is
NOT homogeneous: it has a heavy tail of tight sets.** This is the second-moment face of "one dimension
past Littlewood."

## What this gives the program
1. **The floor is a concentration theorem and the metagraph proves its model:** `H` concentrates with the
   *exact* rate `2/n`; the LRC bulk concentrates likewise (Han-Lee), giving the generic floor rigorously.
2. **The cap is a TAIL theorem, not a concentration theorem:** the residual `CV²_S ≈ 0.083` plateau says
   the tight sets are a fixed-size outlier population. Bounding the cap requires a **tail / large-deviation
   / extremal** argument over the singular series `L(S)`, NOT a variance bound — precisely the `SL(4)`
   quadruple-resonance extremality (HYP-2823's consec-max, the AP as the `L`-minimizer).
3. **Target the `L(S)` minimum, not its variance:** the floor lives in `mean ± √var`; the cap lives at
   `min_S L(S)` (the tail edge). The proof of LRC(14) is an extremal `min_S L(S) > 0` statement, and the
   metagraph (no tail) is silent there — the right rehearsal for the tail is the **resonance-arity /
   additive-relation extremality** (the AP maximizes `S₄`), not the variance.

## Status
- **PROVED (opus):** `CV²(H) ~ 2/n` (`n·CV²→2`); dispersion = succession average
  `E_σ[2^{asc}·1[no desc]]`; `H` fully concentrates.
- **Computed:** LRC `CV²_S(D/q)` plateaus at ≈0.083 (persistent tight tail), NOT → 0.
- **Diagnosis:** metagraph (full concentration) = the generic floor model; the LRC tight tail = the cap,
  beyond any homogeneous second-moment bound; the plateau = the floor↔cap gap = the `L(S)` spread.
- **Refined target:** LRC(14) = `min_S L(S) > 0` (a tail/extremal statement), not a variance bound; the
  rehearsal for it is the additive-relation (`S₄`) extremality, not `Var(H)`.

Related: the metagraph-variance-closed-form reflection (`G(n)`, the succession count), the Siegel–Rogers
moment hierarchy (floor=SL2 / cap=SL4), the Han-Lee generic-floor computation, the razor-thin (measure
vs peak), HYP-2823 (consec variance extremality), THM-501 (singular series `L(S)`), Riordan successions
([A002464](https://oeis.org/A002464)), OPEN-Q-108.
