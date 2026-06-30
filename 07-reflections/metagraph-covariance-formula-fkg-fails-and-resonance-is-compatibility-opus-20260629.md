# The metagraph Hamiltonian-path covariance formula (proved), FKG fails for both H and the LRC, and "resonance = compatibility = positive covariance" — refining HYP-2823's variance extremality, localized at the high-H end

*opus-2026-06-29. Owner: prove the FKG/variance inequality on the finite metagraph as a rehearsal for
the LRC `S₂`/HYP-2823, and work it back and forth with idea 6 (H as Lyapunov / the H-gradient) when
stuck. The rehearsal proved a clean covariance formula, revealed that FKG FAILS, and back-transferred to
refine the LRC variance mechanism; idea 6 localized where the variance lives.*

## The proved result (metagraph, verified n=3..6)
`H = Σ_π X_π`, `X_π = 1[π is a Hamiltonian path]` over the uniform random tournament — a sum of
**dependent** indicators, exactly HYP-2823's shape but finite. For orderings `π, π'` of `[n]`:
> **`E[X_π] = 2^{−(n−1)}`**, and
> **`Cov(X_π, X_π') = 2^{−2(n−1)}(2^{c} − 1)`** if `π, π'` are **arc-compatible**, `= −2^{−2(n−1)}` if
> they **conflict**, where `c =` # common directed consecutive-arcs and *conflict* = some unordered
> pair is a consecutive arc in both with opposite orientation.
> Hence **`Var(H) = 2^{−2(n−1)}[ Σ_{compatible(π,π')} 2^{c} − (n!)² ]`.**

*Proof.* `E[X_π]=P(\text{all }n−1\text{ consecutive arcs of }π\text{ present})=2^{−(n−1)}`.
`E[X_πX_{π'}]=P(\text{both})=2^{−|arcs(π)∪arcs(π')|}` if the union is consistent, else `0`;
`|arcs(π)∪arcs(π')|=2(n−1)−c`, so `E[X_πX_{π'}]=2^{c−2(n−1)}` (compatible) or `0` (conflict);
subtract `E[X_π]E[X_{π'}]=2^{−2(n−1)}`. ∎ (Matches direct `Var(H)=3/4,3,285/16,585/4`.)

## FKG FAILS — and it fails for the LRC too (the back-transfer)
The covariance is **mixed-sign**: positive on compatible pairs (shared arcs), negative on conflicts. And
the conflicts *outnumber* the compatibles (n=6: 295920 conflict vs 222480 compatible). **So `H` is NOT a
positively-associated indicator sum — FKG does not apply.** Back-transferring the test to the LRC: the
danger events `1[‖s_i t‖<1/14]` over `t` are **also mixed-sign** (AP: 49 positive, 28 negative
covariances), and the sign is arithmetic:
> **Resonant pairs (`s_i | s_j`, large `gcd`) have POSITIVE covariance** — the danger arcs nest/align,
> the "compatible" pairs: `(1,2),(2,4),(3,6),(5,10),(6,12),(4,8)` top the list. **Coprime/anti-resonant
> pairs have NEGATIVE covariance** — danger arcs avoid each other, the "conflicts": `(1,13),(1,11),(2,11)`.

So the metagraph's compatible/conflict dichotomy IS the LRC's resonant/coprime dichotomy. **Resonance =
compatibility = positive covariance**; coprimality = conflict = negative covariance — the same
divisibility law as the cut/resonance-killing/Farey structure, now seen in the *second moment*.

## This REFINES HYP-2823 (not pure FKG, but a net balance)
HYP-2823 reads "consec maximizes `Var(N)` because clustered runners empty sectors together (`Cov>0`)" —
a positive-association story. The rehearsal corrects it: **FKG fails; the variance is a NET balance of
resonant (`+`) over coprime (`−`) covariances.** The AP/staircase wins not by all-positive association
but by being **divisibility-rich** — it packs the most resonant pairs (`(1,2),(2,4),(2,6),(3,6),…`, the
positive-covariance "compatible" pairs) while the forced coprimality of spread sets adds negative terms.
This is the second-moment face of the additive-relation-richness (the `SL(4)` resonance count) from the
moment-hierarchy reflection: **consec maximizes the variance because it maximizes the resonant
(divisible) pair-overlaps.**

## Idea 6 worked back in: the variance lives at the high-H end
Where is `Var(H)` carried? By H-value (n=6): **`H≥23` (top half of the range) holds 89.9% of `Σ H²`** —
the variance is concentrated at the **high-H (regular/Paley) end** of the gradient. This is idea 1 ∩
idea 6: the **compatible (positive-covariance) HP-pairs cluster in the high-H tournaments** (a regular
tournament packs many mutually-compatible Hamiltonian paths), so `H` as a Lyapunov potential *localizes*
the variance — it is an excursion phenomenon at the top of the H-gradient, not spread uniformly. The
"almost-DAG" structure (H increasing along edges, MISTAKE-035) means the second moment is an upper-tail
object: to bound `Var(H)` (or the LRC variance), look only near the regular end.

## Two clean takeaways for the proof program
1. **Drop FKG; use the signed covariance.** Both `H` and the LRC danger sum have negative covariances; a
   positive-association (FKG/Holley) argument is unavailable. The right object is the *signed* sum
   `Σ_{compatible}2^c − (conflicts)`, with sign = resonance (`gcd`) — a number-theoretic, not a
   monotone-lattice, statement.
2. **The variance is an upper-tail (high-H / regular) phenomenon.** For the LRC, the analog is that the
   miss-count variance is carried by the *most resonant* configurations (the AP/regular end), so the
   extremality `consec = argmax Var` is a statement about the top of the resonance gradient — exactly
   where the `SL(4)` quadruple-resonance density (and the tight cap) lives.

## Status
- **PROVED (opus):** the metagraph HP pair-covariance formula and `Var(H)` reduction (verified n=3..6).
- **Computed:** FKG fails for `H` (conflicts outnumber compatibles) AND for the LRC danger sum
  (49+/28−); resonant pairs = positive covariance, coprime = negative; variance carried by `H≥23`
  (89.9%).
- **Refined target:** HYP-2823's consec-max-variance is a NET resonant-over-coprime balance (not FKG),
  the second-moment face of additive-relation-richness; it is an upper-tail (regular-end) extremality.
- **Open:** a closed form for `Σ_{compatible}2^c` (a directed-linear-forest sum, messy); the signed-
  covariance proof of consec-extremality.

Related: HYP-2823 (variance extremality), the Siegel–Rogers moment-hierarchy + Han–Lee reflections,
metagraph-siegel-moments (E[H], E[H²]), cuts-as-Farey-geodesics (the resonance=divisibility law),
THM-534 (moment-LP), MISTAKE-035 (H-gradient not strict DAG), eigenvalues-of-the-merged-metagraph
(H = 2nd eigenvector), OPEN-Q-108.
