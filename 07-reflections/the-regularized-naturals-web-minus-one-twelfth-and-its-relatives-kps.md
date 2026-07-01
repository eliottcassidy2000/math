# The regularized-naturals web: −1/12 and its relatives across the project

*kind-pasteur-2026-06-30. Creative exploration of `1+2+3+⋯ = ζ(−1) = −1/12` and the wider family of regularized infinite sums/products/fractions/roots, keeping only the ties that are **structural** (verified, `lrc14_regularized_naturals_zoo_kps.py`). Three of them are genuinely new; one unifies the project's deepest recurring fact.*

> **Concurrent convergence (2026-06-30).** klein/mac-mini **HYP-3768** and opus **HYP-3770** landed the `−1/12 = ζ(−1)` node as `s(n,Φ₆) → −1/12` (the covering-min `E₂`/Dedekind margin), and opus **HYP-3773** ("the sum of natural numbers IS the covering-min": `Φ₆ = 2·Σspeeds + 1`, `killer = 2Σ ≡ −1`, and — matching **C1** here — the `ζ(−2)=0`/squares/hexagonal even node plus the η-weight anomaly = my Casimir→η²⁴ **C4**) covers most of this web on the covering-min object. They pushed first; I defer. What remains distinct here: **C2** (the two Dedekind reciprocity *constants* are the plain & alternating regularized naturals — `η(−1)`, `−ζ(−1)`) and **C6** (the `χ₇`-twisted even-vanishing `B_{2,χ₇}=0` over ℚ(√−7)).

## −1/12 is one node, not the whole story

`−1/12` felt central (the resonance-killing limit, `B₂/2`, the runner-12 coincidence). But it is one value of a whole web, and reading the web shows which appearances are deep.

## C1 — the even-annihilation IS `ζ(−2n) = 0` (the unifying one)

`ζ(−k) = −B_{k+1}/(k+1)` gives the tower `−1/2, −1/12, 0, 1/120, 0, −1/252, 0, …` (and `ζ(−11)=691/32760`, the famous 691). The pattern: **`ζ(negative even) = 0`; `ζ(negative odd) = Bernoulli ≠ 0`.** That parity is the project's single deepest recurring fact, seen four ways — all "even index → 0":

- **Burnside:** `Fix(σ) = 0` whenever `σ` has an even cycle (the A000568 engine).
- **Homology:** `β₂ = 0` for all tournaments.
- **Regularized power sums:** `ζ(−2n) = 0` (the trivial zeros).
- **Twisted (C6 below):** `B_{2,χ₇} = 0`.

and the **odd** side survives everywhere: `ζ(−odd) = Bernoulli`, odd cycles (OCF / Rédei / the odd Hamiltonian paths), `B_{1,χ₇}` (a class number). **The odd/even duality that runs through the entire project is the parity of Bernoulli/`L`-values at negative integers**, and `−1/12 = ζ(−1)` is simply its first surviving (odd) term.

## C2 — the Dedekind reciprocity constants ARE the two regularized natural sums (new)

The barrier residual is a Dedekind sum (THM-563), governed by reciprocity
`s(h,k)+s(k,h) = −1/4 + (1/12)(h/k+k/h+1/hk)`. Its two universal constants are exactly the two ζ-regularizations of the naturals (verified):

> `1/12 = −ζ(−1) = −(1+2+3+⋯)` — the **plain** (unsigned) sum;
> `1/4 = η(−1) = 1−2+3−4+⋯` — the **alternating** (signed) sum.

So the residual's reciprocity is a *regularized-naturals machine*, with the **signed/unsigned split built in as its two constants** — the same split that is the whole difficulty of LRC (the signed cancellation positivity can't reach). The `−1/4` is the alternating/signed constant; the `1/12` the plain/unsigned one. (And `η(0)=1−1+1−⋯=1/2` is Grandi, `ζ(0)=1+1+1+⋯=−1/2=B₁`.)

## C4 — the Casimir road upgrades the "two twelves" coincidence (real bridge)

The earlier `two-twelves` reflection called the runner-12 = weight-12 alignment a coincidence. There is, however, a *real* bridge from `−1/12` to weight 12, and it is not the `n−2` numerology — it is the **string Casimir**:

> `−½·ζ(−1) = 1/24` = the exponent in `η(τ) = q^{1/24}∏(1−qⁿ)`; and `η²⁴ = Δ`, the **weight-12** cusp form; `24 = 2·12` = the GW doubling.

So `−1/12 → Casimir 1/24 → η²⁴ → weight 12` is a genuine classical chain (ground-state energy of the bosonic string). This is the structural skeleton the modular link needs; what is still missing is a theorem tying `η²⁴`/weight-12 to the LRC(14) or tournament objects (the "`η²⁴` = code discriminant" memory is an annotation, not a proof). But the road from the regularized natural sum to weight 12 is real.

## C5 — the product-side dual: `∏ n = √(2π)`

The functional determinant `∏_{n≥1} n := exp(−ζ′(0)) = √(2π)` (since `ζ′(0) = −½ln 2π`). So the regularized *product* of all naturals is `√(2π)` — the Gamma/Stirling `e, √(2π)` of the triangle foundation, **dual to the −1/12 of the sum**. Sum ↔ product, `−1/12 ↔ √(2π)`, two faces of the same Gamma/ζ analytic continuation.

## C6 — the apex-7 twisted −1/12, and why it vanishes (new, ties to ℚ(√−7))

The LRC(14) obstruction lives over `ℚ(√−7)` (the imaginary Gauss sum `i√7`, because `7 ≡ 3 mod 4`). The apex analogue of `ζ(−1)=−1/12` is `L(1−n,χ₇) = −B_{n,χ₇}/n` via the generalized Bernoulli. Computed:

> `χ₇` is **odd** (`χ₇(−1)=−1`, because `7≡3 mod 4`), so `B_{2,χ₇} = 0` — the **twisted even-vanishing** — hence `L(−1,χ₇)=0` (a trivial zero). The surviving datum is the **odd**-index `B_{1,χ₇} = −1`, i.e. `L(0,χ₇)=1`, carrying the class number `h(ℚ(√−7))=1`.

The punchline: **the apex oddness `7≡3 mod 4` that makes the Gauss sum imaginary (`i√7`) is the *same* parity that zeroes the even twisted Bernoulli `B_{2,χ₇}`.** The even-annihilation of C1 repeats over the signed-obstruction field. So if the residual is reorganized over `ℚ(√−7)` (the flagged route to make the signed obstruction real), the relevant regularized constant is not the plain `−1/12` but the `χ₇`-twisted one — whose *even* part vanishes and whose *odd* part is the class number.

## The methodological payoff (why this isn't just numerology)

`−1/12` is the archetype of "**a divergent sum made finite by analytic continuation**." The LRC's central difficulty is exactly that: the singular series' `|T|≥3` tail diverges *absolutely* (`A₃=∞`) but the *signed* sum is conditionally summable (THM-504). Zeta/Abel/Ramanujan regularization is the theory of assigning finite values to such sums — and the repo already uses its working cousin (Abel/Dirichlet summation, HYP-3129 signed SPEC). So the `−1/12` is the emblem of the **regularization method the signed content needs**, not merely a number that happens to appear.

## Net — the web, sorted

- **Structural & verified:** `ζ(−even)=0` = the project's even-annihilation (Burnside/β₂/twisted); Dedekind reciprocity constants = the plain & alternating regularized naturals (`−ζ(−1)`, `η(−1)`); Casimir `1/24 → η²⁴` = a real road to weight 12; `∏n=√(2π)` the product dual; `B_{2,χ₇}=0` the apex-7 twisted even-vanishing.
- **Creative directions worth pursuing:** (i) regularize the residual/singular-series over `ℚ(√−7)` with the `χ₇`-Bernoulli (even vanishes, odd = class number); (ii) chase the `η²⁴`/weight-12 ↔ LRC(14) theorem the Casimir road makes plausible; (iii) the sum/product duality `−1/12 ↔ √(2π)` as two faces of the ζ continuation.
- **Honest:** the modular (weight-12/`η²⁴`) ↔ LRC theorem is still absent; the Casimir chain makes it structurally plausible, no more.

— Related: `two-twelves-...-kps.md`, `the-resonance-killing-game-and-the-zeta-duality-of-the-lonely-runner.md`, `zeta2-governs-the-lonely-runner-floor.md`, `the-barrier-residual-is-a-dedekind-sum-...-kps.md`, `everything-is-the-triangle.md` (the √2,π,e,γ pantheon), THM-563, HYP-3129 (Abel/signed SPEC), HYP-3252/3254 (ℚ(√−7)). Artifact: `04-computation/lrc14_regularized_naturals_zoo_kps.py`.
