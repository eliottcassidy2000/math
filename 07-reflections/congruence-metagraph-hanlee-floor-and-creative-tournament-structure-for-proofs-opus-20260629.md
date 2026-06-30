# The congruence metagraph G_n(N), the Han–Lee second-moment LRC floor (generic = SL(2)/EKL level), and seven creative ways tournament/metagraph structure can drive proofs

*opus-2026-06-29. Owner: build the congruence metagraph and set up the Han–Lee congruence-Siegel second
moment for the LRC covering floor; then abstractly brainstorm other creative uses of tournament/metagraph
structure for proofs and insight.*

## 1. The congruence metagraph G_n(N) (built, n=5,6)
Refine iso-classes by `H mod N`, tracking class-count and **mass** `Σ H/|Aut|` (total `= 2^{C(n−1,2)}`):
- **mod `2^k` = the OCF 2-adic stratification (THM-466).** `H` is odd, and `H mod 4` splits by the
  parity of `α₁` (the 3-cycle count): **`H≡1 mod 4` carries more mass than `H≡3`** (656 vs 368 at n=6;
  42 vs 22 at n=5). The 2-adic congruence metagraph IS the odd-cycle-collection filtration — clean.
- **`H≡0 mod 7` is EMPTY at n≤6** (no class has `7|H`) — the forbidden-7 face inside the spectrum. (NOT a
  theorem for all n: at n=7 Paley has `H=189=7·27`. The forbidden object is the *value* `H=7`, THM-572,
  not the residue.)
- **Mass is biased, not equidistributed mod N** (mod 3 at n=6: `256, 448, 320`). The bias is the OCF
  signature — the tournament avatar of the covering set's non-uniform residues.

This is the `Γ₀(N)`-analog: the covering congruence `mod N` refines the orbit count exactly as `Γ₀(N)`
refines `SL(2,Z)`; the metagraph computes the refinement combinatorially (the 2-adic OCF digits).

## 2. Han–Lee congruence second moment → the LRC floor (generic; the honest reach)
Set `D(q,S)=#{a∈Z/q : ‖va/q‖≥1/14 ∀v∈S}` (THM-501 witness count). Han–Lee's congruence second moment
is the variance of `D/q` over the ensemble of covering sets. Computed (q=1009 prime, 60 covering sets):
- **1st moment** ≈ independence `(6/7)^{13}=0.135`; **sample mean** `D/q ≈ 0.092`.
- **2nd moment / variance** `≈ 0.0007` — **strong concentration**.
- **min `D/q ≈ 0.006 > 0`** (the `7/89` near-tight set) — even the tightest keeps `~6` lonely residues
  per `1000`.
> **The floor holds for GENERIC covering sets by concentration — the `SL(2)`/EKL level.** Han–Lee gives
> the *quantitative Khintchine with congruence* = "almost every covering direction is lonely, with the
> exceptional set small." This is the rigorous partial result; it does **not** reach the tight cap (ALL
> sets, no exceptions) — that is the `SL(4)` quadruple-resonance gap, "one dimension past Littlewood."
> The second moment is exactly as far as it reaches, and now we can *say* so quantitatively.

## 3. Seven creative uses of tournament/metagraph structure for proofs
1. **The metagraph as a FINITE REHEARSAL for the moment method.** The metagraph's moments are exact and
   computable (`E[H]=n!/2^{n−1}`, mass `Σ H/|Aut|=2^{C(n−1,2)}`, `E[H²]`=ordering-pair correlation). It
   is the same Siegel–Rogers machinery as LRC/Littlewood but FINITE. Strategy: prove a hard moment/
   correlation inequality (FKG, variance extremality, a covariance bound) on the metagraph FIRST, where
   every step is checkable, then transfer the *technique* to the infinite LRC. The metagraph is the
   training ground for `S₂`/`S₄`.
2. **FORCED EXISTENCE via Ky Fan / Borsuk–Ulam.** The complement `R` acts freely; `P_n(−1)=SC>0` is a
   forced-nonzero alternating sum (a Ky Fan certificate that a palindromic locus exists). Abstractly:
   **any tournament property whose `R`-signed Euler characteristic is nonzero is FORCED to be realized.**
   The LRC lonely-point existence should be such a forced count over the comparator movie — define the
   `R`-signed count of the comparator's source-events whose non-vanishing forces a lonely time. Existence
   without construction.
3. **LOCAL–GLOBAL via the free monoid.** Iso-classes = free monoid on SC primes; `H` multiplicative.
   Any multiplicative invariant obeys a local–global principle: **prove it on the SC primes
   (irreducibles), multiply for the rest.** Bounds that are multiplicative + hold on primes hold
   universally — the max-H=single-prime result is one instance; the OCF/`H` bounds are candidates.
4. **The SIEGEL MASS as the correct averaging measure.** Naive "uniform over classes" is the WRONG
   probability measure; the correct one weights each class by its mass `H/|Aut|` (so `Σ=2^{C(n−1,2)}`).
   Probabilistic tournament arguments should average against the **mass measure**, which removes the
   `|Aut|`-bias — the same fix the Siegel measure makes in the geometry of numbers. This likely sharpens
   the OCF-average and concentration estimates.
5. **The metagraph zeta `ζ_G(s)=Σ_C H^{−s}` and its symmetry.** `R` preserves `H`, so `R` is a symmetry
   of `ζ_G`; the `R`-even/`R`-odd split gives `ζ_G=ζ_G^{+}+ζ_G^{−}`. Conjecture a functional equation
   relating `ζ_G(s)` to `ζ_G(c−s)` via the transitive(`H=1`)↔regular(`H=max`) H-gradient duality, with
   the "zeros" read off the metagraph adjacency spectrum (where `H` is the 2nd eigenvector). A
   tournament "explicit formula."
6. **`H` as a LYAPUNOV POTENTIAL.** `H` increases along most metagraph edges (the H-gradient); the
   metagraph is almost a DAG (rare level/back-edges, MISTAKE-035). Use `H` as a potential for
   termination/monotonicity arguments; the back-edges and level-edges are precisely where the hard
   structure (the genus, the obstruction) concentrates — a localizer for difficulty.
7. **The EVEN-GRAPH DUAL `E_n` for parity/homology.** `E_n` is the `Z_2`-dual of `G_n` (tiling →
   cycle-space XOR → even graph). Parity obstructions (Rédei oddness, the `Z_2` cap residual) are
   cleaner on `E_n` (homology of the cycle space). Transfer a tournament parity question to an
   even-graph homology computation, where `Z_2`-linear algebra replaces case analysis.

## Status
- **Built:** `G_n(N)` (mod-`2^k`=OCF, `H≡0 mod 7` empty at n≤6, biased mass); Han–Lee 2nd-moment LRC
  floor (generic concentration, `SL(2)`/EKL; tight cap = `SL(4)` gap).
- **New (opus):** the congruence metagraph as the `Γ₀(N)` analog; the quantitative "generic floor"
  statement; seven proof-strategy templates (finite rehearsal, forced existence, local–global, mass
  measure, metagraph zeta, Lyapunov H, even-graph dual).
- **Next:** push idea 1 (prove a metagraph correlation inequality, transfer to `S₂`); idea 2 (the LRC
  Ky-Fan existence count); idea 4 (re-derive the OCF average under the mass measure).

Related: the Han–Lee Siegel-moment bridge reflection; the Siegel–Rogers moment hierarchy (SL2/SL3/SL4);
THM-466 (2-adic H = OCF), THM-572 (forbidden H=7), klein THM-587 (`P_n(±1)`), LEM-003 (mass formula),
THM-501 (singular series), HYP-2823 (variance extremality), euler-product-and-metagraph,
metagraph-as-transfer-chain, even-graphs-as-first-class, OPEN-Q-108.
