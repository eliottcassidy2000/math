# The covering is a congruence subgroup

*mac-mini-2026-06-29-S20. Extending the owner's synthesis: Han–Lee's congruence Siegel moments are the LRC floor with the covering built in as Γ₀(N), the metagraph mass formula is the missing first moment, and G_n(N) is the tournament image of the congruence. New: HYP-3553.*

## The move the whole program turns on

For sessions the LRC floor has been a sum. The surviving lonely mass is `Σ_b φ(b)δ_b`, a totient-weighted Farey sum whose density is the `ζ(2)` Euler product (HYP-3550). It is correct and it is positive, but it has a flaw the proof keeps tripping on: it is bolted on. The covering constraint — which resonances are killed, which residues survive — enters the sum from outside, set by set, by inclusion–exclusion over the particular speeds. That set-dependence is exactly why the uniform floor (the gatekeeper, OPEN-Q-108) has stayed open: you have to control the overlaps for *every* covering set.

Han–Lee's paper makes a different move, and it is the right one. Counting primitive `(p,q) ≡ (p₀,q₀) (mod N)` is not "count primitive pairs, then restrict" — it is an average over a **congruence subgroup** `Γ₀(N) ⊂ SL(2,ℤ)`, with the restriction living inside the group, governed by the index `[SL(2,ℤ):Γ₀(N)] = ψ(N) = N∏_{p\mid N}(1+1/p)` and the density `1/φ(N)`, all normalized by `ζ(2)`. The covering is not a side condition on a sum; it *is* the subgroup. And the consequence is set-independence: the moment depends on `N` — the covering modulus — through `φ(N), ψ(N), J₂(N)`, and on nothing else. For `N=14`: `φ=6`, `ψ(14)=24`, `J₂(14)=144`. Three numbers, no speed set. The covering floor `c_q ≥ 1/(2ζ(2))` becomes a `Γ₀(N)` congruence second moment, and the thing that blocked the uniform bound — the per-set overlaps — is absorbed into the index of a subgroup.

That is the first bullet, and it reframes the gatekeeper. The open piece of THM-579 was "bound `CV(N_R)²` uniformly over all covering sets." Through `Γ₀(N)` it becomes "evaluate the congruence second moment at modulus `N`" — a single arithmetic computation, not a quantifier over sets.

## The mass formula is the first moment the union bound never had

The union bound is an inequality with no equality underneath it. The metagraph has the equality: Burnside's mass formula, `#classes = \frac{1}{n!}\sum_\sigma \mathrm{Fix}(\sigma)`, an exact count of orbits as an average over the group. That is a *first moment*, and the Siegel mass formula (Han–Lee's first moment, with congruence) is its continuous twin: the exact expected number of primitive lattice points in a region, in a congruence class. The LRC floor wants exactly this — the expected number of surviving lonely points — and has only ever had `≥`. Give it the mass formula and the union bound becomes a first-moment equality, the second moment gives concentration, and Chebyshev does the rest. The metagraph isn't an analogy here; it is the finite worked example of the formula the floor is missing.

## Same engine, three voices

The second moment is a sum over pairs, always. On the metagraph it is `μ₂`, the count of ordered pairs of arcs — the two-point correlation, finite and exactly computable, and it does not blow up (`CV(H) ≈ 0.5`–`0.6`). On the lattice it is `∫\hat f^2`, the Rogers–Schmidt pair correlation of primitive points. On the runner it is `Var(N_R)`, the variance of the 14-sheet count that THM-579's coefficient of variation is built from. These are one object spoken in three registers, and the metagraph register is the one where you can hold the whole thing in your hand and watch the variance stay bounded. It is the rehearsal; the congruence Siegel formula is the performance; the LRC is the audience that needs to hear it.

## G_n(N): the tournament at level N

If `Γ₀(N)` is the covering as a subgroup, what is the tournament image? It is `G_n(N)`, the **level-N congruence metagraph** — the iso classes of `Z/N`-circulant tournaments, a tournament carrying a marked `Z/N` structure exactly as a point of `X_0(N)` carries a marked cyclic `N`-subgroup. A circulant tournament on `N` odd vertices is a sign-choice over the `(N-1)/2` antipodal pairs of residues; the multiplier group `(Z/N)^*` acts; the mass is the dihedral Burnside count — `1,1,2,4,4,6` for `N=3,5,7,9,11,13`. A tiny, finite-index, structured shadow of the full metagraph, just as `Γ₀(N)` is finite-index in `SL(2,ℤ)`. And it has a distinguished point: the **Paley class**, the quadratic-residue connection set, the tournament's CM point — the cusp of `Γ₀(N)`, where kps's OCR-genus law already lives (the OCR denominators factor into primes whose modular curves `X_0(p)` have genus 0 or 1). The missing dictionary entry was never missing; it was the modular tournament all along, read at level `N`.

## Where this points

Three creative extensions fall out once the dictionary is in place. Vertex-addition is the modular `T`; as a raising operator `G_n \to G_{n+1}` it is **Hecke-like**, and diagonalizing it would give H-recursions that are the tournament analog of Hecke eigenforms — with the `X_0(p)`-genus predicting which eigenvalues are rational. That same **genus law is a free closed-form predictor**: the LRC floor constant at covering modulus `N` should be "modular/rational" exactly when `X_0(N)` has small genus, telling you in advance which `N` admit a clean first moment. And the metagraph's **spectral gap** (HYP-3552's moments) is the expander statement behind the second-moment method — the formal content of "most tournaments, most covering sets, behave like the mean."

The covering was never a constraint to be enforced. It was a subgroup to be quotiented by. Han–Lee built the door; the metagraph is the room on the other side; the LRC is what the room was for.

See [[the-metagraph-is-a-finite-siegel-transform]] (HYP-3552, the no-congruence version), [[the-modular-tournament]] (tournaments as PSL(2,ℤ), the X_0(p)-genus OCR law), [[the-cut-prefix-is-a-continued-fraction-times-an-euler-product]] (HYP-3550, the totient-sum floor this replaces), [[zeta2-governs-the-lonely-runner-floor]]. New: HYP-3553. External: arXiv:2507.05905.
