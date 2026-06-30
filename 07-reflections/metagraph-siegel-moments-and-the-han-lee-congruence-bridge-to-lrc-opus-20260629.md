# The metagraph has Siegel moment formulas; Han–Lee's dim-2 congruence Siegel transform is the SL(2)/ζ(2) tool for the LRC floor with the covering congruence — and both are "moment formulas for arithmetic-orbit counts"

*opus-2026-06-29. Owner: extend the tournament-metagraph work and define new small useful objects;
read arXiv:2507.05905 (Han–Lee, "Moment formulas of Siegel transforms with congruence conditions in
dimension 2"); think about how metagraph ideas help LRC proofs. The three meet at one frame: Siegel
(first/second) moment formulas for counting orbits of an arithmetic group, with a congruence subgroup =
the covering constraint.*

## The paper (arXiv:2507.05905, Han–Lee)
First and second **moment formulas for the Siegel transform with congruence conditions in dimension 2**:
counting primitive `(p,q)∈Z²` with `(p,q)≡(p₀,q₀) (mod N)`. Density scales with `1/φ(N)`; the counting
is governed by the index `[SL(2,Z):Γ₀(N)] = N∏_{p|N}(1+1/p)`; `ζ(2)` is the normalization. They prove a
congruence **Schmidt random-counting** analog and a **quantitative Khintchine** theorem for `p/q`
approximations under the congruence. This is exactly the `SL(2)/ζ(2)` second-moment machinery the
moment-hierarchy reflection placed at the LRC **floor** — now WITH a congruence subgroup.

## New small useful metagraph objects (defined + verified, n≤6)
The tournament metagraph turns out to carry the SAME Siegel structure:
1. **The metagraph zeta** `ζ_G(s) := Σ_{iso-classes C} H(C)^{-s}`. Special values: `ζ_G(0)=V_n=A000568`,
   `ζ_G(−1)=Σ_C H = 12, 92, 1224` (n=4,5,6 — a new sequence), `ζ_G(1)=Σ1/H = 1.867, 2.923, 5.380`.
2. **The metagraph Siegel MASS FORMULA** (reframe of LEM-003's tiling identity): `Σ_C H(C)/|Aut(C)| =
   2^{C(n−1,2)}` (verified 8, 64, 1024). This is a **mass formula**: `H/|Aut|` is the class mass, the
   tiling-count is the total — the tournament analog of the Siegel `Σ 1/|Aut| = vol` mass formula.
3. **The metagraph moment sequence** `μ_k(n) := E_T[H^k]` over labeled tournaments. **`μ₁ = n!/2^{n−1}`
   EXACTLY** (each of the `n!` orderings is a Hamiltonian path with probability `2^{−(n−1)}` — the
   first-moment/Siegel formula). And:
   > **`μ₂ = E[H²] = Σ_{π,π'} 2^{−|arcs(π)∪arcs(π')|}`** — a pair-correlation over **ordering pairs**,
   > graded by how many consecutive-arcs they share. (`E[H²] = 3, 12, 1185/16, 1305/2`; `Var(H) = 3/4,
   > 3, 285/16, 585/4`.) This is **structurally the paper's dim-2 Siegel second moment** (a sum over
   > pairs graded by coincidence) and the LRC's `S₂` pair moment — the same object three times.
4. **The forbidden value `H=7`** is absent from every spectrum (n=4,5,6) — the known forbidden-`7` state
   (THM-572/HYP-2908), the apex-7 face inside the H-spectrum.
5. **The congruence metagraph** `G_n(N)`: iso-classes refined by `H mod N` (or score-sequence `mod N`).
   This is the **`Γ₀(N)` analog**: the covering congruence `mod N` refines the orbit count exactly as
   `Γ₀(N)` refines `SL(2,Z)`. (To compute next.)

## The unifying frame: moment formulas for arithmetic-orbit counts
| | group action | orbit count | mass / 1st moment | 2nd moment | congruence |
|---|---|---|---|---|---|
| **metagraph** | `S_n ↷` tournaments | `V_n=A000568` | `Σ H/|Aut|=2^{C(n−1,2)}`, `E[H]=n!/2^{n−1}` | `E[H²]` = ordering-pair correlation | `H mod N` (`G_n(N)`) |
| **LRC floor** | `Γ₀(N) ↷` Farey/modular | witnesses `p/q` | `c_q` mean (union bound, vacuous) | `S₂` = `ζ(2)` Farey pair density | covering = `mult-of-q` |
| **Han–Lee** | `Γ₀(N) ⊂ SL(2,Z)` | primitive `(p,q)` cong. | `1/φ(N)` density, `N∏(1+1/p)` | the paper's 2nd moment | `(p,q)≡(p₀,q₀) mod N` |

> **All three are Siegel first/second moment formulas for an arithmetic-group orbit count, and the
> covering/tight constraint is a CONGRUENCE SUBGROUP.** The metagraph computes them by Burnside/cycle-
> index (klein THM-587, `P_n(±1)`); the LRC and the paper compute them by the geometry of numbers on
> `SL(2)`. The `R`-even/`R`-odd split of the metagraph (`±1` of the complement) is the principal-vs-
> congruence (Eisenstein-vs-cusp) split of the Siegel transform.

## How metagraph ideas help the LRC proof (concrete)
1. **Han–Lee is the floor tool, with the covering built in.** The LRC floor `c_q ≥ 1/(2ζ(2))` (HYP-2856)
   is the `SL(2)` Siegel pair density; the paper supplies the **congruence** second moment, so the
   covering constraint (`q ≡ 0 mod the speeds`) enters as `Γ₀(N)` rather than being bolted on. This is a
   cleaner, rigorous floor proof than the elementary totient sum — and it generalizes off the specific
   speed set.
2. **The metagraph mass formula models the LRC mass formula.** `Σ_C H/|Aut| = 2^{C(n−1,2)}` is a
   *closed-form* first moment; its proof (orbit–stabilizer on tilings) is the template for a closed-form
   LRC first moment under congruence — the thing the union bound lacks.
3. **The second-moment correlation is the shared engine.** `μ₂` (ordering pairs), the paper's Siegel 2nd
   moment (lattice-point pairs), and `S₂` (resonance pairs) are one structure; a variance bound proven on
   the metagraph (clean, finite) is a rehearsal for the LRC `S₂` variance bound (HYP-2823's `Var(N)`).
4. **The congruence metagraph `G_n(N)` is the missing dictionary entry.** Computing how `H mod N`
   distributes over classes (the `Γ₀(N)` coset structure) gives the tournament-side image of the
   covering congruence — and may expose which residues are forbidden (the `H=7` obstruction generalized),
   i.e. the tournament avatar of "covering forces `M>1/14`."

## Status
- **Verified:** `E[H]=n!/2^{n−1}`; mass formula `Σ H/|Aut|=2^{C(n−1,2)}`; `E[H²]`, `Var(H)` (n≤6);
  `ζ_G(−1)=12,92,1224`; `H=7` absent.
- **New (opus):** the metagraph zeta `ζ_G`; the Siegel mass-formula reframe; the metagraph moment
  sequence `μ_k` with `μ₁=n!/2^{n−1}` and `μ₂` = ordering-pair correlation; the congruence metagraph
  `G_n(N)`; the unifying "Siegel moments for arithmetic-orbit counts" frame; Han–Lee as the covering-
  congruence floor tool.
- **Open / next:** compute `G_n(N)` (`H mod N` distribution); apply Han–Lee's congruence 2nd moment to
  the LRC covering floor explicitly; catalog `Σ_C H` and `μ_k` sequences (OEIS).

Related: the Siegel–Rogers moment-hierarchy reflection (floor=SL2/ζ2, Littlewood=SL3, cap=SL4),
zeta-duality+Littlewood, cuts-as-Farey+hyperoctahedral-metagraph; THM-501 (singular series), HYP-2856
(ζ(2) floor), HYP-2823 (variance extremality), klein THM-587 (`P_n(±1)`), LEM-003 (tilings·|Aut|=H),
euler-product-and-metagraph, metagraph-as-transfer-chain, eigenvalues-of-the-merged-metagraph, THM-572
(forbidden H=7), OPEN-Q-108.
