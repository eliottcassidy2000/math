# Wild structures: ear decompositions, Bravais-14, Savitch, and the φ⁴ proof lever

*kind-pasteur-2026-06-27-S31ah. The owner asked to be bold and even make wild guesses about recurring
structures, weaving in Savitch, the Bravais lattice, Lee-Yang/φ⁴, and the ear-decomposition theorems,
toward "what maximizes LRC values relating to tournaments." Here are the leads, with the verified
skeleton facts and the bold extrapolations marked.*

## 0. The verified skeleton facts (the hooks)
- **#ears `= C(n,2)−n+1 = C(n-1,2) = #tiles`** EXACTLY for all `n` (verified). A directed ear decomposition
  of a tournament has exactly as many ears as the tiling model has free tiles = the cycle-space dimension.
- **Bravais-14 `= 7 crystal systems × centering = 2·7`**; the count of Bravais lattices by dimension is
  **`1, 5, 14, 64`** (dim 1,2,3,4) — **dim-3 lands on 14 = LRC(14)**, and `5 = ` the last classically-easy
  LRC case (`n=5`). Both factor through the apex `2·7`.
- **Factor-critical ⟺ odd ear decomposition**; `K_n` is factor-critical iff `n` is **odd**.

## 1. EAR DECOMPOSITION = the cycle-space, and the ODD ears = Ω (a new H-recursion lead)
Since `#ears = #tiles = dim(cycle space)`, the ear decomposition is the project's **wiggly/tile/cycle-space
structure** under a new name (CLAUDE.md: "wiggly arcs = cycle-space generators"). The NEW content is the
**ODD** ear decomposition (⟺ factor-critical, available at every odd `n`):
- An odd ear decomposition builds the (odd-`n`) tournament from a base odd cycle by adding **odd ears**
  (odd-length paths). Each odd ear, closing against the existing structure, creates **directed odd
  cycles** — the vertices of `Ω`. So the odd ear decomposition is a **canonical incremental construction
  of `Ω`**, hence of `H = I(Ω,2)`.
- **BOLD (HYP-3100 candidate): `H` has an odd-ear recursion.** Adding an odd ear `ε` multiplies the
  partition function by a local factor depending on how `ε`'s new odd cycles attach to `Ω` (independent
  ear ⟹ `×I(ε-cycles,2)`; overlapping ⟹ a contraction). This is a *structured* deletion-contraction
  whose tree is the ear order — potentially the clean recursion the H-spectrum gaps `{7,21}` live in
  (`7,21 ∉` the odd-ear monoid, cf. THM-115's `{1+2^k}` monoid). **Test:** compute odd ear decompositions
  of small tournaments, track the H-factor per ear, look for the monoid.
- **LRC tie-in:** the fine-scale winding tournament is SC (HYP-3093) ⟺ has an ear decomposition; its ears
  = the **relation lattice `Λ(E)`** generators (the cycle space that controls `meas = Σ_{Λ} K(n)`). So the
  LRC lonely measure is an **ear-decomposition functional** of the winding tournament.

## 2. BRAVAIS-14 = the apex 2·7 (the crystallography of the tight locus)
`14 = 2·7` is the apex of BOTH the LRC and the 3-D Bravais classification, and both split as
**(7 = apex prime / Fano / crystal systems) × (2 = 2-adic / centering)** — the project's exact Cut⊕Cycle /
`14→7→2` descent. The dimension sequence `1,5,14,64` is the count of Bravais lattices; LRC is trivial at
`n≤5` and first-open at `n=14` — the **same two landmarks**.
- **WILD (HYP-3101 candidate): the LRC tight locus ↔ a Bravais/lattice classification.** The tight configs
  are EXACT TILERS (`d·{1..13}` tiles, THM-560) — i.e. **lattices** — and the census `{AP, GW}` is a
  *finite* classification of tilings, like the 14 Bravais types are a finite classification of 3-D
  lattices. The AP = the **most symmetric (cubic-analog)** tiler; GW = a **centered** variant
  (the Jacobsthal `×2` doubling = a 2-adic centering, cf. HYP-2918's "`x→2x mod 14` doubling gate"). The
  "7 crystal systems" ↔ the 7 sectors / Fano lines. **Test:** match the `{AP, GW}`-orbit structure (the
  Jacobsthal-gated doubling operad, HYP-2919) to the centering operations `{P,I,F,C}`; does the census
  count `a(2q)` follow the Bravais-by-dimension growth?

## 3. SAVITCH = the two-scale witness (coarse-guess / fine-refine)
Savitch's theorem proves `NSPACE(f)⊆DSPACE(f²)` by the **midpoint recursion**
`reach(u,v,2ℓ)=∃w[reach(u,w,ℓ)∧reach(w,v,ℓ)]` — guess a midpoint, recurse on halves. The LRC paper's
witness `t = s/14 + r/p` (coarse mod-14 **midpoint guess** + fine mod-`p` **refinement**) is exactly this
shape, and the verification cost `p^{(k+1)/2}` has the tell-tale **`/2` (the Savitch square-root of the
search)**. 
- **WILD (HYP-3102 candidate): LRC witness search is a `Savitch`/space-bounded computation, and the
  `(k+1)/2` exponent is the midpoint-recursion depth.** The project's two reductions — **Mode A** (`n→n-1`)
  and **Mode B** (`n→n-2`, Cayley–Dickson) — are the two halves of the Savitch split; the apex-periodicity
  / gamma-trick is the **memoized** midpoint. If LRC verification is `DSPACE(log²)`, the witness has a
  **bounded-width certificate** (the metagraph walk), which is the structural form of Conjecture 7.1.

## 4. THE φ⁴ PROOF LEVER (the proof-relevant gold, from HYP-3099)
The LRC coverage extremality is a **Lee-Yang / φ⁴** statement: the coverage PGF `Q(z)=Σ q_t z^t` has its 6
zeros on a near-circle `|z|≈R`, the **AP minimizes the off-circle variance `λ = Var(|roots|/R)`** (the
quartic `φ⁴` coupling — verified global argmin at `k=8,9`), and `λ` controls the coverage gap
(`corr(λ, cap−q_0) = +0.70→+0.85`, strengthening with `k`).
> **The lever:** the cap (`C(k+1,2)/91`, mac-mini's Pascal mass) is the **`λ=0` / on-circle / Gaussian
> value** (binomial coefficients = zeros on a circle, de Moivre–Laplace); the binding **dip** at `k=8,9`
> is the **`λ>0` off-circle correction**. So the cover bound `q_0 ≤ cap` would follow from a single
> **Lee-Yang/Asano monotonicity: deforming the coverage zeros off the circle (`λ>0`) does not raise `q_0`
> above its on-circle value.** This replaces the stalled moment-LP with a **zero-locus** argument (Asano
> contraction on the product form `Q(z)=∏(z−r_i)`), which is the natural tool when the object is a
> partition function. Unifies my order-3/4 cover-bound ladder + mac-mini's Pascal-cap + Hankel-dip.

## New signals (the owner's "create new signals to measure")
1. `min|root|` of `I(Ω,x)` — tournament Lee-Yang edge (tracks H-max, `→0` condensation).
2. `λ = Var(|roots|/R)` of the coverage PGF — φ⁴ coupling, minimized by AP (the extremality functional).
3. circle radius `R(k)` and `q_0 = q_6 R^6`.
4. **odd-ear H-factor** per ear (the monoid where `{7,21}` are forbidden).
5. **Savitch depth** = midpoint-recursion depth of the witness search.
6. **Bravais/centering type** of a tight tiler.

Honest: §4 is a concrete, possibly-provable lever (Asano monotonicity) — the priority. §1 (odd-ear
H-recursion) is a real, testable structural lead. §2 (Bravais) and §3 (Savitch) are wild but the
number-coincidences (`14=2·7`, `1,5,14,64`, the `(k+1)/2` exponent) are the kind the project has
repeatedly found load-bearing — flagged for a skeptical test, not asserted.

→ HYP-3099 (the two maps / Lee-Yang), THM-554/485 (partition functions), THM-560 (tilers=lattices),
THM-115 (`{1+2^k}` monoid, `{7,21}` gaps), HYP-2918/2919 (Jacobsthal doubling = centering?), HYP-3093
(fine-scale winding tournament), the paper's `t=s/14+r/p` two-scale witness, [[lrc14-thread]].
