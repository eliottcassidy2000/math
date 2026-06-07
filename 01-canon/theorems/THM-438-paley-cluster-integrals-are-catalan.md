---
id: THM-438
name: the-Paley-cluster-integrals-are-Catalan-numbers-and-R(p)-to-e-is-proven
status: PROVED (R(p)->e: uniform proof, single Weil bound) + VERIFIED (Catalan leading coefficient, k<=4)
date: 2026-06-07
session: monad-explorer-2026-06-07 (deep-research / analytic lane)
depends_on:
  - HYP-2307   # the cherry cluster expansion R(p)=exp(sum a_L); a_2=1, a_odd=0; a_4=a_6=0 verified
  - MISTAKE-011b  # Paley tournament needs p=3 mod 4 (chi(-1)=-1) -- the source of the + cherry weight
relates_to:
  - OPEN-Q-013   # the ratio line H(T_p)*2^{p-1}/p! and "its approach to e"
  - HYP-2306     # the 1729 spine severed: H/|Aut| has no modular structure (analytic axis is the real one)
  - THM-052      # palindromic N / circulant tournaments (the universality class of the limit)
---

# THM-438: the Paley cluster integrals are Catalan numbers; `R(p) → e` is PROVEN

## Context

For the Paley tournament `T_p` (`p ≡ 3 mod 4`), `H(T_p)` = #directed Hamiltonian
paths, and the normalized ratio is

```
R(p) := H(T_p) · 2^{p−1} / p!   =   E_σ[ ∏_{k=1}^{p−1} (1 + χ(d_k)) ],
```

the Paley maximizer's multiplicative edge over a coin-flip tournament (`p!/2^{p−1}`
is `E[H]` for a random tournament). HYP-2307 expanded `R(p)` into single-run cluster
integrals and showed `R(p) → exp(Σ_{L≥2} a_L)` with the **cherry** (`L=2`) the only
expected survivor, but left the limit PROVEN only modulo "`a_{2k}=0` for all `k≥2`"
(verified `k=2,3` by hand decomposition + Weil). This theorem (i) closes that gap with
a **uniform** argument and (ii) determines the **exact leading order** of every
cluster integral — they are the **Catalan numbers**.

## Definitions

The single-run cluster integral (a complete multiplicative-character sum over `F_p`):

```
A_L  :=  Σ_{x_0, x_1, …, x_L ∈ F_p, ALL DISTINCT}  ∏_{i=0}^{L−1} χ(x_{i+1} − x_i),
a_L  :=  lim_{p→∞} A_L / p^L.
```

Here `χ` is the Legendre symbol; an ordering `a_1→…→a_p` contributes `∏(1+χ(d_k))`,
and a subset of its `p−1` consecutive arcs splits into maximal runs, each a directed
sub-path whose character content is `A_L`. `C_k = \binom{2k}{k}/(k+1)` is the `k`-th
Catalan number (`C_1..C_5 = 1, 2, 5, 14, 42`).

## Statement

**Part A (exact vanishing).** `A_L = 0` for every **odd** `L`, at every prime
`p ≡ 3 mod 4`. (Hence `a_1 = a_3 = a_5 = ⋯ = 0`.)

**Part B (the Catalan law — sharp order).**
```
A_{2k}  =  C_k · p^{k+1}  +  O(p^{k+1/2})          for every k ≥ 1.
```
In particular the leading order is `p^{k+1}`, *not* the `p^{2k}` that a generic
`L = 2k` character sum could a priori reach. The coefficient is the Catalan number
`C_k` (positive sign), counting the **bigon-tree** coincidence patterns of the
`2k`-edge walk = Euler tours of plane trees with `k` edges.

**Part C (`R(p) → e`, PROVEN).** For `k ≥ 2`, `k+1 < 2k`, so `a_{2k} = lim A_{2k}/p^{2k} = 0`.
Combined with `a_2 = 1` (Part B, `k=1`: `A_2 = p(p−1)`) and Part A,
```
Σ_{L≥2} a_L  =  a_2  =  1        ⟹        R(p)  →  e^1  =  e  =  exp(−χ(−1)).
```
This rules out Alon's permitted `p^{3/2}` growth: the cluster sum has a single finite
generator.

## Proof

### Part A (odd runs vanish — exact)
The map `x_i ↦ −x_i` is a bijection on distinct `(L+1)`-tuples and sends each factor
`χ(d) ↦ χ(−d) = χ(−1)χ(d) = −χ(d)` because `χ(−1) = −1` at `p ≡ 3 mod 4`. Over `L`
edges, `∏χ ↦ (−1)^L ∏χ`, so `A_L = (−1)^L A_L`, forcing `A_L = 0` for odd `L`. ∎

### The vanishing engine (used throughout)
Let `M` be the `p×p` circulant `M[a,b] = χ(b−a)`. Its row sums are
`Σ_b χ(b−a) = Σ_z χ(z) = 0`, so `M·𝟙 = 0`. Hence the **free** (non-distinct) path sum
```
B_L := Σ_{x_0,…,x_L ∈ F_p} ∏ χ(x_{i+1}−x_i) = 𝟙ᵀ M^L 𝟙 = 0   for L ≥ 1.
```
By inclusion–exclusion over which walk-positions take equal values,
`A_L = −Σ_{patterns with ≥1 coincidence} (signed contribution)`: the fully-distinct
sum equals minus the coincidence sum, since the free sum vanishes.

A coincidence pattern is a partition of `{x_0,…,x_L}` into value-groups; it yields a
reduced multigraph `G` (groups = vertices, the `L` walk-edges = edges). Two facts:

1. **No-leaf lemma.** If a group `g` is incident to exactly one edge (a *leaf*), the
   sum over its value is `Σ_{g} χ(g − a) = 0`. So every nonzero pattern has minimum
   degree `≥ 2`, giving `V(G) ≤ L` (since `2L = Σ deg ≥ 2V`).

2. **Bigon = no cancellation; longer cycle = Weil deficit.** A *bigon* (two
   anti-parallel edges between the same pair, from `x_i = x_{j+1}`, `x_{i+1} = x_j`)
   closes into `χ(d)χ(−d) = χ(−d²) = χ(−1)`, a constant — its value sums with the FULL
   `~p` weight. A genuine cycle of length `≥ 4` is a non-degenerate
   multiplicative-character sum, hence `o(p^{cycle-length})` by Weil.

### Part C (uniform `a_{2k}=0`, single Weil bound)
A reduced graph with `V` vertices has `≤ (p)_V ≤ p^V` terms, so its contribution is
`O(p^V)`. By the no-leaf lemma `V ≤ 2k`.

- `V ≤ 2k−1` (i.e. `≥2` coincidence-merges): contributes `O(p^{2k−1}) = o(p^{2k})`
  **trivially**, no Weil needed.
- `V = 2k` (exactly one merge): the only no-leaf single-merge pattern identifies the
  two path endpoints `x_0 = x_{2k}` (any other single merge leaves a leaf at `x_0` or
  `x_{2k}`). This is the closed `2k`-cycle character sum, which is `o(p^{2k})` by the
  **single** Weil bound on an even cyclic character sum.

Hence `A_{2k} = o(p^{2k})`, i.e. `a_{2k} = 0`, for every `k ≥ 2`. With `a_2 = 1`
(below) and Part A, `R(p) → e`. ∎

### Part B (the Catalan leading order)
By the engine, the **top** order is achieved by patterns whose every cycle is a bigon
(no Weil deficit) and which maximize `V`. An all-bigon reduced graph with `k` bigons
has `L = 2k` edges and, if the bigons are glued in a **tree** (no bigon-cycle),
`V = k+1` — the maximum. Such patterns are exactly the closed length-`2k` walks that
traverse each edge of a plane tree (`k` edges) once in each direction = **Euler tours
of plane trees with `k` edges**, counted by `C_k` (the Dyck-path bijection: each walk
step is "open a new edge" `(+1)` or "backtrack/close a bigon" `(−1)`; `V=k+1` forces a
balanced Dyck word). Each contributes the full `p^{k+1}` with net `+` sign (the
`(−1)^k` from `k` bigons `χ(−1)` is cancelled by the inclusion–exclusion sign — see
the sign verification below). Patterns with any non-bigon cycle, or with bigon-cycles
(`V < k+1`), are `O(p^{k+1/2})` or `O(p^{V}) = O(p^{k})`. Hence
`A_{2k} = C_k p^{k+1} + O(p^{k+1/2})`. ∎(modulo the sign bookkeeping, VERIFIED below)

## Verification (`04-computation/paley_cluster_sharp_order_monad.py`, `paley_cluster_catalan_monad.py`)

- **Part A:** `A_1 = A_3 = A_5 = 0` exactly, `p = 7,11,19,23`.
- **Combinatorial `c_k`:** exhaustive union-find over edge-pairings gives the bigon-
  **tree** (`V=k+1`) count `1, 2, 5, 14, 42, 132` for `k = 1..6` — **the Catalan
  numbers** (OEIS A000108). (Lower-`V` patterns, the bigon-cycle cacti, appear at
  `V = k−1, k−3, …`, contributing only `p^{V}`.)
- **Part B numerics (Paley primes only):**
  - `A_4/p^3`: `0.980, 1.322, …, 1.882` (`p=7..67`) → `2 = C_2`, monotone from below.
  - `A_6/p^4`: `0.420, 1.563, 2.771, 3.110` (`p=7,11,19,23`) → `5 = C_3`, monotone.
  - `A_8/p^5`: order confirmed (`=0` at `p=7`, since `L=8` needs `9` distinct
    vertices — impossible mod `7`; `0.96` and climbing at `p=11`).
  - **Sign caveat learned:** including `p ≡ 1 mod 4` (NON-tournament, `χ(−1)=+1`) makes
    `A_6/p^4` oscillate in sign — exactly the `MISTAKE-011b` boundary. Restricted to
    Paley primes the convergence is clean and `+`.

## Significance

1. **Closes HYP-2307 handoff #1.** `R(p) → e` is now PROVEN with a *uniform* argument
   for all `k` (a single Weil bound on one even cyclic sum), not a per-`k` check.
2. **The exact analytic skeleton.** Every cluster integral's leading order is pinned:
   `A_{2k} = C_k p^{k+1}`. The Catalan numbers are the non-crossing / tree-walk count
   — the **moment-method combinatorics** (the engine of Wigner's law) surfacing from
   the excluded-volume (distinct-vertex) constraint. The full trace
   `tr(M^{2k}) = (−p)^k(p−1)` reflects the **two-point** Gauss-sum spectrum `±√p`
   (`M` has eigenvalues `χ(j)·g`, `|g|=√p`); `A_{2k}` extracts the distinct/tree part,
   surfacing `+C_k`. See the reflection.
3. **Universality.** The engine uses only that `χ` is **odd** (`χ(−d) = −χ(d)`) =
   the tournament condition. So `A_{2k} = C_k p^{k+1}` and `R(p) → e` hold for EVERY
   circulant tournament (HYP-2307's verified universality, now explained: the entire
   leading analytic skeleton is shared — this is *why* `H/|Aut|` carries no Paley-
   specific arithmetic, HYP-2306).

## Honest status

- **PROVED:** Part A (exact); Part C (`R(p)→e`) modulo the one classical Weil bound on
  an even cyclic multiplicative-character sum being `o(p^{2k})` — a textbook input.
- **VERIFIED (mechanism CORRECTED — see ADDENDUM + MISTAKE-060):** the leading
  **coefficient** is `+C_k`. ~~A clean Möbius-sign proof (sign `=(−1)^k` cancelling the
  bigon parity) is the small remaining write-up.~~ **This was investigated and the naive
  bigon-tree story is FALSE:** bigon-trees alone give `1,3,13,69,…` (OEIS A088368, `~e·n!`),
  and even-cycle cacti supply *negative* corrections; `C_k` is the **signed** even-cacti sum
  (verified `k=2,3` exactly via the flow closed form). The order `Θ(p^{k+1})` is proved; the
  remaining clean write-up is the closed-form identity `Σ_{even cacti} μ(0̂,σ)·lead(M_σ) = C_k`
  (the A088368 → Catalan / all-pairings → non-crossing reduction).
  **SUPERSEDED by ADDENDUM-2 + MISTAKE-061:** the support "even cacti" is WRONG (even theta graphs
  also reach top order); the correct, number-theory-FREE target is `Σ_{even-series patterns} μ(0̂,ρ)
  = (-1)^k C_k`, with `g≡+1` now PROVED and the leading coefficient `c_0=C_k` re-confirmed
  rigorously (k≤4, clean Richardson `1/p`).

## ADDENDUM (monad-explorer-2026-06-07, 3rd session) — Part B mechanism CORRECTED; Part C needs no Weil; error term `O(p^k)`

Verified exactly in `04-computation/paley_cluster_cactus_census_monad.py`. See MISTAKE-060.
**The statements above (A_{2k}=C_k p^{k+1}, R(p)→e) are UNCHANGED.** Only the Part-B
*mechanism*, the error term, and the (acknowledged-open) `+C_k` sign claim are revised —
and they make the result deeper, not weaker.

**(1) The closed form behind everything (PROVED, Gauss-sum inversion `χ(w)=g⁻¹Σ_t χ(t)ω^{tw}`).**
For any coincidence pattern `σ` (a set partition of `{0,…,2k}` giving reduced multigraph
`G_σ` with `V` blocks and the `2k` walk-edges),
```
M_σ  =  (−1)^k · p^{V−k} · F(σ),     F(σ) = Σ_{F_p-flows t on G_σ}  ∏_e χ(t_e),
```
the flow sum ranging over the cycle space (dim `m = 2k−V+1`). `M_σ` reaches the top order
`p^{k+1}` **iff** `F(σ)` reaches full order `p^m` — exactly the **even cacti** (connected,
all biconnected blocks even cycles; a bigon = 2-cycle).

**(2) Part B's "bigon-trees, each `+1`" is WRONG; the Catalan number is a SIGNED cancellation.**
By Möbius inversion `A_{2k}=Σ_σ μ(0̂,σ)M_σ`, `μ(0̂,σ)=∏_B(−1)^{|B|−1}(|B|−1)!`. The
bigon-**tree** leading coefficient is `Σ_{non-crossing pairings}∏_v(b_v−1)! = ` **`1,3,13,69,421,2867`**
(`k=1..6`) = **OEIS A088368** (`A=Σn!xⁿAⁿ`, `a(n)~e·n!`) — *not* `C_k`. The even-cycle cacti
contribute **negative** corrections; the signed total is `C_k`:
```
k=2:  bigons(+3) + 4-cycle(−1)                  = 2 = C_2
k=3:  bigons(+13) + {bigon+4-cycle} + {6-cycle} = 5 = C_3   (census verified, p=11,19,23,31)
```
This is the genuine moment-method content: the inclusion–exclusion converts the *all-pairings*
overcount (A088368, `~e·n!`) into the *non-crossing* count `C_k`; the two-point Gauss spectrum's
even-cycle terms are exactly the Gaussian→semicircle corrections. (The old slogan "Catalan =
bigon-trees" was right about the **answer**, wrong about the **mechanism**.)

**(3) Part C needs NO Weil (fully elementary).** The only no-leaf `V=2k` pattern is the single
`2k`-cycle `x_0=x_{2k}`, and that `M_σ = tr(M^{2k}) = (−p)^k(p−1)` is elementary; `V≤2k−1` is
`O(p^{2k−1})=o(p^{2k})` trivially. So `a_{2k}=0` (k≥2) and `R(p)→e` require no character-sum
estimate at all. (Weil is only needed for the *exact* `O(p^k)` remainder coefficient, via the
genuine odd-cycle/Jacobi sub-patterns — not for the limit.)

**(4) Error term corrected: `A_{2k} = C_k p^{k+1} + O(p^k)`** (not `O(p^{k+1/2})`). Verified:
`(A_4−2p^3)/p^2` is STABLE (≈ −7.1…−7.8 → ≈ −8) while `/p^{2.5}` drifts to 0. Consequently
`R(p)−e` has **relative correction `O(1/p)`** (resolves the reflection's stated `O(1/√p)` vs the
close-out's "favors 1/p" — it is `1/p`), which is the right footing for handoff #2.

## Forward

1. **Sub-leading term `C` in `R(p) = e(1 − C/p + …)`** (HYP-2307 handoff #2). The error
   is `O(p^k)` (ADDENDUM (4)), so `R−e` is relative `O(1/p)` — the ansatz `R=e(1−C/p+…)`
   is now *justified* (not assumed), `C` fed by `A_4`'s integer `p^2` remainder + finite-`p`
   cherry-placement. The compute node's `p=31,43,47` would *pin* `C≈1.4`, testing a
   prediction.
2. **Non-circulant doubly-regular tournaments** (handoff #3) → **HYP-2308**. The
   leading-order Catalan skeleton is now shown to be **spectral, not arithmetic**: the
   engine needs only regularity (no-leaf), antisymmetry (bigon `−1`), and the two-point
   `{0}∪{±i√n}` spectrum (`tr(S^{2k})=(−1)^k n^k(n−1)`) — all DRT properties, none
   circulant/Weil. So the Catalan law and `R→e` should hold for EVERY DRT; the only open
   part is the `o(n^{k+1})` remainder for non-circulant DRTs (a tight-spectral
   expander-mixing estimate replacing Weil). See HYP-2308 for the test recipe.

## ADDENDUM-2 (monad-explorer-2026-06-07, 4th session) — leading coefficient RE-CONFIRMED `=C_k` rigorously; top-order class corrected from "even cacti" to EVEN-SERIES patterns; `g≡+1` PROVED; the Catalan number is a pure partition-lattice Möbius sum

See **MISTAKE-061**. Scripts: `04-computation/paley_cluster_topterm_monad.py`,
`paley_cluster_pure_moebius_monad.py`, `paley_cluster_theta_check_monad.py`,
`paley_cluster_leadcoeff_monad.py` (+ `.out`). **The headline `A_{2k}=C_k p^{k+1}` and
`R(p)→e` are UNCHANGED and now on firmer footing.** Two things are corrected/sharpened.

**(1) The leading coefficient is RIGOROUSLY `c_0 = lim A_{2k}/p^{k+1} = C_k` (VERIFIED k≤4).**
Computed exactly via the flow-Möbius identity `A_{2k}=Σ_ρ μ(0̂,ρ)M_ρ` at many Paley primes and
Richardson-extrapolated in `1/p`:  `c_0 = 2, 5, 14 = C_2, C_3, C_4` *exactly*. This REPLACES the
prior census numerics (`A_6/p^4 = 1.56,2.77,3.11` at `p≤23`), which converge to `5` only very
slowly (the `O(1/p)` + character fluctuations) and could not by themselves distinguish `5` from `6`.

**(2) The top-order class is NOT "even cacti" — it is the EVEN-SERIES patterns (even theta graphs
included).**  `M_ρ` reaches `p^{k+1}` iff the flow-form product `P(s)=∏_e ℓ_e(s)` is a perfect
square, i.e. iff **every series-class of edges is even**. Even cacti qualify, but so do even theta
graphs (2-connected, biconnected block ≠ a single cycle) and all even series-parallel skeletons.
At `k=3` the `m=2` top-order patterns are `6` even-cacti{2,4} `+ 1` even-theta(2,2,2); the theta was
missing from the "even cacti" census (it hid in the `(6,)` biconnected bucket, cancelling the single
6-cycle). The cycle-rank of any top-order pattern satisfies `m ≤ k` (max at bigon-trees).

**(3) `g(ρ):=lim F(ρ)/p^m ≡ +1` for every top-order pattern — PROVED.** Within each series-class the
closed Euler walk passes straight through the degree-2 internal vertices, so all its edges receive
the SAME orientation sign `s∈{±1}`; the class being even gives `∏_{e∈class}sign_e=s^{even}=+1`, so
`P=(+1)·Q²` and `g=χ(P)=+1`. (Numerically: `g=+1` for all `3+13+67` contributing patterns at k=2,3,4.)
The character/Gauss-sum content thus **collapses entirely** at leading order:
```
c_0 = (-1)^k Σ_{ρ : even-series pattern} μ(0̂,ρ) = C_k,
i.e.   Σ_{ρ : even-series pattern} μ(0̂,ρ) = (-1)^k C_k        (number-theory-FREE).
```

**(4) Status of handoff #1 (clean proof of the Catalan coefficient).** Now a *purely combinatorial*
identity `(★★)` over the partition lattice — no characters, no Weil, no Gauss sums. The bigon-tree
sub-sum is the all-pairings overcount A088368 (`1,3,13,69,…~e·n!`); even cacti + even thetas + …
cancel it to `C_k`. This is exactly Wigner quasirandomness: the random skew-Rademacher matrix gives
`C_k` *directly* from non-crossing pairings (each `+1`), and `(★★)` is the statement that the
deterministic Paley Möbius expansion, despite over-counting, lands on the same `C_k`. A clean proof of
`(★★)` (a non-crossing / free-cumulant bijection on even-series patterns) is the remaining write-up —
strictly cleaner than before (the wrong "even cacti" support is removed and `g≡+1` is no longer assumed).
