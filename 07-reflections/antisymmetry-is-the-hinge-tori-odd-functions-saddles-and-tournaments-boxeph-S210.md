# Antisymmetry is the hinge: tori, odd functions, saddle points, and tournaments are one oddness

> **SCOPE CORRECTION — MISTAKE-224 (2026-07-21).** The independent skew-rank,
> standard-RPS first-integral, torus Poincare, inversion-fixed-point,
> Bernoulli-sawtooth, and Vandermonde-sign facts below survive. Their claimed
> equivalence does not. A Condorcet winner need not make the subtournament
> below it transitive; a torus admits a nowhere-zero field and needs no
> saddles; RPS level orbits are circles rather than tori; and general
> antisymmetry does not force zero column sums. For LRC, inversion only pairs
> the signed THM-2047 walls/components. Treat it as a symmetry reduction, not
> a positivity or AP-extraction mechanism.

*boxeph-2026-07-21-S210. Owner: see how tori relate to odd-valued functions, saddle points, and
tournaments. Builds on boxeph-S207 (bagel/torus, deficit-1 = reduced Euler), S209 (LRC = toric
arrangement), THM-1820 (relation lattice, t→−t mirror pairs), THM-473 (skew tournament matrix = Hermite),
THM-1830 (the 3-cycle atom). All four pillars verified in
`04-computation/tori_odd_saddles_tournaments_boxeph_S210.py`.*

## The one hinge

A tournament's payoff matrix `M = A − Aᵀ` (`M_{ij}=+1` if `i` beats `j`, `−1` if beaten, `0` on the
diagonal) is **antisymmetric** — it is the *odd* part of the adjacency. Every relationship the question
asks about follows from that single fact, because **antisymmetry forces even rank, and even rank forces
odd kernels, saddles, and toroidal recurrence.** Tori, odd functions, saddle points, and tournaments are
four faces of the same antisymmetry.

## Pillar 1 — odd ⟹ even rank ⟹ odd support (the game-theoretic saddle point)

`M` antisymmetric ⟹ `rank(M)` is **even** (verified for *all* tournaments `n≤5`) ⟹ every odd-size
principal block is singular (its determinant is a Pfaffian of odd order `= 0`). Read `M` as a **symmetric
zero-sum game** (the "tournament game"): by antisymmetry its value is `0`, and an interior equilibrium
needs `Mp = 0` on its support — which requires the support block to be **singular**, i.e. of **odd
cardinality**. This is the Fisher–Ryan theorem, and the census confirms it: over all `n≤5` tournaments,
**every optimal support has odd size** (`1, 3, 5`).

The two poles are the repo's own dichotomy:
- **Transitive ⟹ support 1**: a Condorcet winner (source), a **pure saddle point** of the game
  (minimax = maximin at a single strategy).
- **Cyclic ⟹ support ≥ 3 (odd)**: no pure saddle; the smallest case is the **3-cycle** = rock–paper–
  scissors, uniform on 3 (the 3-cycle atom, THM-1830). Regular odd-`n` tournaments: uniform optimal,
  support `n`.

So **a pure saddle point of the tournament game exists iff the tournament is transitive** — the exact
transitive/intransitive boundary, now read as the presence/absence of a game equilibrium.

## Pillar 2 — the torus *needs* saddles (χ = 0)

A torus has vanishing Euler characteristic, and Morse/Poincaré–Hopf forces **saddles to balance the
extrema**. Verified: `Poincaré(Tⁿ) = (1+t)ⁿ`, Betti `= C(n,k)`, `χ = Σ(−1)^k C(n,k) = 0`. The standing
**bagel** `T²` under the height function has exactly **1 max, 2 saddles, 1 min** (`χ = 1−2+1 = 0`), and
the **two saddles are `b₁ = 2` — the signature of the handle** (`H₁(T²)=ℤ²`). This is the same `χ=0` that
appears in S207 as the **deficit-1** (`bagel − cake = Tₙ − 1`, the reduced-Euler / handle term): the
torus's hole *is* its saddle pair.

## Pillar 3 — transitive = gradient (no torus); cyclic = recurrent (an invariant torus)

Run the **replicator dynamics** `ẋᵢ = xᵢ(Mx)ᵢ` (since `xᵀMx = 0` for antisymmetric `M`):
- **Transitive** (`0` beats `1,2`; `1` beats `2`): the flow runs **downhill to the Condorcet winner** —
  a **sink / gradient flow**, converging to `(1,0,0)` (verified). No recurrence, no torus.
- **3-cycle** (rock–paper–scissors): the product `H = x₀x₁x₂` is **exactly conserved** —
  `dH/dt = H·Σᵢ(Mx)ᵢ` and `Σᵢ(Mx)ᵢ = Σⱼ xⱼ·(column sum of M) = 0` because `M` is antisymmetric
  (column sums `= [0,0,0]`, verified). So the orbits are **closed level curves of `H` around the center
  `(⅓,⅓,⅓)`** — a **recurrent flow foliating an invariant torus** (verified: RK4 `H`-drift `≈10⁻¹⁶`, the
  orbit returns exactly to its start).

So the **intransitive part of a tournament IS the toroidal / recurrent set**, while the transitive part
is gradient-like — the dynamical-systems (Conley) reading of the reify ladder: transitive = the cold
gradient vertex, cyclic = the recurrent torus. The saddle point (equilibrium) is *pure* for the gradient
case and *mixed-on-odd-support* for the toroidal case.

## Pillar 4 — odd functions live on the torus through the involution θ ↦ −θ

The reversal involution `σ: θ ↦ −θ` on `Tⁿ=(ℝ/ℤ)ⁿ` has exactly **`2ⁿ` fixed points** — the **2-torsion**
`{0,½}ⁿ` (verified) — and **odd functions vanish at all of them**. This is the geometric home of the
oddness:
- The LRC good-set **measure is EVEN** under `t↦−t` (the far-set weight `ĝ(k)=−sin(2πkδ)/(πk)` is a
  **sinc = odd/odd = even** in `k`, verified) — so `|G_δ|` is `t↦−t`-invariant, and THM-1820's "mirror
  pairs" are `σ` acting on the *fold structure*.
- The genuinely **odd** object is the **signed-discrepancy sawtooth** `B₁({x})={x}−½`, Fourier coefficient
  `c_k = 1/(2πik)` with `c_{−k} = −c_k` (odd in `k`) — the antisymmetric sector.
- The **transitivity Vandermonde** `∏(aᵢ−aⱼ)` is the **sign character** (odd under `Sₙ`: one transposition
  flips its sign, verified `V=540 ↦ −540`) — literally the same oddness as the tournament payoff.

So the odd sector of the torus (sawtooth, sign character, `σ`-antisymmetric part) is exactly the
antisymmetric "beats" structure of a tournament.

## The synthesis

```
        ANTISYMMETRY  (M = A − Aᵀ,  the odd "beats" relation)
        │
        ├─ even rank ⇒ ODD support (Fisher–Ryan) ⇒ pure game SADDLE POINT ⇔ transitive
        ├─ TORUS χ=0 ⇒ SADDLES balance extrema (bagel: 1 max/2 saddle/1 min = the handle = deficit-1)
        ├─ replicator: transitive = GRADIENT→sink (no torus);  3-cycle = conserved H ⇒ invariant TORUS
        └─ involution θ↦−θ (2ⁿ 2-torsion fixed pts): ODD sector = sawtooth B₁ = sign character = Vandermonde
```

Tori, odd-valued functions, and saddle points are not three separate neighbours of tournament theory —
they are the geometric (`χ=0` needs saddles), analytic (odd functions on the `θ↦−θ` torus), and
game-theoretic (odd support, pure-saddle ⇔ transitive) shadows of the single fact that a tournament is an
**antisymmetric** object. The transitive/cyclic dichotomy is gradient-vs-recurrent, saddle-present-vs-
absent, and even-vs-odd all at once.

## Leverage

- The **3-cycle atom = the elementary invariant torus** (recurrent set) of tournament replicator dynamics;
  the repo's "unstable non-transitive = one 3-cycle atom" (THM-1830) is a Conley-index statement — the
  atom is the minimal non-gradient recurrent block. A Conley/Morse decomposition of tournament space is
  available.
- **Odd support is forced by even rank** — this gives the tournament-game solution structure (Fisher–Ryan)
  a Pfaffian/skew-spectral proof, connecting to THM-473 (skew matrix = Hermite, purely imaginary `±iλ`
  spectrum = the toroidal frequencies).
- The **odd sector on the `θ↦−θ` torus** (sawtooth, sign character) is where the LRC reality-symmetry
  (THM-1820 mirror pairs) and the tournament antisymmetry meet — a single involution governs both.

Links: HYP-8835, THM-1820, THM-473, THM-1830,
[[cake-bagel-and-fibonacci-are-one-pascal-triangle-boxeph-S207]],
[[orlik-solomon-is-a-repo-wide-pattern-toric-arrangements-are-the-lrc-lens-boxeph-S209]].
