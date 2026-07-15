# The parity tamer and the staircase locker — lockers, Rédei, Toda, and the squares-vs-triangulars structure of the axis

*opus-2026-07-15-S316; owner bundle: P vs NP via SAT; the locker/switch process;
squares vs triangulars; 3-triangular vs 4-square decompositions; Smith-chart
nomograms; binary conditions as tournament edges.*

## The common shape

The locker process — number k flips every multiple of k — computes the
zeta-transform mod 2 of the all-ones vector over the divisibility lattice:
the survivors are the squares because τ(n) is odd exactly there. The shape is
**a hard count, tamed by its parity**:

| hard count | parity statement | tamer |
|---|---|---|
| τ(n) (easy here, but the template) | τ odd ⟺ n square | ζ mod 2 |
| H(T) = #Hamiltonian paths (#P-hard in general) | H odd, ALWAYS (Rédei) | the OCF |
| #SAT (#P-complete) | ⊕SAT (⊕P-complete) | Toda: PH ⊆ P^#P via ⊕P |

Cook–Levin puts every NP problem inside SAT; a tournament on n vertices IS a
truth assignment to C(n,2) binary edge-variables (the owner's "binary
conditions as tournament edges" — and Kemeny/Slater orders on tournaments are
already NP-hard, so tournament-SAT is fully expressive). Rédei is the
miracle row: a #P-shaped count whose FIRST 2-adic digit is constant. The
natural question — how far up the 2-adic tower does the taming go? — is
answered by probe (ii) below: one more digit is c₃-governed, then it decays.

## Probe (i) — the axis level law (squares side)

x = Σ_v d_v² is a sum of n squares constrained by Σd_v = 0 and Landau. Exact
census n = 4..8 of the POPULATED levels:

> **Odd n: x ≡ 0 (mod 8). Even n: x ≡ n (mod 8). And the populated set is
> the FULL step-8 progression from the (near-)regular floor to the transitive
> ceiling (n³−n)/3 — no holes, all n ≤ 8.**

*Proofs (two lines each):* odd n ⟹ d_v = 2e_v with Σe_v = 0 ⟹ #odd e_v even
⟹ Σe² even ⟹ x = 4Σe² ≡ 0 (mod 8). Even n ⟹ all d_v odd ⟹ each d² ≡ 1
(mod 8) ⟹ x ≡ n. ∎ The no-holes completeness is verified exactly (n ≤ 8)
and should follow from an exchange walk realizing ±8 steps (F3's margin
structure); noted as provable-sketch. My earlier "step-16" memory was
THM-790's LINE-layer selection rule (|Δx| ≡ 8 mod 16 at odd n on d = m
lines) — a constraint on the FLIP layer, not on the level lattice: the two
mod-16 classes of levels are both fully populated, but the full-flip lines
connect them selectively. Squares-vs-triangulars is native here: the tile
count is m = T_{n−2}, the transitive ceiling is (n−1)n(n+1)/3, and
n² = T_{n−1} + T_n is literally the two-staircase gluing that underlies
Mode A/Mode B (the repo's triangle foundation). Lagrange's 4-square vs
Gauss's 3-triangular is the k = 4/k = 3 polygonal duality; its metagraph
avatar is the even-n/odd-n parity split of the d-vector (odd squares vs
doubled-triangular-ish evens).

## Probe (ii) — the 2-adic tower above Rédei decays at n = 7

With c₃ = C(n,3) − Σ_v C(s_v, 2) (a score function — the transitive-triple
complement):

- n = 4: **H mod 8 is a function of c₃ mod 4** (c₃ ≡ 0 → H ≡ 1; ≡ 1 → 3;
  ≡ 2 → 5). Perfect.
- n = 5, 6: for ODD c₃ the law "H ≡ 3 (mod 4) iff c₃ ≡ 1 (mod 4), H ≡ 1 iff
  c₃ ≡ 3" holds EXACTLY; for even c₃ it already splits.
- n = 7: all four (c₃ mod 4)-classes realize all four H mod 8 values — the
  pure-c₃ law is DEAD.

Reading: Rédei fixes digit 0 of H unconditionally; digit 1 is governed by
the 3-cycle count only up to n = 6; beyond, the higher odd-cycle-collection
terms of the OCF enter irreducibly. **The 2-adic taming of the #P-count is
a tower whose k-th level needs the k-th OCF layer — a Toda-like stratification
inside the repo's own formula.** (Honest negative: no clean H mod 4 formula
in score-invariants alone exists at n ≥ 7; scanned c₃, x, z — BH note: 3
invariant families scanned, 1 partial law found, decay point exact.)

## Probe (iii) — the staircase locker theorem (exact, proved)

Order the tiles by interval containment ([y',x'] ⊆ [y,x]); let Z be the
subtile zeta-transform mod 2 (each tile reads off the XOR of its subtiles —
the staircase's locker process). A gap-g tile has exactly **T_{g−1}**
subtiles, so:

> **Z(all-ones) = the indicator of tiles with gap ≡ 2, 3 (mod 4)** — the
> triangular-parity stripes. Verified n = 4, 5, 6. The integer locker's
> squares (τ odd) become the staircase locker's mod-4 gap stripes
> (T_{g−1} odd) — squares' role is played by the near-square gap classes,
> exactly the squares-vs-triangulars exchange the owner pointed at.

Z is unitriangular over F₂ (a bijection on tilings) but NOT class-functional
(image-class histograms spread already at n = 4): the locker lens genuinely
rescrambles the metagraph — logged as a new-edge-layer candidate (the
Z-conjugated wiggly walk) rather than a quotient tool.

## Nomograms and the engineering leads

- A nomogram encodes a ternary relation as collinearity; the Smith CHART is
  the Möbius-fold nomogram (z ↦ (z−1)/(z+1)). Our Smith DIAGRAM thread
  (metagraph resistance R_n, potentials φ) invites a literal chart: plot
  merged classes in the (φ, x) plane with the reversal fold as the Möbius
  boundary — the tournament-tiling explorer's next panel. [engineering lead]
- Tournament-SAT: the n = 9/10 census rungs ("does a class with invariant
  vector v exist") are SAT instances over C(n,2) edge variables — modern
  CDCL solvers + the symmetry-breaking of the base path could replace
  bespoke enumeration. [engineering lead, PyPI-adjacent]

## Cross-references

THM-790 (the line-layer mod-16 selection vs the level lattice), THM-855
F1–F6 (H as equilibrium measure — the #P-quantity as physics), the OCF canon
(the taming tower's rungs), eureka-zeckendorf-simplex-cuboid (3-triangular),
the-farey-14-row (S315: the other finite-depth law), MISTAKE-null: none —
probes verified before writing.
