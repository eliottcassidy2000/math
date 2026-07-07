# LRC as rotational coloring, and the dihedral ledger

*mac-mini-2026-07-07-S48 (HYP-5047). Owner: think colorings, dihedral groups; LRC is the
assertion that a distance graph which may genuinely require n colors can always be
n-colored by a single rigid rotation of the color circle — with D = {1,…,n−1} as the
tight roots-of-unity case.*

## The dictionary, stated precisely

For a speed set D (|D| = n−1), the distance graph G(ℤ, D) joins x, y iff |x−y| ∈ D. A
witness time t with min_d ||d·t|| ≥ 1/n induces the **rotation coloring**
χ_t(x) = ⌊n·frac(x t)⌋: adjacent vertices differ by d ∈ D, and ||d t|| ≥ 1/n walks the
phase at least one full color cell, so χ_t is a proper n-coloring. LRC(n) is exactly:
*every* such distance graph admits a proper n-coloring of this maximally rigid form —
one number t, colors = arcs of the circle, the coloring equivariant under ℤ-translation
by construction. The classical bridge (Zhu: circular chromatic numbers of distance
graphs) says χ_c(G(ℤ,D)) = 1/κ(D) with κ the lonely-runner constant; LRC(n) ⟺
χ_c ≤ n on this class.

Three project facts snap into place in this frame:
- **The AP is the roots-of-unity case.** D = {1,…,n−1} contains K_n (vertices 0..n−1
  pairwise adjacent): χ = n genuinely. The unique witness t = 1/n places the phases at
  the n-th roots of unity — the coloring exists but with ZERO slack: every color cell
  boundary touched. Tightness of the AP = rigidity of the regular n-gon.
- **The k=8 criticality is a coloring threshold.** klein-S155's mean-out-degree-1
  criticality at (k,θ) = (8, 1/7) reads: 8 phases, 7 color cells — the pigeonhole edge
  where a proper "cell-injective" placement first fails on average; the μ-floor
  measures how often rotation still finds a free cell.
- **The dihedral group is the ledger's symmetry.** The color circle carries D_n, not
  just ℤ_n: rotations = the t-shifts (τ downstairs), the reflection = x ↦ −x =
  step-reversal (THM-639) = the tent evenness Λ(−ψ) = Λ(ψ) that klein-S165 identifies
  as the mechanism halving the H-sums. Every exact law we have proved lives on
  dihedral-symmetric data: THM-641's four Bernoulli terms pair under the reflection;
  THM-643/644/648's parities come from the order-2 elements (ρ, flip) of the same
  dihedral action upstairs. **One dihedral action, three floors of the project.**

## The gate work this frame guided (session results)

- **THM-648 (proved):** blue self-loops only at even n — the odd-n obstruction is the
  score-complement vector's two defective entries at the base path's endpoints; in
  coloring language, the color-circle reflection fixes a cell center only when the
  number of cells is even. The mod-2 program on the strict blue/black structure is now
  closed (parities THM-643/klein-S161, fibers THM-644, Anti-Rédei THM-647, rigidity
  klein-S163, self-loops THM-648).
- **The no-separated-cherry structure lemma (proved, 3 lines):** a k=8 shape with no
  triple e_c ≥ 50(e_a + e_b) satisfies e_max < 50(e_2 + e_3) — the whole shape is one
  bounded band, possibly plus one isolated small speed. Census on 4000 shapes at
  diameter ≥ 27: every no-cherry shape has μ_{1/7} ≥ 0.998 (bar 0.675 — 48% headroom;
  0/600 below). klein-S165's moderate-spread residual is thus never binding
  empirically, and the lemma names the mechanism: no separation ⟹ bounded bands ⟹
  8-point decorrelation. The remaining proof gap is now one shaped statement:
  *single-band 8-shapes at diameter ≥ 27 have μ ≥ 0.675* (observed ≥ 0.998).

## The reflection proper

The owner's coloring sentence is the project's best one-line self-description: LRC
claims that a chromatic problem whose lower bound can be tight (the AP forces n) is
always solved by the most symmetric possible coloring — a group rotation. All our hard
cases are the configurations that try to break equivariance (the moat's signed
cancellation, the composite-14 defect, the two hard cores' transitive tails), and all
our proved tools are statements that equivariance survives averaging (CV floors,
Bernoulli laws, parity theorems). The dihedral reflection is not decoration: it is the
one symmetry both projects share at every level, and it keeps paying — palindromes,
tent evenness, ρ-anti-automorphisms, and now the even-n self-loop law are the same
coin, spent four ways.
