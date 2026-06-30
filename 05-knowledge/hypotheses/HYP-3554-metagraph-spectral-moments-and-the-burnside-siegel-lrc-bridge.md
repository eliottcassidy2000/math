---
id: HYP-3554
title: Metagraph SPECTRAL MOMENTS -- a new cheap invariant extending klein's signed cycle index P_n (THM-587): the derivatives P_n^{(r)}(1) give the mean level kbar=P_n'(1)/P_n(1), the spectral variance Var(k) (1,2.19,4.06,5.25,6.20 for n=3..7), and the mean metagraph eigenvalue d-2kbar (1,1.5,1.33,1.07,0.72 -> 0, the signed-action asymmetry shadow), computable from n! past the 2^{C(n,2)} wall; AND the BURNSIDE<->SIEGEL<->LRC MOMENT BRIDGE -- P_n (Burnside average over S_n) is the finite analog of the Siegel transform (average over SL_n(Z)), arXiv:2507.05905 (Han-Lee) computes Siegel 1st/2nd MOMENTS with CONGRUENCE conditions, and the LRC floor (THM-579) needs the 2nd moment CV(N_R) of the sheet count with the COVERING (=congruence) restriction; so the paper's congruence-Siegel 2nd-moment is the natural tool for THM-579's gatekeeper, with the metagraph's CV(H) over iso classes its finite testbed
status: VERIFIED (P_n reconstruction matches A000568/SC n=3..7; spectral moments + CV(H) computed). The moment-bridge + the THM-579 application are SYNTHESIS / a proposed direction, not a proof.
source: mac-mini-2026-06-29-S19
related:
  - THM-587   # the signed cycle index P_n -- this extends it to its moments
  - THM-584   # complement = antipodal; eigenvalues d-2k
  - THM-579   # the LRC floor as a sheet-count CV (the 2nd moment the paper computes)
  - HYP-3544  # klein equivariant homology (the per-level structure these moments summarize)
  - HYP-3550  # the resonance Euler product (the floor); CV(N_R) controls the channel coupling
external: arXiv:2507.05905 (Han-Lee, Siegel-transform 1st/2nd moments with congruence conditions); Siegel/Rogers/Schmidt mean-value theory
results:
  - 04-computation/metagraph_spectral_moments_macmini_20260629.py
  - 05-knowledge/results/metagraph_spectral_moments_macmini_20260629.out
---

# HYP-3554 -- metagraph spectral moments and the Burnside<->Siegel<->LRC bridge

## (1) New invariant: the metagraph spectral moments
klein's THM-587 gives `P_n(x) = sum_k mult(k) x^k`, the level-multiplicity generating function (eigenvalue
`d-2k`, `d=C(n,2)`). Its EVALUATIONS were known (`P_n(1)=A000568`, `P_n(-1)=SC`); its DERIVATIVES are a new
cheap invariant -- the moments of the spectrum:
- mean level `kbar = P_n'(1)/P_n(1)`; spectral variance `Var(k) = P_n''(1)/P_n(1) + kbar - kbar^2`;
- mean metagraph eigenvalue `= d - 2 kbar`; spectral variance of eigenvalues `= 4 Var(k)`.
VERIFIED (reconstructed `P_n`, n=3..7): `Var(k) = 1, 2.19, 4.06, 5.25, 6.20`; **mean eigenvalue `=
1, 1.5, 1.33, 1.07, 0.72 -> 0`**. The positive-but-vanishing mean is the SIGNED-ACTION ASYMMETRY: the
hyperoctahedral (bit-flip) action breaks the hypercube `k <-> d-k` symmetry, so `mult(k) != mult(d-k)` and
the spectrum is slightly right-shifted, washing out as `n -> infty`. All computable from `n!` permutations
-- past the `2^{C(n,2)}` enumeration wall. Companion: `CV(H) = sd(H)/mean(H)` over iso classes
(`0.50, 0.47, 0.62, 0.57` for n=3..6), the dispersion of the Hamiltonian-path count.

## (2) The Burnside <-> Siegel moment analogy
`P_n` is a **Burnside average over `S_n`** of a per-permutation product -- structurally the FINITE analog
of the **Siegel transform** `f -> hat f(L) = sum_{v in L} f(v)`, averaged over the space of unimodular
lattices `SL_n(Z)\SL_n(R)`:
| | metagraph (Burnside / `S_n`) | Siegel transform (`SL_n(Z)`) |
|---|---|---|
| 1st moment | `P_n(1) = A000568` (count of classes) | Siegel's formula `int hat f = int f` (mean count) |
| signed/2nd | `P_n(-1) = SC` (antipodal Euler #); `Var(k)` (spectral 2nd moment) | Rogers/Schmidt 2nd moment `int hat f^2` (variance) |
| restriction | the bit-flip (hyperoctahedral) signs | **CONGRUENCE conditions** (arXiv:2507.05905) |
Both are "average over a group-quotient, read off moments." The paper (Han-Lee, July 2025) is exactly the
**second-moment Siegel transform WITH CONGRUENCE conditions** in dim 2, giving a Schmidt counting theorem
and a quantitative Khintchine theorem under congruence restrictions on the lattice vectors `(p,q)`.

## (3) The concrete LRC application
The LRC covering floor (THM-579, the gatekeeper) is `R' >= 1 - CV(N_R) sqrt((1-m_Q)/m_Q)`, holding when
`CV(N_R)^2 < m_Q/(1-m_Q)`, where `N_R` is the **14-sheet count** of the small part `R`. The exact identity
`sum_{N!=0} |chat(14N)|^2 = Var(N_R)/196` makes the gatekeeper a **second moment** (variance) of a
**lattice-point count restricted to `14Z`** -- i.e. a CONGRUENCE (covering / divisibility) restriction.
This is precisely the object arXiv:2507.05905 computes: a congruence-conditioned second-moment Siegel
formula. **So the paper supplies the natural machinery to bound `CV(N_R)^2` uniformly** (the open piece of
THM-579), the LRC sheet count being a congruence-restricted primitive-lattice-point count. The metagraph's
`CV(H)` (bounded `~0.5-0.6` over iso classes) is the finite combinatorial testbed -- a Burnside variance
that does not blow up, the analog of the Siegel/Khintchine bound staying finite.

## Extension menu (metagraph, for future sessions)
From the S19 survey, cheap+useful new metagraph objects: (a) the **per-level metagraph** `Pi_n(k)` (the
induced subgraph at level `k`, `V=mult(k)`); (b) the **H-gradient flow poset** (classes ordered by `H`, the
metagraph analog of the LRC `M`-gradient on covering sets); (c) the **bipartite bridge graph** `G_n--E_n`
(tournament<->even-graph duality as one graph); (d) the **multivariate GF** `G_n(x,y,z)=sum x^k y^{c3} z^H`
(moment extraction). The H-gradient poset is the direct metagraph<->LRC link: `M(S)` on covering sets is
the LRC's `H`, the tightest covering set (HYP-3551, `14/183`) its extremal, and the spectral-moment / CV
machinery transfers.

## What it buys
A new computable metagraph invariant (spectral moments, mean-eigenvalue asymmetry) that extends THM-587 for
free; a clean analogy placing the metagraph's Burnside averaging as the finite shadow of the Siegel
transform; and a concrete, actionable LRC direction -- use the paper's congruence-Siegel 2nd moment to
attack THM-579's `CV(N_R)^2` gatekeeper, the covering floor's open piece.
