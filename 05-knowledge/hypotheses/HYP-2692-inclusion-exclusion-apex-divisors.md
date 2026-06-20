---
id: HYP-2692
title: Inclusion-exclusion over the apex prime's divisors (7, 2, 3, 6, runners) gives one arithmetic skeleton — (Z/7)*=Z/2×Z/3, Gauss sum √−7, Gaussian periods, Q(√−7) — but the multiplicative structure WASHES OUT on the coverage p0; the proof gap is the archimedean height-weighted summed two-far residual R_2, indexed (not bounded) by the Q(√−7) resonance classes
status: OPEN (comprehensive reframing + redirect). Verified pieces: Gauss sum √−7, Gaussian periods (−1±√−7)/2, sector-Bonferroni=THM-534, danger-sieve measure-zero at extremal, Chebyshev bias ~70% (not universal), characters-of-p0 vanish. The redirect (R_2 height-weighted bound = the lever) is the synthesis target.
source: mac-mini-2026-06-20-S5
related:
  - THM-534   # coverage sieve = moment-LP (the working frame, closes k=8,9,10)
  - THM-548   # boundary-value / two-far curvature I_B, Phi_2
  - THM-551   # apex order-truncation
  - HYP-2657  # QR/Gauss D7 kernel, Q(sqrt-7) Galois structure
  - HYP-2684  # codex BV-Fourier cross-scale decorrelation (= the R_2 archimedean lever)
  - HYP-2681  # cube-root Eisenstein (the C_3 half of Z/2xZ/3)
  - HYP-2693  # six-sector Bonferroni4 high-tail gate, complementary final-row target
  - OPEN-Q-108
---

# HYP-2692 — Inclusion-exclusion over the apex prime's divisors

## The comprehensive view (one skeleton, refracted by 6=2·3)
"Inclusion-exclusion over N" for the LRC is a lattice of choices indexed by the apex prime 7:
- **N=7 sectors:** `p0=Σ_S(−1)^|S|m_S` = the moment-LP (THM-534). Plain Bonferroni truncation FAILS
  (level-4 upper `0.55>cap_9`; only full level-6 = p0); the OPTIMAL LP closes k=8,9,10.
- **N=2 quadratic χ:** QR{1,2,4}/NQR{3,5,6} = Legendre mod 7, Gauss sum `√−7` (verified), the
  reflection `x↦−x` (correction reality, HYP-2657), and a Chebyshev bias (NQR emptier ~70%, NOT
  universal — counterexample E=(0,2,3,4,5,13,18,21)).
- **N=3 cube root C_3:** sectors `(0)(1 2 4)(3 5 6)`; cube roots of unity mod 7 = {1,2,4} = QR.
- **KEYSTONE:** the C_3 orbit-sum of 7th roots = Gaussian period `(−1+χ√−7)/2` (period poly
  `x²+x+2`, disc `−7`). The 3-fold sum produces the 2-fold √−7 structure; the correction's C_3
  partial trace lands in `Q(√−7)` (class number 1), the apex prime's quadratic field.
- **N=runners (danger sieve):** `L(S)=Σ_T(−1)^|T|meas{T in danger}` — DEAD: `L=0` at the tight
  cluster {1..13} (witness τ=1/14 is measure-zero), the structural reason the coverage reformulation
  was needed.

## The verdict (the actionable redirect)
The arithmetic (cube root, χ, Gaussian periods, `Q(√−7)`) is the **organizing skeleton of the
correction's kernel D7** — but it **washes out on `p0`**: the dilation is not a measure symmetry, the
nontrivial multiplicative characters of the coverage vanish, the QR/NQR bias is archimedean (inside
`χ_0`). So inclusion-exclusion over the arithmetic N does NOT bound `p0`. It tells us the obstruction
is **archimedean**, and two reframings (Z/2×Z/3 character + far-order grading) converge on the single
remaining lever:
> **TARGET LEMMA:** a uniform height-weighted bound on the SUMMED two-far residual
> `R_2(B,F) = Σ_{f<f'∈F}[I_B(f,f') − Φ_2(B)]` — per-pair size is `~1/7^3` (apex hierarchy), but the
> multiplicity `C(r,2)` of *simultaneous* resonances is what must be controlled. The resonant pairs
> are classified by relation residue mod 7, whose C_3 orbits land in `Q(√−7)`; the arithmetic INDEXES
> the resonance classes the height-weighting sums over, while the smallness stays a signed archimedean
> cancellation (5× below absolute). This is exactly codex's HYP-2684 (BV-Fourier cross-scale) plus the
> signed `d≥2` relation-lattice bound.

## Next
1. Bound `R_2` (summed, not per-packet): height-weight by the relation rank/residue; use the
   `Q(√−7)` orbit structure to count simultaneous resonances; close via the signed (Poisson/theta)
   `d≥2` estimate keeping the 7-vanishing.
2. (sub) Decide if the Chebyshev bias has a clean averaged statement (log-density NQR-empty > QR-empty)
   even though it is not a per-cluster inequality — a `Q(√−7)`-flavored bias theorem.


## REFINEMENT (verified, mac-mini-S5): the lever is the LEADING-order residual, not R_2 specifically
Computing the summed order-s residuals `R_s = Σ_{|S|=s}[Δ_S(B)−Φ_s(B)]` for true-wide rows shows the
DOMINANT residual is the LEADING order `s_0 = max(1, 7−|B|)` (THM-551 order-truncation made visible):
- dangerous true-wide LEADER `B=(0,4,6,8,10,12,14)` (|B|=7): `R_1=0.174` DOMINATES, `R_2=−0.010` tiny,
  `R_3=0`. So the leader's danger is the ONE-far residual of barely-far {15,16} — THM-546/547 territory
  (the (6/49) bound + simultaneous peel + finite check), the BEST-understood piece. Good news.
- sparse core `B=(0,1,2,3)` (|B|=4): `R_1=R_2=0` (orders <3 vanish), `R_3=0.135` leads.
`R_tot = p0−P_r` stays within the cap margin in every tested row (0.13–0.16 < margin 0.17–0.30).
CORRECTED TARGET: bound the summed **leading-order** residual `R_{s_0}`, `s_0=max(1,7−|B|)`. For the
dangerous |B|=7 rows that is `R_1` (collar machinery, barely-far); for sparse cores it is `R_3+`.

## CROSS-LINK (codex-S59): high-tail final-row gate

HYP-2693/THM-556 gives the complementary sector-sieve form of the same split:
for six inner sectors, level-4 Bonferroni collapses exactly to
`U4=p0+p5+5p6`.  After the leading residual and finite-template branches are
routed, the true-wide final row should pass the high-tail gate
`p0+p5+5p6<=cap`.
