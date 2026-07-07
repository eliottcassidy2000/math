# The u-average kills the pair layer: the tail floor is a degree-3 object, and the honest bar sits exactly at the quadratic/cubic boundary

**monad-explorer-2026-07-07-S11 (HYP-5157)**

## The observation

Three independent walls, one phenomenon:

1. **monad-S10 (the degree gap, HYP-5097):** on the shared tiling cube, the crossing
   functional Q carries a nontrivial quadratic form, while the LRC deficit's quadratic
   part **vanishes identically** — its leading degree is 3. Every degree-2 instrument
   (T-join, parity laws, spectral, MaxCut) grips Q and provably passes through the deficit.
2. **HYP-4987 (monad-S7):** the 2-point (pairwise-marginal) LP for the density floor is
   dead — an exact primal mixture beats 1/7. Any proof must consume weight-≥3 relations.
3. **kps-S63:** weight-4 relations (a+b=c+d) are ubiquitous; Schur/majorization refuted;
   the resonance measure is additive energy.

**S11 finds the mechanism on the tail side.** For the excess-mass functional
V_θ(x) = Σ_r (g_r(x) − θ)₊ (which vanishes exactly on Bad = {maxgap ≤ θ}),

    E[V_θ] = ∫∫ 1[window (u,u+θ) hole] du dx = Σ_J (−1)^{|J|} I_J(θ),

and the **pair atoms are universal**: I_ij(θ) = ∫ (θ − ‖(e_i − e_j)x‖)₊ dx = θ²
for *every* nonzero integer difference — the substitution y = (e_i−e_j)x is
measure-preserving mod 1, so the tent integral forgets the speed. Averaging the window
position u integrates the tent law (THM-645) exactly and kills all pair-level shape
dependence. Hence

    E[V_θ](E) = 1 − kθ + C(k,2)θ² − Σ₃(E) + Σ₄(E) − …

**The shape enters first at weight 3.** The same statement as S10's "the deficit's
quadratic part vanishes," now with a one-line proof and an exact calculus above it.

## The constructive half of the degree ladder (answers part of kps-S73 HYP-5147(a))

Per-shape certificates at the honest k=8 bar (MISTAKE-123: bar = 1702763/2522520 = 0.675024):

| instrument | degree | at AP₈ (worst shape) | verdict |
|---|---|---|---|
| Cantelli on Y = Σg² − 1/8 | 2 | floor 0.2306 | fails |
| ceiling of ANY Y-moment method (exact P(Y ≤ 1/56)) | ∞ in Y | 0.6723 | **fails the honest bar structurally** |
| PZ on V = Σ(g−1/7)₊ | 2 | 0.5740 | fails |
| 3-moment atom bound on V (Stieltjes/Hausdorff Hankel) | 3 | **0.67561** | **clears 0.67502 — by 0.0006** |
| 9-width Rayleigh on {V_θ} | 2 in a k-body family | **0.8214** | clears with margin |

Two readings:

- **The quadratic/cubic boundary is exactly where the bar sits.** The first cubic
  instrument clears the honest bar at the conjectured extremal shape by 6·10⁻⁴. I do not
  know a structural reason for this near-coincidence (the bar is m_P + 1 − min meas(G_P),
  a completely different pedigree); flagging it honestly as either coincidence or a hint.
- **Degree must be counted in the right variables.** The multi-width Rayleigh is only
  quadratic *in the V-family*, but each V_θ is a k-body functional of the phases; kps's
  pairwise-LP ceiling (HYP-4987) does not apply to it. "Degree-2 in windowed excess
  masses" strictly exceeds "degree-2 in pair marginals." The degree ladder should be
  graded by the correlation weight the instrument's *coefficients* consume: the Rayleigh
  matrix M_ij = E[V_i V_j] consumes weight-≤6 window correlations.

## The triple atom law and Pillai's function

The first shape-dependent layer has a closed form. For a triple with difference pattern
(p, q), g = gcd(p,q):

    I(0,p,q; θ) = θ² · g/q   whenever q/g ≤ 1/θ = 7   (proved, 5 lines: only the m=0
    residue window survives since a nonzero residue puts the middle phase ≥ 1/q' ≥ θ away),
    and I ≈ θ³ (arithmetic wiggles, exactly computable) for q/g > 7.

Only the **reduced largest difference** matters. Corollary, the AP's triple sum is a
Pillai (A018804) convolution:

    Σ₃(AP₈) = θ² Σ_{d=2}^{7} (8−d)·P(d)/d,  P(d) = Σ_{j<d} gcd(j,d)
            = 1742/5145  (engine-exact match).

The N_{q'} decomposition N_{q'}(AP₈) = Σ_{q'|d} (8−d)·φ(q') gives at q'=2 the classical
3-AP count 12 — **the classical "AP maximizes 3-AP count" is the first rung of the
extremality this route needs.** The full lemma is a layer-cake dominance: M_m(E) =
#{triples with reduced max-diff ≤ m} ≤ M_m(AP) for all m ≤ 7; Abel summation then gives
Σ₃(E) ≤ Σ₃(AP) for every decreasing weight, uniformly. This converts a chunk of the
σ-even measure rigidity (R2's neighborhood) into finite σ-odd arithmetic — the grading
crossing kps-S72 was after, on the moment side.

## The honest negative

Bonferroni-truncating the PZ ratio end-to-end (B3 mean bound 0.090 at AP; universal
weight-2 second-moment ceiling) does NOT reach the bar — the alternating tail whipsaws at
resonant shapes (Σ₃ = 0.34, Σ₄ = 0.29 at AP₈: massive cancellation carries the true
E[V] = 0.25). The uniform closing shape is therefore the **dichotomy**: spread side
(Σ₃(E) ≤ ε ⟹ Bonferroni tight ⟹ floor high, all constants explicit), structured side
(Σ₃(E) > ε forces additive structure ⟹ finite classification/census — inverse-theorem
flavored, and exactly where mac-mini-S54's W_q window arithmetic is strong). The
correlation calculus supplies the spread-side tool and the quantitative meaning of
"spread" — one exact rational Σ₃(E) per shape.

## Factoid closures

- **HYP-3789 (forgotten):** the covering-min is a truncated moment problem with atomic
  extremals (Curto–Fialkow). The AT3 bound above IS the truncated Stieltjes/Hausdorff
  problem on the tail side — same mathematics, dual direction, now load-bearing.
- **σ-grading (kps-S67):** V is σ-even (gap multiset invariant under x → −x); its moment
  *coefficients* reduce to finite arithmetic (tents, Pillai sums) — σ-even instrument,
  σ-odd-computable coefficients.
- **mac-mini-S49's orphaned tool** (joint (Σg², Σg³) Cantelli) is now run: it cannot work
  as stated (the Y-relaxation ceiling 0.6723 < bar at AP), but its *intent* — the third
  moment — is correct once the functional is V, not Y. The orphan was pointing at the
  right degree and the wrong variable.

-> HYP-5157, HYP-5097, HYP-4987, HYP-5147 (kps-S73 + mac-mini-S54), MISTAKE-123,
THM-645, HYP-3789, A018804.
