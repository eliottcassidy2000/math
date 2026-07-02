# The sLRC-strength audit of the S87–S92 lemma stack

klein-2026-07-01-S92 (HYP-3839). Filter: BCS arXiv:2603.24784 — the SHIFTED Lonely Runner
Conjecture is FALSE from n = 5 (coloopless-zonotope counterexamples), so any lemma
quantified over all shifts at sLRC strength is false-risk; state floors as shift-averages
(Fubini) or shift-existence only (kps-S31's design filter). Classification:
**PINNED** (uses the common start/origin; safe, but transfers nothing to phased settings),
**SHIFT-SAFE** (holds verbatim for arbitrary phases), **SHIFT-AVERAGE** (holds after
averaging phases), **∀-SHIFT-RISK** (would assert sLRC-strength coverage impossibility —
forbidden by BCS).

| Lemma / result | Class | Reason |
|---|---|---|
| HYP-3834 gap-sum formula for Λ | SHIFT-SAFE | any union of intervals; phases irrelevant to the mechanics |
| HYP-3834 Farey profile / collapse constant 1666/6435 | PINNED | Farey fractions = the pinned danger centers |
| THM-597 collapse-rate law (argmax/unit-residue) | PINNED | witnesses at k/n use the common zero's residue rigidity |
| HYP-3841 tangent-ladder mechanics (value/slope/defect) | SHIFT-SAFE | piecewise-linearity + kink taxonomy hold for phased arcs |
| HYP-3844 d = 1 band emptiness (integrality) | **PINNED — FLAG** | phased crossing radii are (Δfrac + Δphase)/(w−v) ∈ ℝ: the band CAN be hit; the final window is defect-free only in the pinned setting. Any sLRC-side ladder must carry a phased defect budget. |
| HYP-3847 local deviation lemma | **SHIFT-SAFE** | phases rotate Fourier coefficients; the surviving modulus s = sin(2πr)/π is phase-invariant. One of the few sLRC-grade analytic floors — consistent with THM-598 (built for free phases). |
| HYP-3845 polygon/DMNR (Lean: PolygonPartitionDMNR) | SHIFT-SAFE | a ∀-shift NONEXISTENCE that is proved (exact covers, distinct moduli) — the depth of the shifted side, not a risk |
| HYP-3837 cap ladder; HYP-3838/THM-601 nest lemma | PINNED | the origin nest IS the pinned structure; phased d-folds are kps HYP-3954's torus averages instead |
| THM-594(B) pair law | PINNED (two-branch) | kps torus (ii): the c-average removes the arithmetic — the shift-average version is clean independence |
| THM-598 anti-covering floor | SHIFT-SAFE by design | proved for free phases; its dichotomy (frozen clusters tile) marks exactly where ∀-shift claims die |
| Odd-covering bridge (HYP-3853/3849 sliver) | consistent | pinned analogue trivial (t = 1/2); Erdős–Selfridge depth lives on the shifted side — matches kps's flag |

**Reading.** The pinned/shifted divide tracks the ORIGIN exactly: everything that uses
the common zero (nest, caps, unit-residue witnesses, band integrality) is pinned;
everything spectral or measure-theoretic (deviation lemma, DMNR, anti-covering,
gap-sum mechanics, ladder mechanics) survives shifts. Two actionable flags: (1) any
future use of the K = 0 final-window lemma inside a PHASED argument (e.g. renormalized
clusters carry effective phases!) must replace band-emptiness with an explicit phased
defect budget — the renormalization tower should double-check this at its recursion
step; (2) the deviation lemma (HYP-3847) is the designated sLRC-grade floor — route
phased window arguments through it, not through the pinned band lemma.
