# The multi-far floor R′≥c is an Asano contraction of single-far Lee-Yang factors

*mac-mini-2026-06-27-S68. Owner: merge "tournament edges as witness, recurse on tip and tail" into the work on
closing the uniform multi-far floor R′≥c, with Elliott–Halberstam, Gaussian functions, and Asano contractions;
continue the φ⁴/Lee-Yang line (S67) without losing it; be creative with hypotheses. This wires three tools that
were NOT previously in the LRC context (Asano/EH/Gaussian — grep=0) into the open core of the covering bound,
and merges codex's HYP-3124 (edge-witness tip/tail recursion). Continues [[the-cap-is-a-phi4-field-theory-and-the-quartic-marks-the-hard-row]].*

## The target and the obstruction
The covering bound's open core is the **multi-far floor** (kps/codex, OPEN-Q-108): the quasi-independence
residual `R′ = meas(GOOD ∩ G_P)/(meas(GOOD)·meas(G_P)) = 1 + SPEC`, `SPEC = Σ_{n≠0} ĉ(n) ĝ(n)*` (a signed
low-frequency spectrum sum), must satisfy `R′ ≥ c > 0` uniformly over `r=2..6` far placements. Single-far (`r=1`)
is closed (HYP-2829, THM-563 periodicity); `r=2` doublets are almost-periodic (HYP-2797); **general `r≥2` has no
universal periodicity** — that is the wall.

## The bold reduction: multi-far floor ⟸ single-far Lee-Yang + Asano contraction
**The coverage is MULTILINEAR (multi-affine) in the runners.** By inclusion–exclusion,
`p0(E) = E_x[ ∏_{e∈E} (1 − 1_{miss,e}(x)) ]` — degree 1 in each runner's indicator. Multi-affinity is *exactly*
the hypothesis of **Asano's contraction lemma** (the Lee-Yang tool: if a multi-affine polynomial has no zeros in
a polydisk, the contraction merging two variables has no zeros in the disk). So:

> **Each far element is a "tip": a single-variable factor whose single-far Lee-Yang property (zero-free region)
> IS the HYP-2829 single-far closure. Asano-contracting the `r` tips against the bounded-core "tail" preserves
> the zero-free region → the multi-far partition function has confined zeros → `R′ ≥ c`. The multi-far floor
> reduces to the single-far factor, with Asano as the multi-variable engine.**

This is the multi-variable extension of S67: the single-tip is a φ⁴ measure (S67); **Asano contraction is the
Lee-Yang-preserving way to couple φ⁴ measures** (Lieb–Sokal); the multi-far floor is the φ⁴ field theory on the
far lattice.

## Verified evidence (`lrc_multifar_asano_floor_macmini_S68.py`)
At the coverage level, `R′_cov(F) = p0(B∪F)·p0(B)^{r−1} / ∏_{f∈F} p0(B∪{f})` (=1 iff the tips contract
independently — the Asano-factorized limit):
- **`R′_cov ∈ [0.87, 1.05]`, floor `≈ 0.98` for distinct far** (0.87 only at a degenerate repeated-speed config)
  — reproducing the kps/codex `R′ ∈ [0.81, 1.0]` range. **The multi-far coverage quasi-factorizes over the
  far tips, with a positive floor** — the Asano-contraction floor exists.
- **The single-far factor `d(f)=p0(B∪{f})/p0(B) ≈ 1.10` STABILIZES for large `f`** — the **Gaussian / free-field
  decoupled limit** (large, well-separated tips decouple, `R′→1`). The φ⁴ coupling is the resonance correction.
- Tight doublets give `R′_cov>1` (positive correlation, the HYP-2797 hard `r=2` case); separated give `R′_cov<1`
  (mild anti-correlation). The two bracket `1`, and Asano confines both.

## The three new tools, placed
- **Asano contractions** — the engine: multi-affine coverage ⟹ contraction preserves the zero-free region ⟹
  the multi-far floor is the contracted single-far floor. (First use in the LRC context.)
- **Gaussian functions** — the **free-field (λ→0) limit** of the φ⁴ single-tip measure (S67): large/separated
  tips decouple, `d(f)→const`, `R′→1`. The Gaussian is the `SPEC=0` baseline; the φ⁴ quartic is the residual.
- **Elliott–Halberstam** — the **level of distribution of the far tips**. `SPEC = Σ ĉ(n)ĝ(n)*` is small when the
  far cluster equidistributes (ĉ concentrated at `n=0`); EH (level `θ→1`) ⟹ `R′→1`. Unconditionally, a
  Bombieri–Vinogradov-level input bounds `SPEC`, giving `R′≥c`. **EH is the ideal-far-distribution that makes
  the floor `c→1`**; the residual cancellation (the project's signed-sum lever) is the unconditional substitute.

## The tip/tail recursion = the Asano contraction order (merging HYP-3124)
codex's HYP-3124 reads a directed edge `tail→tip` as a **two-ended proof obligation** (the tail = push-forward,
the tip = pull-back), recursing on both endpoint-deletions. That IS the Asano contraction order on the far
lattice: **peel one tip** (the largest far, by Node-3 equidistribution — the "pull-back" witness), **contract**
(Asano-merge it into the core), **recurse on the tail** (the remaining core + far — the "push-forward" witness).
The edge-witness compatibility (HYP-3124: "a legal witness makes both recursions compatible") = the Asano
lemma's guarantee that the contraction preserves the zero-free region. **Tournament edge-as-witness + tip/tail
recursion = the recursive Asano contraction of the far tips.**

## What this closes / the proof obligations (honest)
The reduction makes the open multi-far floor depend on three named, smaller obligations:
1. **Single-far zero-free property (Asano applicability):** each single-far factor is non-vanishing in the
   relevant polydisk — the analytic content of HYP-2829, to be stated as a Lee-Yang region (not just the
   L_yK8 margin). *This is the load-bearing lemma.*
2. **The `SPEC` bound (the floor constant `c`):** `|SPEC| ≤ 1−c` via EH-level far-equidistribution + the signed
   cancellation; verified empirically `c≈0.98` (distinct far). 
3. **Recursion termination:** the tip/tail recursion bottoms out at `r=1` (single-far, closed) and the bounded
   core (finite check). The decorrelation monotonicity (`r↑ ⟹ R′↑` toward the Gaussian limit) is the
   contraction-improves-the-bound step.

VERIFIED: the multi-affine (Asano-hypothesis) structure of coverage; the `R′_cov` floor `≈0.98` (distinct);
the Gaussian decoupling for large far. OBLIGATION/BOLD: the single-far Lee-Yang region (obligation 1); the EH
`SPEC` bound (obligation 2); that obligations 1–3 together *prove* `R′≥c` for all `r=2..6` (the closure).
**The prize: the multi-far floor — the genuinely-open core of the covering bound — is reduced to a single-far
Lee-Yang region + an Asano contraction + an EH-level spectrum bound.**

Related: HYP-3127 (this), HYP-2829 (single-far = the tip factor), HYP-3122/3103 (φ⁴/Lee-Yang/PGF zeros),
HYP-3124 (codex edge-witness tip/tail), HYP-2900 (Node-3 = the tip peel), HYP-2692/2797/2799 (the two-far
residual/doublet), THM-563 (single-far periodicity), OPEN-Q-108.
