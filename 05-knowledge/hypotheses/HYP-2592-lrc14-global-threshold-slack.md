# HYP-2592: LRC(14) S3 Global-Threshold Slack Ladder

**Status:** OPEN proof program / exact case-bank evidence.

The density floor in the S3 residual should be parameterized:

`rho_alpha(P,E) = meas(G_P cap {x : maxgap({e*x:e in E}) > alpha})`.

The recent via-max counterexamples show that `alpha=2/7` has no uniform
positive floor.  The global LRC witness only needs `alpha=1/7`.  The new
threshold-ladder scout suggests a stronger intermediate target may exist:
prove a uniform floor for some `alpha` with

`1/7 < alpha < 2/7`,

possibly as high as `alpha=1/4` or `alpha=4/15`, and then use the resulting
gap slack to absorb the finite Weyl/CRT placement error.

## Evidence

The exact runner `04-computation/lrc14_global_threshold_ladder_codex.py`
scanned `132005` cases: all consecutive clusters, a structured adversarial
bank of perforated near-AP and relation-lattice stress shapes, and `500`
random bounded-spread cases.

- `alpha=1/7`, `1/6`, and `3/14`: zero exact failures; aggregate minimum
  `14249/252252`, attained at the trivial `k=3` G_P floor.
- `alpha=1/4`: zero exact failures; aggregate minimum `1/140`, at
  `P=(1,2,3)`, `E=(0,2,3,4,5,6,7,8,9,10)`.
- `alpha=4/15`: zero exact failures; aggregate minimum `2/525`, at
  `P=(1,2,3,11)`, `E=(0,2,3,4,5,6,7,8,10)`.
- `alpha=2/7`: `54` exact zeros, starting with
  `P=(1,2,3,6,12,13)`, `E=(0,2,3,4,5,6,8)`.

The explicit anti-correlation family is therefore not global-threshold
dangerous.  It binds at or very near the via-max boundary:

- The four named `rho_{2/7}=0` cases keep positive mass at `alpha=1/4` and
  `alpha=4/15`.
- The exact critical probe has floor `51/182 ~= 0.28022`, still above `4/15`,
  among the minima and zero examples checked.

## Proof Target

Replace the ambiguous statement `rho*(P,E) >= c0` by a two-level target.

1. **Slack reservoir:** find an explicit `alpha0 > 1/7` and `c0 > 0` such that
   `rho_alpha0(P,E) >= c0` for all admissible S3 shapes.
2. **Finite placement:** prove Weyl/CRT grid placement cannot miss that
   reservoir once the cluster maximum is large; finite small maximum remains a
   bounded exact check.

If `alpha0=4/15` works, the available gap margin over the LRC threshold is
`4/15 - 1/7 = 13/105`, which is large enough to make the finite-placement
lemma much less knife-edge than the old via-max `2/7` route.

## Tournament Analysis

Vertices are proof obligations / thresholds:

`1/7`, `1/6`, `3/14`, `1/4`, `4/15`, `2/7`.

The pairwise observable is robustness over the exact case bank: fewer zeros,
then larger minimum density.  The tournament is transitive, with leader
`1/7` and sink `2/7` (`score_hist=0..5`, no directed 3-cycles, one
Hamiltonian path).  This quotient preserves the implication chain among
threshold proof obligations, but destroys interval topology and relation
lattice data.  The next quotient should use safe components or relation
patterns if the `4/15` target stalls.

## Caveat

This is not a proof of LRC(14), nor a proof of a uniform `4/15` floor.  It is
a sharper conjectural target created by separating the false via-max threshold
from the true global-witness threshold.

See also `HYP-2591`, `HYP-2590`, `OPEN-Q-108`, THM-527/528, and the stored
result `05-knowledge/results/lrc14_global_threshold_ladder_codex.out`.
