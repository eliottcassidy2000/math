# Handoff: regular-parent symmetry growth

**Current status:** R1--R6 in [PROOFS.md](PROOFS.md) are proved deductions;
the seed and small-order premises have an exact executable audit.
35,661 explicit gates pass in the standard-library verifier.
The old palette recurrences and weighted transitive formula are inherited.

## Do not restart

Completed:

- Regular-card set has size at most two and is pointwise fixed.
- Aut(T+v)=Aut(T)_S for every regular T.
- Exact unmarked extension inventory: subset orbits minus vertex orbits.
- Canonical augmentation rejects precisely N_T^+(w) union {w}.
- Labelled count (n+1)(2^n-n)L_n and both deletion multiplicities.
- Ternary unmarked stabilizer polynomial and child counts through parent81.
- Growing-color Poisson residual groups for C3, R5, P7 towers, with a
  rational higher-node bound and explicit triangle total-variation bound.
- 81/27 hostile pair: exact regular-color formulas and 1620/1728 binary cuts.
- Explicit parity attachments yield rigid P7-free and P11-free families
  at sizes 5^k+1 and 7^k+1, both with transitive maximum 3^k+1.
- Guard hostile: a rigid ten-vertex extension of C3[C3] creates R5.

## Next decisive tests, in priority order

1. **Integrate R2 into a real enumerator.** Identify the current fastest
   canonical augmentation path (Session87 is the inherited comparison
   point, not the old projective prototype). Use regular-parent detection
   to bypass child isomorphism work only in the proved sector. Measure
   wall time and duplicate counts against the same engine and universe.
   Keep nonregular parents on the existing path. R2 alone is not a
   benchmarked full-A000568 speedup.

2. **Beyond regularity: classify intrinsic deletion-root sets.**
   R1 uses a score constraint to force at most two roots, then odd group
   order fixes them. For another parent property P, calculate the possible
   Aut(U)-action on {v: U-v has P}. A set of size three may admit a C3
   action, so cardinality alone stops working there. Find a minimal
   nonregular parent whose attachment moves the new vertex; the transitive
   two-parent/C3 child is the compulsory hostile.

3. **Nonuniform palettes and arithmetic ports.**
   Replace uniform Q colors by a measured feasible palette. The upper
   proof uses the probability of child-orbit coincidences, not merely
   entropy or rank. Derive the relevant collision moments with source and
   clock retained, then combine with an actual simultaneous finite-window
   arithmetic bound. Do not infer a growing-n theorem from fixed-pattern
   equidistribution.

4. **Weighted EH prime-atom certificate.**
   R6 shows that removing automorphisms cannot remove the upper exponent
   barriers. The remaining positive target is the exact weighted
   transitive-set inequality over prime atoms, as specified in
   tournament_open_problem_routes_20261029/PROOFS.md. Start with a complete
   rational polytope certificate for one base; equal weights alone fail
   to certify arbitrary substitution weights (the inherited R5 triangle
   support is the hostile).

5. **Full tails rather than a bounded-observable Poisson limit.**
   Establish whether expected residual group order, or a specified
   return-weight functional, is uniformly integrable at the critical
   scale. Rare upper symmetries can carry large weight. Preserve the
   81/27 control and the inherited LRC star-family Gram-rank obstruction.

   Incoming commit `11f85d93cc` adds
   [source-weighted finite-depth reset tails](../pinned_discrepancy_reset_tails_20261029/FINITE_DEPTH.md).
   It reports a lower tail coefficient retained at every finite reset depth,
   and the actual source27 as a finite-scale hostile. This is a follow-up
   source, not a proof dependency audited here. Compare its two-parameter
   limit issue with the separate growing-palette and weighted-tail limits
   in R4; do not replace conditioned source data by a fresh random law.

## Scope to preserve

Regular parents here have odd order. The one-port theorem is unmarked;
the multiport probability theorem fixes every port individually.
Literal tournament, converse-extended, and switching-class groups differ.
R5's 81/27 counts concern the core stabilizer, not unmarked13 children.
The C3 polynomial is inherited, not a new canon claim.
No Collatz convergence, LRC(14), Sha-finiteness, Fermat-prime finiteness,
or general Erdős–Hajnal theorem has been proved.
