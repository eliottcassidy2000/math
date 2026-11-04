# Independent review of the information packet

Session 176. Reviewed by the sibling spectrum worker without importing author counting code.

## Mathematical review

No substantive defect found.

1. **Context expansion:** maximal module runs determine an unambiguous skeleton and an ordered path cover. Every interior return segment starts at an out-neighbor of the hole and ends at an in-neighbor, which are disjoint; its length is at least two. This proves the formal bandwidth bound for arbitrary exterior orientations.
2. **Triangular probes:** at maximum run count every exterior vertex belongs to a length-two return segment. The condition i≤π(i) forces π to be the identity; ordering the resulting r disjoint pairs gives exactly r!. The inversion divides by this nonzero diagonal, with exact divisibility on actual observations.
3. **Observational equivalence and nesting:** equal q implies equal Hamilton counts in every uniform context by the expansion. The actual probes give the converse. A uniform context around a composite remains uniform relative to its original module, so applying probes after substitution proves equality of the entire resulting profile. This argument does not extend to nonuniform attachments; the proof correctly excludes them.
4. **Formal versus realized rank:** the full formal response map is injective, hence its restriction to W_n has rank dim W_n. The paper does not incorrectly identify this with n. Distinct lowest powers in the ordered triangle/singleton family prove actual linear independence and unbounded realized dimension.
5. **Examples:** the equal-H/different-q five-vertex examples and the labelled nonuniform counterexample are correctly oriented. The cyclic two-vertex wrapper is strong for every nonempty input module, even if the input is not strong.

The lower bound is explicitly for linear sufficient interface states, not arbitrary packed or nonlinear encodings. The probe-size optimality is explicitly about visibility of a formal coordinate, not an unsupported lower bound on the realizable variety. No spectrum completeness conclusion is drawn.

## Replays

The author's `check.py` and `audit.py` were each rerun normally and with `python3 -O`. Each normal/optimized pair agreed exactly.

The independent `reviewer.py` uses a different path-cover method: for each vertex permutation with d backward adjacent edges, every valid cut set must include those d positions, so its contribution to q(k) is binom(n−1−d,k−1−d). This counts all profiles through order five without subset path-cover DP. It computes realized ranks using evaluations of F at 1,...,n rather than the q coordinates.

Additional independent checks cover 48 arbitrary actual uniform contexts, twelve two-level nested equal-profile substitutions, and direct exterior-permutation/gap enumeration of all five displayed probe rows through r=4. The latter checks the row (117,684,908,336,24), which the original literal audit stopped short of checking.

The reviewer normal and optimized results are saved as `reviewer.json` and `reviewer.optimized.json`. Finite checks support the proofs and do not replace the all-order arguments.
