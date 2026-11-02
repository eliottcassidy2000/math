# Boolean tournament inverses and the seventeen seventy one connection

**2026-11-02. PROVED scoped lemmas / FINITE-EXACT certificates;
same-agent audit with independent computational paths.**

Weighted tournament values are expectations. A regular algebraic point
with value 7 or 21 therefore does not produce a tournament with that
path count. The useful repair is to reduce the squared error using the
Boolean relations before evaluating it:
\[
 {\cal E}_h(x)=B[(H-h)^2](x)
             =\mathbb E(H-h)^2.
\]
If h is odd and this energy is below four, successive coordinate
fixings construct an actual tournament with count h, preserving any
fixed path arcs. Finding such a point is a separate problem.

The session gives these exact results:

- Single-edge means 7 and 21, with derivative four and energy four,
  arise between counts 5/9 and 19/23. Exact rational examples also
  lie in the full open probability cube. Four is the sharp universal
  energy lower bound for these forbidden targets.
- Exact energy has strict non-solution local minima. For target 17,
  the complete six-vertex sublevel graph at energy at most four has
  768 isolated non-solutions. Mask 72 supplies a short witness.
- The complete seven-vertex fixed-path chart has twelve literal
  complement lines between H17 and H71. The same ordinary tournament
  at tiles 24511 and 31 has complementary counts 71 and 33, proving
  that its selected path port is indispensable.
- Under the restricted command C alone, the pair (H,H after C) is the
  minimal exact observer. At order seven it has 3,008 states.
  Re-marking and attachments are additional obligations.
- Independent copies of the 5/9 source give expected ordered-join
  count 49; shared copies give 53. Their difference is covariance.
- The old false Boolean partition and actual 81/27 symmetry pair
  survive the new tests. The arithmetic words (1,7) and (7,1) retain
  their distinct source ports 483 and 469 despite the same target 17.

Read [THEORY.md](THEORY.md) for proofs and explicit witnesses;
[RECOVERY.md](RECOVERY.md) for inherited results, sources, the concept
board and scope boundaries; and [HANDOFF.md](HANDOFF.md) for the next
decisive research question.

## Reproduction

From the repository root:

    python -B 05-knowledge/results/boolean_realization_71_17_20261102/verify.py
    python -O -B 05-knowledge/results/boolean_realization_71_17_20261102/verify.py --output 05-knowledge/results/boolean_realization_71_17_20261102/verification_optimized.json
    python -B 05-knowledge/results/boolean_realization_71_17_20261102/seal.py --check

The complete universes are all labelled tournaments at orders 3 through 6
and all 32,768 seven-vertex fixed-path charts. Permutation path counts
independently check the entire labelled census through order five and
each key witness. Weighted path products independently check the exact
interior means; Boolean polynomial reduction and direct expectations
agree. Normal and optimized Python outputs agree bytewise.
There are 1,173 explicit gates, in addition to the enumerated cases.

The two scouting programs and their outputs preserve the discovery
trail; verify.py reconstructs its own universes. seal.py records normalized
text hashes of local artifacts and the actual external code dependency.

The unrestricted spectrum and the smallest carrier closed under all
re-markings and attachments remain open. No new A000568 term is claimed.

