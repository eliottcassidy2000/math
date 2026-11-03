# Audit and scope receipt

## Mathematical review

Two independent reviewers accepted the signed-root arithmetic-progression
proof, the centered-moment identity, the nine-cycle digit map, and the
no-affine-order-nine argument. The audit strengthened the initial obstruction
to the complete real degree-at-most-three numerical classification:

    f(k)=Ak^3+Bk^2+Ck+D, k=-4,...,4,
    magic iff B=0 and (A=0 or C=-15A).

The root supplied the positive cubic witness `k^3-15k`, preventing an invalid
extrapolation of the degree-two obstruction to every nonlinear assignment.

The lane audit produced a 945-system exact Gaussian classifier with independent
minor ranks and a separate normalized-`B` pair-graph elimination. A root-level
auditor independently wrote and ran another 945-system checker, without
importing either producer, obtaining the same nullity census and four nonodd
rays. The odd-plane 384-branch classifier was independently replayed as well.
All are finite exact algebraic classifications of **real coefficient spaces**,
not bounded integer coefficient searches.

## Relevant corrections incorporated

- Integer-square entries plus an integer-square common total is the impossible
  variant. Real-square entries would not support that arithmetic obstruction.
- A central **value** of zero need not have input **index** zero; all 945
  indexed systems are retained, including repetitions.
- The graph is a nine-cycle, not eight states; the final4 is the return.
- Carry suppression changes the permutation order and cannot be used as a
  faithful affine model on the ternary plane.
- Bare bicliques and unequal square line totals do not solve the distinct
  integer-square magic problem.

## Computation receipts

All checkers use standard-library exact integers/Fractions and explicit
exceptions rather than optimization-removable assertions. Positive controls
include ordinary shifted Lo Shu, a constant square grid, linear assignments,
and the cubic witness. Hostile controls include changed Lo Shu entries,
entrywise squared Lo Shu, suppressed carries, negative quadratic coefficients,
and the inherited unequal-square-total witness.

The main checker, 945-system checker, and normalized-nonodd checker have
byte-identical normal and optimized outputs. The portable manifest and all
seven bounded repository source pins are verified. A copied bundle is rerun
from outside both its directory and the repository before handoff.

`verify_hostile_controls.py` makes three mutations in temporary copies only:
source bytes, a duplicated cycle state, and a damaged cubic witness. For each,
both the affected exact checker and the portable manifest reject under normal
and optimized execution: twelve intended failures, with the failure reasons
recorded in `hostile_controls.json`. The unmodified bundle is not changed.

## Nonclaims

No unrestricted magic-square-of-squares solution or nonexistence proof is
claimed. The cubic witness is an ordinary numerical magic square, not a grid
of integer squares. No degree-four-or-higher classification, all-literature
novelty claim, formal proof-assistant build, or whole-repository re-audit is
claimed. The inherited carry principle is explicitly credited, not presented
as a new general extension theorem. The literal prior cycle's historical file
was not located in the bounded source recovery; its exact supplied word is
checked directly instead.

No Git edits or pushes were performed by this lane; publication belongs to
the parent session. The source repository remained clean at lane closeout.
