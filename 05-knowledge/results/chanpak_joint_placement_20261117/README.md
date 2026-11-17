# Chan–Pak joint placement: the source count and its physical repair

The problem is to choose the slot, absolute count and signed correction
from the same target while retaining the exact reflected +2.
[Zehavi's Chabauty–Kim preprint](https://arxiv.org/pdf/2611.18218)
supplied a methodological analogy (especially §§1.2, 2.3 and 3):
change the constraint to an operator whose residual is uniformly
small, and read its surviving fibers at the correct filtration.
No result from arithmetic geometry is invoked in the tournament proofs.

## What advances, and what is still open

- **Exact joint membership.** Given a prescribed absolute source
  and a signed displacement, the [fixed-core decoder](inverse/PROOF.md)
  finds all witnesses in polynomial bit time. On each dyadic
  scale, the entire reversal term is smaller than a lattice
  spacing in the eliminated variable, leaving at most one
  physical source candidate. The [changing-flank inverse](physical/INVERSE.md)
  retains this property with a second independently chosen physical
  scale. On the original near-balanced branch a
  [target-only inverse](graded/TARGET.md) can recover its correction
  as well, rather than being given its value.
- **Graded geometry.** The [fixed-core image](graded/PROOF.md)
  has an exact limiting band law and a right inverse in every
  admissible suffix cylinder. Across levels its leading band
  classifier is an octal digit test. With a fixed near-balance
  constraint it has recurring relative holes, even though all
  fine shifts are available separately.
- **A physical escape from those holes.** The [variable-height
  Farey flank](physical/PROOF.md) independently varies a strong
  physical block and the dyadic suffix, keeping the reflected
  +2 exact. It permits an exact correction in a calculated
  residue class. Its leading coefficient rotates irrationally
  on logarithmic scales and [breaks the recurring relative
  gaps](graded/ESCAPE.md), allowing an arbitrarily close relative
  source with bounded discrepancy from the requested even shift
  (constants depend on the tolerance). An elementary quantitative
  version has error O(N/log log N) with shift discrepancy
  exp(O(sqrt(log N))). These are joint promises, not a fine
  absolute bound.
- **A genuine hidden continuation coordinate.** [Exact collision
  fibers](arithmetic/PROOF.md) preserve the count, order, and sum
  but respond differently to a common one-arc continuation. The
  response difference is exactly 2053G; an unbounded dyadic
  subsequence has an exact three-term linear recurrence for
  the shared counts. Therefore quotienting away all but the
  scalar count and address sum is not sufficient to model repairs.

The larger strong-flank bank still has only O(sqrt(X) log X)
endpoints up to X in its full natural repair width. Removing
a macroscopic hole is **not** enough to hit each integer.
A promising next step is to let the exact decoder choose between
genuinely different continuation states at the same coarse
count, or to add a new physical degree of freedom rather than
make this already-invertible sparse grammar search faster.

## Replays

Run from this directory:

    python3 inverse/decoder.py
    python3 graded/verify.py
    python3 graded/target_decode.py
    python3 graded/escape_verify.py
    python3 arithmetic/collisions.py
    python3 arithmetic/dyadic.py
    python3 physical/verify.py
    python3 physical/inverse.py

All arithmetic comparisons are exact. The changed flank replay
also includes independent literal subset-DP checks on small
tournaments. Finite checks validate the formulas and
implementation; they do not substitute for the all-order proofs.
