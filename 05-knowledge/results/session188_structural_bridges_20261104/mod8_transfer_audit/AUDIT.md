# Independent audit: even-order modulo-eight ambient transfer

PASS. The all-order derivation, physical lift, odd-order correction, and
primary citation check out. The exact proof and producer pins are in both
independent receipts.

## Mathematical review

The unordered-partition identity for q_k/k! is exact. Expanding each part's
Hamilton count modulo four retains at most one odd-cycle contribution.
Contracting all vertices of a fixed odd directed cycle to one atom counts
the relevant partitions by S(N-ell+1,k). This correctly counts cycles that
share vertices across different summands; no disjointness of two cycles is
assumed because only one-cycle terms survive at that stage.

For k=2, a proper odd cycle contributes an odd partition count, while a
spanning odd cycle contributes zero. For k=3, the collapsed atom count has
the parity of N; its Stirling number is even for even collapsed order and
odd for odd order at least three. The spanning-cycle boundary again gives
zero. These facts prove both the even formulas and the exact odd-order
Hamilton-cycle correction. The singleton conventions produce q_2=q_3=0.
Multiplication of the modulo-four expansion by 2 or 6 gives a valid
modulo-eight conclusion, and k! supplies the higher-cover divisibility.

For same-size even cores, only Delta q_1=Delta H and
Delta q_2=2 Delta H survive modulo eight. Thus the ambient multiplier is
c_1+2c_2, exactly H(Q[TT2]), an odd integer. This verifies preservation of
every nonzero modulo-eight residue, not only residue four. The physical
port lift makes the extended core a common module with the same exterior
quotient, so the actual response conclusion is justified. It does not
require coefficientwise integer monotonicity.

For odd cores, Delta q_2=Delta q_3=2 Delta H-4 Delta z, giving the stated
correction. Even Delta z suffices for transfer; without it a residue can
vanish while the exact integer gap remains nonzero. The supplied examples
are such residue failures, not full-response collisions.

## Primary source verification

Independently opened the linked detailed Grinberg--Stanley PDF. Its title
page is dated June 13, 2026. Theorem 7.1 on printed page119 states the
Hamilton-path count modulo-four refinement in terms of nontrivial odd
directed cycles. The source and version used by the producer are accurate:
https://www.cip.ifi.lmu.de/~grinberg/algebra/redeiberge-long.pdf .

## Independent exact controls

verify.py imports no producer code. It uses literal permutations and cut
choices for ordered path covers, literal cyclic permutations with physical
vertex0 fixed for Hamilton cycles, and a canonical unordered-partition
recurrence as a separate cross-check. Its Hamilton counts use independent
recursive endpoint memoization.

Normal and optimized receipts agree. The controls comprise:

* 1,143 profile cases: every labeled tournament through five vertices,
  plus declared seeded samples at orders six, seven and eight.
* 1,041 literal odd-order Hamilton-cycle counts.
* 18,240 literal ambient comparisons: every four-name core against TT4 in
  every rooted quotient through order four, including residues zero, two
  and four.
* Four physical cyclic port-core controls, including odd core orders with
  Hamilton-cycle differences -13 and -15 and exact ambient gaps -600 and
  -720. The even controls retain residue four.

physical_controls.json preserves the exact rows and profiles used in the
last set. Producer replay is checked separately. The finite controls do
not replace the all-order partition and substitution arguments.

Reproduce from this directory with the sibling mod8_transfer packet present:

    python verify.py > normal.json
    python -O verify.py > optimized.json
    cmp normal.json optimized.json
