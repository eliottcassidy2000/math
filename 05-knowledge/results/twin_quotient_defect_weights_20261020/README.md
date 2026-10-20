# Twin-quotient defect weights

**Status: PROVED local formulas and sharp boundary; FINITE-EXACT controls.**
20 October 2026. All work remains on the graph/tournament carrier of the cited
repository results. No unrestricted odd-group equality is asserted.

Start with [`TWIN_QUOTIENT_DEFECTS.md`](TWIN_QUOTIENT_DEFECTS.md) for the
odd-module sign identity, weighted false-twin cycle-index formula, and failure
of the proposed unit-weight split/join cancellation. Continue to
[`SHARP_DEFECT_BOUNDARY.md`](SHARP_DEFECT_BOUNDARY.md) for the new sharp theorem.

## Principal conclusions

- Contracting independent modules of odd size preserves the edge-reversal
  character on partition-preserving automorphisms.
- The alternating twin-triple packing weight is exactly b=E_0-D, where D is a
  signed cycle-index evaluation on the weighted false-twin quotient.
- For every graph of order at most eleven, its defect b-e is -1, 0, or 1;
  globally odd graphs in that range have no negative defect.
- At twelve vertices precisely K_(6,6) and the independent 3-blow-up of P4
  have defect magnitude greater than one: respectively -2 and +2.
- The unique first odd graph with negative defect is 2K_(3,3), with defect -1.
  Odd graphs have both positive and negative defects of unbounded magnitude.

These results correct an all-order extrapolation from the valid small census;
they do not retract its data. The inherited prime-local theorem already makes
the global defect sum zero. The remaining target is an explicit structural
involution on decorated packing orbits, not a new proof that the sum vanishes
or a claim that no decorated pairing exists.

## Files and reproduction

The detailed reproduction commands and exact scopes are in section 8 of the
main note. Python 3's standard library is sufficient. Three paths check the
weighted formula, literal labelled packings, and the sharp-boundary controls;
all normal and optimized outputs agree. See
[`VALIDATION_RECEIPT.md`](VALIDATION_RECEIPT.md) and `SHA256SUMS.txt`.

No claim of literature priority is made. The b-carrier is inherited from
`01-canon/theorems/THM-4723-sylow-profile-form-of-hyp-9340-and-the-defect-reduction-for-the-three-cycle-type.md`;
the formula and counterexamples here are proved directly.
