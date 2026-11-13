# The first next-target neighborhood has no smaller positive defect

Status: FINITE-EXACT for the stated radius-two parent neighborhood, with exact all-height module certificates for its positive stratum. This is not a global obstruction to growing-defect returns.

Start with Q=(14,12,8,64,175,79,535,367,127,447), root9, from MODULE_MIGRATION.md. Its base count is293 and defect1. Append x beaten only by root9, leaving an 11-vertex parent. For every subset of at most two of its55 arcs (1+55+1485=1541 parents), test every sole predecessor of a final newborn y and retain every strong child with total H881. No upper bound on defect is used in this census.

There are exactly11 retained labeled outputs: seven have defect213 and four have defect0. In particular the minimum positive defect is213, identical to the unmodified two-birth output; there are no defects1 through212 in this specified neighborhood. Thus extending the former b1/b2 search to growing defects does not expose a smaller positive off-axis return here.

For all seven positive outputs, independently compile every single-vertex transitive inflation other than the final root or its unique predecessor:70 module choices. Exactly seven lifts match the desired full target H=405*2^k+71, one per positive output, all at module vertex5. Every one has the same defect

    b(k)=108*2^k-3.

The raw two-birth output is one of these lifts. The other six have different parent arc choices but no improvement in the defect or the target law. These are labeled outputs, not a claim of seven ordinary isomorphism classes. The defect-zero outputs are recorded separately in the receipt and were not subjected to a positive-defect module-lift claim.

The strongest survivor is therefore the already known growing-defect target family, with several neutral local realizations. A future repair must use a larger parent neighborhood, another source realization, a different module interface, or another operation. The finite result excludes none of these possibilities.

Reproduce:

    python -B 05-knowledge/results/cycle_opening_six_return_20261112/returns/growing_successor_probe.py

The receipt growing_successor_probe.json records every retained row tuple, flipped arc set, entry port, defect, and exact exponential coefficient dictionary. Runtime in the first run was about13 seconds; it is not used as a mathematical bound.
