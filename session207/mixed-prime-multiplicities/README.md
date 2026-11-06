# Infinite mixed-prime multiplicity and metagraph laws

For fixed b,c>=0 with k=b+c>0, this packet gives the first two eventual polynomial coefficients for ordinary tournament-class counts with H=3^b5^c, the first-order probability of a non-prime-stack core, and concentration on degree5^c in the specifically inherited paired-black converse-merged metagraph Q. See PROOF.md for exact definitions and attribution. It does not give an exact first-order inverse-size formula for this infinite family or assert global count completeness.

The leading coefficient, fixed-count rationality, and qualitative concentration are inherited. The explicit next coefficient, first-order structural defect, and mixed-prime degree deduction are separated from those inputs.

Reproduce author checks, using ordinary Python with assertions enabled:

    python3 check.py > /tmp/session207-mixed-check.out
    cmp check.out /tmp/session207-mixed-check.out

inherited_controls.json is a selected exact extract of the reviewed complete small strong inventories and full core-word formulas; PROOF.md records the repository source and remotely verified blob. Independent mathematical and packaging reviews and independent inventory verification are included separately. manifest.json records every other release file's SHA-256 digest.

Reproduce the separate exhaustive strong-inventory audit through order seven:

    g++ -O3 -std=c++17 independent_inventory.cpp -o /tmp/session207-mixed-inventory
    /tmp/session207-mixed-inventory > /tmp/session207-mixed-inventory.out
    cmp independent_inventory.out /tmp/session207-mixed-inventory.out

This independent program has no external file dependencies. Its full permutation-orbit representatives and endpoint DP verify the five required small strong inventories; the inherited Busch cutoff proves no larger orders contribute to them.
