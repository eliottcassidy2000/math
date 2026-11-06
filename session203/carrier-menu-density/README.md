# Density boundary for the four paired carrier orbits

**PROVED in the stated scope.** This packet bounds only the four selected O/Z carrier orbits from `../staggered-paired-carriers/`. Their upper relative density among positive odd bases is at most 693609731/46994283721, strictly below 1.5 percent. Consequently they cannot contain an odd tail.

This is not a density bound for all physical profiles, all seed banks, or grammars using J or other constructions. F(8,137) and the eventual threshold 23 remain OPEN.

Read PROOF.md. The reviewed proof, its carrier dependency and its separate audit are retained under audit/. The promoted proof changes the status line and corrects two relative citations to the sibling carrier proof; see STATUS_PROMOTION.diff.

Reproduce:

    python3 check_density.py > replay.json
    cmp replay.json DENSITY_NORMAL.json
    cd audit
    python3 check.py > replay.json
    cmp replay.json NORMAL.json

The geometric word-count argument is the all-order proof. The finite checks validate its exact arithmetic and boundaries.
