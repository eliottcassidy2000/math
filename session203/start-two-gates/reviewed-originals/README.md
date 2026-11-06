# Start-two physical profile compiler

**Candidate awaiting independent audit.** Two fixed graph gates, together with the inherited O/Z/J grammar, give a complete physical normal form for every rooted tournament having exactly two root-starting Hamilton paths. The nine-coordinate state is closed, and every represented B has a witness of order at most 3+floor(log2 B).

The universal gate calculation rests on an all-input one-port path-algebra lemma. Exact evaluation on four linearly independent physical basis graphs determines a context matrix; finite samples are not used as an infinite-coverage premise.

Useful exact target distinction, subject to the proof audit:

- B=137 has a seven-vertex start-two source.
- B=139 has no start-two source at any order.
- There are 302 attained start-two bases through 997.

These are statements about Gamma(B), not exclusions from the ordinary Hamilton-path spectrum. The unit-start fiber has the opposite behavior at 137 and 139. The staggered condition F(8,137), and every unconditional eventual odd threshold, remain OPEN.

## Reproduce

    python3 verify.py > NORMAL.json
    python3 -O verify.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json
    python3 compiler.py 137
    python3 compiler.py 139
    python3 compiler.py 997 --inventory

The verifier checks 5,595 literal children of 1,119 rooted inputs, reconstructs the universal matrices by exact rational elimination, compares an independent canonical start-two graph inventory through order 8, and checks positive, absence, singleton, and resource-limit controls. See PROOF.md for all-order justifications and scope.

Only one K1/K2 gate is allowed in the complete start-two compiler. Repeated gates are valid physical constructions with larger start counts, but no all-start classification or all-target coverage is claimed for them. All computations use exact integers. A resource-limit result is INCOMPLETE, never an absence conclusion.
