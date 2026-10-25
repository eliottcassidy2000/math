# Session 55 density bridge

- PROVED: uniformly over odd prime powers, labelled distinct square-value magic arrays satisfy F(q)=q³/512+O(q^(5/2)); an explicit bound is in PROOF.md.
- PROVED: every prime has the square-content recurrence N_p(k)=P_p(k)+N_p(k−2), where P means at least one unit entry. Boundary levels are treated explicitly.
- SYNTHESIS: the already-proved dyadic and ternary ring densities, 1/504 and 27/728, have ratios 64/63 and 1728/91 to the odd-field limit. The source is a two-level square-content descent in a three-dimensional magic parameter space, not a field/ring identification.
- HOSTILE CONTROL: primitive need not mean all-unit beyond p=2,3; a valid mixed example modulo 29 is included.

Zero is allowed, entries are labelled and distinct, and square-root choices are not counted. No integer-square witness, MRB/modular-CM relation, or first-asymptotic novelty is claimed.

The sources are frozen copies with SHA256 provenance in source-pins.json. The packet is portable and standard-library-only. It reuses the source field arithmetic and does not claim a new independent field-engine audit.

Reproduce from any working directory:

    python3 -B path/to/session55-density-bridge/verify.py
    python3 -O -B path/to/session55-density-bridge/verify.py

To deliberately regenerate the certificate, add --write. The exact certificate includes all odd prime powers through 125, two-sided conic-slice error controls, bounded ring recurrence checks at primes 2,3,5,7,29, exact rational density ratios, and finite convergence checks. Analytic proofs establish the unbounded results separately.
