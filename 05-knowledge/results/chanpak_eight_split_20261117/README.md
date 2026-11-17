# The split 8n = n + 7n: paired constructions, factor descent, and phase

This session treats the split as a request for **two actual tournaments**, not as a Hamiltonian count of 8n. Rédei parity rules out the latter. It proves a literal two-vertex sevenfold lift with an arbitrary variable tournament inside it, bounds exactly how that lift can contact a prime-count source in the one-exit architecture, and gives a separate exact Lonely Runner sheet theorem. Selection at every prescribed prime remains open.

## What composite quotients give us

Write S for the Hamiltonian-path spectrum. Ordinal composition gives multiplication: each path visits the ordered strong components in sequence. In particular 7ab is in S whenever 7a and b are in S. This permits descent on the quotient but never cancellation of a physical witness.

Assume that all ordinary odd primes other than 7 are in S. Then for odd n coprime to 7, to prove 7n in S for all n other than 1 and 3 it suffices to prove 7p in S for primes p≥5, p≠7: if n has such a factor, send it to that branch and put the other factors in ordinal components. Otherwise n is a power of 3; 63 gives the start at exponent 2, and ordinal triangles give every higher exponent. The missing quotients are 1 (the unit, not prime) and 3, corresponding to 7 and 21. The quotient 7 is separate from this exactly-one-factor-of-seven discussion; 49 and 343, which already have witnesses, supply all valuations ≥2 by adding exponents 2 and 3.

In particular if the conjecture fails on the 7-coprime quotient lane, under the ordinary prime assumption its **least** unexpected missing quotient is prime. Finding more composite examples doesn't change that obligation. For prime p≠3,7 any witness for 7p has one strong component of count 7p and only singleton SCCs outside it: separating the two prime factors would require a component of count 7. A p-witness has the analogous unique count-p component.

- [Literal sevenfold lift and exact contact constraints](LITERAL_LIFT.md)
- [Why all ordinal combinations of one-apex sources in the old readers still miss almost all primes](THIN_READERS.md)
- [A paired Lonely Runner sheet theorem with sharp 3/16 margin](RUNNER_SHEETS.md)
- [Literal subset-DP replay](replay.py) and [one-apex formula replay](replay_apex.py)

The structural next question is now a concrete landing test. Can every still-uncovered p be realized by a strong C with marked set Q, |Q|≤4, satisfying the cover identity below, or by a different architecture where that identity is not required? One-exit TT padding of the return block offers no extra degrees of freedom. An analogous count and residue are not enough without the marked cover data in the same witness.
