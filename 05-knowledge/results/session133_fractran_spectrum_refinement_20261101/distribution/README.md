# A 397-fraction stream of eleven tournament multiplicity sequences

Start `bounded21.fractran` at **5438**. At states divisible by **2719**, the
exponent of **2729** is the tournament order n; the eleven output register
primes are `program.json`'s x entries 0,6,...,60. They count ordinary
isomorphism classes with H=1,3,...,21 paths, respectively.

This is a running stream, not a halting computation. See [THEORY](THEORY.md)
for the exact all-order scope, classical and inherited dependencies, generic
positive matrix compiler, first-eligible proof, and exact microstep formula.
Run `python3 check.py` (also works under `python3 -O`). No third-party packages.

The inherited rational formulas are recorded in `inherited_formulas.json` as
an explicit prior-work control, not as a new discovery. Full graph-theoretic
inventory certification is in the sibling spectrum and core_reduction lanes.

**Refinement:** `tail21.fractran` uses only **67 fractions** for the same eleven
sequences starting at n=6. Its seed and output mapping are in
`tail_program.json`; all counts at this seed are proved initial values, not
an oracle. The full-prefix program remains available above.

**Further refinement:** `compact_tail21.fractran` uses **43 fractions**, sharing
registers for the two pairs of identical sequences and using an absent prime
for the two zero channels. Its eleven-entry `output_primes` map and exact seed
are in `compact_tail_program.json`. Readout remains a prime exponent.
