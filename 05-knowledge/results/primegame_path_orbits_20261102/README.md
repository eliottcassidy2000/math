# A sharp depth boundary and an all-size Boolean response inverse

**PROVED elementary extension / FINITE-EXACT base certificate / CITED
inheritance. Same-agent audit with independent algorithms.** The assertion
that every positive odd integer other than7 and21 is a tournament Hamilton
count remains **OPEN**. No extension of the existing99,999,999 coverage is
claimed.

This session turns the incoming order-seven legality experiment into a sharp
boundary and an infinite source family:

- Two **strong asymmetric eight-vertex tournaments**, differing in one arc,
  have H171 and H173 but give identical legality answers to **every word of
  at most two commands**, where commands are actual path complement and all
  named single-vertex relocations. A three-command word separates them.
  The inherited complete census through seven makes eight the first order
  for this failure of depth-two legality reconstruction.
- Ordinal sums of r such atoms give **2^r nonisomorphic asymmetric sources
  on8r vertices**, all identical to the depth-two observer. Interior
  three-command probes read every individual atom. Exact reversible words
  test every conjunction, in6|S|-3 commands for nonempty S.
- Consequently any exact deterministic carrier for this ambient command
  language needs at least2^r distinct states on these sources alone. This
  is a lower bound, not a closed2^r-state model: complements leave the family.
- Arbitrary integer source weights have a **unimodular Boolean inverse**:
  conjunction responses form the Boolean zeta matrix, invertible over every
  modulus. There is no torsion cost, but arbitrary laws require2^r linear
  coordinates (one is the total mass).
- Uniform even/odd parity source laws agree on every proper subset of atom
  choices yet have mean Hamilton counts **172^r+1 and172^r-1**. These are
  expectations, not asserted Hamilton counts of individual tournaments.

See [THEORY.md](THEORY.md) for the proof, exact source rows, and the distinction
between Boolean inversion of a response law and the still-open spectrum
inverse.

The user's numerical prompts also have concrete, separately typed meanings:

- The H15 mixed and pure-blue five-vertex cores have automorphism groups3
  and5. Adding a source or sink produces the inherited480/288 trap split
  by2(6!/3) and2(6!/5). This is a genuine map from the five-vertex metagraph
  census to the six-vertex traps.
- PRIMEGAME has12 distinct numerators and7 denominator-only entries, giving
  **19=12+7**. Its entry-incidence multigraph is a forest plus the reciprocal
  two-cycle. Actual execution has the inherited six-cycle4301→…→4301;
  divisibility, simultaneous prime tokens and priority are indispensable.
- Three pure-blue plus five mixed classes are the eight self-converse
  five-vertex classes. Two further converse-pair nodes make ten merged
  nodes. Neither supplied radical list is a loopless adjacency spectrum;
  one is not even a spectrum of an integer matrix.

[CONNECTIONS.md](CONNECTIONS.md) gives the maps, losses, sources and hostile
tests, including the actual PRIMEGAME quartet height character. Numerical
equality alone is not used to transfer a theorem.

## Reproduce

From the repository root:

```text
python -B -X utf8 05-knowledge/results/primegame_path_orbits_20261102/verify.py
python -B -O -X utf8 05-knowledge/results/primegame_path_orbits_20261102/verify.py
```

The standard-library verifier checks literal command transitions, every
depth-two word on the base pair and all word cubes through r4, every subset
conjunction through r4, exact source resets, both base Hamilton counts by
subset DP and permutation enumeration, automorphisms, context examples, the
n5 census, and the PRIMEGAME identities. The separate optional
`probe_relocations.py 6 7 8 9 10` requires python-sat and reproduces discovery;
the proof does not depend on SAT UNSAT answers.

[AUDIT.md](AUDIT.md) records scopes, inheritance and reproduction. No theorem
ID was reserved; no global priority or formal-proof claim is made.
