# Sun's common-witness and height separation, tested on Chan–Pak

[Zhi-Wei Sun's preprint](https://arxiv.org/html/2611.14605v1)
gives an existential definition of the integers in the rationals
using 28 unknowns. The useful analogy is not that Diophantine
definability proves a tournament spectrum statement. It is the
organization of the proof:

* In [§6](https://arxiv.org/html/2611.14605v1#S6), the roots and
  scales are fixed before choosing the contact primes. The
  simultaneous congruences are solved with **one** comparison
  point. A separate real argument keeps the needed quantities
  bounded independently of those moduli.
* In [§7](https://arxiv.org/html/2611.14605v1#S7), height estimates
  turn the congruence constraints into exact relations. Merely
  solving the congruences would not give that conclusion.
* In [§8](https://arxiv.org/html/2611.14605v1#S8), finite choices
  are compressed without losing the exceptional branches or
  the requirement that the same witness satisfies all guards.

Here are two proved transfers of that architecture.

1. **[One common word for every compatible residue triple](JOINT_WORD.md).**
   The strong one-apex family with a transitive spine simultaneously
   realizes the ordinary count and both converse-related rooted
   responses modulo every modulus, subject only to the explicit
   parity and mod-3 conditions. The proof uses literal two-bit
   switches preserving both address residues. At powers of two
   the construction uses O(k²) vertices for precision 2^k.
   Yet some compatible triple at modulus M needs H≥M³/96
   in this class, regardless of implementation. In particular
   the locally admissible pair H=63, R=65 does not become an
   exact pair just by improving precision. This explains why
   a successful local selector needs an independent height bound.
2. **[One exact scalar reader for a full common-source profile](ONE_READER_COMPILER.md).**
   At fixed source order a cyclic reader with two transitive blocks
   encodes all ordered path-cover counts in one ordinary Hamiltonian
   count by greedy exact division. Acceptance **in that reader
   class** certifies one shared physical source for all coordinates.
   Increasing the depth of strong substitutions forces genuinely
   higher path-cover channels; repeated triangular wraps show
   the bound is sharp. A fixed short linear list cannot silently
   substitute for the nonlinear fixed-order encoding.

The testable next requirement is an **absolute-height lift** outside
this transitive-spine chart: choose the quotient and comparison
interface from the prescribed target, and prove that a simultaneous
residue match exists in an interval narrower than its modulus.
It would suffice to bound the possible absolute residual strictly
below its proven divisibility modulus for the same physical source.
The [exact factor and quotient selectors](../chanpak_family_selection_20261114/README.md)
still produce genuine witnesses in their own charts, but neither
they nor the common-witness theorem establish this height step
for every admissible odd target. Chan–Pak remains open here.

Run **node replay_joint_word.js** and **python3 replay_one_reader.py**
for the portable checks.
