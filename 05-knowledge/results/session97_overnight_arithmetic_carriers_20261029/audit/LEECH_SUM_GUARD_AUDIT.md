# Independent audit: Golay-to-Leech sum guard

**PASS** `leech/GOLAY_LEECH_SUM_GUARD.md` and its stated classical scope.
The primary [Allcock manuscript](https://web.ma.utexas.edu/users/allcock/research/algor.pdf),
§2, printed page4, was opened directly: its common parity, mod4 Golay
residue and sum-mod8 conditions agree exactly, with inner product scaled
by1/8. The [errata](https://web.ma.utexas.edu/users/allcock/research/algor_errata.html)
concern other table/exponent entries and do not alter this construction.

The proofs are valid:

- Addition of two odd vectors introduces the all-ones codeword, retained by
  the code. The sum condition and negation are compatible, giving closure.
- Common parity and code-residue indices are2^23 and2^12; the final sum
  parity is genuinely independent, toggled by4e_i, giving total index2^36.
- An independent modular count gives the same result: among8^24 coordinate
  residues there are2*4096*2^23 admissible vectors, hence index2^36.
- Squared integer lengths are multiples of16 for both parity sectors.
  Closure and polarization imply integral pairings after division by8;
  the scaled covolume is1, establishing even unimodularity.
- Minimum-length case analysis is exhaustive. Even vectors have either an
  octad of entries at least2 in absolute value, or multiples of4 whose sum
  guard excludes a lone nonzero entry. Odd all-unit vectors fail the sum
  guard, so one entry of absolute value3 is required.

`check_leech_independent.py` imports no producer code. It reconstructs the
4096-word Golay code and explicitly enumerates every minimum vector using
compact coordinate encodings, checking uniqueness, all guards and squared
integer length32. Counts are1104,97152,98304, totaling196560. The shell
exhaustiveness is a proof from the preceding cases, not a sampled claim.
Normal and optimized outputs agree.

The omission hostile is correct:4e_i/sqrt8 would have norm2, while the
all-ones vector/sqrt8 would have norm3; both fail precisely the sum guard.
Plain Golay Construction A instead has48 roots and is the A1^24 Niemeier
lattice. The note does not mistake it for Leech, infer a Monster action,
or claim new discovery of the classical lattice or its kissing number.
