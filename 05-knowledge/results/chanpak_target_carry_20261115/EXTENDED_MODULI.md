# A second bank of arbitrary-sum certificates

The [arbitrary-high-sum certificate](ABSOLUTE_CARRY_DECODER.md)
has no theorem limiting it to the first six moduli. The independent
[C++ certificate](certify_extended_moduli.cpp) checks these additional
pairs, with state masks implemented as unsigned 128-bit integers:

| Modulus M | Word length L |
|---:|---:|
| 7 | 12 |
| 9 | 12 |
| 10 | 11 |
| 12 | 12 |
| 16 | 10 |

For each required row it exhausts **all** 2^L prescribed sums
and both high residues modulo 2M, on the larger of the two
carry branches. Both adjacent-digit channels are checked. It
iterates until the two masks contain every (address mod oddpart(M),
half-checksum mod M) pair. The actual high prefix only adds a fixed
translation to the half-checksum and may be arbitrarily large.
There is no assertion that the listed lengths are minimal.

The proof of the physical bracket therefore applies unchanged at
these moduli: O_M(N^(3/8)) error, counts congruent to N mod 2M, for
every sufficiently large odd N. In particular moduli divisible by
9, 5, or powers of 2 are not limited to all-ones high sums. This
does not prove a uniform result as M grows.

Replay with a C++17 compiler supporting `__uint128_t` (GCC or Clang):
`g++ -std=c++17 -O3 certify_extended_moduli.cpp -o /tmp/cert_extended && /tmp/cert_extended`.
Compare [the receipt](extended_receipt.txt).
