# Harmonic phase and the exact signed residual

- [Exact Fourier steering](FOURIER.md). Pair free middle bits
  at the conductor's interaction distance. The earliest changed
  column forces a half-period sign; exact orthogonality proves
  mixing at a growing modulus even with arbitrary exterior bits.
  This independently prescribes the count modulo 2^(k+1) and
  the address modulo an odd q (including powers of three).
  The exact dynamic decoder is decoder.py.
- [Sharp one-sided correction](RESIDUAL.md). The product
  approximation is always below the actual count. Alternating
  words uniquely maximize the error (up to the parity phase);
  the exact bound is O(log N) for these readers. A diagonal
  cut in product space and an exact rounding certificate follow.
  Neither locates a suitable integer address in that cut.
- [Riesel, the literal multiset, and Golay](WHEELS.md). The two
  wheel load functions have identical Fourier magnitudes and
  distinct phases. The repeated-7 period list leaves at least
  240 holes per 1820; repeated phase support is not extra load.
  Both overlap sets cannot be Golay octads in one marking, and a
  literal 24-cycle cannot act on Golay coordinates.
- [Research atlas](ATLAS.md). The archived 2/3 questions and their
  extensions are organized into 28 typed problems or subproblems,
  with 11 lines of harmonic work and the distinctions needed
  to use them without treating analogy as a theorem.

## What would turn these pieces into a proof

There are three separate gates: a target must have an admissible
address, its phase must be correct, and the physical count must
fit the target absolutely. For a Riesel cover the last step is
exact division by a certified prime, not a statistical estimate.
For the tournament construction the last step is still open.
The new Fourier bound needs roughly five free bits per dyadic
residue bit; a frozen-exterior repair would need fewer than two
to improve the normalized square-root gap. That makes the next
question concrete: change the reader or retain compensating
prefix pairs across the product curve rather than buying more
phases in one fixed prefix bin. Chan–Pak remains open.

Replay: Python 3 check_fourier.py, decoder.py, and check_residual.py.
