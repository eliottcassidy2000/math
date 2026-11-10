# Sequence generators, recurrences, and the Chan–Pak bridge

Session249, 10 November 2026. **COMPLETE.** New statements are proved
in the lane notes with small exact implementation checks; no independent
canon audit or literature-priority claim. Chan–Pak remains OPEN.

Read the [final synthesis](SYNTHESIS.md) first. Run `python3 run_all.py`
from this directory to replay the four compilers. All use the standard library.

The goal is to derive equations from operations, not to fit longer prefixes.
The four lanes deliberately use different sequence indices:

1. **[Burnside](burnside/README.md):** ordinary tournament classes by vertex order (A000568),
   odd cycle partitions, connected cumulants, divisor and lcm transforms.
2. **[Grammar](grammar/PROOF.md):** Hamiltonian-path counts in repeatedly inflated tournaments;
   rationality versus factorial growth from the growing quotient.
3. **[Heat](heat/README.md):** path-count covariance by vertex order, with edge-noise and
   vertex-support gradings retained together.
4. **[Rays](rays/PROOF.md):** certified path-count values by numerical value;
   positive Mahler equations, binary languages, and the exceptional point11.

## Type contract and inheritance

A000568(n) counts ordinary tournament isomorphism classes on n vertices,
including the empty class at n=0. It is not the set of Hamiltonian-path
counts, the maximum number of paths, or the number of converse pairs.
H(T) is the number of ordered directed Hamiltonian paths. Chan–Pak asks
whether every positive odd integer except7 and21 occurs as H(T).

Closest mechanisms: THM-5097/5098's recent Burnside engines; THM-3121/3181's
path-cover substitution; THM-4761's even-path forest Fourier expansion;
S248's birth/resampling oscillator; THM-5126's certified dyadic lines.
All are inherited rather than re-audited. Their theorem files and exact
paths are cited by the lane notes. A000568's formal series has radius zero;
fixed H fibres are already rational in vertex order (Session24/239).
Neither fact should be rediscovered by indiscriminate recurrence guessing.

Hostile controls: C3 followed by C3 multiplies H to9, not6; three cyclically
interacting growing blocks have factorial growth; the all-odd target with
holes7,21 is not closed under N -> 2N-1 because11 maps to21. Near miss:
finite-order moments and exact residue distributions do not prove that
individual values occur. The underused sidecars are cycle cumulants,
ordered boundary runs, Fourier vertex support, and dyadic ray sources.

The live concept board is those four carriers plus coefficient positivity.
A proposed bridge must retain the actual target value and state what has
been lost. The [final synthesis](SYNTHESIS.md) lists the new equations, runnable
compilers, strongest examples, and the remaining coverage obligation.
