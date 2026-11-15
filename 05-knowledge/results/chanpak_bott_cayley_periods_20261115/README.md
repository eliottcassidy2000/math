# When a pattern survives composition: Bott, Cayley–Dickson, and the Chan–Pak clock

**Main result.** A nonlinear three-packet decoder replaces the old
quadratic-length dyadic guard with a **linear-length physical guard**.
Every sufficiently large odd N has a literal three-count tournament
window on each side, with both exact +2 steps, within O(N^(3/8)). The
same designated count in each window is N modulo 3·2^k with
2^k ≫ N^(1/160). The prior result at that error budget had only
2^(c sqrt(log N)) precision. [Decoder](THREE_PACKETS.md) ·
[Placement theorem](POLYNOMIAL_PLACEMENT.md). Chan–Pak remains open:
this reduces an integer quotient, not to zero.

## Four exact tests for a shared pattern

1. **Is it compatible with every allowed move?** The complementary
   address words factor exactly through three invariants. A hidden signed
   coordinate is invisible to *every* such continuation. Two actual
   terminal probes recover it (also at 7 and 13 where the first probe
   alone is blind). A repeating word gives a contracting direction and
   a translating direction, not a full return. [Reflection and return
   theorem](REFLECTION.md).
2. **Is the hidden defect observable at the next layer?** The old
   Cayley–Dickson doubling continues after the octonions, but the
   associator enters the next product's norm. The tournament cyclic
   join has an exact parallel: the scalar norm is rank two, the
   optional-slot reader is rank three, and two physical sources can
   agree at the first wrap and disagree at the second. The odd-grade
   mod-8 residue still survives because it follows from each visible
   coordinate's parity, not from the rank. In every odd slot grade the
   full 3 mod 8 sheet survives; grades 2 mod 4 give 2 mod 8 and grades
   0 mod 4 give exactly 2 mod 16. [Tower theorem](CAYLEY.md) ·
   [Complete dyadic grades](DYADIC_GRADES.md).
3. **Can the equivalence be carried with a marked inverse?** Bott's
   period-eight theorem is a theorem in its stable category, not a
   claim of eight-step return for integer orbits. The marked Clifford
   matrix corner in the old repo note preserves the Collatz affine
   numerator that similarity can forget. On the tournament side,
   literal probes characterize which path-cover grades a reader sees.
   At equal orders, preserving the first two grades forces changes by
   unmarked cyclic readers to be multiples of 960 (odd orders) or
   2880 (even orders). An order-nine example makes the odd-order
   power of two sharp. [Factorial obstruction](BOTT.md).
4. **Does the lifted state remain physical?** Adding the missing
   imaginary coordinate repairs the algebraic norm identity, but
   products of two physical ordered-path profiles have even
   imaginary coordinate, which no physical input profile does.
   Conversely, two off-diagonal readers give an exact cover of all
   dyadic residues with only O(log² X) actual counts up to X.
   [Norm closure versus placement](NORM_CLOSURE.md).

The construction and obstruction use the **same** mixed carry.
Independent triangular controls force exact valuations equal to
their spatial separation, hence quadratic total width. Correlated
dense packets need no such independence: their mixed interaction
is a perfect finite Fourier pairing, even with arbitrary internal
phases. Three packets give positive support at every target residue.
That is why an apparent obstruction became an ingredient once the
allowed composition changed. [Separation theorem](SEPARATION.md).

The next absolute step is to let a source or reader respond to the
remaining signed quotient without losing marked continuation data.
The factorial obstruction rules out hiding the repair solely in
higher unmarked path-cover grades.

## Replay

Run `python three_packets_check.py`, `python verify_placement.py`,
`python verify_reflection.py`, `python verify_cayley.py`, and
`python bott_check.py` in this directory. The constructions and
obstructions are proved above; these scripts are finite controls,
not evidence substituted for an all-order argument.
