# Norm-seven arithmetic, tournament taxonomy, and finite certificates

2026-10-26. **PROVED deductions with cited inputs; independent audits and
FINITE-EXACT controls in the linked scopes.** No new theorem ID is reserved.

The arithmetic advance is a complete local test for the exceptional conic in the
[inherited degree2/4 fibres](../normseven_descent_moments_20261026/README.md)
over every quadratic node field. The tournament advance was independently
recovered during the session and has now been integrated with stronger concurrent
work. The four proposed TT7 motives remain unidentified.

| Result | Exact scope |
|---|---|
| [Quadratic fibre atlas](arithmetic/QUADRATIC_FIBRE_ATLAS.md) | For every nonzero rational D,δ, the specified resolved fibre has minimum point-field degree 1, 2, or 4. Over Q(√D), the conic splits iff no ramified place of (−1,δ) splits in that field. |
| [Prime-power transport](motives/EXTENSION_TWIST.md) | The six literal cone sums extend over every odd finite field with ε=χ(−1); C04=ε C02 and C06=C05. This fixes the arithmetic carrier before any Euler reconstruction. |
| [Ramification support](arithmetic/CONDUCTOR_SUPPORT.md) | The compact-support cohomology of the specified coordinate-difference double covers is unramified at odd p≠ℓ. Actual constituents inherit this support bound; their existence, purity and exponent at2 still require proof. |
| [Complete tournament code and necklace family](taxonomy/TAXONOMY.md) | Decorated strong-module codes classify every finite tournament. For opposite nonisomorphic blocks on an odd carousel, isomorphism is binary-necklace equivalence and period determines the full automorphism group. |
| [Finite certificates](certificates/ZERO_AND_FINITE_CERTIFICATES.md) | A thirteen-prime integer determinant certifies four independent raw TT7 directions modulo the nine specified functions. Total Euler-stream equality has a Π₁ form; total-zero program behavior has a Π₂ form. |
| [Galois comparison obstructions](certificates/GALOIS_COMPARISON_GATES.md) | Actual rank-eight representations agree on a noncubic bank of fourteen primes but differ at43. Correct GL₂ comparison hypotheses cannot be discarded. Euler factors also forget nonsplit extensions and integral lattices. |

The twelve-vertex 81/27 pair independently confirmed a false all-order clause
of [THM-4803, cyclic-triple fibres](../../../01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md).
While we worked, commits `d761e4cf0f` and `c9a8823436` repaired it and proved
more: the sharp threshold4p, uniqueness of the bad order12 fibre, full
path-cover invariance, weighted necklace spectra and pointed growth. Their
repairs are retained verbatim. See [INTEGRATION](INTEGRATION.md); our independent
checks are confirmation, not another claimed discovery of those results.
The final pull also supplies a literal number-theory bridge: choosing a normal
basis identifies binary-word rotation with Frobenius, and word period with
field degree. [INTEGRATION](INTEGRATION.md) records its source and scope.

The support theorem concerns freely varying labelled configuration coordinates.
It does not transfer to the norm-seven hourglass arrangement without a new
integral-model argument: that arrangement has extra triple collisions at3.

[The concept board](CONCEPT_BOARD.md) records the six/seven probes, exact maps,
lost information and stopping conditions. [AUDIT](AUDIT.md) records proof review
and finite universes. Run `python3 -B verify.py --replay` from this directory
to check hashes and replay all frozen Python receipts in normal and optimized
modes. `--pari` additionally rebuilds the modular columns with PARI/GP.

**OPEN:** the pure subquotients, ranks and conductor exponents of the four TT7
directions; interior rational/integer-square hourglasses; whole-surface F1
models; an inexpensive unrestricted prime-tournament registry; and all-order
black-core connectivity. The work establishes no reduction between RH,
Collatz, P versus NP, LRC(14), or the planar Jacobian conjecture.
