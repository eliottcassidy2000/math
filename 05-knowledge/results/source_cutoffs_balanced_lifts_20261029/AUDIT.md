# Scope and verification audit

**Same-agent review, with alternate computational paths.** No independent
subagent audit is claimed. The mathematical session spans October29--30,
2026; the directory keeps its starting date.

## Mathematical gates checked

1. The Collatz operator is forward, killed on1, with original weights
   1/[j(j+1)]. The clock is an integer count of actual odd steps. The
   older inverse-base Poisson theorem has different weights and was used
   only to recover its finite-observer obstruction.
2. Finite signed LP duality retains one shared endpoint variable. Feasible
   true values make the dual nonempty; discount z<1 makes the stopped
   occupation inverse finite even in cyclic controls. The no-gain result
   assumes a supersolution upper bound, not arbitrary authentic data.
3. Reverse-clock score improvement requires positive residual at both
   endpoints before spending, and a budget no larger than the endpoint
   residual. Negative overspending switches the price to the lower bound.
   The3->5 example is an imported known suffix, not new convergence proof.
4. Ramified primality tests EVERY ideal. The last delta-level is handled
   by a valuation-one parent orbit; the simpler p-reduction test misses
   composite examples. k=1 is excluded from the prime theorem.
5. The Fourier dual action uses the symmetric multiplication matrix in
   basis(1,mu). The coefficient-of-one pairing is unimodular. The row
   energy counts antipodal pairs once and divides by orbit stabilizers.
6. The Gauss-sum proof separates even and odd product valuations. Its
   stationary phase has a unit quadratic coefficient because p>=5,
   a is a unit, and v_p(a²+4)=1. The zero character/product case is
   handled separately. Formula(4) in LOGFREE.md is checked by exact
   cyclotomic arithmetic, not fitted from magnitudes.
7. Legal-menu rounding uses at most p choices per group, an arbitrary
   fractional starting point and M>=G. Every menu leaf has the exact
   parent signed sum. Kernel choices are included, and the inherited
   modes are a lower constraint rather than a new independent draw.
8. The regular-composition kernel is signed, nonzero and one-dimensional.
   Even skeletons are impossible. Canonical strong-module decomposition
   and the actual decorated action justify the class recurrence. A
   switched candidate is re-tested for primality; switching is not called
   a native metagraph operation.
9. The explicit regular11 graph has a rigid prime7 quotient, so its full
   automorphism group is exactly the product of its two triangle groups.
   The incoming regular-parent formula gives505 order12 children, not
   the total number of order12 tournaments.
10. A nonidentity automorphism common to a source and its native path
    complement would fix their agreement Hamilton path, impossible.
    Arithmetic symmetry cannot replace the path or its changed embedding.

No demonstrated error in a maintained prior theorem was found. The
incoming correction to the naive noncoherence test is preserved as the
hostile motivating L1. Newly derived formulas were reviewed before their
PROVED status; no reserved theorem stub is used as a dependency.

## Reproduction

From the repository root, with Python3, NumPy, SciPy and SymPy installed:

```bash
python3 05-knowledge/results/source_cutoffs_balanced_lifts_20261029/check_cutoffs.py
python3 05-knowledge/results/source_cutoffs_balanced_lifts_20261029/check_lifts.py
python3 05-knowledge/results/source_cutoffs_balanced_lifts_20261029/check_composition.py
```

Repeat each with `python3 -O`; the JSON outputs are identical to
`cutoffs.json`, `lifts.json` and `composition.json`. No subprocess downloads
or external solver certificates are needed. The LP solver only proposes
flows; exact rational dual vertex enumeration certifies their scores.
The third box stores the exact upper value at one coordinate and a
coarse value at the others:30 cases fail the supersolution test on the
chosen support. A negative flow returned by the solver is not by itself
evidence of strict gain; the guarded native3->5 example verifies that
gain separately.

For the lifts, integer character polynomials certify parent transport;
cyclotomic reduction certifies row variances and individual orbit norms;
ideal cuts certify the full finite family; direct native pair closures
certify the best625-vertex candidate independently. Spectral intervals
use the inherited rational sine enclosures. Array integer bounds in these
declared universes are below signed64-bit capacity. The proof itself
does not depend on NumPy arithmetic at untested larger levels.

The composition checker tests1024 exhaustive order5 cases and seeded
samples of2000 order7 and300 order9 cases. It verifies full quotient
automorphisms by all5040 permutations only for the displayed order7
witness. It then directly enumerates the2048 binary subsets under the
nine child rotations. Canonical decomposition supplies fullness of that
action; it is not inferred by inspecting only these nine permutations.

## Literature and provenance

Classical modular decomposition is checked against Boussaïri--Chergui--
Ille--Zaidi, [Theorem16](https://arxiv.org/html/1805.04607v1). The arbitrary-
start partial-colouring input is checked against Lovett--Meka,
[Theorem4](https://arxiv.org/html/1203.5747#S2). The finite LP, Pfaffian
kernel and quadratic-sum algebra are proved locally. No priority claim
is made for those classical mechanisms or for their general uses.

Read incoming revisions: `ffb722baad` and, mid-session, `3dea39a85d`.
In particular Session110's prime-affine classes are not substituted for
full graph types, and its nilpotent extra-symmetry hostile is retained.
The primary checkout and all unrelated work remain untouched.
