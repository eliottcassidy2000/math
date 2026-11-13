# Independent review: sixfold positive return

Status: PASS for sections 1–7, after a reproducibility correction.
Reviewed the mathematical proof independently; the source audit does not
assert novelty over classical positive-matrix/cyclic-representation theory.

## Proof inspection

The map `g(x)=-3/(x+3)` is on the projective rational line, including infinity.
The two complex fixed points remain the only fixed points for powers 1–5,
because the eigenvalue ratio has exact order six. The quadratic map in the
historical quarter-integer example is explicitly a different map. The old
preperiodic classification is inherited, not re-proved by this audit.

The symmetric-square calculation and positive functional are correct.
In particular the representing form has determinant `3/4>0`; the observer
coordinate `v` can be negative, whereas the six ray coefficients are the
nonnegative state variables. This distinction is essential.

For the lower bound, the surjection makes the sixth primitive eigenvalue a
quotient eigenvalue, hence an eigenvalue of M. Conjugating by the diagonal
of its positive conserved mass turns M into a column-stochastic matrix.
Transpose and take a maximal-modulus eigenvector component. Every outgoing
positive edge from such a component has equal modulus and advances its
phase by the eigenvalue. The maximal-modulus support is nonempty and closed,
so it contains a simple directed cycle. Its length is divisible by six,
forcing at least six vertices. No irreducibility assumption is missing.

If the total dimension is exactly six, the simple cycle exhausts all states
and gives distinct phases to them. The row-stochastic transition is forced
to the unique next phase, so the original M is a weighted six-cycle. Positive
mass forces the product of weights to one; integer positive weights force
each to one. There is no unexamined reducible extra component at dimension
six. Dropping conserved positive mass or the primitive observer mode invalidates
this argument, as the text correctly warns.

The five-state 2+3 permutation example has group order six but no primitive
sixth eigenvalue. It is a valid hostile to conflating action order and observer
spectrum. No physical tournament transfer is inferred from this construction.

## Independent controls

`verify_six_independent.py` imports no primary lane code and uses explicit
raises, hence remains active under `python -O`. It derives the symmetric-square
matrix from its definition, checks all six rays, positive mass and exact
period, checks 120 rational starting points (counting repeated representations),
and enumerates permutations at dimensions 1–6.

Both normal and optimized runs give 154 checks and the same output:

- dimension5:20 permutations of order6, zero with a primitive sixth mode;
- dimension6:240 permutations of order6,120 with a primitive sixth mode.

The primary `verify.py` initially used bare asserts. The author was asked to
replace them before claiming an optimized replay. This was a test-harness
issue, not a failed mathematical identity.

## General-q extension and branching

Sections6–7 also pass. The condition q>=2 is present and essential: it keeps
the trivial mass representation separate from the primitive cyclotomic one
and makes Phi_q(1) nonzero. The constructive signed image proof is valid:
adding the stated constant multiple of the monic cyclotomic polynomial fits
within q source coordinates. The kernel proof uses exact divisibility, not
only a rational-vector-space dimension count. Thus the stated congruence and
integral index are correct. Nonnegative image membership is explicitly not
inferred from this lattice statement.

The q=6 balanced token identities hold in the actual displayed six rays.
Their applicability requires inventory on the subtractive side, correctly
retained in the proof. Completeness of those moves is not claimed.

The irreducible stochastic cyclic-class iff is sound: irreducibility
propagates equal maximal modulus to every state, and an eigenphase gives
exactly q classes. In the reducible case only an appropriate closed class
is asserted, not a global phase partition. This is credited as classical
peripheral-spectrum structure. The transition normalization multiplier is
separately retained, so the result does not confuse probability conservation
with raw Hamiltonian-count conservation.

The independent verifier now generates cyclotomic polynomials recursively
for q=2,...,30 and checks source-column congruences, constructive signed
image preimages and all degree-bounded kernel generators. It executes 2,199
explicit checks in total. Normal and optimized output are byte-identical.
The primary verifier's assertions have been replaced by explicit require
checks; the initial reproducibility issue is resolved.
