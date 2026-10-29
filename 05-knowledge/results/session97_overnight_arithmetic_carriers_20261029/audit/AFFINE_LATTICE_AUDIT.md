# Independent audit: affine reconstruction and equal-theta lattices

**PASS** `supersingular/AFFINE_RECONSTRUCTION_AND_ISOSPECTRAL_LATTICES.md`.
The new independent checker reconstructs all six zero-free code additions
from response-seven four-subsets, verifies all ordered associativity triples,
checks every octad is an affine three-flat, and recovers the same six Singer
minimal polynomials. No producer module is imported.

Construction A is normalized as L={x/sqrt(2):x mod2 in C}; therefore the
formal variable z indeed records z^(2||v||²)=z^(sum x_i²). The integral,
even and unimodular assertions follow from self-duality/double-evenness and
the stated scaling. With minimum code weight8, norm-two vectors are exactly
the coordinate roots ±sqrt(2)e_i. Their axes are intrinsic, so any lattice
isometry is a signed coordinate permutation and induces a code equivalence.
The two inequivalent code classes therefore yield nonisometric lattices.
Equal entire theta series follows symbolically from the common full weight
enumerator, not from a truncated numerical check.

At lattice norm4, the unscaled sum of squares is8. The permitted types are
two ±2 coordinates (zero parity) or eight ±1 coordinates (an octad); a
±2 coordinate with four odd coordinates would require a forbidden codeword
of weight4. Consequently every odd root-port response counts exactly256
sign choices per incident octad, validating256*R4 and its lower-port analogues.

The independent series expansion recovers norm0,2,4,6,8 coefficients
1,64,160704,64543488,4845725632. Normal and `python -O` outputs agree.
The Golay construction retains48 norm-two roots and hence is not Leech.
The Hamming case has extra roots and is correctly outside the root-axis
criterion's minimum-weight-eight hypothesis.

An exact further classical identification is available without numerical
matching: the six-code affine class is RM(2,5). Indicators of affine
codimension-two flats are products of two affine linear Boolean functions.
Conversely their span includes each x_i*x_j, each x_i by adding two parallel
flat indicators, and the constant by summing four parallel flat indicators.
Thus they span all sixteen degree-at-most-two Boolean monomials. This is an
identification of this code class, not an inferred tournament operation.
