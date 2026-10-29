# Independent audit: all-dimension affine port reconstruction

**PASS** `supersingular/ALL_M_AFFINE_PORT_RECONSTRUCTION.md`.

For a named subset with affine span dimension d, translating an element to
zero and quotienting that span gives exactly [m-d choose2]_2 containing
codimension-two flats. This also proves the zero convention when the span is
too large. Over F2 every three distinct points span a plane, so responses
through order3 are constant. For m>=4, four points distinguish plane from
three-space because the two Gaussian counts differ. With a named origin,
the unique completing plane point reconstructs translated XOR. An affine
chart conversely answers all subsequent incidence queries, proving the
specified faithful-carrier theorem rather than merely a finite example.

The stabilizer of the plane relation is exactly AGL(m,2): after translating
a preserving permutation to fix the origin, reconstructed addition forces
linearity. Thus the labelled-geometry count and fixed-length information
lower bound are justified. The query upper bound N(N-m-1) is a valid simple
construction, not an optimality assertion. The m3 boundary is indispensable:
its codimension-two flats are all pairs, whose incidence stabilizer is S8
and whose entire oracle cannot distinguish thirty source geometries.

The RM(2,m) minimum-word proof via quadratic polar rank is correct. Radical
linear terms give balanced words; otherwise the character-sum magnitude is
2^(m-r), giving the displayed weights. The smallest nonzero weight occurs
at rank2 with the smaller sign and exactly the two-affine-form products.
No self-duality is extrapolated beyond the exceptional m5 case among m>=4.
For larger alphabets the stated result is only a triple-collinearity detector;
semilinear ambiguity is properly retained.

`check_all_m_ports_independent.py` imports no producer code. It constructs
flats by enumerating linear spans and then translating them, unlike the
producer's orthogonal-equation construction. It verifies all named subsets
through size4 for m3,m4,m5, addition recovery at every origin for m4,m5,
and all minimum supports in the full RM(2,4)/RM(2,5) word banks. It confirms
labelled-geometry counts30,64864800,822336494953469303808000000 and the
corresponding fixed-length source counts5,26,90 bits. The first is not an
incidence-observer lower bound at m3 because that observer erases geometry.
Normal and optimized Python outputs agree byte-for-byte.
