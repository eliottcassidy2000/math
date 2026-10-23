# Same cubic tunneling, different fourth-order splitting

Status: PROVED finite-model comparison, using the exact sector polynomials and independently audited edge-model Schur expansion. Neither parameter limit is a continuum Yang–Mills limit.

The models have related zero-syndrome/code ground geometry but DIFFERENT configurations, gauge quotients and kinetic operations.

## Line-flip code model

The128-state Hamiltonian from [FINITE_CODE_SPECTRAL_BRIDGE.md](FINITE_CODE_SPECTRAL_BRIDGE.md) is h times the seven-cube Laplacian plus J times the nonzero-syndrome penalty. Its16 low representation branches have cubic coefficients−42,−6,+6,+42 and fourth coefficients−203,−11,−11,−203, with sector counts1,7,7,1:

    E=7h−7h²/J+c3 h³/J²+c4 h⁴/J³+O(h⁵/J⁴).

Consequently its FULL128-state gap is

    36h³/J²+192h⁴/J³+O(h⁵/J⁴).

The extra coefficient192 is derived exactly from the degree-two/three sector characteristic polynomials, not fitted to numerical eigenvalues.

## Vertex-gauge-projected edge-flip model

The [edge model](gauge/EDGE_KINETIC_FOURTH_ORDER_ADDENDUM.md) acts on the32768-dimensional space invariant under the64 vertex cuts of the21-edge field. Its potential is J(21−R)/8, where R is the sum of21 normalized rainbow Wilson loops. It has16 ground sectors at h=0 and21 single-edge kinetic generators.

Its low energies are

    E_mu=21h−21h²/J−6mu h³/J²−77h⁴/J³+O(h⁵/J⁴),

where mu=7,1,−1,−7 with multiplicities1,7,7,1. Hence its gap on this B-INVARIANT SPACE is

    36h³/J²+O(h⁵/J⁴),

with zero fourth-order coefficient. This is not a statement about the unrestricted2^21 edge Hilbert space, whose additional vertex-gauge sectors can carry different splittings.

## Mechanism and significance

Both ground-sector graphs are the crown graph K_(8,8) minus a perfect matching. Seven minimal three-flip moves, each with six orderings through equal-energy intermediate states, produce the common cubic splitting. In the actual edge model the curvature kernel contains no weight-two or weight-four word, so every four-flip return to the ground space is diagonal. Its fourth-order Schur term is therefore scalar:518 from excited intermediate returns and the energy-dependent correction441 combine to−77. Scalar corrections do not split the low levels.

The line model has different four-step logical structure; its fourth-order term is not scalar. Matching the cubic coefficient therefore did not establish equivalence of the dynamics. The fourth-order calculation is the first displayed discriminator.

This supplies a useful positive and negative result simultaneously: a finite low-order effective description can genuinely transfer through a shared logical graph, but higher-order virtual paths must be retained before claiming a stronger identification.

Reproduce the line-side addition with `python3 check_line_model_fourth.py`; the edge-side executable and independent path audit are linked in its note. The independent code-model polynomial expansion and edge-model path expansion are different validation routes.
