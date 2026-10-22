# Sharply scoped next proof targets

These are OPEN research directions, not new theorem claims or tasks already running.

## 1. Does blue reflection increase every forest coefficient in expectation?

Let Q_n be the blue fixed-path forest mean, P_n the all-presentation mean from the exact Hamiltonian-size-biased formula. The exact n<=8 calculations give nonnegative coefficients of Q_n−P_n. Low degrees are proved: the quadratic difference is delta/8 and the cubic difference is delta*n*(n−3)/16, with delta=n−3 for odd n,0 otherwise. Even n6 first differs at degree4.

Candidate statement: Q_n−P_n has nonnegative coefficients for all n. This is NOT proved. Marginal arc averaging is invalid; simply invoking positive correlation or reflection positivity is not a proof. A direct route would have to compare forests with consistent versus inconsistent orientation demands in reflected edge orbits, preserving roots, the fixed path and acyclicity.

A bounded algebraic next step is to derive the fourth coefficient without a huge class census. With K=diag(d)−S and tr K=0, its centered coefficient is

c2 = (tr(K²)^2/8−tr(K⁴)/4)/16.

On the blue cube K=K0+Σ x_O B_O. Exact Rademacher second/fourth moments reduce its expectation to traces of at most two orbit matrices. For even n, each B_O is a sum of disjoint rank-one nilpotents and B_O²=0. This may give an explicit fourth-coefficient formula after a finite intersection-pattern classification. This reduction itself needs careful implementation/audit before any new formula is claimed.

## 2. Endpoint-preserving approximations to the MRB moment functional

The full measure reconstructs the kernel, but a fixed finite positive atomic quadrature fails at the divergent Mellin endpoint and gives nonvanishing/growing MRB summands. A continuous endpoint model plus a controlled signed correction is the appropriate approximation class. The current model proves absolute summability of the correction, not a faster evaluation algorithm or a modular defect law.

A useful next result would supply uniform, computable remainder bounds for approximating that signed correction while preserving the endpoint expansion. A modular candidate should then specify the actual transformation of this full functional, not merely a fitted normalization at the CM point.

## 3. Forest-refined metagraph ports

The n6 blue witness has equal score energy/triangle count but differing fourth forest coefficients, while the n5 witness has identical full in/out forest polynomials but different H. Thus forest profiles refine some ports but cannot be complete. The next missing coordinate should be path incidence or a marked-endpoint forest tensor, not another scalar spectral slogan.

Any proposed sufficient realizability criterion must account for the complement involution, reflection, orbit stabilizers, colors, endpoint pairing and the independent inventory mass constraints. Matching all present aggregate identities alone is not sufficient.
