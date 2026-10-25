# Parent independent mathematical controls

Status: FINITE-EXACT computations plus explicit proof review. These controls supplement, not replace, packet proofs.

## Ordered arithmetic candidates

`check_ordered_counts.py` constructs projective points by first-nonzero-coordinate normalization, enumerates unordered distinct triples, rejects any coordinate present in exactly one point, and multiplies by six. It independently agrees with the arithmetic packet C_k(q) formula for q=2,3 and k=1..4, and q=5 and k=1..3. All tests use exceptions, not removable assertions.

The transposition Frobenius correction was also checked by a separate argument: a twisted triple consists of one rational point P and a pair Q,F(Q) with Q defined over F_(q²) but not F_q. Their common exact support S has size s; P must be supported in S and is otherwise arbitrary. Thus the product is ((q²−1)^(s−1)−(q−1)^(s−1))·#P^(s−1)(F_q), summed over S. A three-cycle has all three supports equal and automatically satisfies the repeated-support condition, giving #P^(k−1)(F_(q³))−#P^(k−1)(F_q). These are different from three individually rational points.

## Signed flag polynomials

`check_signed_flags.py` uses a last-vertex subset recurrence. Its full coefficient vectors agree with extremal-insertion factors for (r,j)=(3,2),(3,4),(3,6),(2,2),(4,3), including n=14. It checks all jets through j using binomial coefficients, independently of factorial-scaled derivatives. The two lexicographic six-vertex masks255,2015 yield second binomial jets180,204, hence second derivatives360,408. Output: `signed_flags.out`.

Proof review checked: extremal insertion needs an extremal label; signed automorphism convention is orientation inversion, not edge-permutation sign; uniform reference averaging is a finite double-sum identity and does not import same-latent higher moments. The effective torus/affine model supplies a Grothendieck-class equality, not a locus isomorphism. Unbounded jet order does not prove individual profile collisions.

## Geometry and deformation proof review

Reviewed the explicit Laurent semigroup normal form of the shared-root toric variety, odd-field split/nonsplit counts, characteristic-two exception, four-node Galois obstruction, interior branch point and limited blueprint pushout claim. Separate specialist audit supplies independent exact controls.

Independently computed the Smith normal form of the seven integral magic-line differences against the top row: rank6 and six unit invariant factors, so the nine-entry additive group scheme has no hidden torsion factor. The direct group derivation of common total3E avoids division by2 or3. The deformation uses polynomial divided differences BEFORE specialization; undivided product equations would be vacuous at t=0. Ordinary square roots remain different from group doubling.

No broad F1 descent, arithmetic solution, or RH conclusion is certified by these checks.
