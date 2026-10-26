# Two fifteen-element carriers have different exact tournament defects

Session65. **PROVED elementary corollary** of the independently audited
involution-defect polynomial. The two carriers and their distinct fixed-point
characters are recovered from Session52's
`session52_saturation_residues_metagraph_20261024/fifteen-kernels/README.md`,
sections2–3. This is not an identification with the fifteen Ogg primes.

Let Omega be six labelled elements. Its fifteen two-element subsets are
duads, the nonzero vectors of the usual binary symplectic four-space model.
Its fifteen perfect matchings are synthemes, corresponding to Lagrangian
planes and hence the cited Richelot kernel choices. They are different
S6-sets despite having the same cardinality.

Consider the transposition(01). On duads it fixes the pair01 and the six
pairs among the other four labels: f=7, b=4 swapped pairs. On synthemes a
fixed matching must contain01; the other four elements can be matched in
three ways. Hence f=3, b=6 there. A matching without01 pairs0 and1 with
different other labels, and swapping0/1 changes those pairs, proving there
are no other fixed matchings.

The involution-defect theorem therefore gives the exact two polynomials:

    duads:     2^65 z^4 (1+z^2)^40,
    synthemes: 2^57 z^6 (1+z^2)^48.

Both enumerate all2^105 tournament orientations. Their minimum defects are
4 and6, respectively; literal minimizers are supplied in the receipt. Thus
the old fixed-point-character obstruction acquires a quantitative orientation
meaning. No ordinary tournament on either carrier can preserve this entire
S6 action, and no identity-S6-equivariant bijection identifies the carriers.

The classical outer automorphism exchanges the relevant symmetry types.
Our checker confirms that a triple transposition on duads has f=3,b=6,
matching the transposition character on synthemes. This one check does not
prove the whole outer-action equivalence; the explicit inherited Session52
construction is the source for that stronger fact.

This is precisely why retaining the action matters more than counting
fifteen objects. A bare bijection between two finite lists does not preserve
Frobenius, symplectic incidence, tournament defects or computational meaning.

## Verification

    python -B check.py > normal.json
    python -O -B check.py > optimized.json
    cmp normal.json optimized.json

The independent small checker constructs both complete carriers, their actual
permutations, all unordered-pair orbits covering the105 pairs, and explicit fifteen-vertex
orientations attaining the respective minima. It checks the complete
polynomial totals, not by enumerating2^105 tournaments but by the proved
orbit factorization. No use of an arbitrary cardinality-matching bijection
or a claimed new arithmetic Richelot theorem is made.
