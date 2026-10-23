# Independent audit: equivariant residual-certificate reuse

23 October 2026. **PASS, scoped.** Read [RESIDUAL.md](../RESIDUAL.md) and
[check_residual.py](../check_residual.py), checked the group action and
baseline quantifiers, and replayed the revised checker normally and with
Python optimization. Both runs pass19,721 checks. No native tournament-port
census is claimed by this toy transport experiment.

## The proof obligations are sufficient

For every transported branch, the retained maps preserve row/column
incidence, channel identities, bounds, and margins. Therefore they map an
entire feasible integer-table set bijectively to the other. Preserving
each channel's vector contribution in the same observation coordinates
then preserves both every table's baseline and every signed residual arc's
cost. Positive residual capacities, SCCs, and the potential equations
transport. These statements prove the claimed reuse for all branches in
an authenticated orbit; they do not depend on the finite toy enumeration.

The claim is conditional on supplying feasible tables and valid potential
certificates on representative nonempty branches. Empty branches cannot
be certified by simply omitting their feasibility obligation. The report's
explicit feasible-table premise and authenticated-native-table condition
retain this boundary.

Different orbits still need a common baseline. If observation coordinates
are permuted, covariance gives a permuted baseline; invariance requires
the additional fixed-vector condition stated in the report. A merely equal
number of channels, an anonymous signature collision, or an automorphism
of the small cells without a compatible action on the port construction
does not meet these hypotheses. Source marks, path order, converse data,
and colour multiplicities remain part of that construction.

The torsor reduction is also conditional in the correct direction: actual
comparison isomorphisms act transitively on their own set, but reusing one
transport certificate requires covariance of the actual branch models
under that action. The report does not infer the latter automatically.

## Independent check of the finite census

For K_(3,3) with unit margins, all feasible tables are the six permutation
matrices. All have the same scalar cost exactly when every2-by-2 alternating
cost sum vanishes: every such switch extends to two permutation matrices.
Equivalently C_ij=r_i+c_j. Normalize r_0=0.

For entries in {-1,0,1}, let d=max(r)-min(r). The ordered pair(r_1,r_2)
has range d=0,1,2 in respectively1,6,12 ways, and no larger range is
admissible. Each of the three column offsets then has3-d choices. Thus
the number of constant cost matrices is independently

    1*3^3 + 6*2^3 + 12*1^3 = 87.

The remaining3^9-87=19,596 cost matrices are nonconstant. This matches
the checker's complete census without reusing its residual-SCC procedure.
The transported positive control preserves labelled costs by construction;
it is an abstract branch-model test, not evidence that a particular native
port construction has the necessary covariance.

## Hostiles and the small repair made during audit

The paired odd-phase hostile has exactly three admissible shared branches.
Unioning them gives two independent3-by-3 permutation tables, hence36
tables, and off-diagonal-count differences ranging from-3 to3. It exposes
the loss of a common chart choice. It is not itself used as an example of
channelwise observation-preserving symmetry: shifting a matching can
change the off-diagonal cost of its individual channels.

The first draft's second hostile merely stored the numbers0 and3. The
parent replaced it after review with two literal permitted matching
graphs, j=i and j=i+1, each with unit margins and its unique feasible
table. Their residual tests now actually run and both pass; the SAME
off-diagonal-count observation gives baselines0 and3. This directly
verifies that separate within-branch constancy does not establish common
coverage. The revision is included in the normal and optimized replays
reported above.

No further correction is required for the theorem or its stated scope.
