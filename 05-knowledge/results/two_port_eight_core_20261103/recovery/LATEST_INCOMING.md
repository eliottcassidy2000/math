# Incoming proof transfers at the close of this session

The baseline integration in [INCOMING.md](INCOMING.md) is retained unchanged
with its original receipt pins. This addendum incorporates the actual
proofs arriving through main checkpoint `7a64b40972`, including the separate
user-requested two-port research session. Status and operation scope are
kept separate from numerical acceptance.

## A stronger ambient test and a smaller exact converse battery

**PROVED (incoming):** [T1–T4, two_port_twelve_principles_20261103/PROOF.md](../../two_port_twelve_principles_20261103/PROOF.md).
For same-size nonempty modules M,N and any quotient Q,

    H(Q[M])-H(Q[N]) = H(Q)(H(M)-H(N)) modulo4.

The difference of each higher ordered k-path-cover term is divisible by
`2 k!`: its set-partition summands are differences of odd products. Ports
can be adjoined inside the enlarged uniform module. Consequently equality
of ambient full responses implies equality of the intrinsic full responses
modulo four, for any number of ports. This applies to arbitrary module
pairs, not only converse pairs. It is a necessary descent condition; no
exact cancellation of the higher-run terms is claimed.

For a core and its converse the exact integer response difference has
total degree at most three. Quartics cancel under reversal and port exchange;
both pure quadratic families also cancel. Thus the difference is determined
by `1+2n+n^3` queries: **1,753** at twelve. The general full response still
has a 6,241-entry battery there. The smaller battery does not compare
arbitrary pairs of tournaments or replace the full three-path tensor.

The incoming shared-label cubic coefficient is especially concrete. For
distinct u,v at order twelve, the coefficient of `a_u a_v c_v` in
`(F_T-F_Top)/2 modulo2` is the XOR of `S_u(T[W])` over supports W excluding
v with sizes **2,4,6,8**. The middle paired path must be the singleton v;
dropping singleton middle paths would delete this information. Independently,
`(C_uv-C_vu)/2 modulo2` counts odd-support directed paths between u,v in
either direction. This last test is equivalent to the old cover-skew
condition, not a new independent obstruction.

[converse_transfer.py](converse_transfer.py) applies these maps to our saved
full named twelve-core data without importing the incoming implementation.
Its matching [normal](converse_transfer.json) and
[optimized](converse_transfer.optimized.json) receipts verify:

- Every one of the **6,241** exact contrast coefficients and reconstruction
  of the full battery from the **1,753** allowed query values.
- Both shared-label cubic families: **264** ordered checks, retaining the
  complete twelve-by-twelve matrix with **60** nonzero entries.
- **132** directed cover-skew checks against independent path-parity DP.

All **1,584** possible exact mixed cubic coefficients are nonzero for this
particular control. This does not make the control a modular survivor: it
already fails the old cover-skew condition. The independent
[incoming proof review](TWELVE_INCOMING_REVIEW.md) also retains the precise
hostile examples. Independence of the cubic test from the conjunction of
endpoint parity and cover symmetry remains **OPEN**.

Combined with this packet's [C2–C4](../algebra/PROOF.md), the next invisible
twelve-core must be prime and pass the full intrinsic modulo-four converse
contrast. The twin branch is already excluded by the actual ambient sign.
Neither statement proves that a prime survivor exists or that all are
excluded.

## The physical missing44 return is now present, and still isolated

**FINITE-EXACT / VERIFIED:** [Session183 missing44](../../session183_two_branch_frontier_20261104/missing44/PROOF.md)
and its [independent native/port audit](../../session183_two_branch_frontier_20261104/missing44_audit/REPORT.md)
give the actual same-name eighteen-vertex continuation
`1,537,279 -> 4,611,841 = 3H+4`. Three suffix reversals and seven root-edge
reversals are compiled into twenty guarded native complements, deleting
only fresh helpers. The retained defect is **499,420**. All directed joints
and physical ports are recomputed; an equal-count replacement is not used.

**PROVED (incoming):** this specific rewrite has digit four for an arbitrary
nonempty prefix A exactly when A is the transitive three-vertex tournament.
Its coefficient equation is
`62991(P_A(2)-8H(A))=4(1-H(A))`. The cutting injection bounds the prefix
order by seven, and coprimality then forces `H(A)=1`, `P_A(2)=8`.
Filling the finite missing return does not give reusable two-digit closure.

The separately constructed [24-vertex representatives and endpoint
compiler](../../session183_two_branch_frontier_20261104/ports/PROOF.md)
realize the same count with four distinct named graphs and six chart
presentations. A six-state scan computes complete physical endpoint and
directed-joint data when deleting the root leaves a transitive order.
Its preserved sidecar is the entire order and root cut; scalar H and field
phases do not suffice. These representatives are not that eighteen-name
lineage. The [independent audit](../../session183_two_branch_frontier_20261104/ports_audit/REPORT.md)
states the exact family-completeness boundary.

## Exact boundaries for uniform arithmetic return languages

The incoming [uniform native implementation](../../recursive_seed_switching_20261103/transport/UNIFORM.md)
supplies a guarded program for the inherited `7*2^t-5 -> 21*2^t-11`
+4 family at every clock `t>=1`, with defects 6 and 12. This fills its
native-implementation obligation; the target family is still not shown
to return to a class supporting both digits.

**PROVED (incoming):** [Session183 all-mask classification](../../session183_two_branch_frontier_20261104/other/PROOF.md)
shows that the transitive-suffix chart has no uniform +2 map, and its only
uniform +4 maps have masks `7*2^a -> 21*2^b` with clock shift `a-b`.
Its [independent review](../../session183_two_branch_frontier_20261104/other_audit/REVIEW.md)
checks the binary carry proof. This replaces the formerly finite scan's
open all-mask question; it is not an obstruction to nontransitive outputs.

The [quantitative escape theorem](../../session183_two_branch_frontier_20261104/induction/PROOF.md)
and [audit](../../session183_two_branch_frontier_20261104/induction_audit/AUDIT.md)
prove that no infinite lineage with a bounded digit alphabet stays in this
chart, even with varying legal clock shifts. For odd mask u the lower bound
is `|J(3u)-3J(u)| >= 3u/8-4 bitlength(u)-2`; the odd parts along a lineage
would be `3^k u`, contradicting the resulting finite bound.

The [terminal-context law](../../session183_two_branch_frontier_20261104/mining/PROOF.md)
and its [audit](../../session183_two_branch_frontier_20261104/mining_audit/REPORT.md)
expand the positive +2 family: a transitive suffix with its final two
vertices marked admits the prescribed nontransitive output at every outer
transitive clock, regardless of earlier markings. For the particular
paired +4 recipe, arbitrary preceding marked context introduces a slope
term `-2 P_A(2) f`. Both prescribed digits are uniform exactly when that
context is unmarked and transitive. Alternative +4 recipes remain open.

## Orthogonal arithmetic arrivals and the common missing-coordinate test

The fetched [Session184 fixed-support certificate](../../session184_proof_frontiers_20261104/oddperfect/PROOF.md)
gives an exact one-row-per-prime finite selection problem. For its integer
valuation bounds `[L_r,U_r]`, the modulus
`m_r=1+max(|L_r|,|U_r|)` makes modular zero equivalent to integer zero.
Smaller moduli only give rejection filters. This is a prescribed-support
decision certificate, not exclusion of every odd-perfect support. Its
positive integer solution counts must not themselves be reduced modulo a
prime: that would erase witnesses. We retain this precise analogy to the
tournament sieve, not an asserted map between the arithmetic objects.

The same audit question explains the useful connections:

| Source to target | Preserved predicate | Lost information / needed sidecar |
|---|---|---|
| Ambient module responses to core modulo four | necessary congruence equality | exact higher-run counts remain in full endpoint data |
| Core/converse response to 1,753 queries | exact equality within converse pairs | arbitrary unrelated pairs need the general battery |
| Ordered cover skew to odd-support paths | exact vanishing modulo four | integer scale and even supports are not retained |
| Named physical source to scalar H | the count | reentry needs actual names, ports and joints |
| Private coordinates to list constraints | every required unit edge and allowed colour | simultaneous injective geometry needs rigid-contact tests |
| Bounded valuation sum to finite residues | exact integer zero with no-wrap bounds | without those bounds, residue acceptance is inconclusive |

These maps do not resolve universal reconstruction, the global Hamilton
spectrum, arbitrary odd-perfect supports, or six-chromatic forcing. They
identify concrete next tests and the exact information those tests require.
