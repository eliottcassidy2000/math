# Recovered work, correction precedence, and incoming integration

The user supplied `Pasted text.txt`, beginning with the unit-ball formula
and a proposed MRB/projective/Baby-Monster connection. That same mathematical
proposal was already audited in the November2 source below. Recovery started
from clean main `9db342375914`, after integrating102 incoming commits since
the previous response packet. This ledger uses current proofs and explicit
corrections rather than treating the historical reflections as current truth.

## Authoritative bijection and operation routes

|Source|What is carried forward|
|---|---|
|[THM-4531, twisted Mallows–Sloane principle](../../../01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md)|Affine tournament action, invariant reversal quotients, and the correctly scoped naturality question|
|[THM-4721, twisted Burnside and Walsh basis](../../../01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md)|Relative-even graph orbit sums, integral Walsh basis and weighted orthogonality; its explicit October20 corrections outrank older prose|
|[THM-4722, Sylow Hall obstruction](../../../01-canon/theorems/THM-4722-no-equivariant-bijection-tournaments-even-graphs-for-every-n-sylow-hall-violator.md) and [audited uniform obstruction](../uniform_equivariance_obstruction_20261020/README.md)|No S_n-equivariant orbit-bijective tournament-to-even-graph map from n=5; not a ban on all type bijections|
|[October11 SCC/Lyndon inverse](../tournament_even_bijection_20261011c.md)|Specified all-order type codec; canonical ranking is part of the map|
|[Converse and conditional fusion](../tournament_converse_context_20261011d.md)|Graph-side converse is transported, not ordinary complementation; retained leaf-pair frame is necessary for reversible star/path completion|
|[October20 session status](../tournament_even_bijections_session_20261020/RESULT_STATUS.md), [worked maps](../tournament_even_bijections_session_20261020/WORKED_MAPS.md)|Several incompatible but valid maps, explicit failures, and remaining local-map targets|
|[Prime-local matching](../prime_local_walsh_bijection_20261020/README.md), [Sylow action](../prime_local_walsh_bijection_20261020/SYLOW_ACTION_ADDENDUM.md)|One chosen prime's embedded action is preserved; no all-prime coupling|
|[Converse/parity refinement](../prime_local_walsh_bijection_20261020/CONVERSE_PARITY_ADDENDUM.md), [positive coupling](../prime_local_walsh_bijection_20261020/SUPPORT_COUPLING_ADDENDUM.md)|Signed trace, compatible block matchings, and a canonical doubly stochastic coupling|
|[Circulant derivative/integration inverse](../circulant_necklace_bijection_20261020/CIRCULANT_NECKLACE_BIJECTION.md), [completed review](../circulant_necklace_bijection_20261020/ROOT_REVIEW.md)|Atlas-free local map on every odd-prime circulant sector; labelled forward map is a surjection|
|[Primitive Hopf/Lie matching](../primitive_hopf_lie_bijection_20261020/PRIMITIVE_HOPF_LIE_BIJECTION.md)|Connected graphs correspond linearly to primitive combinations; set matching still makes choices|
|[C15 counterexample](../c15_symmetry_counterexample_20261020/C15_ORDER_TEN_COUNTEREXAMPLE.md)|32 versus31 at order ten; hidden odd symmetry appears after forgetting the original orbit partition|
|[Sign-aware augmentation](../sign_aware_augmentation_20261020/signed_root_augmentation.md)|Retain current stabilizer and transporter sign, including odd parent states|
|[Session72 operation envelope](../session72_even_graph_unifying_carriers_20261026/session72-even/EVEN_GRAPH_OPERATION_ENVELOPE.md)|Source-independent reversal kernels diagonalize; quotient averaging removes dual modes|
|[Guarded Fourier response](../clock_even_guarded_composition_20261026/FOURIER.md)|Named ports use relative-evenness; legal Hamilton-word guards mix Fourier supports and require their markings|

The historical `even-graph-equinumerosity-one-cube-four-faces` and
`a-bijection-of-sets-is-not-a-bijection-of-structure` reflections helped locate
provenance. Their Eulerian cycle-space discussion is a different object from
the sign-even graph family above. We do not inherit their speculative
Paley/Ising/LRC identifications or merge the two meanings of even.

The relevant correction lineage is the October20 tournament/even-graph
symmetry entry in [MISTAKES.md](../../../01-canon/MISTAKES.md): it records
the C15 failure, the five-vertex joint-prime obstruction, labelled-surjection
versus type-bijection distinction, and the earlier SCC/Lyndon precedence.
No new error in an accepted theorem was found here; these existing repairs
are carried into this merger rather than duplicated as new discoveries.

## The unit-ball and five-vertex ideas, with actual surviving maps

**INHERITED PROVED:** [Session127 VOL-1](../session127_gsma_norm_collapse_20261031/volume/THEORY.md)
gives V_n=pi^(n/2)/Gamma(n/2+1) and V_(n+2)=2*pi*V_n/(n+2).
Its integer-dimensional maximum is5. Its legitimate arithmetic application
specifies the Euclidean metric, measure and lattice; a rank-two Thue slice,
a degree-five polynomial and a five-dimensional full lattice are different.

**INHERITED REFUTATION / repaired statement:**
[November2 inspiration audit](../two_selector_boolean_20261102/INSPIRATION.md)
and [parallel arithmetic/source audit](../two_edge_transport_20261102/ARITHMETIC_AND_SOURCE_AUDIT.md)
already separate that geometric maximum from tournament self-converse decay.
The latter proportion is1 at order3,1/2 at4 and2/3 at5, so5 is not its
global maximum. A fixed reversing involution leaves floor(n^2/4) free edge
bits, and the elementary upper/lower class bounds in the former audit give

    log(S_n/T_n)=-(log2)*n^2/4+O(n log n).

The pasted fixed-constant MRB expression instead has a leading
-n^2*log(n)/4 term and is incompatible with that bound. No common decay law
follows from the shared small dimension. The pasted anti-automorphism must
also not be called an ordinary automorphism.

The same audit gives the corrected adjacency polynomials
t^5-5t^2-5t-2 and t^5-4t^2-3t-3. Its projective check finds P2,P4,P5
collinear, and every proposed anchor leaves a zero determinant. An external
anchor in PG(2,3) has only four lines through it, so five selected points
cannot all give nonzero pair determinants. Representative signs also need
a gauge. These are failures of the proposed map, not failures of finite
projective geometry or of the valid Conway moving-hole construction.

There is a useful POSITIVE repair, now stronger than the original audit's
snapshot. [Session143 checkpoint3](../session143_magic_square_cyclic_bridges_20261103/CHECKPOINT3.md)
uses the regular tournament A=P+P^2, P^5=I. It proves Q[A]=Q[P] but the
additive lattice index [Z[P]:Z[A]]=11, and repairs integral recovery using
full adjacency operations, including P=A-A circ A^2. Its earlier C5 operator
bridge uses rho(g)+rho(g)^2 in a specified representation; its separately
audited final corollary constructs an ABSTRACT F20 embedding into the cited
Baby Monster normalizer selecting5A. It does not supply explicit Monster
matrices, generator words, or a representation intertwiner. We inherit this
scope instead of repeating the stale claim that no embedding exists.

The relevant connection to our response work is exact: restricted polynomial
or spectral observations can lose data that named adjacency or a marked
operation restores. It is a connection between specified observation maps,
not an identification of the unit ball, Baby Monster and twelve tournament
types. No new Baby Monster character-table claim is made in this packet.

## Incoming named-response work

The earlier [E1-E4 parity separator](../two_port_parity_separation_20261104/PROOF.md)
and [S1-S3 disjoint isolation](../simultaneous_module_reversal_20261104/PROOF.md)
remain valid. The [full quartic normal form](../session188_structural_bridges_20261104/response/PROOF.md)
provides the independent graph-cover computation used here. The
[three-state seed formula](../homogeneous_pair_response_20261103/PROOF.md)
retains both3^r and2^r channels and every physical port label.

The intervening Session188/189 work materially sharpens the old nested
boundary. The reviewed
[order-at-least-eight frontier](../session189_response_frontier_20261104/frontier_eight/PROOF.md)
first excluded prime-seven quotients with arbitrary children. Incoming
commits `c5725b3a7f` and `c87ca88cff` then added the reviewed prime-eight
closure and
[order-at-least-nine frontier](../session189_response_frontier_20261104/frontier_nine/PROOF.md),
with its [independent integration review](../session189_response_frontier_20261104/nine_integration_audit/REVIEW.md).
For a nonidentity XOR family of modules of the ORIGINAL source having the
same full named response, every highest changed original strong quotient
must now be prime of order at least nine. Arbitrary current-module sequences
and all-pairs injectivity among outputs of one source family are not covered.

The proof still uses exact equal-profile restoration of other highest nodes
at one physical query. Thus it genuinely extends the mechanism from this
chat; it does not deduce noncancellation from an unmarked graph bijection.
The original full cubics and quartics remain in their frozen packets; the
present Fourier replay adds a new complete order-five view without replacing
those named data by a graph-type count.

### Closing fetch through d617652d2e

The reviewed [dyadic thickening theorem](../session189_response_frontier_20261104/thickening/PROOF.md),
with its [independent review](../session189_response_frontier_20261104/thickening_audit/REVIEW.md),
now gives every nonempty source T a computable power-of-two divisibility
condition on N such that T[TT_N] is separated from every distinct
original-module XOR target of that inflated source. This changes the
physical object and its vertex set; it does not reconstruct T or provide
an H-preserving graph bijection. The contrast illustrates precisely why
the operation, names and source family must accompany a correspondence.

The [all-even-child transfer](../session189_response_frontier_20261104/even_transfer/README.md)
proves the complete contrast (Delta H,2 Delta H,0,...) modulo16 and its
odd ambient multiplier. Child profiles need match only between source and
target at each position, not across different children. This is ambient
transfer of a nonzero local gap, not a theorem that every local gap is
nonzero, and it does not advance the uninflated frontier beyond nine.
These incoming results reinforce the marked-operation interpretation of
the old even-graph basis without altering the bijection guarantees.
