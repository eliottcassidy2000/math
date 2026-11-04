# Independent audit: Session190 tournament/even-graph correspondence

**Status: VERIFIED proofs in the stated scope; FINITE-EXACT independent checks through n=6.**

Audited packet: `../session190-bijection/`. Audit date: 2026-11-04. This audit did not import or execute either producer and made no changes to the audited packet. Its independent standard-library implementation is `audit.py`. Publication is the coordinating session's responsibility.

## Verdict

No substantive mathematical error found in the requested claims. The automorphism-even definition, orthostochastic normalization, full-flip return law, n=4 loss example, path-Walsh matrix, and pair-stabilizer argument all pass. Existing Walsh and free-marking results are credited, and the packet does not establish a rank-free natural bijection or reconstruction theorem.

One precision edit is recommended: replace “Its first degree moment gives the exit rate N-A(T,T)” with the explicit equation

    N - A(T,T) = 2 sum_G |E(G)| B(T,G).

The current wording “gives” need not mean equality, so this is clarification rather than a disproved theorem. The unscaled mean graph-edge count is half the exit rate.

A useful additional witness makes the off-diagonal loss claim concrete: at n=5, canonical tournament masks 2 and 11 are converses and have identical B rows, but, for the distinct target mask 4,

    A(2,4) = 3,     A(11,4) = 0.

Thus no function of the two endpoint B rows alone can recover all off-diagonal flip counts. This states the lost information precisely without treating an arbitrary finite table with externally supplied tournament names as an information-theoretic encoding problem.

## Primary-source check

Directly checked [Royle–Praeger–Glasby–Freedman–Devillers, arXiv:2204.01947v2](https://arxiv.org/html/2204.01947), §1 and §4. The parity counts edges whose assigned directions disagree after applying an automorphism. It is independent of the orientation; the sign is a homomorphism. Theorem 1.1 gives equality of unlabelled tournament and even-graph counts, and §1 identifies finding a natural unlabelled bijection as open in that paper. The packet accurately distinguishes this from even automorphism-group order and does not claim a comprehensive 2026 status search.

The journal reference and DOI agree with the [author-deposited record](https://arxiv.org/abs/2204.01947): J. Algebraic Combin. 57 (2023), 515–524, DOI 10.1007/s10801-022-01197-0.

Independent counts at n=1,...,6:

    tournament types          1, 1, 2, 4, 12, 56
    automorphism-even types   1, 1, 2, 4, 12, 56
    even group-order types    0, 2, 4, 11, 34, 148

K2 and P3 provide the claimed sign controls. For K2, edge-direction reversal has sign -1 but the edge-set permutation has sign +1. For reflection of P3, these signs are +1 and -1. The n=2 even-group-order interpretation therefore fails.

## Proof audit

### Signed basis and canonical stochastic kernel

The signed graph-copy polynomial is transporter-independent exactly for automorphism-even graph types. Distinct supports are orthogonal under uniform labelled tournaments, giving squared norm r_G. Summing over tournament orbits gives the weighted identity in the packet. The invariant Fourier subspaces are one-dimensional for even graph orbits and zero for odd graph orbits. Together with the tournament-orbit indicator basis, this makes Q square and orthogonal.

Changing the representative, its initial edge directions, or the reference orientation gauge only changes the normalized signed column by an overall sign after the coordinate changes are applied consistently. Tournament relabelling leaves the class polynomial unchanged. Hence B=Q squared entrywise is independent of those choices. Both stochastic marginals, the empty-graph entry B(T,empty)=r_T/2^N, and the Hall-support matching argument follow. The exact Birkhoff decomposition supplied at each order was independently reconstructed from all its permutation matrices and rational weights.

The canonical n=3 choice gives Phi rows (1,1) and (1,-3). An inherited description of the P3 polynomial may use the opposite column sign; this is harmless and is precisely the convention dependence already declared in the packet.

### Full-flip response and its boundary

A acts on the G character by N-2|E(G)|. Its reversible weights are pi_T=r_T/2^N. Inverting the weighted Walsh basis gives the factor pi_U in the off-diagonal expansion, with no missing orbit-size or automorphism factor. The generator for independent rate-one arc flips is A-NI, hence the decay exp(-2|E(G)|t). The exact return spectral measure is B, aggregated by graph-edge count. The m=0 convention and N=0 discrete normalization exception are correctly stated.

The independent test compares every class-transition entry for powers 0,...,4 and at exp(-2t)=1/3 against a separate labelled-target distribution. In the latter, a target at Hamming distance h has probability 2^(N-h)/3^N. This is not a comparison between two copies of the spectral calculation.

At n=4, masks (0,2,4,5) and graph masks (0,3,7,30) give exactly the displayed B matrix. Masks 0 and 4 are self-converse with trivial automorphism groups; mask 4 is strongly connected. Their common continuous-time return law is

    3/8 + (1/2) exp(-4t) + (1/8) exp(-8t).

Their Hamiltonian-path counts are respectively 1 and 5. Independent path lists are saved in `audit_results.json`. The corresponding path-complement rows are (0,0,1,0) and (1,1,2,1), so the row sum and actual incidence response are genuinely different. This is a finite counterexample, not a reconstruction claim.

Independent B ranks, numbers of distinct B rows, and numbers of distinct degree-aggregated return rows are all

    1, 1, 2, 2, 10, 34.

At n=5 and n=6, equality of B rows is exactly equality up to converse, as a finite census statement only.

### Path-complement formula and stabilizers

For each undirected spanning path U at n>=2, exactly one traversal can be directed in a tournament. Adding the two traversal indicators cancels their odd Fourier terms and leaves coefficient 2^(2-n). Complementing outside U changes character H by (-1)^|H\U|. Orthogonality then forces the difference support G xor H to lie inside U with even cardinality, giving the stated matrix formula. No factorial or path-reversal factor is missing.

The formula is symmetric because G\U=H\U whenever G xor H is inside U. It preserves Fourier parity, not Fourier degree. Its n=3 invariant-basis matrix is exactly

    [ 3/2, -3/2 ]
    [-1/2, -3/2 ].

This basis is not orthonormal, so the displayed matrix need not be ordinarily symmetric. Weighted symmetry is consistent with labelled self-adjointness.

The common directed arcs of T and T_p form the directed spanning path p. Any common automorphism fixes that path pointwise, so the stabilizer is trivial and p is recovered from the pair. Dividing ordered incidences by n! gives D(T,U)/|Aut(T)|, explaining both integrality and reversal symmetry. For n>=3, endpoints are distinct, and a permutation interchanging them would preserve the common directed path and thus be the identity. Hence the diagonal quotient counts are even. The n=2 exception is necessary and correctly retained.

Stationary weights for path-rerandomizing D/H are r_T H(T). They are not the full-flip weights r_T, and they do not turn repeated rerandomized steps into the original fixed-path involution.

## Inheritance and nonclaims

Targeted inspection confirmed the cited prior work exists and supplies the claimed mechanisms:

- `math/01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`, Theorem W: the signed invariant basis, both orthogonality identities, and dyadic determinant.
- `math/05-knowledge/results/tournament_even_bijection_20261011c.md`: SCC/Lyndon/component codec, with explicit rank choices at irreducible atoms and an express non-naturality boundary.
- `math/01-canon/theorems/THM-4725-no-uniform-dilation-rule-for-prime-circulants-signed-set-burnside-and-the-p-mod-4-dichotomy-with-the-paley-rule.md`: scoped prime-circulant constructions/obstructions.
- `session34-representations/NOTE.md`, §§4–5: full-flip Fourier degree and free marked-path action.
- `math/01-canon/theorems/THM-848-hamiltonian-path-walsh-needle-stalk.md`: full-arc-cube H response and missing-state distinctions.
- `session33-metagraph-memory/NOTE.md`: literal fixed-path complement versus intermediate rerandomization.

This audit does not re-prove all cited prime-circulant or codec theorems, validate a natural-bijection claim, or establish asymptotic efficiency. The elementary species obstruction has the correct direction: no equivariant map from all labelled even graphs to tournaments exists for n>=2 because of the empty graph's full symmetric stabilizer. It alone says nothing analogous about every possible reverse-direction map.

The homological determinant identity and the tree-center characterization were also checked at proof level. This audit did not independently rebuild the producer's full homology computation and makes no such claim.

## Reproduction and exact scope

Run from the enclosing workspace:

    python3 session190-bijection-audit/audit.py --max-n 6
    python3 -O session190-bijection-audit/audit.py --max-n 6 --out session190-bijection-audit/audit_results_optimized.json
    cmp session190-bijection-audit/audit_results.json session190-bijection-audit/audit_results_optimized.json

Both executions pass and are byte-identical. No software was installed. No producer functions or repository modules were imported. Explicit checks remain active under optimization. Complete graph/tournament orbit enumeration uses independent endpoint transport; signed basis values on every labelled tournament use a coefficient-vector Walsh transform. The path matrix through n=4 is checked by direct double sums, not the producer's butterfly method.

Checks include 1,847,571 labelled polynomial-invariance entries, 19,812 exact full-transition entries, and 1,333 representative path-pair stabilizers at n=2,...,6. An extra n=1 stabilizer check is saved separately in the total per-order output. All 4,164 full labelled path-Walsh entries at n=2,3,4 are checked. No full labelled path-Walsh matrix at n=5 or n=6 is claimed. The audit also excludes endpoint-swapping permutations on every representative incidence at n>=3.

The packet's eight listed SHA-256 hashes were checked successfully before this audit's closure. `SOURCE_SHA256SUMS` pins the source files actually audited; this audit's own files are pinned by `SHA256SUMS`.
