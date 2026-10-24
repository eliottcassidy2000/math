# Independent audit: marked parent growth and unpointing

24 October 2026. **PASS**, with the scope and one repaired wording issue below.
This is an audit of the Session50 marked-incidence transfer, not a claim of a
closed A000568 recurrence from the bare blue/black metagraph.

## Audited sources and inheritance

The tested producer is archived as `source/transfer.py.txt`, SHA-256
`fd3db75c8425ff4fc2946835e7e4e55ffba905ae4d66fb1a6efd5dda61c9c4ba`.
The final proof is `source/THEOREM.md.txt`. `source-pins.json` records its hash
and all inherited source pins. The inherited sources are byte-identical to
commit `f26bc37f4e77f72a7ce4e4f9eabf2b271acc3390` of `eliottcassidy2000/math`.
All copied source documents and code are `.txt` archives, not new active canon.

THM-781, `01-canon/theorems/THM-781-hamiltonian-path-inverse-metagraph-fibre.md`,
already identifies native-word fibres with Hamiltonian-path orbits.
THM-796, `01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md`,
already supplies compatible observer growth and non-lumpability boundaries.
The October20 `metagraph_growth_20261020/growth.md` and `enumeration.md`
already establish native growth, transitive inflation, and explicit failures
of specified unmarked summaries. The elementary orbit/groupoid counting
mechanism is not asserted to be new, and no literature-priority or
asymptotic-complexity claim has been audited or established.

## All-order mathematical review

1. **Marked extension states are correctly typed.** An attachment subset of an
   order-n parent, modulo its actual automorphism permutation action, is
   precisely an order-(n+1) vertex-pointed tournament. Its subset stabilizer is
   the pointed child's stabilizer, not its full automorphism group.
2. **The transport is genuinely constructible from parents.** To change the
   root, delete the new candidate root, identify the remaining order-n
   tournament with the parent atlas, and transport its attachment mask.
   Two identifications differ by a parent automorphism, so the resulting
   subset orbit is independent of that choice. No unpointed child canonical
   label is needed.
3. **The block proof is correct.** In the proof only, fix a child whose vertex
   orbit sizes are r_1,...,r_k. Its D block has every row equal to
   (r_1,...,r_k). Hence D^2=N D, each block has trace N and rank one, and
   A(N)=tr(D)/N. Orbit-stabilizer gives h_p D_pp=|Aut(U)| and
   D_pq/h_p=D_qp/h_q. The weights for self-adjointness are 1/h_p.
4. **Converse is correctly handled.** The full pointed involution C commutes
   with D. A self-converse block contributes N to tr(CD), all other blocks
   contribute zero. Thus SC(N)=tr(CD)/N. The plus/minus projectors
   (D/N +/- CD/N)/2 are mutually orthogonal in the 1/h weighted inner product.
5. **The endpoint quotient is correct.** The new root is a Hamiltonian end
   precisely when L_T(S)>0. Its native observer multiplicity is L_T(S)/h_p
   because the pointed stabilizer acts freely on ordered Hamiltonian paths.
   Restricting root choices to eligible endpoints and normalizing gives the
   stochastic rank-one blocks of R, so tr(R)=A(N).
6. **The anti-endpoint correction is essential.** C does not preserve eligible
   ends, already at N=2. Reversing the child and choosing an eligible end in
   its converse gives K, with K^2=R, RK=KR=K and tr(K)=SC(N). The implementation
   stores extra rows for all full pointed states; the square endpoint domain
   used by the theorem is their restriction to admissible rows and columns.
7. **The symmetry defect and bounds are correct after precise typing.**
   Delta=sum_p(D_pp-1)>=0 and A=(|P|+Delta)/N. Tournament automorphism groups
   have odd order, so all root orbit sizes are odd and Delta is even. Delta
   vanishes exactly when all children are asymmetric. A root bound b_p is
   valid when its cell contains the entire root orbit.

### Repaired issue: an equivariant partition is not enough

The initial bound described an arbitrary “isomorphism-invariant vertex
partition.” An unlabelled equivariant partition may have its blocks permuted:
the discrete partition on the directed triangle is an immediate counterexample
to using singleton cell size as an orbit bound. The final proof now explicitly
requires every cell to be setwise invariant under every automorphism,
equivalently a coarsening of the orbit partition. Stable canonical degree
coloring meets this condition. No predictor count or algorithm changed.

## Independent finite validation

`independent_audit.py` does not use the producer's canonicalizer, automorphism
finder, path dynamic program, deletion function, or sparse multiplication in
its controls. It calls only the producer predictor under a hard runtime guard
forbidding child-order calls to its canonicalizer. Independent canonical keys
are minima of upper-triangle arc codes over literal permutations.

Both normal and optimized Python pass **10,503 explicit checks**, including:

- complete labelled parent universes through order 5;
- complete labelled child universes through order 5;
- all normalized Hamiltonian-path word presentations through order 6;
- every pointed state's full rooted-permutation comparison through order 6;
- every child automorphism group and Hamiltonian path through order 6;
- exact D and C, pointed isotropy, L/h fibres, endpoint support, R and K;
- weighted plus/minus projector identities, balance and traces;
- the N=2 converse restriction and N=3 pointed-division hostile controls.

The predictor receives independent parent atlases, not published target
counts. The expected sequence constants in the producer occur only in its
outer post-computation control routine, which this independent audit never
calls. No child canonicalization or target-count input enters `predict`.
Child degree refinement computes only an invariant upper bound, not a
canonical child label.

| N | Full pointed | Endpoint pointed | A(N) | SC(N) | tr(D) | Delta |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 1 | 1 | 0 |
| 2 | 2 | 1 | 1 | 1 | 2 | 0 |
| 3 | 4 | 2 | 2 | 2 | 6 | 2 |
| 4 | 12 | 7 | 4 | 2 | 16 | 4 |
| 5 | 48 | 33 | 12 | 8 | 60 | 12 |
| 6 | 296 | 236 | 56 | 12 | 336 | 40 |

## Full order-eight replay, not an independent order-eight census

The frozen producer was replayed through order 8 in normal and optimized
Python. Both complete JSON receipts and displayed outputs agree exactly;
both pass **2,641,571 explicit checks**. At N=8:

- 456 parent classes; 54,256 full pointed states; 50,709 endpoint states.
- tr(D)=55,040, tr(CD)=1,408, hence A(8)=6,880, SC(8)=176 and merged(8)=3,528.
- Delta_8=784. Naively dividing the pointed count by eight gives 6,782,
  undercounting by 98.
- Maximum order passed to prediction canonicalization: 7.
- Kernel SHA-256: `dd9d54cd3bd0f1573d4968a315429a7f6d4d627b5db1f5cee41ec7c5ea772775`.

The producer's separate controls exhaust native word cubes through 7 and full
permutations through 6. Its order-eight check does not independently enumerate
all 2,097,152 native words or all labelled eight-vertex tournaments. The
reported order-eight native-word total is an exact mass identity, not a claim
that the full cube was separately classified.

The earlier `pre-witness-*` replay used the archived earlier predictor whose
only difference from the final source is omission of a four-line
`refinement_gap_witness` reporting field. It is preserved for provenance; the
normal/optimized final receipts are the authoritative replay evidence.

## Refinement hostile and its infinite family

`refinement_witness_audit.py` independently checks all 7! permutations of
(74,84,97,22,37,11,56). Its automorphism group has order 3, with orbits
{0,1,2}, {3,4,5}, {6}; degree refinement leaves one cell of size 7. Its three
pointed states therefore contribute 7 to tr(D), not the refinement candidate
21. The excess 14 exactly explains the global order-seven candidate 3206
versus the true 3192. At order eight, the producer's upper trace 55068 exceeds
55040 by 28.

The inherited uniform transitive inflation theorem proves
Aut(T[K_k])=Aut(T). Clone position i has degree 3k+k-1-i, so degree refinement
has k cells of size 7. Counts toward position j are 3+1_(i<j), independent of
the base vertex, so no further refinement occurs. The k clones above fixed
base vertex 6 remain fixed. The refinement replacement therefore fails at
every order 7k, k>=1; the contribution from all this child's pointed states
has excess 14k. Independent finite clone-cell and homogeneous-pair checks at
k=1,...,6 agree. The all-k conclusion uses the displayed formulas and the
inherited automorphism theorem, not extrapolation from those checks.

## Exact scope boundary

The sufficient input is a complete labelled parent atlas with actual parent
automorphism actions and the attachment/deletion incidence it supports. The
algorithm computes cross-parent overlap transport from this information.
It does not establish that naked colored metagraph adjacency, scalar counts,
automorphism orders alone, or a bounded parent signature determine that
transport. The full D still retains classification-strength information. The implementation
constructs every D row: the diagonal-only counting reduction and singleton
root certification are proved diagnostics, not measured test-skipping
optimizations or a demonstrated speedup.
The remaining problem is an all-order closed compression preserving the
required return multiplicities, not an omitted division or trace trick.

## Reproduction and packaging

From any working directory:

    python /path/to/independent_audit.py --receipt /tmp/independent.json
    python -O /path/to/independent_audit.py --receipt /tmp/independent-O.json
    python /path/to/source/transfer.py.txt --max-order 8 --receipt /tmp/order8.json
    python -O /path/to/source/transfer.py.txt --max-order 8 --receipt /tmp/order8-O.json
    python /path/to/verify.py --deep
    python -O /path/to/verify.py --deep

`manifest.json` closes the packet: every regular payload file is listed and
hashed, and unexpected files or symlinks are rejected. The verifier needs no
repository, network, external source path or third-party package. Its deep
mode reruns the independent audit into a temporary external directory. The
recorded relocation test copies the packet to a new temporary tree, changes
working directory, and reruns normal and optimized checks there.
