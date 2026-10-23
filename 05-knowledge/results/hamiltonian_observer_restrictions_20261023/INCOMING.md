# Late incoming work: distinguish the carrier, then retain the seam

Status: **PROVED** elementary comparisons and phase-cycle extension below; **INHERITED / FINITE-EXACT** where explicitly attributed. Integrated commits `1eb9c41c51` and `9928396ee6`, which arrived while the first local checkpoint was being published.

## 1. Two observer graph partitions are incomparable

Incoming [Session43 restricted coupling](../session43_observers_phase_guards_20261023/coupling/THEOREM.md) uses the backward graph `B_x` on all `n` observer positions. Our [native incidence descent](coupling/NATIVE_CONNECTEDNESS.md) uses the compressed graph `G_x` on `n-1` positions, with edge `(i,j)` corresponding to backward native arc `(i,j+1)`.

Both retain exactly the same labelled free bits. After ordinary graph isomorphism, however, they are different observations. In zero-based positions:

| Compressed edge set | Full backward edge set | Shapes |
|---|---|---|
| `{01,02}` | `{02,03}` | `P3` becomes `P3` |
| `{01,12}` | `{02,13}` | `P3` becomes two disjoint edges |
| `{01,23}` | `{02,24}` | two disjoint edges become `P3` |

The first two rows show that compressed graph class does not determine full backward graph class, for every `n>=4`. The first and third show the reverse failure for every `n>=5`. Additional vertices are isolates on the relevant side. Thus from order five the two quotient partitions are incomparable, although the labelled bit transport is bijective.

This resolves the apparent complement discrepancy. The incoming *full-graph* complement of the empty graph is the complete `n`-vertex graph, which violates every forced path bit and is outside the observer subset. The isolated transitive block survives that pairing. Our *compressed-graph* complement flips only free bits, stays inside the observer cube, and pairs the empty source with a nonempty complete compressed source. It joins the two components of our incidence relation.

The identical order-four nontransitive eigenvalue `8/15` is therefore not evidence that these operators agree. Their graph-state masses already differ: `(3,2,1,1)` incoming versus `(3,3,1)` here. Our all-size two-component theorem is not silently transferred to the incoming full-graph partition.

## 2. A genuine native assembly cycle and its representative cut

Incoming [observer bundles and native ports](../session43_observers_phase_guards_20261023/native-ports/OBSERVER_BUNDLES_AND_NATIVE_PORTS.md) separates complete tournament automorphisms from assemblies of internal charts. For two actual directed triangles, its candidate paths are genuine exactly when

    (a,b) in F3^2,  b-a in {0,1}.

There are six such paths, a union of two diagonal cosets. These internal charts need not preserve the entire tournament, so this example does not contradict our single-orbit coset compiler. The incoming branch-aware compiler retains each difference branch and every simultaneous phase lift.

**Phase-cycle proposition.** For an integer `p>=3`, set

    S_p={(a,b) in (Z/p)^2 : b-a in {0,1}}.

Join two distinct elements when just one coordinate changes. The induced graph is precisely `C_(2p)`, with cyclic order

    (0,0),(0,1),(1,1),(1,2),...,(p-1,p-1),(p-1,0).

Proof: every fixed row has exactly its diagonal and next-diagonal vertex, as does every column. These are the displayed two neighbors, and the list visits all `2p` states. For odd `p`, the incoming circulant-cell construction with translated boundary palette `{0,1}` realizes these as actual Hamiltonian path assemblies; the combinatorial proposition also holds for even `p` without that tournament-realization claim.

Choose integer representatives `a,b in {0,...,p-1}` and impose

    a-b in {0,-1}.

Exactly one state, `(p-1,0)`, is deleted. The graph becomes `P_(2p-1)`. Equivalently, for any modulus `L>p`, the residue guard `a-b mod L in {0,L-1}` has this effect: the known difference range prevents aliases. This is a bounded-range assertion; for `L=p` the wrap state survives and the cycle remains.

At `p=3,L=8`, this is the incoming synthetic mod-8 example: six genuine assemblies become a five-state path. The two difference branches now have respectively three and two surviving lifts. The quotient difference alone no longer determines guard acceptance. The remaining state is the actual simultaneous phase, with its chosen representatives.

There is an exact all-modulus refinement. The only possible representative differences on `S_p` are `0,-1,p-1`. Therefore the wrap survives the modulus-`L` guard **iff `L` divides `p` or `p-1`**. A fixed mod-8 guard cuts the odd sizes 3, 5 and 7, but the wrap returns at size 9, and then whenever `p` is 0 or 1 modulo 8. This is an explicit mod-8/mod-9 interaction on the declared phase carrier, not an explanation of the unrelated numerical spectral gap `8/9`.

The cover also retains a small holonomy. Write `D_a=(a,a)` and `U_a=(a,a+1)`. Quotient by simultaneous translation of both phases. There are two vertex types and **two parallel edge types**, one changing `b` and one changing `a`. The base circuit lifts as `D_a -> U_a -> D_(a+1)`, advancing the phase by one. It takes `p` such circuits, or `2p` actual moves, to return. Collapsing the parallel edges to a simple `K2` loses that phase transport. These are valid path-assembly moves across different difference branches, not automorphisms preserving the entire two-cell tournament.

This supplies a precise finite cut mechanism. It is not identified with the electrical mod-8 guard, nor with the five-state valuation quotient. Its connection to the [gain seam](history/VALUATION_CYCLE_GAIN_LIFT.md) is the exact information operation: a quotient identifies representatives, while an actual lifted condition can distinguish them.

## 3. Finite modulus plus a range bound can be sufficient

The incoming [finite-degree guard proof](../session43_observers_phase_guards_20261023/FINITE_GUARDS_AND_DEGREE.md) gives an independent algebraic instance of the same bounded-range repair. For a polynomial map of degree `D`, exact cyclic equivariance modulo `L>D+1` forces its integer monomial charges to equal the required charges, not merely agree modulo `L`. For a Keller map of weights `(1,-1)`, that symmetry forces a diagonal linear map. Simultaneous exact `C8,C9` equivariance has modulus 72, so the degree bound is 70; the shear of degree 71 shows the threshold cannot be dropped.

This is inherited scoped rigidity, not a proof of the planar Jacobian conjecture. The same incoming Gaussian analysis also shows why a single finite angular modulus cannot certify all polynomial powers: charge bandwidth grows with the power.

For our compiler, the useful obligation is now sharper. To replace an actual integer guard by residues, prove either a sufficient range bound or an exact coordinate-factorization theorem. Without one, keep the integer or kernel coordinate. In the phase-cycle example the difference lies in `[-(p-1),p-1]`; in the gain example an integer valuation remains unbounded and cannot be replaced by the same finite residue logic.

## 4. Restricted spectral evidence agrees once the operations are typed

The incoming coupling proof independently supplies the same conditional-variance criterion and positive-weight comparison as our restriction lane. It additionally computes actual full-cube filters by Hamiltonian count modulo eight: at order five the residue-seven filter has gap `226/315`, whereas residue five has gap `1657897/1801800`. One is below and one above `8/9`. These are inherited exact finite calculations, not new computations here and not electrical-mod-8 statements.

Together the packets separate five operations: retaining an observer, changing its graph quotient, coarsening a statistic, conditioning on a hard guard, and sampling only at guarded return times. The proof mechanism and the clock depend on which operation is used.

The new `check_incoming.py` verifies the phase-cycle/cut proposition for `p=3..40`, every modulus `L=2..46`, the labelled two-edge return clock, both quotient counterexamples and all local source/target degree-sequence distinctions. A sibling independently audited the general alias criterion and passed 1,786 checks across these ranges. This does not revalidate the incoming Gaussian or Jacobian results.
