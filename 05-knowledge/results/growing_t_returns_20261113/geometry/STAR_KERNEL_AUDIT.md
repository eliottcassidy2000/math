# Independent selective-star kernel review

**Accepted, scoped.** The reviewed files and exact hashes are in
`star_kernel_audit.json`. The review covers the final unique-minimum and
quantitative triangle additions as well as the original cut compiler.

The one-ear H polynomial partitions paths by whether the new vertex is
first, last, or internal. The F response has a separate singleton-new-vertex
sector H(U); deleting the final vertex of every other first path gives
the stated sum over incoming neighbours. Thus no empty second-path sector
is lost. The canonical second birth cancels the first-start field and
has readiness exactly the selected parent-start sum. Positivity of those
integer starts makes the stated disjoint unique-start cube faces both
necessary and sufficient.

The cut-network directions match source-side x=1: the source field charges
x_i=0, and arc j->i charges x_j=1,x_i=0. On each face, intersection and
union remain feasible and the mixed difference is exactly the stated
sum Q_ij+Q_ji. The inherited pair-deletion bound makes it strictly positive
for incomparable masks. The minimum-chain assertion is valid before the
stronger unique-minimum structure is used.

For a feasible vertex w, S_U(w)=1 forces U=C->TT_a with w in its initial
strong component and S_C(w)=1. The all-tail-out mask leaves the completed
core followed by that transitive tail. Reversing any chosen tail incidence
keeps every old odd-cycle pack and creates v->w->i->v. Therefore the
all-tail-out mask is the unique minimum; this reasoning does not assume
successive individual reversals remain monotone after cycles appear.

For each reversed tail i, the triangles(v,w,i) and(v,y,i) have complementary
Hamilton counts H(C-w) and H(C), respectively. Their pack families are
distinct because they prescribe different cycles through v. All baseline
packs remain. Thus the final quantitative bound

    N-N_min >=2r(H(C)+H(C-w))

is valid, including singleton C. It is an application of the inherited
odd-cycle formula, with that attribution retained in the reviewed proof.

The independent checker uses Python unbounded endpoint subset counts on
the literal A2 source, its three-arc parent, and the returned1691 graph.
It verifies the exact three changed old incidences, one actual birth,
parent state(153,1,27,1539), child state(1691,152,153,1692), and the
target/readiness cancellation. The exhaustive all-star and no-next-star
claims remain the root checker's declared finite universe, not a claim
that this audit re-enumerates broader repair menus. Normal and optimized
receipts match.

Reproduce `python -B .../geometry/star_kernel_audit.py`, then with `-O`.
