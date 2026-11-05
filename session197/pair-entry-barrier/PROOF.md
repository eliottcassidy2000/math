# Two physical pair-entry states with no ordinary predecessor

Session 197. This is a separate corollary to the reviewed pair-entry bridge. It does not alter that snapshot or assert a scalar-spectrum obstruction.

An **ordinary move** means exactly an allowed one-seam gluing with transitive gadget-root deletion, one gadget path ending at its root, and the distinguished root moved to that gadget root. Its three-coordinate transfer is

    (h,e,R) -> (g h−t e+L R, h−[L=1]e, Qh+BR).

The input may be any actual nonempty rooted tournament. The statement below is not about every port operation or every possible rooted graph construction.

## 1. The universal ordinary-output cone

The response-cone theorem in `session197-lowcycle/ARITHMETIC_INVERSE.md` proves that every ordinary output satisfies

    R <= 2(H−e).                                    (1)

For self-contained use, its short calculation is recalled. The gadget recurrences give Q≥g, B≥L, and B≤2L−2 for L≥2. Set D=2g−Q. Even L preserves D; odd L=2l−1 changes it by 2l−2(l mod 2)−B(l). Induction gives D≥1 always and D≥3 at odd L≥3. For L≥2,

    R'−2(h'−e')=(Q−2g+2)h+2t e+(B−2L)R.

For even L this is at most h−2R≤0. For odd L it is at most −h+2e−2R≤0. For L=1 it is h−R≤0. Here every physical rooted graph has R≥h and 0≤e≤h. The singleton seed need not itself satisfy (1), but every ordinary output does.

## 2. Every first-clone switch crosses that cone

Let T be strong C1, let a have sole out-neighbor b, and expose a,b. Add f beating both and retain the two clones a,f. Write the original responses as h,e_a,e_b,R_a,R_b. The reviewed pair-entry identities give the first-clone profile

    H_*=h−e_b+R_b,
    e_*=h−e_a−e_b,
    R_*=h+R_a+R_b.                                  (2)

There is one further elementary adjacency identity:

    R_b=2(h−e_a)−e_b.                               (3)

To prove it, remove a and call the remaining graph V. Every Hamilton path not ending at a splits immediately before a→b into an arbitrary prefix and a b-start suffix, so h−e_a=R_b(V). In a b-start two-cover of T, a must end one of the two paths, since its only out-neighbor b is already the starting vertex of the reserved path. Appending a to the unreserved slot gives R_b(V) choices. Appending it to the reserved slot gives another R_b(V), except that a cannot follow the singleton b; this excludes H(V−b)=e_b. This proves (3).

Substitution into (2) gives the exact excess

    R_*−2(H_*−e_*) = R_a−h+e_b >= 1.               (4)

The first inequality uses R_a≥h, by cutting each Hamilton path just before a. Also e_b=H(T−{a,b})≥1: a must be penultimate in a b-ending Hamilton path, and every path on the remaining vertices can be followed by a,b.

Therefore **no first-clone state in (2) is an ordinary output at any order**, even if its hypothetical input is an arbitrary rooted tournament rather than a C1 graph. In particular it cannot be the endpoint of a nonempty word of ordinary moves.

On the same actual graph, the last clone a has start count one. It is an ordinary output from T by a singleton birth beating b, under the explicit clone-label exchange. It therefore lies in the cone. Both clones are physically exposed, so changing the marked root genuinely crosses the response barrier without changing H or the graph.

This is a state-level separation. The first-clone scalar count is of course attainable, and is the same count as the last clone. No scalar impossibility is inferred.

## 3. Being inside the cone still does not give a physical predecessor

The reviewed nonadjacent-pair construction gives the actual graph

    rows [18,20,25,3,8], root 0,
    (H,e,R)=(13,2,18).

It lies strictly inside (1), since 18<2(13−2)=22. Nevertheless it has no ordinary predecessor.

For L≥2, the output ending count equals the old H. Here that would require old H=2, impossible because every tournament has an odd Hamilton-path count. For L=1 the inverse is unique:

    old (h,e,R)=(R'−h'+e', R'−h', h'−e')=(7,5,11).

No tournament has H=7. One short finite reduction is enough: a nontrivial strong H=7 component has order at most four by 5^(n−1)≤7^3. A tournament with H=7 has exactly one such component, since 7 is prime and SCC path counts multiply. The only strong tournaments of orders three and four have H=3 and H=5, respectively. The verifier independently enumerates all labelled tournaments through four to check this last finite statement. Thus no old graph exists at any order.

The familiar parity theorem and SCC product are used only for this concrete obstruction, not as sufficient conditions for realization.

## 4. Scope and verification

The universal output cone is due to the separate lowcycle arithmetic packet. The additional conclusions here are the exact first-clone excess (4), the all-order no-predecessor statement for that entire physical family, and the inside-cone H13 example with no ordinary predecessor.

The two-port root switch and nonadjacent graph are already independently reviewed in `session197-hybrid-constructions/pair-entry-bridge/PROOF.md`. Its snapshot remains unchanged. The verifier rechecks the adjacent identity and both cone sides on physical C1 controls, the explicit H13 graph, and the finite strong-order cutoff for H7. These controls support the formulas; the all-order statements follow from the displayed identities and reductions.

Neither (1) nor this obstruction solves the uniform inverse problem. A successful inverse procedure must include physically verified pair-entry or rerooting branches, rather than treating the cone or a formal integer inverse as a realization theorem.

## Exact dependency snapshots

- `session197-hybrid-constructions/pair-entry-bridge/PROOF.md`: `fcff7ca192c4055243aeab2344582686214aa7aa25504d1233f3b6b33465a856`.
- `session197-lowcycle/ARITHMETIC_INVERSE.md`: `4ac1727eb27c85f8a39f5598eb9f0e1b020bc4020cbca0a94e6288cb0676c3b7`.
