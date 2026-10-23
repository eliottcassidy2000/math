# Beyond cyclic triples: two regular-five cells need a second shared-incidence coordinate

Status: PROVED finite family classification, with complete120-permutation controls and independently verified Hamiltonian counts. This does not classify arbitrary five-by-five regular blocks or all order10 tournaments.

## Carrier and gauge

Let A and B each be the regular tournament R5 on F5, with i→j when j−i belongs to{1,2}. For a permutation p of F5, orient the cross block A→B exactly on the five pairs (i,p(i)); orient every other cross pair from B to A. Call the resulting10-vertex tournament T_p.

The A vertices have outdegree3 and the B vertices outdegree6, so the two cells are intrinsically distinguished. Each cell's automorphism group consists exactly of translations C5. Therefore T_p and T_q are isomorphic exactly when q(i)=b+p(i−a) for some a,b in F5. This is the rotational double-orbit problem C5\S5/C5, not an affine-F5 classification: most permutations of five symbols are not affine.

There are eight ordinary isomorphism types. Burnside gives(120+16·5)/25=8: the identity pair fixes120permutations; a pair of nonidentity translations fixes five conjugating permutations; a pair with exactly one identity fixes none. Equivalently four affine-slope orbits have size5, and the remaining four orbits have size25. All are strongly connected, since both cells are strongly connected and there are cross arcs in both directions.

Converse exchanges the cells and reverses their internal orientations. Using i↦−i to restore the R5 charts sends p to q(i)=−p^(-1)(−i). Six types are self-converse; the affine slopes2and3 are exchanged. Thus there are seven converse-merged types, not eight.

## Two actual shared-incidence contractions

For a directed internal arc i→j in A, compare whether its matching-owned images p(i),p(j) have the same orientation in B. Define

J(p)=sum_(i→j in R5) sign_R5(p(i),p(j)).

There are ten terms, so J is even. This is the same actual-owner rectangle contraction used for the triple-cell repair, now on five-cell data. It is invariant under the allowed rotations and under converse.

The internal arc i→i+1 is intrinsic: among the two outneighbors of i in R5, its head shares exactly one outneighbor with i; the other outgoing arc shares none. These five short arcs form a directed5-cycle. Define the second, still shared-incidence contraction

K(p)=sum_(i in F5) sign_R5(p(i),p(i+1)).

K is odd. No arbitrary root or order is needed to distinguish the short arcs. K is not asserted to be globally converse-odd: its transformation depends on the orbit. Its purpose is to retain the orientation information that J loses on the one nonself-converse pair.

## Complete table and repair

|p in one-line notation|orbit size|parity|J|K|H(T_p)|
|---|---:|---:|---:|---:|---:|
|01234|5|+|10|5|8265|
|01243|25|−|4|3|8245|
|01342|25|+|−2|1|8221|
|01432|25|−|−4|−3|8201|
|02143|25|+|2|−1|8193|
|02413|5|−|0|5|8175|
|03142|5|−|0|−5|8175|
|04321|5|+|−10|−5|8245|

The table exhausts all120permutations under the proved gauge action, so its conclusions are complete within this declared family:

1. One permutation-parity bit leaves four ordinary types in each value.
2. Exact J classifies the seven merged types. J modulo8 does not:10collides with2, −10with−2, and4with−4. J modulo16 repairs all seven merged types. Sixteen is the first sufficient power-of-two modulus for this scalar; it is not the smallest arbitrary modulus (nine also works).
3. Retaining mod8 without changing modulus is possible: the pair(J mod8,K mod4) distinguishes all eight ordinary types. Because K is odd, K mod4 is only one additional binary coordinate. Each of the four J mod8 buckets contains two ordinary types, so one extra bit is necessary and sufficient given that first coordinate. Thus this actual shared-incidence bit repairs precisely the information the first guard forgets.
4. No single Z/8-linear combination aJ+bK can distinguish these eight types: J is even and K is odd, so aJ+bK always has the fixed parity b and at most four residues. This is a minimality statement only within this declared linear-compression family; an arbitrary nonlinear eight-symbol encoding is not excluded.
5. Hamiltonian count alone also loses information: H=8245 occurs in two distinct self-converse types, while H=8175 occurs in the exchanged pair.

The exact Jmod16 repair and the two-coordinate mod8 repair answer different questions: merged classification versus ordinary classification. Neither authenticates an arbitrary path-labelled complement partner without the actual path-position sidecar.

## Why this is the next boundary of the triple rule

The marked quotient theorem for at most six unmarked vertices works because only two C3cells can interact nonuniformly; a remaining regular5cell is unique and its attachments to coprime-sized triple/singleton cells are uniform. Once two five-cells interact, eight rotation-double-orbits appear even in the restricted matching-block case. The uniform-five-five hypothesis in the adaptive extension is therefore substantive.

For larger regular five-by-five blocks with row sum2or3, the present table does not apply. Their additional bipartite incidence structure must be retained and tested. The new rule is to enlarge the actual incidence bank when the carrier changes, not to extrapolate a successful triple-cell bit to every odd cell.

## Reproduction

Run python3 check_five_cell_barrier.py and python3 -O check_five_cell_barrier.py. The checker exhausts120permutations, verifies Aut(R5), orbit sizes, degree separation, converse pairing, both proposed guards and all64linear mod8compressions. A subset dynamic program computes the displayed directed Hamiltonian counts; the independent audit uses path DFS instead.
