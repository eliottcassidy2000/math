# Sharp labelled odd girth in the remaining residue class

Independent verification and proof, Session25 continuation, 2026-10-21.

**Theorem.** For every r>=1, the labelled tournament nonpath-complement graph L_(4r+1) has odd girth 2r+1. Together with the existing even-order construction and parity obstruction, this completes the labelled odd-girth classification for n>=3: L_n is bipartite for n=0,3 mod4, and has odd girth 2 floor(n/4)+1 for n=1,2 mod4. This says nothing analogous about the unmarked isomorphism quotient.

## 1. Hamilton-path odd cover

Put k=2r+1 and work first on {infinity} union Z_(2k). For i=0,...,k-1 use the Walecki cycle

 C_i = (infinity, i, i-1, i+1, ..., i-(k-1), i+(k-1), i-k, infinity).

The finite edges of C_i have endpoint sums 2i-1 or 2i modulo 2k, and are precisely all non-loop pairs with those sums. Its two infinity edges have other endpoints i and i-k. Therefore the k cycles partition K_(2k+1).

Delete a=0 and b=k-1=2r. The edge ab belongs to C_r, so this cycle becomes a Hamilton path on the n=2k-1 remaining vertices. In every other cycle, a and b are nonadjacent; deletion produces two nonempty path components, allowing an isolated vertex as a path. Its two former neighbors of a lie in different components, as do its two former neighbors of b.

For j=r+1,...,2r define

 A_j=j, B_j=(j+r) mod k=j-r-1, D_j=2j mod k,
 e_j={2j-1,2j}.

Both endpoints of e_j are neighbors of a in C_(A_j), and neighbors of b in C_(B_j). Indeed the finite neighbors of x in C_i, when neither is a formal self-pair, are 2i-1-x and 2i-x modulo 2k. Thus adding e_j to each of these two deleted cycles joins their two components into a Hamilton path.

The pairs {A_j,B_j} partition all cycle indices except r. The edges e_j are mutually vertex-disjoint and avoid a,b. Moreover the original occurrence of e_j lies in C_(D_j), by its odd endpoint sum 4j-1. The three indices A_j,B_j,D_j are distinct: D_j=A_j would require j=0 mod k, and D_j=B_j would require j=r mod k, both impossible. Hence the repair never duplicates an edge within one path.

We obtain k undirected Hamilton paths whose edge multiset is K_n+2M, where M={e_j} is a matching of size r. The three occurrences of e_j are exactly in paths A_j,B_j,D_j. It remains essential, and is not automatic for an arbitrary odd cover, to make them orientation-compatible.

## 2. Ordering lemma for paired triples

**Lemma.** Let X consist of one distinguished element and disjoint pairs P_h={a_h,b_h}. For each pair choose c_h outside P_h and attach any prescribed binary incidence labels d_(h,v) to the three elements v in P_h union {c_h}. There is a total order on X and a binary variable x_v for every v such that, in the ordered triple, the three values d_(h,v)+x_v alternate 0,1,0 or 1,0,1.

All arithmetic in this section is over F2. If q_h(v) is the rank 0,1,2 of v in its ordered triple, the requirement is

 x_u+x_v = d_(h,u)+d_(h,v)+q_h(u)+q_h(v).       (1)

Only two of the three pairwise equations in each triple are independent.

**Proof.** Contract each pair P_h to a node h; retain the distinguished element as a base node. Draw one directed edge h -> owner(c_h), or h -> base when c_h is the distinguished element. Each nonbase node has outdegree one, the base has outdegree zero, and there are no self-loops. Each weak component is consequently either a tree ending at the base or a directed cycle with trees feeding into it. A directed 2-cycle is allowed and corresponds to two parallel edges in the underlying undirected multigraph.

For each directed cycle, call the endpoint of P_h selected by the preceding cycle node's c-value the active element v_h. Call the other endpoint its spare s_h. The outgoing c_h is precisely the active element of the next cycle pair. In every cycle triple use the two equations (1) between v_h and s_h, and between v_h and c_h. Elsewhere use any two equations including the pair's internal equation. The resulting undirected constraint graph is a pseudoforest: contracting its disjoint pair edges recovers the functional graph. Its unique cycle in each cyclic component consists of the edges v_h--c_h. Consequently all equations are solvable exactly when, for every such cycle,

 sum_h [d_(h,v_h)+d_(h,c_h)+q_h(v_h)+q_h(c_h)] = 0.       (2)

There are no compatibility conditions in tree components.

For two selected elements of a triple, their rank parities agree exactly when the third element lies between them. Thus the summand in (2) is toggled precisely when its spare changes from outside the interval between v_h and c_h to inside, or conversely.

Begin with every nonactive element before every active element, in arbitrary orders within these two groups. Each spare is then outside its cycle interval. For a cycle whose parity in (2) is wrong, take one of its spares and reinsert it strictly between its two active neighbors in the total order; for example insert immediately after whichever of those two occurs first. Those neighbors are distinct because there is no functional self-loop. This toggles exactly that cycle obstruction. The moved spare occurs in no other cycle triple: it belongs to just its own pair, and if it were the c-value of a cycle node it would have been declared active. It may occur in triples on feeding trees, but those have no cycle conditions. Relative order of all active elements is unchanged, and unmoved spares remain before all active elements. Therefore repairs for distinct cycles do not interfere. After repairing each obstructed cycle, solve (1) by propagation on a spanning forest. This proves the lemma. QED.

## 3. From the lemma to an actual directed closed walk

Give every constructed undirected path an arbitrary initial traversal direction. For a repeated edge e_h, let d_(h,v) record its direction in path v relative to one fixed orientation of that edge. Apply the lemma to the paired triples {A_j,B_j,D_j}.

Suppose path v occupies position p(v), zero-based, in the resulting order. Reverse that path if and only if x_v+p(v)=1. At the occurrence of e_h with triple rank q_h(v), its displayed direction becomes d_(h,v)+x_v+p(v). To make it a directed edge at that time, its required direction in the initial tournament is

 d_(h,v)+x_v+p(v) + (p(v)-q_h(v))
 = d_(h,v)+x_v+q_h(v).

The subtraction counts earlier moves in which this edge was absent and thus toggled. Equation (1) says exactly that the displayed expression is independent of which of the three occurrences is used. Every singly occurring edge independently determines its initial direction by the same rule. We thus obtain a well-defined initial tournament, and every path is directed at its claimed step.

Every edge occurs an odd number of times in this odd-length sequence, so it is toggled an even number of times and the walk closes. The existing lower bound says every odd closed walk in L_n has length at least ceil(n/2)=k. A repetition of an intermediate tournament would split the constructed odd closed walk into two shorter closed walks, one odd, contradicting the bound. Hence this is a simple k-cycle, proving the theorem.

## 4. Independent exact verification

Run `python3 check.py`. No third-party dependencies or repository imports are used. For every r=1,...,50, hence every n=5,9,...,201, the verifier:

- Checks the complete Walecki edge partition.
- Deletes the named vertices, adds repairs, and traverses each actual path, checking all vertices occur once.
- Checks the exact edge multiset is K_n+2M.
- Builds the functional graph, finds all its cycles, and performs the spare-order repairs.
- Independently solves all F2 equations by graph propagation, checking every equation.
- Infers the initial tournament from every occurrence, checking consistency.
- Executes every directed Hamilton-path/nonpath-complement move and checks return to the initial tournament with all intermediate tournaments distinct.

All 50 cases pass; see `check.out`. This includes functional components with 2-cycles and multiple independently repaired cycles. The all-order assertion rests on the proof above, not this finite range. The previous frozen packet has not been modified.
