# Full n=5 alphabet and resolution of all eight Session29 boundary bridges

## Result

All 324 edges in the original selected-alphabet order-10 wiring family are global nonbridges. The earlier 316 certificates are supplemented here by explicit cycles for the eight remaining edges. No enumeration of all order-10 tournaments was performed.

The full alphabet has only 64 Hamiltonian-spine tiles. Its base incidence graph has 12 ordinary source types and 33 exit-rooted complemented types. The two-block boundary has 78 source nodes, 561 target nodes, and 4096 distinct tile edges. It is a weighted symmetric-square incidence graph, with ordered off-diagonal products retained as separate parallel lines.

## Encoding and dual normalization

A five-vertex tournament is a 10-bit integer: bit k corresponds to the k-th lexicographic pair i<j, and bit 1 means i->j. Spine mask 657 fixes 0->1->2->3->4. Complementation off that spine is XOR 366. Let tau be converse followed by reversing the vertex order. A normalized pair (e,f) means the actual word (e,tau(f)). Its source address is the unordered pair of ordinary canonical types of e,f. Its target address is the unordered pair of exit-rooted canonical types of e XOR 366, f XOR 366. Thus the second source is dual-normalized, not left unchanged.

## General symmetric-square bridge theorem

Let B be any bipartite multigraph. Form its ordered-edge symmetric-square incidence multigraph: a pair of base edges i-j and k-l gives a line between unordered source pair {i,k} and target pair {j,l}. Distinct ordered products are retained.

Every off-diagonal product has its distinct reversed parallel mate and cannot be a bridge. For a diagonal product e^2, with e=i-j:

* If either base endpoint has degree 1, the corresponding doubled vertex has degree 1, so e^2 is a leaf bridge.
* Otherwise choose other incident edges f at i and g at j. If either is parallel to e, another diagonal product already gives a parallel line. Otherwise write f=i-l and g=k-j with l!=j and k!=i. The alternate path is {i,i} --(e,f)--> {j,l} --(f,g)--> {i,k} --(e,g)--> {j,j}. This is a genuine four-cycle after adding e^2.

Consequently every bridge of this boundary graph is a leaf bridge, and the criterion is exactly that one endpoint of the base edge has degree 1. This concerns the raw multigraph, not its underlying simple graph.

The 4096-edge boundary has exactly 19 such bridges. Four are genuine global source-leaf bridges, at diagonal tiles 767,991,1021,1023. The other 15 are only boundary target leaves. This work resolves two of those globally (689,721); it does not classify the other 13 globally.

## Six direct boundary cycles

Each row below gives original actual tile word, then three actual words forming an alternate path between its endpoints. These are literal tile codes, not canonical-type abbreviations.

* Original (723,917): (693,753), (691,945), (915,915).
* Original (753,693): (693,753), (691,945), (945,691).
* Original (727,981): (661,721), (661,945), (669,981).
* Original (761,701): (661,721), (661,945), (761,699).
* Original (979,919): (661,721), (661,977), (923,919).
* Original (981,727): (661,721), (661,977), (765,727).

## Two cross-boundary 4+6 repairs

The last two original blue targets have rooted addresses {842,842} and {335,335}, and Hamiltonian-path counts 585 and 221. Their rooted columns have degree 1 in the entire 64-letter alphabet. A balanced 5+5 repair therefore cannot exist. Their intrinsic cuts include sizes 4,5,6 with the same unique outgoing arc, which supplies the necessary movable cut.

Both admit the following four-cycle in the full metagraph:

BB -- original blue edge --> T -- path flip --> U -- path flip --> V -- existing boundary edge --> BB.

Here U is the SAME ordinary source isoclass in both cases, the ordered join of the strong four-vertex class with canonical mask 4 and the strong six-vertex class with canonical mask 81. This is not BB, whose SCC sizes are 5,5. V has mixed rooted address {335,842} and is already reached from BB by actual tile word (689,661). Its address differs from either original target, so all four cycle vertices are distinct.

Literal certificates, with zero-based vertices and lexicographic 45-bit tournament masks:

1. Original word (689,689), H(T)=585:
   * T mask: 34051591637519.
   * Hamiltonian path P in T: [1,2,3,4,5,6,7,8,9,0].
   * Flip off P to U mask 23708219209456.
   * Hamiltonian path Q in U: [1,2,3,4,7,5,6,8,9,0].
   * Flip off Q to V mask 30137802689039.
   * Converse V, then old-to-canonical permutation [8,5,7,6,9,1,2,4,3,0], yields canonical mixed target {335,842}.

2. Original word (721,661), H(T)=221:
   * T mask: 35013664310799.
   * Hamiltonian path P in T: [1,2,3,4,5,6,7,8,9,0].
   * Flip off P to U mask 22746146536176.
   * Hamiltonian path Q in U: [1,2,3,4,6,7,5,8,9,0].
   * Flip off Q to V mask 30651051279887.
   * Old-to-canonical permutation [0,3,2,1,4,5,9,7,6,8], without converse, yields the same canonical mixed target {335,842}.

The canonical mixed target places the first rooted block on 0..4 with exit root 4; the converse of the second rooted block on 5..9 with entry root 9; all cross arcs point second-to-first except 4->9. This convention makes the isomorphism certificates fully specified.

Only the 585 and 221 Hamiltonian paths of these two targets were inspected to find alternative 4+6 sources. Each yielded just two distinct eligible source classes. A successful return to the known boundary was found after inspecting 14 and 13 neighbor paths, respectively. These counts refer to the deterministic search order, not a claimed minimal search requirement.

## Verification and scope

check.py constructs the full64 base graph and its4096-edge boundary; detours.py finds literal 4+6 certificates. Both read the frozen sibling session29-target-wiring JSON files without modifying them. verify.py independently checks every literal Hamiltonian path, both nonpath-flip operations, canonical relabelings, source SCC component masks, six boundary path incidences, and the general bridge criterion against graph-search results. All checks use explicit exceptions rather than assertions and pass identically in normal and python -O modes.

Files in the frozen sibling directory were not modified. No repository or Git operation was performed. The decisive structural lesson is that a target leaf in a fixed block-size boundary is not an obstruction in the full metagraph: an intrinsic neighboring cut can connect it through a different SCC-size stratum.
