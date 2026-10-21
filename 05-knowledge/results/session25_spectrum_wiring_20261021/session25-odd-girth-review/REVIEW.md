# Independent audit: all-order labelled odd girth

Verdict: PASS, with no mathematical correction required. This audit concerns the labelled tournament graph only, never its isomorphism or converse-merged quotient.

Reviewed frozen producer files:
- RESULT.md SHA256 1c3796377b65b2d91fc6ab885c338ed4b5259fb7d54a11ac8e9bd55e217c7239
- check.py SHA256 ebc819c9f7bf1a98851b4aabb49e416451e979b33b5908dc5b435c913f93c9d9
- check.out SHA256 cee7d3e4e55bcc184df467b98f220e4e1568168c8cecf372be30b3ba451be83a

## Symbolic checks, valid for every r >= 1

Write k=2r+1, m=2k. In the Walecki word, the k finite odd-sum edges are the pairs (i+h-1,i-h), h=1,...,k; the k-1 even-sum edges are (i-h,i+h), h=1,...,k-1. Their sums are 2i-1 and 2i, respectively, modulo m. The word visits every finite label exactly once, so these are distinct and exhaust the non-loop pairs having those sums. Infinity is adjacent to i and i-k. Thus the cycles partition the complete graph. Equivalently the unique source of a finite edge {u,v}, with ordinary representatives, is floor((u+v+1)/2) modulo k; an infinity edge has source u modulo k. The independent checker starts with this edge-source formula, rather than generating the producer's cycles first.

The two deleted labels 0 and k-1 have an even sum 2r, so they are adjacent only in C_r. Deleting them there leaves a Hamilton path. In every other cycle their deletion produces exactly two nonempty path components; this includes the possible one-vertex component. Joining the former neighbors of either deleted vertex reconnects those components.

For r+1 <= j <= 2r, the neighbors of 0 in C_j are 2j-1 and 2j. The neighbors of k-1 in C_(j-r-1) are the same pair modulo 2k. Neither computation encounters a formal self-pair or infinity-neighbor exception: both resulting labels lie between 2r+1 and 4r. Those two index ranges partition all indices except r. Hence exactly the required non-base paths receive one repair. The repair pairs are mutually vertex-disjoint, avoid both deleted labels, and their source is 2j modulo k. This source differs from both repaired indices, since equality would force j=0 or j=r modulo k. Therefore each repair edge appears exactly three times and every other surviving edge once. There are no within-path duplicate edges.

## Ordering lemma and its obstruction

After contracting the disjoint pair edges, the constraint multigraph is the stated functional graph: each pair-node has one outgoing edge, the base has none, and there are no self-loops. Every connected component is a tree ending at the base or a unicyclic component. In a cyclic component, each pair has exactly one active element selected by the predecessor's c-value; its outgoing c-value is the next pair's active element. Choosing the internal pair constraint and the active-to-c constraint makes the unique expanded graph cycle consist exactly of active-to-active edges.

There is therefore precisely one independent F2 obstruction per cyclic component: the XOR of its edge constraints. This statement explicitly uses multigraph cycles. A functional directed 2-cycle gives two parallel active-to-active equations; consistency requires their right-hand sides to agree, exactly the same XOR criterion. It must not be simplified to one undirected edge.

For selected endpoints in a triple, their rank parity XOR is zero if the spare lies between them and one otherwise. Initially all spares precede all active elements. Moving one spare between its own two active neighbors toggles that summand. The spare cannot occur in a different cycle triple: pair membership is unique, and occurrence as the c-value of any cycle node would make it active. It may occur in feeding-tree constraints, which impose no compatibility condition. Thus the operation toggles exactly one cycle obstruction, even in the parallel-edge case and with many cyclic components. Subsequent moves preserve all active relative orders and leave unmoved spares before them. All obstructions can be independently removed, and spanning-forest propagation then solves the complete system. The third equation in each triple is the sum of the other two.

This is an all-order proof, not an inference from finite testing. Arbitrary incidence labels are allowed, so no unstated Walecki traversal-orientation identity is required.

## Actual directed transition, not only an undirected cover

At position p, an edge at occurrence rank q was absent from precisely p-q earlier paths and was toggled that many times. Reversal by x_v+p changes its displayed orientation to d+x_v+p. The necessary initial direction is consequently d+x_v+q over F2. The lemma makes this independent of occurrence for every tripled edge. Single occurrences impose no additional compatibility. All displayed paths are therefore genuinely directed at their prescribed times. Each edge is retained an odd number of times in an odd-length sequence, so it is toggled an even number of times; the final tournament equals the initial one.

The prior odd-closed-walk lower bound is ceil(n/2)=k. Any intermediate repetition in this length-k odd closed walk splits it into two positive shorter closed walks, one odd, contradicting that bound. Thus the construction yields a simple cycle and proves the exact odd girth. Combined with the existing n=2 modulo 4 family and n=0,3 modulo 4 parity obstruction, this establishes the claimed labelled classification for n>=3.

## Separate executable validation

The independent audit.py has no imports from producer code or the repository. It builds cycles by the edge-source formula, checks the Walecki words against them, reconstructs repaired paths, verifies exact edge multiplicities and matching structure, uses different choices of path traversal, active ordering, and cycle-spare selection, then checks all three pairwise constraints of every triple. During each spare move it checks that exactly the selected obstruction toggles. Actual tournaments and transitions are checked using integer bitsets, including every directed path arc and uniqueness of all intermediate tournament states.

All r=1,...,60 plus r=75,100 pass, corresponding to every n=5,9,...,241 plus n=301,401. See audit.out. These include 2-cycles, longer functional cycles, multiple cyclic components, and multiple independently required spare flips. Finite validation supplements, rather than replaces, the proof above.
