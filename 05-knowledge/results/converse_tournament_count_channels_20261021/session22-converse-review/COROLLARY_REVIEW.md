# Tournament and disagreement-graph corollary

Verdict: PASS, including the stated limitation on self-converse correspondence.

For an ordered a-tuple of tournaments on a common vertex set, retain T1 and let Gj, for each j=2,...,a, contain exactly those unordered pairs whose orientations differ in Tj and T1. Each Gj is an arbitrary simple graph. The inverse orients Tj by flipping T1 on precisely the edges of Gj. These constructions are mutually inverse, including a=1 and the empty vertex set.

The construction is equivariant for simultaneous vertex relabeling: two orientations agree before relabeling if and only if they agree afterward. Equivalently, in edge-bit coordinates a relabeling applies the same endpoint-order correction to both tournament bits, which cancels in their XOR. Vertex colors are unchanged by the construction, so their ordinary relabeling action is also respected. Consequently this bijection descends to isomorphism classes and proves that raw F_n(a,b) counts a tournament decorated by a-1 indexed simple graphs and one of 2^b colors at each vertex. The graphs are indexed, not unordered.

Complementing one Gj toggles every disagreement bit and hence reverses Tj while leaving T1 and every other tournament fixed. Reversing all tournaments leaves every Gj unchanged and reverses T1. More generally reversing T1 alone reverses the retained tournament and complements every Gj; this follows from the same XOR identity.

For a=2,b=0 the already-proved numerical identity C_(2n)=F_n(2,0) therefore equates the number of unmarked self-converse tournament classes on 2n vertices with the number of tournament-plus-simple-graph classes on n vertices. The equivariant construction above is ONLY a bijection between same-order tournament tuples and graph-decorated tournaments. It does not construct or establish a canonical bijection with the unmarked self-converse classes on twice as many vertices. No stronger structural identification is claimed.

Executable audit: corollary_audit.py exhausts all tournament tuples of sizes a=1,2,3 on n=0,1,2,3, and every vertex permutation. It directly checks the inverse, permutation equivariance, individual graph-complement operation, and simultaneous reversal. Both normal and python -O runs pass 5,860 explicit checks; receipts corollary-normal.json and corollary-optimized.json are identical. Colors require no extra census because the map is exactly the identity on the color coordinate. No source edits, new branch, or remote mutations were made.
