# Independent path-complement feedback audit

Root's direct theorem is correct: for every Hamilton path P of an n-vertex tournament T, and T'=kappa_P(T), one has tau(T)+tau(T')>=floor(n/3), where tau is the minimum feedback-vertex deletion number. Partition P into disjoint consecutive triples. Exactly one of T,T' is cyclic on each triple, since the two path edges stay and the chord reverses. A feedback deletion set on that side must meet the triple. Disjoint triples imply the sum bound.

For transitive T and its unique Hamilton path, T' has a directed triangle exactly on three consecutive positions: triples with zero or one adjacent gap are transitive. Deleting every third vertex kills all cyclic triples. Thus tau(T')=floor(n/3), establishing sharpness. In particular the transitive node's unique neighbor has exactly this feedback number.

Both endpoint feedback numbers are invariant under isomorphism and converse. The inequality therefore descends to either edge color and all merged loops. Nodes with tau<=d form an independent set when2d<floor(n/3), in particular the one-apex stratum for n>=9. A loop needs tau>=ceil(floor(n/3)/2).

A stronger path-specific certificate colors every consecutive triple by its cyclic owner T orT'. Let nu_c be the maximum number of vertex-disjoint intervals of colorc. Then tau on each side is at least its nu_c. These are ordinary interval packing numbers, computable greedily. This can improve the sum bound (n8 alternating colors give2+2 versus2). It is not a characterization of feedback number: nonconsecutive triangles are discarded.

Independent exact controls: `node verify_wiring.js > wiring-output.json`. All labelled tournaments at orders3,4,5; feedback number by exhaustive triangle-free induced subsets; every directed Hamilton path; independent path-complement operation. Output counts12,192,7680 directed instances; all satisfy the theorem. This is an audit of root's proof, not a claim of discovery or a larger census.
