# Independent audit: six- and seven-pair fixed-half-path maxima

## Verdict and scope

The finite claims F_6=1185 and F_7=8295 pass the audit. Here F_m is the maximum number of fixed Hamiltonian paths for a specified normalized anti-involution in the signed-controller model. It is not the ordinary tournament Hamiltonian-path maximum, the older undefined M_k, or the maximum number of blue tiles.

The audit is intentionally asymmetric:

* Six pairs: a different signed-adjacency/pair-subset dynamic program recomputed every class optimum, over all 34,816 normalized label assignments. All 56 class maxima agree, not only the global maximum.
* Seven pairs: complete producer-source replay gives byte-for-byte equality with all 456 published class rows. Independently written orbit coverage checks cover every controller mask without overlap. The different adjacency DP verifies the emitted optimal witness of every class, but does not recompute its upper bound over all its label assignments. Thus the full seven-pair upper-bound calculation is source-audited and replayed, not independently reimplemented end-to-end.
* The two literal global witnesses pass tournament completeness, anti-involution, fixed-path-count, and centralizer tests.

No proof of F_m=C_m for all m, no universal bound F_m<=C_m, and no asymptotic for F_m/m! follows. The report correctly distinguishes these conjectures from the proved formula and asymptotic for C_m alone.

## Coverage proof

Represent controller arcs on unordered pairs i<j by q_e=0 for i->j and q_e=1 for j->i. Under a permutation p, its transported bit at {p(i),p(j)} is q_e XOR [p(i)>p(j)]. This is exactly the producer transformation. Scanning masks in increasing order, processing an unseen mask, and marking every permutation image therefore selects exactly one representative per ordinary isomorphism class: an orbit is fully marked upon its first encounter and later orbits cannot intersect it.

The independent implementation instead explicitly extracts each directed source and destination, transports both, and re-encodes whether transported source exceeds destination. It checks that every marked mask belongs to just one supplied class, and that all 2^(m choose 2) masks are covered. Results:

* m=6: 56 classes, 32,768 masks, no overlap.
* m=7: 456 classes, 2,097,152 masks, no overlap.

These counts do not rely on assuming the known sequence of unlabelled tournament numbers.

For each nonsink controller, simultaneous flipping of all outgoing labels is a free action. Its first outgoing label can uniquely be normalized to +1. Stars partition the controller edges, so all such normalizations can be made independently. Enumerating all remaining labels visits every tuple of outgoing bipartitions exactly once, even if controller automorphisms create additional redundant assignments.

A useful independent checksum is the total number of assignments. A tournament has at most one sink. The number of m-vertex classes with a sink equals a_(m-1), by deleting or adjoining its unique sink. Thus the exact total over ordinary classes is

N_m = 2^(binom(m,2)-m) (a_m + a_(m-1)).

For m=6, this is 512(56+12)=34,816. For m=7, it is 16,384(456+56)=8,388,608. The implementation's direct degree-based totals agree.

## Phase proof and independent optimization

Each half-path permutation contributes the product of 1-b_vu b_vw over internal vertices v controlling both neighbors, followed by a factor 1+b_vu if its last vertex controls its predecessor. The wedge product is invariant under each outgoing-star flip. Consequently F=A+sum_v B_v, and flipping star v changes only B_v to -B_v. Since these flips are independent, max over phases equals A+sum_v |B_v|. Zero coefficients merely create ties. The producer prints coefficients before optimization and the negative-label mask after optimization; this convention is consistent.

The independent algorithm never evaluates these wedge products. It constructs the signed inter-pair adjacency: if c controls w with label b, then c_s -> w_t exactly when s=-b, independent of t. It counts directed half paths using D[S][v], where S is the set of visited pairs and v is the current signed endpoint. Initialization permits either sign in any pair. Transitions use literal directed arcs to unused pairs. Summing states on all pairs with positive endpoint enforces the normalized fold. Previously chosen signs need not be retained because those vertices are never revisited.

For a baseline normalized assignment it evaluates f, then evaluates f_v after each single star flip. Since B_v=(f-f_v)/2, the exact phase optimum is f+sum_v max(0,f_v-f). It also builds the maximizing combined flip and directly recounts its paths, guarding against a mistaken phase reconstruction. This different evaluator matches every six-pair upper bound and every supplied literal class witness at both orders.

## Literal witnesses and centralizers

The JSON witnesses are read as directed tournament outmasks without rebuilding them from the producer's label formula. They are complete tournaments, have sigma(v)=v XOR 1 as an anti-automorphism, and have positive-to-negative mate arcs. Their counts are 1185 and 8295 by the adjacency subset DP.

An automorphism commuting with sigma maps mate-pairs to mate-pairs. In these even-order normalized witnesses the directed mate edge forces preservation of plus/minus signs. Thus testing all pair permutations exhausts the centralizer; no hidden within-pair flips remain.

Both centralizers are cyclic of order five, rotating pairs 0 through 4 and fixing the tail. Their own even-order blue-tile orbit counts are therefore 1185/5=237 and 8295/5=1659, respectively, under the established free action on ordered Hamiltonian paths. These are witness values, not claims about global blue maxima. The optimal controllers are the cyclic regular five-vertex tournament followed by one or two transitively ordered terminal vertices.

## Candidate and bounded negative tests

C_m=m! sum_{j=0}^{floor(m/2)} 1/(2^j j!) has EGF exp(z^2/2)/(1-z), the partial-matching/ordered-unmatched interpretation, and C_m/m! -> sqrt(e). These algebraic facts are valid for C_m independently of tournament extrema. Agreement through seven pairs remains finite evidence. This audit does not establish a matching injection, an all-order attaining family, or the needed upper bound.

The root's separate 32 translation-invariant eight-pair label tests on regular-seven-core-plus-sink controllers are not part of this audit. Failure of that bounded search to attain 66,465 cannot refute the candidate maximum.

## Reproduction and files

Compile independent_dp.cpp with g++ -O3 -std=c++17 to a temporary executable. Run it with arguments 6 or 7 followed by the corresponding producer m6_classes.txt or m7_classes.txt path. The six-pair run performs complete independent optimization. The seven-pair run performs coverage and all supplied class-witness verification only; it explicitly prints independent_max=not_run.

Compile producer_replay.cpp to a temporary executable; run with argument 7. Compare stdout byte-for-byte with producer m7_classes.txt. Run witness_audit.py for the literal tournament and centralizer checks. witness_audit.py accepts the producer packet directory as its optional first argument; by default it looks for the sibling session29-halfpath-asymptotics directory.

Only source text, reports, JSON, and stdout text belong in the publish scope. Executables were written to /tmp and are excluded. No repository or Git edits were performed.
