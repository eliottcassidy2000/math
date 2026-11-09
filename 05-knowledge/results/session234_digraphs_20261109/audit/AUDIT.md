# Independent audit: semicomplete Hamiltonian-path counts 7 and 21

## Verdict

The conjecture is false under the usual Hamiltonian-path convention. The corrected finite-core classification passes the independent audit, provided its count-three factor includes **both** strong order-three cores. The original assertion that this factor must be a directed triangle was false: a directed triangle with one reverse arc also has exactly three Hamiltonian paths. The correction was independently reproduced during the audit.

There are three strong count-seven cores (order four), and thirty-two strong count-21 cores (sixteen each at orders five and six). All other semicomplete examples are their ordered SCC compositions described below. Diagonal bits have no effect.

## Independent computation

`audit.cpp` uses a deliberately different algorithm from the scanner's subset dynamic program. For each vertex permutation it enumerates every compatible semicomplete adjacency matrix, incrementing that matrix's counter exactly once. Each of the n-1 path pairs permits the required arc alone or both arcs; every other pair permits all three states. Thus the full count array is computed without using a Hamiltonian-path recurrence or the scanner's counts. The tournament version permits only the required orientation on path pairs and either orientation elsewhere.

A uint16 counter is safe at these orders: the absolute bounds are 6!=720 and 7!=5040. All target isomorphism classes are canonicalized by trying every vertex permutation; strongness is independently tested by reachability from every vertex.

Complete labelled semicomplete censuses at orders 1 through 6 have totals 1,3,27,729,59049,14348907. Count-seven frequencies are 0,0,0,72,720,6480. Count-21 frequencies are 0,0,0,0,1760,32640. Every histogram entry, not merely the target frequencies, agrees with the separate scanner. Every order-six target canonical representative, strongness flag, and reciprocal-pair count agrees exactly.

Independent tournament scans inspect **every labelled tournament** through order seven. Strong tournament numbers at orders 3 through 7 are 2,24,544,22320,1677488, and their minimum path counts are 3,5,9,15,25. No tournament in these exhaustive scans has count 7 or 21.

Controls include total graphs = b^(n choose 2), and total path incidences = n!*2^(n-1)*3^((n choose 2)-n+1) for semicomplete graphs, or n!*2^((n choose 2)-n+1) for tournaments. `verify.py` checks these and comparisons with scanner files. All controls pass.

## Structural proof audit

1. **SCC ordering and multiplicativity.** Between distinct strong components every arc goes in one direction: opposing cross-component arcs would create a cycle in the condensation. Semicompleteness makes the condensation a transitive tournament. Its order is unique. A directed path cannot leave a component and return, so every Hamiltonian path visits the components as consecutive blocks in this order. All choices of within-component paths concatenate. Consequently H(D)=product H(C).

2. **Strong semicomplete Hamiltonicity.** Strongness provides a directed cycle, allowing length two. Choose a longest cycle C. Any outside vertex not insertable between two consecutive cycle vertices must strictly dominate all of C or be strictly dominated by all of C. Indeed if c_i->v then noninsertability forbids v->c_(i+1), and semicompleteness propagates c_(i+1)->v around the cycle. If no c_i->v exists, every v->c_i exists. If outside vertices remain, strongness forces a path from the dominated class to C. It must first reach the dominating class, giving an arc b->a between the two classes. Replacing a cycle edge by c_i->b->a->c_(i+1) makes a longer cycle. Contradiction.

3. **Spanning strong tournament.** For n>=3, the spanning directed cycle uses distinct unordered pairs. Retain those arcs and select one orientation from each remaining pair. The resulting spanning tournament is strong, so its path count is at most H(D). This argument explicitly does not work at n=2, where a two-cycle needs both arcs of its only pair.

4. **All-order bound.** The classical Moon strong-subtournament theorem supplies a strong induced (n-1)-vertex tournament in any strong n-vertex tournament, n>=4. Every path in T-v can be extended by inserting v; fix one feasible position for each path. Deleting v recovers the path, so this extension map is injective. Hence the strong tournament minimum cannot decrease beyond order three. The independently verified minima f(5)=9 and f(7)=25 exclude strong semicomplete count-seven cores of order >=5 and count-21 cores of order >=7. This audit relies on the cited classical Moon theorem; it does not claim to independently prove that classical theorem. The cited repo source is THM-1370-h-spectrum-omits-7-21-all-n.md, citing Moon, Canadian Mathematical Bulletin 9(3) (1966), 297–301, Theorem 2, DOI 10.4153/CMB-1966-038-7.

5. **Small factors.** Strong order-two semicomplete graphs have H=2 and cannot be factors of 7 or 21. Every strong component of order >=3 has H>=3. H=3 requires order three because f(4)=5. Exactly two strong order-three types have H=3: rows 010/001/100, and rows 011/001/100. The latter's added arc 0->2 occurs in no Hamiltonian path; the three paths in each are 012,120,201. Larger strong components cannot have H=1. Thus all H=1 components are singletons.

## Complete corrected generators

For H=7: choose one of the three strong order-four cores and put any number of singleton SCCs before and after it.

For H=21: either choose one of the thirty-two strong count-21 cores and put singleton SCCs before/after; or choose one of the three H=7 cores and **one of the two H=3 cores**, choose their order, and put arbitrary singleton chains before, between, and after them. Relabel arbitrarily. There are no other factorizations, because 21=3*7 and all nontrivial strong factors have count at least two. Uniqueness is uniqueness of the ordered SCC decomposition, up to the isomorphism type of each core.

At n>=4 the H=7 isomorphism count is 3(n-3), and labelled count is 3(n-3)n!.

For H=21, sum terms only when their minimum order is met:
- order-five single cores, n>=5: 16(n-4) classes; (44/3)(n-4)n! labelled graphs;
- order-six single cores, n>=6: 16(n-5) classes; 16(n-5)n! labelled graphs;
- mixed order-four/order-three cores, n>=7: 6(n-5)(n-6) classes; 4(n-5)(n-6)n! labelled graphs.

The mixed term uses six choices of core types, two orders and binomial(n-5,2) distributions of singleton vertices. Its orbit weight uses automorphism orders 3 and 1 for the two order-three cores. The order-five orbit weight uses fourteen asymmetric cores and two with automorphism order three. The scanner provides these automorphism details; the independent class counts and labelled totals are consistent with them.

## Edge cases and interpretation

- A Hamiltonian path is an ordered sequence visiting each vertex once, with an arc between consecutive vertices. Reversal is a different path when valid. It is not a Hamiltonian walk, cycle, or diagonal-dependent convention.
- Self-loops are never used. Every result persists under any of the 2^n independent diagonal choices; labelled matrix counts with arbitrary diagonals are multiplied by 2^n.
- At n=1, H=1 regardless of its diagonal. At n=2, the two one-way graphs have H=1 and the reciprocal graph has H=2.
- At n=0 the conventional count is one empty path. Defining it as zero instead changes only this empty case and none of the assertions for 7 or 21.
- A source/sink singleton does not change H; repeated additions form a transitive prefix/suffix. For two nontrivial components, singleton SCCs may also occur between them. Inserting extra arcs *within* or backward between prescribed components is not such padding and need not preserve H.
- The results classify semicomplete digraphs, i.e. no absent unordered pair. They do not classify all unrestricted binary matrices with missing pairs. They completely classify the loopless counterexamples to the user's claimed implication.
- Sharp replacement implications are H=7 => at least one absent unordered pair OR at least two reciprocal pairs; H=21 => at least one absent unordered pair OR at least one reciprocal pair. No diagonal condition belongs in them.

## Reproduction

Run `g++ -O3 -std=c++17 audit.cpp -o audit`. For n=1,...,6 run `./audit 3 n > semicomplete-n.txt`; for n=3,...,7 run `./audit 2 n > tournament-n.txt`. Run `python verify.py`. Plain-text outputs and complete histograms accompany this report. All accompanying outputs are reproducible from the supplied sources.

## Supplement: independently checked loop-decorated classes

`verify.py` independently computes each strong representative's full automorphism group, then enumerates all vertex-loop subsets modulo that group. The numbers of loop-decorated strong core types are 12 for H=3 at n=3; 48 for H=7 at n=4; 480 for H=21 at n=5; and 1024 for H=21 at n=6. These are exact orbit enumerations, not division by automorphism-group order. In particular the directed triangle has four loop orbits, the doubled triangle has eight, and each order-five core with automorphism order three has sixteen loop orbits.

The SCC ordering makes every singleton chain an ordered binary word (loop absent/present), so its ordinary generating function is 1/(1-2x). Thus the ordinary generating functions counting loop-decorated isomorphism classes are:

H=7: 48x^4/(1-2x)^2.

H=21: (480x^5+1024x^6)/(1-2x)^2 + 1152x^7/(1-2x)^3.

The mixed coefficient is 2*48*12=1152. For comparison, the loopless ordinary generating function for H=21 is (16x^5+16x^6)/(1-x)^2 + 12x^7/(1-x)^3, and the loopless labelled exponential generating function is ((44/3)x^5+16x^6)/(1-x)^2 + 8x^7/(1-x)^3. All match the corrected formulas. For labelled matrices with arbitrary diagonal bits, substitute 2x into the corresponding loopless labelled exponential generating function, equivalently multiply each order-n count by 2^n.

## General fixed-count finite-core theorem: independent audit

**Verdict:** The proposed theorem is valid. One qualification is necessary: the stated *exact degree* is meaningful only when at least one core word with product h is realizable. If there is no such word, the generating functions are identically zero. No assumption that every positive h is realizable is needed.

Fix h>=2 and work with loopless semicomplete digraphs first. For any nontrivial strong component C of order m>=2, the spanning-cycle argument above supplies a directed Hamiltonian cycle. Its m cyclic rotations give m distinct ordered Hamiltonian paths, even at m=2. Hence H(C)>=m>=2. This avoids the stronger classical tournament-floor theorem entirely.

By SCC multiplicativity, every nontrivial component has path count d dividing h and order m<=d<=h. Thus there are finitely many possible component isomorphism types, obtainable by finite enumeration through order h. Let W_h consist of all finite ordered words (C_1,...,C_k) of those nontrivial strong core types with product H(C_i)=h. Every d_i=H(C_i)>=2; counting prime factors with multiplicity gives k<=Omega(h). The alphabet and permitted word lengths are finite, so W_h is finite.

For w in W_h write s(w)=sum |C_i| and a(w)=product |Aut(C_i)|. The exact loopless ordinary generating function for isomorphism types and exponential generating function for labelled digraphs, respectively, are

U_h(x) = sum over w in W_h of x^s(w)/(1-x)^(k(w)+1),

L_h(x) = sum over w in W_h of x^s(w)/(a(w)(1-x)^(k(w)+1)).

The ordered SCC decomposition is unique, so distinct words do not overlap. Its k+1 gaps contain arbitrary ordered singleton chains. An order-preserving automorphism of the SCC chain fixes each component position; hence there are no extra permutations of equal cores or factorial denominators for repeated core types. The automorphism group of the entire graph is the product of those of its cores, justifying the labelled weight.

For integers d_i>=2 one has sum d_i<=product d_i when k>=2; this follows from ab>=a+b for a,b>=2 and induction. At k=1 equality holds. Consequently s(w)<=sum d_i<=h. For all n>=h the coefficient formulas are therefore actual polynomials:

[x^n]U_h(x) = sum_w binomial(n-s(w)+k(w),k(w)),

[x^n]L_h(x) = sum_w binomial(n-s(w)+k(w),k(w))/a(w).

The labelled graph count is n! times the latter polynomial. If W_h is nonempty, both polynomials have exact degree K=max_w k(w)<=Omega(h), since all leading coefficients are positive. Thus there is no cancellation. This is an eventual polynomial theorem for isomorphism counts and for labelled counts *divided by n!*; raw labelled counts are not themselves eventually polynomial.

For arbitrary diagonal bits, replace each strong core type by its finitely many loop decorations modulo its automorphism group. Every singleton gap is now a binary word, so replace the denominator (1-x)^(k+1) by (1-2x)^(k+1). The coefficients of a term with total core order s are 2^(n-s) binomial(n-s+k,k) once n>=s. Thus, for n>=h, isomorphism counts are 2^n times a polynomial of degree K (when realizable). The raw labelled counts are n!*2^n times their loopless polynomial, because every labelled off-diagonal matrix admits exactly 2^n diagonal choices. Equivalently the labelled EGF is L_h(2x). Distinguish these normalized and raw labelled statements when presenting the theorem.

For h=1, there can be no nontrivial strong component, since every such component has at least two paths. The graph is a transitive total order. Loopless isomorphism counts are one at every n, and labelled counts are n!. With arbitrary diagonals the isomorphism counts are 2^n and labelled counts are 2^n n!, because a transitive tournament has no nontrivial automorphisms. Corresponding ordinary and normalized labelled generating functions are 1/(1-x) and 1/(1-2x).

For h=0, there are no semicomplete graphs: each singleton component has a path, each nontrivial strong component has a spanning cycle and hence a path, and these concatenate. This includes order zero under the conventional one-empty-path convention; adopting H(empty)=0 introduces that sole exceptional object.

No new external theorem is needed beyond the self-contained spanning-cycle proof already audited, and no novelty claim is made.

## Unrestricted binary digraphs: exact generator and minimum orders

The proposed union/hypergraph generator is correct, with the optional-edge ground set made explicit. Fix the vertex set V and a set P of exactly h distinct vertex permutations. Let U be the union of their consecutive directed edges. Every permutation in P is automatically a Hamiltonian path of U. Retain this P only if H(U)=h, equivalently the Hamiltonian-path set of U is exactly P. Let F consist of all off-diagonal arcs not in U. For every permutation Q outside P, put B_Q=E(Q) minus U, a subset of F. These hyperedges are nonempty because H(U)=h. Choose any X subset F such that B_Q is not a subset of X for every Q outside P. Then D has arc set U union X and Hamiltonian-path set exactly P: the old paths persist, and a new permutation Q is present precisely when all its missing edges B_Q were added. These X are the independent vertex subsets of the hypergraph on the optional arcs F with hyperedges B_Q (duplicate or nonminimal hyperedges may be removed without changing the condition).

Conversely, for any digraph D with exactly h Hamiltonian paths, choose its actual path set P. Its union U is a subgraph of D and has exactly P as its path set: inclusion gives both containments. The unique remainder X=E(D) minus U satisfies the hyperedge condition. Thus every off-diagonal adjacency matrix appears exactly once, indexed by its actual P and X. Diagonal bits can be chosen independently afterward. This is an exact finite generator, not an efficiency claim. Without the condition H(U)=h, an empty forbidden hyperedge would correctly make that proposed P infeasible.

For any unrestricted order-four digraph, either all twelve off-diagonal arcs are present, giving H=24, or at least one directed arc a->b is missing. Exactly 3!=6 of the 4!=24 permutations contain a immediately followed by b, and all six are excluded. Hence H<=18 in the latter case, ruling out H=21 at order four. At orders at most three the trivial bound H<=3!=6 suffices. Therefore the displayed order-five H=21 example is globally minimum even among unrestricted binary matrices. Likewise the order-four H=7 example is globally minimum, because order-three digraphs have at most six paths. Loops do not affect either argument.

The independent permutation-scatter census was extended to all unrestricted loopless digraphs at orders 1 through 5 in `audit_all.cpp`. Its fourth pair state is absent (state 3), while state 2 is reciprocal; this differs from the scanner's encoding but enumerates the same matrices. Counts at order four are H7=120, H21=0; at order five H7=49440, H21=7760. Every target stratum by order, H, absent pairs, reciprocal pairs and strongness agrees exactly with the scanner's separate DP/direct enumeration. Full graph-count and path-incidence controls pass at every order; complete independent histograms are saved as `hist-base4-n*.csv`. The scanner did not initially export complete unrestricted histograms, so only target-stratum equality, rather than equality of those entire histograms, was checked. `verify_all.py` and `verification-all.txt` reproduce the comparisons.
