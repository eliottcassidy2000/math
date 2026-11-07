# Complete universal-source search through order ten

Status: FINITE-EXACT, independently audited using the approved rooted-reciprocity theorem. The unrestricted Hamiltonian-path spectrum and F(8,137) remain OPEN. This does not extend scalar spectrum coverage.

## Result and decisive scope

There is no pair of physical rooted tournaments, BOTH of order at most10, with all eight nonstart coordinates (h,e,L,R,d,C,b,T) identical and starts s1=s2+2 satisfying the width8 residue labels. Equivalently, no universal all-O/Z/J paired seed in the Session203 criterion exists with both source orders at most10.

The complete filtered counts by source order are: zero through order7,41 at order8,301 at order9, and1412 at order10. There are1549 distinct full-eight-plus-label keys and no key containing both residues. These are rooted-cut presentations over canonical unrooted deletion cores, not a claimed enumeration of rooted isomorphism classes.

The reciprocity dependency shows every possible pair must have opposite order parity; starts must belong to (4,2),(10,8),(16,14),(22,20); all eight shared coordinates must be divisible by3; L must equal6 modulo12; C,b must equal0 modulo12; T must equal0 modulo6. This justifies every congruence filter below. The copied theorem in dependencies retains its original candidate label as a frozen input; its separate independent audit has discharged that dependency.

## Complete canonical augmentation, without an external generator

`generate.cpp` starts from the singleton tournament. At step n it takes every retained order-(n-1) representative and every one of its2^(n-1) orientations to a new vertex. It canonically labels the resulting tournament and keeps one representative of each code.

The canonical label uses an ordered partition, initially one cell. Refinement splits each cell by the vector of outneighbor counts to the current ordered cells, lexicographically ordering the resulting cells. When refinement is not discrete, the algorithm selects the first smallest nonsingleton cell and branches over EVERY possible distinguished vertex in that cell. It inserts that vertex as a singleton followed by the remainder, refines again, and continues. At discrete leaves it takes the minimum adjacency code. There is no heuristic pruning.

Isomorphisms preserve every refinement split and induce a bijection of the exhaustive individualization trees, so isomorphic tournaments have the same minimum leaf code. Conversely every leaf code encodes an actual relabeling of the input, so equal codes imply isomorphism. Completeness follows by induction: deleting any vertex from an arbitrary order-n tournament gives a previous representative up to relabeling, and its incidence to that deleted vertex occurs among the enumerated extensions.

The resulting class counts are1,1,2,4,12,56,456,6880,191536 at orders1,...,9. The final stage processes1,761,280 extensions in9.34 seconds in this run. The counts are controls, not substitutes for the exhaustive augmentation argument.

`verify_generation.cpp` compares the canonical equivalence relation against brute-force labeling on all1099 labeled tournaments through order5, checks equality of the fresh order8 class set with the6880 inherited classes, checks canonical idempotence on every one of191536 order9 representatives, and checks2000 deterministic random relabelings. Every tournament row is validated for the tournament axioms while loading.

## Complete root-cut enumeration and safe early filters

Every rooted tournament U of order at most10 is obtained by choosing its deletion Q=U-root from the complete bank at order at most9 and choosing the root's outgoing incidence subset. Relabeling Q transports that cut, so no rooted realization is lost by using unrooted representatives for Q.

For each deletion core, `search.cpp` first obtains exact subset Hamilton-path counts H_A and endpoint counts E_A(i),S_A(i). It rejects H_Q not divisible by3. It also rejects a core if every positive full-core start count exceeds22: any nonzero selected start sum would then exceed the admissible budget. Next it obtains B_A(2) by exact subset convolution and rejects B_Q(2) not divisible by12. These rejection tests are necessary, not heuristic.

The linear root-incidence formulas are the proved formulas in dependencies/SUFFIX_RIGIDITY.md, Section5. They produce s,e,L,R,T for every cut. The allowed starts and the necessary e,L,R,T congruences are tested before any quadratic calculation.

For a surviving cut define
E_A=sum_(i in incoming-root-cut) E_A(i),
S_A=sum_(i in outgoing-root-cut) S_A(i).

Splitting the path containing the root before and after that root gives

h=s+e+sum_(nonempty A subset Q) E_A S_(Q\A),
C=2[L+R-H_Q+sum_(nonempty disjoint A,B) E_A S_B H_(Q\(A union B))].

Empty prefixes/suffixes are counted by the s,e,L,R terms; the displayed double sum uses both portions nonempty. The factor2 selects the labeled bin containing the root. These are exact disjoint decompositions of literal paths/covers. The search applies the remaining h,C congruences and writes EVERY surviving profile.

For order9 deletion cores, the exact run processes191536 graphs, with6365 surviving early core filters,3,258,880 linear cut checks,3722 quadratic cut checks, and1412 surviving rooted order10 profiles. Runtime was12.19 seconds. The number of all theoretically possible cuts is98,066,432; early rejection avoids explicitly forming most of them, using only necessary conditions.

No fixed count cap or floating arithmetic is used. At source order10 every relevant path-cover count is bounded by 10! times the number of weak compositions into at most3 bins, well below signed64-bit range.

## Independent controls and exact matching

The new implementation, replayed on the fresh deletion order8 bank, exactly reproduces the multiset of301 pairable prior order9 profiles from the earlier independent root-incidence tensor compiler. The older raw list also had80 start26 presentations, which cannot pair with a permitted start28 and are excluded by the sharper admissible-start list. This explains the difference from the earlier381-source comparison list; it is not an omitted valid case.

`verify_candidates.cpp` independently reconstructs and directly counts ALL1412 literal rooted order10 graphs, computing their full subset path/cover tables rather than the factored root-incidence formulas. Every nine-coordinate profile agrees.

The key is (h,e,L,R,d,C,b,T,s+2r), with r=s mod3. When the eight nonstart coordinates agree, equality of the last coordinate between residues1 and2 is exactly s1=s2+2. Retaining one representative per such key and residue cannot discard a possible cross-residue match. `combine.py` combines all complete candidate files from source orders2 through10. No key contains both residues. The singleton source is excluded directly because its start is odd.

## Broader but finite positive search

A separate search combines the new order10 modular sources with the prior modular root-cut bank, removes unpairable starts and duplicate full profiles, and retains1874 literal sources. It applies every O/Z/J word of length at most7 independently on each source. There are6,146,720 presentations and2,335,127 full-profile-plus-label keys, again no universal pair.

This last result is NOT an all-depth theorem. It permits source orders up to17 but covers only the stated construction menu. It does not exclude arbitrary order11 sources or deeper words. At order17 the elementary bound17!*binomial(19,2) on ordered three-path covers remains below signed64-bit range, so the exact search does not overflow.

## Reproduction and files

The compressed BANK9.txt.gz is the exact fresh191536-class bank; decompress it to BANK9.txt before using the supplied search and check programs. Alternatively compile generate.cpp and run it with argument9; this regenerates BANK2 through BANK9 deterministically; BANK1 is the supplied trivial singleton. Timing fields in receipts are descriptive and may differ on replay.

Compile with a C++17 compiler:
- generate.cpp and run argument9;
- verify_generation.cpp and run without arguments;
- search.cpp and run `BANK9.txt TARGETS.txt` in this directory;
- verify_candidates.cpp and run without arguments;
- deep_search.cpp and run `DEEP_SOURCES.txt 7` for the separate finite grammar probe.

The replay1 through replay8 folders retain the smaller-order candidate inventories and receipts. `combine.py` rebuilds the finite-completeness receipt. TARGETS.txt contains the earlier381-source superset for the initial order10 comparison; completeness itself uses the exact combined inventories, not an assumption about that target list. DEEP_SOURCES.json gives the literal provenance of every source in the broader grammar probe.

Prior packet files were read only. New code, banks, receipts and provenance are confined to this separate directory. No Git publishing was performed.
