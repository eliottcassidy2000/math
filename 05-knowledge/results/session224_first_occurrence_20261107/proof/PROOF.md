# First-birth fibers: audited counterexample and exact SCC recursion

Status: PROVED for the arguments below; independent small-order computation in check.py is a redundant verification. No claim that every odd value other than 7 and 21 is realized is needed or made.

## 1. The earliest failure is H=9, born at order five

Convention: H(T) counts directed Hamiltonian vertex sequences, with no quotient by automorphisms or reversal. Rows are outgoing-neighbor bitmasks, on vertices 0,...,n-1. Let

A = [14,28,24,16,1], B = [10,28,25,16,1].

A is the transitive tournament i -> j for i<j, except that 0 -> 4 is reversed. B is A with the additional reversal 0 -> 2 to 2 -> 0. Both contain the directed spanning cycle 012340, so both are strong.

A has exactly nine Hamiltonian paths. A path not using 4 -> 0 must follow the increasing order, giving 01234. A path using 4 -> 0 has the form L,4,0,R, where L and R partition {1,2,3} and are each increasing; each of the eight choices for L gives a path. This proves completeness without brute-force enumeration.

Changing A to B removes exactly the A-paths 13402 and 14023 (the paths containing 0 -> 2) and adds exactly the B-paths 12034 and 20134 (the paths containing 2 -> 0). To check completeness of the additions, a B-path using 2 -> 0 but not 4 -> 0 must otherwise be increasing, hence is 12034 or 20134. A path cannot use both 2 -> 0 and 4 -> 0, because a vertex has at most one predecessor. Thus H(B)=9 as well.

The maps i -> 4-i for A and (1 4)(2 3), fixing 0, for B reverse every arc. Hence both classes are individually self-converse. They are not isomorphic: A has exactly three directed 4-cycles, 01240, 01340, 02340; B has exactly two, 01240 and 01340. The cycle list is exhaustive by inspection of the five four-vertex subsets. Alternatively A has nontrivial modules {1,2}, {2,3}, {1,2,3}, while B has none. Both have sorted score sequence (1,1,2,3,3).

Both have trivial automorphism group. In either tournament the unique score-2 vertex is fixed. Within each of the two equal-score pairs, the internal oriented edge forces both vertices to be fixed, since an automorphism preserves scores and cannot swap the ends of a directed edge. The two score-1 vertices and the two score-3 vertices are therefore individually fixed.

No tournament of order at most four has H=9. For order at most two, H=1. At order three the transitive class has H=1 and the directed triangle has H=3. At order four, a tournament with no triangle is transitive. Otherwise take a triangle and the remaining vertex. If it is uniform to the triangle this yields a source-plus-triangle or triangle-plus-sink, each H=3. In either mixed case it yields the unique strong four-vertex class; the one-outneighbor and two-outneighbor cases are isomorphic (check their four-vertex arc sets). This class is the transitive four-vertex tournament with its source-sink edge reversed, so the same partition argument gives H=1+2^2=5. Thus the order-four spectrum is {1,3,5}. Each of 1,3,5 has exactly one class at its own first occurrence (orders 1,3,4). Therefore order five is the earliest possible failure of the asserted birth dichotomy, and H=9 supplies it.

The independent script enumerates all 2^(n choose 2) labelled tournaments for n=1,...,5 and all n! possible paths. It finds precisely 240 labelled H=9 tournaments at order five. Since A and B each have 5!/Aut=120 labellings, these exhaust the birth fiber. Consequently its full symmetry signature is N=2 ordinary classes, S=2 self-converse classes, P=0 non-self-converse pairs, M=2 converse orbits.

## 2. What the two-kind description correctly classifies

Converse is an involution on the isomorphism classes in each fixed (n,H) fiber: reverse the order of a path to obtain a path in the converse tournament. Every converse ORBIT therefore consists of one self-converse class or two mutually converse classes. A fiber can contain multiple orbits of either kind. If N is the number of ordinary classes and S its self-converse classes, the number of two-element orbits is (N-S)/2 and the total orbit count is (N+S)/2. Thus the proposed fiber dichotomy is precisely the extra assertion that this last quantity equals one. H=9 disproves that extra assertion.

## 3. Exact minimal-order recursion from strong components

Every tournament is uniquely an ordered sum of its strongly connected components. All arcs go from each earlier component to each later component. A Hamiltonian path traverses these components in their order, concatenating independently chosen Hamiltonian paths in them. Consequently H is multiplicative. Ordered component isomorphism classes determine the tournament class, and converse reverses the component word and converses each letter. These are inherited mechanisms, not new priority claims.

For q>1 define s(q) to be the minimum order of a nontrivial strong tournament with H=q, with infinity if absent. Let a(q) count its strong classes at that minimum order, and c(q) the self-converse ones there; take both zero when s(q)=infinity. Let b(q) be the minimum order of any tournament with H=q, with infinity if absent, and let N(q),S(q) count ordinary and self-converse classes in that birth fiber.

For recursion ONLY put b(1)=0 and N(1)=1, representing the empty component word. The actual positive-order birth of H=1 is the singleton at order one. All formulas below for q>1 use this auxiliary convention.

A minimal realization of q>1 has no singleton SCC: removing one would preserve H and lower its order. Each remaining component must have the least possible strong order for its own path count, since replacing it with a smaller strong realization again preserves H. Likewise any prefix or suffix component word must minimize its own order. Thus

b(q) = min over divisors d>1 of q [s(d)+b(q/d)],

N(q) = sum over divisors d>1 of q with s(d)+b(q/d)=b(q) [a(d) N(q/d)].

Only finite-cost terms are included. The N formula uses the FIRST component, so it does not overcount different parenthesizations of a product. Equivalently b(q) is the minimum sum of s(d_i) over ordered factorizations q=product d_i with d_i>1; the birth classes are precisely the minimizing words with each letter chosen from its minimal strong fiber. No unordered-factorization quotient is allowed: SCC order matters.

For the symmetry formula extend s(1)=0 and c(1)=1 as an EMPTY center, not a singleton. Then for q>1,

S(q) = sum over x^2 dividing q, with y=q/x^2 and 2b(x)+s(y)=b(q), of N(x)c(y).

Again keep finite-cost terms only. A self-converse component word has a uniquely determined ordinary first half U, followed by an empty center or one self-converse strong component, followed by the reverse-converse of U. A singleton center would violate minimality. The equation on orders forces the first half and center to be individually minimal, since otherwise replacement would lower the total order. Conversely every summand constructs exactly such a minimal self-converse tournament. This proves the formula. The converse-orbit birth count is M(q)=(N(q)+S(q))/2.

This is an exact min-plus recursion with multiplicity and involution data, not a claim to compute s,a,c without strong classification. It cleanly identifies the remaining research target: inventories of minimal strong realizations, followed by an automatic ordered-word synthesis.

## 4. Consequences useful for further research

(a) If q is prime and attained, every minimal realization is strong, because a nontrivial component word would factor q into at least two integers >1. More generally, if q has no factorization into two attained values >1, every minimal realization is strong. Arithmetic primality is sufficient but not necessary.

(b) If q>1 is squarefree, every self-converse birth class is strong: the symmetry formula permits only x=1, so S(q)=c(q) when s(q)=b(q), and S(q)=0 otherwise. In particular, if a squarefree value is first attained by decomposable tournaments strictly before any strong realization, all its birth classes form two-element converse orbits. This statement is a conditional structural consequence, not a claim that such a value has been found.

(c) H=9 illustrates a composite value with a strictly cheaper strong birth: b(9)=s(9)=5, whereas C3+C3 costs six vertices. The equality N(9)=S(9)=2 comes from the strong inventory, not from SCC-word permutations. Strong-component composition alone cannot explain away its two birth orbits.

(d) Scalar H is insufficient as a lossless state for generation. Already the two H=9 birth atoms have the same order, score sequence, automorphism order and converse type, but different 4-cycle and module structures. Any structure-sensitive generation rule must retain the atom identity or an adequate sidecar, rather than using a single two-valued symmetry label for each odd H.

## Provenance and reproduction

Inherited SCC, generating-series and converse formulas: 05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md, especially Sections 1 and 2. That source explicitly already distinguishes ordinary classes N, self-converse classes S, and merged count M=(N+S)/2. Its census.json records the order-five H=9 fiber. The present birth recurrence is the minimum-degree/leading-coefficient specialization of that established architecture.

Run: python3 /tmp/session224-proof/check.py

The check.py and check.json files are compact, independent witness verification. No repository files were edited and no commits were made by this audit lane.

## 5. Independent audit of the all-order palindrome formula

Let F_m(z) count all ordinary tournament classes with H=m by order, including the empty tournament only in F_1, and let Q_m(z) count all self-converse strong components, INCLUDING Q_1(z)=z for the singleton. Let D_m(z) count self-converse tournament classes, also including the empty tournament at m=1. Then

D_m(z) = sum over d^2 e=m of F_d(z^2) [1_{e=1}+Q_e(z)].

This is correct. The empty-center term 1_{e=1} must coexist with Q_1=z; they represent distinct even-length and odd-length SCC words. The first-half coefficient is ordinary F_d, not converse-merged F_d: U and its converse can give distinct complete palindromic SCC words. The center must be a single self-converse strong component, not an arbitrary self-converse tournament, or the decomposition would cease to be unique and be overcounted. These checks audit all three potentially ambiguous choices in the formula.

For squarefree m>1, the only possible d is 1, and every H=1 tournament is transitive. Hence F_1(z^2)=1/(1-z^2) and

D_m(z) = Q_m(z)/(1-z^2).

For m=1, the formula instead gives D_1=F_1(z^2)(1+z)=1/(1-z), as required. Thus the displayed squarefree formula must explicitly exclude m=1. In coefficient form, for squarefree m>1 every self-converse realization is uniquely a self-converse strong center of H=m with the same number of singleton SCCs before and after it.

This palindrome formula and its SCC grammar are inherited from the cited Session 24 PROOF, equation (2), where the two center cases are packaged as C_1^+=1+z. The minimum-order formulas above are its direct minimum-degree and leading-coefficient specialization. Neither the base formula nor the min-plus specialization should be presented as a new independent structural discovery. Their value here is a precise corrected research program and a symmetry refinement of the proposed birth atlas.

Optimization-safe verification: check.py uses explicit require(...) guards raising RuntimeError, never Python assert. Normal and python3 -O executions both print verification=PASS and yield byte-identical check.json and check-optimized.json.
