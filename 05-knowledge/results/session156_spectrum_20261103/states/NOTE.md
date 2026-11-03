# Two-apex automaton and bounded-feedback cutoff

Status: elementary PROVED construction and bound; finite checks are separate from proof. No claim of positive density, exact full spectrum, or novelty.

## Bounded feedback cutoff

Let T be a strong tournament on at least three vertices, S any feedback vertex set of size k, and m=|V(T)\S|. Partition the transitive complement into at most 2^k types by orientations to S. If a type has t vertices, then

    H(T) >= 2^ceil(t/2) + 1.

Consequently, for H(T)<=X with X>=3,

    m <= 2^(k+1) floor(log2(X-1)).

Proof. Choose a median vertex v of the type in the complement's transitive ordering. Since T is strong, v lies in a directed triangle: v has both in- and out-neighbors, and if no out-neighbor beats any in-neighbor, a path from v to an in-neighbor cannot make its first transition from the out-neighbor set to the in-neighbor set. A triangle cannot consist entirely of complement vertices. If its other two vertices are in S, every vertex in the type has the same orientation to both. Otherwise exactly one other vertex w is in the complement and the third is in S. At least ceil(t/2) members of the type lie on the same side of w as v: if w precedes v, all members at or after v qualify; if w follows v, all at or before v qualify. They have identical orientations to both triangle vertices. Hence the induced tournament on these members and the other two vertices is a directed 3-cycle with one vertex replaced by a transitive block of size at least ceil(t/2). Its Hamilton-path count is 2^r+1, as one can count directly by insertion or use the recurrence in S155/primes. Every Hamilton path in a subtournament extends to one in T by successively inserting each remaining vertex, so H(T) is at least this count. Pigeonhole gives t >= ceil(m/2^k), and ceil(t/2) <= floor(log2(X-1)) implies t <= 2 floor(log2(X-1)), completing the cutoff.

Important correction: an S-neighborhood class need not be an externally homogeneous module because its vertices may be interspersed with vertices of different types. The median step repairs this. For a minimal strong example of the failed module claim, take k=2, apex 0 beating apex 1, and type word (1,0,1). Its first and third base vertices have the same S-type but the middle base vertex is beaten by the first and beats the third. This example has H=9. The stronger asserted exponent t without the factor 1/2 is not proved here (failure of the module proof alone does not refute the numerical inequality).

This gives a completely finite exact decision procedure for strong FVS-at-most-k spectrum membership at any fixed X: enumerate k-vertex tournaments S and all 2^(km) type words for 0<=m<=2^(k+1) floor(log2(X-1)), test strong connectivity, and compute H exactly. The one-vertex strong tournament and any tournament with feedback number <=k and n<k can be enumerated separately; otherwise a feedback set of exactly k vertices exists by enlarging one. The bound is not a practical claim: for k=2 its cutoff is 8 floor(log2(X-1)). Nor does it prove asymptotic zero density. This enumeration is compatible with but larger than the exploratory S155 FVS-2 census.

## General finite state transfer bank

Fix an arbitrary tournament S of size k, and list its vertices in a proposed order pi. The remaining m vertices form a transitive tournament x1->...->xm. Each base vertex carries its k-bit orientation type against S. Every Hamilton path restricts to an ordering pi of S and splits the base vertices into k+1 gaps around these marked vertices. Within each gap their order must be increasing in the transitive base order. Conversely, assigning each base vertex to one gap determines a unique candidate path; only the gap boundaries and directly adjacent S vertices must be checked.

For each gap maintain state E (empty), F (nonempty, current last base vertex fails the edge to the right S vertex), or G (nonempty, current last base vertex passes it). A virtual left source/right sink always passes. Process x1,...,xm. For gap i, let L_i mean its left S boundary beats x (true for exterior left), and let R_i mean x beats its right S boundary (true for exterior right). Let M(L,R) be the 3x3 integer matrix, columns indexed E,F,G: the E column is L times the basis vector G if R, F otherwise; the F and G columns are that same basis vector. Then the step matrix on the 3^(k+1) product states is

    A_pi(x) = sum_{i=0}^k I^(tensor i) tensor M(L_i,R_i) tensor I^(tensor(k-i)).

Each summand places x into exactly gap i. The initial state is E^(tensor(k+1)). At termination accept G, reject F; accept E for an outer gap and for an internal gap exactly when the left S vertex beats the right S vertex. Sum the accepted weights over all k! choices of pi. This equals H(T), since it counts each path once. For k=2 there are four input symbols and two fixed 27-state linear representations (one for each apex order). This is an exact richer bank than the two-letter mixed-prefix bank, but no density or coverage assertion follows automatically.

## Finite controls

Reproduce `python3 bank.py`. The two known 7-vertex strong witnesses are described by apex 1 beating apex 0 and base-type words 1,1,3,0,3 (H=35) and 1,1,1,2,3 (H=49), where bit a means base vertex beats apex a. The program compares this gap automaton with independent vertex-subset endpoint DP on all 682 two-apex codes of m=0,...,4. It then scans strong codes for m=2,...,7 and asserts the bound using the largest type class, reporting distinct H counts. The census is finite/exploratory, not an all-order spectrum classification.

## Targeted comparison with inherited work

S155/primes already gives the cyclic-block lower bound 2^r+1 and the stronger k=1 estimate H>=2^(m-2)+1. Session26 CUT_ADDENDUM gives the exact k=1 minimum and a positive lower-density gap. The present median bound and k=2 transfer formalism were not found in these specifically examined sources; no assertion of wider repository or literature novelty is made. The all-order k=2 cutoff is much weaker than the k=1 bound and is not an asymptotic density result.

## Compact finite witness bank through 301

`witnesses.json` records an explicit two-apex code for each of the 148 odd targets H in [3,301] except 7 and 21; all have m<=7. `python3 check_cert.py` (also `python3 -O check_cert.py`) independently builds each tournament, verifies strongness via transitive closure, and counts paths via subset-endpoint DP. It verifies only existence for these targets, not global nonexistence for 7 or 21. H=1 is attained by the singleton tournament. These are independently verifiable certificates, not a claim of novel coverage relative to the broader S155 census.

## Strongness is a regular type-word condition

For fixed apex tournament S, any nontrivial dominator cut A -> V\A has A intersect the transitive base in a prefix. Put R=A intersect S. Necessarily every apex in R beats every apex in S\R; every base vertex in the prefix beats all apexes in S\R; and every base vertex in the suffix is beaten by all apexes in R. Conversely these conditions define a dominator cut. For each R let U_R be the type symbols beating all of S\R and V_R those beaten by all of R. Non-strongness is the finite union U_R* V_R* over dominating R, restricting to a nonempty prefix when R is empty and a nonempty suffix when R=S. Thus strong type words form a regular language. For two apexes with 0->1 and alphabet Sigma={0,1,2,3}, the exact non-strong language is

    3+ Sigma*  union  Sigma* 0+  union  {2,3}* {0,2}*.

`check_language.py` independently compares this regular expression to graph strong connectivity for all 21,845 type words of lengths 0..7, with both normal and optimized Python. This regular language can be intersected with the 27-state weighted representation, without claiming that its scalar image is dense or exhaustive.

For each fixed k, fixed apex tournament S, modulus M and residue r, the set of base lengths m for which some strong type word gives H congruent to r modulo M is ultimately periodic. To see this, augment the strongness DFA state by every apex-order transfer vector reduced modulo M. The augmented state space is finite, transitions depend only on the next type symbol, and acceptance checks both strongness and the sum of final linear functionals modulo M. Replace all type labels by a single letter and determinize. The sequence of states reached by that letter is eventually periodic, proving the claim. This does not imply that any given residue occurs, nor does it imply exact integer coverage.
