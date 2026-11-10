# A physical smooth-sector descent with exact certificate cost

Status: **PROVED elementary application and cost theorem**, using the finite exact graph certificates below. Independently audited; see ../audit/REVIEW.md. Chan–Pak remains OPEN. No new unrestricted coverage interval, smallest tournament order, or historical-priority claim.

## 0. Inheritance, board, and what changes

THM-4094 (`01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`) already gives the prime/7-prime reduction, including carry atoms49,63,343. Session241's `provenance/README.md` recovers it and the earlier Session160/ordered-response-ten {3,5,7} monoid results. These are NOT new here. THM-1862 order-join multiplicativity is the inherited physical operation. THM-5103 (`...faithful-fermat-path-count-arithmetic-and-sharp-coprime-bank-cost.md`) gives a different faithful arithmetic embedding, whose encoded count generally is not the requested input integer.

Session228 `physical_state/PROOF.md` gives an arbitrary four-state physical membership oracle; Session229 `inverse/PROOF.md` gives strictly decreasing proper-graft inverse search without branch coverage. Session237 THM-5092 (`...four-test-tail-content-and-invisible-reroot-escape.md`) warns that even order plus four fixed-root coordinates loses rerootability. Session198 `bilinear-descent/PROOF.md` descends existing representations, not missing targets. None is re-proved or strengthened as a global theorem here.

Board: actual graphs; unavailable factor7; universal join lifting; prime31/41 suppliers; seven-adic carry cost. Anchor: a failure-preserving descent on a declared closed sector. Niche: optimal physical certificate size in a fixed bank. Wildcard: compare genuine preservation with congruent-number descent, without asserting a transfer between problems.

The additional deliverable is explicit certification of the {3,5,7,31,41} sector, with literal217/287 suppliers, and an exact all-order formula and eventual period for the least graph order **inside this specified bank of ordinal sums**. Sector coverage is a concrete application of the inherited atom mechanism, not a newly discovered general principle. The cost theorem is about the declared representatives, not globally optimal representatives of their counts.

## 1. Physical seed bank and independent finite certification

An adjacency row a_i has bit j equal to1 iff i→j. Start a word at the edge0→1 with root0. O appends a sink and keeps the root. Z appends a vertex beating only the old root and keeps the root. J appends that same vertex but marks it as the new root. These are literal operations on tournaments, not operations on numerical tuples. The inherited state updates appear in THM-4888 but are unnecessary for the verification below.

| H | order | word or description |
|---|---:|---|
|3|3|Z|
|5|4|OZ|
|31|7|ZOZOZ|
|41|6|transitive order6 with every distance1 or2 arc reversed|
|35|7|OZJOZ|
|49|8|OJJOOZ|
|63|8|OZZOOZ|
|217|10|OJOOZZOZ|
|287|10|OOOOZJOZ|
|343|11|OOJOZJJOZ|

In particular the two extra7-prime suppliers are not assumed:

    T217=[1014,1020,1016,305,992,960,904,776,512,8],
    T287=[830,1020,1016,1008,992,960,897,257,512,128].

All ten complete row lists appear in CHECK.json. For each graph, `verify.py` checks every unordered pair has exactly one direction, and computes Hamilton counts by the exact recurrence

    P({v},v)=1,
    P(S∪{v},v)=sum_{u∈S,u→v} P(S,u),
    H(T)=sum_v P(V,v).

This partitions paths uniquely by their last vertex, proving the algorithm. An independent permutation enumeration checks all seeds of order at most7. Graph construction plus subset DP, not the O/Z/J numerical recurrence, certifies217 and287. All arithmetic is integer arithmetic. These seeds are existing attainable values, not claims of newly found values.

## 2. The universal lifting operation

For actual tournaments A,B, let A⇒B mean their disjoint union with every cross-edge from A to B. Every Hamilton path first traverses all of A and then all of B: once it enters B it cannot return, and it cannot start in B with A unvisited. Thus concatenation is a bijection and

    H(A⇒B)=H(A)H(B).

Use H(empty)=1 internally; if nonempty tournaments are required, represent the terminal count1 by a singleton. Singleton padding raises an order by1 and never changes the count, but is not included in the minimum-bank cost below, where the empty sum is allowed.

Crucially this lift works for **every** graph B with the desired smaller count. It needs no port, response vector, symmetry, or inverse oracle. A source graph is never replaced by an arithmetic tuple and then presumed physical. Every accepted recipe retains the actual adjacency rows of its factors and the specified join directions. Endpoint data is simply not needed by this universal operation.

## 3. Full admissibility and a failure-preserving descent

Define

    C={3^a5^b31^d41^e7^c : a,b,d,e,c≥0} minus {7,21}.

Every N in C is attained by a sum from the above ten-seed bank. Here is an explicit terminating construction, not merely descent of an already existing representation.

* If c=0, remove any positive ordinary prime exponent using its physical3/5/31/41 seed; finish at1.
* If c≥2 is even, remove49, repeatedly if necessary.
* If c≥3 is odd, remove343 once, then49 as necessary.
* If c=1 and one of b,d,e is positive, remove respectively35,217,287 (choose the first available in that order).
* If c=1 and b=d=e=0, membership in C forces a≥2. Remove63.

The last two steps leave c=0; the c≥2 steps leave either c=0 or c≥2. Every residual stays in C. Each removed block has count at least3, so the residual N' is a strictly smaller positive integer. Reversing the steps by the literal join constructs N. The measure is N (equivalently one can use sum of exponents); termination is well founded, with the unique terminal target1.

For a true infinite-descent formulation set

    Fail(N) := N∈C and there is no finite tournament with H=N.

Each permitted step chooses a **verified physical** block of count b and leaves N'=N/b∈C. Universal lifting gives Fail(N)⇒Fail(N'). If N had no witness but N' had one, adjoining the verified b-block would contradict Fail(N). The terminal target1 has a witness. Thus a hypothetical missing target in C cannot exist. Failure is never asserted of a physical graph; it belongs to an unresolved existential target. The evolving state consists of that residual target together with an actual ordered context of seed graphs, with invariant H(context)·N'=N_initial. This distinction avoids presuming the very witness that an inverse descent must build.

For completeness, the bank cannot produce7 or21: every generator has nonnegative prime valuations; any generator whose prime support is contained in {3,7} is3,49,63,343. Neither7 nor21 has enough7-exponent for49/343 or enough3-exponent for63. This is a bank exclusion. The unrestricted impossibility of7,21 is inherited from THM-1370, not established by this argument.

## 4. Exact least certificate order and a seven-adic period

Let C(a,b,d,e,c) denote the smallest total number of vertices in an ordered sum of the **specified** seed graphs with count3^a5^b31^d41^e7^c, and infinity if no such sum exists. Order of the blocks can change the tournament but not its count or total order.

First define F(c), the minimum cost using only49 of cost8 and343 of cost11. Put F(negative)=F(1)=infinity. For q≥0,

    F(3q)=11q,
    F(3q+2)=11q+8,
    F(3q+1)=11q+5 when q≥1.

Indeed if 2v+3t=c, the cost is8v+11t=(11c+2v)/3. Minimize v≥0 subject to2v≡c mod3 and2v≤c. The least such v is0,1,2 respectively in residues0,2,1; residue1 first becomes feasible at c=4. In particular F(c+3)=F(c)+11 for c=0 and every c≥2, but NOT at c=1.

Every factorization in the full bank is completely specified by the numbers y,x,z,w of blocks63,35,217,287, followed by ordinary prime blocks and a49/343 factorization. Hence the exact formula is

    C(a,b,d,e,c)=3a+4b+7d+6e
       + min [2y+3x+3z+4w+F(c-y-x-z-w)],

where the minimum ranges over

    0≤y≤floor(a/2), 0≤x≤b, 0≤z≤d, 0≤w≤e.

An infinite candidate is discarded. The correction coefficients are literal order differences: 63 costs8 instead of two3 blocks costing6; 35 costs7 instead of a5 block costing4; 217 costs10 instead of a31 block costing7; 287 costs10 instead of a41 block costing6. This accounts for all factorizations because no generator carries an ordinary prime twice, except63 which carries3 twice, and the only remaining generators are49 and343.

Let K=floor(a/2)+b+d+e. For all c≥K+2, every candidate residual c-y-x-z-w is at least2. Applying the F identity to **every** candidate, with exactly the same bounded parameter box, proves

    C(a,b,d,e,c+3)=C(a,b,d,e,c)+11,     c≥K+2.

This is exact eventual periodicity of an optimal physical-certificate cost, not a numerical heuristic or a claim of period3 for the tournament spectrum. It supplies a finite bounded-intercept optimization for an infinite target family. For fixed a,b,d,e, three tail constants determine every sufficiently large seven exponent. The threshold is sufficient, not asserted sharp. All counts in this sector were already implied attainable by the existing finite-prefix and atom mechanisms; the contribution here is the explicit compiler and optimal cost law in its fixed graph bank.

## 5. Hostiles and boundaries

1. **Blind division is not admissible descent.** 35 is realized but35/5=7 is forbidden. Our c=1 rule removes35 as a block, rather than stripping a5 and silently leaving an invalid base case. Similarly63/3=21. The valuation vector alone does not justify every divisor move.
2. **Greedy cost is not optimal cost.** At7^6=117649, three49 blocks cost24 vertices; two343 blocks cost22. The exact minimum is22. The canonical descent of Section3 is sound but not generally size optimal.
3. **Missing physical supplier.** If the217 block is removed, the remaining bank cannot attain217: its prime support is{7,31}, and the only remaining eligible block is31, which leaves7. A valuation-level request for217 does not manufacture a graph. This is the precise obstruction to generalizing the finite prime set to all primes without supplying every required7p atom.
4. **Not a strong-output theorem.** A join of two nonempty blocks is not strong. Atoms may have strong witnesses, but these composite recipes need not. No proof of strong realization of every sector value is offered.
5. **Not a congruent-number transfer.** 31/41/49 enter as actual Hamilton counts and prime/carry suppliers. Their relation as square roots in an arithmetic progression with gap5·12² is a different construction. No map sends an elliptic rational point to these tournament graphs, and reducing an elliptic height supplies neither a217 nor a287 graph. The shared discipline is preserving the exact failure predicate and admissibility, not a mathematical implication between the problems.
6. **Global frontier.** THM-4094 still requires infinitely many prime and7-prime suppliers. This bounded-prime descent neither supplies them nor bypasses S229's physical branch-coverage gap.

## 6. Reproduction, scopes, and connections

From this directory:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O verify.py > optimized.json
    cmp optimized.json CHECK.json

The standard-library verifier checks all ten actual seed tournaments; independent permutations on all seeds of order≤7; eight literal joins up to order14; 1,053 target exponent vectors (a,b,d,e each0..2,c0..12), comparing unrestricted bank dynamic programming against the formula; and2,560 period tests (a,b,d,e each0..3,c=K+2,...,K+11). It also checks each returned recipe product and vertex cost. Its checks use explicit exceptions, not assertions, so -O retains them. Finite tests validate implementations; the unbounded claims are proved in Sections2–4.

Source object→target: a finite list of actual certified seed graphs→an actual ordered-sum tournament. Preserved: multiplication of path counts, the smooth admissibility predicate under the specified inverse rules, and exact sum of vertex costs. Lost if only valuations are retained: which supplier graph exists and its physical size. Required sidecar: seed adjacency rows and ordered recipe. Cheapest decisive test: actual Hamilton DP on217/287, then the35/5 hostile. The inherited Session237 reroot collision is avoided because no reroot operation is used.

No navigation files, theorem IDs, or commits are changed by this lane; parent publication owns those actions.
