# The first base-16 cycle as a shortest labelled arc walk

In the [width-two Kaprekar quotient](KAPREKAR-IMAGES.md), the base-16 gap map is \(g\mapsto |2g-17|\). The first orbit \(1,15,13,9,1\) corresponds to the ordered ordinary Hamiltonian-path heights \(15,225,195,135,15\). A phase is a required checkpoint, not a claim that the arithmetic sends a tournament to a unique representative.

## The literal walk and its optimal length

Start from the outgoing rows
\[
(1022,1020,1016,1008,992,448,640,768,576,32)
\]
on ten vertices numbered from zero. Reverse, in order,
\[
(1,3),(0,4),(0,1),(1,3),(0,1),(0,4).
\]
The actual height itinerary is
\[
15\longrightarrow45\longrightarrow225\longrightarrow195
\longrightarrow75\longrightarrow135\longrightarrow15.
\]
Each arc is reversed twice, so the final labelled tournament is identical to the initial one. The first five vertices carry all edits and have height trace \(1,3,15,13,5,9,1\); they start as a transitive block. The remaining five vertices form an untouched ordinal factor of height 15, with rows \((14,20,24,18,1)\) in its own coordinates. Both detours \(45,75\) are genuine phases of the other base-16 orbit \(45,165,75,105,45\); \(45\) also occurs in the decimal orbit.

Every walk needs at least five reversals, because the [first jump](FIRST-JUMP-SQUARES.md) takes at least two and the other three required transitions each take at least one. A labelled closed walk in the orientation cube has even length. Hence **six is optimal with labels fixed**, among all orders.

**Ten vertices are also necessary at the six-edit optimum.** The order-nine census covers all 30 padded height-15 source shapes. If the first 225 occurs after two edits, it tests each of 44 two-arc arrivals and all remaining schedules satisfying closure. After six labelled edits at most three distinct arcs can have been touched; tracking the symmetric difference of the reversed arcs tests exact closure. If 225 first occurs after three edits, the suffix must be exactly \(195,135,15\); the first three arcs must be distinct, since a single-arc \(15\to225\) is impossible. Their three undos are then tested through those two heights. Ten triples pass the joint test for the heights 225 and 135, counting the choice of the arc left at 135, and none has the required intermediate 195. There are no closed walks in either case. Padding preserves every step, so this also excludes smaller orders.

## What survives if labels are forgotten

The unrestricted optimum up to relabelling is either five or six. A hypothetical five-step walk assigns both extra positions needed in the first jump to its two-step segment; every other segment is one edge. At order 11 the complete height-15 bank, with every ordinal singleton placement, gives 348 two-edit arrivals at 225, 324 subsequent arrivals at 195, and 954 at 135, counting multiplicities. Testing every outgoing arc each time (including undo of earlier flips) gives zero terminal heights 15. The first two arcs can be tried unordered because the resulting orientation is the same, and no intermediate height is prescribed within that first segment. This rules out five through order 11 by padding, even without checking isomorphism of the endpoints.

Five is impossible at every order **if all edits must remain inside the moving block of the ordinal-multiplier construction**. The moving block must then end with exact heights \(13\to9\to1\). Any single reversal from height 1 to height 9 creates a single nontrivial core: the reversed long arc in a transitive five-vertex block. Its rows are \((14,28,24,16,1)\); the ten single-arc neighbor heights in lexicographic order are \(5,9,11,1,9,15,11,9,9,5\), none 13. A flip among pads that does not span the core leaves an untouched factor nine. A cross-component flip that touches or spans the core merges it with at least one other vertex: before reversal there is an alternate forward path between the incident components, and the backward arc closes the whole interval. Height 13 is prime and every strong 13-core has order five, by the order-four maximum and [Busch's lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf). Thus the merged core cannot give 13.

**A finite all-order reduction.** Any five-step return would descend \(225\to195\to135\to15\) by three internal flips, since cross-component flips cannot decrease height. The nontrivial support therefore only shrinks after 225. Delete vertices that are singleton SCCs at both the start and 225. The restriction preserves the heights at those two points; neither of the first two flips can touch a deleted vertex, since it would then leave at most one flip for \(15\to225\). The later three are internal and cannot touch these vertices either. The deleted vertices are thus untouched and have the same outdegree at the start and final 15. Any endpoint isomorphism fixes each of them: it sends singleton SCCs to singleton SCCs, whose outdegrees are distinct (the number of all vertices after that component). Restriction therefore preserves even the relabelled return.

The reduced vertex set is the union of the nontrivial supports at 15 and 225. The 15 support has at most seven vertices. At 225, the product \(q\) of the SCC counts outside the first descending flip divides both 225 and 195, hence \(q\mid15\). For \(q=1,3,5,15\), Busch bounds the changed core of height \(225/q\) by \(11,9,8,6\) vertices respectively; the outside supports contribute at most \(0,3,4,7\). Thus the 225 support has at most thirteen vertices, and the unrestricted five-edit problem reduces to order at most **20**, or exactly 20 by padding. The order-20 census is not done.

For the second base-16 orbit, fixing an ordinal triangle would reduce the first jump \(45\to165\) to \(15\to55\). This normalized jump has \(v_2(55-15)=3\) and ratio different from nine. The [singleton budget](FIRST-JUMP-SQUARES.md#cross-component-gain-and-singleton-budget) therefore caps the source at seven nontrivial vertices plus four pads. The same height-15 bank at order 11 gives \(56\cdot55=3080\) arc tests, with zero hits. This normalization does not rule out unrestricted one-arc \(45\to165\).

## Transporting a literal orbit to larger bases

For an even base \(b\) and odd \(d>0\), put \(B=d(b+1)-1\). The gap map satisfies
\[
|2dg-(B+1)|=d\,|2g-(b+1)|.
\]
Positive gaps \(1\le g\le b-1\) lift to valid gaps, and their literal heights in the new base are \(Cg\) with \(C=d(B-1)\). Thus if a labelled tournament walk has exactly the normalized gap heights \(g\), an ordinal multiplier of height \(C\) transports that physical walk to base \(B\). A walk whose heights are already \((b-1)g\) instead requires multiplier \(C/(b-1)\); lifting the numerical gap map alone does not guarantee that this factor is attainable.

To build \(C\), use \(d=2^k+1\), realized by reversing the long arc of a transitive block of order \(k+2\). A signed shell realizes the other factor. Given core \(A\), tag subset \(J\), one singleton \(b_0\), and nonempty filler \(X\), use \(b_0\to X\to A\), orient \(J\to b_0\), and orient \(b_0\to A\setminus J\). Let \(D=H(A)\). Let \(G\) count ordered splits into two directed paths \(R,S\) covering \(A\) with \(R\) nonempty ending in \(J\) and \(S\) possibly empty. Let \(L\) count those splits where \(S\) starts in \(J\). Then the shell height is
\[
GQ_2(X)+(D-L)H(X).
\]
Indeed, two nonempty filler runs sit on opposite sides of \(b_0\) and give \(G(Q_2-2H(X))\); if a single run lies before \(b_0\) there are \(G-L\) core patterns, and after \(b_0\) there are \(D+G\). The latter cases contribute \((D+2G-L)H(X)\). A transitive filler of positive order \(w\) gives \(G\,2^w+D-L\).

- **From base 16:** take rows \((82,52,121,115,96,65,2)\) with tag mask 66. The exact readings \((D,G,L)=(49,272,34)\) with \(w=k-4\ge1\) give \(B-1=17\cdot2^k+15\). Ordinally join the moving five-vertex block, the \(d\)-block and the shell. It uses \(2k+11\) vertices and the same three edited arcs and six moves. The lifted heights are \(C,3C,15C,13C,5C,9C,C\). At \(k=5\) this is base 560 and the trace is \(18447,55341,276705,239811,92235,166023,18447\). Both detours are again genuine phases of the other lifted orbit. This proves a **six-edit upper bound**, not optimality in each of these bases.
- **From base 4:** the normalized first orbit is \(1,3,1\), physically realized on three vertices by making and unmaking a directed triangle. Take rows \((14,28,24,16,1)\) with tag mask 24. Here \((D,G,L)=(9,40,6)\); \(w=k-3\ge1\) gives \(B-1=5\cdot2^k+3\). For each \(k\ge4\) this yields literal labelled closure in two edits on \(2k+8\) vertices in base \(B=5(2^k+1)-1\). Two is optimal because the two marked heights differ.
  
The shell coefficient and walk checks are in [the runnable replay](kaprekar-base16-cycle.js). All counts are absolute, and no change to an ordinal multiplier is needed during a lifted walk.
