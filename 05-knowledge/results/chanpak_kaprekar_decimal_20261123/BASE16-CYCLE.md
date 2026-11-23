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

## Optimal even up to relabelling

The remaining all-order obstruction is that no two consecutive arcs can realize
\[
195\longrightarrow135\longrightarrow15.
\]
Reverse the proposed chain. The following monotonicity applies to every strictly increasing chain of ordinary Hamiltonian-path counts.

**SCC coarsening and prime valuations.** A cross-component reversal is nondecreasing by the [rooted gain formula](FIRST-JUMP-SQUARES.md#cross-component-gain-and-singleton-budget). Unless the endpoints are two adjacent singleton SCCs, its old forward arc has an alternate forward route, and the new backward arc merges the whole interval. For example, one can route through a middle component, or through an out-neighbor or an in-neighbor in a nontrivial endpoint component. The adjacent singleton reversal just swaps them and does not increase height. An internal reversal which splits an SCC has its endpoints in distinct new components, since otherwise its inverse could not merge the new components. Such a split cannot increase height: its inverse is a nondecreasing cross-component reversal. Thus strict increase only coarsens SCCs. Every terminal singleton was untouched and can be deleted from all stages without changing the heights.

More precisely, each edit lies in exactly one terminal SCC. At all earlier stages the terminal SCC blocks are already ordinal to one another. Each quotient of successive heights is therefore the quotient of the local heights in the block containing that edit. Assign those rational quotients to the corresponding terminal blocks: starting from their integer local heights, all local partial products must be positive integers. This tracks where a prime valuation can be added or removed.

For \(15\to135\to195\) the first quotient is 9. If the second arc were in a different terminal component, the component used by the first arc would keep a factor of nine, contrary to \(9\nmid195\). Hence both edits lie in one terminal SCC. The product \(q\) of the untouched counts divides all three named heights and therefore divides 15. Divide the entire walk by \(q\); the active terminal component has height \(195/q\). The [exact strong minima of Busch](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) give the four small order caps below. The checks allow arbitrary ordinal padding up to the cap, so all smaller realizations are covered.

| \(q\) | Normalized chain | Order cap | First arcs tested | First arrivals | Next arcs tested | Full hits |
|---:|---|---:|---:|---:|---:|---:|
| 1 | \(15\to135\to195\) | 10 | 1890 | 34 | 1530 | 0 |
| 3 | \(5\to45\to65\) | 8 | 140 | 2 | 56 | 0 |
| 5 | \(3\to27\to39\) | 7 | 105 | 10 | 210 | 0 |
| 15 | \(1\to9\to13\) | 5 | 10 | 1 | 10 | 0 |

As a separate check not needing the prime-valuation argument, height 195 has possible nontrivial SCC height lists \((195),(3,65),(5,39),(13,15),(3,5,13)\), up to order. Their total active orders are at most \(10,11,11,11,12\), respectively. Deleting terminal singletons and padding to 12 reduces the full chain to 72 height-15 sources. There are 4752 first-arc tests with 174 arrivals at 135 and zero arrivals at 195 among the 11484 next-arc tests.

The \(15\to225\) jump and the \(195\to135\to15\) suffix therefore force **different** extras. Four marked transitions plus those two extras make six optimal even if no equality or isomorphism of endpoints is required beyond returning to height 15. The ten-vertex witness attains six with exact labelled closure.

**Minimum order.** At six edits the only segment-length schedules are \((2,1,2,1)\) and \((2,1,1,2)\). On all 30 padded order-nine height-15 sources, both give 44 arrivals at 225, then two at 195. With two edits next there are four visits to 135, none with a single-arc return to height 15. With one edit next there are zero visits to 135. These are counts with multiplicity, allowing reuse across segments. Within a two-arc segment the arcs commute, and repeating one would not change the phase height. Padding rules out lower orders as well. Thus **ten vertices are necessary at the six-arc optimum**, also up to relabelling.

Run [`node kaprekar-base16-optimal.js`](kaprekar-base16-optimal.js) for the all-order suffix and minimum-order censuses. It constructs the complete height-15 bank by anchoring one Hamilton path in each possible strong core, enumerating all other arcs, and deduplicating by isomorphism. The bank consists of two strong cores at each of orders five and six and the two orders of a directed triangle and strong height-5 core. The smaller ledger tests use the unique strong \(H=3\) and \(H=5\) cores or just a transitive source. Heights are computed as products of exact subset-DP counts on SCCs. Earlier restricted or low-order checks are retained in [the original walk replay](kaprekar-base16-cycle.js), but the previous order-20 search for a five-edit relabelled walk is no longer needed.

The [second base-16 orbit](BASE16-SECOND-CYCLE.md) also needs six edits under either closure convention, on a minimum of eight vertices.

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

### An infinite lift of the decimal five-phase orbit

The denominator in the transport rule can also be handled exactly. Set \(k\ge7\), \(x=2^k-12\), \(d=9x+1\), and \(B=11d-1=99x+10\). Then \(B-1=9(11x+1)\); the required ordinal factor for the *absolute* decimal walk is
\[
\frac{d(B-1)}9=d(11x+1).
\]
Both factors have signed-shell witnesses with nonempty transitive filler:

| Factor | Core outgoing rows | Tag mask | \((D,G,L)\) | Filler order |
|---|---|---:|---|---:|
| \(d=9\cdot2^k-107\) | \((50,44,9,17,38,76,31)\) | 70 | \((113,576,220)\) | \(k-6\) |
| \(11x+1=11\cdot2^k-131\) | \((6,4,8,19,39,79,31)\) | 54 | \((29,352,160)\) | \(k-5\) |

Ordinally join the [eight-vertex labelled decimal walk](KAPREKAR-METAGRAPH.md) with both shells. The result has \(2k+13\) vertices and realizes the five prescribed phases at eight edits, with every height in the decimal walk multiplied by \(d(11x+1)\). The first such base is 11494. The eight-edit bound is optimal when edits are required to stay within the decimal moving block, by the original labelled lower bound; unrestricted optimality in the larger bases is not claimed.
