# First Kaprekar jumps: the valuation cap and exact two-edit squares

This extends [the decimal walk](KAPREKAR-METAGRAPH.md) and sharpens [the earlier base-16 obstruction](BASE16-OBSTRUCTION.md). An edit means reversing one specified arc in an ordinary tournament. \(H\) is its ordinary Hamiltonian-path count, not a residue or a count on quotient states.

## Cross-component gain and singleton budget

Write the strongly connected components (SCCs) in their ordinal order. Reverse an arc \(u\to v\) from component \(A\) to a later component \(B\). Let \(R\) be the product of the heights of the components outside the closed interval \([A,B]\). Define the following ordered two-slot directed path covers, with either slot allowed to be empty unless rooted:

- \(Q_2(X)\): unrestricted covers of \(X\).
- \(U_A(u)\): covers of \(A\) whose second slot begins at \(u\).
- \(V_B(v)\): covers of \(B\) whose first slot ends at \(v\).

The number of newly gained paths is
\[
R\,U_A(u)V_B(v)\prod_{A<X<B}Q_2(X).
\]
Indeed, cut a gained path at \(v\to u\). Each resulting path visits components in the old order. All exterior components belong to the one slot in which their position permits them to occur. \(A\) and \(B\) impose the displayed rooted constraints, and the middle components are free to distribute across the two slots. Conversely the slot covers concatenate uniquely. If there is a middle component, no old Hamilton path used \(u\to v\), so the gain is the entire change. For adjacent \(A,B\), subtract \(Rxy\), where \(x\) counts spanning paths of \(A\) ending at \(u\) and \(y\) counts spanning paths of \(B\) starting at \(v\).

A singleton has \(Q_2=2\); a nontrivial SCC has \(4\mid Q_2\), with \(8\mid Q_2\) at even order. In the [Rédei–Berge odd-cycle expansion](https://arxiv.org/abs/2307.05569), a packing of \(c\) nontrivial odd cycles with \(s\) uncovered vertices contributes \(2^{s+2c}\) to the two-slot count. This is divisible by four for order at least three, and by eight at even order at least four: a single odd cycle leaves a vertex uncovered, and the empty packing has \(s=|X|\).

Consequently, if \(a=v_2(H_{\rm after}-H_{\rm before})\), each intermediate singleton spends one unit of the budget \(a\), and each intermediate nontrivial SCC spends at least two (three at even order). For a single reversal one may always delete all singletons outside \([A,B]\). The number to retain is at most \(a+2\), or at most \(a+1\) unless the ratio is exactly \(1+2^a\): with two singleton endpoints and only singleton intermediates, the old interval has height one and the new height is \(1+2^a\); with any nontrivial middle component at least two units are spoken for. Adjacent singleton endpoints do not change the height. For an internal flip all ordinal singleton pads can be dropped. In particular when \(a=1\) and the ratio is not three, keep at most two singleton vertices beyond the nontrivial components on the smaller-height side. Singleton padding can then bring a witness at smaller order up to the cap. In general a source of height \(h\) has at most \(\lfloor3\log_3 h\rfloor\) nonsingleton vertices: each triangle has height three, and [Busch's bound](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf) gives \(H(X)\ge 5^{(|X|-1)/3}\ge3^{|X|/3}\) for \(|X|\ge4\). Multiply across components. This gives an order cap even without a core census.

## Exact difference two forces one internal core

The adjacent case has a stronger absolute restriction than the valuation test. Put \(\alpha=U_A(u)-x\), \(\beta=V_B(v)-y\). Then the net gain divided by \(R\) is
\[
\alpha y+\beta x+\alpha\beta.
\]
In a strong tournament every specified vertex starts a spanning path, and every specified vertex ends one. For example, in a longest directed path starting at a specified vertex, an omitted vertex that is beaten by a vertex of the path can either be appended or inserted at its last such position. If no extension exists, all omitted vertices dominate the path, contrary to strength. The same argument with directions reversed gives the ending statement.

The second slot consisting of \(u\) alone gives \(H(A-u)\ge x\), by deleting \(u\) from the paths counted by \(x\). For every slot size \(2,\ldots,|A|\), take a prefix of a spanning path starting at \(u\) and cover the remaining vertices by a Hamilton path in the first slot (empty if needed). Thus \(\alpha\ge |A|-1\) and symmetrically \(\beta\ge |B|-1\). These inequalities include singletons, for which the corresponding difference is zero. Since \(x,y\ge1\), the adjacent gain divided by \(R\) is at least \(|A||B|-1\). This also proves that a cross-component flip is nondecreasing. In the nonadjacent case the gain divided by \(R\) is at least \(|A||B|2^m\), where \(m\) is the total number of intermediate vertices: \(U_A\ge |A|,\ V_B\ge |B|\), and each subset of an intermediate component can be assigned to the first slot, with a path in each part.

**Corollary.** For two odd heights differing by two, a cross-component reversal is possible only for the pairs \(1\leftrightarrow3\) and \(3\leftrightarrow5\). In fact \(R=1\) since it divides both heights. The nonadjacent lower bound forces a singleton, a singleton, and a singleton, giving \(1\to3\). The adjacent lower bound forces a triangle and singleton if the change is nonzero, giving \(3\to5\) (the two-slot excess for a triangle at any chosen vertex is two).

For every other exact \(h\leftrightarrow h+2\), the smaller tournament must use an internal flip; outside its flipped SCC the product divides \(\gcd(h,h+2)=1\), so every outside vertex is an ordinal singleton pad. The flipped SCC must also stay strong: if it split, the reverse edit would be cross-component and would increase, not decrease, height. Thus an exact \(E_p-1\leftrightarrow E_p+1\) seam for \(p\ge3\) is a single edit between **two strong tournaments on the same active vertex set**, up to common ordinal singleton pads. No merger of pre-existing SCCs, or splitting in the opposite direction, can implement that seam. This does not bar rooted or shelled constructions in which the entire construction is already strong.

## Four exact first jumps

In base \(b\) divisible by four, the first width-two Kaprekar jump from digit gap one is \((b-1)\to(b-1)^2\), with valuation one. For \(b=12,16,20,40\) its ratio is not three, so the singleton budget applies. By contrast, a multiplier \(1+2^a\) is physically realizable at every attained height: append a transitive block on \(a+2\) vertices by ordinal join and reverse its first-to-last arc. The block has one old path and \(2^a\) new ones, by distributing its intermediate vertices between the two traversal slots.

For \(H=15\) the source is either a strong core of order five or six, or two cores of heights three and five (orders three and four). For \(H=39\) the source is either a strong core of order six or seven, or two cores of heights three and thirteen (orders three and five). For completeness, the order-five maximum is 15, the order-four maximum is 5, and [Busch's lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf) \(H\ge5^{(n-1)/3}\) gives the upper orders. Anchored enumeration records two strong classes of height 15 at each of orders five and six; one of height 13 at order five; no height 39 at order six, and two at seven. The 15 and 39 sources therefore have at most seven and eight nonsingleton vertices, respectively. Heights 11 and 19 are prime and hence have just one nontrivial SCC, of order five and six respectively, by the same order bounds.

Use all ordinal placements of singleton pads and both orders of any two nontrivial cores. The following census is inclusive; padding reduces all lower orders to the reported order. The one-edit rows are all-order obstructions by the singleton cap. The two-edit rows are minimal-order lower bounds independently of that cap.

| Jump | Edits | Order | Source shapes | Arc sets tested | Hits |
|---|---:|---:|---:|---:|---:|
| \(11\to121\) | 1 | 7 | 33 framed | 693 | 0 |
| \(19\to361\) | 1 | 8 | 114 framed | 3192 | 0 |
| \(15\to225\) | 1 | 9 | 30 | 1080 | 0 |
| \(39\to1521\) | 1 | 10 | 20 | 900 | 0 |
| \(11\to121\) | 2 | 6 | 22 framed | 2310 | 0 |
| \(19\to361\) | 2 | 8 | 114 framed | 43092 | 0 |
| \(15\to225\) | 2 | 7 | 12 | 2520 | 0 |
| \(39\to1521\) | 2 | 9 | 12 | 7560 | 0 |

The framed rows keep duplicates from different Hamilton orders; the other rows quotient the cores by relabelling. For each witness below, row \(i\) is the outgoing bitmask of vertex \(i\), with vertices numbered from zero. The two different arcs commute on the labelled tournament itself.

- **Base 12:** rows \((126,100,72,82,102,76,0)\), with flips \((0,2),(0,6)\). The two routes are \(11\to37\to121\) and \(11\to111\to121\). The optimum is two edits, on at least seven vertices.
- **Base 20:** rows \((510,508,504,304,32,64,152,312,112)\), with flips \((0,4),(2,4)\). The two routes are \(19\to311\to361\) and \(19\to91\to361\). The optimum is two edits, on at least nine vertices.
- **Base 16:** rows \((254,180,248,242,224,192,130,0)\). Flip \((0,1)\) then \((0,7)\): \(15\to33\to225\). Reverse that order: \(15\to231\to225\). The optimum is two edits, on at least eight vertices.
- **Base 40:** rows \((1002,1012,1017,1010,993,960,640,768,576,0)\). Flip \((1,6)\) then \((2,9)\): \(39\to147\to1521\). Reverse that order: \(39\to663\to1521\). The optimum is two edits, on at least ten vertices. The detour through \(663=39\cdot17\) lies on the other base-40 orbit: iterating \(d\mapsto|2d-41|\) from 17 gives \(17,7,27,13,15,11,19,3,35,29,17\), which excludes 1 and 39.

Run `node kaprekar-first-squares.js` for the [complete replay](kaprekar-first-squares.js). It uses anchored enumeration of all nonconsecutive edges with path \(0\to1\to\cdots\to n-1\), exact subset/end-vertex path DP, and relabelling deduplication of the small SCCs. Two-edit tests use unordered pairs of distinct arcs: repeating an arc returns to the source. These squares solve only the first jumps, not the complete physical cycles.
