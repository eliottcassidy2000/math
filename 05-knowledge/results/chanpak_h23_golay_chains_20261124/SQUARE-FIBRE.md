# Exactly which Golay contacts preserve the original shifted square?

**PROVED.** At fixed equal widths \(k\), the unmarked target \(H=(2^{k+1}+1)^2-2\) survives exactly when both expanded modules are transitive. If the 24 marked arcs are internal to the two initially transitive chains, **all 4096 codewords preserve this target iff the marks are disjoint adjacent swaps**. For that construction \(k=24\) is sharp, and all 4096 labelled states are the same ordinary isomorphism class. A less degenerate single-apex placement preserves the target for exactly four words, the maximum in that placement type, at every \(k\ge25\).

## Exact preservation and the adjacency test

Use the [bilinear two-slot law](UNEQUAL-CHAINS.md), with \(D_i\) the ordinary path count and \(Q_i\) the ordered two-cover count with empty slots admitted. \(D_i\ge1\) and \(Q_i\ge2^k\): in the latter, every subset and its complement give at least one path. Both lower bounds are attained by transitive chains. Conversely a directed triangle extends to at least three distinct paths by repeatedly inserting an omitted vertex while keeping the relative order of old vertices, so \(D_i=1\) forces transitivity. The bilinear law is strictly increasing in each of its four coordinates on the feasible set: the two derivatives with respect to \(D_i\) are \(Q_B-D_B\) and \(3Q_A-D_A\), positive since \(Q_i\ge2D_i\). If increasing all coordinates from the transitive minima, first increase the \(Q_i\), then the \(D_i\), so each intermediate profile still satisfies this inequality. Thus the target is retained exactly by transitive modules.

If a marked arc skips an interior chain vertex, flipping just that arc among the at-most-three marked contacts of the resulting triple makes a triangle. If two adjacent marked arcs share an endpoint, flipping both produces a triangle on their three consecutive vertices. Either specified local bit pattern occurs in an extended Golay word, since every projection on at most seven positions is free. Conversely a matching of adjacent swaps simply permutes the order within each chain. Its 24 marks require \(2\lfloor k/2\rfloor\ge24\), with equality possible at \(k=24\). All those words are isomorphic even though the labelled distance between any pair is at least eight.

## Four survivors on an ordered octad trio

Put all 24 marked contacts on a single apex facing 24 consecutive spine neighbours in one chain. For \(k>25\), prefix the same module by \(k-25\) singleton vertices in ordinal order. A bit one means apex-to-neighbour, and transitivity on the star holds iff the contacts along the spine are a threshold word \(0\cdots01\cdots1\). Distinct thresholds differ in precisely the distance between the two cut positions. Since the code has minimum distance eight, at most four cut positions can belong to it.

In the standard generator convention of [the apex-contact chart](APEX-CONTACTS.md), order the positions along the spine by the disjoint octads
\[
 O_1=(0,2,4,5,6,10,11,23),\quad
 O_2=(1,3,7,8,12,13,15,18),\quad
 O_3=(9,14,16,17,19,20,21,22).
\]
The four surviving words are \(0,O_3,O_2+O_3,O_1+O_2+O_3\). The other 4092 strictly raise the count at every \(k\ge25\). These four form a nested flag of supports, not a linear subcode; their count does not say how many attain other shifted squares. At \(k=25\) specifically, an exhaustive exact count gives 4067 distinct heights, with no other value two below a square. The four survivors give \(4503599761588223=67108865^2-2\).

Among \(M_{24}\) coordinate permutations the setwise stabilizer of the four words is precisely the stabilizer of the **ordered** octad trio: the sole weight-8 word determines \(O_3\) and the complement of the sole weight-16 word determines \(O_1\). Compare the [dodecad chart](DODECAD-CHAIN.md), where preserving the full height readout retains \(M_{12}\), and the two-apex chart, where no nonidentity action preserves every exact count.

There is a concrete but limited analogy to [Curtis's symmetric generation using \(M_{24}\)](https://www.cambridge.org/core/books/abs/art-of-working-with-the-mathieu-group-m24/symmetric-generation-using-m24/5875D4E440B8CAAA4238164C2B1BEAB9): a named finite group permutes an indexed set of involutory generators, and additional relations select a structure. Here the physical arc flips *commute*, and the height gate chooses a particular flag of supports or a dodecad; these are not the noncommutative progenitor or a new group-generation construction. The analogy informs which data a quotient must retain.

Replay: run square-fibre-check.py with Python 3.
