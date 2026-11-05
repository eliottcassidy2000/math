# Independent literal check

**VERIFIED:** the standalone counterexample reproduces under two distinct counting methods.

The primary verifier builds the actual source and graft adjacency rows and uses subset dynamic programming. The second verifier reads only those rows, validates each tournament, counts the full five- and six-vertex rooted-source profiles by permutations, then enumerates every directed simple path of the eleven- and twelve-vertex grafts by unmemoized recursion. Hamilton cycles are counted from every starting vertex and divided by the graph order. It imports none of the primary verifier.

Results: u=(9,3,9,3), v=(15,6,18,9), both with exactly one Hamilton cycle and one root-start path; N=417 and M=867, both with exactly two Hamilton cycles. The offset is 33. The identity 69(2k+1)-66U+24=6(23k-11U+15)+3 proves the obstruction for every odd D=2k+1, without finite-search inference or an assumption on U parity.

Both programs pass normally and under Python -O, with byte-identical paired receipts. Run the commands in the parent README. No canonical theorem, proof, index, other researcher's files or shared dirty checkout was changed. This is not an audit of the rest of THM-4883.
