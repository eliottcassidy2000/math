# Order caps for prescribed single-arc itineraries

**PROVED scoped.** A single-arc itinerary is a fixed finite list of labelled tournaments on one vertex set, each obtained from the previous by reversing one arc. Some stages have prescribed Hamilton counts. The following cap gives finite search even if arbitrarily many unmarked singleton pads were originally allowed.

Let \(U(T)\) consist of the vertices in nonsingleton strongly connected components. Put \(H(T)\) for the ordinary Hamilton path count. For each nontrivial component \(C\), [Busch's lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf) gives \(H(C)\ge 3^{|C|/3}\): the 3-vertex case is equality, and \(5^{(|C|-1)/3}\ge3^{|C|/3}\) for \(|C|\ge4\). The ordinal component product therefore gives
\[
 |U(T)|\le \lfloor 3\log_3 H(T)\rfloor .
\]

**Checkpoint cap.** Suppose an itinerary of length \(\ell\) has prescribed positive heights \(h_i\) at a nonempty set of marked stages \(i\). Retain every vertex in \(U(T_i)\) for every marked \(i\), and all the at most \(2\ell\) flipped endpoints. If the set is empty, retain one vertex. At each marked stage, every removed vertex was a singleton SCC; removing any such vertices just discards factors of one and preserves the ordinary count. Every prescribed reversal remains a literal one-arc change after restricting to the common vertex set. Therefore if such a path exists at any order, one exists at order at most
\[
 \max\!\left(1,\;2\ell+\sum_{i\ {\rm marked}}\lfloor3\log_3h_i\rfloor\right).
\]
No height guarantee is needed or asserted at the unmarked stages. Labelled equality at a closed walk is retained; a nontrivial endpoint relabelling need not preserve the chosen subset, so this proof alone does not give the same cap for an unlabelled closing requirement.

**Consecutive different heights.** If every stage is marked and adjacent heights differ, dispense with the extra \(2\ell\) term. For a reversal inside an SCC, both endpoints are already in \(U\) before the edit. For a cross-SCC reversal, the interval of components between the endpoints becomes one strong component whenever there is a retained forward path from the original tail to its head. This always holds if either endpoint component is nonsingleton or an intermediate component exists, by ordinal contacts. A reversal between two adjacent singleton components only swaps two factors of one and cannot change height. Thus in every case of a height change the flipped endpoints belong to the retained union. This variant also yields a finite cap for a prescribed open path in the unlabelled metagraph, by choosing labelled lifts of its edges in succession.

This cap is more general but typically larger than the exact one-core reduction behind the [decimal obstruction](KAPREKAR-METAGRAPH.md). It does not assert existence or prove that the output covers arbitrary integers.
