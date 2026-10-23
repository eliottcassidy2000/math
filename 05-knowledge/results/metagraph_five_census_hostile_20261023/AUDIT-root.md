# Parent review of the all-size boundary and actual-mark repair

23 October 2026. **PASS with the stated observer scope.**

The prime-order multiplier lemma is valid: both regular cyclic groups are
Sylow p-subgroups of the target automorphism group, because this group
embeds in S_p and p occurs just once in p!. Aligning them by a target
automorphism leaves a permutation normalizing the regular cycle; this is
affine. Dividing the number of skew connection sets by p-1 therefore gives
a valid lower bound on distinct tournament isomorphism classes.

The census upper bound and lower-order extension identity count induced
subsets, not labelled embeddings. For a rooted type (H,h), each induced
H contributes |Aut(H).h| root assignments. Vertex transitivity distributes
those assignments uniformly over p actual vertices. This checks the
otherwise important root-orbit factor in the counting proof.

The quantitative corollary uses only t_k<=2^(k(k-1)/2). If
k<=c sqrt(log2 p), c<sqrt(2), the logarithm of the number of possible
censuses is O(p^(c^2/2) (log p)^(3/2))=o(p), smaller than the logarithm
of the number of prime-cyclic isomorphism classes. This proves the stated
liminf lower bound along primes for a complete observer in that family.
It is not a complexity bound for arbitrary algorithms or an obstruction
to global reconstruction from additional port capacities.

Reviewed the common-prefix proof: splitting each selected subset into
prefix and core parts retains the necessary rooted/unrooted core census.
The core stays a single color class; prefix vertices retain corresponding
colors. The transitive-six convolution preserves the explicit difference
92. Padded nonisomorphism also follows independently from the unique
nontrivial strong core.

Reviewed `marked_probe.py`: initial ordered marks have distinct colors,
all other vertices the common uniform seed; every refinement preserves the
old color and both directed neighbor multisets. Strict refinement ends
only when the partition no longer changes. Singleton colors permit reading
the entire adjacency matrix in invariant color order. Exhausting all roots
and ordered pairs therefore covers every actual path-endpoint choice, while
not claiming every pair occurs on a Hamiltonian path. The Paley31 control
blocks promotion of this finite repair to a universal one-mark theorem.

The root replays of the independent census, unique-source transitive-six
audit, and marked probe reproduce their frozen receipts. The primitive
word and the marker are distinct information: counting rooted local motifs
does not preserve a chosen identity in the ambient graph.
