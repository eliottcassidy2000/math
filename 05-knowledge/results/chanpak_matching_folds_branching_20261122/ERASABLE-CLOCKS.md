# Erasable clocks and exponentially expensive clocks

**PROVED (scoped).** A fixed tournament on a distinguished \(k\)-vertex core is joined to a transitive spine. A spine letter records the \(k\) arc directions between its vertex and the core. The core is allowed to be any tournament; it need not be a strong component. Write \(H(W)\) for the number of Hamiltonian paths after installing the attachment word \(W\).

## 1. Long repeat blocks are erasable or exponentially expensive

Consider any consecutive occurrence of \(W^q\) in the spine, with \(L=|W|>0\) and \(q\geq2\), allowing arbitrary exterior parts on both sides. One of the following holds:

1. deleting the whole repeat block leaves \(H\) unchanged;
2. \(H\geq 1+2^{(q-2)L}\).

More explicitly: if any core vertex has both signs on the letters of \(W\), case 2 holds. Otherwise let \(A\) be all exterior vertices that beat the entire block and \(B\) all exterior vertices that lose to it. An arc \(b\to a\) with \(b\in B,a\in A\) implies case 2 with \(H\geq 1+2^{qL}\). If there is no such arc, the block is erasable and can in fact be replaced by a block of the same letters at any length, including length zero. The two alternatives are not logically disjoint as statements about the numerical value of \(H\) alone.

**Proof.** For a mixed-sign vertex \(v\), choose \(x\) in the first copy and \(y\) in the last copy so that \(v\to x\) and \(y\to v\). At least \((q-2)L\) spine vertices lie strictly between \(x,y\). For any even-sized subset of these intermediate vertices, its increasing order along \(x\cdots y\to v\to x\) gives a distinct directed odd cycle. With \(m>0\) intermediate vertices this yields \(2^{m-1}\) odd cycles; for \(m=0\) it yields one. The \(1+2\,(\#\text{odd cycles})\) lower bound follows from the positive [Grinberg–Stanley odd-cycle-packing formula](https://arxiv.org/abs/2307.05569); no disjointness between these individual cycles is required.

In the uniform case, spine vertices outside the block also have uniform arcs to it, by transitivity. If \(b\to a\) exists, the induced tournament ordered \(a,\text{block},b\) is transitive except for \(b\to a\); even allowing only the odd cycles \(a\cdots b\to a\), there are \(2^{qL-1}\) (choose an odd subset of the block). The bound in this case is sharp: with just \(a,b\) outside, all odd cycles use that closing arc and cannot be packed together. If no cross arc \(B\to A\) exists, the entire tournament is the ordinal join \(T[A]\Rightarrow W^q\Rightarrow T[B]\). Before and after erasure the count is \(H(T[A])H(T[B])\), since a forward join concatenates paths and the transitive block has just one path. Empty pieces have count one.

For a single repeated block with fixed exterior, the all-uniform/no-back-arc case gives a constant count for every \(q\), and all other cases grow at least exponentially. Thus an unmarked count that is polynomially bounded as a function of \(q\) must be constant. This is not a claim about rooted or marked observers.

## 2. Boundedly many periodic clocks cannot have a dense spectrum

Fix \(U_0,\ldots,U_r\) and nonempty \(W_1,\ldots,W_r\). Allow all independent nonnegative exponents in

\[
U_0W_1^{q_1}U_1\cdots W_r^{q_r}U_r .
\]

The number of distinct counts at most \(X\) is \(O((1+\log X)^r)\). The same conclusion holds for a finite union, with the largest \(r\).

Indeed for \(X\geq1\), in any realization of count at most \(X\) erase every block with exponent exceeding \(2+\lfloor\log_2 X\rfloor\), one block at a time. Each time the resulting graph is still in the family and has the same count. The exponential alternative is impossible for every remaining overly long exponent. At most \((3+\lfloor\log_2 X\rfloor)^r\) exponent tuples are needed at the end. The activation of a block can depend on other blocks; no independence assumption is made. In particular the erasures stay inside the specified word expression. Simply deleting all singleton strong components would not by itself give this same-language normal form.

The exponent \(r\) cannot be improved uniformly. Use \(r\) core vertices and form an ordinal join of \(r\) near-transitive blocks, each with its first vertex in the core and containing \(q_i\) intervening spine vertices and a last spine vertex. Reverse just the extreme arc in each block. Its two attachment letters are fixed, so the words are of the displayed form (\(W_i\) is a single letter), and the count is \(\prod_i(1+2^{q_i})\). Take \(q_i\in[4^{i-1}m,2\cdot4^{i-1}m]\) independently. Every \(q_i\) exceeds the sum of the preceding maxima, so the expansion in powers of two is a set of distinct subset sums with no carries. Its least positive exponent recovers \(q_1\); cancel \(1+2^{q_1}\) and repeat. There are \(\Omega_r(m^r)\) distinct counts, all at most \(2^{r+2m\sum_{i=1}^r4^{i-1}}\). This gives the matching lower order \((\log X)^r\) for fixed \(r\).

## 3. A finite-automata consequence: a complete compiler needs branching

Let a regular language of attachment words for a fixed core have at most polynomially many words of each length. Then it yields at most polylogarithmically many distinct counts below \(X\). The same holds for finitely many fixed cores with one such language each.

**Proof.** Trim a deterministic finite automaton for the language to states that are reachable from the initial state and can reach acceptance. If a strongly connected state component has two distinct internal outgoing labelled transitions at a vertex, extend both transitions to loops rooted there. Raise the two loops to suitable powers to make their lengths equal; their first letters are still distinct by determinism. Arbitrary concatenations of \(n\) loop blocks give \(2^n\) distinct words of the same length, all accepted after adding a fixed prefix reaching the vertex and a fixed suffix reaching acceptance. This contradicts polynomial growth. Thus a trimmed state component is a singleton without an internal edge or a simple directed cycle, with exactly one internal labelled transition at each vertex. There are finitely many itineraries through the acyclic component graph, including the choice of each labelled exit and the final state. For a fixed itinerary, each cycle's entry and exit are fixed, so the read words form an expression as in §2: the variable part consists of whole turns and the remaining segment is fixed. Taking the finite union proves the claim.

Thus a finite-state fixed-core compiler for a positive-density part of the odd integers must have a branching recurrent component. Branching is necessary, not sufficient; this does not rule out fixed-core compilers with arbitrary interleaving, variable cores, or spines with backward arcs.

## 4. Branching cannot hide only in invisible padding

For an arbitrary attachment word, prune every spine vertex that is a singleton strongly connected component, retaining all named core vertices. The pruned tournament retains the same count by the forward-join product. Its remaining strongly connected components are unchanged by pruning: paths between two vertices of a component cannot leave it and return without absorbing any exterior vertices used. Each surviving component contains a core vertex, since a component made only of spine vertices lies in a transitive tournament and could only be a singleton. There are therefore at most \(k\) surviving components.

Conversely, for a fixed pruned word, all possible pruned-away vertices occur as singleton runs in the at most \(k+1\) ordinal gaps between consecutive surviving components or at the two ends. The attachment letter in each gap is forced by which core vertices lie before or after it. Its position on the spine is likewise determined by the vertices of each existing component. Thus among words of length at most \(m\), a fixed pruned word has at most \(\binom{m+k+1}{k+1}\) preimages (choose nonnegative run lengths with sum at most \(m\)). A fixed unmarked pruned tournament on \(n\leq m+k\) vertices accounts for at most \(n^k\) such pruned words, since naming its ordered core determines the unique transitive order on the complement.

It follows that if a regular compiler has exponentially many accepted words along a subsequence of lengths, it has exponentially many distinct active pruned tournaments there as well, up to an exponential-rate-neutral polynomial factor. It can still have sparse *counts*, but that would require genuine collisions among active tournaments, not just disguised ordinal padding.

## Mechanism, inheritance, and boundary

- **Inherited proved mechanism:** [THM-4878, graded spine grammars](../../../01-canon/theorems/THM-4878-graded-spine-grammars-cut-state-automaton-hankel-and-seam-ranks-cutwidth-censuses.md) and the [apex recurrences](../chanpak_cross_problem_theorems_20261114/APEX_RECURRENCES.md) already give slot automata, the periodic pole bound at integer bases, and a polylog count for *strong* realizations of fixed ordered block powers. This note extends the image bound to arbitrary realizations by erasing long clocks without leaving the expression. It does not repeat the pole proof.
- **Hostiles and sidecar:** a forward join makes its spine length invisible; a single backward outer arc forces the exact \(1+2^{qL}\) growth; one mixed sign at a core vertex already gives the exponential barrier through a very long odd-cycle family. Scalar count forgetting a long run can be legitimate, but only in the first mechanism. Remember the two shores and the direction of their cross arcs.
- **Bridge and next step:** the source is polymer-like positive odd-cycle packing, the target is regular-language growth, and block erasure is the count-preserving map. It destroys the neutral run length but not the active outside incidence. In the opposite direction stored endpoint and matching information survives a homogeneous block (see [twin-block gluing](EULER-GLUING.md)). For direct coverage, try building a selector that uses recurring genuine choices; making a larger fixed list of periodic templates cannot solve the problem.
