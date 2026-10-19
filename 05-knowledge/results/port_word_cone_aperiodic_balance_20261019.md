# A realized port-word cone certificate for aperiodic 7-thin posets

Date: 2026-10-19. Status: proved restricted-class theorem; no claim of solving the general 1/3–2/3 conjecture or of priority for the general transfer method.

## Result first

There is an explicit family of width-three, incomparability-degree-seven posets obtained from arbitrary words in **two non-isomorphic cores** and arbitrary choices from **five two-port seam graphs**. Every block has a named comparison whose probability, under the uniform law on all linear extensions of the whole poset, lies in

\[
[236/525,289/525].                                               \tag{1}
\]

A one-core-buffer version of the same realized cone argument further gives balance at least $5872/12875>9/20$; the small 48-test orbit certificate above is retained as the transparent base theorem.

Subject to a simple endpoint condition, the incomparability graph is connected and has no bridges. The blocks and seams need not be periodic. The first and last blocks supply witnesses among the first and last nine elements, respectively, of **every** reference linear extension.

The proof is a finite boundary-cone certificate realized by legal poset interleavings, rather than a positive-matrix analogy: five shuffle kernels and a permutation symmetry of actual exterior laws give eight left rays and six right rays, and 48 explicit bilinear inequalities per core certify every possible surrounding word. The three-edge seams allow multiple inversions across one interface and genuinely extend the earlier matching-seam construction.

## 1. The cores and permitted seams

Use the original labels throughout. The core \(B_9\) has elements \(0,\ldots,8\) and cover relations

\[
0<2<7<8,\qquad 3<4<5<8,\qquad 0<6,\quad 1<6.
\]

The core \(B_8\) is the induced subposet obtained by deleting element 1. In either core let

\[
V=\{0,3\}\quad\text{(incoming minimal ports)},\qquad
U=\{6,8\}\quad\text{(outgoing maximal ports)}.
\]

The named event is

\[
E=\{5\text{ occurs before }7\}.
\]

Choose any finite nonempty word \(B^{(1)}\cdots B^{(m)}\), each letter being \(B_8\) or \(B_9\). Start with their ordinal sum. At each adjacent seam delete the comparisons in one of the following graphs, where the first coordinate belongs to the earlier block and the second to the later block:

\[
\begin{array}{c|l}
G_0&\{(6,0),(8,3)\}\\
G_1&\{(6,3),(8,0)\}\\
G_2&\{(8,0),(8,3)\}\\
G_3&\{(6,0),(8,0),(8,3)\}\\
G_4&\{(6,3),(8,0),(8,3)\}.
\end{array}                                                     \tag{2}
\]

All other comparisons remain. This is a poset: a deleted comparison goes from a maximal element of one block to a minimal element of the next, and thus cannot be forced by a retained transitive path. Comparisons between nonadjacent blocks are unchanged.

### Theorem 1

Every poset constructed by (2) satisfies (1) for the named comparison in **every** block. Section 5b strengthens its balance to $5872/12875$. It has width exactly three and maximum incomparability degree at most seven.

If the final block is \(B_9\), and every \(B_8\) block has an outgoing seam other than \(G_2\), then the incomparability graph is connected and bridgeless, and its maximum degree is exactly seven. In particular the poset is ordinal-indecomposable. The named comparison in the first block uses elements in the first nine positions of every linear extension, and the named comparison in the last block uses elements in the last nine positions.

There is no restriction on how often either core or any permitted seam occurs. The conditions for bridgelessness are not needed for the probability bound.

## 2. Exact port-word interleaving

For a fixed internal extension \(w\) of a core, write \(\beta(w)\) for its longest initial word consisting only of ports in \(V\), and \(\alpha(w)\) for its longest terminal word consisting only of ports in \(U\). Their possible values lie in

\[
\mathcal B=(\varnothing,0,3,03,30),\qquad
\mathcal A=(\varnothing,6,8,68,86).                              \tag{3}
\]

The core suffix is always nonempty. Its prefix can be empty in \(B_9\). The two port runs are disjoint since \(U\cap V=\varnothing\).

For a seam \(G\), define \(h_G(\alpha,\beta)\) to be the number of shuffles of these two words preserving their internal orders and preserving \(u<v\) whenever \((u,v)\notin G\). This is a count on at most four letters.

### Lemma 2: product of actual seam counts

Fix an internal extension \(w_i\) in every block. The number of global extensions restricting to these internal words is

\[
\prod_{i=1}^{m-1} h_{G_i}(\alpha(w_i),\beta(w_{i+1})).           \tag{4}
\]

**Proof.** A reversed cross-block pair must consist of an outgoing port followed in block order by an incoming port in the next block. If an element of the earlier internal word after that outgoing port were not an outgoing port, it would have to precede every element of the later block, contradicting the inversion. Thus every cross inversion uses the terminal outgoing run. The dual argument restricts its later endpoint to the initial incoming run. There are no inversions between nonadjacent blocks.

Write \(w_i=\beta_i M_i\alpha_i\). Every global extension is uniquely obtained as

\[
\beta_1 M_1\,\operatorname{shuffle}(\alpha_1,\beta_2)\,M_2\,
\operatorname{shuffle}(\alpha_2,\beta_3)\,M_3\cdots M_m\alpha_m.
\]

Only the displayed neighboring port runs can mix. Each chosen shuffle is independently legal precisely when its retained cross comparisons hold. This remains true if a middle word is empty: adjacent seam runs remain ordered by the internal word, by retained comparisons between nonport endpoints, and by all comparisons between nonadjacent blocks. This proves the bijection and (4). \(\square\)

This argument also gives the relevant boundary law. After summing over the whole left exterior, the weight on a core prefix profile is a nonnegative combination of rows of \(h_G\). After summing over the whole right exterior, the weight on its suffix profile is a nonnegative combination of columns of the next kernel. The two sums are independent. Hence a core extension has weight

\[
l_{\beta(w)}r_{\alpha(w)}                                      \tag{5}
\]

for some left and right vectors in these explicitly generated cones. Empty exterior is represented by the all-ones vector. No uniform induced-core marginal is assumed.

## 3. The finite boundary cones

Order profiles as in (3). The five complete shuffle kernels are

\[
K_0=\begin{pmatrix}
1&1&1&1&1\\1&2&1&2&1\\1&1&2&1&2\\1&1&2&1&2\\1&2&1&2&1
\end{pmatrix},\quad
K_1=\begin{pmatrix}
1&1&1&1&1\\1&1&2&1&2\\1&2&1&2&1\\1&2&1&2&1\\1&1&2&1&2
\end{pmatrix},
\]
\[
K_2=\begin{pmatrix}
1&1&1&1&1\\1&1&1&1&1\\1&2&2&3&3\\1&2&2&3&3\\1&1&1&1&1
\end{pmatrix},\quad
K_3=\begin{pmatrix}
1&1&1&1&1\\1&2&1&2&1\\1&2&2&3&3\\1&3&2&5&3\\1&3&1&3&1
\end{pmatrix},
\]
\[
K_4=\begin{pmatrix}
1&1&1&1&1\\1&1&2&1&2\\1&2&2&3&3\\1&2&3&3&5\\1&1&3&1&3
\end{pmatrix}.                                                  \tag{6}
\]

For a direct check of any entry, list the at most six shuffles of a two-letter word with a two-letter word, rejecting those that reverse a cross comparison not deleted in (2).

The distinct row generators are

\[
\begin{array}{c|rrrrr}
i&\multicolumn{5}{c}{l_i}\\\hline
0&1&1&1&1&1\\
1&1&1&2&1&2\\
2&1&1&3&1&3\\
3&1&2&1&2&1\\
4&1&2&2&3&3\\
5&1&2&3&3&5\\
6&1&3&1&3&1\\
7&1&3&2&5&3
\end{array}.
\]

Suppressing the unused empty-suffix coordinate, the distinct column generators are

\[
\begin{array}{c|rrrr}
j&\multicolumn{4}{c}{r_j}\\\hline
0&1&1&1&1\\
1&1&2&2&1\\
2&1&3&3&1\\
3&2&1&1&2\\
4&2&2&3&3\\
5&2&3&5&3
\end{array}.                                                    \tag{7}
\]

Thus every actual boundary pair belongs to the product of these two finitely generated rational cones. We do not claim that every extreme product ray is itself attainable by unconditioned exterior posets.

## 4. The two core count tables

Let \(N^{(s)}_{\beta,\alpha}\) count core extensions with exact profiles \((\beta,\alpha)\), and \(C^{(s)}\) count those also satisfying \(E\). Rows follow \(\mathcal B\); columns are \((6,8,68,86)\). Exact tables are

\[
N^{(8)}=\begin{pmatrix}
0&0&0&0\\0&26&4&4\\0&11&4&4\\0&24&6&6\\0&24&6&6
\end{pmatrix},\qquad
C^{(8)}=\begin{pmatrix}
0&0&0&0\\0&8&1&1\\0&8&3&3\\0&12&3&3\\0&12&3&3
\end{pmatrix};                                                 \tag{8}
\]
\[
N^{(9)}=\begin{pmatrix}
0&85&20&20\\4&90&30&30\\4&69&30&30\\6&60&30&30\\6&60&30&30
\end{pmatrix},\qquad
C^{(9)}=\begin{pmatrix}
0&40&10&10\\1&30&9&9\\3&45&21&21\\3&30&15&15\\3&30&15&15
\end{pmatrix}.                                                 \tag{9}
\]

These have totals \((e(B_8),e_E(B_8))=(125,60)\) and \((e(B_9),e_E(B_9))=(664,325)\).

Here is a short reproducible counting derivation, separate from general-poset enumeration. Shuffle \((0,2,7)\) and \((3,4,5)\) in their \(\binom63=20\) possible ways, then append 8. For \(B_8\), insert 6 in each gap after 0. For \(B_9\), insert 1 and then 6 with \(1<6\) and \(0<6\), equivalently choose their two positions in the final nine-letter word and reject the choices with 6 before 0. These constructions give every extension exactly once. Classify by its first and last port runs to obtain (8)–(9); retain only skeletons with 5 before 7 for the marked table.

As aggregate checks, if \(r\in\{1,2,3,4\}\) is the skeleton position of 0, the numbers of skeletons are \((10,6,3,1)\), of marked skeletons \((4,3,2,1)\), and their insertion multiplicities are \(8-r\) for \(B_8\) and \(36-r(r+1)/2\) for \(B_9\).

## 5. The exact local inequality

By (5), the probability of \(E\) in a block of type \(s\) is

\[
p_E=\frac{lC^{(s)}r}{lN^{(s)}r}.                               \tag{10}
\]

For every pair of rays in (7) and the preceding row list, direct integer multiplication gives

\[
l_i(55C^{(s)}-24N^{(s)})r_j\ge0,\qquad
l_i(31N^{(s)}-55C^{(s)})r_j\ge0.                               \tag{11}
\]

To make the 48 checks per core auditable without long fraction tables, the following vectors give the minimum over the six right rays for each of the eight left rays, in the displayed order:

\[
\begin{array}{c|l}
B_8,\ 55C-24N &(300,740,1180,160,852,1418,0,838)\\
B_8,\ 31N-55C &(575,520,465,1205,1402,1473,1835,2158)\\
B_9,\ 55C-24N &(1939,4138,6337,1379,4460,7100,819,4341)\\
B_9,\ 31N-55C &(2709,2323,1937,5229,5725,5780,7749,8686).
\end{array}                                                     \tag{12}
\]

Every listed denominator \(l_iN^{(s)}r_j\) is positive. Bilinearity extends (11) from generators to arbitrary nonnegative combinations, proving the intermediate bound $24/55\le p_E\le31/55$ for every length and every word of seams. The following actual-law symmetry strengthens it to (1).

For additional precision the exact probability ranges on these **outer product cones** are

\[
B_8:\ [24/55,146/275],\qquad
B_9:\ [1375/3072,265/494].                                     \tag{13}
\]

For example, the first lower endpoint occurs at \((l_6,r_2)\), with marked count 312 and total 715. This establishes sharpness on the stated outer cone, not on the smaller set of all realizable uniform exterior laws. In fact swapping the final two maximal ports, or the initial two minimal ports, imposes additional equal-count constraints on actual exterior distributions.

## 5a. Permutation-orbit cones for the actual exterior law

The preceding individual-profile cones deliberately discard a useful constraint. It has a simple structural proof valid for any number of ports.

### Lemma 3: terminal-antichain orbit equality

Let $U$ be a designated set of maximal elements of a finite poset. Fix $S\subseteq U$. Among uniform linear extensions whose longest terminal run of $U$ consists precisely of $S$, every permutation of $S$ occurs equally often. Dually the same statement holds for initial runs of designated minimal elements.

**Proof.** Permuting those final positions preserves the extension, because the elements of $S$ are maximal and mutually incomparable. It also preserves the longest-run condition and the rest of the word. This gives the required bijections. The initial version is dual. $\square$

Apply the lemma to the **whole exterior**, not just its last or first core. Its designated final ports are globally maximal, and its designated initial ports are globally minimal. Thus terminal-profile weights on $68$ and $86$ are equal, and initial-profile weights on $03$ and $30$ are equal. Valid left cone generators are therefore the two singleton rows of each kernel and the sum of its two doubleton rows, together with the all-ones empty-context vector. Valid right generators are its empty/singleton columns and the sum of its doubleton columns.

After removing repeated rays and dividing by common integer factors, the refined generators are

\[
\begin{array}{c|rrrrr}
i&\multicolumn{5}{c}{\widehat l_i}\\\hline
0&1&1&1&1&1\\1&1&1&2&1&2\\2&1&2&1&2&1\\3&1&2&2&3&3\\
4&2&3&3&3&3\\5&2&3&3&4&4\\6&2&3&6&4&8\\7&2&6&3&8&4
\end{array},\qquad
\begin{array}{c|rrrr}
j&\multicolumn{4}{c}{\widehat r_j}\\\hline
0&1&1&1&1\\1&1&2&2&1\\2&1&3&3&1\\3&2&1&1&2\\4&2&2&3&3\\5&3&6&8&4
\end{array}.                                                     \tag{14}
\]

For these generators, the rowwise minima of the exact deficits are

\[
\begin{array}{c|l}
B_8,\ 525C-236N &(2000,5820,0,5908,6000,7908,20322,3402)\\
B_8,\ 289N-525C &(4625,3720,10155,11158,13875,15783,14022,33327)\\
B_9,\ 525C-236N &(13921,33122,6641,32520,39763,46441,107383,27940)\\
B_9,\ 289N-525C &(21271,15797,43391,44595,59188,65866,52783,135565).
\end{array}                                                     \tag{15}
\]

Bilinearity proves (1) from (15). The sharper one-sided ranges on this refined outer cone are

\[
B_8:\ [236/525,18/35],\qquad B_9:\ [1095/2372,600/1153].         \tag{16}
\]

In particular the common balance constant is $236/525=9/20-1/2100$. The lower extrema use $\widehat l_2=(1,2,1,2,1)$ and $\widehat r_2=(1,3,3,1)$. Even these orbit-summed rays are only an outer description of realizable exteriors; sharpness over actual words is not asserted.

The orbit-sum construction in Lemma 3 applies to arbitrary finite port sets: sum rows or columns over the permutation orbit of each port subset before testing cone inequalities. This can improve a finite transfer certificate without any mixing or asymptotic argument.

## 5b. One whole neighboring core gives an invariant buffer cone

There is a further exact strengthening specific to the stated core-word class. Write $\widetilde K_g$ for $K_g$ with its empty-suffix row removed, so it is $4\times5$. For a central block with a nonempty left exterior, let $s$ be its neighboring core and $x\ge0$ the prefix weights induced by everything further left. Its incoming weight is

\[
l=xN^{(s)}\widetilde K_g.
\]

It therefore lies in the cone generated by the nonzero rows of $N^{(s)}\widetilde K_g$. Union those rays over the two possible core types and five seams, adding the all-ones ray for an empty exterior. Dually use all nonzero columns of $\widetilde K_gN^{(s)}$, and the all-ones ray, for the right cone.

This is an exact, automatically invariant outer-cone construction: a nonnegative row multiplied by any permitted product is a nonnegative combination of that product's rows, and the column statement is dual. It is not an assumed numerical stabilization. It retains every further boundary effect through arbitrary nonnegative $x,y$.

Integer gcd normalization leaves 36 left rays and 15 right rays. The 540 bilinear generator products per central core give the ranges

\[
B_8:\ [5872/12875,58366/115755],\qquad
B_9:\ [759/1622,33199/64887].                                  \tag{17}
\]

Consequently every block has balance at least

\[
\frac{5872}{12875}=\frac9{20}+\frac{313}{51500}.                \tag{18}
\]

For full reproducibility, the rays are defined by (6), (8), (9) and the preceding row/column construction, with no search choices. The exact verifier `verify_port_word_buffer_20261019.py` performs this construction and checks

\[
12875\,lCr-5872\,lNr\ge0,\qquad
7003\,lNr-12875\,lCr\ge0.
\]

Their minimum values over all primitive rays are respectively $(0,102875)$ for $B_8$ and $(285367,465617)$ for $B_9$. All individual integer counts, deficits and ray sources are in `port_word_buffer_cone_20261019.json`; the separate proof is `port_word_one_buffer_20261019.md`. Both lower probability extrema use rays $(17,32,19,32,19)$ and $(5,13,13,5)$. Again these are outer-cone extrema, not a claimed optimum over realizable whole words.

This one-buffer corollary covers empty, one-sided and two-sided exteriors because the all-ones rays are included. It needs only one neighboring block on a nonempty side, not a growing mixing buffer.

## 6. Structural properties and endpoint witnesses

**Width.** Both cores have width three: the displayed chains together with \((1,6)\), or singleton 6 in \(B_8\), are a three-chain cover; \(\{2,4,6\}\) is an antichain. A cross-block antichain uses only two adjacent blocks and only their ports. None of the five deleted-comparison graphs is \(K_{2,2}\), so such an antichain has at most three elements. This proves width exactly three.

**Degree.** In original-label order the internal degrees of \(B_9\) are

\[
(4,7,5,5,5,5,6,5,2).
\]

Those of \(B_8\), in order \((0,2,3,4,5,6,7,8)\), are

\[
(3,4,4,4,4,6,4,1).
\]

Each incoming port gains at most two neighbors, port 6 gains at most one, and port 8 at most two. Incoming and outgoing port sets are disjoint. All resulting degrees are at most seven; the untouched element 1 of any \(B_9\) has degree exactly seven.

**Connectedness and bridges.** The incomparability graph of \(B_9\) contains the spanning cycle

\[
1-8-6-7-3-0-5-2-4-1,
\]

so it is bridgeless. The graph of \(B_8\) has the cycle

\[
0-3-6-7-4-2-5-0,
\]

and its only bridge is \(6-8\). Each seam has at least two edges joining connected core graphs, and every seam edge lies on a cycle through paths in those graphs. A nonstar outgoing seam from \(B_8\) joins both 6 and 8 to the next connected core; together with \(6-8\), these connections give a cycle through that former bridge. Thus the stated endpoint and nonstar conditions remove every bridge. Connectivity is immediate for every permitted word. Connected incomparability is equivalent to absence of a nontrivial ordinal-sum decomposition.

These conditions are useful and honest: the two-block word \(B_8,G_2,B_9\) still has the bridge \(6-8\) in its first block. It satisfies the balance theorem, but not the asserted sufficient condition for bridgelessness.

**Unbounded size and height.** The chain \(2<7<8\) can be concatenated through all blocks, because no seam deletes a comparison into element 2. Thus the family has height at least \(3m\) and arbitrarily many blocks.

**Endpoint control.** Neither 5 nor 7 is an outgoing port, so in the first block both precede every element outside that block. Neither is an incoming port, so in the last block both follow every outside element. Their reference ranks therefore lie within nine of the corresponding endpoint, in every reference extension. The probability guarantee is still computed in the whole poset.

## 7. Exact replay and scope

The standard-library verifier is `verify_port_word_seams_20261019.py`; its saved output is `port_word_seam_controls_20261019.json`.

It independently enumerates the at most six local shuffles for every entry of (6), compares the 20-skeleton insertion construction with explicit core extensions and ideal-DP counts, and checks all 48 exact deficits for each core. Eight fixed assembled posets, including both three-edge seams and the outgoing-star bridge obstruction, agree with unrestricted ideal DP for their total counts and every marked block. A deterministic 48-block mixed word has 409 vertices, degree seven, no bridges, and passes all 48 block probabilities by exact transfer products. No full-poset ideal DP is claimed for that long example.

The refinement is independently replayed by `verify_port_word_refinement_20261019.py`, with all 48 tests per core in `port_word_refined_cone_20261019.json`; `port_word_refinement_20261019.md` records its separate derivation.

The infinite conclusion rests on Lemmas 2–3 and the bilinear certificates (11), (15), not on those finite examples. This is a concrete positive answer for a restricted class to the proposed endpoint-witness and aperiodic-cone route. It does not assert that all 7-thin posets admit these cores, these interfaces, these cones, or a bounded-endpoint witness. The exact terminal/initial port-word kernel is the transferable structural ingredient; the two particular core tables are its present local witness certificate.

Companions: `matching_seam_aperiodic_balance_20261019.md` for the smaller matching-only class and its separate $9/20$ certificate; `finite_cone_balance_route_20261019.md` for the conditional general compiler; `ideal_transfer_positivity_20261019.md` for quantitative ideal-layer contraction. The present proof does not require periodicity, a Perron limit, or an assumed uniform marginal.


Reproduce from the repository root with `python3 05-knowledge/experiments/verify_port_word_seams_20261019.py`. The three standard-library verifiers are colocated for imports; published JSON receipts are in `05-knowledge/results`.
