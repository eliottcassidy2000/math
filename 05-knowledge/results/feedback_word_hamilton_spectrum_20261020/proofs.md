# Proofs: binary cuts, sharp extrema, and the feedback-gap automaton

## 1. Objects, direction conventions, and normalization

Let `v_1,...,v_m` form a transitive tournament, with `v_i -> v_j` iff `i<j`. Add a distinguished vertex `x`. Its binary word is `b_1...b_m`, where

`b_i=1` iff `x -> v_i`.

Write `T(b)` for this tournament and `H(b)` for its ordinary Hamilton-path count. Thus the ordinary tournament order is `n=m+1`.

A leading zero is a source vertex and a trailing one is a sink vertex; deleting either preserves the path count. Repeatedly delete these letters. If no `1` precedes a `0`, the tournament is transitive and `H=1`. Otherwise the reduced word begins with `1` and ends with `0`, and the directed cycle

`x -> v_1 -> ... -> v_m -> x`

shows that the reduced tournament is strong. Conversely, a nontrivial strong tournament in this class must have these endpoint letters. Therefore every nontransitive one-apex tournament has just one nontrivial SCC, itself represented by a reduced word. This observation supplies the all-order quantifier in the obstruction below.

## 2. Exact weighted inversion formula, without a parity black box

Define

`w(1)=1`, and `w(d)=2^(d-2)` for `d>=2`.

**Theorem 1.** For every binary word `b`, including the empty word,

`H(b)=1+2 sum_{i<j} b_i(1-b_j) w(j-i)`.                         (1)

### Direct path proof

Deleting `x` from any Hamilton path leaves an increasing block `A` before `x` and an increasing block `B` after it, with `A,B` partitioning the chain vertices.

- If `A` is empty, exactly one path occurs when `b_1=1`.
- If `B` is empty, exactly one path occurs when `b_m=0`.
- If `max A < min B`, the two blocks form a single ordinary cut. It is valid exactly at an adjacent `01` in the word.
- Otherwise let `i=min B < j=max A`. Validity requires `b_i=1,b_j=0`. Every index below `i` is forced into `A`, every index above `j` into `B`, and each of the `j-i-1` intermediate indices is freely assigned. This gives `2^(j-i-1)` paths.

Hence, for `m>=1`,

`H=b_1+(1-b_m)+N_01+sum_{i<j} b_i(1-b_j)2^(j-i-1)`.

The elementary telescoping identity `N_01-N_10=b_m-b_1` changes the first three terms to `1+N_10`. For adjacent inversions, `N_10` supplies the missing second copy; for all other inversions the displayed power equals `2w(j-i)`. This proves (1). The empty word is the singleton tournament and is immediate.

An alternative derivation observes that all odd cycles contain `x`, counts them by their extreme chain indices, and uses the odd-cycle-collection formula. The direct proof avoids relying on that theorem and makes every Hamilton path explicit.

### A four-coordinate transfer matrix

Initialize `(h,t,l,c)=(1,0,0,1)`. On appending a letter:

- `0`: `(h,t,l,c) -> (h+2t+2l, 2t+l, 0, c)`;
- `1`: `(h,t,l,c) -> (h, 2t+l, c, c)`.

These are fixed `4 x 4` nonnegative integer matrices. The first coordinate after the word is `H`. Here `l` remembers the previous bit and `t` is the geometrically weighted sum of earlier one-bits. Run lengths can therefore be evaluated by exact matrix powers, without enumerating any tournaments or Hamilton paths.

## 3. Strict rearrangement and exact metagraph jumps

Consider two adjacent letters at positions `k,k+1` with `b_kb_(k+1)=10`. Swap them to `01`, leaving every other letter fixed. Let `b'` denote the result. Equation (1) gives

`(H(b)-H(b'))/2`

`= 1 + sum_{i<k} b_i [w(k+1-i)-w(k-i)]`

`    + sum_{j>k+1} (1-b_j)[w(j-k)-w(j-k-1)]`.                 (2)

All bracketed differences are nonnegative, so the swap strictly decreases `H` by at least two. There is no unspecified sign or discarded multiplicity in this identity.

For reduced words, keep positions 1 and `m` fixed and allow only interior swaps, so `2<=k<=m-2`. Define `d(s)=w(s+1)-w(s)`, namely

`d(1)=0`; `d(s)=2^(s-2)` for `s>=2`.

The first and last letters contribute the lower bound

`H(b)-H(b') >= 2[1+d(k-1)+d(m-k-1)]`.                      (3)

If `m>=5`, at least one of `k-1,m-k-1` is at least two. Consequently every such edge changes `H` by at least four.

**Metagraph consequence.** The adjacent-exchange graph of reduced words with fixed numbers of zeroes and ones is connected, but for `n=m+1>=6` none of its edges has ordinary-H increment two. Its connectivity cannot prove a solid odd interval by the unit-step method of THM-4118. This is a theorem about this natural, explicitly defined word metagraph, not about every conceivable tournament metagraph or all single-arc flips.

The move exchanges the orientations of two arcs incident to `x`. It preserves `x`'s score but need not preserve the other scores or regularity. Its vertices are pointed tournaments with a chosen transitive order, not isomorphism classes; forgetting the pointing may identify vertices but does not alter the H-values or create a unit jump on any surviving edge.

## 4. Sharp strong one-apex extrema and an infinite construction

Fix a reduced word of length `m` containing `r` ones, where `1<=r<=m-1`.

Repeatedly apply the decreasing interior swap `10 -> 01`. The unique sorted interior is

`1 0^a 1^b 0`, with `a=m-r-1`, `b=r-1`.

Reversing this sorting gives the maximum word `1^r 0^(m-r)`. The strict inequality in (2) proves uniqueness of these extremal **words for this fixed content**, apart from the fact that the same tournament can have different pointed presentations.

### Minimum evaluation

For `a,b>=1`, equation (1), or the four-coordinate transfer, gives

`H(1 0^a 1^b 0)=(2^a+1)(2^b+1)`.                          (4)

For completeness, partition the inversions into three kinds:

- initial `1` to the `a` early zeroes: weight sum `2^(a-1)`;
- initial `1` to the final zero: weight `2^(a+b-1)`;
- the `b` final one-bits to the final zero: weight sum `2^(b-1)`.

Substitution in (1) gives (4). If `a=0` or `b=0`, the value is instead `2^(m-1)+1`; this boundary exception is important and is not hidden inside (4).

When `m>=4`, an interior choice with `a,b>=1` exists. With `a+b=m-2`, product (4) equals

`2^(m-2)+2^a+2^b+1`.

It is minimized by balanced `a,b`. At `m=4` its value is nine, equal to both boundary cases; for `m>=5` it is strictly smaller than the boundary cases. Therefore, for `n>=5`,

`min H = (2^floor((n-3)/2)+1)(2^ceil((n-3)/2)+1)`.          (5)

For orders three and four the minima are three and five, respectively.

Equation (4) proves an infinite family of **strong** product witnesses, since its word begins with one and ends with zero. The appearance of a product does not imply the tournament is reducible: here it is strong by an explicit spanning cycle.

### Maximum evaluation

For `r,s>=1`, summing the geometric weights gives

`H(1^r 0^s)=(2^r-1)(2^s-1)+2`.                           (6)

Since `r+s=m`, maximizing (6) is equivalent to minimizing `2^r+2^s`, hence balancing `r,s`. Thus

`max H = (2^floor((n-1)/2)-1)(2^ceil((n-1)/2)-1)+2`.        (7)

### First unavoidable failure of the one-apex spectrum grammar

**Corollary.** No one-apex tournament of any order has `H=13`.

Normalize its word as in Section 1. The resulting strong core has order at most five, since (5) is at least 15 from order six onward. Orders three and four give three and five. At order five the four possible reduced words are

`1000,1010,1100,1110`, with counts `9,9,11,9`.

None gives 13. Transitive tournaments give one. This is an all-order exclusion with only four local words, not an enumeration claim. An explicit unrestricted strong witness has out-neighbour masks `(22,12,24,17,2)` on vertices `0,...,4`. It is the two-feedback word `(1,2,0)` with feedback arc `x_0 -> x_1`. The verifier prints its count 13 by both subset DP and literal enumeration of all 120 vertex permutations.

Any multiplicative closure built solely from one-apex tournaments, order joins, and source/sink padding still omits 13, since a product of positive odd integers equal to prime 13 must have a factor equal to 13. The full spectrum conjecture therefore requires a genuinely larger construction grammar.

## 5. Exact finite-state grammar at arbitrary fixed feedback-vertex number

Let `X={x_1,...,x_k}` be a specified feedback vertex set, meaning that the ordinary vertices `v_1,...,v_m` induce a transitive tournament in that order. Fix any tournament `F` on `X`.

Encode `v_i` by a letter `a_i in {0,1}^k`, where `a_i(r)=1` iff `x_r -> v_i`. Every tournament with this marked feedback set is encoded by exactly one such word after the chain order is chosen; conversely every word gives a tournament.

**Theorem 2 (gap-word automaton).** There are explicit nonnegative integer matrices `M_a`, indexed by the `2^k` letters, and fixed initial/final vectors `u,v`, of total dimension at most `k! 3^(k+1)`, such that

`H(T)=u M_(a_1)...M_(a_m) v`.                              (8)

The construction is uniform in `m` and counts ordinary paths, not automorphism orbits.

### State construction

First fix the permutation `sigma` in which a Hamilton path visits `X`. There are `k+1` gaps: before `x_sigma(1)`, between consecutive feedback vertices, and after `x_sigma(k)`.

All ordinary vertices assigned to a particular gap must occur in increasing chain order. For each gap retain one of three states:

- `E`: it is empty;
- `G`: it is nonempty and its current last ordinary vertex points toward the feedback vertex on its right;
- `B`: it is nonempty but its current last vertex does not point toward that right feedback vertex.

For the last gap, which has no right feedback vertex, every nonempty state is `G`. Using three states there anyway gives the stated simple upper bound.

Process ordinary vertices from left to right. Assign the current vertex to exactly one gap.

- If that gap is empty and has a left feedback vertex, require that left feedback vertex to point to the current vertex. This checks the first vertex of the gap once and for all.
- If the gap was already nonempty, its first-vertex condition was checked earlier; no new first condition is imposed.
- Replace its state by `G` or `B` according to the current vertex's orientation toward the right feedback vertex, or by `G` at the last gap.
- Leave all other gap states unchanged.

The transition depends only on the letter, not on the position. Each legal gap choice contributes one. **If different choices produce the same state, their multiplicities are added, not identified.** Thus the letter acts by a fixed nonnegative integer matrix.

Start with all gaps empty. Accept a final state if no gap is `B`, and if every empty internal gap has the necessary direct arc from its left feedback vertex to its right feedback vertex. Empty outer gaps impose no condition.

There are at most `3^(k+1)` states for this permutation. Take the block diagonal union over all `k!` permutations and sum the acceptances.

### Why it counts exactly

Given a Hamilton path, its feedback-vertex order is unique; each ordinary vertex belongs to a unique gap; and the chain's transitivity forces the within-gap order. Its first/last boundary conditions are precisely the transition and acceptance tests above. Conversely, every accepted gap assignment concatenates increasing ordinary blocks and the feedback vertices into one Hamilton path. These constructions are inverse. This proves (8), including paths which a selected insertion matching would miss.

For `k=0` there is one gap and exactly one path. For `k=1` the automaton reduces further to the four-coordinate representation in Section 2. For `k=2`, the uncomplicated bound is 54 states on four letters, even though redundant states can be removed.

### What this resolves, and what it does not

This is a complete **structural representation** for bounded-feedback tournaments. It supplies an exact proof platform for run families, matrix identities, constructive recurrence witnesses, and hostile numerical-image tests. It is not a proof that the numerical image is all allowed odd integers, even for `k=2`.

A sufficient next theorem would be a finite collection of exact matrix-product identities together with a terminating inverse algorithm taking every allowed odd integer to an accepted word. Any such algorithm must supply words for both prime lanes of THM-4094, not merely closure under multiplication. Exhibiting a recurrence or a connected input-word graph without this inverse/existence step is not enough.

The connection contract is explicit:

- source: tournaments with marked feedback set `X` and chosen chain order;
- target: words over `{0,1}^k`, then weighted accepting gap assignments;
- map: record every cross orientation, retain `F`, and partition a path at its feedback vertices;
- preserved predicate: ordinary Hamilton-path count exactly;
- discarded data: only irrelevant vertex labels after a chain order is chosen, never multiplicity;
- sidecar: `F`, feedback permutation, and empty/good/bad gap states;
- cheapest decisive test: literal permutation counting on small encoded tournaments, independently of both transfer and subset DP.

## 6. Historical poisoning-graph audit: a concrete two-cycle

This section documents a new compact witness to a historical failed step, not a new refutation of the exclusions or of the final dichotomy.

In the pre-repair THM-079 Part R.2, the asserted acyclicity of the poisoning graph overlooks directed two-cycles: the proof uses `w_k != w_2`, which fails when `k=2`.

Take an order join of:

- the strong four-tournament on `0,1,2,3`, transitive except `3 -> 0`;
- a directed triangle `B={4,5,6}`, together with vertices `7,8` and arcs `7 -> 8 -> B -> 7`.

The tournament's out-neighbour masks are

`(502,508,504,497,160,192,144,256,112)`.

Its triangles are exactly

`{0,1,3}, {0,2,3}, {4,5,6}, {4,7,8}, {5,7,8}, {6,7,8}`.

Choose `A={0,1,3}` and `B={4,5,6}`, leaving `R={2,7,8}`. Every triangle containing 7 contains 8, and conversely. Thus the poisoning graph has `7 -> 8 -> 7`. The maximum triangle matching is two, and all nine vertices lie on triangles. Its path count is 75.

The actual safe-deletion conclusion is not falsified: deleting vertex 2 works. What fails is the DAG lemma used in that historical argument. The current Moon–Busch/SCC proof of the permanent omissions is unaffected. Incoming repairs are owned by the parent/publisher; this packet stages the hostile without editing shared canon.
