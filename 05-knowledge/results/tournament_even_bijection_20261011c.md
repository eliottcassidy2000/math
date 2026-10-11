# A compositional tournament/even-graph bijection on isomorphism types

**Status:** PROVED relative to the cited all-order equinumerosity theorem: the explicitly specified inverse pair and factor/component invariant. PROVED elementary: the decomposition, factorization argument, and four-type edge-toggle obstruction below. FINITE-EXACT: complete controls through order six. This is an effective type codec with canonical ranking at irreducible atoms, **not** a relabelling-equivariant or rank-free natural bijection, and not a new proof of the cited counting theorem.

## 1. What is inherited, and what is being constructed

The current reference is [THM-4531, twisted Mallows–Sloane principle and scoped natural-map obstructions](../../01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md). The relevant even predicate is

\[
\epsilon_G(\sigma)=(-1)^{\#\{ij\in E(G):i<j,\ \sigma(i)>\sigma(j)\}}=1
\quad(\sigma\in\operatorname{Aut}(G)).
\]

It is independent of the chosen reference orientation. We call this **automorphism-even**. It is not the Eulerian predicate that all vertex degrees are even. A three-vertex path is automorphism-even and not Eulerian; a triangle is Eulerian and not automorphism-even. The [previous transfer packet](tournament_even_transfer_20261011b.md) contains exact small tables and the separate finite Fourier connection to colouring.

The classical input is equality of the numbers of unlabelled tournament and automorphism-even graph types at each order. [Royle–Praeger–Glasby–Freedman–Devillers, Theorem 1.1](https://arxiv.org/html/2204.01947v2). Equal finite cardinalities already allow arbitrary rank bijections. Our aim is narrower than the original natural-bijection problem: retain intrinsic decompositions and use rank choices only to match the remaining irreducible atoms.

The closest constructive mechanism is the [ordered SCC codec and four-state section](tournament_recursive_four_state_20261004.md), inherited from [HYP-3244](../hypotheses/HYP-3244-tiling-half-tiling-interlocking-recursions.md). Its hostile is losing the presentation before performing a named tile flip. The corrected near miss is treating an equinumerosity theorem or a four-element set as a canonical labelled correspondence. The least-used sidecar is the **ordered list of strongly connected tournament types**.

The live board is: **SCC chronology / connected graph components / Lyndon factors / automorphism sign / marked two-bit section**. No plane-colouring, unit-distance, or abc predicate is asserted to survive the new map.

## 2. Two decomposition laws

Every tournament has a unique ordered list of strongly connected components. Between two distinct components all arcs point in one direction, and their condensation is an acyclic tournament, hence a total order. Conversely any finite word of nonempty strongly connected tournament types produces one tournament type by ordered join. The one-vertex tournament counts as strongly connected. The empty word represents the empty tournament.

For ordinary graphs, the corresponding intrinsic decomposition is an unordered multiset of connected graph types. The needed extra fact is exact:

> A graph is automorphism-even if and only if all its connected components are automorphism-even.

Necessity follows by extending a component automorphism as the identity elsewhere. For sufficiency, orient every copy of an isomorphic connected component using the same chosen oriented model. A pure permutation of identical copies preserves every directed model edge and reverses **zero** edges. Every graph automorphism factors into such a copy permutation followed by internal component automorphisms. Their signs multiply to one. There are no intercomponent edges contributing a hidden swap sign. Orientation independence then removes the temporary model choice.

Let $s_n$ count strong tournament types and $c_n$ count connected automorphism-even graph types, and include the empty type in the full counting series. These decompositions give formal power-series identities

\[
T(z)=\frac{1}{1-S(z)},\qquad S(z)=\sum_{n\ge1}s_nz^n,
\qquad E(z)=\prod_{n\ge1}(1-z^n)^{-c_n}.
\tag{1}
\]

Every coefficient involves finitely many types. The cited equinumerosity theorem supplies $T(z)=E(z)$.

## 3. The ordered-word to multiset bridge

Choose a total order on the alphabet of strong tournament types. Our exact convention is to sort by `(vertex count, canonical mask)`, where the mask is the smallest binary encoding under all vertex relabellings. Each letter's weight is its vertex count. A nonempty word $w$ is **Lyndon** if it is lexicographically smaller than every proper nonempty suffix. Equivalently it is strictly smaller than every nontrivial cyclic rotation. This is the classical standard-sequence notion; see [Chen–Fox–Lyndon, §1 and Theorem 1.4](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/fox4.pdf). The following elementary proof records the precise direction we need, rather than relying on a naming convention for factor order.

**Factorization lemma.** Every word has exactly one factorization

\[
w=L_1L_2\cdots L_k,\qquad L_1\ge L_2\ge\cdots\ge L_k,
\tag{2}
\]

into Lyndon words. The empty word has the empty factorization.

First, if Lyndon words $x<y$, then $xy$ is Lyndon and $xy<y$. For $xy<y$, a first differing symbol settles the comparison unless $x$ is a prefix of $y$. In that case $y=xz$ and the Lyndon suffix property gives $y<z$, hence $xy=x y<x z=y$. A proper suffix lying within $y$ is at least $y$. A proper suffix beginning inside $x$ has the form $x'y$ with $x<x'$; neither can be a prefix of the other, since $x'$ is shorter and cannot precede $x$ in lexicographic order. Hence $xy<x'y$. This proves the claim.

Now take the lexicographically least nonempty suffix $v$ of $w$. It is Lyndon, since all its own proper suffixes are suffixes of $w$. Factor the preceding prefix inductively. Its last factor $u$ cannot satisfy $u<v$, for then $uv<v$ would be a smaller suffix of $w$. Thus appending $v$ produces (2). Conversely, in any factorization (2), the last factor is the least suffix: a suffix beginning strictly inside an earlier factor starts with a proper suffix larger than that factor and hence larger than the last factor; one beginning at a factor boundary is at least the last factor, with equality only for the final occurrence alone. So the last factor is forced, and induction proves uniqueness.

This smallest-suffix recursion is the implementation used here. It is finite and deliberately simple. Faster classical factorization algorithms are available; [Duval, *Factorizing words over an ordered alphabet* (1983)](https://doi.org/10.1016/0196-6774(83)90017-2). No speed claim is made for our exhaustive graph atlas.

Let $\ell_n$ count Lyndon words of total **weight** $n$, not necessarily of word length $n$. Equation (2) identifies all words with multisets of Lyndon factors, sorted in decreasing lexicographic order. Therefore

\[
T(z)=\prod_{n\ge1}(1-z^n)^{-\ell_n}.
\tag{3}
\]

Combining (1)–(3) proves $c_n=\ell_n$ for every $n$. The Euler-product exponents are unique: after cancelling all factors of smaller index, the coefficient of $z^n$ is precisely the exponent at index $n$. This argument uses the cited count equality; it is not an independent bijective proof of that equality.

## 4. An explicit inverse pair

At each weight $n$, list the Lyndon words lexicographically and list the connected automorphism-even graph types by increasing canonical edge mask. The two lists have equal finite length. Pair entries of equal rank; denote this atom map by $\rho_n$.

Define $\Phi$ from tournament types to graph types:

1. Read the ordered SCC-type word.
2. Compute its nonincreasing Lyndon factorization.
3. Replace each factor $L$ by $\rho_{\operatorname{wt}(L)}(L)$.
4. Take their disjoint union.

Define $\Psi$ in the other direction:

1. Read the connected-component types of an automorphism-even graph.
2. Replace each component by its inverse atom word.
3. Sort the resulting Lyndon factors in nonincreasing lexicographic order.
4. Concatenate their strong letters and take the ordered tournament join.

**Inverse theorem.** These maps preserve vertex count and satisfy $\Psi\Phi=\mathrm{id}$ and $\Phi\Psi=\mathrm{id}$ on isomorphism types at every finite order. Moreover,

\[
\{\!\{\text{component sizes of }\Phi(T)\}\!\}
=\{\!\{\text{weights of the Lyndon factors of the SCC word of }T\}\!\}.
\tag{4}
\]

The proof is exactly the uniqueness of the three decompositions and the inverse atom ranks. All necessary finite lists can be computed from the input order: enumerate masks and relabellings, test strong connectivity or the automorphism sign, enumerate words of that total weight, then sort. Positive letter weights and finitely many letters below any weight make this a terminating all-order recipe. The shipped implementation intentionally rejects orders greater than six; this is a declared resource boundary, not the theorem's range.

The ranking is explicit and remains part of the construction. Reordering the alphabet or the connected-atom catalogue changes the correspondence. Arbitrary labelled inputs are first interpreted up to isomorphism, and the output is a canonical representative; the output labels do not covary with the input labels. If future operations need particular vertices or tiles, retain a canonicalization transporter and the original marking. The graph type by itself owes no such transport.

### What is not preserved

Ordered join is **not** generally sent to disjoint union. Let $a$ be the one-vertex type and $b$ the directed triangle, with $a<b$. Then $ab$ is one Lyndon factor, whereas $ba$ factors as $b\mid a$. Under our convention

\[
\Phi(a\to b)=K_{1,3},\qquad
\Phi(b\to a)=P_3\sqcup K_1,
\qquad
\Phi(a)\sqcup\Phi(b)=P_3\sqcup K_1.
\tag{5}
\]

This order effect is unbounded. For every $r\ge1$, the word $a^r b$ is Lyndon: a proper suffix has fewer initial $a$'s and meets $b$ sooner, so it is larger. Thus $\Phi(a^r\to b)$ is connected. In contrast, $ba^r$ factors as $b\mid a\mid\cdots\mid a$, and

\[
\Phi(b\to a^r)=P_3\sqcup rK_1.
\]

Moving the same directed triangle across $r$ singleton SCC letters therefore fuses $r+1$ graph components into one. This conclusion is independent of the atom-ranking choices. Under any identification of the output vertex sets, an edge-edit implementation starting from $P_3\sqcup rK_1$ needs at least $r$ edge additions to reach a connected graph: one added edge reduces the component count by at most one, and deletions cannot reduce it. This is a lower bound on that particular implementation of the type transformation, not a computational lower bound for the codec. The script checks $r=1,2,3$, the cases inside the six-vertex atlas.

Nor is the map relabelling-equivariant. At order three it sends the directed cycle to the path $P_3$. A three-cycle automorphism fixes the tournament but cannot fix any labelled path, whose automorphism group has order two. This concrete failure is compatible with THM-4531's separately scoped existence and obstruction results for other maps.

Periodicity is retained: $aa$ is not Lyndon and factors as $a\mid a$, producing two isolated vertices, not a single size-two atom. Sorting all SCC letters before factorization would destroy the distinction in (5). The preserved invariant is (4), not the unordered SCC multiset, edge count, graph chromatic number, or geometric realizability.

## 5. The requested four-type partial tournament

Use the inherited Hamiltonian path $0\to1\to2\to3$, keep $03$ forward, and let bits $\alpha,\beta$ reverse $02,13$. The four states are a transversal of the four tournament types. The present construction gives the explicit inverse table:

| $(\alpha,\beta)$ | Tournament type | Graph type | Canonical graph mask |
|---|---|---|---:|
| $(0,0)$ | transitive | four isolated vertices | 0 |
| $(1,0)$ | directed triangle over a sink | $P_3\sqcup K_1$ | 3 |
| $(0,1)$ | source over a directed triangle | $K_{1,3}$ | 7 |
| $(1,1)$ | strongly connected | $C_4$ | 30 |

Graph masks use pair order $01,02,03,12,13,23$, with a one denoting an edge. The APIs `section_to_graph` and `graph_to_section` implement this marked inverse pair, including arbitrary labelled representatives of each graph type on the inverse side. Its meaning is a specified section and specified atom ordering; it is not an intrinsic preferred correspondence between unmarked objects.

There is a useful obstruction to treating both sides as the same tile algebra. **No choice of labelled representatives of these four graph types is an affine Boolean square under fixed edge toggles.** Their edge counts are $0,2,3,4$, so exactly one representative has odd edge parity, regardless of its vertex labels. On an affine two-dimensional binary subspace, every linear functional is constant or balanced, and therefore has zero, two, or four ones. Edge-count parity cannot have exactly one. This proves the obstruction for every labelling and every pairing of the four types.

The tournament section does have two fixed independent arc flips. Its transported graph moves are well-defined permutations of the four types but necessarily depend on the current state if implemented by edge changes. The script also exhausts all affine-square candidates among the 20 labelled automorphism-even graphs of order four, independently confirming the parity obstruction.

The minimal interaction coordinate is explicit. If $p(\alpha,\beta)$ is the parity of the image's edge count, the table gives

\[
p(\alpha,\beta)=\beta+\alpha\beta=\beta(1+\alpha)
\quad\text{over }\mathbb F_2.
\tag{6}
\]

Every function on two bits has a unique multilinear polynomial. Here the mixed quadratic coefficient is one, so no affine replacement suffices. The discrete port derivatives are $\Delta_\alpha p=\beta$ and $\Delta_\beta p=1+\alpha$: the parity change at one port depends on the other port. This is a precise interaction required by this correspondence, not a general obstruction to nonlinear graph-type encodings.

## 6. Exact controls and limits

[Script](../../04-computation/experiments/tournament_even_bijection_20261011c.py); [output](tournament_even_bijection_20261011c.out).

The complete finite catalogues are:

| Order | Tournament / automorphism-even types | Strong tournament types | Connected automorphism-even types = weighted Lyndon atoms |
|---:|---:|---:|---:|
| 0 | 1 | 0 | 0 |
| 1 | 1 | 1 | 1 |
| 2 | 1 | 0 | 0 |
| 3 | 2 | 1 | 1 |
| 4 | 4 | 1 | 2 |
| 5 | 12 | 6 | 8 |
| 6 | 56 | 35 | 43 |

The checks enumerate all type orbits through six vertices, test both inverse identities and the factor-weight invariant, and compare the factorization with every possible cut of each SCC word. They also test every relabelling of every representative on both sides: 83,744 such controls, including repeated presentations arising from automorphisms. Those tests establish type invariance, not labelled equivariance. The four-section, order-join, periodicity, Eulerian/automorphism-even and malformed-input hostiles are retained. Exact type guards precede cached catalogue lookups, so boolean/float aliases cannot enter the codec.

```powershell
python -B -X utf8 04-computation/experiments/tournament_even_bijection_20261011c.py
python -B -O -X utf8 04-computation/experiments/tournament_even_bijection_20261011c.py
```

Normal, optimized, and saved output agree at **84,251 checks**; LF SHA256 `b7e2afc5e5141b900b03c80a65ece44a42a9f3a720dfe305bc8711ca9e255c23`. The proof above establishes the general recipe; the census validates the bounded implementation. The construction does not give a rank-free atom map, a natural species bijection, a plane-colouring transfer, or an abc consequence. A concrete next improvement would be an atom pairing preserving an independently specified invariant, with equal fibres proved before replacing the rank choice.
