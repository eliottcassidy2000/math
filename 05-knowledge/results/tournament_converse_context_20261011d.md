# Tournament converse as conditional fusion and a reversible star fan

**Status: PROVED relative to the inherited tournament/even-graph type codec:** the all-order converse formula, conditional fusion product, and exact star-family identification. **PROVED elementary:** the marked labelled fan involution and its sharp sidecar multiplicity. **FINITE-EXACT:** complete type controls through six vertices. No natural labelled tournament/graph bijection, graph-colouring equivalence, or geometric realization is inferred.

## 1. The retained context

The source is the specific inverse pair $\Phi,\Psi$ in
[the weighted Lyndon codec](tournament_even_bijection_20261011c.md), not an unspecified equinumerosity bijection. Its conventions remain fixed: strong tournament letters are ordered by `(vertex count, canonical mask)`; Lyndon words of each weight are ranked lexicographically; connected automorphism-even graph types are ranked by canonical edge mask. The graph predicate is **automorphism-even**, not Eulerian. The count theorem and the distinction from equivariant natural maps are inherited with their sources and scopes there and in [THM-4531, twisted Mallows–Sloane principle](../../01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md).

The closest mechanism is unique ordered SCC decomposition followed by nonincreasing Lyndon factorization. The hostile is already concrete: $P_3\sqcup K_1$ and $K_{1,3}$ correspond to opposite orders of the singleton and directed triangle. The corrected near miss is to apply a local change after forgetting which ports or components the change acts on. The least-used sidecar is a **retained pair of star leaves**.

The working board is: SCC chronology / converse / Lyndon factor boundaries / graph components / retained leaf ports. The [marked four-tournament section](partial_four_tournament_20261011c.md) proves that tournament converse swaps its two completion bits after reanchoring the path. Complementing both bits is a different operation.

## 2. A terminating, confluent fusion rule

Write $\rho(L)$ for the connected graph atom assigned to a Lyndon word $L$ by the inherited ranking. Work first with an **ordered list** of Lyndon words. The following rewrite is allowed at any adjacent increasing pair:

$$
(\ldots,L,M,\ldots)\ \longrightarrow\ (\ldots,LM,\ldots)
\qquad\text{if }L<M. \tag{1}
$$

The inherited elementary concatenation lemma says that $LM$ is Lyndon. Each step removes one boundary, so every schedule terminates. A terminal list is nonincreasing; it is therefore the unique Lyndon factorization of the unchanged concatenated word. Consequently **all choices of legal fusion order have the same normal form**. A starting list of $k$ blocks needs at most $k-1$ fusions when nonempty; this is a block-count bound, not a claim about the cost of constructing the rank atlas.

For graph types $G,H$, decode their component atoms into nonincreasing lists $F(G),F(H)$. Concatenate the lists, apply (1), and replace the resulting factors by their graph atoms. Denote the result by $G\star H$. This is an explicit operation using component/factor data:

$$
G\star H=\Phi\bigl(\Psi(G)\mathbin\to\Psi(H)\bigr). \tag{2}
$$

Here $\to$ is tournament ordered join. Hence $\star$ is associative with the empty graph as identity. More usefully, there is an exact local seam test. For nonempty $G,H$,

$$
G\star H=G\sqcup H
\quad\Longleftrightarrow\quad
\operatorname{last}F(G)\ge\operatorname{first}F(H). \tag{3}
$$

If the inequality holds, the concatenated list is already a normal form and no fusion occurs. If it fails, at least one fusion is compulsory, and the resulting graph has strictly fewer connected components than $G\sqcup H$. Thus equality of graph types is impossible. If either input is empty, equality holds without a seam condition. The comparison in (3) requires the ranked atom words; graph vertex count alone does not supply it.

## 3. Converse and its exact fixed-point criterion

Let $\tau$ send a strong tournament type to its converse. It is a weight-preserving involution on the alphabet. For an SCC word $w=w_1\cdots w_k$, define

$$
w^*=\tau(w_k)\cdots\tau(w_1). \tag{4}
$$

Tournament converse reverses all intercomponent arcs and takes the converse inside every component. Its SCC word is therefore exactly (4).

The transported operation $D=\Phi\circ\mathrm{converse}\circ\Psi$ has a direct normal-form algorithm:

1. Decode and sort the graph's component atoms into $F(G)$.
2. Concatenate their strong letters, reverse their order, and apply $\tau$ to every letter.
3. Treat the resulting letters as singleton Lyndon blocks, apply the fusion rule (1), and re-encode the factors as graph atoms.

This proves, at every finite order,

$$
D^2G=G,\qquad D(G\star H)=D(H)\star D(G). \tag{5}
$$

The identities are on isomorphism types. The code returns canonical representatives, so twice applying it to an arbitrarily labelled graph returns its canonical representative, not necessarily its original labels.

The exact fixed-point criterion is that the decoded SCC word is a $\tau$-palindrome:

$$
D(G)=G\quad\Longleftrightarrow\quad
w_i=\tau(w_{k+1-i})\quad\text{for every }i. \tag{6}
$$

Necessity uses uniqueness of the ordered SCC-type word, not merely equality of the multiset of SCC sizes. Let $s_n$ count strong types, $d_n$ the self-converse strong types, and let $T(z)=1/(1-\sum s_nz^n)$ be the full tournament-type series. A palindrome has an arbitrary left half, its forced right half, and either no central letter or a single self-converse central letter. Thus the number of fixed graph types is recorded by

$$
T(z^2)\left(1+\sum_{n\ge1}d_nz^n\right). \tag{7}
$$

This is a decomposition proof, not an inference from the small census.

Two losses are decisive. First, $D$ does **not** act independently on graph connected components: it fixes $P_3$ and $K_1$ separately, but sends their union to $K_{1,3}$. Second, reversing the order of the SCC letters while leaving each letter untouched is insufficient: through order five all strong letters happen to be self-converse, but at order six canonical tournament masks $18$ and $22$ are a converse pair. The internal involution $\tau$ must be retained.

## 4. The star/path pair is an all-order literal graph operation

Let $a$ be the singleton tournament type and $b$ the directed triangle. Both are self-converse and $a<b$. For every $r\ge0$, the word $a^r b$ is Lyndon: every proper suffix meets $b$ earlier. It is in fact the **first Lyndon word of weight $r+3$** in the specified ordering. A word consisting only of $a$'s is not Lyndon at that weight, and there is no strong letter of weight two. Hence any other Lyndon word has fewer than $r$ initial $a$'s and is lexicographically larger.

The first connected automorphism-even graph of order $n\ge3$ in canonical-mask order is the star $K_{1,n-1}$. Its mask is $2^{n-1}-1$, using the initial pair positions $01,02,\ldots,0(n-1)$. Any smaller mask has fewer than $n-1$ edges and cannot be connected. The star is automorphism-even: orient all edges from its centre to its leaves; every automorphism preserves those directions.

The rank-zero atoms therefore match exactly, and the reverse word factors as $b\mid a\mid\cdots\mid a$:

$$
\Phi(a^r b)=K_{1,r+2},\qquad
\Phi(ba^r)=P_3\sqcup rK_1,
\qquad
\boxed{D(K_{1,r+2})=P_3\sqcup rK_1.} \tag{8}
$$

The first equality depends on the specified atom ordering; the factorization does not. This strengthens the earlier conclusion that $\Phi(a^r b)$ is merely connected. At $r=0$ the path is fixed; at $r=1$ these are precisely the three-leaf star and the path plus an isolate.

For $n=r+3\ge4$, retain a frame $(V,c,L)$ consisting of the actual vertex set, the centre $c$, and an unordered two-element subset $L\subset V\setminus\{c\}$. The two graphs in this frame have edge sets

$$
E_0=\{c\ell:\ell\in L\},\qquad
E_1=\{cv:v\ne c\}. \tag{9}
$$

Toggle every edge $cp$ with $p\notin L\cup\{c\}$. This is a literal labelled involution $E_0\leftrightarrow E_1$, covariant under simultaneous relabelling of graph and frame. It realizes (8) after forgetting the marking. The path-plus-isolates graph determines its centre and retained pair intrinsically. The star determines its centre but not the pair.

This information loss is sharp. On a fixed labelled vertex set there are

$$
n\binom{n-1}{2}\text{ copies of }P_3\sqcup(n-3)K_1,
\qquad n\text{ stars}. \tag{10}
$$

The deterministic completion which connects the path centre to every isolate has exactly $\binom{n-1}{2}$ preimages over each star. Any reversible extension of this completion must distinguish all those preimages; the retained pair does so with exactly that many sidecar values. This is a cardinality statement, not an assertion that a fixed-length binary sidecar uses a noninteger number of bits.

There is also no equivariant way to choose an unmarked pair on a star: its leaf automorphism group acts transitively on all two-element leaf subsets, with no fixed subset when $n\ge4$. At order four, the path-plus-isolate stabilizer has order two and the star stabilizer has order six; the labelled counts are $12$ versus $4$. The missing three-way choice is therefore a structural asymmetry, not a bug in canonical labelling.

At order four, marking the retained pair is equivalent to marking the single leaf whose edge will be deleted. Thus the user's three-leaf star and path-plus-isolate are literally two sides of one edge toggle once that leaf port is retained.

The marked operation needs exactly $r$ edge additions in the completing direction and exactly $r$ deletions in the other direction. No edge-edit procedure on the same vertices can use fewer than $r$ additions to connect $P_3\sqcup rK_1$, since each edge addition reduces the component count by at most one. This is not a lower bound on the runtime of a type-level encoding. Nor is the $n=4$ marked completion a labelled tournament/graph bijection: the four tournament types have labelled orbit sizes $24,8,8,24$, different from the graph-side orbit sizes.

## 5. What the interfaces preserve

| Source and target | Actual map | Preserved predicate | Lost information and required sidecar | Cheapest hostile |
|---|---|---|---|---|
| Ordered Lyndon lists to graph component lists | Conditional fusion (1) | Concatenated SCC word and vertex weight | Input boundaries disappear; retain a factor list when boundaries themselves matter | $ab$ fuses, $ba$ does not |
| Graph types under tournament converse | Reverse/converse letters, then normalize | Involution and reversed ordered product | Ordinary graph components are not independent operations | $D(P_3\sqcup K_1)=K_{1,3}$ |
| Labelled path-plus-isolates to a star | Add centre-to-isolate edges | Centre and whole vertex set | Which two leaves were already present; retain their pair | Three unmarked inverse choices on $K_{1,3}$ |
| Marked fan to its opposite side | Toggle the edges in (9) | Exact inverse and relabelling covariance | No loss while the frame is retained | A different retained pair rejects the same source path |

No colouring or unit-distance conclusion is needed for these proofs. A star and a path with isolates both happen to be bipartite; that coincidence alone would not authenticate the inverse frame, the general operation $D$, or an external gluing. The useful next target is to specify which leaf ports an external context can see, then determine whether a smaller target-specific quotient still suffices. That question is separate from reconstructing the entire graph.

There is already a concrete colouring failure at order six. In the inherited pair-mask convention, $D$ sends mask $110$ to $123$. The first has edges $02,03,04,12,13$ and proper two-colouring $(0,0,1,1,1,0)$. The second has edges $01,02,04,05,12,13$: it contains triangle $012$ and has proper three-colouring $(0,1,2,0,1,1)$. Thus their chromatic numbers are exactly two and three. The selected type involution does not generally preserve bipartiteness or chromatic number, even though the entire fan subfamily does.

## 6. Exact finite controls

[Script](../../04-computation/experiments/tournament_converse_context_20261011d.py) and [saved output](tournament_converse_context_20261011d.out). The all-order theorem uses the inherited finite-enumeration recipe; the shipped API deliberately restricts its rank atlas to order at most six.

| Order | All graph/tournament types | Fixed by $D$ | Self-converse strong letters |
|---:|---:|---:|---:|
| 0 | 1 | 1 | 0 |
| 1 | 1 | 1 | 1 |
| 2 | 1 | 1 | 0 |
| 3 | 2 | 2 | 1 |
| 4 | 4 | 2 | 1 |
| 5 | 12 | 8 | 6 |
| 6 | 56 | 12 | 9 |

Every type is checked against literal reversal of every arc, the inverse identity, and (6). The suite checks all 97 Lyndon-block cuts of the decoded words through six vertices and all possible legal rewrite branches; 209 ordered graph pairs and 435 triples test the seam criterion, reversed product and associativity. All labelled marked fans at orders four through six are tested, including all 47,088 frame relabellings. The sidecar multiplicities are $3,6,10$. Boolean/float aliases and wrong-source frame substitutions are rejected before equality or cached atlas access.

```powershell
python -B -X utf8 04-computation/experiments/tournament_converse_context_20261011d.py
python -B -O -X utf8 04-computation/experiments/tournament_converse_context_20261011d.py
```

Normal, optimized and saved output agree at **48,711 exact checks**; LF SHA256 `98d53da845804c83676a20177310c864f49eb49a21af683dfbe84ac61a21ac52`. No finite census is used to extrapolate a new natural bijection, a graph chromatic invariant, or a plane-colouring result.
