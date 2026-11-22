# Which sets of arc flips can Hamiltonian paths detect?

**PROVED.** A set of arc coordinates has identically zero mixed flip difference precisely when its undirected graph contains a three-edge star, a cycle, or a perfect matching of the whole vertex set. The perfect-matching obstruction gives an exact sampler for Hamiltonian-path-marked tournaments with prescribed row parities. Adjoining one vertex makes every formerly perfect matching detectable and gives its exact mean-square response.

## Inheritance and objects

The prior one-arc potential describes individual steps of the metagraph. The even-path-forest expansion identifies the nonzero Fourier coefficients of that potential. Here the question is which sets of steps can interact after the *other* arcs have been assigned arbitrary orientations. In the quotient by vertex relabeling the experiment retains the labels of the flipped pairs, transports them together, and counts corners with multiplicity when classes coincide. The single-arc quotient alone forgets this information. [THM-5271, flip potential](../../../01-canon/theorems/THM-5271-the-count-is-a-potential-on-the-flip-metagraph-the-seam-ladder-implies-chan-pak-the-microcosm-3-1-6-1-21-the-square-sum-cycle-at-32-and-the-order-10-near-join-reach.md) · [THM-4761, signed forest coefficients](../../../01-canon/theorems/THM-4761-random-tournaments-the-forest-identity-is-the-walsh-chaos-decomposition-of-kappa.md)

The four objects used below are the set of flipped pairs, the forest carrying a coefficient, the matching rows of a complete graph, and a Hamiltonian path marked on a tournament. A perfect matching by itself has zero Fourier coefficient at every order. Its mixed difference has a different boundary, determined by what *larger* forests can contain it. A second possible error is to infer the support of a count from its first moment. Two disjoint supports at order four give a hostile example below.

Let $n\ge2$, $E=E(K_n)$, $d=|E|$. Choose a reference direction on each pair. A tournament is a sign assignment $x\in\{\pm1\}^E$. Write $H(x)$ for its number of directed Hamiltonian paths and $\chi_F(x)=\prod_{e\in F}x_e$. For $D\subseteq E$ fix outside signs $y$ and define the unnormalised mixed difference
$$
\Delta_D H(y)=\sum_{x_D\in\{\pm1\}^D}\chi_D(x_D) H(x_D,y).
$$
Changing a reference direction changes at most the global sign of the corresponding difference after reindexing. A forest is *even* here when each nontrivial path component has an even number of edges. Isolated vertices have length zero.

## 1. The universal equations on flip faces

**Theorem 1.** For every $D\subseteq E$, the following are equivalent.

1. Some choice of outside orientations has $\Delta_D H\ne0$.
2. $D$ is contained in an even linear forest of $K_n$.
3. $D$ is acyclic, has maximum degree at most two, and contains no perfect matching of $K_n$.

Thus all universally vanishing individual mixed differences are generated under inclusion by three-edge stars, cycles and spanning perfect matchings. The list of generators can contain redundancies, for example an even Hamiltonian cycle contains a spanning perfect matching.

**Proof.** Write $H=\sum_F\widehat H(F)\chi_F$. The expansion in THM-4761 gives exactly the even linear forests as the nonzero coefficients, including the empty forest. On a parallel $D$-face, character orthogonality gives
$$
\Delta_D H(y)=2^{|D|}\sum_{F\supseteq D}\widehat H(F)\chi_{F\setminus D}(y).
$$
The outside characters in the sum are distinct, giving the first equivalence.

Every component of an even linear forest has odd order, so a spanning perfect matching contained in it would have to match an odd number of vertices within each component. Conversely let $D$ be a linear forest on all $n$ vertices, isolated vertices included. A path has a perfect matching exactly when it has even order. Hence the absence of a spanning perfect matching supplies a component of odd order. Join one of its endpoints to an endpoint of any even-order component. The combined path still has odd order, and the necessary edge exists in $K_n$. Repeat until the remaining components all have odd order. The resulting even linear forest contains $D$. ∎

**Two-level consequence.** Whenever $D$ is a nonempty forbidden set, the two parity halves of every $D$-face have the same total $H$. If only two counts $h<k$ occur on that face, the number of corners at height $k$ agrees in the two halves, by subtracting their common baseline $2^{|D|-1}h$. In particular the larger count cannot occur just once. Both the fixed outside orientations and the restriction to two levels matter for this multiplicity conclusion.

## 2. The response one vertex past the boundary

**Theorem 2.** Let $M$ consist of $m\ge1$ disjoint pairs. At order $2m$, $\Delta_M H$ vanishes everywhere. At order $n=2m+1$, with uniform independent outside signs,
$$
\mathbb E_y (\Delta_M H(y))^2=2^{1-m}(m+1)!.
$$

**Proof.** The even-order statement follows from Theorem 1. At odd order, an even linear forest containing $M$ must carry the unmatched vertex in *each* component, since every other vertex is paired within its component. It is therefore a single spanning path. Contract the $m$ matched pairs and order the resulting $m+1$ blocks. The internal directions give $2^m$ choices and reversal counts each undirected path twice. There are $2^{m-1}(m+1)!$ completions.

For such a path $F$, the normalized coefficient of $\kappa=2^{n-1}H/n!$ has magnitude $2/(n)_{n-1}=2/n!$ by THM-4761. Thus $|\widehat H(F)|=2^{1-2m}$. Squaring the mixed-difference expansion and averaging kills all cross terms. Its value is $2^{2m} 2^{m-1}(m+1)! (2^{1-2m})^2$, as claimed. ∎

## 3. Matching rows and exact sampling

Let $n=2m$. Choose $k$ pairwise edge-disjoint perfect matchings $M_1,\ldots,M_k$ in $K_n$. A full one-factorization has $k=n-1$. Write $p_i=\chi_{M_i}$. A *marked tournament* is a pair $(T,P)$ with $P$ an ordered directed Hamiltonian path in $T$. Choosing this pair uniformly gives an unmarked tournament law proportional to $H(T)$.

**Theorem 3.** Under the uniform marked law, the vector $(p_1,\ldots,p_k)$ is uniform on $\{\pm1\}^k$, including after fixing arbitrary orientations outside the union of the rows. Without outside pins, each row-parity sector contains exactly $n!\,2^{d-(n-1)-k}$ marked pairs, interpreting the power as a rational factor if the exponent is negative. A uniform marked pair in any specified sector can be generated using one uniform permutation, at most one label exchange, and independent fair orientation bits.

**Proof of uniformity.** For a nonempty selection of rows, the product of their parities is the character of their union $D$. The theorem on faces annihilates this character against $H$ for every assignment off $D$, since $D$ contains a spanning perfect matching. Sum over any remaining free orientations. All nontrivial parity moments vanish, giving the conditional uniform law by Fourier inversion. The total number of unrestricted marked pairs is $n!2^{d-(n-1)}$, by first choosing the path order.

**Sampler.** Given target parities, perform the following operations.

1. Draw a uniform permutation $P$ of the vertices. Regard consecutive pairs as directed edges. Call a row saturated if all its pairs occur in $P$. The underlying undirected path on $2m$ vertices has a unique perfect matching, so at most one chosen row is saturated.
2. For each row designate one of its pairs in advance. If there is a saturated row and its forced parity is wrong, exchange the labels of its designated pair in the permutation. Its endpoints occupy consecutive positions. Its row stays saturated, exactly one forced row sign reverses, and the remaining chosen rows are still unsaturated.
3. For every unsaturated row reserve one unforced pair. Assign fair independent signs to every remaining unforced pair. Set each reserved sign to realize the desired row parity.

For a feasible final path order let $s$ be the number of saturated rows, either zero or one. It has $2^s$ permutation preimages in steps 1–2. It allows exactly $2^{d-(n-1)-k+s}$ completions in step 3. Dividing the first count by the second shows that each feasible marked pair has probability $2^{-d+(n-1)+k}/n!$. Listing all pairs, chosen rows, and the output signs takes $O(n^2)$ word operations with vertex-label words, plus the uniform permutation and fair bits. This sampler takes no outside pins. ∎

**Odd-order boundary.** On a path through $n=2m+1$ vertices, take the odd-position edges $M$ and the even-position edges $N$. These are edge-disjoint near-perfect matchings. Each parity alone is fair under the marked law, since a nonempty matching itself has zero coefficient. With the reference directions along the path,
$$
\mathbb E_{\rm marked}[p_Mp_N] = \widehat\kappa(M\cup N)=2/n!.
$$
The two row parities become correlated because their union is a single even-length spanning path. This is an explicit failure of the near-perfect version of the independence claim.

**Loss of count support (FINITE-EXACT).** Label the four vertices $0,1,2,3$ and take the rows $(03,12),(13,02),(23,01)$, with increasing reference directions. Record parity bit one for a negative product. Sector $000$ has four tournaments of count $1$ and four of count $5$. Sector $010$ has eight tournaments of count $3$. Both have marked mass $24$. An exact first moment shared by all sectors still permits disjoint sets of represented odd integers. Chan–Pak remains **OPEN**.

**Connection to the supplied Latin-square source.** Its row representation decomposes $K_{n,n}$ into ordered perfect matchings, and its proposal corrects for the multiplicity of auxiliary histories. Here disjoint perfect matchings in $K_{2m}$ supply independent constraints on a *different* set of marked objects. The operative transfer is to count admissible completions separately for each history. Uniqueness of the perfect matching of a path reduces the needed correction to a single two-to-one label exchange. Forgetting the distinguished path weights each tournament by its number of paths. The Latin-square work bound is a claim about its own algorithm and plays no role in the proof. [Pinned Mira Cult writeup and verification scope](https://github.com/Mira-Cult/math/blob/30ddaaad27b48df2002225c3e3f40cffd02c5253/notes/combinatorics/exact-latin-square-sampling/README.md)

## Reproduction

From the repository root run `python3 05-knowledge/results/chanpak_matching_parity_faces_20261122/probe.py` and compare to [probe.out](probe.out). The script obtains path counts by directly enumerating directed orders and all their completions, then tests all flip sets at $2\le n\le6$. It computes the exact odd-order mean squares at 3 and 5, all row-sector masses at even orders 2, 4 and 6, and the sampler's exact output law in every order-four sector. SHA-256 is `c0f143fda438ea1c42bdc130d905d96443b662acd97ad65eb6ff00c29823e36d` for [probe.py](probe.py) and `8c73cf430555a306778c5a10225d9a150481f0706e40d20716e981bc4798c799` for its output.
