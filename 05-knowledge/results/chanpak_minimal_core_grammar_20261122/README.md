# A target-capped grammar for strong cores and sharp branching on every contact face

**PROVED.** Let a new vertex meet a fixed tournament with \(r\) freely oriented contacts and all other contacts prescribed. Its Hamiltonian-path count takes at least \(\lfloor r/2\rfloor+1\) different values. This bound is attained for every \(r\). For independent fair contacts any one target has probability \(O(1/\sqrt r)\); that exponent is also attained. A positive weighted sum over perfect matchings gives the first result.

A directed minimum cut computes the cheapest *physical* completion of a partially specified mixed attachment. Combining that fact with canonical deletion gives an executable, duplicate-free recursion for every strong isomorphism class below a target count. A separate calculation resolves multiplicity at symmetric parents. For example, there are exactly **1,091** strong order-20 classes of count **11,212,997,148,405** with a Paley-19 deletion. Of these **1,076** have trivial automorphism group. Other deletion families may supply further realizations.

The first nonunique strong fibre already appears at \(H=9\): two rigid, self-converse order-five classes have different vertex-deletion count decks. Scalar multiplicity, graph symmetry and duplication from different choices of a root are three different observables. All orders here are positive; “strong core” means a strongly connected tournament of order at least three. Chan–Pak remains **OPEN**.

## The exact cost of contacts

For a base \(D\), let \(h=H(D)\), and let the new vertex beat exactly \(S\subseteq V(D)\). Define \(w_{ij}\) as follows. On every nonempty even-order subset count paths from \(i\) to \(j\), weighted by the Hamiltonian-path count of the complementary subtournament; sum these integers. Take \(H(\varnothing)=1\). Then the child count is
\[
F(S)=h+2\sum_{i\in S,\,j\notin S}w_{ij}.                                      \tag{1}
\]
In the odd-cycle packing formula a pack contributing to the surplus has a unique cycle through the added vertex. Removing that vertex exposes exactly one even-order directed path with the specified endpoints, and all remaining packings contribute the complementary path count. The vertex subset consisting of \(i,j\) alone shows \(w_{ij}+w_{ji}\ge H(D-\{i,j\})\ge1\). This gives a directed, nonnegative endpoint realization of the inherited quadratic contact observer. [THM-002, odd-cycle packing](../../../01-canon/theorems/THM-002-ocf.md); [THM-5266, quadratic extension grammar](../../../01-canon/theorems/THM-5266-the-extension-grammar-curried-one-vertex-extension-count-is-a-quadratic-form-the-additive-sub-grammar-never-outputs-7-or-21-extension-graph-and-toothpick-tree-to-order-9-first-creators-are-85-percent-with-degree-blind-genealogy.md).

### The face support theorem

Write \(x_i=[i\in S]\) and restrict (1) to any \(r\)-dimensional coordinate face of the Boolean cube. Every remaining pair has quadratic coefficient \(-2(w_{ij}+w_{ji})<0\), unchanged by the pinned coordinates. Suppose the only outputs were \(a_1,\ldots,a_k\) with \(2k\le r\). The multilinear reduction of \(\prod_{t=1}^k(F-a_t)\) would vanish on the whole cube and hence be the zero multilinear polynomial, by Boolean Möbius inversion. Choose any \(2k\) free coordinates \(U\). Their squarefree monomial of degree \(2k\) has coefficient
\[
(-2)^k k!\sum_{M\text{ perfect matching of }U}\prod_{\{i,j\}\in M}(w_{ij}+w_{ji}).
\]
Every factor must contribute two new variables; every matching appears with all \(k!\) orders of its edges. Lower degrees and repeated-variable monomials cannot contribute. Every summand has the same strict sign. This contradiction proves the bound.

In fact a list of \(k\) forbidden counts can be escaped using **any specified \(2k\) contacts**, with the rest set arbitrarily. At least a \(4^{-k}\) fraction of the completions escape: every choice on the other free contacts has an escaping extension. For a strong \(n\)-vertex parent every mixed contact is strong, and the only masks with \(F=h\) are the two constant ones, since every nontrivial directed cut has an arc leaving it. Thus there are at least \(\lfloor n/2\rfloor\) distinct strictly larger strong child counts. Iteration from the directed triangle gives at least \(\prod_{j=3}^{n-1}\lfloor j/2\rfloor\) distinct strictly increasing **count histories** ending at order \(n\). Distinct histories may have the same final count or unmarked class.

### A concentration bound for a target and its history

For \(r\ge2\), every restriction leaving two contacts unfixed still has a nonzero mixed coefficient. Thus it takes at least \(r-1\) pins to force the scalar to be constant. Converting Boolean variables to independent Rademacher variables preserves this property. Applying Theorem 1.1 of Kwan and Sauermann to \(F-N\) gives the absolute-constant bound \(\Pr(F=N)\le C/\sqrt r\) on every face. Exposing one vertex of a uniform labeled tournament gives the same bound for its unconditional count with \(r=n-1\). Equivalently the total isomorphism mass \(\sum_{[T]:|T|=n,H=N}1/|\operatorname{Aut}T|\) is at most \(O(2^{\binom n2}/(n!\sqrt n))\). For a sequence of extensions whose next target is chosen from the previous history before its fair free contacts are drawn, successive conditional bounds multiply. With \(r_j\ge2\) free contacts at the specified steps, the chance of following every prescribed height is at most \(\prod_j\min(1,C/\sqrt{r_j})\). [Kwan and Sauermann, *Resolution of the quadratic Littlewood–Offord problem*, Theorem 1.1](https://arxiv.org/pdf/2312.13826).

For sharpness take the arc-transitive Paley parent on a prime \(p=3\pmod4\), \(p\ge r\). Its count depends strictly on the distinct values of \(s(p-s)\), where \(s\) is the number of wins by the apex. This is the existing arc-transitive attachment law. Pin \(\lfloor(p-r)/2\rfloor\) wins outside the free coordinates, leaving exactly \(\lfloor r/2\rfloor+1\) distinct values in the centered run of weights. One of those values includes a middle binomial layer and has probability \(\Omega(1/\sqrt r)\). One may take \(p\ge r+3\) to guarantee mixed pinned contacts and strong children throughout. [Arc-transitive attachment law, T4](../uniform_prime_transport_20261103/THEORY.md#t4-arc-transitive-one-vertex-attachments-and-a-sharp-paley-escape).

## An isomorphism recursion that prunes at the target ceiling

Fix a ceiling \(B\) and a strong parent \(D\) of order \(n\). Once a nonempty set of contacts is pinned on each side, contract all pinned wins to a source and all pinned losses to a sink in the network of capacities \(w_{ij}\). The minimum child count under these pins is **exactly** \(h+2\lambda\), where \(\lambda\) is the directed minimum source–sink cut. A minimizing cut specifies literal orientations for every free contact. This is an exact existence test for a child with \(H\le B\). Gaps below the ceiling still require an equality check.

Enumerate mixed masks by their first bit and position of the first opposite bit; these \(2(n-1)\) disjoint starting sectors have both sides pinned. Branch on the next bit, rejecting a prefix as soon as its minimum exceeds \(B\). With \(L\) labeled qualifying masks this takes \(O(n(1+L))\) max-flow calls for a given parent once its capacities are available: every feasible prefix is on a depth-\(n\) path to at least one output and every failed child is tested just once. This cost concerns the mask listing, before isomorphism processing.

To remove multiplicity, retain the least mask under the full parent automorphism group. Canonically label the resulting child, retaining *all* minimizing labelings. In one canonical labeling take the first vertex whose deletion is strong. Its position defines a unique vertex orbit transported by all minimizing labelings. Retain the child exactly when the new vertex belongs to that orbit. Every strong child of order at least four has a strong deletion by Moon's theorem. Its designated deletion gives a parent under the ceiling because (1) is monotone in the transition from parent to child. Conversely, two accepted presentations of an isomorphic child have equivalent designated roots. An isomorphism fixing that root then identifies the two masks under a parent automorphism. This proves existence and uniqueness by induction from \(C_3\). The code stops when the next exact Busch strong floor exceeds \(B\); the floors at 3,4,5 are 3,5,9, then multiply by five per three additional vertices. [Moon, *On Subtournaments of a Tournament*](https://doi.org/10.4153/CMB-1966-038-7); [Busch, *A note on the number of Hamiltonian paths in strong tournaments*](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3); [inherited exact canonical labeling and all minimizing transports](../rank_geometry_metagraph_20261026/enumeration/README.md).

The class-level binary relation in this grammar is strongly connected one-vertex extension. The arc capacities are directional; an arc and its reverse can have different weights. Forgetting the marked vertex collapses a whole vertex orbit, and retaining only the numerical minimum loses the construction cut. [grammar.py](grammar.py) carries both sidecars. Its `target` argument prints all matching canonical outgoing-mask lists over all orders. The enumeration cost can be considerable for large targets.

### Small multiplicities and the first internal minimum cut

The following all-order strong-fibre counts are inherited from the exact low atlas; the target-capped recursion additionally returns each graph itself. For \(H=3,5,9,11,13,15,17,19,23,33\) the counts of classes are respectively \(1,1,2,1,1,4,2,2,4,8\). At 33 the split is five classes on six vertices and three on seven. For 7 and 21 the output is empty. The two strong classes at 9 both have trivial automorphism group and are self-converse. Their complete single-vertex deletion count multisets are \((1,1,5,5,5)\) and \((1,3,3,5,5)\); these distinguish the two structural sources of the same scalar. [THM-4810, strong-fibre and small-spectrum atlas](../../../01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md).

**FINITE-EXACT.** Over all strong parents at each indicated order, the minimum and maximum below range over different parents; the maximum counts *labeled contact masks*, and includes the two constant masks.

| Parent order | Strong parents | Minimum number of distinct child counts, including parent count | Maximum masks at one count | Parents whose cheapest strong child requires at least two wins and two losses |
|---:|---:|---:|---:|---:|
| 3 | 1 | 2 | 6 | 0 |
| 4 | 1 | 5 | 5 | 0 |
| 5 | 6 | 4 | 12 | 0 |
| 6 | 35 | 10 | 19 | 1 |
| 7 | 353 | 4 | 70 | 5 |

The first failure of singleton probing is the order-six parent with rows \((2,4,9,35,15,23)\), \(H=17\). The unique cheapest contact is \(S=\{0,1,2\}\), giving 27. The best singleton or complement gives 29. This is why contracting the actual pinned sides matters in the inverse.

## Exact stabilizer spectrum at a single prescribed Paley extension count

Let \(p=2m+1\equiv3\pmod4\) be prime, \(D=P_p\), and \(1\le k\le m\). By T4 the two weight layers \(k,p-k\) have precisely one common count, denoted \(N_{p,k}\). It equals \(H(D)+C_2(D)k(p-k)/(p(p-1))\) where \(C_2(D)\) counts ordered covers by two nonempty paths. Values at different \(k\) are distinct. The following formula is **exact for isomorphism classes with a Paley deletion**, giving lower bounds in the unrestricted strong fibre. The children are strong since their contacts are mixed.

For each divisor \(d>1\) of \(m\), define in descending order
\[
M_d=p\sum_{e=0}^1 [d\mid k-e]{(p-1)/d\choose(k-e)/d}
        -\sum_{d\mid t,\ t>d,\ t\mid m}M_t,
\qquad
M_1={p\choose k}-\sum_{d>1}M_d.
\]
Then the number of different children of count \(N_{p,k}\) with automorphism group of order \(d\mid m\) is
\[
\frac{2dM_d}{pm}-[k=m][d=m].                                               \tag{2}
\]
For \(m=1\) the same formula gives a single child. Every nontrivial stabilizer in (2) is cyclic.

Indeed the automorphisms of the Paley prime parent are exactly \(x\mapsto ax+b\) for nonzero square \(a\); the group has order \(pm\). A proper, nonempty subset cannot be stabilized by a nonidentity translation. Its stabilizer therefore injects into the cyclic multiplier group. A nontrivial stabilizer has a unique fixed point: its generator is a nontranslation dilation. For each of the \(p\) possible fixed points, its unique order-\(d\) subgroup has a singleton orbit at the center and \((p-1)/d\) other orbits of size \(d\). A mask on which this subgroup acts takes zero or one center and a collection of full size-\(d\) orbits. Uniqueness of the center and subtracting the larger stabilizers give exactly \(M_d\) masks of weight \(k\). Orbit-stabilizer counts their orbits. Complementing the masks doubles this count across the two layers. The regular-parent theorem says the only repeated unmarked child is a clone, and counts every automorphism of a child as a mask stabilizer. The one orbit of clones occurs exactly at the middle layers with point stabilizer order \(m\). This proves (2). [Goldberg, *The Group of the Quadratic Residue Tournament*, Theorem 2](https://doi.org/10.4153/CMB-1970-010-8); [regular-parent theorem R1–R2](../regular_parent_symmetry_20261029/PROOFS.md).

At this parent, converse interchanges the two mask-weight layers, because \(P_p\) is self-converse. An unmarked child can therefore be self-converse only if it is a clone in the middle layer. Conversely this clone is self-converse: any anti-automorphism of the base sending the cloned vertex to itself (for instance \(x\mapsto2v-x\)) extends by interchanging its two copies. Thus every noncentral child comes in a converse pair, and the central family has one self-converse member.

| \(p\) | \(N_{p,m}\) when numerically computed | Exact positive counts of classes by automorphism order \(d\) |
|---:|---:|---|
| 3 | 5 | \(d=1:1\) |
| 7 | 657 | \(d=1:2,\ d=3:3\) |
| 11 | 527,725 | \(d=1:16,\ d=5:3\) |
| 19 | 11,212,997,148,405 | \(d=1:1076,\ d=3:12,\ d=9:3\) |
| 23 | given by T4 | \(d=1:10688,\ d=11:3\) |
| 31 | given by T4 | \(d=1:1292610,\ d=3:100,\ d=5:12,\ d=15:3\) |

In general the central family contains \((1-o(1))2{p\choose m}/(pm)\) rigid classes as \(p\to\infty\). Every nonidentity affine map that can fix a proper mask has one fixed point and all other cycles of length at least three. The union bound gives at most \(pm\,2^{1+(p-1)/3}\) nonrigid masks, exponentially small against the middle binomial layer. This yields unbounded strong-fibre multiplicity along a specified sequence of heights, despite the rigid automorphism group of almost every child.

At \(p=7\) the three attained positive-layer counts are 423, 579 and 657. Their Paley-deletion class counts are 2, 2 and 5. The full order-eight atlas instead has 32, 14 and 15 classes. The prior order-seven maximum 189 shows that eight is the *least possible order* for all three. Thus a parent can explain every realization at 657 with nontrivial automorphism (three by the stored order-eight census), while leaving ten other minimal presentations outside its subtree. [Stored order-eight class table](../../../04-computation/experiments/tiling_line_metagraph_n9_20261020_classes_n8.tsv).

## Reproduction

Run `python3 -B grammar.py 33` for the full eight-representative list, or supply another target. Run `python3 -B replay.py` for all strong parents through order seven, all strong fibres through 33 selected above, all 126 nonconstant Paley-seven cuts with canonical child comparison, and formula (2) through the listed primes. The small spectra reuse the stored exhaustive parent classes and compare each ceiling-filtered mask set to direct enumeration through order five. Compile `c++ -O3 -std=c++17 central_counts.cpp -o /tmp/paley-central && /tmp/paley-central` for the four numerical central heights in the table, by literal Hamiltonian endpoint DP. No theorem about the existence or uniqueness of an arbitrary target follows from the local support bounds.
