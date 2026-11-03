# Geometry and metagraphs share Cartesian powers, but not their quotients

**Status.** Explicit constructions and proofs, with a seven-tournament finite certificate. The external large-girth theorem is credited to Tomon (2026). The new work here is the exact tournament interface, the simultaneous Moser-spindle seed, and the resulting profile-invisible cubes. No Hadwiger–Nelson bound or literature-priority claim is made. “Realized in the plane” permits crossing edges; it does not mean a planar graph.

## 1. Definitions

Let M_n be the simple graph whose vertices are isomorphism classes of n-vertex tournaments. Distinct classes are adjacent when representatives differ by exactly one arc reversal. Isomorphism-preserving flips are discarded as loops. Let S_n be its induced subgraph on strongly connected tournaments. (Here S_n denotes this graph only, not the symmetric group, which is written Sym(n).)

For a digraph T, let p_j(T) count its vertex partitions into j nonempty directed paths, with each path carrying its directed order and with no ordering of the collection of paths. Put

P_T(z)=sum_j p_j(T) z^j,
F_T(r)=sum_j p_j(T)(r)_j.

Thus p_1(T)=H(T). The previous session's matching product satisfies P_(T⊕U)=P_T⋆P_U, or equivalently F_(T⊕U)(r)=F_T(r)F_U(r). The ordinal sum T⊕U points every cross arc from T to U.

## 2. Exact Cartesian-product theorem

**Theorem 1.** For any induced subgraphs H_i of S_(n_i), with n_i≥3, the map

([T_1],...,[T_k]) ↦ [T_1⊕...⊕T_k]

is an induced embedding of H_1 □ ... □ H_k into M_(n_1+...+n_k).

**Proof.** The strongly connected components of the ordinal sum are exactly its blocks, and the condensation is the uniquely ordered transitive tournament. Every isomorphism therefore preserves the positions of the blocks. This proves injectivity.

A within-block reversal that leaves that block strong changes only one factor; if it destroys strong connectivity, it increases the number of SCCs and leaves the displayed family. For a reversal between two blocks, the original direction was from an earlier block to a later block. The new backward arc, together with strong paths inside the endpoint blocks and the unchanged forward arcs, makes every block in that interval mutually reachable. The interval merges into one SCC, so the number of SCCs decreases and the result is outside the family. Consequently the only adjacencies inside the family are single-factor adjacencies. ∎

This is stronger than a homomorphism or a subgraph inclusion: no extra edges survive after passing to unlabelled tournament classes.

### 2.1 An induced strong wrapper

**Theorem 2.** For every n≥1, M_n embeds inducedly into S_(5n+10).

Take q=2n+5. Let A and B be fixed regular tournaments of order q, for example the cyclic tournament with differences 1,...,(q−1)/2. Define W(T) on A,T,B with A→T→B→A. It is strong. Its outdegrees are:

- every A vertex: 2n+2;
- every T vertex: 2n+5+d_T(v), in [2n+5,3n+4];
- every B vertex: 3n+7.

These three ranges are separated by gaps of at least three. An isomorphism of wrappers preserves the three sets, so W(T)≅W(U) iff T≅U. After one reversal every degree changes by at most one, and the three bands remain recognizable. If either endpoint lies in a guard A or B, that guard has a vertex whose degree differs from its required constant, so the result is not any W(U). Reversals inside T give exactly W(T^e). This proves the induced embedding. ∎

Together, Theorems 1 and 2 make the class of finite induced subgraphs of tournament metagraphs closed under Cartesian products, after enlarging tournament order. The same is true using only strong tournaments as metagraph vertices.

## 3. A literal Moser spindle in S_6

Encode a labelled tournament on {0,...,5} by a bit mask over pairs in lexicographic order (01,02,03,04,05,12,...,45), with bit 1 meaning i→j for i<j. The following seven masks are all strong and pairwise nonisomorphic:

| Spindle point | O | A | B | C | D | E | F |
|---|---:|---:|---:|---:|---:|---:|---:|
| Tournament mask | 18 | 22 | 26 | 181 | 81 | 146 | 145 |

Their induced metagraph has exactly these eleven edges:

OA, OB, AB, AC, BC, OD, OE, DE, DF, EF, CF.

The file seed_certificate.json records every minimum quotient Hamming distance, an explicit relabelling and flipped edge for each adjacency, and the full profiles. The verifier checks all 720 relabellings for every pair, so it proves the nonedges as well as the edges.

### 3.1 The same graph as an induced planar unit-distance graph

Use O=(0,0), A=(√3/2,1/2), B=(√3/2,−1/2), C=(√3,0). Let R be rotation with cos θ=5/6 and sin θ=√11/6, and put D=RA, E=RB, F=RC.

Every listed edge has length one, including CF because |C−F|²=6(1−5/6)=1. No other pair has distance one. All 21 distances are checked exactly in Q(√3,√11), with no floating-point tolerances.

The graph is not 3-colorable: the diamond OABC forces O and C to have the same color, and ODEF forces O and F to have the same color, contradicting the edge CF. A proper 4-coloring is included in the certificate. Thus its chromatic number is exactly four.

## 4. Arbitrary girth in both realizations at once

We use the following external theorem, not claim it as a result of this session:

**Tomon's Cartesian-power theorem.** For every finite graph H and integer g, some Cartesian power H^□k has an induced subgraph J with girth at least g and χ(J)=χ(H). [Tomon, Theorem 1.5 / 5.2](https://arxiv.org/html/2610.03517v2).

We also use the established fact that Cartesian products of induced unit-distance graphs in R^d, d≥2, are induced unit-distance graphs in the same R^d. Its geometric proof takes generic independent rotations of the factors and maps a tuple to the sum of its rotated coordinates. Finitely many collision and accidental-unit-distance equations are avoided. See Horvat–Pisanski (2010), Theorem 3.4, or [Tomon, Lemma 5.1](https://arxiv.org/html/2610.03517v2).

**Corollary 3.** For every g there are k and a finite graph J such that:

1. J has girth at least g and chromatic number exactly four;
2. J has an induced unit-distance realization in R²;
3. J is an induced subgraph of M_(6k), represented by ordinal sums of the seven certified strong tournaments;
4. J is also an induced subgraph of S_(30k+10), after applying W.

**Proof.** Apply Tomon's theorem to the single common spindle H. Realize its Cartesian power geometrically by rotated sums and tournament-theoretically by Theorem 1. Restrict both realizations to the same vertex set of J. Apply Theorem 2 for the last statement. ∎

This is a shared abstract construction with two genuine realizations. It is not an assertion that arbitrary metagraphs are unit-distance graphs. There is no explicit small value of k or small high-girth certificate in this packet; the arbitrary-girth part depends on the cited existence theorem.

## 5. Entire geometric cubes invisible to all path-cover coefficients

The adjacent masks 18 and 22 have exactly the same profile:

P(z)=15z+93z²+138z³+73z⁴+15z⁵+z⁶.

For every k, form all 2^k ordered ordinal sums whose blocks are chosen independently from these two strong six-vertex tournaments. By Theorem 1 they induce Q_k in M_(6k), and they are pairwise nonisomorphic. Yet every one has

P_T=P^⋆k,   F_T(r)=F_P(r)^k,   H(T)=15^k.

The cube Q_k has an induced unit-distance realization in the plane: map x∈{0,1}^k to sum_i x_i u_i for generic unit vectors u_i. Choosing the vectors off finitely many real-algebraic exceptional sets makes the map injective and eliminates all accidental unit distances. Thus even the complete path-cover profile is constant on an exponentially large, geometrically realizable induced metagraph cube.

### 5.1 Profile equivalence survives module substitution

**Lemma 4.** If P_X=P_Y and |X|=|Y|, replacing a module X by Y in a fixed digraph template preserves the entire path-cover profile.

**Proof.** Restrict a path cover of the substituted digraph to a module, cutting every path at its external excursions. This gives a local path cover of that module. For each j choose a bijection between the j-path covers of X and Y, and for each matched pair choose a bijection of their constituent paths. Replace each local path segment by its matched path, preserving which external segments enter and leave it. Every outside vertex has uniform adjacency to the module, so all external joins remain valid. The reconstruction preserves the path-incidence graph, hence the number of global paths and absence of cycles. Using the inverse chosen bijections reverses the map. ∎

It follows that W carries the profile-invisible Q_k into S_(30k+10), still with one common complete path-cover profile. Their common Hamilton-path number after wrapping need not equal 15^k. The unwrapped count is 15^k; the wrapped one changes with the guards.

The result is a strong limitation on using a path-cover scalar or the complete FI-interface from S165 as a canonical metagraph label. It does not contradict its usefulness for ordinal-sum composition.

## 6. A second bridge: the labelled cube and parity holonomy

For n≥2, the labelled n-vertex tournament flip graph is Q_m, m=binom(n,2). It therefore admits an induced unit-distance realization in R², despite generally being nonplanar as an abstract graph.

Write ε(T)=(−1)^(number of 1 bits). Under vertex relabelling σ,

ε(σT)=sgn(σ) ε(T).

Indeed, permuting edge coordinates preserves the product of their signs, while precisely the inversions of σ reverse the convention i<j. Every tournament automorphism is even. One elementary proof is that an automorphism cannot have even order: an involution would exchange two vertices and reverse their arc; an even-order automorphism would have an involutive power. Thus the automorphism group has odd order, and every element has even permutation sign.

Consequently the quotient of Q_m by Alt(n) has exactly two vertices over every tournament isomorphism class, one of each cube parity. After removing the edges that project to isomorphism-preserving flips, it is precisely the bipartite double cover of the loopless M_n. If such flips are retained, the matching edge between the two copies must be retained too; forgetting this convention gives the wrong cover.

For a closed walk of length ℓ in M_n, lift its steps to labelled one-flip steps beginning at T and ending at σT. Then

sgn(σ)=(−1)^ℓ.

The sign of the closing relabelling is unambiguous because every automorphism is even. Odd metagraph cycles are therefore exactly odd relabelling holonomy. The plane embedding of Q_m does not descend to a plane unit-distance embedding of M_n: its coordinates need not be constant on relabelling orbits, and even its two-coloring changes under odd relabellings. This is the exact obstruction behind the failed “the labelled geometry proves a chromatic bound on the quotient” inference.

### 6.1 A path-count orientation coordinate

Define the alternating Hamilton-path enumerator

A(T)=sum_(directed Hamilton words π of T) sgn(π).

This is not the repository's arc-product signed Hamiltonian permanent. By Rédei's oddness theorem A(T) is odd and nonzero. It transforms by A(σT)=sgn(σ)A(T). Hence A(T)>0 chooses one of the two alternating-group sheets for each isomorphism class, and ε(T)A(T) is a genuine label-independent odd integer. This links the same parity line to path enumeration, without confusing the ordinary positive profile with its alternating refinement.

For the seven certified masks, A(T) in point order O,A,B,C,D,E,F is −1,1,−5,3,9,−11,1. A signed refinement supplies orientation data; it does not remove the metagraph's odd-cycle obstruction.

## 7. Reproducible scope and next questions

Run `python3 check.py` and `python3 -O check.py`. The checks use only Python's standard library, exact integers and fractions, and explicit failures rather than assert statements. They establish the finite seed and controls; the proofs above establish the arbitrary-k statements.

The interesting remaining question is not whether triangle-free chromatic complexity exists in tournament metagraphs: it can be transported by the common Cartesian construction. More discriminating questions are whether one can keep tournament order small, prescribe a strong Hamilton-path count while retaining a high-chromatic induced subgraph, or find a common induced plane-unit-distance/metagraph seed with chromatic number five. None of those stronger claims is established here.

## Credits and boundaries

- S165 [FI interfaces](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/results/session165_geometry_representation_20261103/fi_interfaces/FI_INTERFACES.md) supplies the full path-cover interface and the shared-slot interpretation of ordinal joins.
- B. Horvat and T. Pisanski, [*Products of unit distance graphs*](https://doi.org/10.1016/j.disc.2009.11.035), Discrete Mathematics 310 (2010), 1783–1792, establishes the geometric Cartesian closure used here.
- I. Tomon, *Geometric triangle-free graphs of large chromatic number*, [arXiv:2610.03517v2](https://arxiv.org/abs/2610.03517v2), supplies the large-girth theorem. Its results are external inputs, not discoveries of this session.
- Rédei's oddness theorem is classical; a modern account is Schweser–Stiebitz–Toft, [The Tournament Theorem of Rédei revisited](https://arxiv.org/abs/2510.10659).
- The Moser spindle and bipartite-double-cover constructions are standard. The literal seven-tournament realization and its profile-invisible edge are verified here; novelty against the full literature has not been checked.

## Appendix: the shared-root chain has a finite order threshold

The companion unit-distance packet defines G_m using vertices O and t_i,a_i,b_i, with diamond edges Oa_i,Ob_i,a_i b_i,a_i t_i,b_i t_i and consecutive tip edges t_i t_(i+1). Its G_2 is the spindle above.

The next member G_3, on ten vertices, is an induced subgraph of S_7. In point order

O,t_0,a_0,b_0,t_1,a_1,b_1,t_2,a_2,b_2

its tournament masks are

801,296,673,4712,825,809,817,1584,1589,1593.

`check_chain.py` verifies all 45 pairs against all 5040 relabellings, confirms strong connectivity, and checks the profiles by two independent algorithms. The resulting `chain7_certificate.json` records the exact distances and edge witnesses.

An additional complete finite search gives:

- G_2 does not occur inducedly in M_5. Since M_4 has only four vertices, tournament order six is minimal for an induced spindle.
- G_3 does not occur inducedly in M_6. Therefore order seven is minimal for G_3, and the displayed minimal realization can already be chosen strong.
- G_4 does not occur inducedly in M_7. Consequently neither does any longer G_m, since G_4 is an induced subgraph of each longer chain. This does not rule out order eight or higher, or non-induced embeddings.

The census generator enumerates every labelled tournament for n=5,6,7, partitions them into the full permutation orbits, and records every class's one-flip neighbors. It finds 12,56,456 classes and, respectively, 30,290,4086 loopless edges. The induced matcher uses only necessary degree bounds and exhaustive domain backtracking, with no probabilistic pruning.

Reproduce these additional finite statements with:

```
g++ -O3 -std=c++17 census.cpp -o census
for n in 5 6 7; do ./census "$n" > "metagraph_n$n.json"; done
python3 search_induced.py
python3 check_chain.py
python3 -O check_chain.py
```

`chain_embedding_search.json` is the resulting search record. The positive seed is independently checkable without trusting the census; the nonexistence statements depend on the exhaustive census and search. No complete enumeration beyond order seven is asserted.
