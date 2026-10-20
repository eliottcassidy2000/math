# Every regular outer endcap of order at least five has a three-black bypass

**PROVED using CITED Thomassen and independent FINITE-EXACT boundary
certificates:** every regular tournament of odd order at least seven is
Hamiltonian-connected. Together with the inherited five-vertex boundary
identities, every blue line at \(S^{op}\to Z\to S\), with \(S\) regular
of order at least five, has a black-only detour of length at most three.
This covers every marked blue presentation, not just standard circulant
ones. It is not a theorem of global black-core connectivity.

## 1. Inheritance and the missing finite boundary

The current source is
[SCC boundary moves, Sections 6 and 9](../scc_basins_20261020d/moves.md).
It proves an outer-endpoint criterion: replace a Hamiltonian path in the
terminal SCC by one ending at the same vertex but in a different
automorphism orbit, and a three-black bypass follows from the inherited
common black hub. Section 9 supplies such a replacement whenever the
outer SCC is Hamiltonian-connected.

The regular connectivity estimate already handles orders at least
eleven. The precise remaining finite question was whether all three
regular order-seven classes and all fifteen regular order-nine classes
are Hamiltonian-connected. Ordinary total path counts alone cannot
answer this: endpoints are the missing coordinate.

The old
[regular-tournament census](../lrc_regular_full_census_chi_s592.out)
records the class counts and orbit masses. The new checker uses the
authenticated order-nine representatives from
[the complete order-nine census](../blue_quotient_20261020b/census.md),
but also independently certifies completeness of its regular subset.
No new large-order census is performed.

## 2. Exact endpoint matrices and independent completeness certificates

For each representative compute
\[
 E_T(u,v)=\#\{\text{Hamiltonian paths starting at }u
                         \text{ and ending at }v\}.
\]
Hamiltonian connectivity is exactly \(E_T(u,v)>0\) for every \(u\ne v\).
One computation uses subset dynamic programming with the starting
vertex fixed. A separate depth-first traversal enumerates the actual
Hamiltonian paths. The two complete matrices agree for every class.

The order-seven representatives and endpoint-count distributions are:

| Representative | Outgoing-neighbor bit rows | Aut order | \(H\) | Off-diagonal count distribution |
|---|---|---:|---:|---|
| Carousel | \(14,28,56,112,97,67,7\) | 7 | 175 | \(2^7,3^{14},5^{14},7^7\) |
| Paley | \(22,44,88,49,98,69,11\) | 21 | 189 | \(1^{21},8^{21}\) |
| Triangle reversal | \(26,28,41,112,100,67,7\) | 3 | 171 | \(1^3,3^9,4^{12},5^{15},6^3\) |

Here \(a^b\) in the last column means “count \(a\) at \(b\) ordered
endpoint pairs”; it is histogram notation, not exponentiation.
All 42 off-diagonal entries in each row are positive.

For the fifteen order-nine classes, all 72 entries are positive.
The minimum endpoint count over all classes and ordered endpoint pairs
is 21. Each class's rows, total count, automorphism order, minimum and
endpoint histogram are printed in regular_check.out.

**Independent coverage check.** For a list of desired outdegrees
\(\mathbf d=(d_0,\ldots,d_{n-1})\), the number \(F(\mathbf d)\) of
labelled tournaments with those degrees obeys
\[
 F(\mathbf d)=
 \sum_{\substack{W\subseteq\{1,\ldots,n-1\}\\|W|=d_0}}
 F\bigl((d_j-[j\notin W])_{j=1}^{n-1}\bigr).             \tag{1}
\]
This fixes precisely which later vertices are beaten by vertex zero.
Invalid demands have count zero. Sorting the remaining demand list
only merges equal subproblems; all labelled choices of \(W\) remain
in the sum. The checker uses integer memoization of (1).

It independently checks pairwise nonisomorphism and enumerates every
automorphism of each selected representative. Their orbit sizes give

| Order | Pairwise nonisomorphic classes | Sum of \(n!/|\mathrm{Aut}|\) | Independent \(F((n-1)/2,\ldots,(n-1)/2)\) |
|---:|---:|---:|---:|
| 7 | 3 | 2,640 | 2,640 |
| 9 | 15 | 3,230,080 | 3,230,080 |

The disjoint class orbits exhaust the independently counted labelled
universe. Thus the finite Hamiltonian-connectivity conclusion does not
rest only on an expected class count or on random discoveries.

The compressed order-nine source is additionally pinned by SHA-256
~~~text
959c2276c57108d55b56d4906edfbd761fb9737e708b14114d66dd4ad618db0c
~~~
Its pair-code convention scans lexicographic pairs from the most
significant bit, with bit one meaning \(i\to j\). The new checker decodes
this directly and verifies the tournament axioms and every degree.

## 3. The all-order regular Hamiltonian-connectivity corollary

Suppose deletion of \(k\) vertices from a regular tournament of order
\(r\) leaves a separating cut \(A\to B\), with sizes \(a,b>0\) and
\(a+b=r-k\). Sum outdegrees over \(A\) and indegrees over \(B\). Since
every original degree is \((r-1)/2\),
\[
 r\geq a+2b,\qquad r\geq2a+b.
\]
Their sum gives \(k\geq r/3\). The same argument with \(k=0\) proves
strong connectivity. Therefore every regular tournament has vertex
connectivity at least \(\lceil r/3\rceil\).

**CITED input.** Thomassen's
[Hamiltonian-connected tournaments](https://www.sciencedirect.com/science/article/pii/0095895680900611),
*Journal of Combinatorial Theory B* 28 (1980), 142--163,
DOI 10.1016/0095-8956(80)90061-1, proves that every 4-strong tournament
is Hamiltonian-connected in the ordered-endpoint sense. Its publisher
abstract was retrieved by search on this pass; direct article access
returned HTTP 403, so this note imports the stated theorem rather than
claiming to have audited its original proof.

For odd \(r\geq11\), the connectivity lower bound is at least four.
The cited theorem handles those orders. Section 2 handles seven and
nine. Hence:

> Every regular tournament of odd order at least seven is
> Hamiltonian-connected.

No literature novelty claim is made for this graph-theoretic corollary.
Its role is to close the finite boundary of the repository's detour
argument with independently replayable endpoint certificates.

The cutoff is substantive. The directed triangle and the regular
five-tournament are not Hamiltonian-connected: the checker prints three
and five missing ordered endpoint pairs, respectively. At order seven
some endpoint counts equal one, so no bound of two per pair is inferred.

## 4. Transport to a three-black detour

Let \(S\) be Hamiltonian-connected and fix a blue Hamiltonian path
presentation at
\[
                         U=S^{op}\to Z\to S.
\]
Let the terminal \(S\)-path start at \(a\) and end at \(w\). Because \(S\)
is strong, \(w\) has both an inneighbor and an outneighbor. Choose
\(a'\) on the opposite side of \(w\) from \(a\). Hamiltonian connectivity
supplies an \(S\)-path from \(a'\) to \(w\).

No automorphism fixing \(w\) can send \(a\) to \(a'\), since it would
reverse their relative adjacency to \(w\). The new terminal path
therefore lies in a different automorphism orbit from the old one.
Change only that component path. If the new whole path were also blue,
composition of the two positional anti-reversals would give precisely
the forbidden automorphism between these terminal paths. Hence the new
whole path is black.

The first, second and last vertices of the whole word are unchanged:
the initial SCC has at least three vertices and the terminal vertex is
still \(w\). The inherited common-black-hub lemma applies to the
original blue partner and this new black partner, supplying the other
two black moves. Their concatenation has length three.

By Section 3 this works for every regular \(S\) of order at least seven.
For regular order five use the explicit three-black boundary identities
of [moves.md Section 3](../scc_basins_20261020d/moves.md), which do not
assume Hamiltonian connectivity. Therefore:

> Every blue line at \(S^{op}\to Z\to S\), for regular \(S\) of order
> at least five, has a black-only bypass of length at most three,
> for every blue marked path and arbitrary middle \(Z\).

The source's blue presentation itself supplies the required converse
matching of the outer SCCs and blue middle presentation. There is no
assumption that \(Z\) is pure blue or nonempty.

For completeness, the inherited order-three reduction in moves.md
Section 6 applies when \(S=C_3\) and the whole source is mixed: then
the middle is mixed, and its black path gives the same length-three
bound. Thus every mixed nonstrong source with nontrivial regular outer
SCCs is covered. Singleton outer SCCs have their own center theorem;
nonregular larger outer SCCs remain a separate frontier.

This maps an endpoint-existence certificate to an actual metagraph
detour. It preserves the selected blue line's endpoints in the
converse-merged metagraph. It does not identify distinct black components
globally, prove the ordinary \(H\)-spectrum conjecture, or replace the
regularity-preserving arc moves considered in spectrum.md.

## Reproduction

~~~text
python 05-knowledge/results/rotation_owners_spectrum_20261020e/regular_check.py
python -O 05-knowledge/results/rotation_owners_spectrum_20261020e/regular_check.py
~~~

NetworkX is used for independent isomorphism and automorphism checks.
All checks use explicit exceptions, and the two stdout receipts agree.
The exact controls concern 18 small classes, their full endpoint matrices,
two independently counted labelled universes, and the two small hostile
regular tournaments. The unbounded conclusion additionally uses the
proved cut inequality and the explicitly cited Thomassen theorem.
