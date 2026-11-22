# Isomorphism fibers of a prescribed count and of one physical arc

**PROVED.** For distinct odd counts \(h<k\), forgetting the orientation of the edited edge gives a canonical factorization of every single-arc pair with these counts. For \(k=h+2\) and \(h\ge5\), removing its two exterior singleton runs leaves a pair of strong tournaments. Every such marked pair is rigid. This gives an exact size generating function for physical flips, and sharply separates the strong search at exponent 61 from scalar ordinal-product constructions.

**OPEN.** The upper and lower neighbors of \(E=2^{60}(2^{61}-1)\) are joined by a flip if and only if a strong pair occurs at one of the orders 39 through 157. The upper tournament in every candidate is rigid even when one forgets the arc. The known 66-vertex congruence supplies a different center. [Exponent-31 mechanism](../chanpak_e31_physical_20261122/README.md) · [Exponent-61 congruence](../chanpak_leech_perfect_phases_20261122/EXPONENT-61.md)

The [concurrent common-cut proof](../chanpak_support_birth_grammar_20261122/PRIMITIVE-SEAMS.md) independently establishes the strong middle and the two low-count exceptions, alongside [reversible stripping by Hamilton-path support](../chanpak_support_birth_grammar_20261122/SUPPORT-RIGIDITY.md). Here the double-arc formulation, a lower bound on the cost of crossing a component cut, and the orbit counts give computable polynomials and the order window above. [Exact lower component at exponent 61](../chanpak_e61_fibonacci_target_20261122/PLATEAU.md)

An exact flip at exponent 61 must contain a common directed four-cycle on which the two counts are equal. The marked core grows from this cycle by one- and two-vertex steps. The two-vertex steps preserve the parity of the half-gap, so every successful growth history must contain a one-vertex parity change. [Square, extension grammar and finite slot test](SQUARE-EARS.md)

Every odd count through 3,423,417 other than the two exceptions already occurs at order 14. All smallest representatives of these nontrivial heights are strong, by a component-replacement theorem that applies at every height. [More than a million consecutive odd heights first attained at order 13](ORDER-13-BIRTH-INTERVAL.md) · [Exact interval at order 14](ORDER-14-INTERVAL.md) · [Small-component exclusion](BIRTH-CUT-COVER.md)

## The original single-count series and its missing edge information

Let \(F_h(z)\) count ordinary isomorphism classes at all orders with exactly \(h\) directed Hamiltonian paths. Use the empty tournament only as the neutral object at \(h=1\). Write \(A_d(z)\) for the order polynomial of nontrivial strong classes of count \(d\). The inherited formula for **each** positive odd target, including an absent one, is
\[
 F_h(z)=\sum_{d_1\cdots d_s=h,\ d_i>1}
          \frac{\prod_i A_{d_i}(z)}{(1-z)^{s+1}} .
\]
The sum is over ordered factorizations. Every tournament is an ordered word of its strong components, all paths traverse these components in order, and the \(s+1\) gaps independently admit strings of singleton components. For \(h=1\) this gives \(1/(1-z)\). For an attained \(h>1\), the initial power is the smallest order, and its coefficient is exactly the number of smallest isomorphism classes. For an omitted height the series is zero. The earlier papers already prove this formula, polynomiality at large orders, and minimum-order recursion. Knowing two such series separately does not decide whether any representatives differ in one arc. [SCC multiplicity and converse](../chanpak_collatz_mrb_fibers_20261116/SCC_MULTIPLICITY.md) · [Minimum-order recurrence](../chanpak_birth_atlas_metagraph_20261117/factor/PROOF.md) · [Previous vertex-deletion grammar](../chanpak_core_offsets_and_fiber_grammar_20261114/METAGRAPH.md)

## Strip one edge before stripping vertices

Take a *physical* pair: the two orientations of one specified unordered pair on the same vertex set, up to one simultaneous relabelling. Since \(h<k\), its lower and upper completions are distinguished. Parallel arcs between the same two abstract tournament classes can give different objects here.

Temporarily put **both** directions on the distinguished pair and call the resulting semicomplete digraph \(U\). The two endpoints belong to a common strong component \(C\). Every strong component of \(U\) has a canonical position in a total order, since all arcs between distinct components point uniformly one way. The resulting factorization
\[
                          L\ \longrightarrow\ C\ \longrightarrow\ R
\]
is unique even up to isomorphism. On completing the double arc, \(L,R\) are unchanged and the two new path counts are \(H(L)H(R)\) times the two completion counts of \(C\). Conversely any two exteriors and a strongly connected \(C\) with exactly one double arc produce one marked pair by this construction. Thus if \(W_{p,q}(z)\) counts these central isotypes with completion counts \(p<q\), the exact marked-edge series is
\[
 U_{h,k}(z)=\sum_{de\mid\gcd(h,k)} F_d(z)F_e(z)\,
                         W_{h/(de),k/(de)}(z).
\]
Here \(d,e\) run over **ordered** positive pairs. A strong central completion exists: if one completion is reducible, the other is obtained by reversing an arc from its first strong component to its last. In the original completion, no cut between components is left across which all arcs still go forward after that reversal. More explicitly, all vertices are reachable from the latter endpoint via the new arc and paths through the consecutive components, except in the two-singleton case where both counts are one. If both extreme components are singletons and there is an intermediate component it supplies a forward bypass. If either extreme component is nontrivial, an internal Hamilton cycle supplies a forward bypass of the removed arc. Hence \(W\) is a polynomial, with order bounded by the strong-order cap of the larger count. The common exterior factor divides the gap as well.

## The only unit flips crossing a component boundary

Suppose the lower completion has initial and terminal active components of orders \(a,b\) and \(m\) vertices strictly between them. The switched arc goes from a vertex \(u\) of the first component to \(v\) of the last before the reversal. Choose a full path starting at \(u\) in the first and one ending at \(v\) in the last. These exist by cutting a Hamilton cycle in a strong component, or by taking the singleton. For each nonempty prefix of the former, nonempty suffix of the latter, and any subset of the middle, choose a forward Hamilton path on that subset. Joining and closing with \(v\to u\) gives a new odd directed cycle whenever the selected total order is odd. Distinct choices of lengths and subsets give distinct vertex sets. The only order-two selection is even.

Every old odd-cycle packing survives because no old cycle crosses the forward cut. Each new odd cycle contributes at least two to the gain in the path count by the weighted disjoint-odd-cycle expansion. If \(m>0\), exactly half the middle subsets give the required parity for each prefix and suffix. If \(m=0\), the prefix–suffix rectangle has \(\lfloor ab/2\rfloor\) pairs of opposite parity. Therefore the **half-gain** is at least
\[
    ab\,2^{m-1}\quad(m>0),\qquad \lfloor ab/2\rfloor\quad(m=0).
\]
This also proves that a genuine cross-component reversal increases the count, fixing the direction used above. [Cycle-packing identity, Irving and Omar](https://arxiv.org/abs/2412.10572) · [Earlier exact rank-two gain at a forward cut](../chanpak_whip_seam_20261115/BACK_CROSS.md)

For gain two, \(m>0\) forces \(a=b=m=1\), the transitive triple becoming a triangle. For \(m=0\) we need \(ab\le3\). A strong component has order one or at least three. Thus the only positive cases are \((a,b)=(1,3),(3,1)\), where a triangle and a singleton become the strong four-vertex tournament. These give \(1\to3\) and \(3\to5\) respectively. If \(h\ge5\), coprimality of \(h,h+2\) forces every exterior component to have count one and hence to be a singleton. After stripping the exteriors, both orientations must be strong. Equivalently, **the partial digraph obtained by removing the edited arc entirely is strong** in this central pair: any source cut in the partial digraph could acquire incoming arcs in at most one of its two completions. The converse follows since both orientations contain the partial digraph.

Write \(P_h(z)\) for the polynomial counting these strong partial digraphs with completion counts \(h,h+2\). We obtain
\[
                      U_{h,h+2}(z)=P_h(z)/(1-z)^2\qquad(h\ge5).
\]
Every monomial \(w_s z^s\) contributes exactly \(w_s(n-s+1)\) marked classes at each order \(n\ge s\). In particular the first nonzero coefficient counts the genuinely smallest physical flips, and the eventual slope counts all central flip types. For the two exceptions the numerators for the full marked-pair series are \(z^3\) at \(1,3\) and \(2z^4\) at \(3,5\). In the latter case the singleton may be on either side and the three triangle ports form one rotation orbit.

## Why the mark removes all automorphisms

Every automorphism of the partial digraph fixing the two completions acts freely on the ordered Hamilton paths of each. An automorphism interchanging the completions would interchange their distinct counts, so this is the entire group. Its order divides \(\gcd(h,k)\). In particular each unit pair has a trivial group. In a given completion the endpoints of an admissible unit arc form a base for its automorphism group: only the identity fixes them both. If equally many paths are possible in both completions, the interchange argument no longer applies.

This gives an isomorphism-counting algorithm without canonical labelling for unit edges. With the lower arc fixed as \(0\to1\), enumerate digraphs on the other pairs and retain exactly the requested completions and central strongness. Divide by \((n-2)!\). The action of permutations fixing 0 and 1 on these accepted encodings is free. One can instead use an existing lower-class catalogue and count acceptable labelled arcs in each class divided by its automorphism order. These are equivalent counts of marked objects.

Under converse the partial digraph's two exterior singleton runs are interchanged. If \(P_h^{*}\) counts self-converse central pairs, the self-converse marked series is \(P_h^{*}(z)/(1-z^2)\). At a rigid upper endpoint any anti-automorphism is unique if it exists and is an involution. For the marked upper class to be self-converse it must exchange the two endpoints of the arc. Thus at most \(\lfloor n/2\rfloor\) of its arcs can have this symmetry.

## The order-61 exponent as a prescribed absolute target

Busch's strong-order bound gives \(B(E-1)=B(E+1)=157\). For the opposite direction, split the positions of a Hamilton path on 38 vertices into four contiguous blocks of sizes 8,10,10,10. For each ordered partition of the vertex set, ignore the three interblock constraints. Using the maxima 661 and 15,745 at orders 8 and 10 gives
\[
 H(T)\le \frac{38!}{8!(10!)^3}\,661\,(15745)^3
       =700390589337447041726395546833990000 < E-1.
\]
This also excludes smaller orders by singleton padding. Thus it suffices, and is necessary, to search strong partial digraphs with **exactly one missing arc** in orders 39 through 157. No exhaustive search over that window is claimed. [Busch's strong minimum](https://doi.org/10.37236/1141) · [THM-4889, maxima by order](../../../01-canon/theorems/THM-4889-the-hamiltonian-path-spectrum-by-order-and-by-value-exact-censuses-to-order-11-resolved-by-cyclic-triangles-least-orders-strong-windows-level-floors-and-the-zone-law.md)

More generally, if \(\gcd(h,B(h)!)=1\) for \(h>1\), every tournament with count \(h\) has trivial automorphism group at **every** order. Each nontrivial component has count \(d\mid h\) and order at most \(B(d)\le B(h)\). Its automorphism order divides both \(d\) and its order factorial. Each component is therefore rigid, and the component order prevents exchanging isomorphic components. In particular \(\gcd(E+1,157!)=1\) by ordinary Euclidean division. Thus the coefficient of \(F_{E+1}\) at each order \(n\) equals the number of labelled tournaments with that height divided by \(n!\), and multiplicity among the smallest upper witnesses, if present, must come from genuinely different rigid types. The same rigidity test applies to all prime heights \(p\ge5\) since \(B(5)=4\) and the strong-order floor exceeds \(n\) for every \(n\ge4\).

## Two marked edits can join two unique smallest classes

**FINITE-EXACT with all-order scope.** At counts 11 and 13 there is a unique strong type of order five for each, and no other strong types at either height. There are nevertheless *two* inequivalent arcs from the former to the latter. Thus both ordinary series are \(z^5/(1-z)^2\), while their marked-flip series is \(2z^5/(1-z)^2\). At counts 25 and 27 both first appear at order six; their sole central order for a physical edge is seven, with two marked edges. The polynomials are \(P_{25}=2z^7\) and \(P_{27}=2z^6+2z^7\). Among the pairs \((h,h+2)\) with \(1\le h\le73\) where both scalar heights are attainable, the **entire marked-edge series vanishes** precisely at lower counts
\[
                    35,\quad39,\quad53,\quad55,\quad59,\quad63.
\]
For instance 63 and 65 have individual witnesses but no single arc in any order joins them. These conclusions follow by taking every strong lower isotype through order eight and every arc, dividing each resulting multiplicity by its group order. The strong minimum at order nine is 75, making the table complete through the stated lower count. [Earlier smallest-order 25/27 obstruction](../chanpak_birth_atlas_metagraph_20261117/metagraph/PROOF.md) · [Underlying class catalogues](../../../session212/cycle-levels)

Run `python3 -B fibers.py --height 11` to print both rational series. Add `--order 100` to return their coefficients at that order. The `--small-table` option computes all ordinary birth multiplicities and unit-edge numerators for odd lower heights through 73. The default checks the formulas against direct enumeration at orders through six, and `--e61` checks both order bounds and the global rigidity certificate. The uncapped general series algorithm terminates for every positive odd height by the strong-order bound. Beyond order eight it generates each strong class under the target ceiling once using the [canonical deletion and minimum-cut grammar](../chanpak_minimal_core_grammar_20261122/README.md). Large ceilings can still require impractically many classes.

## The five Fermat-prime birth fibers

The known Fermat-prime heights are a concrete family in which prime factorization cannot split a minimum witness. Every one of the types is strong. Their smallest orders were already determined in [THM-4861, least orders of the five Fermat primes](../../../01-canon/theorems/THM-4861-the-five-fermat-primes-are-hamiltonian-path-counts-of-tournaments-with-least-orders-3-4-6-8-11-a-vertex-transitive-tournament-has-a-prime-count-only-for-c3-and-fermat-primes-are-free-vertices-of-the-prime-square-sum-graph.md). Resolving the counts by cyclic triangles in the inherited exhaustive census gives the following multiplicities.

| Hamiltonian paths | Least order | Smallest classes | Number of cyclic triangles and corresponding class count |
|---:|---:|---:|---|
| 3 | 3 | 1 | 1 triangle, 1 class |
| 5 | 4 | 1 | 2 triangles, 1 class |
| 17 | 6 | 2 | 4 triangles, 2 classes |
| 257 | 8 | 30 | 12, 13, 14 triangles with 2, 18, 10 classes |
| 65,537 | 11 | 13,602 | 48, 49, 50 triangles with 164, 13,092, 346 classes |

Every class in the last four rows is rigid by the prime-height test. For 65,537 alone, 13,092 equally small rigid tournaments even share their triangle count. Thus the multiplicity persists after fixing the target, the size, the group, and this additional cycle statistic. The jump from 1 to 13,602 among these examples does not give an all-prime growth law. Run `python3 -B fibers.py --fermat-births`. The option `--ordinary --height 65537 --order 11` reads this coefficient directly. [Exhaustive order-11 census](../path_spectrum_two_angles_metagraph_20261106/data/census_11.txt.gz) · [Fermat and Platonic rotation restrictions](../chanpak_fermat_solids_seams_20261122/README.md)

## Multiplicity at the first order can span the edit graph

**PROVED consequence of the cited typical-birth estimate.** Fix \(0<\delta<1/2\) and choose an ordinary isomorphism class uniformly among all order-\(n\) classes. With probability \(1-\exp(-\Omega(\sqrt n))\), its count is born at order \(n\), and its fiber there contains at least
\[
         2^{(1-H_2(\delta))\binom n2-O(n\log n)}
\]
classes pairwise more than \(\lfloor\delta\binom n2\rfloor\) arc reversals apart even after arbitrary relabelling. Here \(H_2(x)=-x\log_2 x-(1-x)\log_2(1-x)\). This probability samples a class uniformly. The induced distribution on odd integers weights a height by the size of its fiber.

Indeed put \(M=\binom n2\), and write \(t_n\) for the number of all order-\(n\) classes. We have \(t_n\ge 2^M/n!\) and at most \(n!\) possible path counts. The fraction of classes in fibers of cardinality below \(t_n e^{-\sqrt n}/n!\) is therefore at most \(e^{-\sqrt n}\). The inherited large-order concentration and maximum estimates give the same exceptional scale for failure to be a first birth. On the intersection the entire fiber consists of smallest representatives. A radius-\(r\) ball about one isotype contains at most \(\sum_{i\le r}\binom Mi\) isotypes: fix a labelling of its center and realize any quotient-distance minimizer by flips in that labelling. Greedily select one type and discard its ball. The ordinary Hamming-ball entropy bound now gives the displayed packing. In particular first-birth multiplicities can grow without bound and can consist of widely separated representatives, even after identifying a type with its converse (this at most doubles the ball bound). The elementary packing argument holds without an asymptotic input for any given finite fiber of known cardinality. [Quantitative first-birth theorem](../chanpak_birth_atlas_metagraph_20261117/families/PROOF.md) · [Divoux–Wang, maximum Hamilton path count](https://arxiv.org/abs/2610.23755)

## Two small components can never coexist in a smallest witness

**PROVED using finite-exact input.** At **every** attained height greater than one, the orders of any two distinct strong components of a minimum-order representative sum to at least **16**. Consequently all smallest representatives of orders at most 15 are strong, and a smallest representative on \(n\) vertices has at most \(\lfloor n/8\rfloor\) strong components if it has two or more. The census argument below reaches sum 13 and explicit smaller witnesses settle sums 14 and 15. [Proof and witness cover for the last two sums](BIRTH-CUT-COVER.md)

Write \(b(h)\) for the least order at height \(h\). Each component of a smallest witness must itself have order \(b(a)\) at its height \(a\), by substitution. For any two components of heights \(a,d\), even when separated in the ordered word, the original count can also be obtained by placing one tournament of height \(ad\) at the earlier position and omitting the later component. All other component heights are unchanged. Thus a necessary condition on **every pair** in a smallest word is
\[
                         b(ad)=b(a)+b(d).
\]
The same reasoning works for any submultiset of components. In particular, a single strict inequality excludes coexistence at every height with those two factors represented by those components. A singleton is already redundant. This strengthens the inherited [individual-component minimality criterion](../chanpak_birth_atlas_metagraph_20261117/factor/PROOF.md).

Read the first orders of all heights in the ordinary census through order 11. For a candidate sum \(s\le14\), the two component orders are at least 3 and at most 11, so all candidate height pairs are in this census. For each unordered pair of heights (including equal heights), look up the product.

| Sum of orders | Candidate height pairs | Products realized below that sum | Least certified saving |
|---:|---:|---:|---:|
| 6 | 1 | 1 | 1 |
| 7 | 1 | 1 | 2 |
| 8 | 5 | 5 | 1 |
| 9 | 16 | 16 | 2 |
| 10 | 80 | 80 | 2 |
| 11 | 349 | 349 | 2 |
| 12 | 1,753 | 1,753 | 2 |
| 13 | 9,011 | 9,011 | 2 |
| 14 | 51,818 | 22,565 | 3 |

The strict savings through sum 13 exclude those pairs everywhere. Adding the two witness covers and summing over the pairs of a word of length \(k\ge2\) gives \((k-1)n\ge16\binom{k}{2}\). At sum 14 the first unclassified product in the stored ordering is \(3\cdot26{,}839=80{,}517\). The separate witness cover realizes this product on at most 12 vertices. Thus 14 is the demonstrated exclusion boundary of this census argument. Run `python3 -B birth_cuts.py` to reconstruct the table. [Order-resolved census](../path_spectrum_two_angles_metagraph_20261106/data)

## When does the first coefficient equal one?

For an attained \(h>1\), list the unordered multisets of **actual strong isomorphism types** whose heights multiply to \(h\) and whose orders sum to \(b(h)\). If one multiset contains \(k\) types counted with repetitions \(m_i\), it contributes exactly \(k!/\prod_i m_i!\) smallest representatives. Every permutation gives the same scalar height. Two words are isomorphic exactly when their component sequences agree term by term. Hence uniqueness occurs precisely when there is one optimal multiset and it consists of copies of a single strong type. In particular, at a height which is not a proper power, a unique smallest representative must be strong. As the Fermat-prime table shows, a single strong component can itself have many different isomorphism types.

There is a divisibility prediction even before counting these types. Let \(g=\gcd_{q\mid h} v_q(h)\), over prime divisors. Cyclically permute the \(k\) positions of a smallest word on \(n\) vertices. If its stabilizer has order \(d\), the word repeats a block \(d\) times, so \(d\mid k,n,g\). Each orbit thus has cardinality divisible by \(k/\gcd(k,n,g)\). That integer divides the number of minimum classes with exactly \(k\) components. For a height that is not a proper power it is \(k\). Minimum words with a primitive period can also be counted by the usual necklace inversion. If \(W(h,n,k)\) counts all ordered words of \(k\) nontrivial strong components of total order \(n\) and height \(h\), then the number of cyclic classes of aperiodic such words is
\[
 \frac1k\sum_{d\mid\gcd(k,n,g)}\mu(d)\, W(h^{1/d},n/d,k/d).
\]
This follows by Möbius inversion on the number of repetitions followed by division by \(k\) on free orbits. Taking \(n=b(h)\) restricts the answer to smallest representatives.

The recursion also assigns an exact integer address to every class. At state \((n,h)\), a specific first component of order \(s\) and height \(d\mid h\) owns an interval of length \([z^{n-s}]F_{h/d}\). The empty word has weight one exactly at \((0,1)\). Place the first components in canonical order and recurse inside their intervals. The ordered component word is unique, so the intervals partition the whole coefficient. A uniformly drawn integer in this coefficient therefore samples each **unlabelled** class with equal probability. This is a concrete link to exact Latin-square sampling: local steps must account for the mass of their possible completions. Here the side information required for an ordinal suffix is only its remaining order and target quotient. Complexity still includes obtaining the strong alphabet. [Mira Cult, exact Latin-square sampling and history weights](https://github.com/Mira-Cult/math/blob/30ddaaad27b48df2002225c3e3f40cffd02c5253/notes/combinatorics/exact-latin-square-sampling/README.md)

Run `python3 -B rank_classes.py --height 17 --unrank 1 --rows` for the second of the two smallest 17-path classes, or add `--order` to work at a specified size. The inverse accepts a representative under any vertex labelling. The program reports when its adjustable strong-component enumeration cap leaves a query unresolved. [Class ranking and construction](rank_classes.py)
