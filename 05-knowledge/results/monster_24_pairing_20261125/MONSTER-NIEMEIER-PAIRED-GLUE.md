# When the 24-dimensional half-swap survives: Monster valuations, Niemeier glue and Golay contacts

The useful question behind the repeated integers is not whether they occur, but **which operation survives the change of object**. Here the same test distinguishes a prime-valuation display from a group action, a root-system half-swap from a lattice isometry, a 26-site incidence polarity from a lattice map, and a Golay permutation from a Hamilton-count symmetry.

## 1. Exact Monster bookkeeping and a different, honest 13×24

For the [Monster](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/M/) and [Baby Monster](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/B), the enlarged common primes \(2,3,5,7,11,13\) have valuation differences \(5,7,3,4,1,2\), respectively. The unchanged common primes are \(17,19,23,31,47\); the new ones are \(29,41,59,71\). Thus
\[
 |M|/|B|=(2^5 3^7 5^3 7^4 11\,13^2)\,(29\cdot41\cdot59\cdot71).
\]
The [old valuation chain](../session32_square_monster_20261023/VALUATION_CYCLE_AND_QUOTIENT.md) is \(13\to2\to5\to3\to7\to4\) with separate \(11\to1\). Closing it by \(13=4\) changes its state space; it is not a congruence or a canonical group action. Moreover the true 2A centralizer is \(2.B\), so the index subtracts one from the exponent of 2 and changes the closed chain to a two-cycle. The distinct 2B centralizer is \(2^{1+24}.\mathrm{Co}_1\): there 24 is a dimension and 25 is the exponent of the order of the extra-special factor, not interchangeable readings. See [the Monster character-table verification](https://arxiv.org/abs/2412.12182).

There **is** an exact 13-by-24 group object on the Baby side: its unique order-13 class has normalizer \(S_4\times(13:12)\), so its centralizer is \(S_4\times C_{13}\), of order 312. This uses the normalizer and faithful 12-part action, not just the equality \(312=13\cdot24\). [Breuer's character-table reconstruction](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/htm/ctblbm.htm) records the structure and the centralizer. Its \(S_4\) has nine involutions. A different 24-group below has precisely one.

The [old level-two prime partition](../session32_square_monster_20261023/monster/MONSTER_BABY_PRIME_PARTITION.md) retains a complementary point: at each of \(29,41,59,71\), the coarse 2-isogeny read has one nonrational Frobenius pair of parallel edges even though their endpoints are rational. The 24 in its eta-quotient exponent and the additive 24 in the normalized 2A coordinate do not identify those edges. The detail that was forgotten is the isogeny, not its endpoint.

## 2. The six values are exactly the pairable pure-A systems, not exactly the liftable ones

In the [Niemeier list](https://www.mi.uni-koeln.de/~mgriffin/PDFs/moonshine.pdf), the pure-A rank-24 root systems are \((A_{h-1})^r\) with \(r(h-1)=24\). Consequently the possible \(h\) are \(2,3,4,5,7,9,13,25\). **Exactly** the archived set \(\{4,2,5,3,7,13\}\) has even \(r\), equivalently \(h-1\mid12\). So these are precisely the systems that admit a fixed-point-free component permutation built from disjoint swaps. The added \(9,25\) have three and one components and fail this test before considering the glue. This is a structural characterization of the *set*, not an explanation for the *order of the valuation arrows*. Replacing 24 by 25 or 26 in the root-rank identity is only a divisor calculation; positive even unimodular lattice ranks are divisible by eight.

The new distinction is exact:

| \(h\) | Components \(r\) | Minimum order of a Niemeier isometry inducing disjoint swaps on all components |
|---:|---:|---:|
| 2 | 24 | 2 |
| 3 | 12 | 4 |
| 4 | 8 | 2 |
| 5 | 6 | 4 |
| 7 | 4 | 4 |
| 9 | 3 | no paired component permutation |
| 13 | 2 | 4 |
| 25 | 1 | no paired component permutation |

**Obstruction at odd prime \(p=h\).** The root lattice's discriminant module is \(\mathbb F_p^{2m}\) for \(r=2m\), with bilinear form a nonzero scalar multiple of the dot product. The unimodular glue is a totally isotropic subspace of dimension \(m\). A root isometry exchanging all components acts by a signed permutation on this module (Weyl operations act trivially). If the isometry has order two, each swap has sign product \(+1\). Its two eigenspaces have dimension \(m\), with forms represented, up to a common scalar, by \(2I_m\). The invariant glue splits into its intersections with those eigenspaces. For odd \(m\) the dimensions of the two isotropic intersections sum to at most \(m-1\). For even \(m\) both must have dimension \(m/2\), so both forms have to be hyperbolic. This requires \((-1)^{m/2}\) to be a square in \(\mathbb F_p\), by the determinant of a hyperbolic form. The test fails for \(p=3,7\) (nonsquare \(-1\) with odd \(m/2\)) and for \(p=5,13\) (odd \(m\)). This rules out every order-two lift, not only the unsigned one.

**Existence and minimality.** The companion constructs each glue explicitly. For the four odd primes it verifies a signed swap squaring to \(-I\) on the whole lattice, hence of order four. The certificates enumerate the glue classes, verify isotropy and the self-dual size, and verify that the least nonzero glue-coset norm is four, so no extra roots appear. The binary Golay involution comes from \(x\mapsto-1/x\) on \(\mathbf P^1(\mathbb F_{23})\); the octacode gives an unsigned pairing at \(h=4\), checked on generators derived from its [quaternary cyclic description](https://arxiv.org/abs/math/0207208). These disjoint component swaps become actual order-two lattice isometries for \(h=2,4\).

**A sharper 24 comparison.** For \(A_6^4\) the glue can be generated by the rows of \([I_2\mid \left(\begin{smallmatrix}2&3\\4&2\end{smallmatrix}\right)]\) over \(\mathbb F_7\). Among the 4! component permutations, exactly the 12 even ones lift to a signed action on this code, with two lifts each and kernel \(\{\pm I\}\). The replay checks all signed permutations. Over the normal Klein four group of \(A_4\), the three nonidentity elements are disjoint double swaps, and their six lifts have order four by the obstruction just proved. The inverse image is \(Q_8\); adjoining an element of order three permutes its three order-four axes and produces the binary tetrahedral group of order 24. The Baby Monster's 13-centralizer instead has \(S_4\) as its 24-factor. **Equal size, different lift:** the two 24-groups have one and nine involutions respectively.

## 3. A 26-site incidence plane really produces the rank-24 root lattice

This is a literal link to the [Singer order-13 Moore block](../moore_fixed_order_spectra_20261124/SINGER-ORDER13-MODULES.md). Take its 13 points and 13 lines, with incidence matrix \(N\). In \(\mathbb Z^{26}\), impose zero sum separately on points and lines. This gives the real integral lattice \(R=A_{12}\oplus A_{12}\), of rank **24** and determinant \(13^2\). The incidence relation \(NN^{\mathsf T}=3I+J\), together with row and column sum 4, gives actual numerical maps on \(R\).

* The phase polarity \((P_x,L_y)\mapsto(L_{-x},P_{-y})\) exchanges the point and line coordinate sets and preserves incidence. On the discriminant module it is the *unsigned* swap. It does **not** preserve the Niemeier glue line \(b=5a\) in \(\mathbb F_{13}^2\), because it changes the slope to \(1/5=-5\). So it is an involution of the 26-site incidence graph, but not an involution of this glued lattice.
* The signed component swap preserves that line and lifts with square \(-I\). The **signed incidence map**
  \[
       S=\begin{pmatrix}0&N\\-N^{\mathsf T}&0\end{pmatrix}
  \]
  also preserves the Niemeier overlattice, and has \(S^2=-3I,\ S^{\mathsf T}S=3I\). On the discriminant it acts by \((a,b)\mapsto(4b,-4a)\), because \(N\) sends \(e_i-\mathbf1/13\) to \(4(e_j-\mathbf1/13)\) modulo the target root lattice. This is an integral similarity of norm multiplier 3, not a graph automorphism or a lattice isometry.

The neighbouring [Singer planes](https://personal.ntu.edu.sg/bernhard/Publications/pub/Singer.pdf) make a useful control. For \(v=q^2+q+1=7,13,21,57\), two zero-sum sides have rank \(2(v-1)\). Among these four, only \(v=13\) has an even unimodular gluing of these two A-components. At a prime \(\ell\equiv3\bmod4\) dividing \(v\) exactly once, the local two-dimensional discriminant form is a scalar multiple of \(x^2+y^2\), anisotropic, obstructing an isotropic subgroup of order \(\ell\). Use \(\ell=7\) for \(v=7\) and \(\ell=3\) for \(v=21,57\). For \(v=13\), the line of slope 5 is isotropic and has the required size. Thus the same data that let 26 split off two constant directions create a distinctive 24-dimensional glue, not a generic property of all projective-plane incidence graphs. No missing Moore graph is needed to build the incidence plane itself.

## 4. A Golay swap that preserves the 24-bit code but not the Hamilton count

In the repository's [physical 24-contact dodecad board](../chanpak_h23_golay_chains_20261124/DODECAD-CHAIN.md), choose the binary Golay involution \(g:x\mapsto-1/x\) just used. Retain all 24 orientations, not merely a weight. Its fixed code has dimension six. Indeed the replay finds 64 fixed words and verifies \(\mathrm{im}(1+g)=\ker(1+g)\) on the code. After merging paired coordinates the fixed code is length 12 and self-dual: the dot product of a fixed word with one representative from each pair of \((1+g)y\) equals the dot product of its full doubled word with \(y\). Its representative weights are \(0,4,6,8,12\) with multiplicities \(1,15,32,15,1\).

The twisted equation \(g(c)=c+\mathbf1\) has 64 solutions, each a transversal of the twelve pairs. Choose one, \(D=\texttt{0x5dfb}\), as the dodecad for the first twelve physical contacts. Then \(g\) exchanges \(D,D^c\). The whole twisted fibre is \(D+C^g\). Writing \((a,b)\) for the two contact weights, its exact profiles in ordered form are
\[
 (0,12):1,\quad(4,8):15,\quad(6,6):32,\quad(8,4):15,\quad(12,0):1.
\]
The physical formula reads \(H=3^{a+b}(4\cdot2^{2k-a-b}+3\cdot2^{k-a}+2^{k-b}-1)\), \(k\ge36\). Under \(g\) the height changes by
\[
 H(c)-H(gc)=2\cdot 3^{a+b}(2^{k-a}-2^{k-b}).
\]
There are exactly **five** actual Hamilton heights on this fibre, and precisely 32 words keep their count under the swap. The permutation of contacts is not an automorphism of any underlying tournament asserted here. Physically relabelling the two halves of the edit register need not preserve the two different ports.

The [older unequal-chain weight projection](../chanpak_h23_golay_chains_20261124/UNEQUAL-CHAINS.md) provides an independent precision test at **24,25,26**. Puncture the same full code to 22, 23 or 24 coordinates and read the weight \(w\) as lengths \((1+w,1+k-w)\) for \(k=22,23,24\). The sums of lengths are respectively **24,25,26**, yielding **11,8,5** distinct exact heights. These are length budgets, not tournament orders (four other vertices remain). The 25 case was already archived. Injectivity of the unequal-chain formula at fixed sum and the full weight enumerators certify the two neighbours. In the unpunctured 26-budget chart all 64 twisted words collapse to weight 12 and hence one height **268,468,223** at lengths (13,13). The five physical heights versus that single projected height is an exact measure of lost placement information.

## 5. Archive index: other meanings worth retaining separately

| Prior lane | Which invariant survives, and which does not |
|---|---|
| [June involution modulus](../../../01-canon/theorems/THM-484-involution-modulus-24.md) | 24 is the largest modulus whose units square to 1; the Paley/Gleason read uses 8 and 24. This is not an order-two lift theorem for the Niemeier glue. |
| [June Hurwitz–Markov–Pell cannonball](../../../07-reflections/lrc14-hurwitz-markov-pell-cannonball-codex-s243.md) | \(29\cdot169-70^2=1\) and \(1^2+\cdots+24^2=70^2\) are different exact equations sharing a coordinate. In the [Borcherds/Conway cusp](https://archive.intlpress.com/site/pub/files/_fulltext/journals/cdm/1996/1996/0001/CDM-1996-1996-0001-a001.pdf) there are 25 spatial slots \(0,\ldots,24\) and one time slot 70; the null reduction removes two directions. The repo's [25-score chart](../../../01-canon/theorems/THM-4752-tournament-scores-in-ii-25-1-backward-length-and-the-a24-deep-hole.md) measures cyclic triangles, not all Hamilton paths. |
| [October deep-hole swap](../leech_deep_hole_swap_theorem_20261022.md) | The local Weyl-point and cusp comparisons use their own theta readouts. Neither a root label nor a shell number replaces the isotropic-glue test. |
| [October March-to-Golay Baer lane](../../../04-computation/golay_baer_bridge_88.py) and [rank-24 Paley/Golay experiment](../../../04-computation/experiments/paley_golay_leech_moonshine_20261018.py) | Numbering sheets or counts is not attaching a physical symmetry. The explicit projective involution and dodecad above supply one actual contact permutation, together with exactly what it fails to preserve. |
| [October 25 ledger](../analogy_reduction_atlas_20261025/LEDGER_15_25_MULTIPERFECT.md) | Cyclic versus elementary-abelian order-25 tournaments, the 25-score cusp and the 25-state Goldbach machine have distinct definitions; a shared construction size is not an isomorphism. |
| [November native birth at 26](../eleven_wedge_honeycomb_20261102/NATIVE_26.md) | New square-sum ports are exactly 10–26 and 23–26; deleting internal 26 exposes a two-path forest, not a 25-path. This is another marked-interface effect, independent of the Golay puncture. |
| [November marked 24:13](../marked_24_13_descent_wheel_20261116/ARCHIVE.md) and [13-power selector](../moore_fixed_order_spectra_20261124/THIRTEEN-ODD-ALPHABET.md) | \(12\mid24\) in the order-13 wheel, Bernoulli denominators and the \(A_{12}^2\) root system. The wheel gives modular selection, not arbitrary absolute counts; the new norm-3 incidence similarity is a separate actual lattice map. |
| [November all eight signed pure-A labels](../chanpak_absolute_liouville_leech_20261122/LEECH-LABELS.md) | The six labels omit 9 and 25 in the valuation display and in the root-half pairing, but a separate port construction realizes their height rails. This separates inability to swap components from inability to realize a count. |

The concrete consequence for future analogy searches is to preserve the right boundary data: **prime centralizer structure, glue subgroup, or ordered two-port profile**. Each can survive a symmetry lost at a different level. The companion is dependency-free and checks all finite certificates above.
