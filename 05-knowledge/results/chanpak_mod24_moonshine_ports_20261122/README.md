# The two endpoints of a flipped arc control its switching table

**PROVED.** The signed change of the Hamiltonian path count under a specified arc reversal has zero mixed difference in the two endpoint switches. Away from those endpoints, marking the arc recovers a path–cycle identity at the first two orders where unrestricted interpolation allows two distinct component counts. The exact rational observation dimension at order nine is 182.

**FINITE-EXACT.** Within one nine-vertex switching class and its one-arc neighbors there are 66 marked gaps of two, at 29 integer centers. The centers cover every even residue modulo 24 using four physical arcs, and four are necessary. The original switching class has trivial relabeling stabilizer.

## Two commuting changes

Let $T$ be a labelled tournament on $n\ge3$ vertices, with skew edge signs $x$ in a fixed orientation gauge. A cut with vertex signs $s_i\in\{1,-1\}$ replaces $x_{ij}$ by $s_i s_j x_{ij}$. Write $T^s$ for the result. Reversing an unordered pair $e=\{a,b\}$ is denoted $T^e$; in the formulas below the superscript is typed by its argument. For the directed Hamiltonian path count $H$ or cycle count $C$ (up to rotation) put
$$
 f_{e,H}(s)=H((T^e)^s)-H(T^s),\qquad f_{e,C}(s)=C((T^e)^s)-C(T^s).
$$
The two changes commute. For a set $U$ of switch variables, $\Delta_U$ is the ordinary mixed Boolean difference about a specified cut. Global negation of all spins describes the same cut. The object is a *signed* difference on a physical arc, before taking its absolute value or identifying isomorphic tournaments.

The closest prior mechanisms are the pair-cover identity for one arc and the even-path forest expansion for the path potential. The source-independent divisibility law for cut switches provides the odd-prime comparison and the radius-four modulo-24 reconstruction. A competing scalar model which permits both endpoints in an independent cut coefficient will fail the integer identity below. [THM-5267, arc-flip formula](../../../01-canon/theorems/THM-5267-seams-at-a-prescribed-centre-by-one-port-contact-the-arc-flip-lemma-deletion-monotonicity-the-order-10-extension-graph-and-the-heaviest-parent-tree-of-first-creators.md) · [THM-4761, Walsh forest expansion](../../../01-canon/theorems/THM-4761-random-tournaments-the-forest-identity-is-the-walsh-chaos-decomposition-of-kappa.md) · [THM-5313, switching derivatives, sections 2 and 6](../../../01-canon/theorems/THM-5313-switching-derivatives-of-hamiltonian-paths-and-cycles.md)

## 1. Endpoint restrictions and the exact odd-prime order

**Theorem.** Let $U\ne\varnothing$, $d=|U|$, $r=\lceil d/2\rceil$, and let $b_0$ be the least integer at least $r$ of the same parity as $n$. If $U$ contains both endpoints of $e$, then $\Delta_U f_{e,H}=\Delta_U f_{e,C}=0$ identically. In all other cases define
$$
 L=n-2r-2\mathbf1_{\{d=2r,\ e\cap U=\varnothing\}}.
$$
If $L<b_0$ the two differences vanish identically. Otherwise, for every odd prime $p$ and any fixed $e,U$ of the indicated type,
$$
 \min_T v_p(\Delta_U f_{e,H})=v_p(b_0!),\qquad
 \min_T v_p(\Delta_U f_{e,C})=v_p((b_0-1)!).
$$
For integer $a_0$ the corresponding minimum on $f_{e,H}-2a_0 f_{e,C}$ is the minimum of $v_p((c-a_0)(c-1)!)$ over all $b_0\le c\le L$ with $c\equiv n\pmod2$. Zero has infinite valuation. In particular for every fixed choice of the other switches the four signed seam values at the two endpoints satisfy $f_{00}+f_{11}=f_{01}+f_{10}$ over the integers.

**Proof.** A nonconstant cut character of either observable comes from a linear forest $F$ with even-length paths, and has switch support the set $\partial F$ of endpoints of its nontrivial components. A spanning cycle of the cycle observable has empty switch support. Flipping $e$ and differencing in $U$ leaves exactly those forests with $e\in F$ and $U\subseteq\partial F$. Their edge-character coefficients acquire an extra factor $-2(-2)^d$. The path and cycle coefficients before this factor are signed $2^{j-n+1}c!$ and $2^{j-n}(c-1)!$ with the same sign, for $j$ nontrivial paths and $c$ total components. Distinct forests retain distinct characters of the independent source signs. [THM-5313, section 6](../../../01-canon/theorems/THM-5313-switching-derivatives-of-hamiltonian-paths-and-cycles.md)

No edge in an even-length path has both ends as leaves. Always $j\ge r$, $c\ge b_0$, $n-c\ge2j$, and $c\equiv n\pmod2$. If $d=2r$ and $U$ misses the marked pair, then either $j>r$ or all forest endpoints are required. In the latter case the marked edge joins two interior vertices, requiring a path of length at least four. In both cases $n-c\ge2r+2$.

Conversely take any listed $c$. Start with $r$ two-edge paths and $c-r$ isolates. Add the even number $n-2r-c$ of remaining vertices to one path, in pairs. If one endpoint of $e$ is required, position $e$ as the first edge of a path. If $d$ is odd and $e$ is disjoint from $U$, place one endpoint of $e$ at the unused leaf and the other at its neighbor. If $d$ is even and $e$ is disjoint, there are at least two extra vertices, so position $e$ between adjacent interior vertices of the lengthened path. In all cases relabel the remaining leaves to include $U$. This realizes every allowed $c$. Walsh inversion over the labelled source signs uses denominators that are powers of two. Over the localization at $p$ the ideals generated by the function values and by its coefficients therefore coincide. This proves the minima.

**A restored cancellation.** If $d=2r$ and the marked arc misses $U$, the orders $n=3r,3r+1$ give zero. At $n=3r+2,3r+3$ exactly one component count, $c=b_0$, remains. Consequently
$$
 \Delta_U f_{e,H}=2b_0\,\Delta_U f_{e,C}
$$
over the integers for every source. At $n=3r+2$ both $c=r$ and $r+2$ contribute to the unmarked derivative. Specifying the flipped edge outside the switched vertices removes the higher component count. At nine vertices this forces the quartic part of the seam of $H-6C$ to be supported on quartets containing exactly one endpoint.

## 2. How many cut samples determine an entire signed seam?

At nine vertices fix a gauge vertex outside the marked pair. On the eight free coordinates, all cuts of weight at most four determine $H$ modulo 24. Omitting cuts containing both endpoints therefore reconstructs the seam modulo 24 from $163-(1+6+15)=141$ signed observations. For integer reconstruction there are 247 cuts of weight at most six. Omit the 57 containing both endpoints, and the single six-cut avoiding both, leaving 189. These deletion rules form downsets, so the retained Boolean coefficients are obtained by ordinary inclusion-exclusion on the retained samples. The missing coefficients vanish in the stated ring.

**Sharp rational-linear dimension.** Over all labelled sources, the dimension of the span of the integer seam tables for one fixed physical edge, as functions of the cut, is
$$
 \sum_{t=1}^{\lfloor n/3\rfloor}
 \left(2\binom{n-2}{2t-1}+
   \mathbf1_{n\ge3t+2}\binom{n-2}{2t}\right).
$$
At $n=9$ the contributions are $35,105,42$, giving 182. Thus a source-independent set of 182 cuts and rational linear interpolation suffice. Fewer observations cannot give a universal *linear* reconstruction from cut evaluations.

Indeed the surviving even spin characters have endpoint sets of size $2t$. The first term counts those touching exactly one marked endpoint. The second counts those avoiding both; they require one longer path as in the proof with $j=t$. Two distinct endpoint sets have disjoint collections of coefficients in the *source* edge characters, each nonempty. Thus the associated coefficient functions of the source are linearly independent. The spin characters themselves are independent on cuts, since exactly one representative of the simultaneous sign quotient is chosen and the supports all have even size. This gives equality in the dimension formula. A maximal independent set of evaluation rows recovers the functions. The probe selects such rows over a finite auxiliary prime, certifying full rank over the rationals. The dimension bound addresses exact linear interpolation; it supplies no divisibility criterion on nonlinear decoders.

## 3. Marked contacts and the Monster comparison

Let $G$ be the group of relabelings preserving a labelled switching class. It acts on pairs $(X,e)$ consisting of a cut output and a marked unordered pair. For a given pair of ordered heights $h=H(X)$ and $h'=H(X^e)$, its stabilizer acts freely on the Hamiltonian paths of each tournament, because an ordered spanning path has no nonidentity label symmetry. Its order therefore divides $\gcd(h,h')$. In particular, the action is free on each exact gap-two center: both heights are odd and differ by two. The number of marked contacts in each such center is a multiple of $|G|$. This includes restrictions to $G$-invariant sets of marks or additional conditions preserved by $G$. Forgetting the mark generally enlarges the stabilizer. [Earlier sheet conservation and orbit weights](../chanpak_moonshine_prime_blocks_20261118/SWITCHING_SHEET_CONSERVATION.md)

In the Monster, retaining an involution of class $2A$ selects its centralizer $2.B$. Retaining a tournament and its gap-two arc selects a trivial joint stabilizer. These examples share an orbit-stabilizer operation and have different stabilizer groups. The Monster and the Baby Monster are also excluded as faithful groups of relabelings of a single tournament switching class at **any** order. Babai and Cameron's classification requires cyclic or dihedral Sylow 2-subgroups. The ATLAS gives extraspecial subgroups of orders $2^{25}$ and $2^{23}$ respectively. Each has exponent at most four, since every square lies in its center of order two, so neither can be cyclic or dihedral of its order. Subgroups of cyclic or dihedral 2-groups are cyclic or dihedral. This obstruction is specific to the stated action on a switching class. [Babai and Cameron, *Automorphisms and Enumeration of Switching Classes of Tournaments*](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v7i1r38) · [ATLAS Monster subgroups](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/M/) · [ATLAS Baby Monster subgroups](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/B/) · [Prior marked-involution and prime-support comparison](../session32_square_monster_20261023/monster/MONSTER_BABY_PRIME_PARTITION.md)

In the earlier Paley construction, 24 is the length of the extended binary Golay code and the rank of the resulting Leech lattice. Its 196560 minimal vectors enter the weight-two lattice vertex-algebra count $196560+24+300=196884$. In the Monster moonshine vertex operator algebra, an Ising vector and its commuting conformal complement split the **central charge** as $24=1/2+47/2$. The commutant is the Baby Monster vertex operator algebra, the even part of Höhn’s shorter moonshine superalgebra. For the marked tournament data the operation is the Chinese remainder isomorphism. A count or a signed difference modulo 24 is uniquely specified by its ordered pair of residues modulo 8 and modulo 3. [Prior lattice construction](../chanpak_moonshine_prime_blocks_20261118/MASS_TRANSPORT.md) · [Höhn, Lam and Yamauchi, section 5.1](https://www.math.ksu.edu/~gerald/papers/e7Aug05.pdf) · [Höhn’s shorter moonshine construction](https://arxiv.org/abs/0706.0236)

## 4. One exact nine-vertex experiment

Use the fixture $\mathtt{0x7103d59dc}$ of THM-5313, with lexicographically ordered pairs on $0,\ldots,8$ and bit one directing the edge from the smaller endpoint to the larger. Inspect all 256 cuts and all 36 physical arcs. Prefix and suffix path dynamic programs compute both path-cover readings across each possible arc. Three direct flipped-tournament recounts control sign conventions; cut interpolation is checked against all 9,216 readings.

The 66 absolute-gap-two contacts span 29 centers, representing all 12 even residues modulo 24. A minimum portfolio of four marked arcs covers every phase.

| Mark | Center residues modulo 24 |
|---|---|
| $(1,2)$ | $4,10,14,16,18,22$ |
| $(1,7)$ | $0,2,6,20$ |
| $(2,3)$ | $2,8$ |
| $(2,4)$ | $12$ |

The exact minimum comes from subset dynamic programming on the 12 residues with all 36 arc supports. The group's free action on each center gives $29|G|\le66$. Its order is odd at odd vertex order: any involution would fix a vertex, which can be normalized to a source by switching, and the induced root stabilizer acts by ordinary tournament automorphisms on its complement. A nonidentity ordinary tournament involution would reverse the arc of a transposed pair. Therefore $G=1$. The finite source realizes all odd phases modulo 70; at modulus 48 it misses precisely 41. A phase of a marked contact certifies the existence of a particular integer center, and lifting it to all integers of that phase remains open.

Run `python3 05-knowledge/results/chanpak_mod24_moonshine_ports_20261122/check.py` from the repository root. The [script](check.py) prints the complete list of centers, residues and successful even moduli, together with the sample counts and rank. [Exact output](check.out)
