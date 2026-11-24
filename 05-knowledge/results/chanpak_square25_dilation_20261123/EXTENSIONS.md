# Twelve marked tournaments suffice at every order from eleven onward

**PROVED.** There are twelve fixed, marked, strongly connected tournaments on eleven vertices. For any prescribed order \(n\ge11\), replace the mark of each by a transitive tournament on \(w=n-10\) vertices. The resulting twelve tournaments have transitive odd Hamilton-path precedence, are rigid, and realize every odd residue modulo 24 exactly once. Both eleven and twelve are sharp here. The inherited lower bound excludes order ten, and each marked template has just one output at a fixed order. This is a residue theorem. The formulas below also give the exact integers its selector produces. [Eleven-vertex theorem and lower bound](STRONG11.md)

## Counting the visits to one vertex

For a tournament \(T\) on \(d\) vertices with mark \(v\), erase the labels of the transitive block that replaces \(v\) in a Hamilton path and contract each maximal run of that block to one \(v\). Let \(c_k\) count the resulting possible words with exactly \(k\) visits. Equivalently, all exterior vertices occur once, the mark occurs \(k\) times, and every adjacent pair is an arc of \(T\). In particular \(1\le k\le d\), and \(c_1=H(T)\).

**Visit formula.** For every integer \(w\ge1\) the expanded tournament has
\[
 H_w=\sum_{k=1}^d c_k\,k!{w\brace k}.
\]
Indeed each surjection of the \(w\) labelled block vertices onto the \(k\) ordered visits gives exactly one filling. Within each visit the labels must appear in transitive order. Every allowed path arises uniquely from its word and this surjection. This is the inherited one-mark interpolation formula. [Exact reader interpolation](../chanpak_moonshine_prime_blocks_20261118/EXACT_READER_INTERPOLATION.md)

**Cycle-packing formula.** Take the set of simple directed cycles through \(v\), each rooted at \(v\). Let \(\mathcal P_s\) be the unordered packs of \(s\) such cycles disjoint outside the mark. Write \(U(P)\) for their combined exterior vertices. Then
\[
 c_{s+1}=s!\sum_{P\in\mathcal P_s} H\bigl(T-U(P)\bigr).
\]
Between consecutive visits a run word makes a simple directed cycle through \(v\). Its exterior vertices occur nowhere else. Delete these interiors and contract the visits to a single \(v\) to obtain a Hamilton path of the residual tournament. Conversely insert the cycles at the mark of that path, in any order. This proves the formula. In particular all \(c_k\) are divisible by \((k-1)!\), and \(c_2\bmod2\) is the parity of the number of simple directed cycles through \(v\), since every residual height is odd.

**The growth base is a matching number.** Build a bipartite graph \(G_v\) between the out-neighbors and the in-neighbors of \(v\), retaining each pair \((a,b)\) with \(v\to a\to b\to v\). Then
\[
 \max\{k:c_k>0\}=1+\nu(G_v),\qquad
 \lim_{w\to\infty}\frac{H_{w+1}}{H_w}=1+\nu(G_v),
\]
where \(\nu\) is the maximum matching size. Indeed each cycle through \(v\) takes a first step in its out-neighbors and a last step in its in-neighbors. The first transition between the two gives a triangle through \(v\) supported on that cycle. Disjoint cycles therefore yield a matching of equal size. A matching gives a pack of triangles, and the residual tournament always has a Hamilton path. The visit formula has leading term \(c_b b^w\) for \(b=1+\nu(G_v)\). Flipping one arc changes the growth base by at most one: for an exterior arc at most one bipartite edge changes, and for an incident arc deleting the other endpoint gives the same bipartite graph before and after. In this selector eight marks have matching number one and four have matching number two. The equality for the maximum number of visits is the existing marked-visit theorem; the pack formula identifies its exponential rate and its response to a single flip. [Marked-visit matching](../chanpak_moonshine_prime_blocks_20261118/MARKED_VISIT_MATCHING.md)

Modulo 24 only \(k=1,2,3\) contribute. Inclusion-exclusion for surjections gives
\[
 H_w\equiv c_1+(2^w-2)c_2+(3^w-3\cdot2^w+3)c_3\pmod{24}.
\]
For every marked tournament the residues are therefore periodic with period dividing two from \(w=3\). The twelve templates below have even \(c_2,c_3\), so \(H_2\equiv H_4\) as well. Their whole order selector uses width one, even width, and odd width at least three. For an arbitrary mark the period-two law already holds from width two exactly when the number of directed cycles through the mark is even. Indeed \(H_4-H_2\equiv12(c_2+3c_3)\equiv12c_2\pmod{24}\). The eventual residues are constant exactly when \(3\mid c_2\) and \(4\mid c_3\), since \(H_4-H_3\equiv8c_2+6c_3\pmod{24}\).

**All moduli.** More generally, write \(M=\prod p^{a_p}>1\), set \(A=\max a_p\), and let \(L=\lambda(M)\) be the Carmichael exponent of its unit group. For every fixed marked tournament,
\[
 H_{w+L}\equiv H_w\pmod M\qquad(w\ge A).
\]
Expand each surjection count into integer multiples of \(j^w\) using inclusion-exclusion. For a given \(p^{a_p}\), a nonunit \(j\) divisible by \(p\) has \(j^w=0\) modulo this power from \(w=a_p\), and a unit repeats after \(L\). The Chinese remainder theorem gives the result. Terms with \(M\mid k!\) may be discarded before forming the sum.

**Exact transfer test.** For any two marked tournament templates of orders at most \(K\), the exact difference \(H'_w-H_w\) equals a prescribed constant \(g\) at every positive width precisely when \(c'_1-c_1=g\) and \(c'_k=c_k\) for \(2\le k\le K\), taking missing coordinates as zero. Checking widths \(1,\ldots,K\) suffices. The matrix of surjection counts is lower triangular with diagonal \(w!\). If any higher coordinate differs, let \(b\) be the largest such index. The leading term of the difference is \((c'_b-c_b)b^w\), so even bounded differences along infinitely many widths force the constant case. For two templates differing in one exterior arc, every inflated pair still differs in precisely that arc. The test distinguishes a persistent exact flip from a modular coincidence.

## Four contact boards and twelve fixed marks

The following four tournaments have vertices \(0,\ldots,9\). Entry \(i\) is the outgoing bitmask at \(i\).

| Board | Outgoing masks |
|---|---|
| A | 766, 1020, 1016, 496, 992, 960, 896, 768, 513, 8 |
| B | 542, 1020, 1016, 1008, 864, 705, 385, 785, 33, 320 |
| C | 750, 924, 1016, 1008, 737, 322, 898, 800, 17, 288 |
| D | 1022, 1004, 216, 880, 930, 964, 912, 776, 516, 4 |

Adjoin vertex 10 with the outgoing mask shown, then mark the indicated old vertex. The three rightmost columns list residues. At the desired width, select the row whose corresponding column contains the target modulo 24. The exact value then follows by using the displayed coefficients in the visit formula. All \(c_k\) with \(k\ge4\) vanish, so the same three-term formula above even holds in the integers on each row.

| Board | Mask | Mark | \(c_1\) | \(c_2\) | \(c_3\) | \(w=1\) | even \(w\) | odd \(w\ge3\) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 512 | 1 | 1489 | 712 | 0 | 1 | 9 | 1 |
| D | 1 | 8 | 15339 | 15976 | 5382 | 3 | 11 | 15 |
| A | 1008 | 1 | 2261 | 1098 | 0 | 5 | 17 | 17 |
| D | 773 | 8 | 7783 | 4210 | 0 | 7 | 3 | 19 |
| A | 962 | 2 | 5625 | 4466 | 1072 | 9 | 13 | 21 |
| D | 1021 | 8 | 1859 | 1026 | 0 | 11 | 23 | 23 |
| A | 128 | 1 | 1525 | 724 | 0 | 13 | 21 | 13 |
| B | 8 | 1 | 9327 | 4982 | 0 | 15 | 19 | 3 |
| A | 1009 | 1 | 2345 | 1230 | 0 | 17 | 5 | 5 |
| C | 32 | 7 | 7147 | 6390 | 1872 | 19 | 7 | 7 |
| A | 993 | 1 | 2541 | 1298 | 0 | 21 | 1 | 9 |
| B | 1014 | 1 | 6191 | 5156 | 1218 | 23 | 15 | 11 |

**Why the graph properties persist.** For the eleven-vertex parent, let \(D(T)\) order a pair by which direction of precedence occurs in an odd number of Hamilton paths. The replay checks that this is a total order, the parent is strong, and all exterior in-neighbors of its mark precede all exterior out-neighbors in that order. The substitution law keeps the exterior order and makes each block vertex see the exterior exactly as the original mark did. A two-vertex block has one precedence arc in one of the two directions. A larger transitive block keeps its internal order. Thus the module inserts at the indicated cut for every width. Its precedence is transitive, so every automorphism fixes each vertex. Strong connectivity lifts by following paths in the original quotient; within a block one can return via an exterior directed cycle through the mark. [Substitution law](STRONG-SELECTOR.md#what-survives-substitution)

**FINITE-EXACT certificates.** The replay verifies the twelve parents and every run coefficient through \(k=11\). It computes precedence from parities of path prefixes and suffixes, and counts run words by their last vertex, exterior subset and number of visits. Each of the three selector columns is a permutation of the odd residues. Hence the finite table, visit formula, and order-cut argument prove the all-order assertion. Run `python3 05-knowledge/results/chanpak_square25_dilation_20261123/extensions.py`. [Replay and integer selector](extensions.py)

## A dihedral group of order 24 on the twelve residue addresses

**PROVED for the displayed bank.** Name each marked template by its width-one residue. The two other width choices give permutations of these names: let \(E\) give the residue at an even width and \(O\) the residue at an odd width at least three. The actual contact table gives
\[
\begin{aligned}
E&=(1\ 9\ 13\ 21)(3\ 11\ 23\ 15\ 19\ 7)(5\ 17),\\
O&=(3\ 15)(5\ 17)(7\ 19)(9\ 21)(11\ 23).
\end{aligned}
\]
Thus \(E\) has order 12, \(O^2=1\), and \(O E O=E^{-1}\). Since \(E\ne E^{-1}\), the involution lies outside its rotation subgroup. The generated group is dihedral of **order 24**, acting faithfully with orbit sizes 4, 6 and 2. In the split/inert terminology for odd integers, the 6-cycle is the entire inert residue class \(3\bmod4\); the split class \(1\bmod4\) divides into the 4-cycle and the 2-cycle.

The reflection has two fixed vertices on the split four-cycle and no fixed vertex on the inert six-cycle. With the displayed cyclic labellings it acts by \(j\mapsto-j\bmod4\) and \(j\mapsto3-j\bmod6\), respectively. Changing either starting vertex changes the reflection constant by an even number. Thus the difference between a vertex-axis and an edge-axis reflection survives any cyclic choice of viewpoint.

There is an arithmetic reason for some of this symmetry. At width two the difference from width one is \(2c_2\), a multiple of four on these marks. At width three it is \(6(c_2+c_3)\), a multiple of twelve. Hence the odd phase must permute each antipodal pair \(\{r,r+12\}\) internally and so is automatically an involution in any complete bank with even \(c_2,c_3\). The stronger reflection identity and the three particular rotation lengths come from these chosen contacts. Composing the permutations means using a resulting residue to choose the next template by its width-one name. It does not posit an order-24 automorphism group of any of the underlying rigid tournaments.
