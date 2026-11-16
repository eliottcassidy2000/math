# 1729 as a target-decoder test, not a numerical analogy

The earlier [cubic lens](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/07-reflections/the-cubic-lens-where-1729-is-a-theorem-and-where-it-is-a-coincidence.md) separates the actual Eisenstein-norm mechanism from the isolated Paley path-orbit count 1729. The [two-cube divisor theorem](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-463-two-cube-representations-are-a-divisor-property-on-the-split-axis.md) works backward: for target $N$, try $d=x+y$ dividing $N$, then require $4N/d-d^2=3(x-y)^2$ and legal recovered summands. At 1729 the two admitted sums are $12+1=13$ and $10+9=19$. The useful lesson is the order of operations: *factor the target first, then reconstruct and verify the actual inputs*. Equal scalar outputs do not identify which construction remains usable.

## A literal divisor test for asymmetric rooted bridges

Use the [single reversed seam](../BRIDGE_RESET.md) from rooted $A$ to the converse of rooted $B$. In each input assume normalized signed response $O_i-E_i=g_i$ and set
$$
B_i=E_i+O_i,\quad \Delta_i=h_i-B_i+2g_i,\quad c_i=2g_i-\Delta_i,
\quad x_i=h_i+c_i/2.
$$
The input counts $h_i$ and deletion counts $g_i$ are odd, so $B_i$ is odd and $c_i$ even. The scalar bridge formula reads
$$
N=h_Ah_B+B_AB_B-g_Ag_B
 =2x_Ax_B+C,\qquad C=c_Ac_B/2-g_Ag_B.
$$
In the first equality the right end-profile is the left start-profile of $B$ read in converse, so the surplus is $2(O_B E_A+E_B O_A)=B_AB_B-g_Ag_B$. In the second equality substitute $B_i=h_i+c_i$. This is an exact target-first factorization for prescribed $(g_i,\Delta_i)$. If $(N-C)/2$ is not a positive integer there are no candidates. Otherwise each divisor pair gives $x_A,x_B$ and exactly one numerical candidate $(h_A,h_B)$; retain it only if actual physical inputs with both prescribed rooted profiles exist. In the ready sector $\Delta_A=\Delta_B=0$ this simplifies to
$$
N-g_Ag_B=2(h_A+g_A)(h_B+g_B),\quad N\equiv g_Ag_B\pmod8.
$$

There is a direct continuation filter too. If the next root is the bridge's right root and $B$ has more than one vertex, its deletion count is $g_{\rm out}=g_Bh_A$. Therefore the desired next deletion count chooses divisor candidates for $g_B$ and fixes $h_A$ before the scalar factorization. To be normalized there, the full start count at the root of $B$ must equal $g_B$: the new signed coordinate is minus that start count times $h_A$. For the unsigned output, keep the weighted root-end sum of $B$ and $P_{B-r}(2)$, as specified in the right-reroot transfer. These checks are the analog of reconstructing the cube summands from $d$ and their difference. Factoring alone does not furnish a tournament or a next fresh digit.

## Mirroring collapses to a norm and fails real targets

For $A=B$ and common defect $\Delta$, the formula completes the square:
$$
N+\Delta^2/2=(g-\Delta)^2+2(h+g-\Delta/2)^2.
$$
A fixed defect can thus be screened by the norm $a^2+2b^2$ and signs, then by physical profile reconstruction. All mirrored outputs have $N\equiv1\pmod8$ since $N=h^2+B^2-g^2$ and all three numbers are odd. A prime $p\equiv5,7\pmod8$ (hence inert for this norm) must divide $N+\Delta^2/2$ to even valuation: a single factor forces both coordinates divisible by $p$, and iteration proves evenness.

At $N=1729=7\cdot13\cdot19$, the ready mirror with $\Delta=0$ is ruled out by both 7 and 13. In fact no normalized mirror, at **any** defect or order, gives 1729. The [terminal-root lemma](../quadruple_splice/DECODER.md) forces the root into the last component and a prefix count $p$ multiplies all four coordinates. A mirror of a prefixed source has output divisible by $p^2$. Squarefreeness of 1729 leaves only removable singleton prefixes. A terminal singleton then yields output one. For a strong nontrivial source, $B\ge g$ gives $h^2\le1729$, so $h\le41$. The [Busch bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) gives a minimum of 45 at order eight, so the strong source has at most seven vertices. **mirror1729.cpp** enumerates every fixed-Hamilton-path orientation to that order with count at most 41, using positive induced monotonicity to prune, verifies each rooted profile by subset DP, and finds none. This is a restriction on a bridge, not on the absolute tournament spectrum.

Nor is the small-defect obstruction peculiar to 1729. On the old three-count row $756\cdot2^s-33+\{0,2,4\}$, the outer two counts are 7 and 3 modulo 8, so cannot be mirrored. The middle $N_s=756\cdot2^s-31$ meets the residue 1 but still fails every fixed mirror defect $0,2,4$ on the infinite progression $s=255\pmod{387660}$. For $p=71,13,5$ respectively, $2^{387660}\equiv1\pmod{p^2}$, and at $s=255$ the residues of $N_s+\Delta^2/2$ modulo those squares are 3124, 26, and 10; each has valuation exactly one at its indicated inert prime. These are fixed modular exponentiations, not prime-distribution assumptions. A universal selector must change the architecture, not just enlarge a bank of those symmetric choices.
