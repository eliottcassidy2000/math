# Cyclic clone growth: visible poles and the minimum literal counterexample

Follow-up to [Clone tests](CLONE_MARKS.md). The previous finite test allowed an exponential total width for cyclic marks. We can use a quadratic width instead, extract a smallest witness directly, and often reject a proposed switch without testing any widths. This is an obstruction and a certificate, not an all-prime construction.

## 1. A switch may only rotate the marked strong components

**Pole obstruction.** If two switching frames have equal q-slot cover counts for every independent tuple of nonempty replacement blocks, then on the marked subtournament the cut must consist of a prefix or suffix of its ordered strongly connected components (possibly empty or full). Here q can be any one fixed positive integer. In particular, if the marked subtournament is strongly connected, a switch reversing even one marked-to-marked arc cannot pass the universal height test. Ordinary, unmarked vertices need not follow this condition.

Here is the proof. Put a variable at every marked color. Before requiring every color, the compressed word series is rational by the stretch resolvent in [Clone tests](CLONE_MARKS.md). Requiring all colors takes inclusion–exclusion over subsets of the marked vertices. A denominator for the full-support series is thus a product of powers of determinants for induced restrictions of the marked tournament. A determinant splits over strong components.

For a strongly connected component C of at least three vertices, write D_C = det(I - A_C X_C). This polynomial is irreducible over the rationals and its variable support is exactly C. Indeed, it is multiaffine with constant 1, so a nontrivial factorization would partition its variables between the factors. But in the formal identity
```
-log det(I - A_C X_C) = sum_(k>=1) tr((A_C X_C)^k)/k
```
all coefficients on the right are nonnegative, and a closed walk crossing any proposed partition contributes a mixed monomial. Strong connectivity provides such a walk, contradicting the absence of mixed monomials in the logarithm of the factored expression. The same argument with a closed walk visiting any chosen vertex shows it really is a variable of D_C.

Every such **maximal** C in the marked quotient gives a genuine pole, even after requiring all marks. Fix a Hamiltonian path in the quotient and insert any closed marked walk starting and ending at a chosen x in C into the position of x. This gives a subset of the full-support words in one slot, with the other slots empty. Set the variables of C to u, the others to a fixed sufficiently small positive number, and let u rise to 1/rho(A_C). The sum over closed walks at x diverges by the positive Perron left and right eigenvectors (the diagonal resolvent has a nonzero residue at the simple positive Perron eigenvalue). All coefficients of the full-support series are nonnegative. At this radial limit, determinants from proper subsets of C stay nonzero, since proper principal submatrices of an irreducible nonnegative matrix have strictly smaller spectral radius. Determinants from the other components are nonzero by the small-variable choice. If the factor D_C canceled, the rational series would therefore remain bounded at the limit, which is impossible.

Universal equality of physical counts implies equality of the full-support rational series: for independent transitive widths the Stirling transform is triangular with nonzero diagonal at each finite total degree. Thus every maximal nontrivial marked component C of one quotient must occur as the variable support of a pole in the other. Every pole in the other has support contained in one of its marked strong components, even after vertex restrictions. Apply the containment both ways. Since the component partitions are disjoint, they coincide, including singleton components. On each common nontrivial component the actual determinant must also coincide, since the only restriction with that full component support includes the whole component.

The coefficient of the monomial on three marks in its determinant is -1 for a directed triangle and 0 otherwise. A triangle crossing the switch cut becomes transitive. Every strongly connected tournament is connected by chains of overlapping directed triangles: if we group vertices by that overlap relation, an outside vertex has uniform orientation to each triangle (nonuniform orientation would form a crossing triangle), hence to the whole group. The groups are themselves strong, have uniform intergroup arcs, and their quotient is triangle-free, hence transitive. Consequently the cut must be constant on each marked strong component. The directed order of the components then has no pattern 010 or 101 in the cut's membership word. Such a pattern would become a new crossing directed triangle, contradicting equality of the component partitions. Equivalently the membership word changes at most once. Switching on one such boundary rotates the two parts of the ordered list and preserves the pole signature. Preserving the pole signature is necessary, not sufficient, for neutral replacement. With k marked strong components and e ordinary vertices, only k·2^e cut masks up to global complement pass this routing test, rather than 2^(s+e-1). In particular, if every vertex of a strong tournament is marked, no nontrivial cut can agree at even one slot count across all independent replacements.

## 2. A quadratic cutoff, sharpened by common poles

Write s for the number of marks, g for the number of ordinary vertices plus the number q of slots, and F_A = P_A / D_A^g for the rational word series before requiring every color. The denominator D_A has degree at most s and constant 1; the numerator has degree at most gs by the stretch count. For two such frames A and B apply the mixed derivative in all s variables to F_B - F_A. It kills missing-color terms while retaining the entire full-support difference, with each surviving exponent reduced by one.

Each differentiation of P/D^j replaces its numerator by D partial(P) - jP partial(D) over D^(j+1); the numerator degree increases by at most deg(D)-1. Let L be the polynomial least common multiple of D_A,D_B, with constant 1, and let d_min be the smaller denominator degree. The derivative difference over L^(g+s) therefore has numerator degree at most
```
(g+s) deg(L) + g(s - d_min) - s.
```
If its numerator is not zero, the denominator cannot remove its lowest-degree nonzero terms. Returning from the differentiated series to the original adds s to total degree. The resulting witness bound is (g+s)deg(L)+g(s-d_min), at most **2s(g+s)**: deg(L) is at most d_A+d_B, with each d at most s. One can lower the bound further by using the actual numerator degrees. If the pole obstruction passes, D_A=D_B, so the bound simplifies to gs+s deg(D_A), at most **s(g+s)**.

For eight marks and seven ordinary vertices at q=2, the previous sufficient total width was 18,432. It is now at most 272 for arbitrary pairs and at most 136 after the pole obstruction passes. The number of tuples in a brute-force test can still be large.

There is no need to enumerate all those tuples if the exact numerator is in hand. Take any nonzero term c times t^u of **minimum total degree** in the differentiated numerator. Put r_i=u_i+1. Then width r in independent transitive blocks is a smallest-total witness, and the literal q-slot difference there is exactly
```
c * product_i (r_i - 1)!.
```
Indeed the original coefficient is c/product r_i. The transitive widths r allow only componentwise smaller or equal run counts, and every strictly smaller one has lower total degree and vanishes. The remaining run vector is filled in product r_i! ways. This gives an actual witness directly from a sparse coefficient, not a residue-class obstruction.

## 3. Physical regression for the shared-block caveat

Take outgoing masks A=(10,4,9,2), B=(2,12,1,5), obtained by switching vertex 3. Mark 0,1,2, a directed triangle in each. Both have the same marked pole. The converse of A is B after exchanging marks 0 and 2, fixing 1 and 3. Substituting an identical block Q at every mark therefore gives equal slot counts at every q: take converse, swap marks, and use the fact that reversing each individual path preserves the ordered nonempty cover counts of Q. This holds for any shared Q, not just transitive blocks, but it does not change H either.

Independent transitive widths (2,1,1), however, give (H,B)=(9,84) and (11,100). Total replacement width 4 is the smallest possible distinguishing width. The exact symbolic replay recovers a width-4 witness from the derivative numerator at q=1,2,3 and independently checks literal path-count differences. Run `python3 clone_symbolic.py` (SymPy required).

For a height-equal base pair that fails the pole test, take the elementary H=2399/2393 switch at common height 38,344 from [Frame repair](FRAME_REPAIR.md), and mark vertices 0, 3, 8. They form a cyclic component in one frame and a transitive triple in the other. Doubling only mark 0 gives heights 225,088 and 222,280, a smallest-total witness because all widths one give equal heights.

This exposes two separate filters for a useful lattice move: the removed cut must not change the marked strong structure, and the remaining boundary reads must agree at q=2 but differ at q=1. It is the second filter, followed by absolute target placement, that remains constructive work. It can sometimes be pursued with a correlated choice of replacement blocks, where the independent obstruction here is not applicable.
