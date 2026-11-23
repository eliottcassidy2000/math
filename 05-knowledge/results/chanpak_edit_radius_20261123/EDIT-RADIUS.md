# E61 two-edit screen and the height cost of bounded arc edits

This extends the [same-strong-core seam restriction](../chanpak_kaprekar_decimal_20261123/FIRST-JUMP-SQUARES.md). Write \(H(T)\) for the ordinary directed Hamiltonian-path count. Edits here are literal reversals of individual arcs; an edit budget counts distinct arcs after cancelling repeated reversals.

## Neither exact E61 endpoint lies at transitive edit distance at most two

Let \(E=2^{60}(2^{61}-1)\). A transitive tournament has height 1. One flip gives height 1 if the two vertices are consecutive, and otherwise \(1+2^s\), where \(s\) is the number of interior vertices. Indeed, the path using the new backward arc cuts the interior into two independent increasing pieces. Neither \(E\) nor \(E-2\) is a power of two: the first has the odd factor \(2^{61}-1\), and the second has valuation one and exceeds two.

For two distinct flipped pairs, mark their endpoints in the transitive order, numbered from zero, and let \(x,y,z\) be the numbers of unmarked vertices in the consecutive gaps (omit \(z\) in the three-mark case). Ignore exterior singleton padding. Put \(A=2^x,B=2^y,C=2^z,T=3^y\), and let \(\delta_w=[w=0]\).

| Flipped pairs | Exact final height |
|---|---|
| \(01,02\) | \(1-\delta_x-\delta_y+A+2AB-B\) |
| \(01,12\) | \(1-\delta_x-\delta_y+\delta_x\delta_y+A+B+AB\) |
| \(02,12\) | \(1-\delta_x-\delta_y+2AB-A+B\) |
| \(01,23\) | \((1-\delta_x+A)(1-\delta_z+C)\) |
| \(02,13\) | \(1+2AB-A\delta_z+2BC-\delta_x C+2ACT\) |
| \(03,12\) | \(1-\delta_y+4ABC-2AC+B+2ACT\) |

**Symbolic proof.** Prescribe the ordered list of backward arcs actually used, of length zero, one, or two. Cut the path into increasing runs at these arcs. Assign each mark to exactly one run. Each used arc fixes the larger mark as the final vertex of one run and the smaller mark as the first of the next. Test these endpoint constraints; they also bound which runs are available to each unmarked gap. A flipped pair assigned to the same run cannot occur as a forward adjacency. If no mark of that run separates the pair, impose the requirement that some unmarked vertex in the enclosed gaps be assigned to the run. Apply inclusion–exclusion on these at most two requirements. In every summand the vertices of each gap choose independently among the surviving runs, and their order within a run is forced. This gives monomials \(\prod_j c_j^{w_j}\), with \(c_j\in\{0,1,2,3\}\), \(0^0=1\), and \(0^w=[w=0]\). The included symbolic compiler compares the resulting integer coefficients with the six displayed formulas, not just their numerical specializations.

Any positive-width gap spanned by a flipped pair lies in one nontrivial final SCC: each interior unmarked vertex makes a directed triangle with the pair. The middle gap for disjoint pairs \(01,23\) is irrelevant to the formula and can be set to zero. On an SCC of order \(n\), [Busch's bound](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf) gives \(H\ge5^{(n-1)/3}\). An SCC height is a factor of the final height, and \(E+1<2^{121}\); hence each relevant width is at most 157. Testing the larger box \(0\le w_j\le157\) is safe.

| Flipped pairs | Width tuples | Survivors mod 65521 for \((E-1,E+1)\) | Survivors mod both 65521 and 65519 |
|---|---:|---:|---:|
| \(01,02\) | 24964 | \((1,1)\) | \((0,0)\) |
| \(01,12\) | 24964 | \((0,0)\) | \((0,0)\) |
| \(02,12\) | 24964 | \((1,1)\) | \((0,0)\) |
| \(01,23\) | 24964 | \((0,0)\) | \((0,0)\) |
| \(02,13\) | 3944312 | \((50,46)\) | \((0,0)\) |
| \(03,12\) | 3944312 | \((0,0)\) | \((0,0)\) |

These are exact filters on an exhaustive finite box, not a density argument. Run [the replay](e61-two-edit.js) with `node e61-two-edit.js`; it checks the symbolic coefficients, independently checks 108 small formula instances by subset DP, and tests the two moduli. This excludes both target heights at transitive edit distance at most two at **every** order. It does not exclude farther strong cores or give either side of the seam.

## Fixed starting height: finite templates and a polylogarithmic exact oracle

Fix an attained starting height \(h\) and let \(t\) be the maximum possible number of nontrivial SCCs for a tournament of height \(h\). The SCC product law and Busch's bound imply that the ordered list of nontrivial starting cores is chosen from a finite set: \(t\le\log_3 h\) and their total order is at most \(\lfloor3\log_3 h\rfloor\). The latter follows by multiplying \(H(C)\ge3^{|C|/3}\) across nontrivial SCCs, which holds for triangles and follows from Busch for orders at least four. Use the same bound with \(X\) for the total order of the nontrivial final SCCs when the final height is at most \(X\).

Mark the at most \(2r\) originally singleton vertices incident to one of the at most \(r\) flips, and remember all starting cores. Unmarked singleton blocks between these units have unchanged, uniform arcs toward all vertices outside the block, and are transitive internally. If any member lies in a nontrivial terminal SCC, a directed return route outside the block connects a vertex beaten by the block to one beating the block. To see this, follow a directed cycle containing that member. Take the first external vertex after exiting the block and the first block return; directed paths internal to a transitive block alone cannot make a cycle. Uniformity means the same outside route works for every member of the block, so they all lie in that SCC. Otherwise they are all terminal singletons and the entire block can be deleted. Exterior unmarked blocks are always deletable: they still precede or follow everything. Deletion of singleton SCCs preserves height on both ends.

There are therefore at most \(d=2r+t-1\) width slots for \(r\ge1\), and their **sum** is at most \(L=\lfloor3\log_3 X\rfloor\). For fixed \(h,r\), only finitely many choices of starting cores, marks relative to the cores, and reversed pairs remain. Consequently the number of attainable final heights in \([1,X]\) from this fixed starting height within \(r\) edits is
\[
O_{h,r}\left(\binom{L+d}{d}\right)=O_{h,r}\big((1+\log X)^d\big).
\]
For \(r=1\) the [two-slot formula](../chanpak_kaprekar_decimal_20261123/FIRST-JUMP-SQUARES.md#cross-component-gain-and-singleton-budget) sharpens this to \(O_h(1+\log X)\): once the rooted endpoints and intervening nontrivial components are fixed, every intervening singleton multiplies the gain by two. The possibilities are finitely many exceptions and finitely many dyadic rays \(h+a\,2^s\). In particular this is a statement about exact heights, not their residues.

**The exponent is sharp for \(h=1,r=2\).** In the crossing pattern \(02,13\), take \(1\le x<y<z\le m\) and let \(D=(H-1)/2=2^{x+y}+2^{y+z}+2^{x+z}3^y\). Then \(S=v_2(D)=x+y\) and \(T_0=v_2(D/2^S-1)=z-y\). The integer \((D/2^S-1)/2^{T_0}=3^y+2^{2y-S}\) is strictly increasing in \(y\) for fixed \(S\), so the height decodes \(y\) and hence \(x,z\). These give \(\binom m3\) distinct heights at most \(7\cdot12^m\), proving \(\Theta((\log X)^3)\) with the upper bound above.

**Exact target oracle.** The no-edit case is immediate. For a given core-and-edit template, the unmarked vertices in each block are exchangeable as ordered sets: only their increasing order, not the actual positions, matters in an induced subtournament. Use a dynamic program tracking the set of remaining fixed vertices, the remaining size of each block, and the last symbol (a fixed vertex or block, plus an initial state). Append an allowed fixed vertex or an allowed nonempty run of \(q\) vertices from a block, with weight \(\binom{w}{q}\) when the block has \(w\) vertices left. Do not append two consecutive runs from the same block; they would be a single maximal run. Edges across distinct symbols are uniform when a block is involved, so compatibility depends only on the symbols. For a fixed selection of the \(q\) vertices, the run is forced to be increasing; different selections leave equivalent ordered sets. This recurrence counts each Hamilton path once by its maximal runs. Tables for all total widths up to \(L\) have \(O_{h,r}(\binom{L+d}{d})\) states, polynomial-size transition lists, and polynomial-bit counts (bounded by a factorial at logarithmic order). Comparing each full-state count with the exact target \(X\) gives a complete oracle polynomial in the bit length of \(X\) for fixed \(h,r\), with potentially large parameter-dependent constants. The oracle is not an existence proof for all \(X\).

## A growing edit budget for broad coverage

For the varying-\(r\) lower bound, a simpler count avoids enumerating templates. Let \(\mathcal H\) be a fixed finite bank of starting heights, and let \(\delta>0\). Any witness with final height at most \(X\) can delete **every** vertex singleton in both endpoints, even if incident to a flipped arc. Both endpoint heights are preserved, and the induced tournaments differ by no more flips than before. The remaining order is at most \(n=\lfloor3\log_3 X\rfloor+O_{\mathcal H}(1)\), the bound on the union of the two nontrivial SCC supports. Pad at one end to a common \(n\) if needed. Up to isomorphism the number of source shapes at order \(n\) is \(O_{\mathcal H}(n^t)\), where \(t=\max_{h\in\mathcal H}t(h)\): choose the finite ordered core list, then distribute singleton padding among the gaps. For each shape, label the vertices once; every at-most-\(r\) edit gives an arc subset, with no labelling multiplier. Thus for \(M=\binom n2\) the number of possible final heights is at most
\[
O_{\mathcal H}(n^t)\sum_{j\le r}\binom Mj .
\]
For \(1\le r\le M/2\), the sum is bounded by \((r+1)(eM/r)^r\).

Write \(Y=\log X\), with natural logarithms. Covering \(\delta X\) distinct heights up to \(X\) requires
\[
r \ \ge\ \frac{Y}{\log Y+\log\log Y+O_{\mathcal H,\delta}(1)}
       \ =\ (1-o(1))\frac{\log X}{\log\log X}.
\]
Indeed the logarithm of the bound is at most \(r\log(C_{\mathcal H}Y^2/r)+O_{\mathcal H}(\log Y+\log(r+1))\). Evaluate this at \(r=Y/(\log Y+\log\log Y+C)\) for a sufficiently large constant \(C\); it is below \(\log(\delta X)\) for large \(X\), and the main bound is increasing up to this threshold. If \(r>M/2\), the claimed lower bound already holds. This improves the coordinatewise-width estimate, whose looser exponent yields a leading constant of \(1/3\). It gives a necessary edit cost for a proof strategy from a fixed bank, not an upper bound on what is possible with more edits or with starting heights that grow with \(X\).
