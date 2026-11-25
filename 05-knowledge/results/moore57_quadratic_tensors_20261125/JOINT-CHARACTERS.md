# Every fourth-source slot has an integral signed completion

**FINITE-EXACT.** Take the 28-trade partial Moore chart with three saturated branch stars on 56 sheets. One can delete any four of the 54 exterior branches and still realize arbitrary integral marginal counts on the three stars, with prescribed integral sums on each retained branch. The necessary equality is that the total on each star equals the sum prescribed across branches. Multiplicities in this conclusion may be negative. [Three-star construction](THREE-STARS.md) · [Exact certificate](joint_characters.py)

For each of the 54 possible fourth branches, there are 53 signed integral matching matrices, each with unit row and column sums, satisfying all 56 source-sheet equations against each of the three stars simultaneously. This removes all additive obstructions over every abelian group, including mixed obstructions involving all three stars. Nonnegativity is the remaining condition in this integer system for a single new star. Compatibility between successive stars remains a further condition for a Moore graph.

## The integer lattice of demands

Write the stars as \(O,A,B\). For each exterior branch \(h\) and sheet \(x\), let \(v_{h,x}\in\mathbb Z^{168}\) have three unit entries, at the sheets of \(O,A,B\) reached directly from \((h,x)\). Changing a branch assignment from sheet 0 to sheet \(x\) changes the marginals by \(v_{h,x}-v_{h,0}\). These differences belong to the rank-165 lattice \(L\) of triples of integral 56-vectors having sum zero in each vector. Dropping sheet 0 in each coordinate identifies \(L\) with \(\mathbb Z^{165}\).

The certificate partitions the 54 branches into five disjoint packets. For each packet, its difference matrix contains square minors of order 165 with the following absolute determinants.

| Packet | Branches | Absolute determinants |
|---|---:|---|
| 1 | 11 | 3307165519, 21343166339 |
| 2 | 11 | 48920247, 20683932191 |
| 3 | 11 | 964971499, 3201780363 |
| 4 | 11 | 7194767353, 1729777383, 11590632807 |
| 5 | 10 | 265755223, 5333293091 |

Each row of determinants has gcd one. Indeed, for each square minor of determinant \(d\), its integer adjugate shows that the span of the chosen rows contains \(dL\). Bézout's identity therefore shows that each packet spans \(L\) over the integers. Any four deletions leave at least one packet intact.

Start with the prescribed branch sum concentrated at sheet 0 in every retained branch. The desired change to its three marginal vectors belongs to \(L\). Express it as an integral combination of the differences just proved to span \(L\). Each difference preserves its branch sum. This realizes all the requested marginals.

## Fourth-source sheets and joint characters

For a source sheet \((b,x)\), the currently specified paths to each of the three stars reach three different target sheets. To fill its slot, omit branch \(b\), prescribe sum one on each of the remaining 53 exterior branches, and prescribe marginals zero on the three used sheets and one on the other 53 at each star. The preceding construction gives an integral signed solution in each of the \(54\cdot56\) instances.

Now fix the entire source branch \(b\). For its first 55 source sheets, use the construction independently, writing the resulting integers as \(m_{x,h,t}\) for source sheet \(x\), intermediate branch \(h\), and intermediate sheet \(t\). On the last source sheet set
\[
m_{55,h,t}=1-\sum_{x=0}^{54}m_{x,h,t}.
\]
For each \(h\), the matrix \((m_{x,h,t})_{x,t}\) has column sums one by construction. Its last row sums to \(56-55=1\). For a fixed target star \(j\) and sheet \(y\), exactly one \(t\) in each of the 53 branches maps to \(y\). Across all source rows the signed coverage of \(y\) is therefore 53. The three pre-existing routes to \(j\) are permutations, so they use \(y\) three times across the 56 sources. Hence the total demanded exterior coverage is also 53. Subtract the 55 satisfied demands to get exactly the last one. This proves the simultaneous assertion. Requiring every entry to be nonnegative would make each integral matrix a permutation matrix.

Alternatively, let three sheet functions \(f_j\) take values in any abelian group. If \(\sum_j f_j((h,x)\to j)\) is independent of \(x\) on each retained branch, its associated homomorphism annihilates the difference lattice. Since this lattice equals \(L\), each \(f_j\) is constant on its own star. The joint test detects genuine failures in the original cyclic chart. There the coordinate differences \(x_O-x_A\) and \(x_O-x_B\) are each constant along every exterior branch. With their residues modulo 2, these two identities obstruct 2,352 of the 3,024 proposed slots. Modulo 7 they obstruct all 3,024. The 28 trades therefore change the integral feasibility of the routing problem. [Single-board obstruction](LOCAL-LATIN.md)

Run `python3 -B joint_characters.py` from this directory. The checker constructs every difference row from the specified matchings, selects the square minors, evaluates their exact integer determinants, and checks the cyclic counts. Its standard-library determinant routine can optionally be accelerated by python-flint.

## Every sufficiently large multiplicity is achievable at source (1,0)

**FINITE-EXACT certificates, PROVED consequence.** A \(d\)-fold completion of this source slot assigns nonnegative integers to its legal triples. Each of its 53 branches and each of its 53 available sheets on each of the three stars must be used exactly \(d\) times. Every sufficiently large integer \(d\) admits such a completion. This statement concerns the single specified source sheet. The threshold and the status of \(d=1\) do not follow from the argument.

There are 2,505 legal triples after rejecting any choice that would use one of the already occupied sheets. In each branch subtract its first legal triple from the others. Drop one reference sheet from each of the three target coordinates. This leaves a \(2452\times156\) difference matrix. Two maximal minors have absolute determinants 47304122656582883 and 16329197407783105, whose gcd is one. As above, this proves the existence of a signed *unit* cover \(s\) entirely on the legal support.

Give every legal triple baseline weight 354158094286698 and add the 209 positive integer corrections in the certificate. Each prescribed sum is now exactly \(D=35415809428669800\). Dividing by \(D\) gives a rational unit cover assigning **at least 1/100 to every legal triple**. Let \(v\) be the unnormalized positive cover. For \(d=qD+r\), \(0\le r<D\), the integral vector \(qv+rs\) has all prescribed sums equal to \(d\). Since \(v\) is strictly positive on this finite support and there are only \(D\) residues, all entries become nonnegative once \(q\) is sufficiently large.

Run `python3 -B joint_slot_covers.py`. The weighted cover was found by a linear program with a positive baseline. Its integer verification and the two lattice minors can be checked with the Python standard library. [Checker](joint_slot_covers.py) · [Integer weights](joint_slot_weights.txt)
