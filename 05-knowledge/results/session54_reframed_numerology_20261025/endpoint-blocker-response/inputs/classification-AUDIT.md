# Independent cycle-potential classification audit

Verdict: PASS. Frozen producer manifest 6b05de37509f3a6f041362b21c8bb98f09ec8e27c4c9470df90bf0c3433a1693.

## Mathematical argument

Arc-transitive orbit averaging applies separately to each Hamiltonian cycle, without any transitivity assumption on cycles. A minimum blocker of size n-1 therefore intersects every Hamiltonian cycle exactly twice. The isolated-vertex star supplies the upper bound, and one Hamiltonian cycle plus arc transitivity supplies the lower bound.

The bipartite support graph has 2n vertices and is connected. Row-column potentials have dimension 2n-1; every cycle difference is orthogonal to them. Independently verified maximal minors give rational cycle-difference rank 8 for P7 and 34 for P11, so their orthogonal complements equal precisely this potential space. The determinant -2 is a rational-rank certificate, not a saturated integral cycle-lattice certificate.

A rational potential representation of an integral edge word becomes integral after fixing one potential to zero and propagating on the connected support tree. Normalization by the minimum gives nested threshold sets. Compatibility follows from nonnegative edge values; each cut has size r times its integer deficiency. Zero deficiency would give no support edges across either direction of the bipartition, contradicting connectedness at a nontrivial active level. Binary edge values force disjoint level cuts, and total size 2r forces either one deficiency-two cut or two deficiency-one cuts. No potential-height case is omitted.

The doubly regular neighborhood-union bound follows by writing h_u=|N^+(u) intersect T|, using sum h_u=rt and sum h_u^2=rt+lambda*t*(t-1), and applying Cauchy-Schwarz. Both polynomial comparisons in the README are correct. Deficiency one has only t=0,n-1. Deficiency two has only t=0,9 at P11, and t=0,1,4,5 at P7. At P7 t=4, S omits one c, so compatibility gives T contained in V minus N^+(c)={c} union N^-(c), and equality of sizes forces that set. This accounts for both exceptional seven-element families.

Two nested deficiency-one cuts cannot be two distinct rows or two distinct columns. The remaining disjoint column/row pair is either a star or a mixed pair with w->v. The verified forward endpoint Hamiltonian path, transported by verified affine automorphisms, excludes every non-star mixed pair. Every Hall cut blocks HPs by predecessor deficiency, and every star isolates a vertex. Thus the analytic argument proves completeness, not just a candidate census.

The general sufficient criterion remains conditional on arc transitivity, a Hamiltonian cycle, connected support, maximal cycle-difference rank, and, for the Hall/star simplification, double regularity and a forward endpoint witness. The deficiency-two expansion gives only the two endpoint Hall families when r>=5; r=3 is the explicitly treated exceptional case. No all-regular or all-Paley maximal-rank assertion is made.

## Independent finite checks

The audit uses Fraction Gaussian elimination rather than the producer's Bareiss implementation. It obtains determinants 1 and -2, checks every cycle and the row/column rank upper bound, checks each supplied support-tree edge, and verifies every affine map preserves all arcs.

An independent direct enumeration of all pairs of vertex subsets S,T gives deficiency-one/two cut classifications. It is not an enumeration of deletion sets. Every resulting binary potential word is independently reconstructed by integral support propagation, normalized, and checked to have positive level deficiencies totaling two. There are 84 candidates at P7 and 176 at P11, respectively partitioned into 56/28 and 110/66 one-level/two-level words.

Independent forward subset HP-DP returns exactly the published 63 and 121 blocker masks. Every surviving path in all 21 and 55 rejected mixed candidates is validated edge by edge; both endpoint witnesses are checked. All candidates have intersection two with all listed cycles, explicitly confirming this condition alone is insufficient for HP blocking.

Normal and optimized audit outputs are identical. Producer replay passes exact-file-set/digest checks, arbitrary-directory normal/optimized replay, and mutation hostiles. The release bundles every producer file and all audit inputs; replay additionally rejects missing files. No producer or repository edits occurred.

Prior P7 counts and P11 Hall counts are explicitly credited as recovered results. Repository search novelty claims are not expanded into literature-priority claims. This packet supersedes earlier P11 total-count uncertainty with the proved total 121.
