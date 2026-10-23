# Independent metagraph audit

Date: 2026-10-23. Status: PASS for the n<=7 full-network audit, blue-only n=8, and both synthetic switch controls. No Git mutation.

## 1. Independent metagraph reconstruction

`check_metagraph_independent.py` reconstructs the fixed-path tournaments without importing the research producer. It uses full outgoing-row words and descending degree cells, rather than the producer's upper-triangle integer and ascending degree encoding. Literal reversal-converse symmetry determines blue presentations. Literal off-path complement determines the line pairing. Black multiplicities are observed to be even, not imposed.

The independently reconstructed results are:

| n | vertices | blue components | binary bicycle dimension | critical 2-rank | half-form radical | v2(tree count) | positive Smith valuations |
|---:|---:|---:|---:|---:|---:|---:|---|
| 3 | 2 | 1 | 0 | 0 | 0 | 0 | none |
| 4 | 3 | 2 | 0 | 1 | 0 | 1 | 1 |
| 5 | 10 | 3 | 1 | 3 | 0 | 3 | 1,1,1 |
| 6 | 34 | 23 | 3 | 25 | 1 | 26 | 24 copies of 1; one 2 |
| 7 | 272 | 185 | 11 | 195 | 3 | 198 | 192 copies of 1; three 2s |

The verifier computes the full integer tree determinant by fraction-free Bareiss elimination, not the producer's modular Smith pivot routine. In particular the n=7 determinant has valuation 198. The binary nullity 195 and half-form radical dimension 3 already force total valuation at least 198; equality proves there are no elementary divisors divisible by 8. This is an independent route to the complete 2-primary answer, not merely a repeated rank computation. Both ordinary and optimized Python runs agree after removal of elapsed-time fields.

## 2. The first dyadic theorem

For a blue incidence matrix C over F2, the map ker(CC^T) -> im(C^T) intersect ker(C), given by x -> C^T x, is onto and has kernel ker(C^T). That kernel consists of the potentials constant on the blue components. Therefore nullity Lap(Blue) = c_Blue + dim bicycle(Blue).

For G=Blue+2 PairedBlack, Lap(G) has that same reduction modulo 2. On each connected G-component, use as an integral potential basis the nonroot coordinate vectors followed by the all-ones vector. This is a unimodular basis, and the full Laplacian becomes the reduced Laplacian plus one zero row and column. Consequently the reduced mod-2 nullity is c_Blue-c_G+dim bicycle(Blue), which counts the even Smith factors of the componentwise critical group.

This proof handles loops and parallel lines correctly. Loops have zero incidence columns; parallel lines remain separate columns. Component counts include isolated vertices. The statement does not require all-order connectedness of the tournament graph.

## 3. The half-Laplacian form

For a symmetric integral reduced Laplacian A and V=ker(A mod 2), define beta(x,y)=x^T A y/2 mod 2 using integral lifts. The numerator is even. Changing a lift by 2z changes the quotient by z^T A y, which is even because y is in V. Addition of binary lifts is therefore compatible, and the form is symmetric and bilinear. In the full graph formulation, adding constant potentials contributes exactly zero, so the construction descends to the unrooted quotient by constants.

A vector x in V lifts to ker(A mod 4) precisely when A z=-A x/2 mod 2 is solvable. Symmetry gives im(A mod 2)=V-perp. Thus liftability is exactly membership in rad(beta). Integral Smith coordinates show that this image of the mod-4 kernel has one dimension for each invariant factor divisible by 4. The root, kernel basis and integer lifts do not affect the result.

The color decomposition beta=half-blue-energy+paired-black-energy follows directly from A=A_Blue+2 A_PairedBlack. The form is not determined by its diagonal, and the ranks of the two summands do not determine the rank of the sum.

Additional controls in `check_algebra_controls.py` verify the statement on 964 nonsingular symmetric integer matrices of dimensions 1 through 7, including 804 direct mod-4 lift enumerations and 275 examples with nonzero radical. These are independent integer Smith controls, not graph-only tests.

## 4. Two-color trees and blue-only torsion

For connected G on r vertices, the blue count in a spanning tree ranges exactly from c_B-1 to r-c_L. Both endpoint constructions are valid: maximize a forest inside the components of one color and connect those components using the opposite color. Connectivity of the spanning-tree exchange graph fills every intervening integer count.

At the lower endpoint, every B-component must contribute a black spanning tree and the quotient contributes a blue spanning tree. At the upper endpoint, reverse the colors. Counting the choices of the doubled black lines gives precisely the factors 2^(r-c_B) and 2^(c_L-1) in the manuscript's two endpoint formulas.

`check_tree_controls.py` independently reconstructs the n<=6 polynomial coefficients using exact Bareiss determinants and forward-difference interpolation, and checks the contraction products. All coefficients match. It independently confirms the blue-only 2-primary valuations through n=7. At n=7, a separate symmetric-congruence algorithm uses odd 1x1 pivots or invertible 2x2 blocks; the exact integer determinant independently checks the total valuation. The blue-only positive valuations are 1,1,3,3,4,4,5,5,6,9,10. This genuinely contrasts with the full graph's maximum valuation 2.

The electrical limits are sound: scaling the full conductances by 1/epsilon turns black conductances into 2/epsilon and contracts their components, giving epsilon R_Gepsilon -> R_quotient across components. For vertices inside the same black component, the weak blue conductances disappear in the unscaled limit, giving one half the paired-black resistance. These are electrical limits, not literal complement-orbit dynamics.


## 5. The additional blue-only order-eight hostile

`check_blue8_independent.py` independently enumerates the 4,096 literal reflection-fixed presentations using descending-degree, full-row canonical words. This is the blue-supported graph only, not the order-eight full line graph. It finds 176 supported vertices, one connected component, bicycle dimension 11, and positive blue Smith valuations 1,1,1,1,2,2,3,3,3,3,4. The independently implemented symmetric dyadic elimination agrees with the exact integer tree determinant, whose 2-adic valuation is 24.

The half-blue-energy form has rank 4 and is not alternating. Exhausting all 2,048 binary harmonic classes finds a minimum nonalternating cut with 914 crossing blue lines; 914 is 2 modulo 4. Thus the cut is Eulerian modulo 2 but its half-energy is odd. This disproves an unrestricted alternating or doubly-even strengthening. It does not create an order-eight full-network torsion census.

## 6. Two synthetic incidence controls

`check_switch_controls.py` independently recomputes exact integer Smith forms and tree determinants for both order-five rewires. Repeated edge rows are accumulated as multigraph multiplicities, including the repeated paired-black edge (4,7) in the second control.

- Blue switch: replace (0,5),(6,9) with (0,6),(5,9), preserving every blue degree and the whole black graph. The tree count changes from 41096 to 37184 and the critical 2-rank changes from 3 to 4.
- Black switch: hold the blue graph fixed and replace paired-black edges (1,4),(3,7) with (1,3),(4,7). Every paired-black degree and the full Laplacian modulo 2 are unchanged. The tree count becomes 37952; positive 2-primary valuations become 1,1,4 rather than 1,1,1; and the half-form rank falls from 3 to 2. In the same binary kernel basis, both forms have diagonal (1,1,1): the off-diagonal data are essential.

Both controls pass and are correctly labeled synthetic, not alternative tournament realizations. In particular the second control shows why the blue graph and all color degrees fail to determine higher dyadic arithmetic: the black wiring changes it while leaving the first dyadic operator completely unchanged. The factor of order 16 also blocks a putative exponent-four theorem based only on abstract involutions and those inventories.

## 7. Uniform odd global conductance rescaling

The first two dyadic counts are invariant under arbitrary positive odd global blue and literal-black weights: the mod-2 matrix stays fixed and both coefficients of the half-energy form stay odd. This general statement does not by itself control higher valuations.

For the five actual full graphs n=3,...,7, `check_odd_weights_independent.py` separately verifies the stronger complete-profile claim. Divide the reduced matrix by the odd literal-black weight in Z_2; its residue modulo 8 depends only on the odd blue/black ratio t in {1,3,5,7}. Independent symmetric elimination at exactly three bits completes for all twenty cases with no leftover coordinate, and gives the original histogram each time. Completion at this precision certifies that no factor divisible by 8 is hiding. Hence the full 2-primary groups for these five fixed networks are unchanged for every such odd positive global pair. This is not an all-order statement or an edge-by-edge weight invariance.

## 8. Reproduction and exact file list

The audit code imports only Python's standard library, NumPy and SymPy. No producer code is imported. Some scripts read producer receipt JSON files under the sibling `session36-metagraph/` directory; retain that relative input location when replaying, or provide an equivalent staged directory layout. The scripts can be run from any current directory because inputs are located relative to the script files.

Stage these audit code and receipt files:

1. `METAGRAPH_AUDIT.md`
2. `check_metagraph_independent.py`
3. `metagraph_independent_receipt.json`
4. `check_algebra_controls.py`
5. `algebra_controls_receipt.json`
6. `symmetric_dyadic.py`
7. `check_tree_controls.py`
8. `tree_controls_receipt.json`
9. `check_blue8_independent.py`
10. `blue8_independent_receipt.json`
11. `check_switch_controls.py`
12. `switch_controls_receipt.json`
13. `check_odd_weights_independent.py`
14. `odd_weights_independent_receipt.json`
15. `metagraph_replay_agreement.json`

Optional successful logs are `metagraph_independent.out`, `metagraph_independent_optimized.out`, `algebra_controls.out`, `algebra_controls_optimized.out`, `tree_controls_symmetric.out`, `tree_controls_optimized.out`, `blue8_independent.out`, `blue8_independent_optimized.out`, `switch_controls.out`, `switch_controls_optimized.out`, `odd_weights_independent.out`, and `odd_weights_independent_optimized.out`. Do not stage `tree_controls.out`: it is an abandoned slow full-integer-Smith attempt superseded by the successful independent symmetric-congruence route and exact determinant check.

Ordinary and optimized Python replay results agree after removing only elapsed-time fields. All checks use explicit exceptions, not assertions. The audit is finite-exact where a finite census is claimed; the general statements rest on the proofs above. It does not certify source-search completeness, literature priority, full-network order-eight data, all-order connectedness, or an all-order exponent-four bound.
