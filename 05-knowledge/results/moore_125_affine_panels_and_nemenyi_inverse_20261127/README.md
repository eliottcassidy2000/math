# Two 125-label panels: the 225+25 split and a Neményi-type inverse fiber

This extends the [fifth-power halves](../moore_fixed_order_spectra_20261124/FIFTH-POWER-HALVES.md), the [450 general-position bound](../general_position_moore_design_and_pollock_separator_20261127/README.md), and [Session 324's Neményi comparison](../session324_busch_inverse_and_gap_structure_20261113/BUSCH_NEMENYI_FERMAT.md). The degree-57 graph is still only a hypothesis. The new obstruction is for one specific crossing template, not every two-panel quotient.

## 1. A literal 1,5,25 progression and what 225+25 can mean

In the Petersen and Hoffman–Singleton models, the crossing degree of the two equal halves is 1 and 5. For degree 57, an equitable equal split would have crossing degree 25: the sign vector of the halves is an odd nonconstant eigenmode, hence has eigenvalue 7, giving internal degree 32 and crossing degree 25. A semiregular order-13 action preserving that split gives \(1625/13=125\) orbit labels per half, as in the prior work (see [Mačaj–Širáň](https://www.sciencedirect.com/science/article/pii/S0024379509003735) for the order-13 restrictions).

There is a concrete common *crossing template* behind the numbers. In dimension \(r=1,2,3\) over \(\mathbb F_5\), use two panels of pairs \((x,z)\) and \((a,c)\), with \(x,a\in\mathbb F_5^{r-1}\), \(z,c\in\mathbb F_5\), and match by \(z=a\cdot x+c\). Panel widths are \(5,25,125\); crossing degrees are \(1,5,25\). The first two are the actual crossing blocks. The third is a candidate for the **orbit quotient**, not an edge list of the full graph.

For \(r=3\), write \(C\) for this \(125\)-by-\(125\) matrix, and let \(D\) be the block diagonal of 25 copies of \(J_5\), the 25 five-point columns with \(x\) fixed. Distinct \(x,x'\) give one nontrivial equation on the two slope coordinates, hence five shared neighbors; distinct points on the same column have none. So
\[
CC^\top=25I+5(J-D).
\]
A panel has 100 independent within-column differences and 25 column means. Across both panels the coordinate identity is literally
\[
2\cdot125=200+25+25=225+25.
\]
These are coordinates, not a canonical split of 250 vertices into a distinguished set of 225 and 25, and the 225 in this identity is not the list of matched edges in the general-position problem.

There is nevertheless a literal intersection with that extremal problem. If a size-450 general-position set \(S\) and an equitable equal split both exist, **exactly 225 selected vertices lie in each half**. The centered \({\bf1}_S\) lies in the \(-8\) eigenspace, orthogonal to the sign vector of eigenvalue 7. The cut between the halves must put an odd number of the 225 selected matching edges across, since the internally matched vertices in each half occur in pairs. An order-13 action cannot fix \(S\) setwise, since 450 is not divisible by 13. Thus even in a two-panel quotient the selected vertices require per-orbit occupancies, not a choice of whole orbit labels. If \(y\) records the \(0,\ldots,13\) selected vertices per orbit and \(B\) is the quotient adjacency, then summing the equality \((A+8I){\bf1}_S=9{\bf1}\) over orbits gives \((B+8I)y=117{\bf1}\), with each panel's occupancies summing to 225. No equality case is assumed to exist.

## 2. The affine crossing template eliminates two of the three order-13 strata

Assume **all three** extra hypotheses: order-13 symmetry, a compatible equitable equal split, and the displayed affine \(0,1\) crossing block \(C\) in the orbit quotient. Let \(U,V\) be the internal adjacency blocks, of row sum 32. The prior note's three possible totals of internally cyclic orbits are 30, 120 and 210. The rank 101 of \(C\) already rules out 210 by its rank table; the column geometry below also rules out 120.

The off-diagonal Moore identity implies \(U\) commutes with \(CC^\top\) (on centered vectors \(UC=-C(I+V)\); extend to constants by regularity). Thus \(U\) preserves both the within-column zero-sum subspace \(W\) of dimension 100 and the column-constant subspace of dimension 25. On \(W\), \(CC^\top=25I\), so \(U^2+U=31I\). Its two roots \((-1\pm\sqrt{125})/2\) have equal multiplicity by rationality; this adds exactly \(-50\) to \(\operatorname{tr} U\).

Let \(E\) be the internal **25-column** quotient on the other subspace. It is symmetric, integral, of row sum 32, and
\[
E^2+E=56I+40J_{25}.
\]
Indeed the 125-label diagonal identity has constant term \(13J_{125}\) because each orbit has size 13. On column constants this acts as \(65J_{25}\), and the crossing Gram subtracts \(25J_{25}\). If \(e_i=E_{ii}\), subtracting the row sum from the diagonal identity yields \(e_i+\sum_j E_{ij}(E_{ij}-1)=64\); hence \(e_i\) is even.

These 25 columns per half are equitable cells of **65 actual vertices** each. The induced regular graph of degree \(e_i\) has girth at least five, so counting its radius-two ball gives \(e_i^2+1\le65\). Equality at \(e_i=8\) would give a degree-eight diameter-two Moore graph on 65 vertices. Its two nonconstant eigenvalues would be \((-1\pm\sqrt{29})/2\), each of multiplicity 32 by rationality, giving trace \(8-32\ne0\). This is impossible. Consequently \(e_i\le6\), and \(\operatorname{tr}E\le150\).

If \(a\) is the multiplicity of eigenvalue 7 on the 24 centered column coordinates, then
\[
\operatorname{tr}E=15a-160,\qquad
2c_U=\operatorname{tr}U=15a-210,
\]
where \(c_U\) is the number of internally cyclic order-13 orbits on that half. The second relation uses the fact that each such orbit contributes diagonal weight 2. It forces \(a\) even, and the trace bound gives \(a\le20\), so \(c_U\le45\). The other half obeys the same bound. Thus the affine template leaves **only total 30**, with possible splits \((0,30),(15,15),(30,0)\). These are survivors of the tests, not realizations. Arbitrary nonaffine crossing blocks are not restricted by this argument.

## 3. Neményi: count the true inverse fiber before choosing marks

Over \(\mathbb Q\), \(C^\top\) recovers all 100 differences on the point panel and one total mean. Its kernel consists exactly of the 25 column-constant values with their sum zero, a 24-dimensional inverse fiber. On \(W\) the decoder is \(C/25\) applied to the centered crossing data. To recover a point signal fully, one marked sample in each of 24 columns is sufficient and 24 scalar linear marks are necessary for arbitrary data; recover the last mean from the total. The same statement holds for the opposite slope panel. The proposed 25 mean records per half thus contain precisely 24 hidden coordinates and one already visible mean.

The [Sáenz formulation of Neményi's problem](https://link.springer.com/article/10.1007/s00707-008-0118-3) is a different inverse problem about stresses and their trajectory nets. Session 324 gives an exact local Cartesian five-constant model. To make the analogy literal at the level of sampling, embed the signed quadratic stress data
\[
a y^2+b y+c,\qquad -a x^2+d x+e
\]
as column-constant signals on the two panels, interpreting the coordinate labels as rational representatives. Centered crossing readouts erase all five parameters, since crossing sends column means to a constant. Three evaluations in the first component at \(y=0,1,2\), and two in the second at \(x=0,1\), recover all five (sampling determinant 2). If **uncentered** crossing readouts are retained, the two total means supply two independent conditions, leaving three hidden parameters instead. This five-dimensional model is a marked subspace of a much larger inverse fiber; it is not a reduction of continuum elasticity or a claim that arbitrary Moore contacts have five degrees of freedom.

Finally the arithmetic field matters: the affine \(125\)-label incidence has rank **101 over \(\mathbb Q\)** but only **35 mod 5**. Over \(\mathbb F_5\), each row is the function \(1-(z-ax-by-c)^4\) of \(a,b,c\). The degree-at-most-four monomials give the 35-dimensional upper space, and their coefficient functions in \(x,y,z\) are scalar multiples of the distinct monomials \(x^i y^j z^{4-i-j-k}\) for output monomials \(a^i b^j c^k\), except that the constant coefficient uses \(1-z^4\). They remain independent on \(\mathbb F_5^3\) since their individual degrees are below 5. The exact Smith replay gives 35 unit factors, 56 factors of 5, 10 factors of 25, and 24 zeros. Counting a modular nullity as if it were a rational family of contacts would dramatically overstate the freedom.

Run `python3 replay.py` with SymPy for Gram, Smith, side-trace, and sampling checks.
