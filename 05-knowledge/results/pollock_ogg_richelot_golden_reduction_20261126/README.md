# The supersingular fifteen and Pollock's charts: a marked Richelot bridge

**Proved and replayed:** an incidence-preserving use of the outer six-set, Clebsch line and point counts in every residue characteristic, and a marked Richelot reduction which selects one of the six golden pairs. **Not asserted:** any natural prime-to-chart bijection, an Ogg-prime characterization from the Clebsch cubic, or a proof of Pollock.

This joins the [Pollock–Clebsch cancellation packet](../pollock_five_clebsch_and_defect_20261126/), the [six-branch genus-two curve](../chanpak_six_five_boundary_atlas_20261123/CURVE.md), the [earlier Richelot/Frobenius and outer-six-set work](../session97_overnight_arithmetic_carriers_20261029/supersingular/RICHELOT_FROBENIUS_AND_FANO.md), and the [supersingular level filtration](../supersingular_divisor_shapes_qseries_20261022.md). We use actual labeled kernels and lines, not an equality of two counts.

## 1. Where the fifteen really match

Mark infinity among six branch points and label the other five by the five Pollock coordinate positions. A rational cancellation line has form
\[
 (y_a,y_b,y_c,y_d,y_e)=(u,-u,v,-v,0).
\]
Send it to the branch partition $(ab)(cd)(e\infty)$. This is a bijection of the fifteen lines with the fifteen geometric Richelot kernels, equivariant under the permutations fixing infinity. Its inverse drops the pair through infinity. The map to kernels uses the usual even-subset model of genus-two two-torsion; it is not an identification of two surfaces or of isogeny codomains. Marked infinity is essential to this particular map.

There is a second six-label description, built from the charts themselves. The fifteen charts are edges of the Petersen graph on duads of the five coordinates. The six perfect matchings of Petersen are the crossing edges between an unordered pair of complementary pentagons in $K_5$. Indeed a pentagon and its complement each give a 5-cycle in Petersen, with exactly one crossing edge at each vertex. Distinct perfect matchings intersect in exactly one edge, and each edge lies in two; the small census in the replay also certifies this. Label a chart by the two matchings which contain it.

**The actual line intersection graph is $KG(6,2)$ in these second labels.** In the original five labels, two distinct lines meet precisely if the charts share a duad, or if they omit the same coordinate. In the first case the common point is the pair-difference direction. With the same omitted coordinate the two pairings make a 4-cycle, giving one alternating common direction. Otherwise a nonzero common coordinate would propagate through the alternating pairs and either hit a forced zero or lie in a common pair; there is no intersection. This gives degree $4+2=6$. Two edges in one Petersen perfect matching share no duad, and omit different coordinates: the 5-cycle translation permutes the five crossing edges and translates their omitted coordinate. So chart intersections imply disjoint second labels; disjointness on two-subsets of six also has degree six, proving equality.

In particular the fifteen lines have intersection matrix $A(KG(6,2))-I$. To see its rank, let $B$ be the $6$ by $15$ incidence matrix of the two-subsets. Then $A-I=J-B^{\mathsf T}B$ and $BB^{\mathsf T}=4I+J$. The eigenvalues of $A-I$ are $5$, five copies of $-4$, and nine zeros. Thus the fifteen rational lines span a nondegenerate rank-six subspace of the rank-seven geometric Picard group. This retained intersection form, not the number fifteen, is what lets the graph compute a Frobenius trace.

## 2. The Ogg primes split 7 / 1 / 7, but so do ordinary residue tests

For arithmetic use the **primitive** Clebsch model $S:\ e_1(y)=e_3(y)=0$ with five projective coordinates. On $e_1=0$, $\sum y_i^3=3e_3$. Using the undivided sum of cubes as the equation at 3 would incorrectly produce the whole hyperplane.

This model is smooth except at 5. On $e_1=0$, the partial derivatives of $e_3$ are $e_2+y_i^2$, so a singularity requires all squares equal. In characteristic two all coordinates would be equal and their sum does not vanish unless they are zero. In other characteristics all coordinates are signed copies of a common nonzero value. The only possible vanishing signed sums occur in characteristics three and five. In characteristic three the four-plus-one sign pattern has $e_3=\pm(4-6)a^3\ne0$. At five the only possibility is the constant vector, which is on the surface.

For $p\ne5$ use the Fourier plane spanned by $(\zeta^i)$ and $(\zeta^{-i})$ in cyclic order, with $\zeta$ a primitive fifth root. The orthogonal projector is defined since 5 is invertible: its off-diagonal entries are $(\zeta+\zeta^{-1})/5$ on edges and $(\zeta^2+\zeta^{-2})/5$ on nonedges. These values differ even in characteristics two and three, since their difference squares to 5. Therefore the twelve planes are distinct; none has a zero coordinate functional, unlike the fifteen paired planes. They satisfy the cubic in characteristic three too: expanding $e_3$ gives a cubic in the two Fourier coefficients invariant under the 5-cycle and no degree-three monomial has zero cyclic weight. The same invariant argument works at every $p\ne5$. These account for all 27 geometric lines on the smooth model, by the classical line theorem (see [Huybrechts](https://www.math.uni-bonn.de/people/huybrech/Notes.pdf) for the characteristic-zero realization). Frobenius $\zeta\mapsto\zeta^p$ keeps all twelve lines when $p\equiv\pm1\bmod5$ and exchanges the two in each golden pair for $p\equiv\pm2\bmod5$. In the latter case the action on Picard is nontrivial: distinct lines cannot have identical divisor class, as self-intersection is minus one but intersection with another line is nonnegative. It is an involution fixing the rank-six space from §1, so it acts by minus one on the remaining dimension. In particular
\[
 |S(\mathbb F_p)|=p^2+(6+\chi_5(p))p+1\quad(p\ne5),
\]
where $\chi_5$ takes $+1$ on $\pm1$, $-1$ on $\pm2$. Here we use the standard [surface point-count formula via Picard Frobenius](https://www.cambridge.org/core/journals/mathematical-proceedings-of-the-cambridge-philosophical-society/article/abs/del-pezzo-surfaces-over-finite-fields-and-their-frobenius-traces/3D70794B5F27B0D655F34987CB2EB6F4). In fact over $\mathbb F_{p^r}$ the corresponding trace is $6+\chi_5(p)^r$.

At 5, the constant vector is an ordinary node. Write a point as $y=t\mathbf1+v$ with $v_0=0,\ \sum v_i=0$. The cubic is $3t\sum v_i^2+\sum v_i^3=0$ before division by the unit three. Its conic $\sum v_i^2=0$ in the direction plane is nonsingular: the bilinear form on the remaining four coordinates restricted to sum zero is nondegenerate as 4 is invertible. There are six common conic/cubic directions given by bijective labelings $v$ of the five coordinates by $\mathbb F_5$, modulo affine changes of labels. Each direction is transverse, since on the four nonzero coordinates no quadratic with leading coefficient one can agree with an affine function at four distinct inputs. Thus these exhaust the intersection by Bézout. They give six lines through the node. The other lines project to lines of the direction plane whose intersection with the conic must be a length-two subscheme of the cubic intersection: on a nonnode line $t$ is a linear function of the direction coordinates. The six points are simple, so the projecting lines must be secants between two distinct ones, with at most one lift each. The fifteen old rational lines are distinct and miss the node, hence attain the upper bound. All 21 are rational. Projection counts the points by adding the node to one point per direction off the conic and five points for each of the six special directions. This gives $|S(\mathbb F_5)|=5^2+6\cdot5+1=56$.

Moreover each complementary pair of golden lines reduces to **one** of the six nodal lines. The cycle adjacency eigenvalues $(\pm\sqrt5-1)/2$ both become 2. Its eigenspace in characteristic five is exactly $\ker(T-I)^2$ for cyclic shift $T$, because $T+T^{-1}-2=T^{-1}(T-I)^2$. This is the plane of affine functions on the five cyclic positions, spanned by $\mathbf1$ and the cyclic labeling $v$. Properness of the Grassmannian gives specialization of each plane into that two-dimensional kernel, so equality holds. Changing the cyclic step by a unit only changes the affine label; this is the same complementary-pair equivalence.

The [Ogg predicate](https://arxiv.org/html/1809.07421v2) requires **all supersingular elliptic $j$-invariants in the given characteristic** to be in the prime field. It does not refer to any of these planes. The exact comparison is:

| Ogg subset | Primes | Clebsch lines over the prime field |
|---|---|---:|
| Golden inert | 2, 3, 7, 13, 17, 23, 47 | 15 |
| Golden ramified | 5 | 21 |
| Golden split | 11, 19, 29, 31, 41, 59, 71 | 27 |

The same inert and split rows occur outside Ogg (37 and 61, respectively). Thus neither the field of definition of the fifteen rational lines nor that of all 27 detects Ogg membership. Also the repo's [square-sum graph on Ogg primes](../outdeg_converse_descent_fifteen_20261116/FIFTEEN.md) is a tree; the natural graph on the fifteen lines here is 6-regular, so that earlier prime graph cannot give an incidence-preserving identification either. The 7/1/7 split is a comparison of two predicates, not a characterization of the finite Ogg set.

## 3. At the central prime, superspecial reduction selects a particular golden pair

Use the actual inherited six-branch curve, with the affine scaling $z=3x$:
\[
 C:\quad Y^2\propto z(z^2-1)(z^2-9),\qquad
 z\in\{\infty,-3,-1,0,1,3\}\ {\rm at\ the\ branch\ points}.
\]
The constant twist does not change branch points or the decomposable-kernel predicate. For the actual model in the inherited $x$ coordinate, the reduction at 5 is $y^2=x-x^5$ and is superspecial (zero Cartier–Manin matrix), as proved there. At every good prime $p\ge5$ the six branches remain distinct and rational, so **all fifteen marked kernels are rational**, including at 37 and 61, even though their Clebsch line counts differ.

For a chart $(ab)(cd)(e\infty)$ write $a,b,c,d,e$ for its actual five branch values. The [Richelot coefficient determinant](https://arxiv.org/abs/2101.00917) of $(Z-a)(Z-b),(Z-c)(Z-d),(Z-e)$ is, up to sign,
\[
 \Delta=(a-e)(b-e)-(c-e)(d-e).
\]
This follows by evaluating the first two quadratics at $e$ in the coefficient determinant; the linear factor represents the binary quadratic with one root at infinity. Vanishing is the [decomposed-product criterion](https://arxiv.org/abs/2003.00633), not vanishing of the isogeny. For the scaled marked values its absolute multiset is $0^4,8^3,10^6,18^2$. This produces four decomposed kernels over characteristic zero and every good prime except 5, but ten in characteristic five. In particular the extra decompositions are triggered by **5 alone** in the good-prime range, not by being one of the fifteen Ogg primes.

The five nondecomposed kernels at 5 have a useful exact description. Relative to the unused finite value $e$, the four remaining residues are $\pm1,\pm2$. Pair products agree for the two mixed sign matchings, but disagree for $(e-1,e+1)(e-2,e+2)$, with products $-1$ and $-4$. As $e$ runs through $\mathbb F_5$, this is exactly the Petersen perfect matching crossing between the additive pentagon (differences $\pm1$) and its complementary pentagon (differences $\pm2$). Thus on the marked curve, the surviving five nondecomposed kernels pick out **one golden conjugate pair of Pollock charts**, equivalently its specialized line through the Clebsch node. The cyclic action is the translation from the earlier Paley–Frobenius result. Scaling all finite labels by 3 only interchanges the two pentagons and leaves this unordered pair unchanged.

There is a useful distinction in what survives: at the same prime the two real golden planes merge on the cubic, all fifteen Richelot kernels stay prime-field rational, and six of their *codomain types* change. That is the bridge through actual reduction data. It gives a distinguished choice of five marked charts, but no guarantee that their quadratic corrections meet a prescribed Pollock target. The prior cubic-defect density barrier still applies.

Replay: run `python replay.py`. It checks all 105 chart-pair incidences, the determinant and the selected pentad, and directly counts the surface points in all fifteen Ogg characteristics and four controls. Line counts in its output are the formulas proved above, not a brute-force enumeration.
