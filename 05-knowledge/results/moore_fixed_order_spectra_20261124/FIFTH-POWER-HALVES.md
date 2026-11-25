# Why 125 is a real Moore index, but not a homogeneous lift

This answers the 5 → 25 → 125 analogy without identifying a tournament isomorphism class with a graph. The degree-57 Moore graph is hypothetical throughout. The prior [249-site marked encoding](PENTAGON-PORT-ATLAS.md) assumes an order-11 action; no 125-site replacement for that encoder is claimed here.

## 1. The fifth powers occur in the links across the two halves

Petersen is a pentagon and a pentagram joined by a perfect matching. Hoffman–Singleton, in the [Robertson–Hafner model](https://www.kurims.kyoto-u.ac.jp/EMIS/journals/JACO/Volume18_1/np3862514r41w675.fulltext.pdf), is two 25-point halves, five pentagons versus five pentagrams; every vertex has one neighbour in each opposite five-set. Put side size and crossing degree in separate columns:

| Graph or conditional split | Side sizes | Internal degree per half | Crossing degree |
|---|---:|---:|---:|
| Petersen | 5 + 5 | 2 | 1 |
| Hoffman–Singleton | 25 + 25 | 2 | 5 |
| Degree 57, **if** it has an equitable equal split | 1625 + 1625 | 32 | 25 |

In a regular graph an equitable split into equal halves has internal degree a on both sides, by double-counting the crossing edges. Its sign vector is a nonconstant eigenvector of eigenvalue \(2a-d\). The Moore identity \(A^2+A=(d-1)I+J\) gives roots \(1,-2\) at degree 3, \(2,-3\) at degree 7, and \(7,-8\) at degree 57. Since \(2a-d\) is odd in each case, the rows follow. The 32/25 conclusion also appears in [Renteln’s corrigendum](https://ajc.maths.uq.edu.au/pdf/79/ajc_v79_p193.pdf).

Thus **1 → 5 → 25 is exact for the crossing degrees**, while the side sizes are **5 → 25 → 1625**, not 125. The missing factor is \(1625/125=13\). In particular a direct copy of the small graphs' internal pentagon rule (internal degree 2) cannot supply an equitable degree-57 half.

There is a literal way to recover 125 on a side. If the missing graph admits an order-13 action, it is semiregular (all orbits have size 13, by [Mačaj–Širáň](https://www.sciencedirect.com/science/article/pii/S0024379509003735)). If it **also** has an equitable equal split preserved by that action, the two sides have **125 orbit labels each**. An order-13 action cannot swap two sides of a split it preserves setwise. This suggests a 125-vertex tournament indexed by one half's orbits, or a paired construction on that index set, with the two internal blocks and crossing records separately carried. Neither existence of the action nor existence of a compatible split follows from the arithmetic.

The cross block cannot be a full-rank affine incidence matrix in this branch. Here is a concrete filter. In the order-13 quotient write \(U,V\) for the two symmetric internal blocks (row sums 32) and \(C\) for the crossing block (row and column sums 25). On centred vectors the off-diagonal Moore equation is \(UC=-C(I+V)\). If the centred crossing rank is h, the image of C and the orthogonal complement of its kernel carry paired eigenvalues of U and V summing to −1. The remaining \(2(124-h)\) centred dimensions have eigenvalues 7 or −8 by the diagonal Moore equation. Consequently
$$
                 \operatorname{tr}(U)+\operatorname{tr}(V)
                 \le 64-h+14(124-h)=1800-15h.
$$
The 250 order-13 orbits have internal degree 0 or 2: two distinct unoriented steps would form a square. In a nontrivial Fourier character the dimension is 250 and the eigenvalues are 7 and −8, so the trace is an integer. Rationality of the trace forces equal counts t of each of the six internal step classes (the only rational relation on primitive powers is the cyclotomic relation). Its value is −t, so \(t\equiv5\pmod{15}\) and \(6t\le250\). This rederives cyclic census \(c=30,120,210\). Applying the trace inequality to \(2c\), and adding the constant mode, gives

| Cyclic orbit count c | Max rank of the 125×125 crossing block | Min nullity |
|---:|---:|---:|
| 30 | 117 | 8 |
| 120 | 105 | 20 |
| 210 | 93 | 32 |

For the smaller actual graphs, the corresponding crossing ranks are **5 of 5** in Petersen (a permutation matrix) and **21 of 25** in Hoffman–Singleton. In the latter, rows are indexed by affine points: a point has five lines, two points in the same vertical column have no common line, and different columns have one. Thus the squared crossing matrix is $(5I+J-\mathrm{diag}(J_5,\ldots,J_5))$, with a four-dimensional kernel. The new order-13 rank restrictions continue the loss of crossing rank and do not just extrapolate sizes.

For comparison, on a full equitable split (no order-13 assumption) the same trace argument gives rank at most 1521 of 1625. The centred crossing block has squared Frobenius norm \(1625\cdot25-25^2=40000\), and each squared singular value is at most \(225/4\) by \(CC^{\mathsf T}=56I-U^2-U\) on that space. Thus the centred rank is at least 712. These rank restrictions concern a graph incidence block, not an ordinary tournament orientation.

## 2. A different interpretation of 125 is provably not homogeneous

Suppose instead we try to partition all 3250 vertices into **125 equitable cells of size 26**, and then put one tournament vertex per cell. Write B for the symmetric weighted quotient. It satisfies \(B^2+B=56I+26J\) and each row sums to 57. If \(a_i=B_{ii}\) is a cell's induced regular degree, the diagonal equation subtracting the row sum gives
$$
             a_i+\sum_j B_{ij}(B_{ij}-1)=25.
$$
Thus \(a_i\) is odd. Girth five bounds \(a_i^2+1\le26\), so \(a_i\le5\). Equality would make the cell a degree-five Moore graph: every vertex reaches all 26 vertices in at most two steps. Such a graph is impossible even by trace: its nonconstant roots \((-1\pm\sqrt{17})/2\) require equal multiplicities by conjugacy, but there are 25 nonconstant dimensions. Hence \(a_i=1\) or 3.

If s cells have degree 3, the quotient trace is \(125+2s=57-8\cdot124+15m\), where m is the multiplicity of 7. Thus \(s\equiv10\bmod15\), giving \(s=10,25,40,55,70,85,100,115\). **Every such quotient must have both kinds of cells and has some multiple-weight links in every row.** No transitive cell action (in particular no translation-invariant weighted quotient on the group \(\mathbb F_5^3\)) can furnish this interpretation. This does not exclude asymmetric 125-cell quotients or arbitrary 125-site encodings.

Conversely, if “125” is the **size** of each of 26 equitable cells, constant internal degree would have to be 2: it is even by the same diagonal equation, at most 11 by girth, and \(26a\equiv57-8\cdot25\pmod{15}\). A cyclic-transitive quotient on those 26 cells still fails. Its integer circulant's nontrivial characters of orders 2, 13 and 26 have multiplicities 1, 12 and 12. Trace 52 requires 13 copies of eigenvalue 7, hence the order-2 character and exactly one of the size-12 classes take value 7. At any nonzero even shift modulo 26, each size-12 character sum is −1 and the order-2 value of the character is +1. Fourier inversion demands quotient weight \((57+7-7+8)/26=5/2\), not an integer. This obstruction again concerns the proposed uniform structure, not all 26-cell partitions.

## 3. Which tournaments are natural analogues?

On five cyclic labels orient differences +1 and +2 positively. The +1 and +2 layers are directed pentagon and pentagram on the same five labels. This specific regular cyclic tournament has **15** Hamilton paths, and is natural for the Petersen matching that pairs just those five labels. It is a marked representation of the matched *index*, not a claim that its Hamilton count reconstructs Petersen.

One honest analogue on 25 and 125 labels is the self-similar lexicographic sequence \(T_5,\ T_5[T_5],\ T_5[T_5[T_5]]\), on \(\mathbb F_5,\mathbb F_5^2,\mathbb F_5^3\). Orient the first coordinate on which a pair differs using the cyclic 5-tournament. This gives a concrete isomorphism class at each size; translation is regular, so the path counts are divisible by 5, 25, 125 respectively (the action on ordered Hamilton paths is free). No Paley orientation is being assumed. The 25-label tournament has **2,398,512,090,439,546,875** Hamilton paths, checked by compressing a Hamilton word to its consecutive block runs. The run word follows the outer 5-tournament and the \(r\) runs assigned to one block can be filled in 15, 90, 210, 240, 120 ways for \(r=1,\ldots,5\). This tower is a qualitative analogue, not an edge encoder: in particular the transitive-on-labels pattern cannot describe an equal 125-cell weighted quotient found above.

**Search recommendation.** Aim at the order-13 plus equitable-half branch: two 125-label panels, a 25-regular crossing block with large nullity, and 32-regular halves. Use the 125-label tournament as a marked address system, not as a reason to assume a uniform graph lift. The two competing ways to read 125 (side-orbit count vs. a uniform whole-graph quotient) now have different, checkable predictions.
