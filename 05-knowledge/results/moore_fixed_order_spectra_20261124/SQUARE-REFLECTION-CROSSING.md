# Square-sum reflections constrain the 125-label Moore crossing block

This continues the [125-label half model](FIFTH-POWER-HALVES.md) and the repo's [square-sum reflection transport](../../../01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md). The Moore claim concerns only the branch with an equitable equal split preserved by a free order-13 symmetry; the 249-site order-11 code is a separate model. We impose the additional test hypothesis that **only the crossing block** on the two 125-label panels is cyclic after relabelling. Neither internal block is assumed cyclic.

## 1. What the small square-sum graphs transfer

Let \(Q_n\) have labels \(1,\ldots,n\) and edge \(xy\) when \(x+y\) is a square. Each square label \(s\) gives the partial involution \(x\mapsto s-x\), hence a matching on that label's domain. Odd square labels cross parity and even labels preserve it. The exact finite replay gives:

| Window | Edges | Across parity | Within parity | Triangles | Four-cycles | Hamilton paths up to reversal | Hamilton cycles up to reversal |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 25 | 32 | 17 | 15 | 0 | 0 | 10 | 0 |
| 32 | 46 | 24 | 22 | 1 | 0 | 392 | 1 |

Thus 25 labels and 32 edges really coexist at the first window, but they do **not** mean 25 crossing and 32 internal neighbours. The triangle at the second window is \(6,19,30\), using sums 25, 49, 36. It consists of two crossing reflections and one internal reflection, exactly the kind of collision that a Moore wiring must forbid. The small graphs illustrate which information to retain: a reflection label and its domain, not just whether a step crosses a cut.

The actual equal Moore split, if present, has internal degree 32 and crossing degree 25 ([Renteln's correction](https://ajc.maths.uq.edu.au/pdf/79/ajc_v79_p193.pdf)). At any vertex the same-side destinations are exhausted by
\[
1+32+32\cdot31+25\cdot24=1625.
\]
The opposite side is exhausted by \(25+32\cdot25+25\cdot32=1625\). This is exact disjoint routing, not a growth estimate. In the block matrices \(U,V,C\) of the order-13 quotient, the opposite-side equation reads \(UC+CV=13J-C\). The two terms are the two orders of applying a crossing and an internal step.

## 2. A cyclic crossing block is forced to be a 25-transversal

If full reflections indexed by residues \(s\) pair left label \(x\) to right label \(s-x\), negating the right labels makes their sum a circulant matrix. This is the direct square-sum-inspired test. More generally suppose \(C\) is any circulant nonnegative **integer** 125×125 crossing matrix with row sum 25. Its coefficient polynomial is \(f\) with \(f(1)=25\). Here are consequences even if the internal blocks \(U,V\) are arbitrary.

**Bounded spectrum and a lower rank.** On centred vectors the Moore identity gives \(CC^{\mathsf T}=56I-U^2-U\), so the squared singular values are at most \(225/4\). By integrality,
\[
\|C\|_F^2-25^2\ \ge\ 125\cdot25-25^2=2500.
\]
Thus at least 45 centred singular values are nonzero, or rank at least 46 in total. On the other hand the [prior balanced-quotient argument](FIFTH-POWER-HALVES.md) gives ranks at most 117, 105 or 93 for cyclic-orbit counts \(c=30,120,210\), respectively.

**Cyclotomic restriction.** The possible nontrivial zero packets for the integer polynomial \(f\) have sizes 4, 20 and 100, for character orders 5, 25 and 125. The rank lower bound prevents the 100-packet. Nullity at least eight then forces the order-25 packet. This already rules out \(c=210\), whose required nullity is at least 32.

To sharpen this, aggregate the coefficients modulo 25. Vanishing at a primitive 25th root implies that this vector of length 25 is constant over each residue class modulo 5: writing the corresponding polynomial of degree at most 24, it is a degree-at-most-four polynomial times \(\Phi_{25}(z)=1+z^5+\cdots+z^{20}\). Let the five repeated values be \(a_i\). They are nonnegative integers with \(\sum a_i=5\). The order-five coefficients are \(5\sum_i a_i\zeta_5^{ki}\). Parseval and the four squared-singular-value bounds yield
\[
 25(5\sum_i a_i^2-25)\ \le\ 4(225/4).
\]
Anything other than \(a_i=1\) for all five raises \(\sum a_i^2\) to at least 7, which fails the inequality. Therefore **the crossing shifts are binary, with exactly one in each residue class modulo 25**. The order-five packet vanishes as well. The crossing rank is exactly **101**, and its centred kernel consists precisely of functions constant on residue classes modulo 25.

This rules out a literal global lift using square shifts, or a common affine image of square shifts: the squares have just 11 residue classes modulo 25, whereas all 25 are compulsory. This argument is about full modular reflection channels; it does not rule out partial reflections with varying domains, as in \(Q_{25}\) and \(Q_{32}\).

## 3. The 25-channel residue rows force coarser internal structure

From \(UC=-C(I+V)\) on centred vectors and its transpose, \(U,V\) preserve the centred kernels of \(C^{\mathsf T},C\), respectively. Including the constant direction, the internal blocks hence respect the same coarse grouping by residues modulo 25. The full graph would have 25 physical cells of 65 on each side: \(5\) order-13 orbits per cell. For the 25×25 coarse internal quotients \(U_0,V_0\) the crossing block is \(J_{25}\), and
\[
 U_0^2+U_0=56I+40J,\qquad V_0^2+V_0=56I+40J.
\]

Each diagonal entry is the internal degree of one 65-site cell. It is even: subtracting the row sum 32 from the diagonal equation gives \(a+\sum_j b_j(b_j-1)=64\). Girth at least five gives \(a^2+1\le65\). If \(a=8\), the cell attains equality and is a degree-eight Moore graph. Its two nonconstant eigenvalues \((-1\pm\sqrt{29})/2\) must occur equally often by rationality, contradicting trace zero on 65 sites. Therefore the diagonals are at most 6 and the two traces total at most 300.

The complementary centred 100-dimensional spaces for \(U,V\) pair eigenvalues adding to \(-1\) via invertible centred crossing transport there. Thus the coarse trace total equals \(2c+100\). **The \(c=120\) case is impossible**, requiring 340. The only surviving cyclic crossing candidate has \(c=30\). In that case the coarse traces total 160 and each is \(20\bmod30\) (trace \(32-8\cdot24+15m\), even). The unordered trace pairs are exactly the remaining numerical candidates **(20,140), (50,110), (80,80)**. No claim that any pair extends.

A concrete next search is now a 25×25 panel of perfect matchings between 65-site cells, and a same-side complement constraint. For any two distinct same-side cells, routing out and back through the 25 opposite cells gives **25 edge-disjoint permutation matchings** between them: otherwise two distinct common neighbours create a square. The other 40 destinations per source in that cell must be supplied by internal adjacency or a two-step internal route. This is the useful version of the small square-sum reflection picture: use individually labelled matchings, but test their overlaps two steps later. It is not justified to assume either internal graph is circulant or to identify its routes with the six small square labels.

The companion script checks the finite square windows and the integer inequalities. The order-13 assumption remains conditional; see [Mačaj–Širáň](https://www.sciencedirect.com/science/article/pii/S0024379509003735) for the symmetry restrictions.
