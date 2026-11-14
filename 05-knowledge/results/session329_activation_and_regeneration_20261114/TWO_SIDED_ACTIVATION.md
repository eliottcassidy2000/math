# Where six- and seven-vertex reservoirs really enter

Use a cyclic triangle with modules \(U\to W\to x\to U\), with \(x\) a singleton, and nonempty reservoir tournaments \(U,W\). This is a two-module case of the repository's [cyclic substitution path-cover calculus](../../../01-canon/theorems/THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant.md). The output is strong. Write \(H_i=H(T_i)\), \(q_i=Q_{T_i}(2)\), and \(K_i=q_i-H_i\).

## The joint term, and why it matters

A Hamiltonian path's module word is a subword of repeated cyclic order; the unique \(x\) can separate at most two visits to a given module, but **can separate two visits to both**. With both visited only once there are three orders. With exactly one module visited twice the only possibilities are \(UWxU\) or \(WxUW\). With both visited twice the order is \(UWxUW\). A two-visit module contributes \(q_i-2H_i\) because an ordered two-slot cover must use both slots. Thus the exact count is
\[
\begin{split}
N&=3H_UH_W+(q_U-2H_U)H_W+H_U(q_W-2H_W)
           +(q_U-2H_U)(q_W-2H_W)\\
 &=K_U K_W+2H_UH_W.
\end{split}
\]
The missing final term in the preliminary calculation would have treated two individually possible splits as mutually exclusive. It predicts the wrong height by a large amount. The replay counts the full output graphs, including paths using both splits.

Also \(q_T\equiv(-1)^{|T|}H_T\pmod3\) by the inherited slot-polynomial parity. Thus an odd-order module has \(K_T\equiv H_T\pmod3\), and an even-order module has \(K_T\equiv0\pmod3\). Two odd orders always give \(3\mid N\); in the other parity patterns \(N\equiv2H_UH_W\pmod3\). This is a necessary constraint on unit targets even after both reservoirs become affordable.

## Exact grade budgets

A reservoir's grade is its largest strong component order, not its total order. The exact joint minima of \((H,K)\) at strong orders \(3,\ldots,8\), each attained simultaneously, are \((3,9),(5,27),(9,75),(15,201),(25,523),(45,1371)\). The [isomorphism census and output replay](replay_triangle_boundary.py) checks all through order eight. For larger orders the inherited [Busch slot floor](../session327_clock_orbits_and_activation_20261114/OPTIONAL_SLOT_FLOOR.md) and the weak-cut inequality \(q_T\ge(n+1)H_T\) imply \(K_T\ge n q_T/(n+1)\); the weak-cut injection takes every full path and cuts it at each of the \(n+1\) positions. These bounds for the next size already exceed each tabulated \(K\) minimum, and the Busch \(H\) floor increases. Adding ordinal components only increases both positive coordinates: \(q_Aq_B-H_AH_B\ge q_A-H_A\) since \(q_B\ge H_B\ge1\). Consequently the minima also hold for reservoirs of grade *at least* the corresponding orders.

| Minimum reservoir grades | First possible count |
|---|---:|
| 4 and 7 | 14,371 |
| 5 and 6 | 15,345 |
| 4 and 8 | 37,467 |
| 5 and 7 | 39,675 |
| 6 and 6 | 40,851 |
| 6 and 7 | 105,873 |
| 7 and 7 | 274,779 |

That gives a genuine shape in the neighborhood of 30,000: seven vertices are affordable opposite a four-core, and six opposite a five-core, but 5+7 and 6+6 are both blocked. The earlier slope-99 chart at the same target size can use six but not seven in its *one* source. These are different budget constraints, not rival bounds on the entire Hamiltonian spectrum. In particular, actual thresholds don't depend on the location of the primorial.

### Literal near-30,000 constructions

In \([29000,34000]\), among the 4+7 and 5+6 grade pairs, the only prime or \(7p\) outputs are
\[
29167,\quad29477=7\cdot4211,\quad30347,\quad31657,\quad32837=7\cdot4691,\quad32911,
\]
all from 4+7. This includes all padding in the stated grade pairs. From the joint minima, two pads in any allocation already cost more than 34,000. An extra nonsingleton component does too: it multiplies \(H\) by at least 3 and \(q\) by at least 12. Thus inspect the exact core pairs with zero or one pad. The replay independently counts all six assembled tournaments and verifies primality by trial division.

In particular the four-core rows \((2,4,9,3)\), with \(H=5,q=32,K=27\), combine with each of two unpadded seven-cores: rows \((2,12,17,53,99,71,15)\) give \(H=83,q=1144,N=29477\); rows \((6,12,24,65,43,15,55)\) give \(H=89,q=1180,N=30347\). The nearest prime below 30,000 in these lanes is \(29167\), using a neutral singleton appended to the minimal seven-core.

## The real neutral-vertex transfer

A neutral singleton put in a reservoir doubles its \(q\) but leaves its \(H\) unchanged. Moving that singleton from one reservoir to the other at fixed total order then changes the cut output by
\[
N_{\text{pad }W}-N_{\text{pad }U}=H_W q_U-H_U q_W .
\]
The determinant survives the joint-split correction, even though the two output baselines change. It is a way to adjust without paying another vertex, but not a coverage theorem. Both arithmetic phases and the positive absolute budget have to be checked, just as one may activate a square-sum reflection hexagon without having a spanning path on which to use it.

Bounded grade still cannot give a cutoff here. Since \(q_i\ge2H_i\), one has \(q_i\le2K_i\le2N\). With a fixed grade cap all four \(H,q\) coordinates are smooth over the same finite set of primes as in [the single-seam argument](BOUNDED_GRADE.md), regardless of the number of components. This again gives only polylogarithmically many possible outputs through \(X\). The extra degree of freedom opens a real chamber but does not, by itself, fill successive chambers.
