# How many Moore branch matchings can quadratic maps explain?

**PROVED, conditional only on the existence of the graph.** Fix any vertex of a degree-57 diameter-two Moore graph. Assign arbitrary distinct real coordinates to the 56 sheets of each of its 57 branches, independently in each branch. Among the 1,596 unordered pairs of branches, at most **386** of their physical matchings can be written, in at least one direction, as \(x\mapsto x^2+c_{ij}\). Thus **at least 1,210 matchings require a different formula in both directions**. Even one-sided formulas on all pairs of any seven branches are impossible. The parameter and the chosen direction may vary on every edge.

| Permitted formula on a selected matching from \(i\) to \(j\) | Maximum selected pairs | Minimum exceptions |
|---|---:|---:|
| \(x^2+c_{ij}\) | 386 | 1,210 |
| \(F_i(x)+c_{ij}\), with an arbitrary function per sending branch | 456 | 1,140 |
| \(a_{ij}F_i(x)+c_{ij}\), with arbitrary nonzero real edge scales as well | 632 | 964 |

The last row covers edgewise quadratics with one shared critical point per sending branch. The statement places no restriction on the inverse formula or the maps on the exceptional pairs. A general quadratic with a separately chosen critical point for each edge remains outside these bounds. The proof uses only squares forbidden in the graph, first and second moments of coordinates, and a count of pairs. The rational boundedness conjecture plays no role.

## A rectangle of formulas already makes a square

Between two different branches the sheets are perfectly matched. Two different intermediate branches cannot send even one starting sheet to the same terminal sheet, since the two routes would make a 4-cycle. Choose one permitted direction on every pair to which the proposed formulas apply. This gives a partial orientation on 57 vertices. Write \(d_i,e_i\) for its outdegree and indegree, \(m\) for its number of arcs, \(S_i\) for its sheets' coordinates, and \(\mu_i\) for their mean. [Rooted transport identity and its equivalence with the Moore graph](../moore57_switching_flags_20261124/TRANSPORT.md)

Suppose \(i,k\) both send to \(j,l\) via formulas \(F_i+c_{ij}\) and \(F_k+c_{kj}\). Transport from \(j\) back through \(i\) and forward to \(l\) sends \(y\) to \(y+c_{il}-c_{ij}=y+\mu_l-\mu_j\). Transport through \(k\) gives this same translation. This proves that two distinct senders have at most one common receiver. The reverse formulation also gives at most one common sender for two receivers. Thus
\[
  \sum_i {e_i\choose2}\le {57\choose2}=1596.
\]
Since \({e\choose2}\ge 7e-28\) and \(57{8\choose2}=1596\), it follows that \(m\le456\).

Equality implies every indegree is eight: the convex bound is tight at degrees seven and eight, with mean eight. Every pair of senders is on a unique eight-point line given by a receiver's incoming neighborhood. Each sender is then on eight such lines, as its 56 partners occur seven per line. Two lines intersect at most once and counting at a point gives the projective plane axioms. The equality case therefore forces a plane of order seven **on the branch labels**. The incidence bound alone is sharp: translate the repository's Singer difference set in \(\mathbb Z/57\) by three to get
\[
 D=\{2,11,18,22,41,51,53,54\},\qquad D\cap(-D)=\varnothing .
\]
The arcs \(i\to i+D\) form an oriented 456-edge witness for the incidence condition. Obtaining actual matching formulas would also require the sheet permutations. [Session 35 Singer coordinates, trace and Frobenius action](../session35_dimensions_dualities_20261023/arithmetic/SINGER_LINE_FLIPS.md)

For nonzero scales \(a_{ij}\) standardize each receiver coordinate by subtracting its mean and dividing by its positive standard deviation. On a fixed sender the standardized output is now its standardized \(F_i(x)\) multiplied by \(\mathrm{sgn}(a_{ij})\). For two common receivers, the two transports between them coincide if the products of the two edge signs through the respective senders agree. Equivalently, for a given pair of senders the relative signs at all their common receivers must be pairwise distinct. Thus there are at most two common receivers, and \(\sum_i {e_i\choose2}\le3192\). At 633 edges the smallest possible left side already equals \(51{11\choose2}+6{12\choose2}=3201\), proving the third row.

## The second moment makes each collision graph bipartite

Return to the monic square formula. Fix a sender \(i\). Join two different outneighbors \(k,l\) when they in turn send to a common receiver \(j\). Sample \(x\) uniformly on \(S_i\) and let \(Z=x^2\). Since \(i\) has an outgoing matching, all \(Z\)-values are distinct. Along \(i\to k\to j\) the final coordinate is
\[
  Z^2+2c_{ik}Z+c_{ik}^2+c_{kj}.
\]
Its mean is \(\mu_j\). If \(c_{il}=c_{ik}\), equality of means gives \(c_{lj}=c_{kj}\), and then the two transports coincide. For distinct first shifts the requirement that both variances equal \(\operatorname{Var}(S_j)\) instead gives
\[
 c_{ik}+c_{il}
 =-\frac{\operatorname{Cov}(Z^2,Z)}{\operatorname{Var}(Z)}.
\]
Indeed the variance at first shift \(c\) is \(\operatorname{Var}(Z^2)+4c\operatorname{Cov}(Z^2,Z)+4c^2\operatorname{Var}(Z)\). The right-hand side of the sum rule depends only on \(i\). Adjacent first shifts are reflections of each other in the same midpoint and lie on opposite sides. The collision graph is bipartite, with at most \(\lfloor d_i^2/4\rfloor\) edges.

By the rectangle result, \(L=\sum_i{d_i\choose2}\) counts distinct unordered pairs having a common sender, and \(R=\sum_i{e_i\choose2}\) counts distinct pairs having a common receiver. Their intersection has at least \(L+R-1596\) members. Each member is one edge of the collision graph at its unique common sender. Consequently
\[
 \sum_i\left[{d_i\choose2}-\lfloor d_i^2/4\rfloor+{e_i\choose2}\right]\le1596 .
\]
For nonnegative integers the first difference is at least \(3d_i-12\) and the last term at least \(6e_i-21\). These follow by successive differences, or by separating the two parities in the floor term. Hence \(9m-33\cdot57\le1596\), or \(m\le386\).

**FINITE-EXACT local boundary.** The test that every ordered pair has at most two directed two-step routes leaves just one of the 456 tournament types on seven vertices: the Paley orientation \(i\to i+\{1,2,4\}\bmod7\). On it, the outneighbors of 0 are \(1,2,4\); each pair has a common receiver, so their collision graph is a triangle. Thus even the 21 matchings on seven branches cannot all have the one-direction quadratic formula. This local constraint applies in the existing Hoffman–Singleton graph as well as the hypothetical degree-57 graph. At six branches seven of 56 orientations survive the bipartite test. Their coordinate realizability is open here. [Enumeration and two independent boundary examples](quadratic_chart.py)

Two distinct routes can genuinely share a target set over \(\mathbb Q\). For \(j=3,\ldots,30\) take both legs \((j^2-1)/(j^2+1),2j/(j^2+1)\) of the rational unit-circle point. The 56 distinct coordinates start with \(4/5,3/5\). Their squares form a set invariant under \(z\mapsto1-z\) without fixed points. The two length-two routes \(x\mapsto x^2\mapsto x^4\) and \(x\mapsto x^2-1\mapsto(x^2-1)^2\) therefore have precisely the same 56 terminal values and disagree at every starting sheet. Here the Pythagorean parametrization witnesses sharpness of the local variance bound.

## Changing the field changes the obstruction

**PROVED over arbitrary fields.** For additive formulas \(F_i(x)+c_{ij}\), replace the real mean by the translation stabilizer \(H=\{t:S_{j_0}+t=S_{j_0}\}\) of one common receiver \(j_0\) of two senders \(i,k\). For every other common receiver \(j\), both translations \(u_j=c_{ij}-c_{ij_0}\) and \(v_j=c_{kj}-c_{kj_0}\) take \(S_{j_0}\) to \(S_j\). Their difference belongs to \(H\). Transport from \(i\) via \(j\) to \(k\) is transport to \(j_0\) from \(i\), then translation there by \(u_j-v_j\), then inverse transport to \(k\). Different receivers require distinct differences; \(j_0\) uses zero. Thus a sender pair has at most \(|H|\) common receivers. The additive group \(H\) acts freely, so \(|H|\mid56\). In characteristic zero \(H\) is trivial; in positive characteristic \(p\) its order is a \(p\)-power. The bound on \(\sum_i {e_i\choose2}\) is therefore \({57\choose2}\) times 1, 7 or 8. Balancing indegrees in integers gives the bounds on \(m\):
\[
\begin{array}{c|ccc}
  \operatorname{char}&0\text{ or }p\notin\{2,7\}&7&2\\ \hline
  m\text{ at most}&456&1157&1235.
\end{array}
\]
For characteristic two this already forces at least 361 exceptions for completely arbitrary \(F_i\). It explains precisely why one must retain translations invisible to the unordered receiver set in small characteristic.

Over a finite field of characteristic two, \(x\mapsto x^2+c\) is bijective and its second iterate with a different shift is \(x\mapsto x^4+c^2+d\). Eight distinct additive shifts over \(\mathbb F_8\) give eight disjoint routes, so the real two-route test cannot be used there. Nevertheless no complete 57-branch matching chart can have all its edges given in one direction by \(a_{ij}x^2+c_{ij}\), \(a_{ij}\ne0\), in a common algebraic closure of \(\mathbb F_2\). All finite coordinates and parameters lie in one finite field. Each map there is an invertible affine-semilinear map and its inverse is in the same group. This group is solvable by its successive additive, multiplicative and cyclic Galois factors. Its closed-walk subgroup stabilizing one sheet set has a solvable permutation image, whereas Moore holonomy on the 56 sheets is primitive and nonsolvable. The same argument permits edge-dependent Frobenius powers. [Moore holonomy nonsolvability](../moore57_switching_flags_20261124/TRANSPORT.md)

## Which twelve belongs to rational dynamics?

- Poonen's conditional classification assumes that no rational orbit of a quadratic polynomial over \(\mathbb Q\) has exact period at least four. It yields **12 possible graph types** and at most **nine rational projective preperiodic points** for one polynomial, including infinity. Exact periods four and five have been excluded unconditionally. Stoll's period-six exclusion assumes BSD for a specific abelian variety. Recent work on exact period three still treats higher periods as open. [Poonen](https://arxiv.org/abs/math/9512217) · [Flynn–Poonen–Schaefer](https://people.maths.ox.ac.uk/flynn/arts/art11.pdf) · [Stoll](https://www.cambridge.org/core/journals/lms-journal-of-computation-and-mathematics/article/rational-6cycles-under-iteration-of-quadratic-polynomials/5824B91D481A3F10535DBB65EE28F7CD) · [Käenmäki–Kosunen, October 2026](https://arxiv.org/abs/2610.39820)
- A different twelve counts **points**. Manes bounds rational preperiodic points with eventual period at most four in the rational-map family \(az+b/z\), \(a,b\in\mathbb Q^\times\). A twelve-point bound for its full rational portrait follows under exclusion of larger periods. A third twelve counts the conjugacy classes of quadratic postcritically finite rational maps defined over \(\mathbb Q\), under conjugacy over the algebraic closure. [Manes](https://doi.org/10.1112/plms/pdm044) · [Lukas–Manes–Yap](https://arxiv.org/abs/1212.1518)
- The repository's earlier twelve **affine candidates across two maps** consist of the four rationals at \(c=-7/4\) and eight at \(c=-29/16\). For the latter, the complete set is \(\{\pm1,\pm3,\pm5,\pm7\}/4\), with unique rational cycle \(-7/4\to5/4\to-1/4\to-7/4\). At \(c=-21/16\) exactly the same eight rational coordinates instead give two fixed points and a 2-cycle. The trace-one matrix of order six lifting the 3-cycle at \(-29/16\) supplies a different map. This is why the Moore construction needs the full matching, even when first and second moments can rule one out. [The two-parameter twelve](../quadratic_cycle_phase_20261015b.md) · [Denominator-16 classification](../hourglass_triangular_dynamics_20261103/DYNAMICS.md) · [THM-4146, rational three-cycle order-six lift](../../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md)

For the finite statements run <code>python3 -B quadratic_chart.py</code> here with nauty's <code>gentourng</code> on the path. It enumerates each unlabelled branch orientation through seven sites, checks the 56 rational Pythagorean coordinates, and develops the 456 Singer incidences.
