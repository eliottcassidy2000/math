# Pointed-monoid slot phases, Frobenius correction, and an exact observer collision

**Proved for the specified support-phase model.** This does not construct a field with one element, and it does not turn geometric counts into a Chan–Pak realization rule.

## 1. The recovered notation and its limits

The main historical source is [Session 59's F1 synthesis](../session59_f1_geometry_and_information_20261025/SYNTHESIS.md), especially §4 on split versus geometric unordered quotient counts and §7 on projective cyclic unit supports. Its [cyclotomic support note](../session59_f1_geometry_and_information_20261025/cyclotomic-support/CYCLOTOMIC_SUPPORT_DYNAMICS.md) fixes the pointed monoid \(A_m=C_m\cup\{0\}\) and proves that nonzero projective vectors of exact support \(S\) have \(m^{|S|-1}\) choices after common scaling. The \(m=1\) monoid has one unit, not field addition. In particular actual pointed \(A_1\) support counts must not be identified with evaluating a finite-field polynomial at \(q=1\); the field comparison uses \(m=q-1\), not \(m=q\).

This \(\mathbb F_1\) notation is independent of repository labels such as F1 and F2 for [physical odd-count families in the finite one-arc-join checker](../../../session203/coverage-recertification/README.md), \(F_1(x)=(3x+1)/2\) for a specified Collatz valuation branch in [signed parameter fill](../collatz_signed_parameter_fill_20261007e.md), and \(F_1,F_2\) for polynomial coordinates in a Keller map. None of those labels denotes a base field or monoid model. Earlier [fair-attachment transforms](../session26_random_bridges_20261022/session26-attachments/proof.md) already use nonempty ordered path-cover counts and transitive surjections; the present lift adds unit phases, exact quotient behavior, and a support-preserving collision.

## 2. Slot occupancy as a native rank-one phase model

For an \(n\)-vertex tournament \(T\), let \(F_T(k)\) count ordered spanning covers by \(k\) nonempty directed path words. Over every nonempty slot \(P\) put the projective unit phase
\[
\Phi_m(P)=C_m^{V(P)}/C_m
\]
with the diagonal scalar action. Vertex labels and path words remain marked. This is a rank-one line on its support, not an additive linear transformation. It has exactly \(m^{|V(P)|-1}\) choices by normalizing any one entry.

Every cover therefore supports exactly \(m^{n-k}\) phase choices. Over a finite field with \(m=q-1\), the corresponding scheme is literally the disjoint union, over covers, of products of split tori. For \(r\) ordered slots allowing empties its point count is
\[
Z_{T,r}(q)=\sum_k\binom rk F_T(k)(q-1)^{n-k}.
\]
The binomial factor chooses the occupied slots and preserves their order. For a transitive tournament \(F_T(k)=k!\,{n\brace k}\) because an onto assignment to the slots fixes the internal path words uniquely. At \(A_1\) there is one phase per occupied support, recovering \(\sum_k\binom rkF_T(k)\). In contrast, evaluating \(Z_{T,r}(q)\) at \(q=1\) keeps only all-singleton covers. For exactly \(k\) occupied slots, the coefficient of \((q-1)^{n-k}\), not its unnormalized value at \(q=1\), recovers \(F_T(k)\).

**Ordinal glue retains a missing clutching phase.** In an ordinal join \(T=T_1\to\cdots\to T_j\), restrict a decorated path to each nonempty segment. Suppose a slot visits \(c\) components. For fixed local projective phases the global phase has a fiber of exactly \(m^{c-1}\) choices: choose a representative on each component, scale independently there, and quotient only the one global scalar. The same fiber is the split torus of dimension \(c-1\) over a finite field. Thus after forgetting phase, the empty-slot count multiplicatively factors over ordinal blocks, but for \(m>1\) the global rank-one phase is **not** the product of the projective phases of its segments. The inter-component units are the extra marks.

## 3. Exact geometric quotient counts

First quotient by permutations of the \(k\) nonempty slots. There is no Frobenius correction. All supports are disjoint and nonempty, so each path word is pinned by its support as an element of the unordered collection. Frobenius does not move support labels, and an invariant unordered collection therefore consists of individually invariant phases. The count is \(F_T(k)(q-1)^{n-k}/k!\). This is a concrete contrast to the degree-two and degree-three conjugate packets in Session 59's unlabeled projective configuration problem.

Now let \(G=\operatorname{Aut}(T)\) act on unordered covers and their phases. For any \(g\in G\) and any underlying cover \(D\) fixed by \(g\), its paths break into orbits \(O\). For each such orbit let \(L_O\) be its number of paths, and \(s_O\) the common vertex count of each path. Define a geometric quotient point to be a \(G\)-orbit on geometric decorated covers that is stable under \(q\)-Frobenius. Its count is
\[
\frac1{|G|}\sum_{g\in G}
 \sum_{D:\,gD=D}\ \prod_{O\text{ path orbit of }g\text{ on }D}
 (q^{L_O}-1)^{s_O-1}.
\]
Indeed the number of stable orbits is \(1/|G|\) times the sum of the numbers of solutions of \(\mathrm{Frob}_q(y)=g y\). In a stable orbit every point has \(|\mathrm{Stab}_G(y)|\) such \(g\), and the product with the orbit size is \(|G|\). For a path orbit of length \(L\), the equation propagates the phase all the way around the orbit. The permutation \(g^L\) fixes each ordered position of the path, hence acts as identity on its support. The initial phase therefore lies over \(\mathbb F_{q^L}\), giving the displayed factor. If one instead counts \(G\)-orbits of individually rational configurations, ordinary Burnside gives the **same formula with every \(q^{L_O}\) replaced by \(q\)**. This describes precisely where a nontrivial Frobenius correction can enter.

For \(k=1\), no nonidentity automorphism fixes an ordered Hamiltonian word (it would fix every vertex). Both quotients therefore give
\[
\frac{H(T)}{|\operatorname{Aut}(T)|}(q-1)^{n-1}.
\]
This one-slot rank-one observer cannot distinguish the transitive triple from the directed triangle: their actual counts are \(1,3\) and their group orders are \(1,3\), respectively. Recording group order repairs the one-slot scalar ambiguity, not the all-slot reconstruction problem.

A concrete genuine twist occurs for the ordinal join of two directed triangles, at \(k=3\). Mark vertices as \(a_i,b_i\) for \(i\bmod3\). The three paths \((a_i,b_i)\) give an invariant cover under the simultaneous rotation. Its one path orbit has \((L,s)=(3,2)\), so the geometric twisted contribution is \(q^3-1\) rather than \(q-1\). Every summand in the geometric formula is at least its individually rational counterpart for \(q\ge2\), hence the overall correction is strictly positive in this example.

For this example both quotient polynomials can be written in closed form:
\[
P_{\rm split}(q)=\frac{114(q-1)^3+12(q-1)}9,\qquad
P_{\rm geom}(q)=\frac{114(q-1)^3+12(q^3-1)}9.
\]
To count identity covers, the ordered empty-slot count for a triangle is \(A(r)=r^3+2r\) (its nonempty ordered cover counts are \(3,6,6\)). On the ordinal double, empty-slot counts square, and inclusion–exclusion gives \((A(3)^2-3A(2)^2+3A(1)^2)/6=114\) unordered covers. A nonidentity automorphism moving just one triangle fixes no cover: any moved vertex would need a three-path orbit, which uses up all slots, and none could then contain a fixed vertex. If it rotates both triangles, again the three paths must form one orbit and all have length two. A path inside one triangle would force all three paths into that triangle, so they all cross. There are precisely three equivariant cross matchings for each of the four such automorphisms. Thus the geometric correction is \(\frac43(q^3-q)\), a native count rather than an unspecified nonsplit contribution.

## 4. All-slot collision and a sharp marked reconstruction

**Structural collision.** Permute the strong components of an ordinal join while retaining the internal tournaments and their original vertex labels. Every directed path word in a cover is a concatenation of directed segments from its visited components. Rearrange these segments into the new component order. This gives a bijection on covers, for every \(k\), with inverse the opposite reorder. The support of each path does not change, so leave the original projective phase coordinates unchanged. The map commutes with slot permutations and Frobenius. It also intertwines the automorphism groups, identified with the product of the component automorphism groups: in either total condensation order every component position is fixed by every automorphism. Consequently all the preceding split and twisted phase counts coincide. The collision is not merely numerical; the decorated support models are equivariantly isomorphic as disjoint unions of tori.

For instance compare \(TT_2\to C_3\to TT_1\) and \(C_3\to TT_3\), each on six vertices. The first has two singleton SCCs preceding its triangle, the second has zero, and its converse has three. They are not isomorphic or converse. More generally one may distribute \(d\) singleton SCCs around a single fixed self-converse strong core to get \(\lfloor d/2\rfloor+1\) distinct converse-merged classes with the same observer. This failure survives the entire rank-one family, rather than being repaired by a different finite field, a new unit conductor, or remembering all supports.

**Precise repair.** For \(k=n-1\) every cover has one two-vertex path and otherwise singletons. Every unordered pair occurs in exactly one underlying unordered cover at this layer. Recording the ordered two-vertex word recovers the orientation of that pair and hence the whole labeled tournament. If only support, phase, and slot occupancy are kept, these words lose their direction and the raw \(k=n-1\) data recover just the complete undirected pair set. Automorphism-quotient data can retain some symmetry information, but the explicit SCC reorder shows they cannot restore the lost order in general.

## 5. Consequence for constructive inverses

The useful distinction is between coefficient mass, support, clutching phase, directed path order, and Frobenius descent. A geometric one-unit specialization keeps actual slot states, whereas a finite-field polynomial specialization can erase positive torus ranks. Adding all unit conductors or all fixed-point traces cannot undo a support-preserving SCC reorder. The needed extra observable for constructions sensitive to order is a genuinely directed boundary word or a comparably strong native interface, not another scalar value of the same counting model. No geometric implication for full odd-count coverage is asserted.

The companion `replay_slot_geometry.py` tests the six occupancy layers of the marked collision, and the split/geometric controls \((14,22),(104,136),(346,426)\) at \(q=2,3,4\) for the double triangle, respectively.
