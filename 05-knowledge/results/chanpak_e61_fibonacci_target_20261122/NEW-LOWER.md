# More than a million strong lower witnesses at exponent 61

**PROVED, with finite-exact arm data.** The height \(E_{61}-1\) has at least **645,120 strongly connected isomorphism classes of order 107** and **737,280 of order 108**. These are different from the two previously identified order-124 classes. In a subfamily of 720 of the order-107 classes, examining every one of 4,083,120 physical arc reversals gives no \(+2\) exit. A separate four-layer obstruction proves that none of the **1,382,400** classes can reach the upper endpoint by a single flip. For every exponent \(p\ge5\), the obstruction also rules out all shorter words that use a transitive arm of width \(p\) to realize the lower target. An unrestricted straddle at exponent 61 remains open. [Original three-port lower](../chanpak_support_birth_grammar_20261122/ODD-PERFECT-PORT.md) · [Its neutral component](PLATEAU.md)

## A product with a recoverable root

For a nonempty rooted tournament \((V,b)\) write \(h(V)=H(V)\), \(Q(V)=\sum_{S\subseteq V}H(S)H(V-S)\) with \(H(\varnothing)=1\), and \(s_b(V)\) for the number of ordered two-slot covers whose first path is nonempty and ends at \(b\). Splitting any Hamiltonian path at \(b\) injects into these covers, so \(s_b(V)\ge h(V)\). Let \(A\) be any nonempty tournament. Orient \(A\to V\), add a root \(a\) beating every old vertex except \(b\), and put \(b\to a\). Denote the result by \(K(A,V,b)\). Then
\[
 H(K)=h(A)h(V)+Q(A)s_b(V),\qquad s_a(K)=H(K).
\]

For the first identity use the odd-cycle packing formula. Packings avoiding \(a\) contribute \(h(A)h(V)\). The cycle through \(a\), when present, has a path on an arbitrary subset of \(A\) followed by a nonempty path ending at \(b\) in \(V\); the combined path must have even order. Conversely, every such even path completes to the cycle. When its \(A\) part is empty, even order forces at least two vertices in \(V\), hence its start is different from \(b\). The path count on its complement splits across the same ordinal cut. Exactly half of the \(Q(A)\) ordered covers have even first-slot size and half have odd first-slot size. Indeed the signed sum is zero: expand \(H(S)H(A-S)\) as two disjoint odd-cycle packs. For fixed cycles and uncovered singleton vertices, each cycle and singleton can independently lie on either side, giving opposite signs since all blocks have odd size and at least one block exists. The weight two of the cycle through \(a\) cancels the half. [THM-002, odd-cycle formula](../../../01-canon/theorems/THM-002-ocf.md)

For the second identity, a cover ending its first path at \(a\) fails to concatenate only if the second path starts at \(b\). In that case the first path must consist just of \(a\), since \(b\) is its only possible predecessor. The second would then start a spanning path of \(A\to V\) inside \(V\), which is impossible. For strong \(V\) (including a singleton), \(K\) is strong: all old vertices can reach \(b\to a\), and \(a\to A\to V\) gives reachability in the other direction.

Starting from a singleton, the equality \(s_b=h\) therefore persists and each new layer multiplies the count by \(R(A)=h(A)+Q(A)\). Although \(R(A)\) is the familiar path count of the triangle wrapper, the marked construction realizes **ordered products inside strong tournaments**. It retains an apex where the mechanism may be iterated. [Prior triangle-wrapper transfer](../radix_three_shields_20261116/TRIANGLE_WINDOW.md)

Whenever \(|A|\ge2\) the unmarked isomorphism class recovers the layer. Every vertex of \(V\) has indegree at least \(|A|\), whereas \(a\) has indegree one. A second indegree-one vertex exists precisely when \(A\) has a source. In that case the one beating the other is \(a\). Its sole inneighbor identifies \(b\), and after deleting \(a\) the terminal strong component is exactly \(V\), leaving \(A\). Induction recovers the ordered word of actual arm classes. In particular permutations of distinct arms give distinct isotypes.

## The two orders and the complete class count in this grammar

Use a transitive arm of width 61 with \(R=2^{61}+1\). The five other example factors satisfy
\[
 1353\cdot3963\cdot5551\cdot4681\cdot8275=2^{60}-1,\qquad
 (2^{61}+1)(2^{60}-1)=E_{61}-1.
\]

| \(R\) | Order | \(h\) | \(Q\) | Outgoing row masks |
|---:|---:|---:|---:|---|
| 1,353 | 7 | 93 | 1,260 | 8, 5, 73, 50, 39, 71, 27 |
| 3,963 | 8 | 219 | 3,744 | 6, 12, 72, 129, 15, 31, 59, 119 |
| 5,551 | 8 | 335 | 5,216 | 6, 20, 136, 67, 45, 79, 151, 59 |
| 4,681 | 8 | 257 | 4,424 | 10, 12, 65, 164, 15, 23, 59, 119 |
| 8,275 | 8 | 563 | 7,712 | 136, 37, 73, 82, 39, 141, 51, 94 |

These six arms give \(6!=720\) classes on \(1+6+61+7+4\cdot8=107\) vertices. To obtain the larger counts, allow any arms of orders two through eight together with exactly one transitive 61-arm. The complete small catalogue contains 243 arm isotypes with \(R(A)\mid 2^{60}-1\), occupying 66 different pairs \((|A|+1,R(A))\). Let \(D(X,z,t)\) count ordered words of the small arms having product \(X\), marking total order cost with \(z\) and length with \(t\). Set \(D(1,z,t)=1\). For \(X>1\) sum \(z^{|A|+1}tD(X/R(A),z,t)\) over all actual types with \(R(A)\mid X\). The exhaustive recurrence gives
\[
 D(2^{60}-1,z,t)=t^5(107520 z^{44}+122880 z^{45}).
\]
Inserting the long arm has six positions. It and the initial singleton add 63 vertices. Unique decoding proves the announced counts. The extra types and orders here give no upper bound on the full count fiber. [Ordinary class catalogue and recursion](../chanpak_arc_fiber_grammar_20261122/README.md)

## Sampling the class fiber and counting its singleton extensions

The word recurrence also gives an exact uniform choice of an **isomorphism class** in either displayed order. At the small-word state \((X,b,j)\), choose a particular next arm \(A\) with probability
\[
 \frac{[z^{b-|A|-1}t^{j-1}]D(X/R(A),z,t)}
      {[z^bt^j]D(X,z,t)}
\]
when its factor divides \(X\), assigning zero to impossible choices. Start at \((2^{60}-1,44,5)\) or \((2^{60}-1,45,5)\), then choose the large arm's position uniformly among six positions. The recurrence normalizes the probabilities; their product telescopes to the reciprocal of the certified class count. Equivalently partition an interval of that many integer ranks in proportion to the continuation counts. The reverse decoding of the unmarked roots makes each completed word one class. This illustrates the role of a normalized reverse history in the [pinned Latin-square sampling writeup](https://github.com/Mira-Cult/math/blob/30ddaaad27b48df2002225c3e3f40cffd02c5253/notes/combinatorics/exact-latin-square-sampling/README.md). For these tournament words the reverse history is unique.

Prepending and appending an arbitrary run of singleton components gives a contribution to the all-size generating function for \(E_{61}-1\). Combine the two orders above with the disjoint [90-vertex](../chanpak_e61_compact_volchkov_20261122/STRONG-90.md) and [enlarged 91- and 93-vertex](PROFILE-FIBERS.md) families. The following is a coefficientwise lower bound, and is the exact generating function of this particular subfamily:
\[
 \frac{2{,}349{,}088{,}560z^{90}+414{,}387{,}118{,}080z^{91}
 +30{,}453{,}618{,}843{,}648{,}000z^{93}
 +645{,}120z^{107}+737{,}280z^{108}}{(1-z)^2}.
\]
In our two six-arm orders, taking the converse reverses the order of the arms and takes each arm to its converse. The unique 61-arm can never remain in the same position of an even-length word. Unique decoding and closure of the small bank under converse therefore give exactly 322,560 and 368,640 pairs of distinct converse classes in these two families. The single strong component in a padded word recovers its order and class, and the two singleton runs are recoverable independently. The inherited strong-middle theorem implies no member of this padded family can reach the upper endpoint with a single edit: any such pair strips to the strong lower class with a corresponding upper exit. [SCC and marked-pair series](../chanpak_arc_fiber_grammar_20261122/README.md)

## Every arc of the 720 explicit permutations

The quotient has one symbol per arm and one singleton per root. In the order outer root, outer arm, next root, next arm, and so on until the innermost singleton, all intersymbol arcs point forward except those between consecutive roots. To return to a given arm after leaving it, a skeleton path must cross back over that arm by the unique backward arc spanning its position. Thus each arm occurs in at most two visits on a single path and at most three across two disjoint paths; a fourth visit would require using this same root arc twice.

For each pair of physical vertices \(u,v\), enumerate two nonempty path slots partitioning the full vertex set, one ending at \(u\), the other starting at \(v\). Concatenation counts exactly the old paths using \(u\to v\) if this arc is present. Subtracting this number from the opposite ordered boundary count gives the exact change on reversal. For the quotient enumeration we use the prior disjoint-mask kernel. To fill \(r\) ordered visits to an arbitrary arm, partition its vertices into \(r\) nonempty sets and multiply the path counts on those sets. The required end and start of boundary visits are pinned by the ordinary last-vertex recurrence or its converse on their respective sets. When both boundaries are in one arm they belong to distinct visits. Its transitive specialization is the increasing-rank filling formula. The replay compares every physical arc of a small two-arm host (one cyclic arm, one transitive) to direct Hamilton-path DP. [Original two-slot quotient kernel](rail_arcs.py)

Across the 720 chosen lower classes the signed changes are negative 96,600 times, zero 45,600 times, and positive 3,940,920 times. The smallest positive difference is 10,648,826,675,531,779,512,729,600. These 4,083,120 exact tests concern all arcs of the 720 words in the table. The general proof below covers the other words without testing individual arcs.

## No four-layer product can exit by two

**PROVED.** Start at a singleton \(V_0=r_0\) and iterate \(V_j=K(A_j,V_{j-1},r_{j-1})\) with root \(r_j\) and arbitrary nonempty tournaments \(A_j\). For **every** word of length \(k\ge4\), of any orders, no arc changes its height by \(+2\) or \(-2\). At length three this fails: taking all arms to be singletons gives a strong seven-vertex host at 27 with two arcs increasing its count to 29. The obstruction thus excludes every word in the 107/108-vertex family above.

Two divisibility tests apply to opposite ends of the word. First, if both flipped endpoints are inside \(V_{k-2}\), the next two outer layers are unchanged. The first of those produces a root with \(s=H\) for either orientation, and the second makes the new and old count divisible by \(R(A_k)>1\), an odd integer. A change of two is impossible.

Second, put \(W=V_2-r_2=A_2\to V_1\), with distinguished vertex \(d=r_1\) inside \(V_1\). Suppose neither flipped endpoint belongs to \(W\). After deleting \(r_2\) the rest still splits \(L\to W\). The vertex \(r_2\) beats all of \(W\) except \(d\); its contacts with \(L\) can be arbitrary. Split a Hamilton path at \(r_2\) into two slots of \(L\to W\). If the first slot uses no vertex in \(W\), the second contains a full path of \(W\). It cannot start that path at \(d\) because \(A_2\) is a nonempty initial segment. Thus its contribution, for each permissible assignment and ordering in \(L\), is exactly \(H(W)\). If the first slot uses \(W\), it must end at \(d\) before entering \(r_2\); the second cannot start at \(d\) because \(d\) has already been used. In this case the \(W\)-contribution is exactly \(s_d(W)\). Both possible contributions are multiples of \(H(V_1)>1\), since
\[
 H(W)=H(A_2)H(V_1),\qquad s_d(W)=Q(A_2)s_d(V_1)=Q(A_2)H(V_1).
\]
This common odd divisor survives any edit avoiding \(W\).

Accordingly one endpoint must lie in \(W\) and one outside \(V_{k-2}\); these are disjoint sets when \(k\ge4\). Place \(r_j\) at level \(j\) and the arm \(A_j\) at level \(j-1/2\). All arcs between distinct levels point from higher to lower, except the unit increases \(r_{j-1}\to r_j\). At \(k\ge5\) a directed path of length two cannot join the inner endpoint to the outer: their level difference is at least two, with equality only between two half-integers, whereas a two-step increase by two requires three consecutive integer roots. Their outside neighborhoods are therefore nested. The classification of nested-contact unit flips permits just four lower heights, the largest 13. Here \(H(V_k)\ge3^k\). [Nested-contact theorem](../chanpak_cech_perfect_seams_20261122/NESTED-CONTACTS.md)

At \(k=4\), the inner levels are at most \(3/2\) and the outer levels at least \(5/2\). An outward two-step path must use two successive root increases to bridge these sets: one increase with a remaining level-preserving step cannot start and end at the requisite half-integers. The only possible pair is \(r_1,r_3\), originally \(r_3\to r_1\). The difference on reversing it has a positive expansion in four-arm nonempty path-run counts. For each arm let \(B_t(A_i)\) be its number of ordered covers by exactly \(t\) nonempty paths. The coefficient of \(\prod_{i=1}^4 B_{t_i}(A_i)\), in lexicographic order on \((t_1,t_2,t_3,t_4)\in\{1,2\}^4\), is
\[
 (60,22,33,11,33,12,18,6,22,8,12,4,11,4,6,2).
\]
All other coefficients vanish. These sixteen integers are counts of two-slot quotient skeletons with ordered boundary \((r_1,r_3)\) minus those with \((r_3,r_1)\), by exactly the disjoint-mask recurrence above on the nine-symbol quotient. The first term alone guarantees a gain of at least 60 for every choice of arms.

## An upper obstruction independent of the inner tournament

Even allowing arbitrary \(V,b\), \(H(K(\mathrm{TT}_p,V,b))=E_p+1\) is impossible for every \(p\ge3\). Indeed it forces \(h(V)\equiv2^{p-1}+1\pmod{2^p}\). But \((2^p+1)h(V)\le E_p+1=(2^p+1)(2^{p-1}-1)+2\). The smallest positive candidate for the residue exceeds this upper bound on \(h(V)\). Consequently a lower word ending in the 61-arm cannot gain two by any modification leaving that outer interface in the stated form. Two of the other example arms, of factors 1,353 and 3,963, have both \(h(A)\) and \(Q(A)\) divisible by three; their entire output is likewise divisible by three, whereas \(E_{61}+1\) is not.

Run `python3 -B NEW-LOWER.py --class-count` for the complete word recurrence and `python3 -B NEW-LOWER.py --all-orders` for all arcs on the 720 specified words. The default checks one 107-vertex lower class. [Replay](NEW-LOWER.py)

## No perfect straddle with a width-\(p\) transitive factor in a product word

**PROVED.** Let \(p\ge5\), with no primality assumption. Suppose a word from the singleton has lower count \(E_p-1\) and some arm is \(\mathrm{TT}_p\). No arc reversal reaches \(E_p+1\), whatever the number, sizes or internal orientations of the remaining arms. This applies to arbitrary positions of the long arm.

The four-layer theorem reduces this to two or three layers. One layer with the specified arm has the wrong count. An internal arc of any arm changes only its factor, and another odd factor greater than one survives. Distinct arm blocks have nested contacts. An arm vertex can have an outward two-step path to or from an outer level at distance at most one: one of the two steps touches the arm and cannot increase level. Thus the only potentially nonnested arm/root edges have the bounding roots as endpoints. At three layers the cuts in the previous proof remove both end arms and leave just the middle one. At two layers either arm is possible. Reversing all arcs and reading the word in the opposite order yields the same construction with the converse of each arm. Hence it suffices to consider flipping the edge from the higher bounding root to a vertex \(u\) of the candidate arm \(B\).

The remaining root/root cases also have large gaps, irrespective of the arms. As before take differences of the two boundary skeleton tables, in lexicographic run order, and write roots as \(r_j\). At two layers the coefficients at the three original arcs \(r_0\to r_1,r_1\to r_2,r_2\to r_0\) are respectively
\[
 (-6,-2,-3,-1),\quad(-6,-3,-2,-1),\quad(4,3,3,2).
\]
At three layers the only cases surviving the two end cuts and nested-contact test are \(r_2\to r_0,r_3\to r_1,r_1\to r_2\). Their respective coefficient arrays are
\[
 (16,6,9,3,11,4,6,2),\quad(16,11,9,6,6,4,3,2),\quad
 (-18,-6,-9,-3,-6,-2,-3,-1).
\]
Each has one strict sign throughout and first term of absolute value at least four.

**The endpoint equation.** Write \(F_L,h_L\) for the factor \(R(A)\) and count \(H(A)\) of the arm immediately inside \(B\), with both equal to one if absent. Use \(F_R,h_R\) for the outer arm in the same way. Only these arms can remain in the cases under consideration. In \(B\) write \(t,e\) for full paths starting or ending at \(u\), and \(s,m\) for all ordered two-slot covers whose first path starts or ends at \(u\), respectively. The signed change on flipping the higher root's arc to \(u\) is
\[
              \Delta=F_L\{F_Rm-h_R(s+t)\}-h_Rh_Le.                 \tag{1}
\]
Indeed the left inner tournament has count and marked ending mass both \(F_L\); its number of full paths starting at its root is \(h_L\), since removal of that root leaves the left arm followed by the initial singleton. Both statements include the no-left case. At the root above \(B\), the new paths acquired by permitting \(u\) to end the first slot number \(F_Lm-h_Le\): when the second slot starts in the inner tournament its sole forbidden start is its root. The paths lost by disallowing \(u\) as a start are \(F_L(s+t)\): either the first slot ends at the inner root or the first slot is empty. The change in marked ending mass at the higher root is \(F_Lm\), since it now also accepts preceding paths ending at \(u\). The outside arm weights the first change by \(h_R\) and the mass change by \(Q(A_R)=F_R-h_R\). This gives (1), including an absent outer arm.

Set \(M=p-1\), \(J=2^M-1\). We have \(F_LF_BF_R=(2^p+1)J\) with \(F_B=H(B)+Q(B)\). If \(F_L\) is the specified transitive factor, then \(h_L=1\), and (1) is an integer multiple of \(2^p+1\) minus \(h_Re\), which lies in \([0,J]\). It cannot equal two. If \(F_R\) is the transitive factor, \(h_R=1\), \(m\ge1\), \(s+t\le F_B\), \(e\le H(B)\), and \(h_L\le F_L\). Therefore \(\Delta\ge F_L(2^p+1-F_B-H(B))>2\) since \(H(B)<F_B\le J\).

It remains to put \(B=\mathrm{TT}_p\). At rank \(0\le j\le M\) we have \(m=2^j,\ s=2^{M-j},\ t=[j=0],\ e=[j=M]\). At the last rank (1) is at least \(F_LF_R(2^M-2)>2\). At all earlier ranks \(F_L\) divides (1), excluding equality if \(F_L>1\). Otherwise \(F_R=J\). At rank zero the difference is \(J-h_R(2^M+1)<0\). At ranks \(2\le j\le M-2\) it is divisible by four. Rank \(M-1\) gives \(2^{M-1}J-2h_R>2\). Equality at rank one would require \(2^M-2=2^{M-2}h_R\), impossible modulo four for \(M\ge4\). This proves the all-length exclusion.

## A full contact face of the older 124-vertex lower

In the eight-symbol source of the original lower, vary all \(p\) arcs from singleton 7 to port 6 of width \(p\), keeping the remainder fixed on the odd-exponent balanced rail. Write \(x_i\) for reversal of its initial arc to rank \(i\), and \(B=2^{p-1}-1,\ Z=2^{(p-3)/2}\). The count difference from the original lower obeys
\[
               \Delta H\equiv (1-3Z)x_{p-1}\pmod B.
\]
Removing singleton 7 leaves the ordinal tournament \(\mathrm{TT}_p\to U\), with \(H(U)=B\) by the [three-port lower formula](../chanpak_support_birth_grammar_20261122/ODD-PERFECT-PORT.md). In \(U\) the six labels of the proof are \(0,\ldots,5\), with 4 and 1 inflated; the numbers of spanning skeletons with run counts \((1,1),(1,2),(2,1),(2,2)\) are \(15,11,5,4\). Expanding gives \(-1-3Y+3Z+4YZ\), which equals \(4Z^2-1=B\) on the diagonal. Every mixed coefficient of this contact cut is divisible by \(B\): the even paths between two vertices of the transitive segment stay there, leaving a complement with count \(H(U)\). The [single-rank reflection formula](PLATEAU.md) gives linear coefficients divisible by \(B\) except at the last rank, where \(D-B=-5-3(Z-2)=1-3Z\). The inherited [quadratic cut formula](../chanpak_minimal_core_grammar_20261122/README.md) finishes the claim. At exponent 61 this excludes the upper height anywhere on the whole \(2^{61}\)-vertex face of contact words. In particular an exact lower word in this face must preserve its last contact.

## Trying to replace the long transitive arm

Use the factorization \(E_{61}-1=32767(2^{15}+1)(2^{30}+1)(2^{61}+1)\). If \(Q(B)=2^sH(B)\), an ordinal prefix of \(m-s\) singleton vertices gives an arm with \(R=H(B)(2^m+1)\) whenever \(m\ge s\). Three such arms with ordinary heights multiplying to 32,767, including a nontransitive choice for \(m=61\), would give a three-layer lower outside the preceding exclusion. This criterion itself does not guarantee a unit exit.

**FINITE-EXACT.** No nontransitive tournament of order at most nine with \(H(B)\mid32767\) has the required dyadic ratio. The order-nine census has 191,536 classes. Of those, 348 have nontrivial height dividing 32,767: there are 14 at 31, 35 at 151, 63 at 217 and 236 at 1,057. For a complete replay from this directory, compile and run the following, with nauty's tournament generator in its default format.

~~~sh
c++ -O3 -std=c++17 NEW-LOWER-DYADIC.cpp -o /tmp/chanpak_dyadic_arm
for n in 2 3 4 5 6 7 8 9
do
    gentourng -q "$n" | /tmp/chanpak_dyadic_arm "$n"
done
~~~

The checker reconstructs every subset's Hamilton count by its first vertex and computes \(Q(B)=\sum_S H(S)H(B-S)\). The ordinary census at each size is asserted. [Exact dyadic arm test](NEW-LOWER-DYADIC.cpp)

This search is also conclusive for **all orders within the one-apex family**. Remove leading zero and trailing one contacts of the transitive spine. These are exterior singleton components; each removal preserves \(H\) and divides \(Q\) by two, thus preserving the dyadic-ratio predicate. For a remaining nontransitive word the binary numeral \(v\) begins with one and ends with zero. The ordinary count is \(1+f(v)\); the last zero alone contributes at least \(\lfloor v/2\rfloor\ge2^{w-2}\) to \(f(v)\) at width \(w\). Thus \(H\le32767\) forces \(w\le16\). Among all positive even addresses below \(2^{16}\), exactly 11 have nontrivial count dividing 32,767: one at 31, two each at 151, 217 and 4,681, and four at 32,767. None has dyadic \(Q/H\). Run `python3 -B NEW-LOWER.py --dyadic-spines` for the complete word sweep. [Inherited binary and two-slot reads](../chanpak_e31_physical_20261122/ALL-WORDS.md)
