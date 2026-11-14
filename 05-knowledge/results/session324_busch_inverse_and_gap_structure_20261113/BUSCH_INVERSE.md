# Inverting an ordered Busch overlap

[Busch's merge proof](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/) injects pairs of Hamiltonian paths of strong tournaments overlapping at one vertex into Hamiltonian paths of their union. The repo's [THM-4745](../../../01-canon/theorems/THM-4745-exact-gluings-of-interval-chains-and-the-reduction-of-busch-equality-to-a-one-way-lemma.md) separates non-restricting paths from excess merges and treats interval-chain equality. Here we impose one-way cross arcs, count the entire excess explicitly, characterize zero and unit defects for arbitrary strong factors, and use near-transitive blocks to invert the exact absolute target.

## 1. The exact excess formula

Let strong A,B meet only at v; A0=A−v, B0=B−v, and all cross arcs go A0→B0. Write H_A=H(A), h_A=H(A0), E_A for the number ending at v. Let L_A count ordered two-path covers (P,Q) of A0, empty paths allowed, such that Q is empty or v→first(Q). Thus this is the Hamiltonian path count through v with its outgoing constraint retained and its incoming constraint dropped. Put D_A=L_A−H_A. On B use the dual R_B, S_B, and D_B=R_B−H_B; R_B retains the incoming constraint, and S_B counts paths starting at v.

Every glued path has shape P_A P_B v Q_A Q_B. Its four cases according to whether P_B and Q_A are nonempty give, in the order (yes,yes), (no,yes), (yes,no), (no,no),

    H(T) = (L_A−h_A)(R_B−h_B) + (H_A−E_A)h_B
           + h_A(H_B−S_B) + E_A S_B.

Indeed, if P_B and Q_A are present, only their joins to v need be legal. With just Q_A present, A is legal on both sides and B is a single trailing path; the case with just P_B present is dual. If both are absent, A must end and B must start at v.

Writing e_A=h_A−E_A and e_B=h_B−S_B, expansion gives

    H(T)−H_AH_B = (H_A−E_A) D_B + (H_B−S_B) D_A
                    +(D_A−e_A)(D_B−e_B).

All factors are nonnegative: D_A counts covers whose first path is nonempty and ends in Y=N_A^+(v), with the second empty or starting in Y. Here e_A is the number with second path empty. D_B is dual. By a spanning directed cycle H_A−E_A and H_B−S_B are positive, so exact multiplicative gluing is equivalent to D_A=D_B=0.

## 2. The bottom two defects for arbitrary strong cores

**Claim.** If |Y|≥2, at least two defective covers have both paths nonempty. If Y={y}, then D_A is exactly the number of Hamiltonian paths of A0 ending at y. Thus D_A=0 iff v has one out-neighbor and A0 is not strong. D_A=1 iff v has one out-neighbor y and A0 has exactly one Hamiltonian path ending at y (necessarily A0 strong).

A path can start at precisely the vertices of the initial SCC, and end at precisely the vertices of the terminal SCC, by spanning cycles inside strong components and their total order. Strongness of A forces Y to meet the initial component C of A0: otherwise v could not reach C. Take y∈C∩Y and z∈Y−{y}.

* If z∈C: for two distinct vertices y,z of a strong tournament C, either y can start a Hamiltonian path of C−z, or z can start one of C−y. If y is not in the initial component K of C−z, then z has an arc into K by strongness; K survives removing y and reaches every later component. This places z in the initial component of C−y. This gives a cover whose first path is a singleton and whose second starts at the other Y. The dual argument gives a cover whose second path is a singleton and whose first ends at the other Y. The two are distinct since |C|≥3. Append the subsequent components to the second path in both covers, preserving distinct restrictions to C.
* If z∉C: the first cover is ({z}, a path from y on A0−z). If |C|≥3, take a Hamiltonian cycle of C and let w be the predecessor of y. Then (w→z, a path from y on A0−{w,z}) is a second cover, using the cycle cut inside C−w. If C={y} and another vertex w≠y,z beats z, the same cover works since y is still a source. If C={y} and no such w exists, use ({y}, a path from z on A0−y).
* If Y is a singleton it cannot occupy both endpoints of disjoint paths. Then only the covers counted by e_A remain. Since y lies in the initial SCC, y is a possible end iff A0 itself is strong.

The result with in-neighbors is obtained by converse. Thus the equality-side condition is expressed by one degree and one SCC test, without assuming either core is an interval chain. A defect of one really forces a single marked exit, not two choices with accidentally matching totals.

## 3. A target-first decoder and its reuse cost

Let B=R_j be the order-(j+2) near-transitive tournament, forward along a chain except for its last-to-first edge, with v its first vertex, j≥1. Then H_B=2^j+1, h_B=S_B=1, D_B=0. The exact law is

    N = (2^j+1)H_A+2^j D_A = H_A+2^j L_A.

For a prescribed N, the source is therefore not an arbitrary smaller odd integer. It has to lie on the finite progression

    H_A ≡ N (mod 2^j),   3 ≤ H_A ≤ N/(2^j+1),
    D_A = (N−(2^j+1)H_A)/2^j.

A candidate passing the integer test is a construction only if the indicated rooted profile exists. For j=1, D=0,1,2 determine a unique H for each residue of N mod 3; if these three rooted defect types were available uniformly on smaller attained H, this would give a ternary induction after checking its initial cases. The present classification proves that the D=1 type cannot be supplied by any core with two marked exits. It does not prove the uniform profiles.

For the old root v, h stays fixed, and L transforms to 2^(j+1)L−h: when the A suffix is present all 2^(j+1) splits of the B chain are allowed; otherwise the only forbidden suffix is the singleton final vertex. In particular repeated j=1 steps at a fixed root obey

    h'=h,  H'=H+2L,  L'=4L−h,  D'=(H−h)+2D.

At the newly appended final vertex d the defect is zero: d has the single exit d→v and removing d leaves A followed by a transitive chain. Thus the same operation gives both a low-defect reset at a changing root and high defect at the original root. In particular any obtained count begins an exact ×3 ray, including 63, without requiring 21. The factor closure itself is inherited; the two roots distinguish physically different continuations.

A source of literal unit-defect cores for every j≥1: take A0=R_j with y its last vertex, add v with v→y and all other old vertices beating v. The unique old path ending at y is the forward chain. For a two-cover whose second path begins y, there is one case with singleton second path; otherwise y goes to the first vertex and each of the j intermediate vertices can independently go to either path, in increasing order. Hence L_A=H(A0)+(1+2^j)=2^(j+1)+2, D_A=1 and H_A=2^(j+1)+1. This supplies exact offsets N=(2^k+1)(2^(j+1)+1)+2^k at every j,k≥1. It is an infinite family, not a classification of unit-defect counts.

## 4. A full elementary menu of positive dyadic defects

Continue with A0=R_j, but allow any choice of the unique old vertex y beaten by v. Use forward indices 0,…,j+1. For y at index i with 1≤i≤j, the count of old paths ending at y is 2^(i−1): use the backward edge, then choose which earlier interior vertices lie on its right. To count two-covers with second path starting at y, distinguish the cases of no backward edge, its use in the first path, and its use in the second path:

    C_y = [2^(j+1−i)−[i=1]] + 2^(i−1)3^(j−i) + 2^(i−1)3^(j−i).

In the first case the second path chooses any subset of later vertices. The only false forward cover has first path exactly the two extreme vertices; this is possible only at i=1. In either backward case the earlier i−1 interior vertices have two choices and each later interior vertex has three. Since D_A is the number ending at y in A0 and L_A=h_A+C_y, this gives

    H_A = 2^j+1 − 2^(i−1) + 2^(j+1−i) + 2^i 3^(j−i) − [i=1],
    D_A = 2^(i−1).

At i=0, D=1 and C_y=2^(j+1)−1 (the second path cannot be just the two extremes), yielding H_A=3·2^j−1. At i=j+1, D=1 and H_A=2^(j+1)+1 as above. In particular D=2 at i=2 for every j≥2, giving

    H_A = 3·2^(j−1) + 4·3^(j−2) − 1,
    N = 9·2^(j−1) + 12·3^(j−2) + 1  (attach a triangle).

The three smallest ternary defect types therefore all have unbounded supplies, but the H coordinates are sparse and not interchangeable. This is an exact constructive example of why realizing a digit is weaker than realizing the digit at the prescribed source value.

## 5. Reading the two exceptional values inside this chart

No strong source can make 7: the least target is the two-triangle product 9. For 21, j≥3 is too large; j=2 would require H_A≤4 and H_A≡1 mod 4, impossible. For j=1 the candidates H=3,5,7 require D=6,3,0 respectively. The first two strong source types are C3 and the four-vertex near-transitive tournament (by Busch's floor and the four-vertex classification); their defect menus are {0} and {0,1,2,4}. The last source 7 is globally absent. Thus the single-block chart sees 21 fail by source-value *and* endpoint-defect constraints, not just residue. These tests do not reprove the global exclusion outside this chart.

## 6. Removing the entire search over degree-one cuts

Let A0 have initial strong component C containing y, and tail D (possibly empty). Reconstruct A by adding v with its only exit v→y; this gives a strong tournament provided A0 has at least two vertices and, when D is empty, C is nontrivial. Let W(C,y) count ordered two-path covers of C whose second path starts at y. Then

    h_A=H(C)H(D),  L_A=H(C)H(D)+W(C,y)Q_D(2),
    D_A=0 if D is nonempty, otherwise # Hamiltonian paths of C ending at y.

The formula uses the fact that each path in the cover follows the SCC order: all choices in the tail contribute independently as an arbitrary two-slot cover, and the second slot already starts at y. It removes the search over all subsets for the low-defect side: defect zero is a pointed strong component plus a nonempty serial reservoir; defect one is a pointed strong component with unit end and no reservoir. They are genuinely different inverse charts.

This also disproves the naive three-digit induction at its first nontrivial boundary. A strong core with H=11 has at most five vertices by Busch. If it had D_A=0, the classification would have |C|+|D|≤4 with D nonempty and C strong. For C a singleton the possible H(D)+Q_D(2) are 3,5,9,15 (D is transitive of size 1–3 or C3). If C=C3 the only possible D is a singleton, giving 3+2·3=9. Thus count 11 has no zero-defect source, even though it is attained. The 33 target cannot descend to 11 via the low-digit j=1 rule, but is itself attainable. A working recursion must change charts or allow more than the three smallest defects; it cannot simply fill missing source integers.

## 7. An unbounded strong ray starting exactly at 63

Use the earlier H=63 witness with outgoing masks (238,220,216,16,225,206,8,72) on vertices 0,…,7. Vertex 3 has unique exit 3→4. Deleting 3 leaves 6 a sink (6's only outgoing arc was to 3), so the defect at 3 is exactly zero. Attach a right triangle as above, then use its new last vertex as the next root. This gives explicit strong tournaments of order 8+2k and count 63·3^k for every k≥0. It explains how the first realized member of the exceptional 7-chain can supply its own indefinite strong continuation. It does not manufacture the missing first two members, or place arbitrary prime cores on this ray. This unlimited two-vertex return is compatible with the separately verified failure of a fourth consecutive one-vertex triple extension from a strong 63 seed.
