# A finite marked-cover reduction for the complete named response modulo eight

## 1. Statement and retained observer

Let Q be any fixed tournament on b quotient vertices, and let A_1,...,A_b be nonempty tournaments on fixed disjoint physical name sets. Let F(a,c) count Hamilton paths after adjoining two physical ports x->y, with a_u=[u->x] and c_u=[u->y]. The complete unique Boolean multilinear polynomial F has bidegree at most (2,2). Shared physical names at the two ports are retained.

**Finite reduction.** Modulo eight, the complete named polynomial F is determined by the child endpoint-constrained ordered-cover tables defined below, with at most four distinguished path slots and at most seven total paths in a child. More precisely a table with s distinguished slots is needed only for h<=s+3 paths. Unmarked child profiles require only q_1,q_2,q_3. The quotient run skeletons that can survive have at most 3b+4 old runs in total.

For a coefficient with old-name support W, the sharper child cutoff is h_i<=3+|W intersect A_i|, and the total run cutoff is 3b+|W|. Thus shared-index cubic and quartic coefficients are not treated as disjoint-name coefficients. In particular a quartic whose two ports share an old name has at most three distinct old names.

This is a sufficient reduction, not a reconstruction or noncancellation theorem. Congruent retained marked tables imply congruent complete named responses for the same quotient and fixed child name sets. No converse, realizable collision, or efficient algorithm for computing the marked tables of an arbitrary large child is claimed. The number of physical-name entries still depends on the original child orders.

## 2. Marked ordered covers and an exact factorial

An ordered h-cover of A is an h-tuple of nonempty vertex-disjoint directed paths covering A. Designate s distinct slots, and prescribe a physical start name, a physical end name, or both, at each designated slot. There are at most four prescribed endpoint roles in the tables needed here. A slot with the same start and end name is a singleton; incompatible prescriptions have count zero. Repeated physical names in different slots also give zero. Denote the number of covers satisfying the prescriptions by M_A(h;P).

The distinguished slots may be moved to the first s by a fixed slot permutation without changing the count. Permuting the remaining h-s slots is a free action of the symmetric group: the paths are nonempty and disjoint, so no two such slots carry the same path. Every endpoint prescription is preserved. Consequently

    (h-s)! divides M_A(h;P).                           (1)

This is an exact integer statement; it does not rely on tournament parity. Since eight divides 4!,

    M_A(h;P)=0 modulo eight whenever h>=s+4.            (2)

The case s=0 is the usual factorial divisibility of ordered-cover profiles. For one, two, three and four distinguished slots, the respective retained path cutoffs are four, five, six and seven.

## 3. Physical-word substitution

A Hamilton word with two ports has one of the forms

    P0, x, P1, y, P2;          P0, y, P1, x, P2.

The P_i are old directed path segments, possibly empty. In the second form P1 must be nonempty because y->x is absent. Every port has at most one old predecessor and one old successor. Thus there are at most four old contact roles in total.

Split each old segment into maximal child runs. Adjacent runs in one segment have distinct child labels and their compatibility is determined by Q. A port also separates runs even when their child labels agree. For each child, its runs in global word order form an ordered cover. Record which runs carry old contacts, and prescribe the physical endpoint names at those contacts. A run can carry two contact roles. In particular an old singleton between ports supplies the same physical name at both roles, not two independently chosen names.

Fixing this run skeleton and all contact names makes every cross-run and port edge independent of the unprescribed child interiors. The number of physical realizations is exactly the product over children of their marked-cover counts. There is no extra slot factorial: global run order already orders the cover slots. Summing over all admissible skeletons and contact assignments gives the exact full physical response. The attachment factors are precisely a_u or 1-a_u at x and c_u or 1-c_u at y, with factors omitted when a port is adjacent to another port or is an endpoint of the whole word.

There are at most four distinguished slots across all children, since each distinguished slot has a contact role. By (2), a skeleton term is zero modulo eight if h_i>=s_i+4 for any child. Otherwise sum_i h_i<=3b+sum_i s_i<=3b+4. This proves the full-response reduction with no cubic or quartic truncation.

## 4. The coefficient-sensitive refinement

Expand each attachment factor 1-a_u or 1-c_u before extracting a named coefficient. A role whose constant term is chosen imposes no endpoint condition in that coefficient's cover count. The remaining selected roles prescribe only physical names occurring in that coefficient. Sum over all unspecified contact endpoints; this precisely removes their endpoint conditions, with no multiplicity added to a fixed cover.

Each distinct physical name can belong to only one child run. Therefore at most |W intersect A_i| slots in child i remain distinguished. The free action in Section2 now applies to every other slot, even if that slot is next to a port: its selected attachment term is constant, so no endpoint condition is imposed there. This gives the refined cutoff h_i<=3+|W intersect A_i|. Shared-name contact roles must occupy the same run; inconsistent placements contribute zero.

A sum of such signed marked-cover products yields every coefficient of the complete bidegree-at-most-(2,2) polynomial. Its Boolean values, including all nonuniform physical attachments, then follow by evaluation. The reduction does not identify names between different children or replace a child by a smaller unnamed model.

## 5. A literal-query cutoff with few exceptional names

There is also a direct cutoff for an actual query. Suppose the joint port-incidence pattern (a_u,c_u) is constant on a child except at a set of s physical names. Fix a run skeleton and distinguish every slot containing any exceptional name, whether that name is an endpoint or internal to its path. There are at most s such slots. Permuting the other slots preserves every internal, quotient and port edge, since all their names have the same joint port-incidence pattern. The action is free. Therefore the number of Hamilton words with this skeleton is divisible by (h-s)! when h>=s, and every h>=s+4 contribution vanishes modulo eight.

This is a literal-query counting corollary, not an additional endpoint table assumption. It applies to whole four-variable faces with at most s exceptional names in a child, including constant all-zero or all-one port backgrounds. It also shows why uniform probes require only the unmarked q1,q2,q3 profile.

## 6. Other moduli and limitations

For an integer modulus m>=2, put tau(m)=min{t:m divides t!}. The same proof retains h<=s+tau(m)-1 in a child and at most b(tau(m)-1)+4 total runs. For example tau(8)=4 and tau(16)=6. This is again only a finite marked-profile sufficient statistic for a fixed quotient.

The argument is a free slot-permutation action and ordinary maximal-run substitution; no novelty claim is made for either underlying device. Its purpose here is to expose a full named finite-state boundary complementary to the five unmarked modulo-eight profile types. The unmarked five-type theorem applies only to uniform ports. This reduction also retains the endpoint data needed for nonuniform cubic and quartic observers, but does not show that those data force a separator in every prime quotient.
