# Cyclic coarse separators with no odd-order ambient correction

**PROVED, subject to the accompanying finite exact certificate:** let three nonempty children all be nontrivial. If at least one child has odd order or is a directed pair, their cyclic quotient converse has a block-uniform physical two-port separator congruent to four modulo eight in every fixed ambient tournament, with no restriction on the total core parity. The selected query is stronger than an arbitrary nonzero local query: all its higher augmented cover-coordinate differences vanish modulo eight.

## 1. All-cover extension of the five-type transfer

Use the five actual canonical representatives and odd units from `../five_profile_types/PROOF.md`. For any fixed quotient, any child-uniform actual port query, and every h>=1, its ordered h-path covers have a maximal-child-run expansion. Order the h path slots, then order the child visits along each path. For each child i the resulting run list is an ordered cover of that child by k_i nonempty paths and contributes exactly q_(k_i)(A_i). Each global skeleton therefore contains exactly one complete-profile coordinate per child.

The global count is consequently multilinear in the complete child profiles, just as the Hamilton count is. If R_i is the actual canonical representative and lambda_i=H(A_i)/H(R_i) modulo8, then

    q_h(augmented Q[A_i])
       = (product_i lambda_i) q_h(augmented Q[R_i]) (mod8)    (1)

for every h. This uses fixed ordered path slots, so paths ending before another path begins are not silently joined, and empty paths are not inserted. It applies to arbitrary nonempty children, although the cyclic theorem below assumes nontriviality.

The scalar in(1) is the same for all h, all queries, and both quotient orientations. It is odd. In particular a canonical contrast with

    (Delta q1, Delta q2)=(4,0) modulo8                       (2)

transfers to exactly the same two residues for the original physical children. Restoring the selected quotient bits onto the whole actual child sets gives the actual query; no physical child is replaced or relabelled in the comparison.

## 2. The complete odd-total certificate

For nontrivial children the possible types are P=TT2, E=TT4, O=TT3, C=directed C3. Exactly32 of their64 ordered triples have odd total order: one or three types are O/C.

For every one of these32 triples, the verifier checks all64 actual coarse queries on both actual canonical quotient orientations. It computes H and

    q2(G)=sum_(nonempty proper S subset V(G)) H(G[S])H(G[V-S]),

the exact number of ordered two-path covers. Every term uses disjoint complementary vertex sets. Thus no product of unrelated marked endpoint counts is substituted for an actual cover.

The complete retained integer table certifies that each of the32 triples has at least one query satisfying(2). It retains every good query index, not merely a sampled witness. The numbers of good queries range from4 to20. The all-O triple has12 good queries; every one of its nonzero modulo-eight coarse queries is good.

The first good query in the fixed order(a,c)=(0,0),(0,1),...,(7,7) is an explicit deterministic witness for each type triple. This finite certificate is exhaustive over all32 odd-total types. The unbounded child orders enter only through the proved transfer(1).

## 3. Vanishing of every higher coordinate

The compared augmented physical cores have the same order. Their order is odd in Section2. The exact odd-order formulas in `../mod8_transfer/PROOF.md` give

    Delta q3 = Delta q2 (mod8),
    Delta qh = 0 (mod8) for h>=4.

Hence the selected query has the full augmented profile contrast

    (Delta q1, Delta q2, Delta q3,...)=(4,0,0,...) (mod8).   (3)

This is why the odd-order Hamilton-cycle correction cannot cancel the separator. It is not enough to start with an arbitrary local residue-four query. The inherited TT2,TT2,TT3 example at query(1,2) has a nonzero q2 contrast and loses its Hamilton residue in a cyclic exterior; that boundary is retained here alongside the selected good query for the same type triple.

For the31 even-total nontrivial type triples other than(E,E,E), the earlier complete five-type certificate supplies a query with Delta H=4. The even-order formulas automatically give Delta q2=2Delta H=0 and Delta q3=0 modulo8, with all higher coordinates zero. Thus(3) holds for those selected queries as well.

## 4. Arbitrary fixed ambient, all total parities

Place the core in the root r of any fixed ambient tournament R. Give both actual ports the root's original incidence to every exterior vertex. The union of the core with its two ports is a module in both actual augmented graphs. Its run substitution gives

    Delta F = sum_(h>=1) c_h(R,r) Delta q_h,
    c_h>=0, c1=H(R),

with the same coefficients on both sides. Applying(3),

    Delta F = 4H(R) = 4 (mod8),                            (4)

because H(R) is odd. This proves the stated arbitrary-ambient theorem for all63 nontrivial type triples except(E,E,E), including the32 odd-total cases.

The excluded type means all three child orders are even and at least four. Its whole coarse contrast bank is zero modulo8 by the separate five-type classification. It may still have nonzero integer coarse contrasts, and named nonuniform queries may still separate it. No exact collision is asserted.

The target may preserve the source child interiors or replace them by tournaments with the same complete child profiles. The two-port query is uniform on each whole child; no individual marked-name reconstruction is claimed.

## 5. Original-module root families

The independently audited `../monotone_dual_queries/PROOF.md` now applies without any common-profile or size assumption beyond the present type condition. Given the selected baseline query z, pair it with its core converse-complement Jz=(1-c,1-a); keep both exterior port incidences fixed at the root incidences. Their entire augmented profile contrasts are exact negatives before reducing modulo8. One has negative baseline source-minus-target contrast.

For every original-module XOR family containing the cyclic root, proper selected modules lie within the three children. Hereditary original-module monotonicity makes each target child cover profile coordinatewise no smaller than the source profile, even with crossing descendants. Positive run substitution therefore preserves the negative direction at that one member of the pair.

Consequently the two fixed actual queries {z,Jz} detect every root-containing family for the stated source, also in any fixed ambient, and one member chosen from the known baseline signs works simultaneously for every such family. No all-zero scalar fallback is needed. This is source-relative family detection, not pairwise reconstruction of arbitrary orbit states. Exterior changes are excluded.

## 6. Verification and scope

The verifier binds the complete frozen local classification, retains the full32×64 integer(H,q2) table, compares every H value to the earlier physical certificate, checks both modes, and includes original larger-child and arbitrary fixed-exterior controls. It retains the bad odd-total local probe and its nonzero integer but zero-modulo-eight ambient gap as a sharp warning against the invalid automatic lift.

The full named cubic/quartic data remain preserved in the pinned `five_profile_types` packet. This extension does not replace them by a low-degree truncation; all new Hamilton values are literal counts of the whole physical query.
