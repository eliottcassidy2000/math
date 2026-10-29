# Signed cuts, rooted charts, and a genuine simple symmetry group

**Status:** elementary PROVED identities and finite exact group certificates. The bijoin/switching framework and abstract switching-group classification are **CLASSICAL**, not discoveries of this session. This is a partial classification layer, not an atlas of every prime tournament.

All cut-matrix ranks in this report are over Q (equivalently R). They are not binary ranks: modulo2 every ±1 cut becomes an all-ones matrix, erasing the signed-cut distinction.

## 1. There is more than one meaning of an atom

Use a skew adjacency matrix S, with S_ij=1 for i→j, S_ji=-S_ij. Vertex switching is S↦DSD for a diagonal sign matrix D. Global converse is S↦-S. Neither operation is the native Hamilton-path complement move.

A proper nonempty subset M is an ordinary module when all vertices outside M see every vertex of M in the same direction. A cut M|M^c is a **signed rank-one cut** when rank(S[M^c,M])=1. These are the classical tournament bijoin cuts. Empty sides are excluded; a nontrivial bijoin has at least two vertices on each side.

**Cut characterization.** A subset M can be made an ordinary module by a vertex switch if and only if its cut matrix has rank one. For a ±1 matrix this is equivalent to every row being equal or opposite to a chosen row. Necessity follows from the cut matrix of a module, and invertible diagonal row/column scaling preserves rank. For sufficiency, multiply each column by the sign of its entry in a fixed row; every resulting row is constant.

**Survival rule.** If M was already an ordinary proper module, it remains one after DSD iff D is constant on M. Indeed an outside row was constant on M, and multiplication by the internal signs preserves that constancy exactly in this case. A cut-rank classification therefore survives switching, while a fixed ordinary modular decomposition need not.

Two exact hostile objects show why the atom vocabulary needs its operation:

- The cyclic5-tournament with differences1,2 is ordinarily prime, but has five nontrivial signed cuts (ten sides).
- A new source vertex joined to Paley7 has an ordinary module consisting of all seven old vertices, but has no nontrivial signed cut.

Thus neither kind of primeness can silently replace the other; their trivial-cut conventions also differ.

## 2. One root converts signed cuts to ordinary modules

Fix r. Put d_r=1 and d_i=S_ri for i≠r, and define S^(r)=DSD. Now r dominates every other vertex. For M not containing r,

    M is a signed rank-one cut in S
    iff M is an ordinary module in S^(r) with r deleted.

The forward implication follows because the normalized cut includes the all-ones r-row, forcing every other rank-one row to be constant. The reverse implication restores that all-ones row and then undoes diagonal scaling.

This gives an exact rooted classification chart: compute ordinary modules of an(n−1)-vertex tournament to recover all signed cut sides avoiding r. It is an elementary matrix version of classical Seidel/bijoin normalization. A different root changes the chart and transports its addresses; it does not justify identifying unmarked attachment positions.

For context, Bui-Xuan, Habib, Limouzy and de Montgolfier's [umodular decomposition work](https://arxiv.org/abs/cs/0604065) develops the corresponding decomposition theory. We use the named framework, while proving the matrix statements needed here directly. No new general decomposition algorithm is claimed.

## 3. Root stabilizers dictate the possible powers of two

For n≥3 define G_sw(S) as the permutations σ for which some signs d_i satisfy

    S_(σi,σj)=d_i d_j S_ij.

Define G_sw^± similarly allowing a common factor ε∈{1,-1}. For a fixed σ, ε is unique when n≥3: otherwise sign ratios e_i would satisfy e_i e_j=-1 for all pairs, impossible on a triple. Thus ε is a genuine group homomorphism on the extended group, with kernel G_sw.

**Root-stabilizer theorem.** The stabilizer of r in G_sw is isomorphic to Aut(S^(r)−r). In the extended group it is Aut^±(S^(r)−r).

Proof: after normalization, a permutation fixing r must satisfy1=εd_r d_i for every i≠r. For ε=1 all d_i equal d_r; for ε=-1 all nonroot signs equal-d_r. In either case the internal map is exactly the asserted ordinary or converse isomorphism. Conversely every such map extends by fixing r with those signs.

Every finite tournament automorphism group has odd order: an involution would interchange a pair and reverse its arc, while an even-order finite group contains an involution. Hence every Sylow2 subgroup of G_sw acts freely on vertices. In particular

    |G_sw|_2 divides n.

For the extended group, a vertex stabilizer has2-part at most2, so every orbit of a Sylow2 subgroup P has size |P| or |P|/2. Consequently

    |G_sw^±|_2 divides 2n.

Here |G|_2 denotes the largest power of2 dividing |G|. These are elementary necessary conditions in a specified action, not a sufficient abstract-group classification.

**Classical stronger result, credited.** Babai and Cameron classify abstract groups realizable as full switching-class automorphism groups by cyclic or dihedral Sylow2 subgroups. Their prescribed-action criterion additionally requires those subgroups to act semiregularly. This is their Theorems5.1/5.2, not a new claim here: [primary paper](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v7i1r38/pdf). Repo THM-4531 already uses this source. The distinction between full abstract realization and containment in a specified permutation action is essential.

## 4. The finite-simple-groups analogy becomes literal in a small example

Ordinary tournament automorphism groups have odd order; their solvability follows from the classical odd-order theorem if that theorem is imported. Switching-class symmetry is a different object and can already be nonabelian simple.

Take P7 on F7 with i→j when j−i is a nonzero square, and adjoin a source∞. This is an eight-vertex tournament S. The checker examines all8! permutations and proves

    G_sw(S)=PSL(2,7), order168,
    G_sw^±(S)=PGL(2,7), order336.

The identification is not inferred from group orders. Independently enumerate every2×2 matrix over F7 of determinant1, respectively nonzero determinant, and its induced Möbius permutation on F7∪{∞}. The resulting permutation sets equal the directly computed signed-automorphism sets. The root stabilizer has order21, exactly Aut(P7). The2-parts8 and16 attain both bounds of §3.

The same test for source+C3 gives PSL(2,3)≅A4 of order12 and PGL(2,3)≅S4 of order24. These are finite exact identifications, not an all-q automorphism theorem.

There is also a standalone finite simplicity certificate for the168-element group. Its conjugacy classes have sizes

    1,21,24,24,42,56.

A normal subgroup is a union of conjugacy classes containing the identity, and its order divides168. Enumerating the32 possible unions shows that the only permitted orders are1 and168. Hence this explicitly constructed group is simple. It is nonabelian, for example because it has nonsingleton conjugacy classes.

This is a genuine bridge to simple groups, but it classifies a symmetry group of a switching class, not the tournament itself, not all switching classes, and not the native metagraph. Those distinctions prevent the analogy from overstating what has been achieved.

## 5. How this fits the recursive program

There are three layers that should be kept side by side:

1. Ordinary modular substitution retains child isomorphism types and quotient symmetries; its signed cycle-index recursion is in the grammar lane.
2. Switching-stable cuts retain rank-one signed interfaces and rooted-chart transporters. Classical bijoin theory controls their decomposition; the root stabilizer makes symmetry computation concrete.
3. Native path complement is a guarded operation on a whole Hamilton path. The metagraph lane recovers the prior obstruction that two such moves cannot imitate a nontrivial vertex switch at order≥8 with labels fixed. Thus switching decomposition is not automatically a transition grammar for that metagraph.

A complete program would have to combine atomic objects, attachment interfaces, automorphism actions and legal transition rules. A list of prime names alone cannot supply this. The examples above show both a precise classical classification that is already available and a structural barrier to importing it unchanged into the native move graph.

## Verification

Run `python -B switching/check_signed_cuts.py` and its optimized counterpart. Both produce identical JSON. The checker covers all1,098 labelled tournaments of orders2–5, all31,668 proper nonempty cut sides, the ordinary-module survival law for the stated deterministic switching controls, signed group actions at orders3–5, and the explicit projective examples at orders4 and8. The matrix arguments prove the unbounded identities; finite enumeration certifies the displayed examples and does not classify all prime tournaments.
