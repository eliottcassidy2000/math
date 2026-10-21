# Gluck regular subsets, Fano symmetry, and tournament sidecars

Status: CITED paper statements; PROVED elementary transfers relative to Gluck's Corollary 1 and the stated repository factorization; FINITE-EXACT receipts. No literature-priority claim. Prepared outside the repository; no Git operations or repository modifications.

## Useful conclusions first

1. Every finite tournament has a binary vertex coloring with trivial ordinary automorphism stabilizer. This is a classical consequence of Gluck, not a new result.
2. Every finite tournament with at least two vertices has a three-color vertex coloring with trivial stabilizer even when anti-automorphisms are admitted. Three are sometimes necessary: the seven-vertex Paley tournament.
3. For any particular ordered partial path, binary colors on its omitted vertices can kill its residual ordinary stabilizer, even when three or more vertices are omitted. This does not make the original uncolored path action free and does not replace its Burnside correction.
4. In the inherited disjoint-action decomposition Aut(T)=H x C3^t, with H faithful on exceptional B and u fixed singleton vertices, the exact number of regular binary subsets is R_H 6^t 2^u. The minimum selected-class size is rho(T)=rho(H)+t <= floor(|B|/2)+t <= delta+t when |B|<=2delta. Thus bounded symmetry slack does not bound the size of a symmetry-breaking subset.
5. On the same seven-point carrier, the ordinary Paley group, its converse extension, and the full Fano group have distinguishing numbers 2, 3, and 4, respectively. The missing coordinate is the allowed group action, not the number seven.

## 1. Actual uploaded paper and exact scope

David Gluck, *Trivial Set-Stabilizers in Finite Permutation Groups*, Canadian Journal of Mathematics 35(1) (1983), 59-67. DOI: https://doi.org/10.4153/CJM-1983-005-2 . The 2018 date on the publisher landing page is an online-publication date, not the original paper year.

The supplied PDF is this group-theory paper, not a Fano/projective-geometry paper. It does not itself formulate a tournament application or mention the Fano plane. The application below is a deduction, and the basic tournament two-color consequence has independent prior literature.

PDF SHA-256: b9b484ab2e28417b3f6f2c696862e10cd73ebe36569c2037eda58a4320d3f01c . All nine pages were extracted; printed pages 65-66 (PDF pages 7-8) were visually checked for theorem hypotheses, exception symbols and powers in Proposition 1.

### Terminology

For a faithful permutation group G on a finite set Omega, a subset S is regular if its **setwise** stabilizer G_S is trivial. This is equivalent to a regular orbit of G on the power set and to a distinguishing coloring with two named colors (S and its complement). It is not the same as a pointwise base. A color swap is not allowed in this definition.

Gluck's (*) means existence of a regular subset. His (**) adds |S| != |Omega|/2, equivalently a regular subset with |S|<|Omega|/2 after complementing.

### Theorem 1, printed page 65: primitive solvable groups only

Let G be finite, primitive, solvable on Omega, with point stabilizer M. Then (*) holds except exactly in these cases:

- degree 3: |M|=2;
- degree 4: |M|=3 or 6;
- degree 5: |M|=2 or 4;
- degree 7: |M|=6;
- degree 8: |M|=21;
- degree 9: |M|=16, 24 or 48, or M is isomorphic to D8 (the order-eight dihedral group).

When (*) holds and the degree is not 2, (**) also holds. This is not a classification of arbitrary primitive groups or arbitrary solvable imprimitive groups. The faithful primitive solvable group has affine form K semidirect M, where K is regular elementary abelian and M acts faithfully irreducibly on K.

For the degree-seven exception, G is AGL(1,7), order 42. Lemma 6, printed page 64, explicitly contrasts its absence of regular subsets with the order-21 group, for which every two-element subset is regular, and the order-14 dihedral group, which has a regular orbit of fourteen three-element subsets.

### Theorem 2, printed pages 65-66: a sufficient imprimitive lifting criterion

Let G be transitive solvable with point stabilizer H. Suppose a chain H=H0 < H1 < ... < Hn=G has each subgroup maximal in its successor; each index [Hi:Hi-1] is not 2; and each induced primitive permutation group of Hi on the cosets of Hi-1 has (*). Then (G,Omega) has (*).

The induced groups are permutation images, not necessarily faithful copies of Hi. The proof uses non-half-size regular subsets in the blocks and complements them according to a regular subset of the block quotient. Cardinality distinguishes chosen versus complemented blocks. The no-index-two hypothesis is essential for this lifting mechanism.

### Corollary 1, printed page 66: odd order

Every finite odd-order permutation group has a regular subset. The proof invokes solvability of finite odd-order groups and Theorems 1-2. The assertion includes intransitive actions; one works on their induced orbit actions and combines the subsets. Do not silently claim the original paper gives an elementary proof avoiding the odd-order theorem.

### Other precise claims

Proposition 1, page 66: if n>=2, m=2^n, and G is a Sylow-2 subgroup of S_m in its natural action, then G has no regular subset. Here |G|=2^(m-1)>binom(m,k) for every k. The PDF text extraction drops some exponent formatting; the visual check confirms m=2^n, not m=2n.

Corollary 2, pages 66-67: let F be a field, G finite, H<=G, L=core_G(H), and [G:L] odd. If the faithful right F[G]-module V is induced from a one-dimensional F[H]-module W, then some vector v has abelian stabilizer C_G(v). The proof uses a regular subset of the permuted one-dimensional summands to force C_G(v)<=L, and faithfulness plus diagonal action makes L abelian. This is not a general assertion that all solvable modules have regular vectors.

## 2. Inheritance and novelty audit

Read repository policy, START-HERE, CURRENT-FRONTIER, ACTIVE-GUARDRAILS, RESEARCH-PROTOCOL and CORE-PAPERS. The relevant current sources are:

- `05-knowledge/results/marked_path_tournament_series_20261021/session23-symmetry/README.md`, Sections 1, 6, 7: ordered path action, reversed converse, and original free-action boundary m=n-k<=2. The m>=3 correction is positive at the population level; this does not say every individual object has nontrivial stabilizer.
- `05-knowledge/results/continuation_width_symmetry_20261020/README.md`, Theorem 5 and Corollary 5C: exceptional B of size <=2delta, independent quiet-triangle factors, and the faithful disjoint-action factorization used here.
- `01-canon/theorems/THM-481-paley-gauge-generates-gleason.md`: the repo already has genuine Paley/Fano/coding bridges, including explicit classical attribution. This report does not treat a known code/design bridge as newly discovered.
- Current mistakes ledger: preserve actual actions and quotients; historical messages contain overstatements such as identifying the two Fano systems with two orientations of one set of seven triples. The finite computation here retains the two distinct systems.

Targeted repository searches for Gluck, Seress, the paper title, regular subsets, setwise-trivial stabilizers and distinguishing number found no prior direct treatment of the imported theorem or these exact transfer formulas in the maintained theorem/reference/result files. This is a scoped search finding, not a claim of global originality.

Primary literature checks:

- Meslem and Sopena, *On the Distinguishing Number of Cyclic Tournaments: Towards the Albertson-Collins Conjecture*, https://arxiv.org/abs/1608.04866 , explicitly derives the cyclic-tournament two-color result from Gluck and proves the named canonical coloring works for Paley tournaments.
- Lozano, *Distinguishing Tournaments with Small Label Classes*, https://arxiv.org/abs/1707.05549 and the author-institution published PDF https://upcommons.upc.edu/server/api/core/bitstreams/e3c8f7ff-1373-4877-8a23-2505f40f9174/content , explicitly states rho(T)<=floor(n/2) and constructs tournaments attaining this bound. Thus neither two colors nor the general minority bound is new. The published article is Acta Math. Univ. Comenianae 88(3) (2019), 923-928.

Targeted searches for distinguishing colorings together with anti-automorphism, converse or duality found no exact three-color statement in the maintained repo or the primary literature results inspected. This does not establish priority; the statement below is presented only as an elementary corollary of the cited deep input.

Closest proved mechanism: the actual residual/quotient stabilizer, not an ambient relabeling group. Canonical hostile: Paley7 when converse is adjoined. Corrected near miss: confusing colors with a uniformly small selected set. Least-used sidecar: the number and minimum size of regular subsets of the bounded exceptional group.

## 3. Binary rigidification of ordered partial paths

**Proposition (PROVED using Gluck).** Let T be a finite tournament, P=(v1,...,vk) an ordered directed path of distinct vertices, k>=1, and U=V(T)\{v1,...,vk}. Put H=Aut(T;v1,...,vk), with every path position individually fixed. There exists S subset U with H_S=1. One can choose |S|<=floor(|U|/2).

**Proof.** A tournament automorphism group has no involution: any transposed vertex pair would have its unique arc reversed. In a finite group, non-self-inverse elements pair with their inverses; thus absence of nonidentity involutions makes its order odd. H is a subgroup and is odd. Its action on U is faithful because it already fixes every marked vertex. Apply Gluck to this permutation action. Complementing within U preserves the stabilizer and gives the size bound. If U is empty the action is already trivial and S is empty. QED.

**Connection contract.** Source: an ordered partial-path object with actual stabilizer H. Target: the same object plus a named binary color on U. Map: append the characteristic function of S. Preserved: all original tournament arcs and all path positions. New predicate: trivial stabilizer of the augmented object. Lost if the sidecar is forgotten: precisely this guarantee of rigidity. This is not an isomorphism-equivariant choice of one S, not a new identity for uncolored path counts, and not an arbitrary representative chosen canonically from a symmetric object.

**Equivariant-selector obstruction.** Suppose a rule assigns a single subset S(T) and commutes with all tournament isomorphisms. For every g in Aut(T), equivariance gives g S(T)=S(T), so the setwise stabilizer is all of Aut(T). Such a rule cannot choose a regular subset when Aut(T) is nontrivial. The identical argument applies to a fixed marked object and its residual stabilizer. The lawful invariant object is the family of all regular colorings, or suitable orbits of that family, with its multiplicity retained. Merely asserting existence does not furnish a canonical symmetry-breaking choice.

If R(T,P) is the number of such regular subsets, the regular colored extensions over the fixed uncolored isomorphism class number R(T,P)/|H|. This quotient is integral because H acts freely on its regular subsets. Replacing it by 2^|U|/|H| would count nonregular colors incorrectly.

The same statement works for any individually marked scaffold. It works for a faithful odd permutation action on a quotient carrier, not for a nonfaithful action with a forgotten kernel. If the kernel is nontrivial, every subset stabilizer contains it and no coloring of that carrier can remove it.

## 4. A sharp three-color remedy for converse

**Proposition (PROVED using Gluck).** For every finite tournament T with at least two vertices, three named vertex colors suffice to kill every nonidentity automorphism and every anti-automorphism. For some tournaments two do not suffice.

Let E=Aut^+(T) union Aut^-(T), where the second set may be empty. When nonempty, E is a group with ordinary subgroup G of index two and sign homomorphism E -> C2. Take a binary regular subset S for G. The stabilizer E_S has trivial intersection with G, so the sign homomorphism injects E_S into C2. Therefore E_S is either trivial or consists of the identity and one involutory anti-automorphism a. In the latter case choose a moved vertex v and give v a new third color, leaving all other colors unchanged. The original two-color partition can be recovered by merging the third color back into v's old color, so every new stabilizer belongs to E_S; fixing v excludes a. QED.

The restriction |V(T)|>=2 ensures the sign of an isomorphism or anti-isomorphism is uniquely determined by its vertex permutation. At order zero or one, the formal signed action can have a color-invisible negative identity; no vertex coloring kills that signed kernel. If one instead quotients to the ordinary permutation image, that image is trivial, but this is a different target.

This proof actually applies to every faithful finite permutation group with an odd normal subgroup of index at most two. No statement is made for a general solvable group or an arbitrary Fano symmetry group.

For reversed-converse partial paths, let E preserve P in order in the positive channel and reverse P in the negative channel. Its ordinary subgroup is H above. If |U|>=2, the action of E on U is faithful: an anti-automorphism cannot fix two tournament vertices, and an ordinary kernel element fixes all vertices. The same argument therefore gives at most three colors **on U**. With zero or one omitted vertex this conclusion need not hold, because an anti-map can be invisible on U. For example, a transitive two-vertex full path has no omitted vertices but has an anti-map reversing the path. Retain the converse/position channel rather than applying Gluck to the wrong image. If k is odd and |U|=1, there is no compatible anti-map, since the middle marked vertex and the omitted vertex would both be fixed.

## 5. Exact fixed-slack regular-subset transfer

Use the inherited disjoint action, not merely an abstract group isomorphism:

    V(T) = B disjoint_union A1 ... disjoint_union At disjoint_union F,
    |Ai|=3, |F|=u,
    Aut(T)=H x C3^t,

where H acts faithfully on B, each C3 rotates its own quiet triangle, and F is pointwise fixed. Define

    R_H = #{S subset B : H_S=1},
    rho(H) = min{|S| : S subset B, H_S=1}.

Gluck gives R_H>0, including R_H=1 and rho(H)=0 for empty B with trivial H.

**Theorem (PROVED relative to the inherited decomposition).**

    R(T) = R_H * 6^t * 2^u,
    regular-colored isomorphism classes above [T]
          = (R_H/|H|) * 2^(t+u),
    rho(T) = rho(H)+t <= floor(|B|/2)+t <= delta+t.

**Proof.** A subset has trivial full stabilizer iff its restriction to B has trivial H-stabilizer and its restriction to each Ai is nonempty and proper. A directed triangle has six such subsets, of minimum size one. Fixed singleton colors are arbitrary. Independence of the disjoint factors gives the first formula and minimum. Divide by |H|3^t for free regular-subset orbits. Apply Gluck and complementing on B for the bound. QED.

A useful polynomial version records the selected-class size exactly. Let

    r_H(z)=sum_(S subset B, H_S=1) z^|S|.

Then

    r_T(z) = r_H(z) * (3z+3z^2)^t * (1+z)^u.

Complementation makes r_H reciprocal about |B|/2, but does not select a preferred color class. This is a bounded exceptional polynomial multiplied by the exact unlimited-word factors. In a canonical compiler the coefficients must be attached to its actual quotient stabilizer; using the full automorphism group of an uncolored core in place of the stabilizer of its retained word can change r_H.

**Boundary.** Arbitrarily long ordered joins of triangles have delta=0, Aut=C3^t, and rho=t. Thus two available colors are enough but the number of vertices assigned the minority color is unbounded at fixed slack. The result gives a faithful marked cover and exact fiber multiplicities, not a canonical unmarked encoding with those extra symmetries erased for free.

## 6. Fano/Paley/converse: three genuinely different actions

Use V=Z/7Z and Q={1,2,4}. Define Paley7 by i->j iff j-i is in Q. Let

    L+ = {b+Q : b in V},
    L- = {b-Q : b in V}.

Each is a Fano triple system: every pair lies in exactly one of its seven triples. They are disjoint systems. Their union is exactly the fourteen directed cyclic triples of Paley7. A cyclic triple has one orientation in this fixed tournament; the fourteen triples are not 'the same seven triples with both orientations.'

The exact actions are:

- Aut(Paley7) = {x -> a x+b : a in Q, b in F7}, order 21.
- Aut^pm(Paley7) = AGL(1,7), order 42. Nonresidue multipliers are the anti-maps, and they exchange L+ with L-.
- Aut(V,L+) has order 168 and is GL(3,2) on its seven nonzero vectors. An explicit identification of cyclic labels 0,...,6 with nonzero vector bitmasks is (1,2,4,3,6,7,5).

The ordinary order-21 group preserves both Fano systems. The order-42 group is not a subgroup of the automorphism group of the individually selected system L+: its anti-coset swaps the two systems. In particular one must not describe 21 < 42 < 168 as a chain of these specified actions.

### Exact distinguishing data

All 5,040 permutations and 128 subsets were checked; a separate implementation constructs the groups arithmetically and counts unlabeled set partitions before restoring named color labels.

- Order 21: 84 regular subsets, with counts by size 0,...,7 equal to (0,0,21,21,21,21,0,0). Its distinguishing number is 2, minimum selected-class size 2, and it has four regular binary-subset orbits.
- Order 42: no regular subset, distinguishing number 3, and 1,512 regular named three-colorings (36 regular orbits).
- Order 168 Fano action: no regular two- or three-coloring, distinguishing number 4, and 4,032 regular named four-colorings (24 regular orbits).

The ordinary binary-regular subsets all have full signed stabilizer of size exactly two, generated by a unique anti-involution. This follows without enumeration: they kill the ordinary subgroup, hence leave at most C2, while the order-42 group has no regular subsets. Thus ordinary binary rigidity at this exact degree is always accompanied by one surviving converse symmetry.

The no-regular-subset statement for order 42 is immediate from constant-cardinality orbit sizes: any subset orbit has size at most max_k binom(7,k)=35<42. The same argument rules out regular subsets for the degree-seven Fano group, but that group is not covered by Gluck's solvable classification.

For three colors in AGL(1,7), give two distinct vertices separate singleton colors and all other vertices the third. An affine map fixing two distinct points is the identity. Together with the binary obstruction, this proves the sharp value 3 without exhaustive computation.

### An elementary proof that the full Fano action requires four colors

Four suffice: uniquely color three linearly independent nonzero vectors and give all others a fourth color. Any color-preserving linear map fixes a basis.

For three colors, consider the nonempty color-class sizes. Two colors already fail by the orbit-size bound. With three nonempty classes the possibilities, after sorting, are (5,1,1), (4,2,1), (3,3,1), (3,2,2).

- (5,1,1): fixing two distinct nonzero vectors still leaves a nontrivial shear fixing their span pointwise.
- (4,2,1): let x be the singleton and {y,z} the pair. If x=y+z, a nontrivial shear fixes their plane pointwise. Otherwise x,y,z form a basis; interchange y,z while fixing x. Either map preserves all classes.
- (3,3,1): fix singleton x. The six remaining vectors form three pairs {v,v+x}. If a three-class A picks one from each pair, either its vectors sum to 0 or to x; cyclically permuting those three vectors extends to a linear map fixing x and preserving A and its complement. Otherwise its pair occupancy is (2,1,0). Fix x and the selected vector from the split pair; adding x to one basis vector in the fully selected pair gives a nontrivial shear preserving all classes.
- (3,2,2): if the triple is a line, the group fixing that line pointwise acts by the four translations of the remaining affine plane. The difference of either pair is the same nonzero translation, preserving both pairs. If the triple is an independent basis, its stabilizer is S3, fixes the sum of the basis vectors, and permutes their three pair-sums. In the two-two partition of these remaining four vectors, fix the pair containing the basis sum; the transposition of the two basis elements in its other vector preserves both pairs.

Thus no three-coloring is regular. This proof agrees with the two differently represented exact enumerations.

### The second order-168 trap

Gluck's degree-eight exception is the solvable affine-semilinear group AΓL(1,8), of order 8*7*3=168, acting on eight field elements. It is not the degree-seven full Fano group. Its lack of regular subsets also follows from max_k binom(8,k)=70<168. The independent audit constructs its 168 maps over F8 using x^3+x+1 and checks all 256 subsets. Equal orders do not imply isomorphic groups, equal degrees, or transferable stabilizer data. No Fano conclusion is imported from that degree-eight exception.

## 7. What this does and does not suggest for open problems

Distinguishing colors and acyclic tournament colors are different predicates. The former prohibit color-preserving symmetries; the latter require each color class to induce a transitive tournament. An acyclic multicoloring or its covering number is not a regular-subset count, and none of the palette counts here should be substituted into a Hamiltonian-path or dichromatic-coloring formula.

This is a genuine symmetry/quotient transfer to the marked-path and fixed-slack grammar lanes. It is not a proof mechanism for LRC(14), JC(2), or another flagship conjecture by itself. To transfer there, a concrete faithful action on a retained carrier and a target-preserving marking operation are needed. A seven-window carrier is not automatically the Fano action, and adding a color that is unavailable to the physical problem would be a new sidecar, not a proof of its original predicate.

Cheap next tests:

1. Attach r_H(z), not just |H|, to the compiler's existing actual-stabilizer states and verify which word/port refinements change it. The direct product formulas above give controls.
2. Keep three symmetry channels in Fano-inspired experiments: selected Fano incidence, oriented tournament, and converse-merged tournament. A claimed equivalence must specify how its action crosses these channels.
3. For a reversed-converse partial path, inspect the action kernel on omitted vertices before seeking a regular coloring. The zero-omission transitive path is the cheapest hostile.
4. For a desired canonical rule, test equivariance separately from existence. No nontrivial group can equivariantly choose one individually labeled element of a free orbit from an unchanged symmetric object. A regular subset orbit can be an invariant output, but choosing one subset requires breaking symmetry or retaining alternatives.

## 8. Reproduction and scope of checks

Run from this packet directory:

    python3 check_groups.py --output receipt.json
    python3 -O check_groups.py --output receipt-O.json
    cmp receipt.json receipt-O.json
    python3 audit_independent.py --output audit-receipt.json
    python3 -O audit_independent.py --output audit-receipt-O.json
    cmp audit-receipt.json audit-receipt-O.json

The first implementation scans every permutation to discover the Paley and Fano groups. The second uses affine arithmetic and GL(3,2) matrices, counts set partitions rather than named color words, and constructs the degree-eight affine-semilinear hostile. It additionally checks the disjoint-product transfer for three product actions and checks all 1,283 ordered nonempty directed-path objects in all labelled tournaments of orders 1 through 4 for faithfulness, odd residual stabilizer, existence of a regular omitted subset, and the minority-size bound.

Both scripts use explicit exceptions, not optimization-disabled assertions. Their receipts contain source hashes. The implementations are independent representations written by the same agent, not an independent-agent review. The general statements are proved above; the finite checks are positive/hostile controls, not induction substitutes.
