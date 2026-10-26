# Session68 — decorated composition and a complete ternary chirality stratum

**Status:** CITED canonical modular decomposition and independent C3 node choices; PROVED elementary classification/group and full-realization-stratum consequences below; FINITE-EXACT controls in `receipt.normal.json`. No global novelty claim, CFSG analogue, closed catalog of all prime tournaments, or unrestricted metagraph classification is asserted.

## 0. Evidence and inheritance

Read-only repository pin: `fc3768ac5fcd87022fde9e6d7e3f06c0e5bc22e8`, checkout `session67-repo`. Publisher-owned working changes were not used as canonical truth. The explicitly supplied Session67 correction `session67-incoming/INCOMING_INTEGRATION.md` is an additional, uncommitted incoming review input; its 12-vertex witness is independently recomputed here. Source integrity and source correctness are distinct.

**Concurrent prior integrated before freeze:** commit `d761e4cf0ff92264b043586a6160b274dae0b29c`, `05-knowledge/results/tournament_composition_classification_20261026/`. Its `CLASSIFICATION.md`, `fibre/RESULT.md`, and `arithmetic/RESULT.md` already establish the decorated classifier and split extension, the sharp 4p valuation threshold, the unique bad order12 fibre with four ordinary types, full path-cover invariance, and the Sylow/carry deficit laws. These overlaps are inherited/independently verified here, not claimed as new. The new continuation is the complete all-depth ternary realization stratum and its full stabilizer-polynomial spectrum.

Closest proved mechanisms:

- `05-knowledge/results/session56_rmq_lca_connections_20261025/tournament-hierarchy/REPORT.md`: canonical strong-module tree; exact leaf-LCA arc recovery; ordered SCC species; prime cycle-index substitution.
- `01-canon/theorems/THM-1960-tournaments-compose-from-regular-seeds-the-spectral-substitution-law.md`: modular-prime seeds differ from strongly connected SCC atoms.
- `05-knowledge/results/weighted_core_transport_20261020/core_count.md` and `README.md`: the canonical **pair-free transitive-chain core**, positive weights, weight stabilizers, and the insertion/splitting growth dichotomy. Pair-free does not mean modular-prime; gcd-primitive is a third notion.
- `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`: complete path-cover profiles compose, scalar H does not.
- `01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md`: iterated C3 wreath subgroups and ordered-path-cover divisibility already exist. Those facts alone are not new Session68 results.

Corrected near miss: THM-4803's global assertion that automorphism order is determined by labelled cyclic triples is false. Its local formula is valid, but changing a descendant's ordinary isomorphism type changes an ancestor's decoration stabilizer. The supplied Session67 correction proves the first failure is at order12, conditional on cited BCIZ structure. This note independently tests the order12 witness, but does not claim a new proof of the minimality cutoff.

Primary structural sources checked online on 2026-10-26:

- [Boussaïri–Chergui–Ille–Zaidi, 2018](https://arxiv.org/html/1805.04607v1), modular decomposition and common strong modules of C3 realizations; labelled node-choice theorem. This is a cited input, not rediscovered here.
- [McConnell–de Montgolfier, Linear-time modular decomposition of directed graphs](https://www.cs.colostate.edu/pubserv/pubs/McConnell-rmm-digraphDecomp.pdf), canonical quotient representation.

The board is **prime quotient / ordinary child type / stabilizer gauge / chain-core scale / marked-path runs / ternary necklace arithmetic**. Anchor: a rigorous classification carrier. Niche: automorphism jumps invisible to C3. Wildcard: a complete restricted family with a stabilizer-polynomial recursion. The actual new deliverable relative to the inherited repo route is the entire ternary chirality stratum, its automorphism spectrum, and its normalized marked-path fibre spectrum.

## 1. The exact finite recursive classification carrier

Let D(T) be the inclusion tree of nonempty strong modules. Its leaves are the original vertices; an internal node has either a transitive quotient (linear node) or a prime quotient of order at least3. A **reduced decorated tree** retains:

1. each linear node's ordered list of children, with no linear child of a linear node;
2. each prime node's actual oriented prime quotient, and the identification of its vertices with the child subtrees;
3. recursively, the ordinary oriented isomorphism type of every child.

At a prime node, equivalent presentations are exactly quotient isomorphisms carrying each child to an isomorphic child. At a linear node, only positionwise child isomorphisms are allowed. Equal children of a linear node do **not** acquire a symmetric-group permutation factor.

**Classification theorem (CITED structure + elementary consequence).** Evaluation of reduced decorated trees gives a bijection with finite nonempty tournament isomorphism classes, after the indicated local quotient identifications. Labelled evaluation is recovered by the LCA rule: the arc between two leaves is the arc between their two child ports at their least common ancestor.

**Proof.** Canonical strong modules supply a reduced decorated tree. Conversely, substitution into a prime quotient makes its nonempty blocks the maximal proper strong modules (Session56's proof). Substitution of strongly connected factors into a linear quotient makes them the SCCs, in their unique order. The no-adjacent-linear convention precisely flattens repeated ordered sums. Induction identifies the evaluated tree as the canonical tree. Any tournament isomorphism preserves strong modules and quotient arcs, hence induces exactly the stated tree isomorphism; the converse follows by evaluation. This proves both existence and uniqueness.

This is an effective finite recursion when the actual prime quotients through the requested size are supplied or enumerated. It is **not** a finite list of all prime families and not a closed formula for A000568. A prime tournament can itself be the root quotient of order n, so the difficult part is not removed.

The useful comparison with simple groups is recursive factorization plus action/extension data, not identity of the notions of prime. For example, the odd carousel R_9 is modular-prime but has automorphism group C_9, which is not simple; the companion arithmetic lane proves the carousel statement. R_9 and C3[C3,C3,C3] are both regular order9 tournaments with30 cyclic triples, but the latter is modular-composite and has automorphism order81 rather than9. Neither regularity nor the cyclic-triple count determines the composition hierarchy.

## 2. Full automorphism groups, including the splitting gauge

Consider a canonical node T=Q[T_i : i in I]. Put

    K = product_i Aut(T_i),
    G = {sigma in Aut(Q) : T_(sigma(i)) is isomorphic to T_i for all i}.

Then there is an exact sequence

    1 -> K -> Aut(T) -> G -> 1,                       (2.1)

and it splits. Consequently

    Aut(T) is isomorphic to K semidirect G,
    |Aut(T)| = |G| product_i |Aut(T_i)|.             (2.2)

**Proof and gauge.** Every automorphism permutes the intrinsic child blocks and induces a member of G; the kernel is K. For each ordinary child isomorphism type c choose a representative R_c and, for each child of that type, an isomorphism phi_i:R_c -> T_i. Define

    s(sigma)|_(T_i) = phi_(sigma(i)) phi_i^(-1).

It preserves all internal and interblock arcs, and s(sigma)s(tau)=s(sigma tau). This proves surjectivity and splitting. In the chosen gauges, G permutes equal-type factors of K. Replacing the phi_i by other such gauges conjugates these coordinated sections by an element of K. There is no canonical choice of phi_i; we do not claim that every abstract complement of K is kernel-conjugate to every other complement.

For a linear node G=1. For a prime node G is the **ordinary decorated stabilizer**, not Aut(Q) without decoration, not the stabilizer of sizes, and not the stabilizer of C3 types. This is precisely the missing coordinate in the 12-vertex example.

**Presentation-partition hostile.** For an arbitrary substitution partition the same formula computes only its block-preserving subgroup. Set

    T=C3[C3,C3,C3],
    Q=C3[C3,1,1],
    T=Q[1,1,1,C3,C3].

The second presentation is not canonical: its group preserving the five displayed blocks has order27, while full Aut(T) has order81. Outer rotations send a singleton presentation block into the interior of a size3 block. The independent vertex-level checker verifies this example.

## 3. Growth and the ancestor-stabilizer response

If a fixed **valid reduced tree** is retained while one subtree is changed, only its internal group and the decorated stabilizers along the ancestor route can change. The exact order response is the ratio of the new to old internal group order times the product of the new to old stabilizer orders at those ancestors. This ratio need not be monotone: decorations can split or merge isomorphism types.

Substituting arbitrary nonempty factors into a fixed prime quotient always preserves those top canonical blocks. Replacing a child of a linear node by a nonstrong factor instead requires SCC flattening; retaining the old node as if it remained canonical is invalid. For arbitrary vertex addition/deletion, recompute the canonical tree rather than assuming a local tree update. The Oct20 chain-core insertion/splitting theorem is a different, already proved representation of growth, not a modular-prime insertion theorem.

## 4. A complete family inside one labelled C3 fibre

Let A be `source -> C3`, and let B=A^op. The source A has scores (3,1,1,1); B has scores (0,2,2,2), so they are not isomorphic. Both have automorphism group C3 and the same single labelled cyclic triple.

Let m_d=3^d and I_d=(3^d-1)/2. Take the complete rooted ternary tree of height d, orient every internal quotient as the same C3, and replace its m_d terminal leaves by A or B. Fix all leaf blocks and labels throughout. Write T_w for the resulting order4m_d tournament, where w is a binary word of length m_d.

**Common-C3 theorem.** Every T_w has exactly the same labelled cyclic triples, with

    c_d = (8m_d^3 - 5m_d)/3.                         (4.1)

**Proof.** Reversing a module changes no cyclic triples: triples wholly inside retain cyclicity, those with two module vertices and one exterior vertex are transitive, and other triples are unchanged. Here only the terminal A/B modules are reversed. The count satisfies c_0=1 and c_d=3c_(d-1)+(4·3^(d-1))^3, which solves to (4.1).

All internal C3 substitution blocks are canonical maximal strong modules. Therefore T_w and T_v are isomorphic exactly when their binary-colored cyclic ternary trees are isomorphic; only cyclic rotations at internal nodes are allowed. Write W_0=1 and W_(d+1)=W_d^3 semidirect C3, in its natural action on the m_d leaves. Then

    Aut(T_w) is isomorphic to C3^(m_d) semidirect Stab_(W_d)(w).   (4.2)

The action of the stabilizer on the terminal factors is its actual leaf permutation. In particular,

    |Aut(T_w)| = 3^(m_d+e(w)),   0 <= e(w) <= I_d.     (4.3)

The group, not just its order, is obtained recursively: if the three ordinary colored child types are not all equal, use the direct product of their groups; if they are all equal, use their common group cubed semidirect C3. Two different colored types may still have isomorphic abstract automorphism groups; the group is an invariant, not the classification code.

### 4.1 The palette meets every ordinary realization class

Write C_d for this common labelled cyclic-triple hypergraph, and R(C_d) for **all** labelled tournament realizations of C_d. The palette has only2^m_d labelled tournaments. Nevertheless:

**Full ordinary-fibre theorem (with CITED node-choice input).** Every ordinary isomorphism class in R(C_d) has a representative T_w in the palette, and two such representatives are isomorphic exactly by the cyclic-tree criterion above. Therefore the polynomial in Section5 classifies all ordinary realization types of C_d, not merely a subclass of them.

**Proof.** BCIZ gives the common strong-module tree and independent quotient choices. At every upper C3 node the only quotient choices are C3 and its converse. A reflection of its three child slots changes the latter to our standard C3, carrying the whole child subtrees with their slots. Recursively normalize each child. At a terminal four-vertex node the linear quotient has children of sizes1 and3, so its two orders give exactly A or B as ordinary types; either orientation of the internal triangle is absorbed by a triangle reflection. All these normalizing permutations preserve C_d. Thus every realization has a palette representative. Conversely an ordinary tournament isomorphism must preserve its canonical child blocks and the oriented quotient; after normalization only C3 rotations, not arbitrary slot reflections, can identify two palette presentations. This proves the stated criterion and exhaustiveness.

For this particular hypergraph the full automorphism and realization counts are

    |Aut(C_d)| = 6^(m_d+I_d),
    |R(C_d)| = 2^(2m_d+I_d).                          (4.4)

At each terminal triangle, its single hyperedge permits S3, while the unique size1 and size3 children cannot be interchanged. At every upper node all three child hypergraphs are identical and every S3 permutation preserves the complete transversal triples. Canonicity of the shared strong-module tree allows no further automorphisms, proving the first formula. There are m_d+I_d prime C3 nodes and m_d linear nodes of arity2; independent node choices prove the second. An ordinary type with exponent e has exactly

    |Aut(C_d)| / |Aut(T_w)| = 2^(m_d+I_d) 3^(I_d-e)   (4.5)

labelled realizations. This quotient is valid because two realizations of the fixed C_d are ordinarily isomorphic precisely when related by an automorphism of C_d.

## 5. The stabilizer polynomial and full spectrum

Let P_d(z) be the sum of z^e over the distinct ordinary tournament types in this family, where e is the quotient-stabilizer exponent in (4.3). Then

    P_0(z)=2,
    P_(d+1)(z) = (P_d(z)^3 - P_d(z^3))/3 + z P_d(z^3).   (5.1)

For an ordered triple of child types, the exponent is the sum of the three child exponents. Nonconstant triples form C3-orbits of size3. A constant triple of child type with exponent e contributes instead exponent3e+1. Removing the constant triples, dividing by3, and restoring them with the additional rotation proves (5.1), including coefficient integrality.

Examples:

    P_1 = 2+2z,
    P_2 = 2+10z+8z^2+2z^3+2z^4,
    P_d(1) = 2,4,24,4624,32955830624,... .             (5.2)

Thus r_d=P_d(1) satisfies r_(d+1)=(r_d^3+2r_d)/3. This is the C3 necklace formula used recursively on child **types**, not an ordinary binary necklace of length3^d. In particular the d=2 family has24 types, not the60 ordinary binary necklaces of length9.

**Full spectrum theorem.** For every d, every coefficient of P_d from degree0 through degree I_d is a positive even integer; both endpoint coefficients are exactly2. Hence every order

    3^(m_d), 3^(m_d+1), ..., 3^(m_d+I_d)              (5.3)

occurs inside this one labelled C3 fibre, with the exact number of types given by P_d.

**Proof.** Complement every terminal color. This preserves stabilizer exponent and has no fixed isomorphism type because the number m_d of leaves is odd, so the two color counts cannot agree. Every nonempty exponent class therefore has at least2 types. Inductively suppose the exponent interval is0..I. Every integer0..3I is a sum a+b+c from that interval. If the three exponents are not all equal, any choice of types is a nonconstant triple. If they are equal, use two distinct types with that exponent, again obtaining a nonconstant triple. These give the entire interval0..3I. A constant triple at exponent I gives3I+1=I_(d+1). Formula (5.1) gives endpoint coefficient2 from endpoint coefficient2; the base is immediate.

The two maximal types are the all-A and all-B colorings. The two minimal types can be constructed recursively by

    X_0=A, Y_0=B,
    X_(d+1)=C3[X_d,X_d,Y_d],
    Y_(d+1)=C3[Y_d,Y_d,X_d].                           (5.4)

Their leaf-B counts are (3^d-1)/2 and (3^d+1)/2, so X_d and Y_d are distinct. Their child decorations break every quotient rotation. Therefore their full groups are elementary abelian C3^(m_d). The maximal group is the full iterated cyclic wreath group W_(d+1), of order3^(m_d+I_d) and exponent3^(d+1), whereas the minimal group has exponent3. This is a group-structure separation, not merely unequal orders.

The maximal/minimal automorphism-order ratio is

    3^I_d = 3^((n-4)/8),     n=4·3^d.                 (5.5)

Thus C3-blind symmetry variation grows exponentially along this explicitly classified size sequence. A single changed terminal A->B has stabilizer index m_d in W_d, because W_d is transitive on its leaves; its ratio to the uniform automorphism order is only3^d. Breaking every local rotation requires a distributed decoration, not one exceptional leaf.

Orbit-stabilizer supplies an exact independent checksum:

    3^I_d P_d(1/3) = 2^(m_d).                        (5.6)

The two rotation-rigid types together contain exactly2·3^I_d labelled colorings, a fraction (2/sqrt(3))(sqrt(3)/2)^m_d of the palette. The palette does not exhaust the whole **labelled** C3 fibre: independent node reversals yield additional labelled realizations. Section4.1 proves that it does meet every **ordinary isomorphism orbit**. Combining (4.4)–(4.5) with (5.6) gives the full-fibre mass identity

    sum_e [z^e]P_d(z) · 2^(m_d+I_d) 3^(I_d-e)
      = 2^(2m_d+I_d).                               (5.7)

### 5.1 The precise base-three carry interpretation

The concurrent `tournament_composition_classification_20261026/arithmetic/RESULT.md` defines delta_3(T)=v_3(n!)-v_3(|Aut T|) and proves its local carry-minus-stabilizer recursion. Applying that law, rather than claiming a new general carry theorem, gives a concrete interpretation of our full spectrum. Since n=4·3^d has base-three expansion11 followed by d zeroes,

    v_3(n!) = 2m_d-1,
    delta_3(T_w) = m_d-1-e = 2I_d-e.                 (5.8)

Thus the exact deficit spectrum is the interval I_d..2I_d, with multiplicities from P_d in reverse order. Each upper C3 node combines three children of size4·3^h and contributes two carry units. One unit is cancelled exactly when its three ordinary child types agree and retain the quotient rotation. All terminal A/B nodes have zero deficit. Hence delta_3 is **I_d unavoidable carry units plus I_d-e broken-rotation units**. The uniform family attains the global strong maximum 3-part at d=1, but not at d>=2: the concurrent strong-max theorem has deficit1 whenever the base-three digit sum is2, while our uniform deficit is I_d>1.

## 6. What transports to path and metagraph incidence

Every unordered path-cover coefficient pc_T(k) is identical across this palette. For a module reversal, split each cover path into maximal runs inside the module and reverse each such run. All exterior transitions remain directed because the module is uniformly seen from outside; applying the operation again gives the inverse. The same bijection works for ordered covers and Hamilton paths. The quotient walk-content kernel of THM-3121 is the recursive numerical transport law for the common profiles.

Let H_d be the common Hamilton-path count. The automorphism action on directed Hamilton paths is free: a fixed ordered vertex path fixes every vertex. Likewise the action on ordered spanning k-path covers is free. The largest group in the palette therefore gives

    3^(m_d+I_d) divides H_d,
    3^(m_d+I_d) divides k! pc_T(k) for every k.         (6.1)

The factorial is indispensable for unordered covers. For the ordinary, unmerged marked-Hamilton-path fibre,

    t_c(T_w) = H_d / 3^(m_d+e(w)).                   (6.2)

Consequently every value in this full geometric interval occurs, and the largest/smallest fibre ratio is3^I_d. At d=1, H_1=462915, so the two fibre values are5715 and17145 even though all67 labelled cyclic triples agree.

This is **not** a classification of the full metagraph or its port adjacency. A path-cover count forgets which actual vertices form each run, their endpoints, and how repeated block visits allocate internal vertices. A lossless marked-path lift of the decorated tree must keep the actual ordered run components and their occurrence assignments at every quotient node. That data reconstructs the full labelled path and hence its complement partner. Counting runs, keeping H, keeping a cycle index, or retaining only a node-level tree does not do so. The Oct20 allocation/shuffle hostiles remain decisive. Neither merged-converse factors nor loop/edge-degree conventions have been silently identified with (6.2).

## 7. Verification and limits

Run in this directory:

    python3 -B check_classification.py > receipt.normal.json
    python3 -B -O check_classification.py > receipt.optimized.json
    cmp receipt.normal.json receipt.optimized.json

The script has no external input and uses explicit exceptions, not assert. Universes:

1. All1,099 labelled tournaments of orders1..5: brute subset strong modules, decorated-tree classification against full-permutation canonicalization, and recursive automorphism orders against independent degree-pruned arc-preserving permutation backtracking.
2. All8 binary palettes at order12: exact labelled C3 sets, full unordered path-cover profiles via independent subset recurrence, Hamilton count462915, and complete vertex-level automorphism lists.
3. All128 independent node choices in the full labelled order12 C3 fibre: exact C3 equality, explicit hypergraph-preserving gauge into the eight-member palette, ordinary orbit sizes16,48,48,16, and independent complete hypergraph-automorphism backtracking giving1296.
4. All512 binary colorings of the nine-terminal depth2 tree: literal action of all81 W_2 permutations versus recursive codes and stabilizers.
5. Explicit cyclic-child-type enumeration through depth3 (4,624 types), independent of the polynomial operations; polynomial recursion and mass/integrality/support controls through depth5.
6. Noncanonical nine-vertex partition hostile, full group81 versus presentation-preserving subgroup27.

These finite checks validate implementations and examples; all-size conclusions come from the proofs. No enumeration through order36, through all prime tournaments, or through higher-depth full labelled C3 fibres is claimed.
