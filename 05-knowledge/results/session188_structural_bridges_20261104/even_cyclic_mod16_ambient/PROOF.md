# A special modulo-sixteen cyclic query survives every fixed ambient

**PROVED:** for arbitrary even-order children A,B,C, each of order at least four, the fixed core query a=A,c=A union B distinguishes their cyclic quotient from its opposite modulo sixteen in every fixed ambient tournament. The residue is eight. This is a special all-cover identity, not a general modulo-sixteen transfer theorem.

## 1. Exact profile domain and physical convention

The source is T=C3[A,B,C], with A->B->C->A. The baseline target U=C3^op[A,B,C] keeps every physical child interior. Child interiors on the target may instead be replaced by tournaments with identical complete cover profiles. The two actual ports have x->y. The fixed query is

    a=1_A, c=1_(A union B),

where a_v=[v->x], c_v=[v->y]. The source and target use the same physical names and bits.

For every even child of order at least four, `../even_cyclic_mod16/PROOF.md` proves the complete normalized modulo-sixteen upper domain

    q/H=(0,1,b,c,8,d,0,...),
    b in{6,14}, c in{4,12}, d in{0,8}.

Every H is odd. These eight formal residue profiles cover every actual child; no assertion that every formal vector is realizable is needed. All coordinates of index at least six vanish modulo sixteen.

## 2. Complete ordered-cover run expansion

For any h>=1, put the h paths of an ordered spanning cover in order. Within each path compress every maximal visit to a child into one run. The two ports remain single occurrences. For each child i, list its run slots in the resulting global path order. Filling those slots is exactly one ordered q_(k_i) cover of that child. Thus every skeleton contributes exactly one profile coordinate per child, and

    q_h(T+)-q_h(U+)
      = H(A)H(B)H(C) D_h(r_A,r_B,r_C) (mod16),          (1)

where r_i is its normalized profile. Path boundaries are retained: two consecutive appearances of the same child separated by a path cut are distinct runs and must not be merged. Every global path is nonempty.

Since q_k of a child vanishes modulo sixteen for k>=6, each child needs at most five runs in(1). Consequently there are at most seventeen symbols in a skeleton, including the two ports. For h<=5 the finite automaton below is exhaustive, not a bound inferred from sampled tournaments.

## 3. Exact finite cover-word certificate

A state records the three child visit counts in0..5, which ports have appeared, the last symbol, and the number of path cuts so far. Append an unused port or a child with fewer than five visits.

- Without a cut, the directed arc from the last symbol is required. Equal consecutive child symbols are forbidden because they would be one maximal run.
- With a cut, no arc is required; equal child symbols are allowed. The number of cuts increases by one.
- The first symbol starts the first path. Cuts are inserted only between two symbols, so no empty path is created.

Accept states in which both ports have appeared and all three child counts are positive. A state with h−1 cuts contributes to the ordered h-cover kernel. Sum over all possible final symbols and all lengths through seventeen. The resulting exact integer tensors count every skeleton once.

The source and opposite-quotient tensors for the fixed query are retained completely. Evaluating them on all512 ordered triples of the normalized upper-domain profiles proves

    (D_1,D_2,D_3,D_4,D_5)=(8,0,0,0,0) (mod16)          (2)

for every triple. The source also checks D_1 against the independently audited rational Hamilton-word functional, coefficient conventions and all512 exact integer values. A separate independently implemented cover-word automaton reproduces(2).

The odd scalar in(1) fixes eight modulo sixteen. For every actual augmented tournament,16 divides h! for h>=6, so every such higher ordered-cover count vanishes modulo sixteen. Thus(2) proves the entire physical augmented-profile identity

    (Delta q1,Delta q2,Delta q3,...)=(8,0,0,...) (mod16). (3)

## 4. Arbitrary fixed ambient lift

Let R be any fixed tournament with a distinguished root r, replaced by T or U. For every old exterior vertex z, give both ports the original root incidence: z->r means z->x,z->y, and r->z means x->z,y->z. These bits are the same on both sides. The union of the core with the two ports is a common module, and its collapsed ambient is R.

Maximal-module-run substitution gives fixed nonnegative integer coefficients c_h with c1=H(R), such that the actual full-query contrast is

    Delta F=sum_(h>=1)c_h Delta q_h
           =8H(R)=8 (mod16).                            (4)

The last equality uses the oddness of every tournament Hamilton count. The exterior's order, orientation and strong connectivity are unrestricted. No bound on the sizes, profiles or orientations of the even children is used beyond the stated order conditions.

The known failure of a GENERAL modulo-sixteen ambient transfer is not contradicted: a general module pair need not satisfy the all-coordinate vanishing(3). This fixed cyclic query does, by a separate complete certificate.

## 5. Scope

This proves baseline cyclic quotient-converse separation in arbitrary fixed ambient for precisely the all-even-order-at-least-four child domain. The companion modulo-eight and modulo-four results handle the remaining child domains. The generic negative-dual monotonicity lemma can then extend baseline separation to original-module root-containing families; that synthesis is stated separately.

The query remains child-uniform, so this is not a full named-profile reconstruction theorem for arbitrary quotients. Every recorded physical count is the complete Hamilton count of the actual two-port extension, including all cubic and quartic response contributions. No exact collision is asserted.

## 6. Controls

The verifier stores the exact integer source/target cover tensors and every512 formal profile evaluation. It computes complete first-five ordered-cover coordinates on actual port cores by an independent physical subset/partition method, and literal Hamilton counts in fixed strong and ordered exterior contexts. Larger unequal child profiles and a profile-increasing original-module example are included. Full named cubic/quartic data remain preserved in the separately pinned local `even_cyclic_mod16` packet.
