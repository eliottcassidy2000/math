# Every cyclic quotient converse is visible in arbitrary fixed ambient

## 1. Unified statement with physical names

Let A,B,C be arbitrary nonempty tournaments on disjoint fixed name sets. Put

    T=C3[A,B,C], with A->B->C->A,
    U=C3^op[A,B,C],

retaining every child interior. More generally, the target interiors may be replaced on their same name sets by tournaments with identical corresponding complete ordered-path-cover profiles.

Let R be any fixed rooted tournament. Replace its root r by T or U, keeping every exterior name and arc fixed. Attach actual ports x->y. A block-uniform core query is specified by three-bit masks(a,c), with a_v=[v->x], c_v=[v->y] constant on each whole child. For every exterior vertex z give BOTH ports the original root incidence: z->r gives z->x,z->y, while r->z gives x->z,y->z. The bits and physical names are identical in the two comparisons.

**Theorem 1.** Some such block-uniform query distinguishes the two full physical ambient responses, for every choice of A,B,C and R. No child-size, parity, common-profile, transitivity, or strong-connectivity restriction remains.

The proof has three separate residue branches. Their exact hypotheses and the reason their higher cover coordinates vanish are retained below. A nonzero local integer gap alone is never used to infer an ambient gap.

## 2. Singleton-containing children: modulo four

The complete five-type transfer in `../five_profile_types/PROOF.md` assigns every nonempty actual child to S=TT1, P=TT2, E=TT4, O=TT3 or C=directed C3, with the complete block-uniform response multiplied by a query-independent odd unit modulo eight. The C representative is scaled by H(actual)/3; normalized formal entries are not passed off as literal tournaments.

Exactly61 of the125 ordered type triples contain S. The retained complete canonical table proves that EVERY one of these61 triples has a query whose contrast is2 modulo4. The verifier extracts the first such query for every triple and retains the full mod4 row, so completeness is explicit. These are all actual canonical graphs, and their complete banks were independently counted in the five-type audit.

Restoring the selected block bits to the actual children and multiplying by the odd unit keeps the Hamilton contrast2 modulo4. For any two same-order nonempty tournaments, all ordered-cover coordinates q_h with h>=2 have universal residues modulo4:

    q2=2 for order>=2,
    q3=2 at odd order>=3 and0 at even order,
    qh=0 for h>=4.

These follow from q_h/h!=S(N,h) modulo2 and factorial divisibility, as proved in `../parity_profile_transfer/PROOF.md`. Thus the compared augmented cores have complete difference

    (Delta q1,Delta q2,...)=(2,0,...) modulo4.             (1)

No exposed-laminar assumption is used. The singleton branch is a complete type classification and a same-order cover-residue argument.

## 3. Nontrivial children outside the all-even-large case: modulo eight

Assume every child has order at least two, and at least one child has odd order or is a pair. This is precisely the63 nontrivial type triples other than(E,E,E).

The theorem in `../cyclic_odd_profile_lift/PROOF.md` supplies a block query with complete augmented profile difference

    (Delta q1,Delta q2,...)=(4,0,...) modulo8.             (2)

For even total core order the higher-coordinate conclusion follows from the even profile formulas. For odd total order, a separate exhaustive32-type certificate selects a query with Delta H=4 and Delta q2=0; the odd-order formula then gives Delta q3=0. Every higher coordinate vanishes. The proof retains the unsafe odd-total query whose residue disappears in an ambient and the stronger replacement query; it does not ignore the odd-order correction.

## 4. Three even children of order at least four: modulo sixteen

The only remaining case is(E,E,E). Each actual child profile, normalized by its odd H, lies in the complete upper domain

    (0,1,b,c,8,d,0,...), b in{6,14}, c in{4,12}, d in{0,8}.

The local packet `../even_cyclic_mod16/PROOF.md` proves the fixed query a=A,c=A union B has Hamilton contrast8 modulo16 for every one of the512 formal profile triples. Its separate companion `../even_cyclic_mod16_ambient/PROOF.md` tracks all ordered path cuts and proves the stronger entire-profile identity

    (Delta q1,Delta q2,...)=(8,0,...) modulo16.            (3)

This is a special cyclic query identity. The known failure of general modulo-sixteen ambient transfer remains valid and is not contradicted.

## 5. The common ambient step

Under the prescribed physical exterior lift, the augmented core is a module on both sides, with the same contracted rooted ambient R. The exact module-run substitution formula gives

    Delta F=sum_(h>=1)c_h(R,r) Delta q_h,
    c_h>=0 and c1=H(R).

Every tournament H(R) is odd. Applying(1), (2) or(3), the actual ambient contrast is respectively2 modulo4,4 modulo8, or8 modulo16. It is nonzero in each case. This proves Theorem1 with all quantifiers as stated.

The selected core query depends only on its finite profile-type branch and witness table; the exterior bits depend only on the fixed root incidence. No large child is silently replaced by its canonical representative in the actual experiment.

## 6. Every root-containing original-module XOR family

Let the source core be the displayed T, with arbitrary nonempty children. Select any family of modules of this ORIGINAL core, with repetitions reduced modulo two, and include the whole root oddly. Toggle each edge according to the parity of selected modules containing its two ends. Proper selected modules may overlap or nest and need not remain modules during a sequential ordering of operations.

Every proper original module lies within one child, since the root quotient is a directed triangle. The target therefore has opposite root quotient and child interiors B_i satisfying

    q_h(B_i)>=q_h(A_i) for every i,h,

by hereditary original-module monotonicity in `../module_decomposition/PROOF.md`. The internal converse from the selected ancestor does not affect q.

Choose the baseline separator z from Theorem1 and pair it with its CORE converse-complement Jz=(1-c,1-a). Keep the original exterior incidences for both ports in both queries. The exact all-cover duality of `../monotone_dual_queries/PROOF.md` makes their baseline ambient contrasts exact negatives. At least one is negative. At that member, replacing the baseline opposite-quotient children by the profile-increased target children cannot decrease the target Hamilton query value, because the children remain modules with uniform port incidences and nonnegative run coefficients.

**Theorem 2.** The fixed pair{z,Jz} separates the stated source from EVERY root-containing original-module XOR family, in every fixed ambient with all selected operations inside the core. For a given source and fixed ambient, one member chosen by its negative baseline sign works simultaneously for all family targets. That choice may depend on the ambient; a single member with the correct sign uniformly over all ambients is not asserted.

There is no scalar-query fallback, no common-profile requirement on the source, and no equality requirement on the target profiles. This is source-relative detection, not pairwise injectivity among arbitrary family outputs. Exterior changes are excluded from this root-family statement.

## 7. Sharpened highest-node residual

For a general original-module XOR comparison, the earlier `../highest_quotient_reduction/PROOF.md` shows that a possible full-response collision has scalar equality, hereditary cover-profile equality, and only prime quotient reversals at its highest changed strong nodes; a highest changed linear quotient is already ruled out.

A highest changed C3 quotient is now ruled out as well. Its source/target child profiles agree under scalar equality. Restore every other disjoint highest node using profile equality, and restore the changed child interiors inside the chosen node to their original profiles. The lifted ports see each disjoint outside module uniformly, so these restorations preserve the actual queried Hamilton count exactly. The remaining comparison is precisely a baseline cyclic quotient converse in a fixed ambient, contradicted by Theorem1.

Hence every highest changed quotient in any still-possible collision must be prime of order at least five. There is no prime strongly connected tournament of order four, also checked directly in the verifier. In particular, if every prime quotient of the original strong-module tree is a directed triangle, every original-module XOR target different from the source is distinguishable by the full named response. This conclusion concerns the stated original-module action, not arbitrary tournament reconstruction.

## 8. Reproduction and retained evidence

The verifier checks all61 singleton canonical rows, directly counts complete first-three augmented cover coordinates for their chosen queries, retains larger-child ambient controls, and enumerates all original child-module families in a singleton-containing cyclic source in whole-core and cyclic-exterior contexts. A multiple-highest-node control verifies exact restoration of another disjoint C3 change. The nontrivial branches bind their separately frozen proofs and complete certificates. Full named cubic/quartic data remain preserved in those component packets; no truncated observer replaces a physical query here.
