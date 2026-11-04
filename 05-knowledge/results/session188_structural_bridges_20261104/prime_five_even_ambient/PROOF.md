# Prime five-name even-child converses are visible in every fixed ambient

## 1. Fixed physical setting

Let Q be a prime five-name quotient and let A_0,...,A_4 be arbitrary tournaments of even orders at least four on fixed disjoint named sets. Put T=Q[A_i] and U=Q^op[A_i]. Target child interiors may also be replaced while preserving their complete source profiles. Let R be any fixed rooted ambient tournament whose root is replaced by T or U.

Use the quotient block masks selected below inside the core. At every exterior vertex v, both physical ports imitate the original root: v->root means v->x,v->y, and root->v means x->v,y->v. Keep x->y. Thus the entire port-augmented core is a common module in both actual extensions. All exterior vertices and arcs are fixed.

**Theorem.** At most two fixed block-uniform queries, lifted by these incidences, distinguish the source and opposite quotient in every fixed ambient. Child orders and profiles need not agree with one another. The theorem uses the independently certified local packet `../prime_five_even_profiles/PROOF.md` and the separate all-cover certificate below; it does not use an unrestricted modulo-sixteen module-transfer claim.

## 2. Two quotient types already survive modulo eight

Use the complete prime representatives and labeling conventions in the local packet. For Q10=(20,1,2,7,14), the fixed query(0,2) has normalized contrast12+8(epsilon0+epsilon2+epsilon4), hence four modulo eight. For Q41=(18,8,3,5,14), query(0,1) has normalized contrast6+8(epsilon0+epsilon1+epsilon2), hence six modulo eight. Multiplication by the product of the actual odd child Hamilton counts leaves both contrasts nonzero modulo eight.

The port-augmented core has even order. The previously proved even-order modulo-eight transfer in `../mod8_transfer/PROOF.md` therefore multiplies each core contrast by an odd unit. Both fixed queries remain separators in every fixed ambient. This step does not require the new modulo-sixteen certificate.

## 3. The regular quotient: all necessary cover coordinates

Let Q76=(24,17,3,6,12). Normalize each actual child modulo sixteen as

    r_i=(0,1,14+8e_i,4+8f_i,8,8g_i,0,...),
    e_i,f_i,g_i in {0,1}.

This is the complete upper domain from the local packet and `../even_cyclic_mod16/PROOF.md`; its fifteen bits need not all be realized. Let Lambda be the product of the odd actual child Hamilton counts. Every ordered-cover count of the two-port extension is multilinear in the child profiles, so normalization multiplies every coordinate by the same Lambda.

For z0=(a,c)=(0,1) and z1=(2,18), the separate explicit-cut run calculation proves, for all2^15 formal bit assignments,

    (Delta q1,Delta q2,Delta q3)(z0)
         = (8(e2+e3),0,0) modulo16,
    (Delta q1,Delta q2,Delta q3)(z1)
         = (8+8(e2+e3),0,0) modulo16.                    (1)

The fifteen first-order coefficient vectors and both baseline vectors are retained for each orientation. All products of two perturbations vanish modulo sixteen. Hence these affine coefficient identities prove the entire domain, in addition to the exhaustive state replay. In particular f_i and g_i have not been discarded by a Hamilton-only involution. The ordered-cover calculation retains them and verifies their cancellation in the stated differences.

At exactly one of these queries, the actual first-three-coordinate contrast is (8,0,0): multiplying by the odd Lambda preserves eight and zero. The two augmented actual cores have the same order N. For k=4,5, forgetting slot order gives

    q_k/k! = S(N,k) modulo2,

because every induced nonempty tournament has odd Hamilton count. Since4! and5! are both eight modulo sixteen, q4 and q5 depend only on N modulo sixteen, so their source-target differences are zero. For every k>=6, sixteen divides k!, so the difference is zero too. Therefore the selected query has the COMPLETE augmented-profile contrast

    (8,0,0,0,...) modulo16.                            (2)

No tail coordinate has been omitted.

## 4. The ambient lift and source-relative families

Exact common-module substitution expresses the ambient Hamilton contrast as sum_k b_k Delta q_k, where b_k are nonnegative integers determined solely by the fixed exterior and b_1=H(R), which is odd. By(2), the selected regular-quotient query has ambient contrast8H(R)=8 modulo16. Together with Section2 and the complete prime classification, this proves the theorem.

The same result yields a source-relative original-module-family consequence. Restrict a simultaneous XOR family of original modules to the core and include the core itself after reducing repeated selections modulo two. Since Q is prime, every proper original module lies within one child. Hereditary original-module monotonicity from `../module_decomposition/PROOF.md` makes each target child cover profile coordinatewise no smaller, while the quotient reverses.

For a baseline separating query z, pair its core masks with Jz=(1-c,1-a), keeping the same lifted exterior incidences. The complete profile duality proved in `../monotone_dual_queries/PROOF.md` gives opposite baseline ambient contrasts. Choose the member whose source-minus-baseline contrast is negative. Uniformity on every child and nonnegative run substitution imply that increasing target child profiles cannot cancel that negative difference. Thus the baseline battery together with its two duals, at most four fixed physical queries, detects every such root-containing original-module family. For a fixed source and ambient, one sign-selected baseline query works for all these families.

This is source-relative detection with a fixed exterior, not pairwise injectivity of arbitrary orbit outputs. Operations involving exterior vertices or changing the ambient quotient are outside scope.

## 5. Explicit ordered-cover calculation and reproduction

A directed ordered h-cover is encoded by a global run word together with h-1 cuts between consecutive symbols. A transition without a cut must follow a quotient/port arc and cannot repeat a child symbol. A cut permits any transition, including repeated child symbols. Each path is nonempty, each port appears once, and each child's runs across all path slots form one ordered cover. Thus this calculation retains every slot placement, both port orders, all cuts, and all full cubic/quartic interactions.

The dynamic program keeps cut counts0,1,2 for the first three coordinates. It separately retains all fifteen perturbation channels; it does not reuse the Hamilton-only q3/q5 vanishing argument. Baseline terms whose total child-weight valuation is at least four vanish modulo sixteen exactly. Every first-order term has all other children at one run. All2^15 profile states are checked; their two profile vectors are((0,0,0),(8,0,0)) or the reverse, each in16,384 states.

Two direct physical subset-DP counts on24 names independently verify the nonzero regular battery query in a cyclic ambient with two exterior names. The physical ports use the required root incidences. The counts are exact residues modulo sixteen, not full-integer Hamilton counts. A separate independent audit must verify the all-cover identities and physical lift before publication.
