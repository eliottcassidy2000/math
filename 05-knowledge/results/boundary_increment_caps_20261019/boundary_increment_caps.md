# Exterior degree caps resolve the second mixed (3,3,1) boundary obstruction

Date: 2026-10-19. Reviewed fixed-core result. This proves a fixed-core, fixed-projected-order exclusion and a general boundary-increment inequality. It does not prove the entire canonical (3,3,1) exclusion or the seven-thin conjecture.

## 1. Setting and boundary counts

Let P be a finite poset of maximum incomparability degree at most k. Let x,y be a critical incomparable root pair, and Q its full core: the roots and all their incomparable neighbors. Every vertex outside Q is either in A below both roots or B above both roots, and A<B. No q in Q can be below any member of A, or above any member of B. Thus every A/core relation is a<q or incomparability; every core/B relation is q<b or incomparability.

Write s_A(q)=|{a in A:a incomparable with q}| and define s_B(q) dually. Necessarily

    s_A(q)+s_B(q) <= k-d_Q(q).

For a fixed core word u which can appear before both roots, L(u) is the number of interleavings of A and u respecting P, divided by e(A). Define U_empty=1 and U_(wq)=L(wq)-L(w). Appending q at the end of an interleaving injects the smaller set into the larger. The extra interleavings counted by U_(wq) are exactly those with at least one A letter after the last core letter q. Define R(v), V_empty=1 and V_(qv)=R(qv)-R(v) dually on suffixes after both roots. The increments are nonnegative.

If a core extension has boundary words pre and suf before its first and after its last root, its number of ambient lifts equals

    e(A)e(B)L(pre)R(suf).

Expanding both telescoping sums gives generator coefficients lambda_(u,v)=U_u V_v. Each generator counts core extensions beginning with u and ending with v, with both roots between the words. Only boundary pairs supported by at least one such extension contribute. Unsupported pairs need not have zero coefficients; we simply do not use them.

## 2. Quantitative ancestor cap lemma

For every nonempty prefix word wq,

    0 <= U_(wq) <= s_A(q) U_w.

This includes w empty, with U_empty=1. Dually,

    0 <= V_(qv) <= s_B(q) V_v.

Proof. Fix a linear extension a_1,...,a_m of A. Delete q from an extra interleaving of wq. When w is nonempty, at least one A letter remains after the last letter of w, so the image is an extra interleaving counted by U_w. When w is empty, the image is the fixed A extension. For any such image, a legal reinsertion must place q after every letter of w, after every A predecessor of q, and before the final A letter. Let a_j be the last A predecessor of q in this A extension, taking j=0 if none exists. Every a_(j+1),...,a_m is incomparable with q, since q cannot lie below A. Hence m-j <= s_A(q). There are at most m-j possible A gaps after a_j but before a_m. Imposing the requirement to follow w only removes choices, and within any A gap there is at most one position after all of w. Thus every image has at most s_A(q) preimages. Sum over all A extensions and divide by e(A). If A is empty, every nonempty increment is zero and the inequality still holds. Order duality proves the suffix case.

Consequences include hereditary support: a positive increment requires every ancestor increment to be positive. Iteration gives

    U_(q_1...q_r) <= product_i s_A(q_i).

A separate direct gap argument gives L(u) <= product_(q in u)(s_A(q)+1): for each fixed A extension, each q can occupy at most s_A(q)+1 A gaps, and a gap assignment uniquely determines the merged word because the core order is fixed. This latter bound is not needed below.

For every supported pair (u,v), multiplication of the ancestor cap by the opposite nonnegative increment gives

    lambda_(wq,v) - (k-d_Q(q))lambda_(w,v) <= 0,
    lambda_(u,qv) - (k-d_Q(q))lambda_(u,v) <= 0.

The parent of every supported child is itself supported. These are valid linear constraints on the supported generator coefficients alone. We relax s_A and s_B independently to k-d_Q; we do not need the stronger joint-slack restriction or any allocation branching.

## 3. Fixed second obstruction and exact certificate

Set k=7, x=0,y=1, and use predecessor masks

    (0,28,0,0,12,12,13,125,125)

with projected strict majority order

    (2,3,0,4,5,6,1,7,8).

The core has 408 linear extensions, 133 supported boundary generators, and degrees

    (5,5,2,2,3,4,3,2,2).

Thus common vertex5 has exterior slack3; there is no degree6 vertex, so the one-slack split does not apply. No equality of the lower twins2,3 or upper twins7,8 in the ambient poset is assumed.

For each generator g let E_g count its core extensions and N_g(a<b) count the extensions satisfying a<b. Write M_g(a,b)=3N_g(a<b)-2E_g. Define the weighted comparison row

    W_g = 4M_g(0,4) + 6M_g(6,1) + 5M_g(4,5)
          + 9M_g(5,6) + M_g(7,8).

All five comparisons are incomparable and forward in the projected order. Each cap row C has +1 at its child, minus the stated degree-slack bound at its parent, and zero elsewhere. The fixed JSON certificate supplies 18 nonnegative integer cap multipliers t_C. Exact enumeration verifies

    -540 <= W_g - sum_C t_C C_g <= -1

for every one of the 133 generators.

For actual ambient coefficients lambda>=0, every cap row satisfies C dot lambda <=0. Therefore

    W dot lambda <= sum_C t_C(C dot lambda) - sum_g lambda_g < 0.

The strict inequality holds because lambda_(empty,empty)=1. But if every listed forward comparison had probability strictly above2/3 in P, each individual comparison margin after summing actual lift counts would be positive, forcing W dot lambda>0. Contradiction. Thus this full core and projected order cannot occur in an actual seven-thin ambient poset with all forward incomparable probabilities strictly above2/3.

## 4. Why this closes a real relaxation gap

The second obstacle admits a nonnegative factorized formal array with U_empty=V_empty=1, U_(2)=V_(8)=1, U_(2,3,4,5)=V_(5,6,7,8)=6 and all other nonempty increments zero. Its formal weighted core distribution has total mass948 and every forward comparison strictly exceeds2/3. This does not arise from any exterior: for instance U_(2,3,4,5)>0 but its immediate parent U_(2,3,4)=0, contradicting the ancestor cap. More importantly, the certificate above excludes every array satisfying the necessary degree caps, not merely that one witness.

This establishes why free factorization alone is too weak and why hereditary support alone would also be an incomplete fix: strict feasible formal arrays can be perturbed by small positive entries, while the degree caps impose quantitative ratios.

## 5. Replay, dependencies, and scope

verify_caps.py is standard-library only. It enumerates all9! permutations independently, retains the 408 core extensions, reconstructs the133 generator columns and checks all18 cap bounds and the five comparison rows. It uses explicit exceptions rather than assert statements; normal and optimized-mode verification outputs match. Its only input is fixed_certificate.json. Numerical discovery files are not needed for verification.

The proof repeats the root-core factorization and supplies the ancestor-cap argument in full. It needs neither a completed morphology census nor the six-thin theorem, canonical tie-breaking, rank-one polynomial certificates, or any external computation. Maximum-span or canonical assumptions can imply the stated root setup but are not used beyond it. No ambient-poset counterexample, whole-morphology exclusion, or novelty claim is made.
