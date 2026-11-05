# Stronger insertion and compression bounds for the inverse bridge

Session 197. This short supplement strengthens, without changing, the earlier valid bounds in `PROOF.md`, `ADDENDUM.md`, and `FACTOR_BRIDGE.md`. **All-order exact inequalities and profile compression; complete odd prime unit-amplitude classification through 29. The full spectrum remains open.**

## 1. Insert every deletion path, not just endpoint paths

Let (T,r) have at least one out-neighbor, with the usual profile (h,e,R,d), slack f and U=R-d. For **every** Hamilton path P of T-r, insert r immediately before the first out-neighbor encountered along P. Every preceding vertex beats r, so insertion is legal. The result never ends at r, and deleting r recovers P. Hence

    h-e >= d,
    U=h-e+f >= d+f.                                 (1)

This is stronger than the endpoint-only injection h-e>=e used in the first packets. The root must have an out-neighbor; the universal-sink case has U=0 and is excluded from positive amplification.

For singleton-partner unit synthesis, beta=d-f, so (1) yields

    2d <= U+beta.                                   (2)

For the general two-piece determinant beta=dg-ft, the source dichotomy in `FACTOR_BRIDGE.md` still applies. At beta=-1 one has f>=d+1, and at beta=+1 either f=d-1 with g=t=1 or f>=d+1. Thus **the same bound (2) holds for every physical unit synthesis**, even with varying partners. In the nontrivial cancellation case f>=d+1, the stronger sign-independent bound is2d<=U-1.

The earlier bound3d<=2U+2beta-1 remains valid, but is weaker.

## 2. Compression controlled by U alone

Let S=T-r, and let C_j be its first strong component containing an out-neighbor of r. Delete every singleton component of S before C_j, keeping all other vertices and r. Call the resulting rooted tournament (Q,q).

All components before C_j beat r and every later S vertex. They are an ordered prefix of T as well as of S. Removing singleton components from this prefix changes none of its Hamilton counts or the root-starting cover response. Thus

    profile(Q,q)=profile(T,r).                       (3)

Choose x in C_j beaten by r. The exact SCC cover transport gives

    U >= R_x(S)
      = (product_(i<j) H(C_i)) R_x(C_j)
          (product_(i>j) P_Ci(2))
      >= 2^s,                                       (4)

where s is the total number of vertices after C_j, since P_C(2)>=2^|C|. Therefore s<=floor(log2 U).

The total number of vertices in nontrivial S components is at most d, as proved in the first packet. Before C_j all remaining components are nontrivial. C_j itself adds at most one singleton vertex, and all vertices after it total at most floor(log2 U). Adding r proves

    |Q| <= d+floor(log2 U)+2.                        (5)

This is exact profile-preserving compression at every positive U, with **no bounded-slack assumption and no rearrangement of the retained vertices**. It does not supersede the bounded-(d,f) dyadic-family theorem, which permits U to grow unboundedly and proves a different sparsity conclusion.

## 3. A simpler finite factorwise certificate

For a unit amplifier alpha=U L, beta=±1, (2) and (5) allow the source to be replaced by an actual piece satisfying

    |X| <= floor((U+beta)/2)+floor(log2 U)+2.          (6)

A factor pair is impossible if floor((U+beta)/2)<1.

For the partner, apply (5) to its converse W. The deletion count is t=H(Y-c), starting-cover response is L, and U_W=L-t. If U_W=0, the unit determinant forces L=g=t=1, so a singleton partner suffices. Otherwise, writing v=L-t>=1,

    |Y| <= t+floor(log2 v)+2 <= t+v+1=L+1,           (7)

using floor(log2 v)<=v-1. The normalization preserves g,t,L exactly. Hence every unit ray has a physical star certificate with the same factorization alpha=U L, the same sign, and source/partner bounds (6)-(7). The compiler remains the literal S196 graph construction.

This strengthens the earlier factorwise bounds and keeps the exact finite inverse procedure valid: enumerate factor pairs, then actual rooted pieces within these sharper bounds. Arithmetic solutions without those physical pieces remain insufficient.

## 4. Complete odd prime unit amplitudes through 29

For odd prime alpha=p, S196's complete star classification forces L=1 and U=p. Equation (2) gives d<=(p+beta)/2. The all-order d<=13 table in `ADDENDUM.md` therefore decides every negative case below, and every absent positive case through p=23. The positive p=29 case is exhibited directly.

**Theorem 6.** For odd primes p<=29, a physical binary cloning amplifier with ray p*2^m+beta, beta in {+1,-1}, exists exactly for

    (p,beta)=(3,-1),(5,-1),(11,+1),(29,+1).           (8)

The first three witnesses are already in the first packet. For the last, take the actual source

    X=[2,44,49,53,3,17], root5,
    (h,e,R,d)=(23,4,40,11), U=29, f=10,

and a singleton partner. The actual K-family has counts59,117,233,465,..., i.e.29*2^m+1.

To check completeness, the relevant primes are3,5,7,11,13,17,19,23,29. For negative offset the largest allowed deletion count at p=29 is14, hence odd d<=13. For an absent positive offset the largest p is23, whose bound gives d<=12, hence odd d<=11. The complete table applies, and none of those omitted pairs has its required U=p. Positive offset at p=29 is already a genuine construction, so no enumeration at d=15 is required for this existence classification.

These are amplitude statements, not forbidden Hamilton counts. For example p=23,beta=-1 would give a particular H45 clone presentation, and p=23,beta=+1 a particular H47 presentation. The theorem does not exclude H45 or H47 via other amplitudes or operations.

## 5. A shorter route to the H15 two-sign obstruction

The stronger inequality also removes the need for the extra f=6 test in the factorwise proof, although that earlier test is correct. For alpha=8,beta=-1, U divides8 and2d<=U-1. Only d=1 is possible when U=4; only d=1 or3 when U=8. A deletion-one source with U a power of two has f=0, impossible for a negative determinant. At d=3,U=8, (1) gives f<=5, while f>d and even force f=4. Its exact U bases are{9,15}, so U=8 is impossible.

Together with the excluded alpha7,beta=+1 this again proves that H15 has no unit binary-cloning root. `H15_CLASSIFICATION.md` independently proves the stronger one-step statement by classifying every H15 tournament and its complete cloning response.

## 6. Verification

`verify_sharpening.py` checks (1) for every positive-U root in every labelled tournament through order6, checks profile compression and (4)-(5) on a complete smaller universe and adversarial/larger controls, rebuilds all small unit-pair certificates with (6)-(7), verifies the prime-amplitude budget arithmetic against the exact earlier tables, and constructs the p=29 positive family. Normal and optimized runs must agree. No earlier frozen file is modified.
