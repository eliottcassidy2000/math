# Independent higher-adic audit

Status: **PASS: all-k class projection and all-valuation closed form, with independent exact controls at moduli3,9,27.** Nothing here derives a height-weighted law.

Independent verifier `audit_higher.py` groups normalized states by

    (d,g,D mod g), where g=gcd(d,U,q).

This differs from the worker's forward-component discovery. For every group at q=3,9,27 it verifies invariant membership, indegree2, and forward and reverse connectivity. It then projects arbitrary input mass by one O/Z step followed by uniform measure on the group, and calculates exact state probabilities through zero, one and two Js. Every residue probability matches the worker's HIGHER.txt.

In particular the two-J valuation-local factors are:
- q=9:47/54,173/144,11/8, for valuations0,1,and the residue0 bucket.
- q=27:47/54,173/144,1751/1296,461/324, for valuations0,1,2,and the residue0 bucket.

These are exact rational calculations, not finite-length numerical fits. Their fixed-modulus limiting interpretation uses the valid resolvent argument: initial d=1 class is aperiodic; all peripheral eigenvalues of a finite Markov transition are semisimple; only eigenvalue1 can contribute a pole of order3 in the three-block initial-row generating function. Other roots of unity contribute at most order2, yielding O_q(1/L) after normalization. No assertion uniform in q follows.

Normal and python3 -O verifier outputs agree. Artifacts:HIGHER_VERIFY.txt and HIGHER_OPT_VERIFY.txt. The q81 calculation has not been independently recomputed by this audit.

## Independent all-valuation derivation

The proposed closed form has now been independently derived, rather than fitted. This argument uses the all-k class-projection rule; its connectivity proof is audited separately below.

Use Haar-uniform3-adic variables to abbreviate compatible finite residue distributions. After the initial unit-d stationary block, x=D and y=D+U are independent Haar variables and E is0 or1 equiprobably. After J and one O/Z step, h=min(v3(x),v3(y)) and the frozen D coset modulo3^h is0 with weight3/4 or−1 with weight1/4. If h>=1 the latter branch has D unit, so the second J and projection give completely uniform F. Its total mass is(1/4)Pr(h>=1)=1/36. Put w_0=1 and w_h=3/4 for h>=1.

For the remaining mass there are two disjoint branches:
- A_h: v3(x)=h<=v3(y), mass(2/3)9^(−h)w_h. The projected D and U are independent Haar multiples of3^h. Their minimum valuation is j=h+t with probability(8/9)9^(−t). If t=0 the final F is uniform on multiples of3^h. If t>=1 it is3/4 uniform on multiples of3^j plus1/4 uniform on the valuation-h stratum, because the second frozen coset is−x in that branch and x has valuationh.
- B_h: v3(y)=h<v3(x), mass(2/9)9^(−h)w_h. The projected U has exact valuationh, so the next gcd has valuationh and F is uniform on multiples of3^h.

These statements include the x=0 or y=0 finite-modulus boundaries by taking compatible residue limits; the zero3-adic values have Haar mass zero. The final residue-zero buckets are recovered by summing all valuations at leastk, not by discarding finite-modulus zero states.

Consequently the density relative to uniform residue mass on the valuation-v stratum is exactly

    f(v)=1/36 +(22/27) sum_(h=0)^v w_h 3^(−h)
         +(4/9) sum_(h=0)^(v−1) w_h 3^(−h) sum_(t=1)^(v−h)3^(−t)
         +w_v 3^(−v)/36.

The first sum combines branchB and branchA,t=0; the second is the deeper-zero-coset contribution; the last is the nonzero coset of exact valuationv. Evaluating the geometric sums gives

    f(0)=47/54,
    f(v)=157/108−3^(−v)(v/6+85/144), v>=1.

Indeed S_v=sum_(h=0)^v w_h3^(−h)=11/8−(3/8)3^(−v) for v>=1, and sum_(h=0)^(v−1)w_h=(3v+1)/4. Substitution yields the displayed constant and correction exactly.

The total probability normalizes:

    sum_(v>=0) [2/3^(v+1)] f(v)=1.

The zero-residue factor at modulus3^k is forced by the tail sum:

    T_k=3^k sum_(v>=k) [2/3^(v+1)]f(v)
       =157/108−3^(−k)(k/8+11/24).

It gives34/27,11/8,461/324,2803/1944 at k=1,2,3,4, matching the independently checked first three cases and the worker's fourth. The limiting factor is157/108 and its ratio to the unit factor is157/94. These are fixed-length-ensemble statements only.

## General class connectivity: independent short proof

For normalized states E in{0,d}, fix g=gcd(d,U,3^k) and D modulo g. Divide all changing coordinates by g after subtracting the frozen D coset. It suffices to treat gcd(d,U,3^k)=1. On the E=0 slice, the allowed words generate permutations

    S(D,U)=(D,2U),
    A(D,U)=(D+U,U+h), h=−d/2,

where A is ZO followed by the inverse of S². All inverses are realizable by positive powers of permutations on this finite slice. Since2 is primitive modulo3^k, R=S^(3^(k−1)) sends U to−U. Then(R composed with A)² sends(D,U) to(D−h,U), giving every D-translation in the ideal(d).

If d is a unit, all D shifts are available; combining A with a D correction gives arbitrary U shifts. If d is nonunit, U must be a unit. Powers of S reach every unit U. To adjust D while fixing a chosen U=u, apply A^m and then a power of S restoring u+mh exactly to u. This power exists because both are units and2 generates the unit group. The D increment is mu+h*m*(m−1)/2, congruent to mu modulo the ideal(d). Choose m to set D modulo that ideal, then use the pure D translations to correct the remaining difference. If d=0, A already fixes U and translates D by the unit U; handle this case directly. Applying Z reaches the E=d slice.

Thus every proposed class is strongly connected. Each has indegree2, so its stationary law is uniform. The initial d=1 class has the O self-loop at E=U=0 and is aperiodic. For arbitrary initial E, the first O/Z step reaches normalized E and freezes the D coset as D or D−E modulo g. This proves the all-k projection and validates the fixed-q resolvent extension, including all first-step boundaries.

Primitive-root input:2 generates units modulo3^k because its order modulo3 is2 and v3(2^(2*3^r)−1)=r+1 by induction (or LTE). Hence its order modulo3^k is2*3^(k−1). No unproved number-theoretic conjecture enters.

## Final comparison with ALL_K.md

Reviewed the worker's completed `../threeadic/ALL_K.md`: **PASS**. Its alternative density-integral derivation agrees with the disjoint-branch derivation above. In particular, T_j(a) correctly distinguishes a<=j from a>j; in the latter case only v3(U)=j and v3(D)>=j survive the density. The derived masses m_j=(24j+41)/(54*9^j) and Z_j=(11/18)9^(−j) are correct. Its nonzero-coset tail is1/(48*3^v) for v>=1, while v=0 has the separate1/36 contribution. All sums converge absolutely.

The Haar formulation is justified by compatible cylinder probabilities, not an interchange of k and L limits. Finite Markov Cesaro projections commute with reduction modulo lower powers of3. Equivalently, construct the law directly by independent Haar d,U and normalized Haar measure on each finite-index D coset; its reduction gives exactly the finite-q projector, including the bucket where both d,U vanish modulo q. The asymptotic statement always fixes k before L tends to infinity.
