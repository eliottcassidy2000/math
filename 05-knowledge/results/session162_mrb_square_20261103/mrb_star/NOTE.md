# MRB filters on actual tournament star cubes

Status: **proved conditional on the inherited odd-cycle packing identity**, with exact finite controls. This is a generative connection, not a formula for the MRB constant or a proof of the Hamilton-path spectrum conjecture. No novelty claim against the literature.

## The all-order strong tournament / weighted-cut construction

Start with the transitive tournament X on 0,1,...,n−1, n≥2. Independently reverse 0→i iff x_i=1, for 1≤i<n. Replace the node X in a directed triangle by X to get G(X), so X→a→b→X. G(X) is strong. Define

    c_1=1,  c_d=2^(d−2) for d≥2,
    q_d = c_d + 2^(n−d−2) [3^(d−1)+1],  1≤d≤n−2.

The exponent n−d−2 is nonnegative. The number of Hamilton paths is exactly

    H(G(X)) = 1 + 2^n + 2 sum_{1≤i<j<n} q_(j−i) (1−x_i)x_j.      (1)

This realizes every value of the indicated weighted directed cut as the exact count of an explicit strong tournament. It is a *source-level* construction, not an assertion that its range covers all required odd integers. The full Boolean word matters: two words with the same cut weight need not have equal responses to future flips.

Proof. Every directed cycle in X uses 0. It is uniquely 0→i→u_1→...→u_s→j→0, where 1≤i<j<n, the u's are an increasing subset of the strict gap, x_i=0 and x_j=1. The cycle is odd iff s is even. All directed cycles meet at 0, so an odd-cycle packing contains at most one cycle. The inherited odd-cycle packing identity used in [S161 FAN](https://github.com/eliottcassidy2000/math/blob/ab860360b8dcae14454ea150a849ac8616233039/05-knowledge/results/session161_source_matching_20261103/matched_sources/FAN.md) and [profile formula](https://github.com/eliottcassidy2000/math/blob/ab860360b8dcae14454ea150a849ac8616233039/05-knowledge/results/apex_wedge_constructor_20261103/PROFILE.md) give F_r(X)=r^n+2Σ_C r^(n−|C|+1) and H(G(X))=F_1(X)+F_2(X). For a fixed gap d=j−i, the contribution to half the excess equals

    q_d = Σ_{s even} binom(d−1,s)(1+2^(n−s−2)).

The first summand is c_d; the second is 2^(n−d−2)[3^(d−1)+1] by the even binomial identity. This proves (1).

Direct subset-endpoint DP checks (1) for every star word for n=2,...,7, 126 cases. The controls do not replace the all-order proof.

## An MRB finite-difference projector

The actual MRB source profile is f(m)=exp(log(m)/m)−1, m≥1, and C=Σ_{m≥1}(−1)^m f(m), not the divergent ordinary sum with unshifted m^(1/m) terms. The previous [MRB research](https://github.com/eliottcassidy2000/math/blob/ab860360b8dcae14454ea150a849ac8616233039/05-knowledge/results/session144_jacobian_eth_mrb_20261103/mrb/RESEARCH.md) treats zero-mean periodic filters and Euler finite differences while also showing that finite samples, asymptotics and branch jumps do not determine C. We retain the actual finite profile and make no claim about recovering C.

For ANY tournament T and any k separately reversible edges all incident to one fixed vertex, write H(S) for the path count after reversing exactly S. As a Boolean polynomial H(S)=Σ_{U⊆S} h_U, it has degree at most two: each Hamilton path has at most two edges incident to the center, and its orientation indicator is a product of at most two binary edge literals. Thus h_U=0 for |U|>2, and an unweighted alternating projection over any such k≥3 edge cube vanishes.

For any function g:{0,...,k}→R, define forward differences Δg(t)=g(t+1)−g(t). Exact grouping of supersets gives the **filtered formula**

    Σ_{S⊆[k]} (−1)^|S| g(|S|) H(S)
       = (−1)^k Σ_{|U|≤min(2,k)} h_U Δ^(k−|U|)g(|U|).       (2)

Indeed, for each U with |U|=d, its coefficient is (−1)^d Σ_{j=0}^{k−d} (−1)^j binom(k−d,j)g(d+j)=(−1)^k Δ^(k−d)g(d). The natural MRB choice is g(t)=f(t+1). Formula (2) is a finite, exact analogue of the Euler-transform difference operation, not an equality to the infinite MRB constant. It provides a sharp audit: low-order observed jumps need their mixed h_U sidecars, even when a cancellation makes an unweighted high-order alternating sum zero.

## Explicit seven-vertex certificate and S161 source matching

Take n=5 in (1), fix x_1=0, and write x=x_2,y=x_3,z=x_4. Then q_1=q_2=9, q_3=12, and

    H(x,y,z)=33+18x+36y+60z−18(xy+xz+yz).                 (3)

Bit-mask order 000,100,010,110,001,101,011,111 yields 33,51,69,69,93,93,111,93. Two independent exact algorithms, subset-endpoint DP and all 7! vertex permutations, confirm all eight values. The specific chain 000→100→110→111 yields 33→51→69→93; its equal first two +18 jumps have a nonzero pair interaction −18. Thus this gives an explicit explanation of why source state and mixed differences matter for the S161 matching problem.

The aggregate Hamilton-path counts on the subsets of sizes 0,1,2,3 are 33,213,273,93. Therefore for g(t)=f(t+1),

    Σ_S (−1)^|S| f(|S|+1)H(S)
      = −213 f(2)+273 f(3)−93 f(4)
      = 33+273*3^(1/3)−306*sqrt(2)
      ≈ −6.015217392244597.

Here f(1)=0 and f(2)=f(4)=sqrt(2)−1. By contrast the same cube with constant weight has zero third alternating difference. The nonzero MRB-weighted response is a finite algebraic expression, not C; this cleanly exhibits what the nonconstant finite-difference kernel sees that the unweighted projector destroys.

Run `python3 certificate.py` and `python3 -O certificate.py` for independent seven-vertex enumeration plus all 126 weighted-cut checks. The ordinary finite controls do not establish any claims about all tournaments, the exceptional values 7 and 21, or existence of arbitrary prescribed strong counts.

## A sharp limitation, with an explicit strong H=35 witness

For n≥3, q_1=q_2=1+2^(n−2) (where q_2 exists). For 2≤d≤n−3,

    q_(d+1)−q_d = 2^(d−2) + 2^(n−d−3)(3^(d−1)−1) > 0.

Thus q_d≥q_1 for every d. A nonzero directed cut contributes at least 2+2^(n−1) to (1), and this bound is attained by the word with only x_2=1. In particular 35 does not occur anywhere in the family: n≥6 have baseline ≥65, n=5 has 33 and then ≥51, while the exact n=2,3,4 value sets are {5}, {9,15}, and {17,27,37}.

This is a real restriction, not an unrealizable-spectrum exception: the seven-vertex tournament with adjacency rows

    0010000
    1010011
    0000100
    1110100
    1100000
    1011101
    1011100

is strong and has H=35. The attached certificate checks strongness by transitive closure, independently counts paths by both methods, and checks the small family ranges. Hence even allowing all star words and all n, this explicit strong carrier alone cannot prove the complete spectrum.
