# Exact mixed profile transfers and a Mersenne-divisor inverse

Write H(T) for its Hamilton path count, including H(empty)=1, and
    Z_k(T)=sum_{f:V(T)->[k]} product_{i=1}^k H(T[f^{-1}(i)]).
It counts k ordered (possibly empty) disjoint directed paths partitioning the vertex set. Thus Z_1=H; for a nonempty tournament Z_0=0. Ordered joins multiply the full Z-profile: in each color the forced part order makes the Hamilton paths multiply.

This profile is the principal specialization of the [Grinberg–Stanley Rédei–Berge symmetric function](https://arxiv.org/abs/2307.05569). In particular its exact odd-cycle packing expansion is
    Z_k(T)=sum_{P} 2^{c(P)} k^{n-v(P)+c(P)},
where P ranges over packs of vertex-disjoint nontrivial odd directed cycles, c is cycle count and v covered order. One may apply their H expansion in each color and count colorings of packs; v-c labels are identified. Thus Z is a polynomial with nonnegative coefficients, degree n, leading coefficient one, and no constant if n>0. The empty packing gives Z_k>=k^n.

## Physical mixed combiner and exact inverse
Let U,V point around a cyclic triangle with a new singleton b, in the order U→V→b→U. If both U,V are nonempty, a Hamilton path with b at an endpoint contributes 2H(U)H(V). An interior b permits a path of the form U_before,V_before,b,U_after,V_after, with U_after and V_before nonempty, and either other slot possibly empty. Its contribution is (Z_2(U)-H(U))(Z_2(V)-H(V)). Thus
    H(C3[U,V,1]) = 2H(U)H(V)+(Z_2(U)-H(U))(Z_2(V)-H(V)).

In fact, for every k>=1,
    Z_k(C3[U,V,1]) = k [Z_k(U)Z_k(V)
       + (Z_{k+1}(U)-Z_k(U))(Z_{k+1}(V)-Z_k(V))
       + (Z_k(U)-Z_{k-1}(U))(Z_k(V)-Z_{k-1}(V))].
Fix the color of b. Other colors concatenate their U and V parts normally. In b's color apply the preceding expression, correcting when a part is empty: subtract 1_{U empty}H(V) and 1_{V empty}H(U), and add 1_{both empty}. Summing the endpoint term over all colorings gives 2Z_k(U)Z_k(V), and the corrections give the three products with Z_{k-1}. The W before/after split gives the Z_{k+1}-Z_k terms (one previously available color is split in two; the mandatory subslot excludes the unsplit count). Collect to get the displayed formula and multiply by k choices of b's color. For V singleton, the expression becomes k² Z_k(U)+k Z_{k+1}(U)-k Z_{k-1}(U).

If a fixed profile has H=H(U), B=Z_2(U), append a t-vertex transitive ordinal pad and take V a chain of m>=1 vertices. Then the honest tournament height is
    N = 2H + (2^m-1)(2^t B-H),   t>=0.
Exact inversion: if N<=2H, reject; otherwise enumerate Mersenne divisors d=2^m-1 of N-2H, and test whether ((N-2H)/d+H)/B is an integral nonnegative power of 2. Enumerating m with 2^m-1<=N-2H avoids any factoring assumption. Concrete (H,B): regular C5=(15,120); direct Petersen chord table=(5213,104408); original doubled C5=(6473,134288); compressed Petersen=(8265,155040); first native direct-Petersen descent neighbor=(1747,50656). For direct Petersen with m=1, the first four t produce 109621,214029,422845,840477. The regular C5 profile has gcd 15 and gives only multiples of 15. The direct Petersen and doubled C5 profiles have gcd 1, but are not universal charts: a finite bank of K such two-parameter families covers at most K(r+1)^2 residues mod 2^r, far fewer than the 2^{r-1} odd residues.


See [the general sparsity theorem](bounded-templates.md) and [literal recurrence replay](moore-chanpak-profile-check.py).
