# Classifying interleaved two-site switches: exact cancellations and a sharp scope limit

**Proved.** In the one-apex transitive-spine family, let each independent binary control toggle two distinct sites, the supports of different controls disjoint but potentially interleaved. Two controls have zero mixed difference exactly in the two types below. Consequently additive legal subsets admit only polynomial, not Fibonacci/exponential, capacity. In fact no pair of these controls, even allowing nonzero interaction, yields four consecutive odd counts. This extends the [contiguous-control obstruction](../chanpak_graceful_ternary_20261115/RECTANGLE.md) to interleaved two-site supports; it says nothing about interleaved controls on larger supports.

## 1. Classify zero mixed interaction

If the apex beats the spine site at bit 1, write (H=1+2F), where

    F(w) = Σ_{i<j} c(j−i) w_i(1−w_j),
    c(1)=c(2)=1, c(d)=2^(d−2) (d≥2).

This is the [odd-cycle expansion](../chanpak_graceful_ternary_20261115/GRACEFUL.md). Orient each two-site control to put the first site from 0 to 1. The second changes either 0→1 (symmetric, sign +1) or 1→0 (antisymmetric, sign −1). The mixed half-count difference is minus the kernel bilinear pairing: for disjoint supports X,Y with change vectors δ,η it is −Σ_(i∈X,j∈Y) c(|i−j|) δ_i η_j. The exterior choices do not modify this difference.

Let a<b<d<e be the four sites, with gaps r=b−a, h=d−b, t=e−d, and span L=r+h+t. Write s,z for the second-site signs. Three cases exhaust interleaving:

- Disjoint intervals (a,b),(d,e) give
  (c(r+h)+z c(L)+s c(h)+sz c(h+t)). Replacing c(1) by 1/2 gives exact product (2^{h-2}(2^r+s)(1+z2^t)) of absolute value at least (2^{L-4}). If h=1 the only correction is (s/2). With L≥4 this is nonzero, and the L=3 values are given directly by (1+2z+s+sz), none zero.
- Crossing (a,d),(b,e) gives (c(r)+z c(L)+s c(h)+sz c(t)). For L≥4 the largest term has magnitude (2^{L-2}), while the other three together are at most (3·2^{L-4}); nonzero. At L=3 check (1+2z+s+sz), again nonzero.
- Nested (a,e),(b,d) gives (c(r)+z c(r+h)+s c(h+t)+sz c(t)). If both signs positive all terms positive. If both negative it is (-f_h(r)-f_h(t)), where (f_h(u)=c(u+h)-c(u)). It vanishes only for r=h=t=1. For opposite signs it is up to sign a difference of values of (f_h) or of (g_h(u)=c(u+h)+c(u)). Both are strictly increasing for integer u≥1, by direct evaluation at 1,2 and doubling thereafter. It vanishes exactly for r=t.

Thus the vanishing pairs are **concentric opposite types**, and **two antisymmetric controls on the opposite pairs of a consecutive four-site block**. For a bank of m pair controls the zero graph decomposes by midpoint. Each midpoint component is complete bipartite across signs, plus at most one extra antisymmetric edge (the unique consecutive quartet at that midpoint). Every additive legal subset is a clique. At most three controls can combine, and the number of such subsets is (O(m^2)): at a midpoint with S symmetric and A antisymmetric there are at most (SA+1) edges and S triangles. Three compatible controls are possible (symmetric on sites 1,6 and antisymmetric on 2,5 and 3,4). But disjoint two-site interleaving alone cannot implement an exponentially long exactly additive matching interval. This is a statement about additivity on each allowed subset as its genuine simultaneous value, not about other nonlinear codebooks.

## 2. The zero rectangles already miss four consecutive odds

Four consecutive odds in an additive rectangle would require the absolute half-step sizes to be 1 and 2.

For the nested consecutive antisymmetric exception the core words are (a b (1-b)(1-a)), with internal half-steps 4 and 1. If u is the ordinary binary left exterior and v the reversed complement of the right exterior, inner half-step is (1+u+v) and outer half-step is (4+(7(u+v)-[(u mod 2)+(v mod 2)])/2). Indeed each exterior bit contributes the difference of the two distance weights; the nearest contributes 3 outside the outer pair and 1 outside the inner pair, and all farther terms double. In particular the outer is at least 4.

For opposite-type concentric nesting, consider the antisymmetric pair with separation h0. Its half-step under 01→10 is at least (c(h0)+Σ_(i=1,…,h0−1) min(c(i),c(h0-i))). The inside contribution has one term at every intervening site and the exterior contributions are nonnegative by monotonicity. Thus h0≥3 gives at least 4. If h0≤2, the concentric symmetric partner must be exterior with equal distances r0 on the two sides. Its additional contribution to the antisymmetric half-step is (c(h0+r0)-c(r0)). So h0=2 gives ≥3. For h0=1, only r0=1,2 could give half-step ≤2. At r0=1 the symmetric half-step from the two exteriors is (9(floor(v/2)-floor(u/2))+5(v mod 2-u mod 2)), with zero internal contribution, never ±1 or ±2. At r0=2 the antisymmetric step is at least 2 and exceeds 2 if any active exterior bit exists. Without active exterior, the two fixed interior bits give symmetric half-steps 5,0,0,−5, by the four evaluations of the displayed F formula. None is ±1 as needed.

## 3. Nonzero rectangles reduce to small internal spans and two affine equations

For four consecutive odd H, after ignoring the common offset the four F values are a permutation of 0,1,2,3, whose mixed differences belong to {0,±2,±4}. Thus for nonzero interaction we need |K|≤4. The disjoint and crossing bounds of §1 force L≤6. For nested equal positive types, the bound (c(h+max(r,t))≤4) gives (h+max(r,t)≤4). Both negative types, with (f_h(u)≥c(h+u)/2) when h+u≥3, give (h+max(r,t)≤5). Opposite types with max v>min, the gap of the two (f_h) or (g_h) values is at least (2^{h+v-4}): for v≥3 both functions double from v−1 to v, and for v=2 direct evaluation gives (2^{h-1}). So h+v≤6. In all nonzero cases L≤11.

This bounds the internal word between the first and last toggled sites but **not** the exterior. The latter is handled over all integers by separating its two low parities from the remaining addresses. With internal sites indexed 0,…,L, a control δ has external contribution to its half-step

    − U Σ_i δ_i 2^i + V Σ_i δ_i 2^(L−i)
    − u0 Σ_i δ_i c(i+1) + v0 Σ_i δ_i c(L−i+1),

where the left exterior binary address is (2U+u0), and the right reversed-complement address (2V+v0). For the two controls the desired shifts therefore give two linear equations in nonnegative integers U,V. The [exact verifier](verify_two_site.py) enumerates internal words and the two parity bits, matches each admissible permutation of the four F values, and solves these equations by Cramer's rule or extended gcd with exact bounds. None of the 7,920 residual systems has a solution. This finishes the all-order no-four-consecutive theorem. The inverse problem has not been truncated in absolute address size or spine length.

**Design implication.** The useful invariant is not merely separation or “put the switches in a clever order,” but the entire graph of mixed differences. Here, reflection creates true independent repair channels; the exponential weight kernel limits them to a very special nested geometry. An exact polynomial clock that needs exponentially many *additive* choices must use larger or overlapping supports, a wider reader, or a nonlinear target-selected decoder. This says nothing about the known nonlinear modular packets, which intentionally exploit nonzero interference.
