# A sharper growth bound for the exact ±7 one-exit route

This strengthens the [contact passport](CONTACT_PASSPORT.md). Let R be the deletion, g=H(R), and Q its tagged vertices. The one-exit reader is WP₂(B)+KH(B) with K=g−U, where U counts ordered two-covers of R whose end of the first path and start of the second path are both tagged. For an exact additive offset ±7, H(B)=1 because H(B) is odd and cannot be 7; hence B is transitive and K=±7.

**New bound.** For p≡3 or 5 (mod 8), a witness for 7p in this exact offset architecture requires deletion count g at least a constant times `p^η`, where

    η = 1/(1+3 log₅ 2) = 0.4362945258…

This replaces the exponent about 0.257 on this route. It is not a lower bound on other representations of 7p.

## The offset confines the marked components

For any nonempty marked component C, its weighted start count W_s(C) is at least `|Q∩C| H(C)`. For each full path cut immediately before the tagged vertex; the tag and the path are recoverable from the ordered cover by concatenation. Cutting immediately after a tagged vertex similarly gives W_e(C)≥|Q∩C|H(C).

If F and L are distinct first and last tagged strong components of R, use a weighted start cover of F and a weighted end cover of L. The counted global cover has its second path start in F and first path end in L. All intermediate components split arbitrarily between the slots; the earlier prefix is forced to the first path and the later suffix to the second. These choices are independent and the ordered slot covers can be recovered blockwise. Hence they give at least

    g · |Q∩F| · |Q∩L| · 2^(number of intermediate components)

covers counted by U, using P₂(C)≥2H(C) on every intermediate component C.

For K=+7, U=g−7<g, so every tag lies in a single component. For K=−7 and g>7, U=g+7<2g; tags either lie in one component, or occupy precisely two adjacent components with one tag in each. The finite small-g table in the earlier note already rules out g≤7 on the sevenfold route.

A wholly untagged prefix or suffix multiplies both g and U by its full path count, so it multiplies K by the same count. Since K=±7 and no tournament has count 7, the combined outer full count is 1. All those exterior components are thus singletons. Trimming them changes W only by a power of two. The trimmed deletion has either one nonsingleton strong component, or two adjacent components with at most one singleton when g>7.

## Bounding the odd part of the slope

For a nontrivial strong component C on n vertices with h=H(C), the [odd-pack formula](https://arxiv.org/abs/2307.05569) gives

    P₂(C) ≤ 2^n + (h−1)·2^(n−2) = (h+3)·2^(n−2).

Indeed the empty pack contributes 2^n. Any nonempty pack has at most n−2 atoms (its cycles and uncovered singletons), as even one directed odd cycle replaces at least three singletons by one atom. Each pack contributes its H-weight multiplied by two raised to that number of atoms. Since h≥3, the ratio P₂(C)/h is at most 2^(n−1), and [Busch's bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) makes this at most `h^(3 log₅ 2)`.

On the trimmed deletion R*, W≤P₂(R*). Two-slot counts multiply over ordinal components; the possible singleton costs only a factor 2. Thus

    oddpart(W) ≤ g^(1+3 log₅ 2)      for K=+7,
    oddpart(W) ≤ 2g^(1+3 log₅ 2)     for K=−7 and g>7.

For an output 7p the slope satisfies `oddpart(W)=7 oddpart(p∓1)`, since P₂(B) is a power of two. For p≡3 or 5 (mod 8) this is at least 7(p−1)/4 in either sign. Combining the inequalities proves the stated exponent; for example the common bound is `g ≥ [7(p−1)/8]^η`.

## In the two-component case, only one can grow

Suppose K=−7 and the marked interval comprises adjacent components F,L with full counts h_F,h_L and one tag each. Write W_s(F)=h_F+x, W_e(L)=h_L+y, and let E(F),S(L) count full paths ending and starting at their respective tags. Both same-component U terms vanish because endpoints of disjoint paths cannot be the same tag. The contact ledger is therefore the exact budget

    7 = h_L x + h_F y + xy + E(F)S(L).

If both h_F and h_L exceed 7 this forces x=y=0. Then the split injections are bijections. A weighted start cover of F whose second path is just its tagged vertex must concatenate, showing E(F)=H(F−tag). Dually S(L)=H(L−tag). Their product would have to be 7, impossible because a factor 7 is forbidden. Thus at most one tagged strong component can grow; an adjacent partner has count 1,3 or 5.

There is a further one-sided restriction when F is the large component and L the small one. Here y=0 and x≤7/h_L. The actual trimmed slope is

    W = P₂(L)(h_F+x)+h_F W_s(L).

The small component has order at most four by Busch (or is a singleton), so P₂(L)≤16h_L and W_s(L)≤P₂(L). Hence `W≤32g+112`. In that orientation, for p≡3 or 5 (mod 8), the same odd-part test forces **linear**, not just polynomial, growth of g. In the special case that L is a tagged singleton, the formula is `W=3g+2x≤3g+14`. This leaves a single strong core, or a large last tagged component with a small preceding partner, as the less restricted configurations. None of these inequalities selects the core or proves a new coverage theorem.
