# An actual bounded +2 deposit at arbitrary height, with a different multiplier

There is a positive unit test outside the rigid count-three
companion class. Use the reverse bridge with the old root
retained. Take the two right readers below, both rooted at
vertex 6 (zero-indexed rows). Outgoing masks are given in
check_unit_reader.py.

| Reader code | H | P₂ optional | E_s,O_s | E_e,O_e |
|--:|--:|--:|--:|--:|
| 4578838542729556224 |109|2848|158,158|59,72|
| 4578836343708259584 |109|2848|158,158|60,71|

The count of either right reader and the complete retained-root
transfer match. The right end totals are both 131, but the
signed parity changes from −13 to −11. Consequently, for
an arbitrary incoming state (H,X,D,g) the two outputs satisfy

    H'_0=109H+131X+13D,
    H'_1=109H+131X+11D,
    X'=2848X−316g, D'=0, g'=109g.

Thus **if D=−1, the new deposit is exactly +2** without
changing the full continuation. This is a change of reader,
not a large determinant with merely the right residue.
It remains exact as the entering H grows without bound.

A literal ready family: reverse the terminal arc across m
intervening sites of a transitive tournament, and root at
the original final vertex x. The old order contains only
that backward arc. Root deletion is transitive, and x
beats only the first spine vertex, so
(H,X,D,g)=(1+2^m,1+2^m,−1,1) for m≥1.
This gives the exact unbounded pairs

    (240(1+2^m)−13, 240(1+2^m)−11).

For m=4 the source H is 17 and the pair is 4067/4069.
On the separate singleton left input the same readers give
227/229, both prime. Tests use literal subset DP on joined
tournaments as well as an atlas scan that matches the
full signatures, rather than inferring a +2 from the
desirable right end difference alone.

**Scope:** the multiplier is **109, not 3**, and afterwards
D'=0. The unit reaction is therefore not already a closed
ternary lift or an independently stackable register. In
the fixed-root language a subsequent operation producing
nonzero D from g uses a multiple of 109. A complete
iterative rule must refresh the old unit from another
coordinate or choose a target-dependent splice.

The search uses every root of every order-eight atlas
representative, with keys (H,P₂,E_s,O_s,E_e+O_e). It finds
four (key, adjacent end-even levels) pairs, with right
Hamilton counts 109,197,319,345. This is a finite discovery
inventory, not a classification at all orders. The displayed
unit deposit requires only the two explicitly verified readers.

There is also an **all-strong version**. Codes
8871896820243039496 and 8007222596720985376, both at root 6,
are strong and have right H=197, P₂=3392 and starts (318,333).
Their ends are (111,156) and (112,155), so on the same
incoming D=−1 family they give (464H−45,464H−43). For
H=17 the joined tournaments are strong and give 7843/7845.
Their common continuation has D'=−15g; it still does not
re-arm a unit at the next step. This strong construction
avoids relying on a disconnecting or isolated shield.
