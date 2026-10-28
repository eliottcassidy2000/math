# Self-converse lift: exact inheritance and complexity status

## Repository proof recovered

At source pin `03408f84ac43fcdcfb0887e5d6bc1fc8792415cf`, the exact files are:

- `05-knowledge/results/converse_tournament_count_channels_20261021/session22-converse-channels/README.md`;
- `05-knowledge/results/converse_tournament_count_channels_20261021/session22-converse-review/COROLLARY_REVIEW.md`;
- `01-canon/theorems/THM-283-sc-burnside-formula.md`.

Their bytes are copied under `sources/` here. Session 22 proves the raw simultaneous tuple/color Burnside interpretation of `F_n(a,b)`, and proves `C_(2n)=F_n(2,0)`, `C_(2n+1)=F_n(2,1)` by twisted Burnside. Its review separately proves the equivariant codec from a tournament tuple to one retained tournament and the indexed disagreement graphs. The review **expressly disclaims a structural/canonical bijection from that codec to unmarked self-converse tournaments**. Consequently a new explicit isoclass lift needs its own proof of surjectivity and independence from the chosen anti-involution; it is not already established by the inherited count equality.

The S22 source also corrects an old explanatory error in THM-283: allowed anti-permutations can have odd-length unordered-pair orbits (a transposition's antipodal edge is fixed). The consistency condition is equality between orbit-length parity and endpoint-swap parity. This matters when constructing the lift rather than merely summing a formula.

## A simple independent odd-order TI reduction

For an arbitrary tournament T on n vertices define `W(T)=T ⊕ {v} ⊕ T^op`, with every earlier block beating every later block. This is self-converse under exchanging the first and last copies and fixing v. Moreover

```
T ≅ S  iff  W(T) ≅ W(S).
```

Forward is immediate. Conversely an isomorphism preserves the uniquely ordered strong components. In W(T), the initial union of strong components of total size n exists and is exactly T. An isomorphism sends it to the equally positioned union in W(S), recovering S. This gives an explicit polynomial-time many-one reduction from ordinary TI to isomorphism of odd-order self-converse tournaments. The reverse direction is the identity inclusion of that promised class into TI.

**Source / target / map:** tournament isomorphism pairs map to the ordinal-sum construction. **Preserved:** yes/no isomorphism and an explicit recoverable prefix. **Lost:** no input isomorphism information; the construction does not cover every self-converse tournament and is not the general tuple-to-self-converse class bijection.

This elementary complexity fact is not claimed novel, and does not give GI-hardness, a polynomial TI algorithm, or a canonical decoder for every self-converse tournament. It also does not establish anything about *recognition* of self-converse tournaments merely by proving isomorphism hardness on the promised class.

## Literature-priority check, not an asserted new complexity theorem

A bounded literature search located the historical title *Isomorphism problems involving self-complementary graphs and tournaments*, C. J. Colbourn and M. J. Colbourn, Proceedings of the Eighth Manitoba Conference / Congressus Numerantium XXII (1979), 153–164. The author's [publication list](https://www.public.asu.edu/~ccolbou/src/cv.pdf), item 106, verifies the bibliographic item. The original paper was not recovered in this lane. Therefore it is a priority lead, not a claimed primary-source verification of every variant.

A secondary survey explicitly discusses TI versus self-converse restrictions, and its updates flag a problem in some additional regular/digraph construction claims. Those leads justify avoiding novelty language; the short ordinal-sum proof above is the mathematical justification used here. No result about the new full lift is inferred from the secondary source.
