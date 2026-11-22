# The exact gap around the odd-exponent lower rail under backbone reorientation

**PROVED by a finite symbolic enumeration.** Fix an odd integer `p=2t+3≥5` and write `E=2^(p−1)(2^p−1)`, `L=E−1`. Start with the eight-vertex lower-rail source with outgoing masks `(2,8,3,37,15,151,63,95)`. Inflate the ports `(6,4,1)` to transitive blocks of widths `(p,t,t)`. Permit arbitrary reorientation of the ten arcs among the five remaining vertices, keeping all port interfaces fixed. This gives 1,024 labelled sources. [Inherited strong lower rail](../chanpak_support_birth_grammar_20261122/ODD-PERFECT-PORT.md)

For **every** such exponent no resulting tournament has count `E+1`. For every odd `p≥9` the immediate distinct neighbors of `L` among all these counts are

\[
\begin{aligned}
 H_{\rm below}&=L+2-2^t+16\cdot4^t-8\cdot8^t,\\
 H_{\rm above}&=L+2+16\cdot4^t=E+2^{p+1}+1.
\end{aligned}
\]

At `p=5` the two changes from `L` are `−16,+18`; at `p=7` they are `−258,+252`. At exponent 61 the nearest larger count in this entire family exceeds the desired upper endpoint by exactly `2^62`.

The upper boundary has a factored form. With `U=2^(p−1)`, the two attained counts are `L=(U−1)(2U+1)` and `H_above=(U+1)(2U+1)`. To interchange the first factors, reverse the two arcs joining vertex 3 to vertices 0 and 2 in the original skeleton. Reversing only the arc between 2 and 3 gives the lower boundary. Each altered skeleton is strong. Indeed the original, the singly altered and the doubly altered skeletons all contain the common directed cycle `0→1→3→5→7→6→4→2→0`. Their nonempty port inflations remain strong.

## Symbolic calculation for every width

The three ports induce a transitive tournament (`6→4→1`, `6→1`). After compressing each maximal visit to a port into one symbol, a path visits each ordinary vertex exactly once. Between consecutive ordinary vertices, and before the first and after the last, the port symbols form a subsequence of this three-term order. Thus every port has at most six runs, independently of the ten changed arcs.

For each oriented backbone the replay enumerates compressed paths by their last symbol, visited ordinary mask, and three run counts. A port of width `w` appearing in `r` runs has as many fillings as surjections onto `r` labelled run positions. Every assigned set is read in increasing order. There are exactly

\[
             \sum_{s=1}^r(-1)^{r-s}{r\choose s}s^w
\]

fillings, including zero when `r>w`. Multiply these three counts and sum against the compressed path multiplicities. Because the widths are `(2t+3,t,t)`, each product summand is an integer coefficient times `(i²jk)^t`; its coefficient includes `i³`. The resulting polynomial in exponential terms gives the exact path count at every positive `t`. Subtracting `L=32·16^t−4·4^t−1` leaves 491 distinct coefficient vectors across the 1,024 inputs.

For every nonzero expression `D(t)=Σ c_b b^t` with largest base `B`, its sign is that of `c_B` once `|c_B|B^t>Σ_(b<B)|c_b|b^t`. The inequality persists at every greater exponent by division by `B^t`. Applied to the changes, to their differences from two, and to their differences from each stated boundary, every required strict inequality starts by `t=4`. The replay checks every earlier exponent directly. It classifies the zero expressions separately. This proves the uniform ordering and exclusions with finite integer arithmetic.

In lexicographic order on the ten unordered pairs of `(0,2,3,5,7)`, let bit zero reverse the first pair. Masks 0 and 1 give `L` for all widths on the rail. The original vertices 0 and 2 have identical contacts with all other vertices, so reversing their mutual edge amounts to interchanging their labels. The only additional occurrences of `L` have `p=5` and masks `3,16,148,149,162,163`. The displayed lower and upper formulas themselves have masks `3,16` and `18,19`, respectively. The source implementation prints the separate formula and exact exponent-61 value for every one of the ten single-arc reversals.

More generally, with the port identities fixed, the backbone vertices have contact words `001,001,000,011,111` in the port order. Thus only 0 and 2 can be interchanged, and this relabeling acts freely on the 1,024 sources because it reverses their mutual edge. Of the 512 resulting isomorphism classes, the full three-port grammar distinguishes 505 signatures. Exactly 499 signatures have a unique class, five have two classes and one has three. On the perfect-number width line this coarsens to 491 formulas, with respectively 471, 19 and one of multiplicity one, two and three. The loss from isomorphism to grammar and the further loss from evaluating on one width line can be counted separately.

## The exponent-61 inverse with three independently varying widths

**FINITE-EXACT.** Keep the same port interfaces and any of the 1,024 oriented backbones, and now allow **every** positive width triple `(a,b,c)`. The upper target `E_61+1` still has no representation. The lower target `E_61−1` has precisely the two twin-relabelled sources, both with the original widths `(61,29,29)`.

The 505 signatures have run counts bounded coordinatewise by `(2,3,2)`. For fixed `b`, a middle port with one, two or three runs contributes `1`, `2^b−2` or `3^b−3·2^b+3` fillings. Consequently every exact count has the form

\[
       d_{00}+d_{10}(2^a-2)+d_{01}(2^c-2)
                  +d_{11}(2^a-2)(2^c-2),
\]

with four nonnegative coefficients determined by `b` and the source signature. For a positive target `N`, a slot wider than `floor(log₂(N+2))` can only contribute one-run terms in a solution if its width is at least three: a two-run filling alone costs `2^w−2>N`, and three runs cost at least that much for `w≥3`. Such a slot could be reduced to width one without changing the count. For both targets the width limit is 120. Loop over the middle and final widths up to this limit, then solve the displayed linear equation for `2^a−2` and check that adding two gives a power of two. A zero coefficient of `2^a−2` simply permits taking `a=1`. No enumerated hit at either target has a width equal to one. This excludes additional solutions hidden by normalization. Altering the port interfaces or the internal tournament of a port gives a different family.

Run `python3 check_perfect_backbone.py` from this directory.
