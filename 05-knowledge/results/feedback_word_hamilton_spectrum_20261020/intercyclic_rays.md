# Intercyclic spectrum: one-apex values plus two dyadic rays

**Status: PROVED, using the elementary reduced-word classification below and the cited two-vertex transversal theorem.** Ordinary Hamilton-path counts only. No global spectrum completeness or priority claim is made.

Let `S_1` be the all-order spectrum of tournaments which become transitive after deleting at most one vertex. Let `S_ic` be the spectrum of intercyclic tournaments: those with no two vertex-disjoint directed cycles, equivalently no two vertex-disjoint directed triangles.

**Theorem.**

`S_ic = S_1 union {3*2^t+1:t>=2} union {9*2^t-3:t>=1}`.       (1)

Both dyadic rays have explicit strong witnesses. In particular the values in `S_ic \ S_1` have counting function `O(log X)` up to `X`.

**Corollary.** Every tournament with `H=29` contains two vertex-disjoint directed triangles. The value 29 itself is attained by a strong six-tournament.

## 1. Exact word condition

By Bang-Jensen–Bessy, [*Small cycle transversals in tournaments*](https://www.lirmm.fr/~bessy/publis/SmallTransversTournament.pdf), Theorem 1.8, every intercyclic tournament has a feedback set of size at most two. The zero- and one-vertex cases are already in `S_1`. Otherwise choose feedback vertices `x,y` with `x->y`, and encode the transitive remainder by letters `0,1,2,3` as in the previous packet: bit one means `x` beats the chain vertex and bit two means `y` beats it.

Let `A` be the pairs of chain positions making a `10` inversion in the x-bit and `B` the corresponding pairs for the y-bit. These correspond exactly to directed triangles using `x` alone and `y` alone, respectively. Every triangle uses at least one feedback vertex. Any triangle using both meets every other triangle. Thus the tournament is intercyclic **iff every edge of A intersects every edge of B**.

If `A` is empty, deleting `y` leaves a transitive tournament. If `B` is empty, deleting `x` does. We may therefore assume both are nonempty. Delete initial `0` letters (actual sources) and terminal `3` letters (actual sinks); this preserves H and neither inversion family is changed.

## 2. Complete reduced grammar

Every reduced word with nonempty cross-intersecting inversion families is one of the following, where `a,b>=1`:

- `1^a20`, `2^a10`;
- `1^a2 1^b`, `2^a1 2^b`;
- `1^a3^b0`, `2^a3^b0`;
- `3 0^a1^b`, `3 0^a2^b`;
- `31 2^a`, `32 1^a`;
- `3 0^a`, `3^a0`;
- the two short words `310`, `320`.

The word `30` appears twice in the displayed parametrization; that harmless duplicate is not a second tournament count.

### Endpoint proof

The first letter lies in `{1,2,3}` and the last in `{0,1,2}`.

**First 1, last 2.** The x-bit has the outer inversion. Every y-inversion would have to meet an endpoint, but the initial y-bit is zero and the terminal y-bit is one, so neither endpoint can occur in a y-inversion. This contradicts `B` being nonempty. The `2,...,1` case is symmetric.

**First 1, last 1.** Every internal x-zero gives the A-edge from the first position, and every internal y-one gives the B-edge to the last position. Their chain endpoints must coincide. Since both families are nonempty, there is exactly one position of either type, it is a letter 2, and all remaining positions are 1. This gives `1^a2 1^b`. Swapping the bits gives `2^a1 2^b`.

**First 1, last 0.** The outer edge is in A. A B-edge cannot use the first position, so every B-edge must end at the last position. Before the last position, the y-bits are nondecreasing. If at least two positions have y-bit one, then any earlier x-zero would create an A-edge from the first position disjoint from a B-edge ending at the last. Thus the initial block is all 1 and the later block all 3, giving `1^a3^b0`. If there is just one y-one position, it must be penultimate; it is either 3 (the preceding family) or 2, giving `1^a20`. Swapping bits gives the two families starting with 2 and ending with 0.

**First 3, last 1 or 2.** Reverse the word and complement both coordinates. This preserves the cross-intersection condition and reduces to the preceding case, giving `3 0^a1^b`, `3 0^a2^b`, `31 2^a`, and `32 1^a`.

**First 3, last 0.** Both inversion families contain the outer edge, so there are no internal inversions in either coordinate. If the interior contains a 0, cross-intersection forces every interior letter to be 0; similarly, if it contains a 3, every interior letter is 3. Otherwise the interior consists of 1s and 2s. Two 1s would give disjoint inversions using the outer endpoints, as would two 2s; one 1 and one 2 create an internal inversion disjoint from the other outer edge. The remaining possibilities are a single 1, a single 2, or an empty interior. These give `310`, `320`, `3 0^a`, and `3^a0`.

Conversely every listed word satisfies the cross-intersection condition by inspection of its two inversion families. This proves the grammar without an unbounded search.

## 3. Remove the families which are actually one-apex

The following deletions leave a transitive tournament:

- `1^a2 1^b`: delete the unique letter 2;
- `1^a3^b0`: delete the terminal 0;
- `3 0^a1^b`: delete the initial 3;
- `3 0^a` or `3^a0`: delete the singleton endpoint.

One can check this directly from the remaining coordinate words: there are no single-feedback inversions, and no letter 2 to form a triangle with both feedback vertices. Equivalently every triangle in the original word uses the deleted vertex.

The remaining nontrivial templates reduce by tournament duality to four cases. To be precise, dualizing all arcs, reversing the chain, and relabelling the two feedback vertices so that their arc is again `x->y` sends the word to its reverse under the letter map

`0 <-> 3`, `1 -> 1`, `2 -> 2`.

It pairs `1^a20` with `32 1^a`, `2^a10` with `31 2^a`, and `2^a3^b0` with `3 0^b2^a`. This operation preserves H.

## 4. Direct odd-cycle count for the four templates

All these tournaments are intercyclic, so no two odd directed cycles are vertex-disjoint. The odd-cycle-collection formula is cited from John Irving and Mohamed Omar, *Revisiting the Redei-Berge Symmetric Functions via Matrix Algebra*, arXiv:2412.10572v3, Corollary 20, PDF page 14, which attributes the formula to Grinberg–Stanley, Theorem 1.39 and Lemma 6.5 ([primary PDF](https://arxiv.org/pdf/2412.10572)). Applying it here gives

`H=1+2(c_x+c_y+c_xy)`,                                    (2)

where the three counts classify odd directed cycles according as they use only `x`, only `y`, or both. A directed cycle is counted once modulo rotation, never once per vertex set.

The exact counts are:

| Template | `c_x` | `c_y` | `c_xy` | `H` |
|---|---:|---:|---:|---:|
| `1^a20` | `3*2^(a-1)-1` | `1` | `3*2^(a-1)` | `6*2^a+1` |
| `2^a10` | `1` | `3*2^(a-1)-1` | `3*2^a-2` | `9*2^a-3` |
| `2^a1 2^b` | `2^(b-1)` | `2^(a-1)` | `3*2^(a+b-1)-2^(a-1)-2^(b-1)` | `3*2^(a+b)+1` |
| `2^a3^b0` | `2^(b-1)` | `2^(a+b-1)` | `2^(a+b)-2^(b-1)` | `3*2^(a+b)+1` |

The single-feedback entries follow directly from the weighted-inversion formula. Here are explicit derivations of the less immediate both-feedback entries. Every cycle through both feedback vertices either uses the arc `x->y`, or is uniquely of the form `x -> A -> y -> B -> x`, with both A and B nonempty increasing chain subsequences. Odd cycle length means an odd total number of chain vertices.

**`1^a20`.** Write `U` for the a letters 1, `z` for 2, and `t` for 0. The only odd cycle using arc `x->y` is `x,y,z`. Without that arc, the choices are:

- A contains a nonempty subset of U followed by t, and B is z: `2^(a-1)` choices;
- A is a nonempty even subset of U, and B is z: `2^(a-1)-1` choices;
- A is an odd subset of U, and B is z,t: `2^(a-1)` choices.

Together with the first triangle this gives `c_xy=3*2^(a-1)`.

**`2^a10`.** Write U for the a letters 2, z for 1, and t for 0. Cycles using arc `x->y` give `2^(a-1)` choices ending in U, and `2^a-1` choices ending at t (whether z is included is fixed by parity). Thus their count is `3*2^(a-1)-1`. Cycles without that arc have A=z or A=z,t; splitting according to whether B includes t gives the same total. Hence `c_xy=3*2^a-2`.

**`2^a1 2^b`.** Write U,z,V for the three blocks. For cycles using `x->y`, omitting z gives `2^(a+b-1)` choices. Including z requires nonempty subsets of both U and V with even combined size; this gives `2^(a+b-1)-2^(a-1)-2^(b-1)+1`. Without the arc, A must be the singleton z and B is a nonempty even subset of U union V, giving `2^(a+b-1)-1`. Their sum is the stated entry.

**`2^a3^b0`.** Write U,V,t for the blocks. Cycles using `x->y` either end in U (`2^(a-1)` choices), or end at t after a nonempty even subset of U union V (`2^(a+b-1)-1` choices). Without the arc, A is a nonempty subset of V followed by t and B a nonempty subset of U. The combined subset size must be even, giving `2^(a+b-1)-2^(a-1)-2^(b-1)+1`. Summation gives the stated entry.

The exceptional short words `310,320` have counts 13 and 15, respectively. These are already the first values of the two rays. They can be checked by the same cycle count or by the explicit six- and seven-cycle counts in the verifier.

This proves the inclusion of `S_ic` in the right side of (1). The reverse inclusion is constructive: one-apex tournaments are intercyclic, `1^(t-1)20` realizes the first ray, and `2^t10` realizes the second. These two word families are strong. For the first, a spanning cycle is `x -> U -> t -> y -> z -> x`; for the second it is `x -> y -> U -> z -> t -> x`.

## 5. The concrete obstruction H=29

The one-apex sharp extrema give:

- at order at most six, H is at most 23;
- at order at least eight, every strong one-apex core has H at least 45;
- at order seven, fixed-content minima for reduced length-six words are 33 (one or five ones), 27 (two or four ones), and 25 (three ones). The strict rearrangement theorem says every nonminimal word lies at least six above its own fixed-content minimum, because (3) of `proofs.md` gives a six-unit minimum exchange at chain length six.

Thus 29 is not a one-apex value at any order. It is also not in either ray of (1), since neither `3*2^t+1=29` nor `9*2^t-3=29` has an integer solution. Therefore it is not intercyclic.

An ordinary strong witness with H=29 has out-neighbour masks

`(50,16,3,7,12,30)`.

The independent verifier checks all 720 permutations, subset DP, strong connectivity, and a pair of disjoint directed triangles. This establishes a genuinely new boundary for the proposed additive/intersecting-cycle strategy: it already fails at the globally attainable value 29, so it cannot discharge the two prime lanes.

## Reproduction / provenance

`verify_intercyclic.py` checks the reduced grammar through chain length eight, both ray formulas by direct directed-odd-cycle enumeration and independent Hamilton counts for small parameters, the exceptional words, and the H=29 witness. Its finite checks validate an all-order proof rather than replacing the endpoint classification or the exact cycle counts.

The structural kernel was developed jointly during independent proof review: the cross-inversion interpretation and endpoint classification were combined with a two-vertex inversion-cover refinement, the exact dyadic-ray reduction, direct cycle counting and the H=29 obstruction. No external mathematical novelty claim is asserted pending broader review.

## 6. Exact triangle-packing number at H=29, and a divisibility hostile

The permanent oddness theorem and SCC product law imply that a tournament with prime count 29 has exactly one nontrivial strong component, whose whole count is 29. By the Moon–Busch minimum theorem, a strong component of order at least eight has count at least 45. Thus the nontrivial component has at most seven vertices and can contain at most two vertex-disjoint triangles. Section 5 supplies the matching lower bound. Consequently

`H(T)=29 => maximum number of vertex-disjoint triangles in T is exactly 2`.

Singleton SCCs do not change this packing number. The cutoff is imported from the correctly cited Moon–Busch theorem in the inherited spectrum packet, not from a fitted recurrence.

The whole-count hypothesis matters. The strong one-apex word `1100100` has `H=87=3*29`; its out-neighbour masks are

`(126,124,248,240,96,192,128,19)`.

All its cycles contain the apex, so its triangle-packing number is one. The verifier checks its count by the binary formula, subset DP and all 40,320 vertex permutations. Therefore divisibility by 29 does **not** force two disjoint triangles.

## 7. Positive-density obstruction for the entire intercyclic strategy

This argument was developed jointly during independent proof review.

Let `X=2^L`, `L>=5`. Every nontransitive one-apex count at most X has a reduced word of length at most `L+1`, since the sharp minimum is at least `2^(m-2)+1` for reduced length `m`. There are

`sum_{m=2}^{L+1} 2^(m-2)=X-1`

such candidate words before the upper-count filter.

The involution `b -> 1-reverse(b)` preserves H, by tournament duality. It has no fixed word at odd length, and `2^(m/2-1)` fixed words at even length m, so at most

`2^floor((L+1)/2)-1`

fixed words occur in this whole bank.

At the maximal length `m=L+1`, every word starting in `11` and ending in `00` has

`H >= 1+9*2^(m-4) > X`.

The lower bound counts only the four inversions between those two initial and two terminal letters. At `m=4` the adjacent-pair correction only strengthens it. This excluded family is stable under the involution and has `2^(m-4)=X/8` words.

After excluding it, the number of involution orbits is at most half the sum of the remaining word count and the original fixed-word count. Add the transitive value one. This proves the explicit bound

`|S_1 intersect [1,X]| <= 7X/16 + 2^(floor((L+1)/2)-1)`.

For `L>=5`, the two dyadic rays in (1) contribute at most `L-3` and `L-4` additional values, respectively. Therefore

`|S_ic intersect [1,2^L]|`

`<= 7*2^L/16 + 2^(floor((L+1)/2)-1) + 2L-7`.                (3)

There are exactly `2^(L-1)` positive odd integers at most `2^L`. At least

`2^L/16 - 2^(floor((L+1)/2)-1) - 2L + 7`                    (4)

of them lie outside the intercyclic spectrum. In relative terms, at least `1/8-o(1)` of the positive odd integers below these dyadic cutoffs cannot be realized by intersecting-cycle tournaments.

This is **not** a positive-density claim about forbidden ordinary H-values. A value outside `S_ic` might be unattainable, or it might require at least two vertex-disjoint triangles. The theorem proves that the proposed `H=1+2*(number of odd cycles)` intercyclic strategy cannot establish full odd-spectrum coverage; it does not refute that coverage.
