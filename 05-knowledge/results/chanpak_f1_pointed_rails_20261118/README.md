# The transitive point, six uncolored sheets, and seven binary-colored sheets

This develops the [two-rail switch construction](../chanpak_two_rail_pell_switches_20261118/README.md) in the pointed-set model of the field with one element. The repo's [older ground-state picture](../../../07-reflections/the-fundamental-truths-in-every-language.md#f₁-geometry-hyp-1539) identifies H−1 as departure from the transitive count. The precise development below uses finite pointed sets, as in [Szczesny](https://arxiv.org/abs/1006.0912). It neither posits an ordinary field of cardinality one nor says that every tournament has a distinguished Hamilton path.

## A natural transitive basepoint on the packing side

Let Z_T(z) count vertex-disjoint directed odd-cycle packings, marking each cycle by z; unused vertices require no mark. There is a canonical basepoint, the empty packing, and [Grinberg–Stanley](https://arxiv.org/abs/2307.05569) gives H(T)=Z_T(2). The three specializations are distinct: z=0 keeps just the basepoint, z=1 forgets cycle colors, and z=2 is the Hamilton count. For nonnegative integer z, decorating each cycle with one of z colors gives a literal finite pointed set whose rank (number of nonbase elements) is Z_T(z)−1.

Ordinary ordinal join makes these packing sets Cartesian products, not smash products. The product of pointed sets decomposes into the common basepoint, the two axes (exactly one input nonbase), and the smash part (both inputs nonbase). So for reduced counts x,y the composition is x+y+xy, the multiplicative formal-group law in coordinates shifted by one. This is not the linear map x↦7x. The two-rail construction below is a different operation entirely: replication on the nonbase part.

## Exact colored factor, not just a count coincidence at two

Put the vertices in the order 0,A,B,C,4,5,6, where A is a nonempty transitive block of width a and B,C are single vertices. All arcs are forward except B,C→0. This is the low tournament L. The upper tournament U also reverses 0→4 and 4→6. Put d=3·2^(a−1)−1. Then the exact polynomials are

    Z_L(z) = 1 + d z,
    Z_U(z) = 1 + 5d z + d z²,
    Z_U(z) − 1 = (5+z)(Z_L(z) − 1).

Proof: in L, every cycle uses 0, a nonempty ascending subsequence I of A and a nonempty subsequence J of B,C. If |I| is odd, J has size one (two choices); if |I| is positive even, J must contain both. There are 2^(a−1) odd subsets and 2^(a−1)−1 positive even subsets, giving d odd cycles.

In U without 0 the only nontrivial directed cycle is the tail triangle 4→5→6→4: the other vertices are still in forward order before the triangle. A cycle through 0 using A has the form 0 I J followed by either no tail (J then nonempty to close at 0) or the tail suffix 4,64,564, all ending at 4→0. Count |I|+|J|+|tail| even. For odd I the counts are 2 without tail and 6 with tail; for positive even I they are 1 and 6. Without A the two odd cycles are 064 and the tail triangle. Thus U has 8·2^(a−1)+7·(2^(a−1)−1)+2=5d single cycles. Any pair of disjoint odd cycles must contain the tail triangle and one original cycle avoiding the tail, exactly d choices. No triple exists.

This actually lifts to a uniform, support-controlled pointed bijection. Mark the first vertex of A and fix the set J0 of all other A vertices used in a cycle. If J0 is nonempty, the optional marked vertex gives two parity choices: there are 3 low single cycles and 15 upper single cycles in this fiber. If J0 is empty there are respectively 2 and 10 (the upper two cycles outside A are counted here). Within each fiber, sort by the literal vertex word from its least vertex and match the upper single cycles to five indexed copies of the low cycles. This is explicit and independent of the size of the unused part of A, preserving membership outside the marked vertex. For an upper double pack, remove its tail triangle to read the low cycle; its color specifies one further sheet. The empty packing maps to the common basepoint. The bijection commutes with order-preserving relabelings and with enlarging the ambient ordered block through unselected vertices while retaining the mark. It does not state that the five-copy ordering is canonical without the marking and linear order.

Thus z=1 gives six copies of the uncolored nonempty part, z=2 gives seven copies of its binary-colored counterpart. More generally for a finite field k one may express the extra sheet label as a disjoint union of five points and the affine line over k. Although `5+2=7=|P²(F₂)|`, the polynomial `5+z` is **not** `1+z+z²`; a projective-plane incidence or GL(3,2) symmetry does not follow. The original 7H−6 is seven wedges sharing a single anchor rather than seven disjoint copies with separate anchors.

## On actual Hamilton paths: yes with slots frozen, no with physical deletion

The common path `0 A B C 4 5 6`, with A filled in its increasing order, works in both tournaments at every a and provides a genuine basepoint on each path set. Compress each maximal A run to one A. There are at most two runs (0 is the only external vertex that can precede an A run). For one run, all A vertices go there. For two runs the allocation datum is a nonempty proper subset of A for the first run; the second is its complement. Template enumeration (templates.py) gives:

| Number of A runs | Low | Upper |
|---|---:|---:|
| one | 5 | 29 |
| two | 3 | 21 |

Removing the same base template leaves 4 versus 28, and 3 versus 21, respectively. A pointed sevenfold bijection of literal paths can therefore be chosen by template rank and leaves the entire two-run allocation unchanged. It commutes with a **slot-preserving** restriction: delete named A vertices if all the original slots still have vertices, and send to the basepoint otherwise. But physically removing a slot can leave a perfectly valid surviving path. Is there also a pointed sevenfold lift compatible with *actual* deletion?

No. Look at width two and the two restrictions obtained by keeping only its first or only its second clone. Delete the unused vertex from the path; if the resulting order is not a path, send it to the common basepoint. These are pointed maps (and such restrictions compose throughout the nonempty width family: losing a slot is the only potential break; after it is gone, all further nonempty restrictions have the same surviving one-run template). Count paths with **exactly one nonbase image** under the two restrictions. A one-run template never contributes to this statistic. Each two-run template has two allocations at width two and contributes two if just one of its faces is nonbase.

- Low templates: two with one face, one with both. The statistic is 4.
- Upper templates: sixteen with one face, five with both. The statistic is 32, not seven times four (28).

The check is stronger than using the given anchor. Under any coherent basepoint choice the low possibilities are {2,4,6} and the upper possibilities are {30,32,34}; no element in the second set is seven times an element in the first. This can be read off from the printed face tables: one-run templates map the same way to the two singleton copies, and no two-run template maps to the same valid rendering on both faces. Thus any coherent basepoint must come from a one-run template at width two. A pointed isomorphism commuting with these physical restrictions is impossible even without presuming that it preserves run data.

This gives an explicit repair demand for a stronger sevenfold mechanism, separate from adding six anchors to the count: new paths that behave as persistent one-run anchors do not change the face statistic. With the original basepoints and the same two-run template total, a deletion-compatible product needs to turn two one-face upper templates into two-face templates (a necessary marginal change, not sufficient for a full equivalence). The field-with-one-element viewpoint helps locate which structure the counted replication retains and which it loses.

Run:
```sh
python templates.py | diff -u templates.out -
python packing_replay.py | diff -u packing_replay.out -
```
