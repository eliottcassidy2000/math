# Physical inverse tests: cut the forced component, then repair the route

The [signed shortcut chart](PAIRED_SHORTCUT_WINDOWS.md) constructed both counts at each of its successive heights. This note asks a stronger question: when are the two graphs related by a **literal edit**? A one-arc edit and a merged metagraph line are different operations, and each has its own exact failure test.

## An all-order certificate from strong components

**Single arc.** Suppose reversing one arc changes H from h to N. If the endpoints are inside one old strong component, everything outside that component has a common ordinal multiplier d. If the endpoints belong to different old strong components, take the interval from the first endpoint's component to the last endpoint's component. Reversal makes that entire interval strong unless it just swaps two adjacent singleton components, which leaves H unchanged. Indeed the forward route from the first endpoint to the last survives via an intermediate component or via a nontrivial endpoint component, and the reversed edge returns; the ordinal arcs and paths inside the endpoint components give routes through every other vertex of the interval. Everything outside the interval remains an ordered common factor. Thus any genuine one-arc change has a restricted witness whose counts are h/d,N/d and with at least one endpoint strong. The factor d divides gcd(h,N). Conversely any such witness with a realizable multiplier d can be padded ordinally to give the original pair, with the same single changed arc.

**Actual merged line.** A line complements exactly the nonpath arcs of a chosen directed Hamiltonian frame, then merges isomorphism and converse classes. At order n≥3 at least one end of every line is strong. If the source is not strong, any Hamiltonian path traverses its ordered strong components, so the first-to-last arc points forward across two different components. It is not a frame arc. Complementing it makes the frame close to a Hamiltonian cycle in the other endpoint. This is true for blue and black lines, independent of the chosen frame. Therefore any merged line, and any nonempty same-order path of merged lines whose counts stay at most C, has order at most the strong-component bound at height C.

By [Busch's exact lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3), a strong tournament of order at least 9 has at least 75 Hamiltonian paths. Thus for endpoints at most 73, direct single-arc tests reduce to connected core intervals of order at most 8; scanning all divisor-scaled pairs is unnecessary in the signed case because the parent and target are coprime. Indeed if d divides x and its signed stripped parent, it divides both x and 3x±1 and hence is 1. For the adjacent-odd +2 controls, the endpoints are also coprime. Every capped merged-metagraph test through 73 is likewise all-order exact: no larger order can hide a path below the cap. A path with **uncapped overshoot** has no such exclusion.

## An exact failure atlas and a choice of new parent

Use the complete archived [order-seven representatives and order-eight census](../actual_hamilton_multiplicities_20261021/session24-multiplicities/). Induced restrictions of order-seven representatives supply all smaller orders (adjoin singleton sinks to extend any smaller tournament to order seven). The replay recounts literal endpoints. For the metagraph, it enumerates Hamiltonian frames, complements only their nonpath arcs, and canonicalizes up to relabeling and converse. Edges have an explicit frame certificate. A blue flag means frame reversal is an anti-automorphism, not simply that the scalar endpoint is self-converse.

Among the 31 legal signed shortcut requests through target 73 (ignoring requests whose target or prescribed arithmetic parent is 7 or 21), there are **11 all-order single-arc failures**:
```
1→11, 1→43, 5→27, 13→17, 13→35, 17→45,
23→61, 25→67, 29→39, 35→47, 47→63.
```
In the actual merged metagraph, 25 of the 31 are direct lines. Four more have increasing-count detours not exceeding the target:
```
1→9→11;   13→33→35;   1→17→37→43;   11→51→59.
```
Only **13→17 and 47→63** cannot be routed without overshoot. For 13→17 the smallest possible cap is 19; at order 6 the route 13→19→17 uses two black lines. For 47→63 no route is possible even with cap 73. Rather than let the arithmetic parent dictate a failed move, there are genuine replacements:
- 1→17 is the order-six transitive blue line.
- 29→63 is an order-eight black line: from masks (0,1,11,51,71,87,143,63) use frame (2,3,5,6,7,4,1,0); the partner has masks (252,237,248,224,42,192,144,16).

In fact every admissible target through 73 other than the unit has some direct line from a smaller count. This is a finite exact observation about actual node representatives, not a statement that an arbitrary previously chosen representative at the smaller count has the right continuation. It separates the scalar-parent choice from the graph and frame that actually lift it.

## Constructive larger tests, with the two graph types kept distinct

The exact signed arrows 10,007→7,505 and 10,009→7,507 have literal same-order single-spoke witnesses. The replay supplies their deletions and apex contacts. The first one flips the arc to vertex 5, the second to 6. A Hamiltonian frame avoiding that spoke is common to both graphs. On such a frame, nonpath complementation commutes with flipping the spoke. This gives a real square with two ordinary arc edits and two genuine metalines. For 10,007→7,505 the opposite side has counts 9,849→7,355, and both metalines are black. From these **same** 10,007 and 10,009 source graphs, other frames give direct black lines to counts 7,505 and 7,507. At 10,007, the metaline lands in the very same merged class as the one-spoke edit; at 10,009, the two routes land in different merged classes despite the same count. This is a direct check of why scalar equality is not enough to carry the next move. There are also genuine direct black lines 257→193, 259→97 and 263→197 using different representatives from their single-arc examples. For 30,013→11,255 a literal three-spoke path is 30,013→24,961→13,955→11,255. A different order-12 source reduces that to two spokes with either order strictly decreasing, via 22,797 or 14,321. Another order-12 source similarly gives 30,011→5,627 via 29,325 or 8,315, with the same two spoke indices. These paths are in the arc-edit graph, not asserted to be metagraph lines.

Four signed families with the same variable transitive width admit real one-arc edits at every positive width: the pairs of affine reads (16·2^w+1,12·2^w+1), (64·2^w+5,12·2^w+1), (64·2^w−17,48·2^w−13), and (128·2^w−23,96·2^w−17). They use the [one-back-edge reader](PAIRED_RECURSION.md); flipping one specified exterior core arc changes its two coefficients on **the same** incoming module, verified in the replay. This is stronger than independent endpoint counts on those affine rails.

## Replacing one failed large line by two actual black lines

A general exchange principle: up to taking converse, a merged line reverses the **frame arcs only**. Suppose P starts at x and ends at y with the endpoint arc x→y. Reversing P creates the backwards Hamiltonian cycle using x→y to close it. Choose a second frame by deleting an internal P-arc f from that backwards cycle. Reversing the second frame restores every ordinary arc except that it flips both f and x→y in the original. If x,y are consecutive in a homogeneous transitive block and f is disjoint from them, swapping those two labels removes the x→y flip. The net change is exactly f, attained by two true lines.

Apply it at every width w≥2 with deletion D=(2,8,3,69,47,79,151,63), contact mask I=128, and core arc f=(1,3). Label the transitive module 8,…,w+7 and the external apex b=w+8. The physical frame
```
P = (8,10,11,…,w+7,4,5,2,0,1,3,6,7,b,9)
```
uses the indicated core arc. The frame for the second line is the backwards rim minus f. The swap of 8,9 is valid because their module contacts are identical outside the pair and they are consecutive inside it. The signed selector is therefore realized by two true lines for the exact pair
```
64·2^w+59 → 12·2^w+11.
```
The lines are both black: on P the arc from its first vertex to its third disagrees with its reflected converse (8 beats the third vertex, whereas 7 loses to 9); on the second frame the corresponding comparison is the absent 1→2 and present 7→3. At widths 2,3,4 the literal paths of counts are 315→369→59, 571→711→107 and 1083→1299→203. The complete step shrinks even when its first line overshoots.

The two lines are **necessary** for this pair of graph classes at every w≥32. Every frame in the source uses at most two transitive-module runs, because only the apex can immediately precede a noninitial module run. Reversing its frame arcs creates a directed triangle from each disjoint triple of consecutive vertices within a run, while the third chord is unchanged. This gives at least ceil((w−4)/3) vertex-disjoint triangles in any line partner. In the actual target, every triangle meets one of just nine exterior vertices, so there are at most nine vertex-disjoint triangles, also after isomorphism or converse. At w≥32 this rules out a single line. The explicit two-line construction shows the exact distance is two, not merely that individual endpoint counts are attained.

The working principle is to allow a change of physical representative or a short path when the scalar selector names a parent that cannot make a direct move. A universal selector still requires a theorem supplying one of these routes at every height; the finite atlas and all-width examples do not supply that theorem.

Replays: `atlas.py` and `metalines.py` take the archived `n7_representatives.txt` and optionally `n8_classes.csv`; `physical_edges.py` verifies the larger literal witnesses, squares, and two-line controls.
