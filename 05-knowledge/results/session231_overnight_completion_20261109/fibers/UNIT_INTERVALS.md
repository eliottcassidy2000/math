# Unbounded consecutive-odd repair intervals from one saturated moving core

Status: PROVED, independently audited. Superseded in efficiency, not validity, by SHORT_INTERVALS.md. This does NOT settle unrestricted Chan–Pak completeness.

## Main theorem

For every integer k≥0, there is a single physical rooted full interface supporting k+1 consecutive odd Hamilton-path counts. Every displayed witness is a strong tournament on1188k+12 vertices. Only ONE nontrivial deletion core is used; no repeated-core scalar divisor develops.

The physical seven-vertex core is the graph `adj` in `UNIT_LOOP_WITNESS.json`, cut30. In canonical upper-triangle encoding it is1147012. Its state (h,q,A,G,T,S,E) is

    (37,740,269,422,132,12,21).

It is saturated: h divides q. This is essential to the unit-drift escape; the unsaturated cycle-quantization theorem would forbid it for h>1.

Let X,Y be the length1188 binary words in `UNIT_LOOP_WITNESS.json`. They end/start with X_last=0,Y_first=0. Form

    L_(k,j)=s1 ⇒ s0 ⇒ X^j ⇒ C ⇒ Y^(k−j) ⇒ s0,  0≤j≤k.

A bit denotes a singleton in the deletion graph, selected by the root exactly when the bit is1. All arcs between consecutive listed blocks point forward. Add v,r with L→v→r, and r beating precisely the selected vertices of L. The leading selected singleton gives a path from r to every later vertex; every deletion vertex reaches r through v. Thus every closure is strong.

There is an odd baseline D_(k,0) such that all displayed counts are

    D_(k,0), D_(k,0)+2, ..., D_(k,0)+2k.

The core order is7, there are1188k+3 singletons, and closure adds2 vertices, giving1188k+12 total.

## Common physical interface

Set b=2^1188, R=value(Y), V=value(reverse(Y)), and J_k=(b^k−1)/(b−1). All choices j share

    h=37,
    q=740·2^(1188k+3),
    A=1749·2^(1188k+1)+74R J_k,
    G=1725+2960V J_k,
    S=37, E=0.

These are actual path-count coordinates of the displayed deletion graphs, not formal admissible numbers. The baseline is specified exactly by j=0 and the ordinal state recurrence. For k=0 the closure count is6057. At larger k its binary length is1188k+13 in the checked range; this finite observation is not needed for the theorem.

## Finite identity and all-order proof

The literal core enumerator gives the stated seven-coordinate state. Applying the exact ordinal composition law to the fixed words verifies

    s0 ⇒ C ⇒ Y ⇒ s0,
    s0 ⇒ X ⇒ C ⇒ s0

have equal (h,q,A,G,S,E), and the second has T exactly2 smaller than the first. This is a finite identity with exact integers, not a conjecture based on increasing k.

For blocks differing only in their (T,S,E) coordinates, but agreeing in h,q,A,G, surrounding singleton contexts change the T difference to

    ΔT + last_bit(left) ΔS + first_bit(right) ΔE.

This follows immediately by composing on the left and right. For each movement from position j to j+1 in L_(k,j), the left context ends0 and the right context begins0: the words X,Y and the external zero singletons ensure this. Hence EVERY move changes T by−2, regardless of k or j. Therefore

    T_(k,j)=T_(k,0)−2j.

The closure count is2A+h−T, so it increases by2. The common-interface formulas are obtained at j=0; each movement preserves them. This proves the all-order theorem.

## How the unit word is constructed

The enhanced carry state is (a,g,p), p=previous y, with

    x=g mod2, a'=2a+hy−qx, g'=(g+qy−hx)/2,
    w=a y+h p y−x g'−h x(g' mod2).

`UNIT_LOOP_SEEDS.json` contains TWO actual directed closed walks of length108 through this same physical core interface. Their gains are+34 and−90. Concatenate eight positive walks and three negative walks. They remain a closed walk through the same physical interface, and

    8·34+3·(−90)=2.

Their concatenated x/y words are X,Y; length is11·108=1188. This use of signed cycle gains is legitimate because both walks are themselves forward, physically executable carry walks. No edge is reversed and no negative number of repetitions is used. The positive and negative labels name their scalar gain, not traversal direction.

The discovery code originally supplied several cycles with gcd-two gains. Gcd alone would not suffice. What closes the construction is the displayed NONNEGATIVE repetition identity, together with the common physical basepoint and compatible endpoint bits. This is precisely the missing physical sidecar in a bare gcd argument.

## Verification and scope

`verify_unit_family.py` depends only on the physical core, fixed word witness, and ordinal composition routine. It checks the core by independent literal permutation-and-split enumeration; every edge and the exact closure of the long carry word; its gain2; the finite flanked word identity; and all positions through k8. `UNIT_FAMILY_CHECK.json` and its optimized counterpart agree. The proof above supplies every k. A separate seeded-word checker is included for the+34/−90 decomposition.

This yields arbitrarily long intervals of consecutive odd counts located at particular exponentially growing baselines. It does not show that these intervals overlap as k changes or meet an arbitrary requested odd target. In particular it does not prove Chan–Pak completeness, remove the exceptional counts7,21, or produce a global occupancy theorem. Its concrete advance is unbounded consecutive repair width in one physical full-interface fiber, using one moving saturated core and no growing repeated-core divisor.

The next mathematical obligation is to connect these real intervals to an effective baseline-selection/overlap mechanism. A smaller unit loop would improve the width-to-order ratio, but is not needed for the existence theorem.
