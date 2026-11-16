# Changing the deletion layer and the reader

## An exact fourth-height handoff at a different deletion count

The [deletion-15 scan](ROW_SCAN.md) misses the fourth height. Changing $g$ to 45 works: the old +4 branch of count **4021** can be replaced at its next targets **12063, 12065, 12067**, with exactly common continuation $(8943,8988,45)$. Unlike the fifth-height example this is not asserted to triple a known physical 4021 source's continuation.

Use the base rows below with the last vertex the root. Into each insert a nonempty transitive block between the root deletion and the root: the deletion beats all of the block, which beats the root, and all old arcs stay unchanged. At block length one the outputs in order are the three displayed counts. Their profiles for block length $b\ge1$ are
$$
(8943\cdot2^b-5823+d,\ 8943\cdot2^{b-1},\ 8943\cdot2^{b-1}+45,\ 45),\qquad d=0,2,4 .
$$

| Digit | Row masks before adding block |
|---:|---|
| 0 | 2046,4092,2040,4080,1376,3776,1408,3856,3616,3152,0,1109 |
| 2 | 2046,4092,2040,4080,1248,1216,3968,3840,3632,3120,0,1077 |
| 4 | 2046,4092,2040,4080,1248,3264,3968,1792,3632,3120,0,1173 |

This is another exact single-height alignment, not a lift at all heights of the old row.

## The old +2 branch has all-height alignment for a restricted reader

For $h_s=252\cdot2^s-13$, $s\ge0$, the following base rows followed by a transitive block of length $s+2$ give $3h_s,3h_s+2,3h_s+4$. All three have common $E=189\cdot2^{s+1}$. Their deletion counts are **37, 37, 49**, so their full continuations do not all match.

| Digit | Base row masks |
|---:|---|
| 0 | 958,1020,952,944,992,448,909,768,0,288 |
| 2 | 30,108,24,48,98,69,13 |
| 4 | 190,252,152,240,32,68,21,112 |

The common $E$ suffices exactly for every future word of **retained-root singleton reverse bridges**. The [bridge formula](../BRIDGE_RESET.md) is $(H,E)\mapsto(H+2E,2E)$ in the normalized sector $O=E+g$, with $g$ unchanged. Thus the chosen counts remain consecutive through any number of those doublings, with the same new $E$. That realizes synchronization relative to a specified continuation language, not a new fresh digit or synchronization under arbitrary right readers. It is an all-height alignment on a weaker interface.

## Equal-weight contact swaps at a fixed deletion

In the [terminal decoder](README.md), put $K=A_e-J$, the weighted even-vertex paths with **both** endpoints in the old outcontacts. Then the intercept is $D=g_0-2K$. A swap of one outcontact for an incontact of equal weighted start total preserves the slope $F$; changing $K$ by one gives exactly a +2 or −2 while preserving the full output continuation.

For a literal example use deletion rows $(30,12,56,16,34,11)$, of count 29. The weighted start totals by vertex are $(86,43,50,31,31,50)$. Use old root outcontact sets $\{3,5\},\{4,5\},\{2,4\}$ in order. The two swaps keep the slope 81 and change $K=8,7,6$. Therefore after the terminal block the counts are $81\cdot2^b+13,81\cdot2^b+15,81\cdot2^b+17$ at common $(81\cdot2^{b-1},81\cdot2^{b-1}+29,29)$. At $b=1$ these are 175,177,179 around $3\cdot59$. Unlike swapping labels on an old transitive tournament, these moves change actual root arcs. The nonlinear endpoint sum, not the score, determines the repair.

A useful bound makes this decoder finite at a fixed intercept and specified deletion factor $g_0$. In a strong realization the root contacts the first deletion component. Strip ordinal singleton components after the last old outcontact: each doubles $F$ and leaves $K,g_0$ unchanged. Among $J_0\ge1$ singleton components strictly between the earliest and latest old outcontacts, choose any even subset for the directed path between those endpoints. Those paths contribute at least $2^{J_0-1}$ to $K$, since their omitted path counts are at least one. So there are at most $3+\lfloor\log_2 K\rfloor$ singleton components in the reduced search when $K>0$ (at most the two endpoints when $K=0$). The nonsingleton components are already bounded at fixed $g_0$ by their count factors and Busch's bound. Thus a fixed $(g_0,D)$ has a finite bank of slopes up to powers of two, which tests *unbounded* transitive block heights. This is an alternative to raising the host-order ceiling at every new height.

Run **python3 check_swaps.py** for the literal profiles and suffix controls. It also verifies an independent count-1729 tournament: a transitive deletion with root outcontact word 11101001000 gives $(H,E,O,g)=(1729,932,933,1)$. That is different from both the Paley orbit quotient and a normalized mirror bridge.
