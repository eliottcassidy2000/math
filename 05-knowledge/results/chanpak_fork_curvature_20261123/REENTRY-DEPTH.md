# How much a repeatable fork must remember

**PROVED scoped / FINITE-EXACT.** This extends [the two-read physical fork](PHYSICAL-TOGGLES.md) and [its uniform graft obstruction](UNIFORM-AMPLIFIER.md). The fixed-order arc toggles still work as proved there. Here the operation is different: wrap a filler in new vertices and then use the resulting tournament as the next filler.

## The exterior matching gives the exact number of needed slots

For a nonempty filler \(M\), write \(P_k\) for ordered covers by \(k\) nonempty directed paths and \(Q_s\) for ordered covers by \(s\) paths with empties allowed, so \(Q_s=\sum_k\binom{s}{k}P_k\). Let \(B\) be the exterior vertices beaten by all of \(M\), and \(A\) those beating all of \(M\). Compute the maximum matching size \(\kappa\) in the bipartite graph of the exterior arcs \(B\to A\).

In an output \(s\)-cover, each gap between successive \(M\)-runs within one path starts in \(B\) and ends in \(A\). One \(B\to A\) edge from each gap gives a matching, since the gaps are vertex-disjoint. There are at most \(s\) paths containing a run. Thus there are at most \(s+\kappa\) runs. Conversely, to get any \(1\le k\le\min(|M|,s+\kappa)\), put one module placeholder in each of \(\min(s,k)\) paths and put \(k-\min(s,k)\) matching edges between consecutive placeholders in one path. Insert every unused exterior vertex into a path in turn. One can always prepend, append, or insert at the switch from beating the inserted vertex to losing to it; all existing vertices keep their relative order and the gaps survive. Fill the placeholders with any ordered nonempty \(k\)-cover of \(M\). This proves that the \(s\)-read is a linear combination of \(P_1,\ldots,P_{\min(|M|,s+\kappa)}\) with **positive coefficient on every index**. Those coefficients depend only on the exterior skeletons, not on the filler.

This budget adds exactly under nested uniform wraps, each uniform to the whole previous tournament. When a new layer is added, its incoming vertices beat *all* earlier vertices, while its outgoing vertices lose to all earlier vertices. Thus the \(B\to A\) matching contains no edges between different layers; the final \(\kappa\) is the sum of the individual matching sizes. The linear budget is sharp for unbounded filler order: \(P_k\) on transitive tournaments is \(k!\left\{\begin{smallmatrix}n\\k\end{smallmatrix}\right\}\), and the evaluations at \(n=1,\ldots,L\) give a triangular full-rank test for the first \(L\) coordinates. This proves linear sharpness, not a general nonlinear lower bound at every depth.

## A nonlinear no-go, with actual fork twins

The following strong nine-vertex fillers, given by outgoing row masks, have identical \(D=H=2231\) and \(Q_2=36088\):
\[
\begin{split}
M_A&=(42,124,345,432,481,68,137,295,99),\\
M_B&=(266,412,121,464,417,203,403,261,36).
\end{split}
\]
They give the same factor-seven fork \((30356309,30360771,30365233)\). Their exact nonempty-slot difference \(P(M_B)-P(M_A)\) is \((0,0,240,480,240,0,0,0,0)\).

- If \(\kappa=0\), the extension is ordinal \(A\to M\to B\); \(D,Q_2\) update by separate multiplication. The pair is universally sufficient.
- If \(\kappa=1\), the next \(H\) depends only on the pair, but the next \(Q_2\) depends positively on \(P_3\). These twins prove no deterministic rule, linear or not, can update the pair universally. The coefficient of \(P_3\) is at least two: use a bridge in one slot and a lone placeholder in the other, or exchange slots.
- If \(\kappa\ge2\), \(P_3\) also has positive coefficient in the next ordinary \(H\), at least two by taking the two bridges in either order. The twins separate both reads. Differences from higher slots are nonnegative and cannot cancel.

This tests the exterior first, before counting any paths. Adding more uniform vertices cannot remove an existing matching edge as long as the original arcs and contacts are retained.

## The cyclic sandwich is an exact slot operator

Write \(W(M)\) for \(M\to b\to c\to M\). Its matching size is one. If \(b,c\) occupy different slots, then \(b\) is at its slot's end and \(c\) at the beginning; contribution \(s(s-1)Q_s\). In the same slot, if \(b\) comes first it immediately precedes \(c\), and separates that slot into two independently allowed pieces, giving \(sQ_{s+1}\). If \(c\) comes first the two vertices enclose a nonempty module run, giving \(s(Q_s-Q_{s-1})\). Hence
\[
Q_s(W(M))=sQ_{s+1}(M)+s^2Q_s(M)-sQ_{s-1}(M).
\]
This is an actual reusable transition on a slot polynomial, not a new physical ternary center from the two original integers alone. After \(d\) sandwiches the top input \(Q_{s+d}\) has coefficient \(s(s+1)\cdots(s+d-1)\) (the unique sequence of upward shifts). Thus reading \(H(W^d(M))\) for \(0\le d\le r-1\) is an integer-triangular observation of \(Q_1,\ldots,Q_r\), with diagonal \(d!\).

For the twins the first sandwich gives \((38319,743964)\) and \((38319,744444)\). The second ordinary counts are \(782283\) and \(782763\). These are exact separations of 480, independently recounted by subset DP in the replay.

## A finite modular controller, and what the two triples actually distinguish

For any modulus \(m\ge2\), put \(K=\max\{k:m\nmid k!\}\). Permuting the \(k\) nonempty slots freely gives \(k!\mid P_k\). Thus modulo \(m\) only \(P_1,\ldots,P_K\) matter in **every** uniform continuation, for every output slot count. The binomial-basis expansion remains valid on negative integers. Using the slot-polynomial parity \(Q(-s)=(-1)^{|M|}Q(s)\), the readings \(Q_1,\ldots,Q_{\lceil K/2\rceil}\), plus known order parity, supply \(K+1\) consecutive evaluations including \(Q_0=0\). Newton finite differences recover the truncated coefficients modulo \(m\), using integer addition only and no division (binomial basis shifted to start at the first evaluation point). Thus the same number of readings plus the updated parity suffices at **every** step, for arbitrary successive uniform exteriors, even while exact matching depth grows.

For moduli \(5,6,7,8,9\) this register bound is respectively \(2,1,3,2,3\). In particular, at a known order parity:

| Modulus | Do two reads determine the next two in every uniform graft? |
|---|---|
| 5, 6, 8 | Yes, and the conclusion iterates. |
| 7, 9 | No; three suffice and the twins disprove two already on \(W\). |

The distinction is independent of the [spherical/flat/hyperbolic comparison](README.md). The exact same two readings and two exit directions at the initial fork conceal distinct successor residues mod 7 and 9 but not mod 5, 6 and 8. At 5 the actual repeatable sandwich is particularly small: setting \(\sigma=(-1)^{|M|}\), the two-vector is multiplied by \(\left(\begin{smallmatrix}1&1\\-2&4+2\sigma\end{smallmatrix}\right)\). For even/odd order its vector orders are 24/12 and its projective orders are 6/3, respectively. These finite periods are for the sandwich map, not the projective generators of the triangle group.


## Exact observation lattice: the bridge permutation gain

The factorial cutoff can be sharpened when only the first few *output* reads matter. Fix \(s\ge1\). For the integer lattice generated by every uniform exterior's first \(s\) output reads, the exact generators are
\[
W_k^{(s)}=\frac{(k-1)!}{(\min(k,s)-1)!}\,P_k.
\]
Here “lattice” means taking integer combinations of readouts, not realizing every combination as one exterior.

Indeed, in a skeleton with \(k\) runs distributed across \(r\le\min(k,s)\) occupied slots, fix the occupied positions, each slot's exterior prefix and suffix, the other exterior-only slots, and the unordered set of its \(k-r\) distinct bridge words. All orderings and weak distributions of these bridges among the \(r\) slots are valid, because uniformity makes every bridge attach to any module segment. The number of skeletons with these data is \((k-r)!\binom{k-1}{r-1}=(k-1)!/(r-1)!\). Thus all coefficients are divisible by the asserted weight. Conversely, for \(k\le s\) use bare \(Q_k\) in triangular order. For \(k>s\) use \(Q_s\) of the [ordinal pair ladder](../chanpak_cofactor_contradiction_20261123/UNIFORM-GUARDS.md) with \(k-s\) matching pairs. At the top all \(s\) slots meet \(M\) and all exterior vertices are used by matching bridges. The only possible matching is the identity (the matched second indices must be at least the first and have the same sum), so bridge redistribution gives exactly the asserted coefficient. The preceding terms can be eliminated over the integers.

In particular **\(s=1\) and \(s=2\) have the same weights \((k-1)!P_k\)**. Two fillers agree on every future \(H\) modulo \(m\) if and only if they agree on every future pair \((H,Q_2)\) modulo \(m\), equivalently if their weighted vectors agree. The state is finite: \(W_k^{(s)}\) vanishes automatically once \(m\) divides its weight times \(k!\). It is also stable under continued uniform grafts, since composing a later probe with the current exterior is itself a uniform probe of the original module. To update it explicitly, form the finitely many triangular probe reads around the new filler, express each composite probe as a linear combination of the old weighted reads, and eliminate again. This stronger state needs no parity metadata. In particular, without knowing parity, the transitive three- and four-vertex fillers have equal \(H=1\), \(Q_2=0\bmod8\), but one sandwich separates their next \(Q_2\). Thus the known-parity qualification in the preceding table matters.

## A restricted exact renewal that does work

[The two wrapper towers](RESTRICTED-RENEWAL.md) provide a positive alternative to enlarging the state: two integers recover the marked strong nucleus and the number of source pads. Repeatedly wrap the nucleus or add a source, retaining the original fork exterior. This is closed on a specified family, not on all fillers. The same companion gives a minimum four-vertex exterior whose two matching edges make the twins’ Hamilton counts diverge in one step.

Replay: `python3 reentry-depth-check.py` and `python3 matching-restricted.py`.
