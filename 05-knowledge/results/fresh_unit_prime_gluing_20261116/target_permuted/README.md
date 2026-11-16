# Choosing a fresh continuation from an absolute target

**Exact handoff.** At height 5 of the [old doubling row](../ternary_deposit/README.md), the previous +4 child has 8053 Hamiltonian paths. It can be re-encoded by a rooted tournament with $(H,E,O,g)=(8053,6479,6484,5)$. Its next exact digit gives **24159, 24161, or 24163** with the same continuation **(19437,19452,15)**. All three new outputs are strong. This is a different source graph at the same scalar, not an induced extension rule. At height 5 it actually triples the new source's continuation, unlike the earlier [997 handoff](../changed_continuation_dyadic/HEIGHT_TWO.md).

Among heights 1 through 8 the three prescribed counts share a continuation at $g=15$ exactly at heights **1, 2, 5**, even allowing arbitrary rerooting, nonstrong hosts, and $E=O$. This is specific to deletion count 15; a missing height might use a different $g$. The [exact scan and literal rows](ROW_SCAN.md) spell out the scope.

## Exact inverse for a terminal slot

A general identity explains these infinite bands without extrapolating finite numerical fits. Take a rooted tournament with a marked vertex whose sole outward neighbor is the root $r$. Delete that vertex and $r$ to obtain $L$. Were the mark replaced by a nonempty tournament $M$, the root deletion would be precisely the ordinal sum $L\to M$, and the root would lose to all of $M$. Let $u=H(M)$ and $R=P_M(2)$ (the count of two ordered path slots, each individually allowed empty).

On $L$, let $A_e,A_o$ count weighted old paths starting in an outneighbor of the root, where a path is weighted by $H$ of its complement in $L$. Let $J$ be the even part restricted to paths ending in a vertex beating the root. Write $F=A_e+A_o$, $g_0=H(L)$. Then the four coordinates are
$$
H=FR+Du,\quad E=FR/2,\quad O=FR/2+g_0u,\quad g=g_0u,
\qquad D=g_0+2J-2A_e.
$$
The even and odd optional-slot totals of any nonempty $M$ equal $R/2$ by odd-cycle packing cancellation. Each nonsingleton root-start path begins in $L$ and sees those balanced choices for its optional tail; the singleton root contributes $g_0u$ to $O$. For $H$, the [rooted odd-cycle formula](https://github.com/eliottcassidy2000/math/blob/c3cf550db9ba849e8d1dfb9cf2c6eefc61b550e5/05-knowledge/results/triangular_backbone_chanpak_20261116/CORE.md) adds twice the weighted even old paths from an outcontact to an incontact. The paths stopping in $L$ contribute $uJ$; the ones continuing into $M$ contribute $A_e(R/2-u)+A_oR/2$. This gives the identity. Equivalently, from the singleton case $(H_1,E_1,O_1,g_1)$ the signature is $(F,D,g_0)=(E_1,H_1-2E_1,g_1)$.

**Choosing the continuation from the target.** For $F>0$ and a transitive tail of length $t\ge1$, $u=1,R=2^t$. For given $N$, each *physical* signature is tested by dividing $N-D$ by $F$: it must be a power of two at least two. On success the test gives both the tail length and the exact continuation $((N-D)/2,(N-D)/2+g_0,g_0)$.

For a three-digit target window, use three signatures $(F_i,D_i,g_0)$, $i=0,1,2$. They give $N,N+2,N+4$ with the same continuation **if and only if**:
1. $D_i=D_0+2i$;
2. all $F_i$ and the positive number $N-D_0$ have the same odd part;
3. $v_2(N-D_0)\ge1+\max_i v_2(F_i)$.

Necessity is equality of $F_i2^{t_i}=2E$; sufficiency uses $t_i=v_2(N-D_0)-v_2(F_i)$. This is an exact target selector once there is a physical signature bank, not a theorem that the bank covers every target. It compresses an infinite range of tail lengths without dropping an attainable continuation.

At height 5 the output signatures are $(19437,-14715+d,15)$ for $d=0,2,4$. The independently found input signature is $(6479,-4905,5)$. The arbitrary-module formula therefore extends the tripling and its three digits to all transitive tail sizes, with $H_{
m in}=6479\cdot2^t-4905$. That new family lines up with the desired old-row target only at the base size: at $t\ge2$ its lowest output is 1 modulo 4 while the old lowest output is always 3 modulo 4. This same valuation obstruction applies to a fixed band of any odd $F$ lined up at its base size.

## Permutohedral moves: change $g$ and keep the exact target

The [transitive-flip theorem](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-207-permutohedron-transitive.md) identifies transitive orientations with the vertices of the ordinary permutohedron. More generally the convex hull of labeled score vectors is that permutohedron: each pair contributes either of its coordinate unit vectors, and generic linear maximization orders the coefficients and yields a transitive orientation. A score vector can have many orientation lifts. Lifts with the same labeled scores are connected via directed-triangle reversals. Indeed the disagreement arcs are Eulerian, hence split into directed cycles; each cycle flip decomposes into directed triangle flips by splitting at a chord and ordering its two sub-flips according to the initial chord direction. Thus score preservation is exact, but it does not force $H$ or the continuation.

Here are two explicit same-score recontacts. Rows are outgoing masks, slot 0 is marked, and the root stays fixed. The indicated triples are reversed in order, each a directed triangle at that step.

| Base source rows | Root | Triangle sequence | Endpoints $(H,E,O,g)$ |
|---|---:|---|---|
| 4,17,82,39,105,135,171,31 | 2 | (2,3,4), (2,5,6) | (283,158,187,29) to (283,158,191,33) |
| 8,21,73,194,45,15,179,55 | 3 | (1,2,3) | (309,166,193,27) to (309,166,189,23) |

No reversal touches the slot, so the same moves work after replacing it by any common nonempty module. By the terminal identity the two endpoints preserve $H=158R-33u$ or $166R-23u$ respectively, together with $E$, all labeled scores, and all contacts to the module. The defect $H-2E+g$ changes from $-4u$ or $+4u$ to zero. These are genuine count-preserving endpoint re-encodings at unbounded height. The negative-to-zero example temporarily changes its base count to 291. The positive-to-zero example uses just one directed-triangle reversal, so no intermediate count is spent; its output masks are 8,25,67,196,45,15,179,55. Also zero defect is only a ready *equation*, not proof these particular graphs admit all three fresh digits. Run **python3 check_permutation.py** for the direct DP of both paths, several substitutions, and all four literal height-five profiles.

A permutohedral swap of adjacent singleton contact marks $01\leftrightarrow10$ by itself is not a neutral move. Within the same ordinal singleton gap, the two unsigned-start weights are $2w,w$ for a positive $w$ determined by the preceding component counts and following slot factors. Passing the 1 earlier increases the weighted start sum by $w$. A coordinated swap must compensate that change or intentionally change the continuation. This is why the correct search groups **both** exact target counts and rooted coordinates, not just inversion distance or score.

## Target-pruned search and scope

**selector.cpp** reads an absolute lowest count $N$, a maximum deletion order, and a number of profiles at which to stop. It uses a bank of strong components with product 15, constructs exact weighted even-path cut matrices by parity convolution along the ordinal list, and searches binary root contacts. As soon as both bits of a weighted cut edge are assigned, it charges its nonnegative contribution. A partial cut over $(N+4-15)/2$ is therefore pruned without losing any target; no accumulated count or residue is used as a heuristic. At a surviving exact target it computes $(E,O,15)$ and intersects the three support sets. For example:
```sh
c++ -std=c++20 -O3 selector.cpp -o /tmp/target_selector
/tmp/target_selector 24159 18 1
```
The first positive strong triple at this target has E=19437 and O=19452. The complete all-host-orders scope at $g=15$ for the first eight heights uses the separate [row replay](ROW_SCAN.md), including prefixes and suffixes. Neither search assumes that the source is retained or claims universality across $g$.

The complementary [1729 factor decoder](TAXICAB.md) works on the two-core side instead of the terminal tail. It provides another exact way for a target to prescribe candidate continuation data, and a proof that symmetric cores alone cannot solve the task.

## Follow-up: the fourth height works by changing the deletion count

At $g=45$, the fourth old-row height reopens: 4021 goes to 12063,12065,12067 with continuation (8943,8988,45). The [follow-up](NEW_READERS.md) also gives an all-height alignment on the old +2 branch for the restricted singleton reader, and exact +2 moves using equal-weight contact swaps in one deletion graph.
