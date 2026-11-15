# A two-contact switch which preserves both tail interfaces

The [previous recontact](../chanpak_carry_admission_20261115/OUTLET_WIDTH.md) bought a genuine +2 but did not preserve the full continuation state. A different, smaller switch does both. It also preserves the *incoming* reader: one may grow arbitrary painted ordinal blocks to the left as well as the right. A separate direct witness gives three consecutive counts with this stronger interface.

## 1. Same spine, two changed arcs, exactly two new Hamilton paths

Give the five-vertex spine B outgoing rows (6,24,18,5,9), indexing bits from zero. Attach one apex q with outgoing bits 5 (vertices 0 and 2) or 3 (vertices 0 and 1). Here h=H(B), Q=Q_2(B) counts ordered two path slots including empties, and the marked boundary data are:

- s and e: total and odd-minus-even weighted starts from vertices beaten by q;
- t and f: total and odd-minus-even weighted ends into vertices beating q.

Each nonempty path is weighted by the Hamilton count of its complement. The two masks have the exact common signature
```
(h,Q,s,e,t,f) = (15,120,42,6,63,9).
```
Their apex counts N are 41 and 43, respectively. Both resulting six-vertex tournaments are strong, since B is strong and the apex has both kinds of contact.

The preservation mechanism is not a chance equality. Relabel B in cyclic order (0,2,1,4,3); from each vertex its two successors in this order are its outneighbors. Thus B is vertex-transitive, so every vertex has exactly the same four singleton weighted endpoint statistics. All masks of the same cardinality therefore preserve (s,e,t,f). In general a vertex-transitive tournament of order n gives for a k-mark mask
```
s=k(Q−h)/n, e=k h/n, t=(n−k)(Q−h)/n, f=(n−k)h/n .
```
The unsigned start sum over vertices is Q−h, deleting the empty first path; the signed odd-minus-even sum is h, since the total signed sum including the empty part vanishes by the odd-packing expansion, with the empty first path removed. End weights obey the same argument after taking converse. These are common boundary statistics, not enough to recover N; the replay proves the exact differing counts by literal subset-endpoint DP. One changed apex arc on a fixed spine cannot preserve s, since each vertex's weighted start sum is positive (its singleton path alone contributes). Thus two recontacts are minimal on a fixed spine.

## 2. Why any ordinal attachment preserves the interface

Use the [four-coordinate cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md). Appending a painted block with local data (h,Q,N,s,e,t,f) to a state (a,A,S,d) produces the count and state
```
a' = a h
A' = h A + (N−h)a + t S + f d
S' = Q S + a s
d' = a e.
```
This follows by splitting the single odd cycle containing the apex at the ordinal seam and classifying its two path parts by parity; other odd cycles are left or right and their weights factor. In the start update, the possible right continuation of a path starting on the left doubles into Q ways, and its signed sum vanishes because every cycle and singleton in the [odd-packing expansion](https://arxiv.org/abs/2307.05569) has odd order.

Consequently a local shift Δ between blocks of equal full signature becomes aΔ with the same outgoing (a',S',d'). Every further common ordinal block multiplies the difference by its ordinary Hamilton count without otherwise distinguishing the state. In particular arbitrarily many painted **transitive** vertices on either side preserve the absolute gap +2. General cyclic blocks remain allowed, but scale the gap by the product of their ordinary Hamilton counts. This does not assert preservation under arbitrary nonordinal consumers, or allow iterating +2 on the same marked spine to obtain an unlimited ladder (the two-state switch reverses itself).

## 3. Explicit left/right transfer: two tails interact by a product

Let the left transitive word have ordinary binary value x, with 1 meaning the apex beats that vertex. Let the right transitive word use the canonical reverse-complement address y from the [binary suffix recurrence](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md): its least significant address bit is 1 iff the right word's first vertex beats the apex. Allow empty tails, and strip neutral leading zeros on the left or trailing ones on the right. Write p(z) for the common binary apex correction
```
p(0)=0, p(2u+1)=p(u), p(2u)=p(u)+u+(u mod 2).
```
The exact two-sided output of the switch is
```
N_i(x,y) = 41+2i + 15[p(x)+p(y)] +63x+42y+120xy
                    +9(x mod 2)+6(y mod 2),   i=0,1.
```
Indeed the left word gives (1,1+p(x),x,x mod 2); applying the displayed transfer puts the marked start at 120x+42 and bias at 6. The right word has unmarked end total y and end bias y mod 2, since the even/odd sums of paths to any later transitive vertex cancel. Applying it gives the formula. The exact product 120xy expresses interaction of two bent tails through the same core, rather than two independent additive perturbations. No carry-free condition is used. For nonempty canonical tails the first left mark is 1 and last right mark is 0, so the enlarged graph is strong.

## 4. Three consecutive values with a matched two-sided interface

Here are strong six-vertex spines with one marked apex. Row and mask conventions are the same.

| Rows | Mask | Apex count N |
| --- | ---: | ---: |
| (6,28,40,17,37,11) | 6 | 135 |
| (26,44,17,36,10,21) | 10 | 137 |
| (26,44,17,36,10,21) | 6 | 139 |

Their entire common signature is (45,432,129,15,258,30). The last two witnesses even use the same spine and switch only two apex arcs. For arbitrary transitive left/right tails the whole window is exactly
```
135+{0,2,4} +45[p(x)+p(y)] +258x+129y+432xy
                       +30(x mod 2)+15(y mod 2).
```
As an order-nine literal check, set x=1 and y=1 (one newly appended vertex at each side). Each six-vertex spine plus two tails and the apex then has exactly nine vertices, and their counts are 999,1001,1003. This is a three-count window at order nine, not the previous nine-count window in a single common-interface fiber. The new result changes what a window can survive, not yet its maximal width or where it can be placed.

Replay: `python replay_two_sided_recontact.py` checks the literal profiles, the rotation, and independent whole-tournament path counts for both two-tail formulas through a grid of addresses including the order-nine triple.
