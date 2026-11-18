# Separate the two axes: an exact paired carry and a minimum-seed theorem

The [paired recursion](PAIRED_RECURSION.md) transport rule is linear on (h,e,t,d), where h is the ordinary path count, e is the end-at-root count, t is the distinguished-start two-path count, and d is the ordinary count after deleting the root. It preserves exact ×7 between two physical states. This note finds a constraint explaining why this architecture **cannot start at the usual least seeds**, and gives one way to use its freedom to achieve an exact ±2 carry. It does not establish universal coverage.

## The holes are inherited by the reservoir, not just by the final count

If (h,e,t,d) and 7(h,e,t,d) are both rooted tournament states, then d and 7d are actual ordinary path counts. Hence **d cannot be 1 or 3**, even when h is arbitrarily large. Since d≤h, uniform four-coordinate pairing cannot start at h=1 or h=3 either. This is not a modular effect. A purported root reset to a transitive deletion would require a real tournament with seven paths on the other side.

More strongly, **every paired four-coordinate state with lower h=5 is divisible by 5 in all its coordinates**. In particular a homogeneous paired recursion cannot bootstrap from 5 to a different prime.

Proof. A Hamiltonian path of T−a can always be extended by inserting a. If a lies in a nontrivial strong component, take a Hamiltonian cycle in that component and delete a; its remaining path permits insertion at both ends, yielding strict inequality for the component's counts and therefore for T. Thus d=h exactly when the root is a singleton strong component. The paired state at h=5 has d≥5 by the two forbidden totals, so d=5, and both roots are singleton components.

Counts multiply along the ordered strong components. A count-five tournament has exactly one nontrivial component; the same is true at 35 because factoring 35 nontrivially would use the forbidden value 7. By [Busch's strong-tournament bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) their nontrivial components have at most seven vertices. The archived order-seven class representatives, including all induced restrictions, give the respective (H,B) reads (5,32) and (35,676); the replay verifies these exact finite inputs. A singleton root before its unique core has t=B_core times a power of two, because the core can divide between both slots; a singleton root after the core has t=H_core times a power of two, because the core must be in the optional slot. The odd parts of t on the lower side are thus 1 or 5; on the upper they are 169 or 35. Exact ×7 forces the after-core choice on both sides. So t is divisible by 5, as are h=d=5 and e∈{0,5}. This proves the stronger boundary.

The same census checks the basic holes: a strong tournament with count 7 or 21 would have at most six vertices, where neither occurs. Factoring into distinct ordered strong components cannot produce 7, and cannot produce 21 without a factor 7. These exclusions are separate from any proposed recursive completeness argument.

**Design consequence.** Keep an ordinary lane for h=1,3 and for useful deletion counts 1,3. Do not require proportional rooted states every time a scalar pair exists. In particular 5 and 35 are both actual counts, but pairing them in all four reads forces a nonprimitive seed. A construction that changes the reader, changes the root signature, or independently carries the ordinary and upper witnesses can bypass that constraint.

## An exact two-sided carry that works in actual graphs

The paired attachment rule has a free extra letter if keeping the root: append an ordinary singleton without reversing any edge. Call this 0, and call the previously used singleton with its sole backward arc to the root 1. Both retain d, while neither path heights nor the other reads are guessed:
```
0: (h,e,t,d) → (h,0,2t,d)
1: (h,e,t,d) → (h+t−e,d,2t−d,d).
```
All four formulas are physical applications of the ordinal join or the one-arc seam. A transitive block with a single back arc from its root is exactly a word of these letters with one 1, so varying the boundaries of such blocks cannot give a new mechanism at fixed vertex word.

Take the two choices at a first bit, then apply **the same suffix** of bits b_i, indexed from 0 in physical order. Put v=sum b_i 2^i, epsilon=b_0, and
```
K(v) = sum over i<j with b_i=b_j=1 of (2^(j−i−1)+[j=i+1]).
```
The lower count after first bit 0 is h+2tv−dK(v). Along the suffix the initial difference in t is −d and doubles each step, and the initial difference in e is d and is erased after one suffix letter. Thus
```
H_after_initial_1 − H_after_initial_0 = t−e−d(v+epsilon).
```
This is an **exact carry chooser**, not a limit argument. To prescribe both an endpoint p and a gap delta, a proposed d,e and suffix specify t=e+delta+d(v+epsilon) and h=p−2tv+dK(v). What remains is to realize those seed readings. Conversely from an available seed the gap fixes v+epsilon by one exact division.

Physical example: start from the paired seed (29,11,31,13) and its 7-multiple. Prepend 0001, getting (277,13,483,13) and its 7-multiple. The request delta=2 forces v+epsilon=36. Choose v=35; its marked bits are 0,1,5 and K=26. The endpoints are **33,749 and 33,751**, both prime by trial division; the upper counts are **236,243 and 236,257**. Both pairs use the same suffix after the variable first bit. Independent subset DP checks all four order-19 tournaments. The useful new feature is that the difference is chosen algebraically, and an actual previously available seed supplies the center. The two continuations need not be identical after this carry.

## What an all-target descent would need

The unsigned ternary maps h↦3h,3h+2,3h+4 give a unique-parent tree on positive odd integers, ignoring the loop at 1. Deleting 7,21 leaves a forest with **exactly six roots**: 1,23,25,63,65,67. The other five roots are surviving children of 7 and 21. A physical readiness invariant preserved by those three moves would therefore need to initialize just these roots, rather than a long interval. This statement is conditional on the physical moves at all heights, not a proof that they exist. But the first surviving child 23 already has a physical fork: three different outdegree-one rooted witnesses at H=23 have t−h=0,1,2, respectively. Attach a sink-rooted TT2 by one reversed arc to get h+2t=3h+2(t−h). The replay verifies 69,71,73 by literal counting. Thus the three branches can be supported by different actual sources even at the same height; asking one fixed root to support them all is unnecessary. At H=25 a different architecture is productive: the replay takes one fixed source and three different apex contacts to obtain 75,77,79. At H=67 two sources and three apex choices give 201,203,205. These exact controls use the height to choose the architecture, not just its digit. They do not yet make the resulting continuations ready at the next generation.

A different design allows either signed representative of each nonzero ternary digit: 3h+2 or 3h−4, and 3h+4 or 3h−2, as well as 3h. Every odd target above 5 not divisible by 3 then has two smaller positive odd predecessor choices, differing by 2. They cannot both be holes. Multiples of 3 still use their unique tripling predecessor: 7→21 stays forbidden, while 21→63 is repaired by starting from an actual count-63 seed. Thus **only seeds 1,3,5,63** are needed arithmetically for this overlapping design. The physical obligation is to realize one legal predecessor lift for each target with its continuation intact, not necessarily all signed branches on a single profile. This explains why the pair of holes creates both a need for a separate ordinary lane and a finite restart at 63, without claiming that the signed lifts already cover the whole spectrum.

Replay: `python3 paired_bits.py n7_representatives.txt` with the [archived representatives](../actual_hamilton_multiplicities_20261021/session24-multiplicities/n7_representatives.txt). For the final four literal counts, use `python3 paired_bits.py --fixtures | ./rooted_pair_literal` with the C++ counter from the earlier paired recursion.
