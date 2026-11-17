# The first ready cores for 7p, and what 18 measures

**Results.** The fixed-deletion ready spectra are finite unions of dyadic rays at all orders. A size estimate from the strong-component floor makes low deletion counts quantitatively insufficient. Exact small-deletion classification explains why 1 and 3 cannot rearm eligible 7p targets, and why modulus 18 selects the first deletion-count-five two-run unit family. There is an explicit permanent 1211/1213 pair, with 1211=7·173. This does not prove a ready re-encoding for every target.

## 1. The numbers 1, 2, 3 and 7, 18, 21

The first two odd counts mean different structures: H=1 is transitivity; H=3 is one directed triangle with possible ordinal singleton padding. In the [odd-cycle expansion](https://arxiv.org/abs/2307.05569), 1 is the empty pack, and a single directed triangle adds 2. The two optional ordered path slots appearing below give a second, literal role for 2. It is not itself a permissible Hamiltonian count.

Multiplying the first two odd counts by 7 gives 7 and 21; both are excluded by the [strong-order floor and finite core proof](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md), not by congruence. The subtraction 21−3=18 gives no physical 18-step from the triangle. But **18 is an actual continuation coefficient**: the [15/17 rooted pair](../one_port_neutral_continuations_20261117/README.md) has deletion count 3 and 2E=18. All near-sinks carry it to counts 18·2^k−3 and 18·2^k−1. With this readiness criterion, the lower member is always divisible by 3. The number 18 appears differently, as a genuine phase filter, when the deletion core changes to 5 (§4).

## 2. An exact positive readiness test

Fix a root r, with D=T−r and g=H(D). Write E,O for the weighted even/odd root-start counts (weight of a path is the Hamiltonian count of its complement). Call the root normalized when O−E=g, and put d=H−2E+g. Ready in this note means normalized and d=0, the condition in the existing ternary-ready construction. It is not by itself a certificate that a specific next digit can be deposited.

Let U be the outneighbors of r. Let K count even-vertex paths in D starting and ending in U, weighted by H of their complement in D. Let J count them with start in U and end outside U. The odd-cycle expansion at r gives H=g+2J. Prepending r gives O=g+J+K, with the root singleton contributing g. Hence

    d = 2(O−E−K); in the normalized sector, d=2(g−K).

So a ready root satisfies **K=g**. This is a positive count, not a cancellation target; a search can reject a partial contact assignment as soon as it already exceeds g. The equality d=0 has more bite than its scalar appearance suggests.

## 3. Only a dyadic tail can grow with bounded g

Order the strong components of D. Preceding singleton components before the first component meeting U are global sources in T and may be deleted without changing the profile. Between the first and last outcontact components, fix a vertex of U in each extreme when those components differ. If there are s intervening singleton components, each even-cardinality subset of these vertices gives a distinct even path between the extremes, each of weight at least one. Thus s≤1+floor(log₂ K) for K>0. Only two endpoint singleton components may be additional. If the extreme components coincide there are no intervening singletons. The cases used below have K>0.

Let L end at the last component meeting U and let M be the remaining ordinal suffix. Every member of M beats the root. If M is nonempty, let R=P_M(2) count ordered two-slot spanning covers, slots allowed empty. The even and odd allocations of vertices to the first slot are equal: in the signed sum over subsets use the odd-cycle expansion on each part, then allocate every odd cycle and every uncovered singleton freely to either part. Each contributes a factor 1−1. Let F count the nonempty paths of L starting in U, weighted by their complement in L. The rooted profile has E=FR/2, O=FR/2+g, and K is the corresponding L-count times H(M). Removing singleton components from M divides R by the corresponding powers of two but leaves H(M) unchanged. Retain one singleton whenever needed to keep M nonempty. Thus normalized profiles keep g,d and the odd part of 2E, reducing its exponent only. If M is empty there is no tail to compress. Conversely ordinary near-sink growth preserves g,d at a normalized root and doubles 2E, so every reduced ray really continues. The exact law at fixed g,d is therefore a finite union

    H = A·2^t − g + d,   with odd A and t bounded below per ray.

The same compression holds for d in any fixed set with K>0 and K=O(g); in particular it covers d=0,±2 for g≥3. This is an all-order restriction even allowing arbitrary count-preserving re-encodings, not just an obstruction to one switch.

**Stronger quantitative bound.** The nonsingleton components of D have at most 3 log₃ g vertices in total: multiply the component counts and use [Busch's strong floor](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3), which implies H(C)≥3^(|C|/3) (the triangle is equality). Add the active singleton bound, the retained suffix singleton, and the root. For d=0 the whole template has at most 5+3 log₂ g vertices. Counting even all labeled tournaments on this many vertices yields at most

    2^((9/2)(log₂ G)^2 + O(log(2G))) · O(log(2X+G))

distinct ready values H≤X with g≤G. In particular for G=2^(sqrt(log₂ X)/3) this is at most X^(1/2+o(1)). Only a vanishing fraction of primes p or numbers 7p≤X can use such small g. Thus the proposed landing/rearming separation needs **growing deletion complexity or a different ready criterion**. The original background proof had a weaker exponential-in-G description; the strong floor is what sharpens it here.

## 4. Three complete small-deletion tests; a neutral pair at 1211

**g=1.** D is transitive. K=1 forces exactly two outcontacts (three already give three two-vertex paths) with at most one intervening vertex. Normalization requires the root to lose to the final deletion vertex; see [the one-root bit law](../one_port_neutral_continuations_20261117/README.md). The two words after trimming are 110^t and 1010^t, t≥1. They give exactly 3·2^t−1 and 5·2^t−1, neither divisible by 7: the powers of 2 mod 7 are 1,2,4.

**g=3.** D consists of a triangle padded by ordinal singleton components. Its compressed contact word gives E+O divisible by 3, hence readiness makes H divisible by 3. None can be eligible 7p. The exact ray list is A=9,15,27, all t≥1. The complete finite reduction needs at most five singleton vertices total alongside the triangle; replay.py checks all 2,568 masks on those shapes and verifies the lowest ray witnesses by literal path DP.

**g=5.** Being prime, g forces one nonsingleton component. The strong floor forces order at most four, and the unique strong 4-tournament has H=5. Alongside it retain at most six singleton vertices, sufficient also for normalized d=±2 since K≤6. Replay checks 12,304 contacts, with literal subset checks of the least seeds. The following table lists **every** ready ray and the least exponent of any ±2 companion with the same odd part A of 2E. A companion at a smaller exponent can be grown by near-sinks to the ready exponent. Empty cells mean there is no such profile at any height with g=5.

| odd A | least ready t | least t at d=−2 | least t at d=+2 |
|---:|---:|---:|---:|
| 3 | 5 | | |
| 5 | 5 | 3 | |
| 9 | 2 | | |
| 15 | 1 | | |
| 19 | 2 | | 1 |
| 21 | 1 | | |
| 23 | 1 | | |
| 25 | 1 | | |
| 71 | 1 | | |

Of the ready rays, only A=3,5,19 can be divisible by 7. A=5 is always divisible by 5 and never 35 at its legal exponents. A=3 gives the ready 7·13=91 at t=5 but has no matched neighbors. Thus all eligible g=5 *ready plus matched unit* counts must lie on A=19. Divisibility by 7 requires t≡0 mod 3. When t≡3 mod 6 the count is divisible by 3, so it cannot be eligible 7p. The remaining lane is

    H = 19·2^(6k)−5, k≥1, with H≡5 mod 18 and p=H/7≡11 mod 18.

This is a necessary phase condition within the g=5 matched-ready architecture, not a prime-generating assertion.

The A=19 start has profiles (71,5,38,43) and (73,5,38,43) at the retained root 6, with outgoing masks

    (62,36,104,114,102,64,3)
    (242,244,185,243,160,192,20,64).

After four ordinary near-sink moves the counts are **1211 and 1213**, both at (g,E,O)=(5,608,613). Here 1211=7·173 and 1213 is prime. Any further common sink/near-sink word preserves the +2 gap, and removing the sole live root limits the appended tape to two runs, as in the one-port construction. Only the normalized ready condition has the sparse phase restriction above; arbitrary future words still preserve the gap but need not keep readiness. This is a new core of the prescribed type, not a claim that either scalar count was previously unattainable.

Run python replay.py and python paired.py in this directory. The general reduction proves completeness of the finite mask bounds; the programs check the ray tables and literal pair independently by subset path counts.
