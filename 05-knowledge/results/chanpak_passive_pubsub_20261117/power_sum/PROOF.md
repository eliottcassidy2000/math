# Synchronized changes by matching power sums

The known cyclic two-run reader furnishes an unusually simple compensation mechanism. Substitute a nonempty B in the cyclic order TT_l → B → singleton → TT_l, with l≥1. Exactly two runs of B consume the singleton and then a nonempty subset of TT_l between them, so the two-run skeleton coefficient is 2^l−1. For one run the singleton can be first, last, or between the B-run and a nonempty final segment of the transitive block; the coefficient is 2^l+1. Thus with C₂=P₂−2H the strong lift reads

    (2^l−1)(P₂(B)−H(B))+2H(B).

The exterior subscribes to the difference d=P₂−H as its high-amplitude topic. Changing both inputs together can leave that topic intact. More interestingly, these correlated changes compose across *several* modules, not just on a single pair of blocks.

## Product composition theorem

Suppose we have realizable generators A_i having different odd H-counts h_i but the same positive d=P₂−H. Compose k of them in ordinal fashion. Both H and P₂ multiply, so the new high topic is

    ∏_(i=1)^k (h_i+d) − ∏_(i=1)^k h_i
        = ∑_(j=0)^(k−1) d^(k−j) e_j(h_1,…,h_k),

where e_j are elementary symmetric polynomials. Choose two multisets of k generator counts having the same first k−1 power sums. Newton's identities give the same e_1,…,e_(k−1); the high topic cancels exactly, while the low topic changes by twice the difference of products. The physical gap is constant at every l, including arbitrarily large l. Think of ideal Prouhet–Tarry–Escott solutions, but both sets must be drawn from one *physical* same-d fiber; abstract integer partitions alone do not supply such generators.

This also explains why naively toggling one member doesn't work: its variable part is multiplied by the other members' profiles. The power-sum pairing balances the assembly globally rather than making individual switches pretend to commute without cross terms.

## Two explicit infinite gap families

At order eight there are strong cores with d=4287 and h=241,257,273 (masks in replay). Compare (241,273) and (257,257). Their sum is 514, and their products differ by 256. The lift gives

    20581887·(2^l−1) + 131586     and     +512 relative to it.

For a cubic example order-eight strong cores with d=4803 have h=293,301,317,325. Compare (293,317,317) and (301,301,325). These have equal sum and sum of squares but their products differ by 2048. The strong lift gives

    133559139123·(2^l−1) + 58886554     and     +4096 relative to it.

Both constructions jointly change whole inner modules; no claim is made that their switches are a single arc reversal. The second uses the integer template (0,3,3) versus (1,1,4), scaled by eight and translated by 293. No infinite prime values are claimed. The numerical gaps are fixed; the absolute centers of these two families still run on sparse dyadic strands. To make this avenue select a target, the unknown is a target-dependent same-d fiber whose substituted product polynomial has the required common high topic and sufficiently many product values.

Run `python power_sum/replay.py` from the packet root. It independently checks the full literal strong order-eight cores, the cyclic reader coefficients for small l, the equalities of moments, and the resulting closed formulas. The search for such cores was only a witness finder, not an exhaustive census.
