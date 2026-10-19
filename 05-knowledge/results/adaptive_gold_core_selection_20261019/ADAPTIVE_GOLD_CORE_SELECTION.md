# A common-core first-comparison obstruction, with adaptive Gold repair

19 October 2026. Status: FINITE-EXACT obstruction and segment repair; PROVED two-question selection lemma. Root reviewed and approved. No general balance or Gold theorem, counterexample, or historical novelty is claimed.

## 1. Result and exact scope

Two actual seven-thin posets, on 12 and 13 vertices, have the same induced ten-vertex full critical-root core. Their sets of 1/3–2/3-balanced CORE pairs are disjoint. Consequently no first comparison chosen solely from that common core, independently of the exterior, can begin a two-comparison Gold strategy for both posets. Allowing different second comparisons on either first branch does not fix this obstruction.

There is a positive repair on the entire segment between their projected core laws. The same two questions

    A: 0<2,       B: 5<8

produce strict Gold after selecting which question comes first. Every mixture on that segment has all four Gold margins greater than 1/25 of total mass under the explicit selection rule below. Thus this example isolates the need to adapt the FIRST comparison rather than an inability to partition these laws.

The obstruction concerns core comparisons. It does not exclude a common comparison involving exterior vertices, does not rule out different pairs for individual posets, and does not show that every interior point of the mixture segment is realizable by a poset.

## 2. Two small ambient posets

Use strict predecessor bitmasks on labels 0 through 9:

    Q=(0,0,0,6,0,19,2,18,86,151),       roots=(0,3).

This is the established ten-vertex no-good-pair seed. Its earlier Gold analysis is in [packet136](../conditional_gold_exchange_and_minimal_obstruction_20261019.md); the independent exchange-capacity obstruction is in the [support-three Hall note](../gold_support_three_hall_20261019/GOLD_SUPPORT_THREE_HALL_BARRIER.md).

We use the chain-profile convention proved in the [boundary chain compression theorem](../boundary_chain_compression_20261019/BOUNDARY_CHAIN_COMPRESSION.md). A lower chain of length h has c_i<q exactly when i<=h-alpha(q), with i numbered from 1; an upper chain has q<d_j exactly when j>beta(q). Every lower-chain vertex is below every upper-chain vertex.

The two profile pairs, in core label order, are

    P0: alpha=(0,1,0,0,1,0,0,0,0,0)
        beta =(0,0,0,0,0,0,0,0,1,0)

    P1: alpha=(0,1,1,0,0,0,0,0,0,0)
        beta =(0,0,0,0,0,1,0,1,1,2).

For direct checking without the profile theorem, their full predecessor tuples are

    P0=(1024,0,1024,1030,0,1043,1026,1042,1110,1175,0,1791)

    P1=(1024,0,0,1030,1024,1043,1026,1042,1110,1175,0,1119,3583).

Both have incomparability degree at most 7. The full set consisting of roots 0,3 and all their incomparable neighbors is exactly Q, and D(0) is contained in D(3), while U(3) is contained in U(0), in both ambient posets. The extension totals are

    e(P0)=29193,       e(P1)=84707.

The complete balanced-core-pair lists are

    P0: {0,2}, {1,4}, {3,5}, {5,9}, {6,7}, {8,9}
    P1: {0,4}, {0,6}, {1,2}, {2,4}, {3,7}, {5,8}.

The lists are disjoint. The exact certificate records all 28 incomparable core-pair numerators, not just these positive lists. Reversing a pair preserves balancedness, and comparable pairs cannot be balanced, so these checks exhaust all first core comparisons.

Why balancedness is necessary for Gold: normalize total mass to 1 and let h>=1/2 be the heavy first-branch mass. Its second comparison splits it into c and h-c. Gold requires both h+c<=1 and h+(h-c)<=1. Since max(c,h-c)>=h/2, this forces h<=2/3. Thus the two disjoint lists exclude a universal first core pair even when the second comparison may be selected separately for each ambient poset and each first outcome.

## 3. Exact adaptive repair of the whole mixture segment

Let mu_i be the core-order distribution induced by uniform extensions of P_i, and put

    mu_t=t mu_0+(1-t)mu_1,       0<=t<=1.

This is a normalized mixture of projected laws. Its coefficients are not unnormalized exterior-order multiplicities.

For questions A=0<2 and B=5<8, the cell counts in order

    (A yes/B yes, A yes/B no, A no/B yes, A no/B no)

are

    P0: (10042,3382,9710,6059),       total29193;
    P1: (16922,9696,29630,28459),     total84707.

Divide each row by its OWN total before mixing. These four counts alone certify the following rule:

    if t<=1/2, ask B first and A second;
    if t>=1/2, ask A first and B second.

At t=1/2 either order works. The second question is the same on both branches; no branch-dependent choice is needed here.

For a normalized cell of mass c inside a first branch of mass b, its Gold margin is 1-b-c. Each of the four margins is affine in t for a fixed question order. Therefore checking the endpoints of each half interval proves the rule throughout. The exact verifier checks all four margins at t=0,1/2 with B first and at t=1/2,1 with A first. Their common minimum is

    113042071 / 2472851451 > 1/25.

This gives a strict positive margin for every law on the segment, including any realizable interior laws. It is a segment theorem, not a claim covering the full chain-profile hull of Q.

## 4. A general two-question selection lemma

LEMMA. For any probability distribution and any fixed two binary questions, if one order of asking the questions is weak Gold, then asking the more balanced question first is weak Gold. The analogous statement holds for strict Gold. The same second question is asked on both first branches.

Proof. Orient each question so its affirmative event is a heavy side. Let the heavy probabilities be h<=g, both at least 1/2, and let z be the probability that both heavy events occur. With the h-question first, its heavy branch has cells z and h-z. With the g-question first, its heavy branch has cells z and g-z. A light first branch automatically satisfies weak Gold because its mass is at most 1/2. Therefore the respective nonautomatic conditions are

    h+max(z,h-z)<=1,       g+max(z,g-z)<=1.

The function u+max(z,u-z) is nondecreasing in u, proving the weak assertion. For strict Gold the heavy inequalities are preserved strictly. If h>1/2 the light branch is automatic strictly. If h=1/2, failure of strictness on a light cell would require that cell to have mass 1/2; then one of the g-question's cells has mass at least 1/2 inside a branch of mass at least 1/2, contradicting strict Gold in that order. This proves the strict assertion too.

There is a compact exact scalar criterion for WEAK Gold. For a nonnegative 2-by-2 cell array of total E>0, let M be its largest cell, R its largest row sum and C its largest column sum. Then

    rows-first is weak Gold iff R+M<=E;
    columns-first is weak Gold iff C+M<=E;
    some order is weak Gold iff min(R,C)+M<=E.

If the largest cell lies in a heavy row, the first condition is exactly the heavy-row test. If it lies in a light row, M<=E-R and all heavy-row cells are at most M, so both the scalar condition and Gold hold automatically. The column argument is identical. This is a useful piecewise-linear form of the adaptive selection rule.

For STRICT Gold, replacing <= by < is sufficient but not necessary. For example, cells (3,3;4,0) have E=10, R=6 and M=4, while all rows-first path sums are at most 9. Thus the strict segment certificate above uses the actual four branch margins rather than treating R+M<E as an equivalence.

For the present segment, asking whichever of A and B is more balanced therefore also gives strict Gold. The switch is exactly

    t=336799641/674187622.

Indeed P(A)<1/2<P(B) throughout, so the switch condition is P(A)+P(B)=1. The verifier also checks that this observable more-balanced-first rule retains the stronger uniform margin >1/25. Thus one need not know a hidden mixture decomposition to choose the order: the two comparison probabilities suffice.

## 5. Reproduction and next use

Run `python verify_adaptive_gold_selection.py` or `python -O verify_adaptive_gold_selection.py`. Both use only integer ideal dynamic programming and exact rational arithmetic and produce the same 424-guard receipt. The dynamic program counts each ambient extension directly, with added comparisons for each event; it does not rely on the earlier profile-weight formula, numerical optimization, or discovery arrays. Its guard set checks thinness, the full critical-root core, all first-core-pair counts, all selected joint cells, the affine interval certificates and elementary edge cases for the ordering lemma. The manifest freezes script and receipt hashes.

For a general finite boundary-atom model, every fixed two-comparison tree has an explicit polyhedral validity region in the normalized atom weights: each Gold margin is affine. Adaptive strategies amount to covering the realizable parameter set by these regions. The two models above prove that one fixed first core pair may be insufficient, while their segment admits a two-region cover with a comfortable strict margin. The full profile universe and its nonlinear realizability restrictions remain the next target; no cover of that entire universe is supplied here.
