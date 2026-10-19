# Exact local-screen limitations and an outside-word orbit repair

Date: 2026-10-19. Status: exact scoped examples and a fixed-poset orbit inequality. None is a counterexample to the 1/3–2/3 or Gold Partition conjecture. No arbitrary-exterior conditional migration theorem is asserted.

This packet complements [the positive one-relation twin Gold theorem](one_relation_twin_gold_theorem_20261019.md). It addresses the kind of endpoint-locality issue already illustrated by [the rooted-five Gold obstruction](rooted_five_core_gold_obstruction_20261019.md): a failed restricted certificate must not be mistaken for a failed conjecture.

## 1. Setup and exact reconstruction

Labels0=x,1=y,2=b,3..6 are four common incomparable neighbors of the roots,7..8 are upper slaves, and9=e is lower exterior. The core has b<y and x<7,8. A predecessor mask p[v] has bit u set exactly when u<v. Every list below specifies the full transitive strict order, rather than only Hasse edges.

In each marked example z has exactly two core comparabilities b<z<w and is incomparable with e. The posets are seven-thin. Their size10 is the smallest possible for this required nine-point-core-plus-exterior morphology; no broader minimality claim is made.

The [standalone checker](../experiments/verify_local_screen_orbit_limits_20261019.py) checks closure, cycles and all degrees, independently enumerates extensions and verifies totals by subset DP, and checks all displayed counts, rational inequalities and orbit margins. Its [normal](local_screen_orbit_limits_controls_20261019.json) and [optimized](local_screen_orbit_limits_controls_optimized_20261019.json) receipts agree. All checks remain active under Python's optimization flag.

## 2. Two actual-poset limitations

### Equal-weight root aggregate

    p=(512,516,0,512,512,4,532,513,685,0), z=5,w=8.

There are7,290 extensions, and

    P(x<z)=883/1215,
    P(z<y)=743/1215,
    P(x<z)+P(z<y)=542/405>4/3.

Thus a universal upper bound4/3 for this root sum is false. The pair(z,y) is balanced in this example.

### Both root orientations, even with the lower exterior orientation

    p=(512,516,0,4,512,528,540,513,529,0), z=3,w=6.

Its degrees are(6,7,6,7,5,7,5,7,6,2), and it has8,522 extensions. Exactly,

    P(x<z)=2880/4261>2/3,
    P(z<y)=2865/4261>2/3,
    P(e<z)=4016/4261>2/3.

Therefore no nonzero nonnegative weighted average of these three probabilities has a universal upper bound2/3 across all actual exteriors in this morphology. Since e<x, the third probability is already at least the first. Adding the deterministic orientation z<w cannot repair such an upper bound.

The balanced witness in this example can instead be(b,e):

    P(b<e)=1437/4261.

This does not prove that(b,e) is always balanced, nor that it is forced to be balanced whenever both root probabilities exceed2/3. The latter conditional statement remains unproved here.

## 3. Exact named fibres in a fixed poset

Use the different core

    (0,4,0,0,12,28,4,1,197), z=6,w=8.

Add e=9 below x,y and vertices3,4,5,7,8, incomparable with b,z. Its full masks are obtained by adding bit512 to each core mask except those of b,z and appending zero for e. This core survived both the initial combined analytic screen and the subsequent side-budget screen when this packet was prepared; that is only context for choosing the example, not a prerequisite to any calculation below.

Delete the marked vertices b,e,x,z,y and fix a linear-extension word of the remaining five named core vertices. The exact fibre vector records

    (E,Nxz,Nzy,Nbe).

For the two outside words

    W1=(3,7,4,8,5), W2=(3,7,8,4,5),

the vectors are respectively

    (286,217,188,94), (218,149,158,78).

Both are fibres of the same actual ten-point poset, with the same exterior relations. Artificial raw weights9 and2 produce

    (3010,2251,2008,1002).

The two root probabilities exceed2/3 while P(b<e)=501/1505<1/3. Hence no common nonnegative linear certificate in those three signed deficits can hold on every separate fibre and prove the desired conditional lower bound. This rules out a particular fibrewise linear proof strategy.

It does not give an actual-poset probability vector. In the actual ten-point poset all outside words enter with their ordinary raw coefficient one. For larger exteriors, their weights are constrained by the actual insertion process; independent arbitrary coefficients are a relaxation requiring justification.

## 4. An involution repairs that artificial-mixture defect

The outside induced poset consists of chains3<4<5 and7<8, and has exactly ten extension words. Define τ on every outside word: exchange4 and8 if adjacent; otherwise fix the word. Since4 and8 are incomparable, every exchanged word is valid. The rule is an involution, including its fixed points, and has eight orbits.

Sum fibre count vectors within each orbit. Every orbit satisfies

    6Nzy−15Nbe+E≤0.

Thus the fixed ten-point poset satisfies

    P(b<e)≥[6P(z<y)+1]/15.

In particular P(z<y)>2/3 would imply P(b<e)>1/3. This is a fixed-poset statement proved by all eight exact orbit margins, not an arbitrary-exterior theorem obtained by coefficient fitting.

The previously problematic two fibres belong to one orbit and have sum(504,366,346,172). Their equal raw weighting gives P(b<e)=172/504>1/3. The artificial9:2 weighting violated the involutive orbit's actual equal-weight pairing.

## 5. One additional actual exterior defeats the same coefficients

Add t=10 with t<b, with t below x,y,z,4,5,7,8, and with t incomparable with e and vertex3. The full predecessor masks become

    (1536,1540,1024,512,1548,1564,1028,1537,1733,0,0).

This is seven-thin and has5,095 extensions. Its exact probabilities are

    P(x<z)=760/1019,
    P(z<y)=3164/5095,
    P(b<e)=196/1019,
    6P(z<y)−15P(b<e)+1=9379/5095>0.

Thus the same orbit-summed coefficients are not valid for arbitrary permitted exterior. The new lower vertex changes the source-release conditions and delays b relative to e. The adjacent swap remains a valid involution when its vertices are adjacent in the full outside word; the deficit sign need not remain the same.

This example does not refute conditional migration, because P(z<y)<2/3. The reusable lesson is the need to preserve actual orbit weights and exterior release constraints, not a claim that a general migration theorem is impossible.

## 6. Reproduction and limits

Run from the repository root:

    python3 05-knowledge/experiments/verify_local_screen_orbit_limits_20261019.py

The checker prints a fresh receipt without modifying files. It contains explicit constants for every poset and needs no local exploratory scripts or external libraries. No conjecture-wide conclusion or literature priority claim is made.
