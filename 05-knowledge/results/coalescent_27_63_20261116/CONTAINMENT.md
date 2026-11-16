# Exact prohibition of an induced 27→63 lift, and sharp rewrites

Let H be the directed Hamiltonian path count. The known exclusion of H=7 and H=21, the multiplicative SCC factorization of H, and H≥3 for a nontrivial strongly connected tournament show that a tournament of count 63 has exactly one nontrivial SCC, of count 63. All other SCCs are singleton padding. Any induced source meets that padding only in neutral vertices, so its count equals the count of its restriction to the strong core. The inherited Busch bound places the H63 core at order at most eight, because 63³=250047<5⁸. The complete n8 atlas sees smaller cores too via singleton padding. The script check_27_63.py checks all restrictions and gives exactly two strong H63 classes, both of order eight. They are distinct converses and have representatives (2,16,3,23,165,79,31,111) and (16,5,17,7,74,31,175,63) in outgoing row masks. Neither contains even one induced H27 subtournament. This excludes an induced H27→H63 extension at **any order**, not merely a one-vertex extension. H63 itself is attained.

## Strong H27 sources and exact repair cost

The complete atlas gives two strong H27 classes of order six and three of order seven; no higher strong H27 orders occur (27³<5⁷). The reverse search in rewires27.py exhausts every restriction of both strong H63 classes and all reversals of at most two edges in each restriction. It gives the following cost to a **strong H63 output** when all old vertices must stay:

| H27 class | Minimum reversed old arcs | Added vertices |
|---|---:|---:|
| each of the two order-six classes | 1 | 2 |
| two of the three order-seven classes | 1 | 1 |
| the remaining order-seven class, representative (2,4,17,7,43,79,31) | 2 | 1 |

If neutral padding is permitted in the output, all five classes need exactly one old reversal. To confirm the last exceptional class, take its relabeling (72,1,67,102,79,87,2). Reverse {1,6}, and then add two vertices. In the final target (16,5,17,7,74,31,175,63)→TT1, the old vertices are the vertices 0,1,3,4,5,6,8 in order; the newcomers are vertices 2 and 7. With only one newcomer the eight-vertex H63 output would necessarily be strong, so the extra vertex is essential at one reversal for this class.

## Connected strong runs

Starting from the literal original H27 seed (62,124,120,112,32,64,17), all steps below reverse exactly one old edge and append exactly one vertex. Every output is strongly connected. The first step is globally optimal in both budgets.

| Target | Reverse pair | Append outgoing mask |
|---:|---|---:|
|63|{3,6}|63|
|189|{0,2}|24|
|567|{0,6}|288|
|1701|{4,8}|304|
|5103|{3,6}|16|
|15309|{0,1}|56|
|45927|{0,1}|8056|
|137781|{3,6}|2576|

These are literal witnesses, not a proof that an unlimited strong ternary multiplier exists. The orbital count identity via adjoining C3 applies indefinitely but gives a decomposable target.

## A separate exact no-contact check

For each order-six H27 class take its ordinal sum with C3, in both block orders. Every such source has H81. No outgoing mask on the new vertex gives H189 in any of these four cases. In particular this rules out the potential one-new-vertex H189 bridge from an H27 strong nucleus of order six and one entire directed triangle. The standalone check_constructions.py validates all 2048 masks and uses a full subset path DP to verify the strong chain.
