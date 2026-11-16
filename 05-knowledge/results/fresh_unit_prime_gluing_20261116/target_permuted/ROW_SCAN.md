# Changing the exact interface along the old doubling row

Put \(L_s=756\cdot2^s-33\). Among heights 1 through 8, the three
counts \(L_s,L_s+2,L_s+4\) share a rooted continuation at \(g=15\)
exactly at \(s=1,2,5\). The matching terminal \((E,O,g)\) are
\((1107,1122,15),(2897,2912,15),(19437,19452,15)\), respectively.
In particular the fifth level has three strong witnesses of counts
24159,24161,24163. The absence at s=3,4,6,7,8 allows arbitrary rooting,
rewiring, nonstrong hosts, and signed-zero profiles.

## Exhaustion and arbitrary suffixes

The product \(a\) of Hamilton counts of components outside the root's
strong component divides both \(H\) and \(g=15\). Thus enumerate
\(a=1,3,5,15\). Delete the root from that component to get an ordinal
list of SCCs of product \(15/a\). The only nonsingleton types needed
are the triangle, the unique strong count-5 type, and the four strong
count-15 types. The latter are enumerated by fixed Hamilton order
through order six (Busch's bound), pruning monotone count thresholds
via the positive odd-cycle packing identity. Enumerate every factor
order and all singleton slots up to the strong order bound for H/a.

Within each deletion tile compute weighted paths by parity, with
starting and ending vertices specified as needed, where the weight
is the Hamilton count of its omitted vertices. Ordinal decomposition
makes all intermediate tile path choices optional and independent,
so their pair of parity totals convolves in the group ring of C2.
Before and after the endpoint tiles, skipped components contribute
their whole Hamilton counts. This computes exact weighted starts and
the even weighted path table D[i,j]. For the root's outcontact set M,
the count is deletion count plus twice the cut sum of D[i,j] for i in M,
j outside M, by inserting the odd cycle through the root. Contact
bits are traversed in Gray order. The odd and even weighted start sums
give E and O−g. Strongness requires a contact in the first deletion SCC
and an incontact in the last, and these suffice in an ordinal list.

A nonempty ordinal suffix forces E=O. To see it, in the positive pack
expansion for the signed sum, assign each odd cycle or uncovered
singleton of the suffix to one of the two slots; flipping any one
block cancels the sign. For the unsigned sum U=E+O, divide a into prefix
count a/b and suffix count b. A raw terminal U0 becomes (a/b)·R·U0,
where R counts two optional ordered paths covering the suffix. This
cover count multiplies along ordinal components. Scan all divisors
b of a, and all R from ordered SCC factorizations of b; any ordinal
singleton multiplies R by 2. Thus a common signed-zero continuation
exists iff some three enumerated unsigned sums share an odd part:
one can append an adequate number of ordinal sink singletons to
each host to match all three and enforce signed zero. A singleton
root component would have H=g and is out of these windows. The
three positive heights have shared unsigned odd parts 2229, 5809,
38889; the others have none.

## An exact general inflation law

If a vertex v has sole outward neighbor the named root r, the vertex
is a terminal singleton in the deletion; r does not beat it. Replace
v by u transitive twins. Every nontrivial rooted prefix in the
deletion minus v may append any subset of the u sinks, with the
same omitted-subgraph weight; exactly half the choices have either
parity. The root singleton still contributes weight g and cannot
go directly to these sinks. Thus E(u)=2^(u−1)E(1), and O(u)=E(u)+g.
For the odd-cycle insertion cut, paths ending before this terminal
block are unchanged. To stop at an incontact within the block
requires appending a nonempty subset of the right parity. Compared
with u=1, there are 2^(u−1)−1 more such subsets regardless of prefix
parity. The unsigned weighted prefix count equals E(1) by the same
parity-balance calculation. Hence H(u)=H(1)+(2^u−2)E(1).

All three s=5 witnesses have such a vertex immediately before the
named root. Therefore they give a new all-u strong band
\(19437\cdot2^u-14715+\{0,2,4\}\) with common E(u),O(u),g.
This is not a synchronized placement at higher old heights:
if E(1) is odd and H(1)=L_s, the intercept H(1)−2E(1) is
1 mod 4, whereas every L_t is 3 mod 4. Inflating the terminal
slot further misses the prescribed first count along the old row.

## An exact inverse source at the fifth height

There is a strong tournament with full profile (8053,6479,6484,5).
Its outgoing masks are

    16374 32764 32760 16369 32736 32704 16256 32512
    32256 31744 30720 12288 8192 16384 6217

The root is last and its predecessor only beats the root. The old
first-stage +4 output has scalar count 8053 at s=5, but the physical
source here has a changed continuation. Inflating its terminal
slot gives \(H'=6479\cdot2^u-4905\), \(E'=6479\cdot2^{u-1}\),
\(O'=E'+5\), \(g'=5\). The entire three-output replacement band
above has exactly \(3P'+(d,0,0,0)\), with the same d=0,2,4 all
heights. This is an exact positive re-encoding from the original
scalar 8053 at its one aligned height; it does not use old induced
arcs and is not a universal old-row selector. The exact rows are also checked in **python3 check_permutation.py**.

## Replay

Compile **c++ -std=c++17 -O3 replay_rows.cpp -o /tmp/g15-row**. Run **/tmp/g15-row s** with s = 1 through 8. A support bitmask of 7 in standard output means all three counts share that terminal continuation. Standard error prints the signed-zero screening result and literal rows. The separate C++20 selector in this folder finds strong triples by a target-pruned cut search instead of a Gray traversal.
