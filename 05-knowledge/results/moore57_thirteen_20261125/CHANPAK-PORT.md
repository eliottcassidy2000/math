# Exactly when arbitrary modules preserve the five-vertex 11-to-13 flip

**PROVED.** There are two marked five-vertex isomorphism classes of single-arc flips from 11 to 13, exchanged by converse. Allow *all three* vertices outside the flipped arc to be replaced independently by arbitrary nonempty tournaments, with each block's external arcs following the old vertex. The two new orientations still differ by exactly one arc. Their counts differ by two in absolute value **if and only if** two specified modules are singletons and the third is transitive. For transitive width \(w\ge1\) at that port the exact physical flip is
\[
                   6\cdot2^w-1\quad\longrightarrow\quad6\cdot2^w+1 .
\]
The entire tournament is strong in both orientations. The theorem classifies simultaneous substitutions away from the arc at every order. Other kinds of extension can retain a unit flip. [All-order 11/13 marked census and converse classification](../chanpak_arc_fiber_grammar_20261122/SELF-CONVERSE.md)

Use outgoing masks \((24,1,3,6,14)\) and change the arc on 1 and 3, indexing from zero. Call the modules at 0, 2 and 4 respectively \(A,B,C\). For a nonempty tournament \(X\), let \(P_r(X)\) count its ordered covers by \(r\) nonempty directed paths, so \(P_1=H\). A compressed Hamilton word visits each module in a positive number of runs; multiply each possible word's count by the corresponding \(P_r\)'s. The coefficient for \(r\) visits to a single marked vertex is divisible by \((r-1)!\): the consecutive returns are simple directed cycles disjoint off the mark. Order these cycles arbitrarily and insert them at the mark of the residual Hamilton path. In particular the two-visit coefficient has the parity of the total number of directed simple cycles through the mark, by Rédei on the residual tournament. [Cycle-packing form of the visit formula](../chanpak_square25_dilation_20261123/EXTENSIONS.md)

## A modulo-four test for any prescribed seam

**PROVED.** Given any same-vertex pair \(T,T'\) and modules \(X_i\) common to both, write \(c_2(T,i)\) for the two-visit coefficient. The inflated difference satisfies
\[
 H(T'[X_i])-H(T[X_i])\equiv
 (H(T')-H(T))\prod_i H(X_i)
       +2\sum_{i:\,|X_i|>1}\bigl(c_2(T',i)-c_2(T,i)\bigr)\pmod4 .
\]
In particular, for a base difference of two, an odd number of inflated sites with changed rooted cycle parity makes a unit flip impossible, either sign. The ordinary height mod four also reads one plus twice the number of simple directed odd cycles by the [odd-cycle-packing identity](https://arxiv.org/abs/2307.05569). The correction here sees all directed return cycles through the individual expanded vertices.

Indeed \(P_r\) is divisible by \(r!\) by permutation of its slots. Two multiply visited modules contribute zero modulo four, as does a module visited at least four times. For a lone three-visit module its even return coefficient supplies the second factor two. For every nonsingleton \(X\), the count \(Q_2(X)=P_2(X)+2H(X)\) admits empty slots, so pairing complementary subsets of its vertices writes it as an even number of pairs of equal odd summands. Thus \(Q_2\equiv0\) and \(P_2\equiv2\pmod4\). The other heights are odd, which yields the formula. Singleton modules have zero repeat count.

The changes in the two-visit coefficients at \((0,2,4)\) for this marked seam are \((-5,6,0)\). Any nonsingleton \(A\) therefore forces difference zero modulo four, independently of the other two blocks. This obstruction already treats both signs and arbitrary intransitive blocks.

## The exact surviving integer difference

Set \(A\) to a singleton. The following are the compressed-word counts with \(i\) visits to \(B\) and \(j\) to \(C\). All unlisted entries are zero.

| \((i,j)\) | lower word count | upper word count | change |
|---|---:|---:|---:|
| (1,1) | 11 | 13 | 2 |
| (2,1) | 6 | 12 | 6 |
| (3,1) | 0 | 2 | 2 |
| (1,2) | 6 | 6 | 0 |
| (2,2) | 4 | 8 | 4 |
| (3,2) | 0 | 2 | 2 |

All changes are nonnegative. The first contributes \(2H(B)H(C)\). Every nonsingleton \(B\) has \(P_2(B)>0\) (cut any Hamilton word), producing a further positive contribution from the second row. Thus exact unit distance requires \(B\) a singleton and \(H(C)=1\), and these conditions suffice. A tournament has one Hamilton path exactly when it is transitive: a cyclic triangle has three distinct paths, each of which can be extended to a spanning path while preserving its old order.

With both earlier modules singletons, \(C\) is entirely arbitrary and the two counts are exactly
\[
       11H(C)+6P_2(C)=6Q_2(C)-H(C),\qquad
       13H(C)+6P_2(C)=6Q_2(C)+H(C).
\]
For a chain, \(Q_2(C)=2^w\), proving the claimed absolute inverse. For arbitrary \(C\), the last display also realizes **every attainable height \(h\)** as the half-difference of an actual single-arc flip between two strong tournaments on just \(|C|+4\) vertices. Each possible filler has a strong expanded host: all original intervertex paths lift, and a directed cycle through a port links any two of its members. Thus for either seed the set of all strictly positive half-differences *at this specific port* is the ordinary tournament Hamilton spectrum, with an explicitly prescribed midpoint. This transfers the original integer-target problem to strong marked pairs without requiring the filler to be strong. The exponent-61 perfect number is 1 modulo 3, so the midpoints in this last display exclude that specific crossing.

For the finite coefficient tables, run the following from the repository root. The compressed-word enumerator is exhaustive because the ports 2, 4 induce an acyclic digraph; it can reuse a port but must consume an ordinary vertex to return.

~~~sh
cd 05-knowledge/results/chanpak_square25_dilation_20261123
python3 -B - <<'PY'
from selector import skeleton
lo=(24,1,3,6,14)
hi=(24,9,3,4,14)
for v in (0,2,4):
    a,b=(skeleton(t,(v,)) for t in (lo,hi))
    print(v,a,b,b.get((2,),0)-a.get((2,),0))
for t in (lo,hi):
    print(skeleton(t,(2,4)))
PY
~~~
