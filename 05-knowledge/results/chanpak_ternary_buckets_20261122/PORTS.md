# Nine vertices are necessary for a one-reversal ternary fork

Let \(H\) count directed Hamilton paths. Call a center a tournament of height \(3h\) with individual arc reversals to \(h,3h-2,3h+2\). For \(h>1\) the three target heights are distinct. All statements in the first table allow **every** choice of tournament and of the three arcs. There is no apex, distinguished vertex, regularity, strong-connectivity or flank assumption.

**FINITE-EXACT.** The smallest order of a center with \(h>1\) is nine. At that order there are exactly four ordinary center isomorphism classes, two at each of the following heights.

| Source \(h\) | Middle \(3h\) | Other one-arc neighbors | Ordinary centers at order nine |
|---:|---:|---:|---:|
| 29 | 87 | 85, 89 | 2 |
| 279 | 837 | 835, 839 | 2 |

To ask which of these can propagate the same three *physical* edits, choose a vertex of the center and expand it into a uniform transitive tournament of arbitrary positive size. Each outside vertex retains its orientation toward every vertex of the replacement. The edits are fixed outside the replacement. A growing solution has unbounded source height as the replacement grows. We include the possibilities with \(h=1\) at width one in this search.

**FINITE-EXACT classification with PROVED propagation.** Through nine vertices there is precisely one growing source function in this model,
\[
                   H_S(T_s)=21\cdot2^s-13.
\]
The source functional on a general nonempty tournament \(M\) replacing the same vertex is \(29H(M)+21P_2(M)\), where \(P_2\) counts ordered covers by two nonempty directed paths. Its center is three times this functional and its side branches differ from the center by \(\pm2H(M)\). The two height-837 classes do not give a second growing function at this order by expanding one vertex away from the three edits. A canonical center from the replay has adjacency rows \((494,508,504,368,481,448,128,264,64)\), with port 0. Reversing \((6,7),(3,5),(7,8)\) gives the source, lower and upper neighbors, respectively.

## Enumerating all arcs and all ports

For a tournament on vertex set \(V\), use the usual subset path recurrence in both orientations. For distinct \(i,j\) define \(D(i,j)\) as the sum over \(S\ni i,\ S\not\ni j\) of the number of paths ending at \(i\) on \(S\), times the number starting at \(j\) on \(V\setminus S\). Cutting at a proposed arc gives a bijection with these ordered pairs. Thus reversing \(i\to j\) changes the height by \(D(j,i)-D(i,j)\). This evaluates all one-arc neighbors directly on each ordinary representative.

For each possible common center retain its sets of qualifying arcs. For each proposed port \(v\), discard incident arcs. If its outdegree and indegree are \(d^+,d^-\), at most \(K=1+\min(d^+,d^-)\) occupied visits can pass through the expanded block. In a Hamilton path all visits except possibly the first need distinct outside predecessors, and all except possibly the last need distinct successors. None of the retained edits changes these bounds. For each width \(1\le s\le K\), count the expanded graph and its neighbors by independent subset recurrences and intersect the sets of surviving edges. The triangular ordered-cover interpolation in [the block substitution lemma](PHYSICAL-FORKS.md#what-survives-replacing-a-transitive-block-by-any-tournament) proves that exactly these arcs obey the required identities at every width. Positive source coefficients in slots \(k\ge2\) distinguish growing solutions from constant ones.

The numbers of ordinary tournaments checked at orders 2 through 9 are \(1,2,4,12,56,456,6880,191536\). At orders 4 through 9 there are respectively \(2,3,4,5,6,7\) ordinary centers with \(h=1\) if the source and lower neighbor may coincide. None supplies a growing solution. Order one has no arc. [McKay and Piperno's nauty](https://pallini.di.uniroma1.it/) supplies one representative per ordinary isomorphism class. Run from this directory.

~~~sh
c++ -std=c++17 -O3 -Wall -Wextra -Werror port_forks.cpp -o /tmp/cp-port-forks
for n in 2 3 4 5 6 7 8 9; do gentourng -q "$n" | /tmp/cp-port-forks "$n"; done
~~~

## A classification at every tournament order

**PROVED with finite-exact strong cores.** For odd \(1<h\le71\), an individual tournament of height \(3h\) has reversals to all three heights \(h,3h-2,3h+2\) at some order if and only if \(h\in\{29,71\}\). The ordinary generating functions for the numbers of possible *center isomorphism classes*, by order, are respectively
\[
                  \frac{2z^9}{(1-z)^2},\qquad
                  \frac{2z^{10}}{(1-z)^2}.
\]
The counted object is the center, with the three moves required to exist. Additional exits or a choice of labels do not weight it.

By [core locality for a difference of two](../chanpak_collatz_edge_buckets_20261123/README.md#3-a-general-no-go-for-the-cheapest-ternary-lift), both side reversals take place inside a unique nonsingleton strong component. The remaining reversal **lowers** the height from \(3h\) to \(h\), so it too is internal. Indeed a cross-component reversal can never lower the height. Every old path that uses the old cross arc \(u\to v\) visits these two vertices consecutively. Swapping them gives a new path using \(v\to u\). The predecessor of \(u\), if present, belongs to its own component or an earlier one and beats \(v\). The successor of \(v\), if present, similarly loses to \(u\). This injects all lost paths into gained paths, while paths avoiding the arc survive. The three moves therefore descend simultaneously to the strong core, and any singleton padding before and after that core preserves all three.

[Busch's strong-tournament bound](https://doi.org/10.37236/1141) gives \(H(S)\ge5^{(|S|-1)/3}\). When \(|S|\ge11\) this exceeds 213, since \(5^{10}=9,765,625>213^3=9,663,597\). Thus a complete census through order ten decides every source through 71 at *all* orders. The order-ten census examined all 9,733,056 ordinary classes. The table gives all of its nontrivial source heights (the other eight centers have source 1).

| Source | Ordinary centers at order ten | Strong among them | Self-converse strong centers |
|---:|---:|---:|---:|
|29|4|0|0|
|71|2|2|0|
|117|2|2|0|
|139|2|2|0|
|279|4|0|0|
|431|2|2|0|

At each of 29 and 71 there is precisely one pair of strong centers exchanged by reversing every arc. This operation preserves path counts by reading a path backward. The two centers in each pair are nonisomorphic, and ordered singleton padding before and after a given center has series \((1-z)^{-2}\). This proves both generating functions. Among growing one-port functions the order-ten sweep yields only the order-nine function and its shift \(h=42\cdot2^s-13\).

The same converse action illustrates the difference between asking for a count and asking for a fork. For heights 5, 11 and 13 there is a unique smallest ordinary class, necessarily self-converse. The smallest order at height 29 is six, with four classes, two self-converse and one pair. At height 71 the smallest order is seven, with seven classes, one self-converse and three pairs. The two center generating functions above impose the three arc reversals. The [canonical replay](center_classes.py) computes the counts for both questions.

To extend the command above through order ten and summarize the strong cores, save the actual center rows using the optional flag. The 10-vertex stream can also be processed by the disjoint residue/modulus shards of gentourng.

~~~sh
for n in 2 3 4 5 6 7 8 9 10; do gentourng -q "$n" | /tmp/cp-port-forks "$n" --bases; done > /tmp/cp-centers
python3 -B center_classes.py < /tmp/cp-centers
~~~

## Two independently growing vertices still give one effective width

**FINITE-EXACT classification, PROVED for unbounded widths.** Take *any* tournament through ten vertices and choose two distinct vertices to replace by transitive blocks of sizes \(s,t\). Keep the three candidate flips fixed outside both blocks. Require a growing center of height \(3h\) with reversals to \(h,3h-2,3h+2\) for every \(s,t\ge1\). Up to exchanging variables, the only possible source functions are
\[
         h(s,t)=21\cdot2^s-13,\qquad
         h(s,t)=21\cdot2^{s+t}-13 .
\]
The second option first occurs at order ten. Its two variables enter only through their sum. Even allowing any fixed positive integer rates \(a,b\) and requiring the equalities only at \((s,t)=(1+ak,1+bk)\) for every \(k\ge0\) produces exactly the restrictions of these same two functions. Starting a synchronized progression away from \((1,1)\) is a different problem.

For the proof, remove both blocks from a Hamilton path. In each gap of the exterior permutation the occupied portion can be empty, from the first block, from the second, or one visit to each block in the direction of their mutual arc. Each surviving exterior pattern therefore gives an ordered pair of visit counts \(i,j\). Its contribution at arbitrary nonempty input tournaments \(A,B\) is \(P_i(A)P_j(B)\). Counting all patterns by a subset recurrence gives the coefficient table \(b_{ij}\). The product of the triangular matrices for the two inputs makes this table unique. Fixed exterior edits can be tested by comparing the tables. On fully independent widths the additive shift of two changes only the coefficient \(b_{11}\).

For synchronized widths, expand \(i!{s\brace i}=\sum_{x=1}^i(-1)^{i-x}\binom ix x^s\). Each term of the resulting height has exponential base \(x^a y^b\) as a function of \(k\), with coefficient including the initial factor \(xy\). Group equal bases and compare coefficients. Distinct bases are linearly independent by Vandermonde. A common divisor of \(a,b\) does not change which bases coincide. Suppose the rates are coprime and \(x^a y^b=X^a Y^b\) with different ordered pairs. At a prime where \(x\) differs in valuation from \(X\), the relation \(a(v(x)-v(X))=b(v(Y)-v(y))\) forces both rates to be at most three. Indeed the gap count bounds \(1\le x,y,X,Y\le9\), so each valuation difference is at most three in absolute value. At larger reduced rates every base is distinct and the independent-width test suffices. The only remaining coprime rate pairs are \((1,1),(1,2),(2,1),(1,3),(3,1),(2,3),(3,2)\).

All 55 candidate ordinary centers through order ten at width \((1,1)\) were included, even the centers of height three. There are 983 possible port pairs surviving the ordinary three-neighbor filter, of which 668 have growing central height. For each rate test the same six rooted presentations succeed, all at order ten. Four have source coefficient table \((b_{11},b_{21})=(29,21)\) up to swapping the ports. Two have \((b_{11},b_{21},b_{12},b_{22})=(71,42,42,21)\). The second table gives \(21Q_2(A)Q_2(B)-13H(A)H(B)\) on arbitrary inputs. The [two-port replay](two_ports.py) checks a planted positive case and prints every successful pair.

~~~sh
python3 -B two_ports.py < /tmp/cp-centers
~~~
