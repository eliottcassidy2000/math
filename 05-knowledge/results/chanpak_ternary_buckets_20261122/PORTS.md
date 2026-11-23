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
