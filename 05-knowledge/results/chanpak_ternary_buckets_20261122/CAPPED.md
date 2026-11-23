# All single-arc ternary forks through source height 501

**PROVED with a finite-exact core census.** Let \(1<h\le501\) be odd. A tournament of height \(3h\) has individual arc reversals to each of \(h,3h-2,3h+2\) if and only if
\[
                  h\in\{29,71,117,139,155,279,323,431\}.
\]
For each listed height exactly two strongly connected center isomorphism classes exist. They are exchanged by reversing every arc. The table gives the order \(m\) of both. At every order, the ordinary generating function counting *all center isomorphism classes* for each such source is \(2z^m/(1-z)^2\); at every other odd source in the stated range it is zero. The [locality argument](PORTS.md#a-classification-at-every-tournament-order) puts the three moves inside a unique strong component, and any number of singleton components can lie before and after it in the ordered decomposition.

| Strong center order \(m\) | Source heights \(h\) |
|---:|---|
| 9 | 29, 279 |
| 10 | 71, 117, 139, 431 |
| 11 | 155 |
| 12 | 323 |

The two new higher-order sources, 155 and 323, belong to the [direct family](DIRECT-FORK.md) with \(s=3,4\). Each has precisely two strong centers even when arbitrary orders and unrelated constructions are permitted.

## Why a height cap closes the search

Start with the directed triangle. Every strong tournament of order \(n\ge4\) contains a directed cycle of length \(n-1\); its vertex set induces a strong parent. Every Hamilton path on the parent admits an insertion of the deleted vertex, and choosing the first possible insertion gives an injection since deletion recovers the old path. Thus every parent obtained in this way has height at most that of its child. Conversely a one-vertex extension of a strong parent is strong exactly when the new vertex has both an incoming and an outgoing arc to the parent. This gives an exhaustive, duplicate-tolerant generation using just those nonconstant masks. Once a generation is empty, all higher generations are empty. [Moon, *On subtournaments of a tournament* (1966)](https://doi.org/10.4153/CMB-1966-038-7)

For efficiency, compute the old path starts \(S_i\) and ends \(E_i\), and let \(D_{ij}\) count ordered path pairs covering the old vertices whose first path ends at \(i\) and whose second path starts at \(j\). If \(M\) is the set of old vertices beaten by the new one, the extended height is
\[
            \sum_{j\in M}S_j+\sum_{i\notin M}E_i+
                        \sum_{i\notin M,\ j\in M}D_{ij}.
\]
This follows by splitting at the new vertex, including the possibilities that it comes first or last. Adding one element of \(M\) updates this reading in linear time. The same matrix evaluates each individual reversal at the parent, since reversing \(i\to j\) changes the count by \(D_{ji}-D_{ij}\). The [replay](capped_cores.cpp) computes both endpoint tables by subset recurrence and canonicalizes each surviving child as a directed graph with [nauty](https://pallini.di.uniroma1.it/). It checks each augmented height against a fresh endpoint sum when processing that child.

| Order | Strong classes with height at most 1503 | Strong fork centers |
|---:|---:|---:|
|3|1|0|
|4|1|0|
|5|6|0|
|6|35|0|
|7|353|0|
|8|6,008|0|
|9|91,986|4|
|10|108,185|8|
|11|65,636|2|
|12|27,875|2|
|13|6,941|0|
|14|490|0|
|15 and greater|0|0|

Run from this directory with a C++17 compiler and nauty. Each success prints canonical directed adjacency rows, one arc for each of the three heights, and whether the class is self-converse. The optional final flag also prints each order and height's strong class counts, followed by the counts of classes with a negative step, positive step, and both steps of size two.

~~~sh
# With Homebrew nauty, expose its pkg-config metadata.
export PKG_CONFIG_PATH="$(brew --prefix nauty)/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
c++ -std=c++17 -O3 -Wall -Wextra -Werror capped_cores.cpp $(pkg-config --cflags --libs libnauty) -o /tmp/capped-cores
/tmp/capped-cores 1503 15 --fibres > /tmp/capped-fibres
~~~

## Every ordinary height series through 1503

The 307,517 strong classes in the table simultaneously determine the ordinary class series for all these heights. If \(S_d(z)\) counts strong isomorphism classes of height \(d>1\), marking their orders, the exact series for all ordinary classes of height \(h\) is
\[
 \sum_{k\ge0}\frac{1}{(1-z)^{k+1}}
       \sum_{\substack{d_1\cdots d_k=h\\d_i>1}}
            S_{d_1}(z)\cdots S_{d_k}(z).
\]
The factorization is ordered; singleton components occupy the \(k+1\) gaps independently. At height one the empty word contributes \(1/(1-z)\), including the empty tournament at order zero. Removing its constant term gives the convention requiring a nonempty tournament.

**FINITE-EXACT, all orders.** Among attainable odd heights at most 1503, precisely the following have a single smallest ordinary isomorphism class:
\[
  1,3,5,11,13,41,61,85,87,125,127,147,155,157,171,175,189,629,637.
\]
The distinction from fork multiplicity is visible at 155. There is one ordinary tournament of smallest order seven with 155 paths, while the strong fork centers of height 465 occur in a converse pair at order eleven. For any input height the [class-series replay](capped_cores.py) prints the ordered-factorization numerator terms and denominators. At 117, for example, the complete ordinary series is
\[
  \frac{7z^7+15z^8+5z^9}{(1-z)^2}+\frac{8z^{10}}{(1-z)^3}
       +\frac{3z^{11}}{(1-z)^4}.
\]

**FINITE-EXACT, all orders.** Among odd \(5\le h\le1501\) for which both \(h,h+2\) are attainable, there is a one-arc reversal between them except precisely at
\[
        h=35,39,53,55,59,63,83,93,107,147,187,195,261 .
\]
In particular all 620 consecutive steps from 263 to 1503 exist individually. Locality puts every such step into a single strong component, so the capped strong census decides the ordinary all-order question. Its positive-step counts also give the series for ordinary classes of height \(h\) with some reversal to \(h+2\), by multiplying their strong polynomial by \((1-z)^{-2}\). For example the census gives eight strong centers of height 263 with an upper exit. No strong class of height 263 has both a lower and an upper exit.

~~~sh
python3 -B capped_cores.py 1503 < /tmp/capped-fibres
python3 -B capped_cores.py 1503 117 < /tmp/capped-fibres
~~~

## A tournament chosen from the even target

For every even \(264\le N\le1502\) the [certificate table](target_flips.rows) gives the lower tournament and an arc \(u\to v\) whose reversal takes its height from \(N-1\) to \(N+1\). Each row is \(N,u,v\) followed by all outgoing adjacency masks in zero-based vertex order. For the perfect number 496 the order-eight lower rows are \((6,76,104,193,15,27,177,55)\). Reversing \(4\to2\) changes rows 4 and 2 to 11 and 120 and takes the count from 495 to 497. To select a target, the [standalone Python reader](target_flip.py) needs only the checked-in row table. It validates the tournament encoding and the arc, constructs both row matrices, and counts the Hamilton paths in each one.

~~~sh
python3 -B target_flip.py 496
python3 -B target_flip.py --all
# Regenerate the complete table from the height-capped augmentation.
# Progress goes to stderr and the sorted certificate rows go to stdout.
/tmp/capped-cores 1503 15 --edges 264 1502 > /tmp/target_flips.rows
cmp target_flips.rows /tmp/target_flips.rows
~~~

Of the 620 midpoints, 150 first occur at order eight, 469 at order nine, and only 312 first occurs at order ten. The generator considers strong centers in increasing order; core locality ensures an ordinary witness cannot appear sooner. Adding a singleton source preserves either endpoint count and their marked arc. Thus ten is the least common order at which every midpoint in the interval has an exact single-arc witness.
