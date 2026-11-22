# A fixed ordinary height has a finite catalogue of auxiliary rays

The second-slot read \(Q=Z_T(2)\) is a different auxiliary count
from the preceding-prefix register of the backchain. Its possible
values at a fixed ordinary count can nevertheless be described
exactly using the smaller strong parts of a tournament.

For a positive integer \(N\), normalize a tournament of height
\(N\) by deleting every singleton strongly connected component.
The surviving strong components are ordered and of order at least
three. If their heights are \(h_i\), then \(\prod_i h_i=N\).
By [Busch's strong minimum](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3),
a strong component of order \(n\ge3\) has at least \(3^{n/3}\)
paths: for \(n=3\) directly, and for \(n\ge4\) use
\(5^{(n-1)/3}\ge3^{n/3}\).
Thus every normalized host has total order at most
\(\lfloor3\log_3 N\rfloor\).

Ordered \(q\)-slot path covers multiply over the ordinal sum.
Indeed, each slot splits uniquely into its ascending sequence of
component subpaths, and conversely those internal slots concatenate.
The deleted singleton components each have \(q\) choices. In
particular for every \(q\ge1\), with a common padding exponent \(k\),
\[
 Z_T(q)=q^k\prod_i Z_{C_i}(q).
\]
Conversely all \(k\ge0\) are realized by adding ordinal singletons
to any given nonempty normalized host. If \(N=1\), the normalized
host is empty and a nonempty tournament needs \(k\ge1\).

**Exact fibre consequence.** For fixed \(N\), the set of compatible
second reads is a finite union of rays \(Q_0 2^k\); its joint
second- and third-read set is the union of the synchronized rays
\((Q_0 2^k,R_0 3^k)\). Having a permissive congruence at height \(N\)
does not make an arbitrary prescribed second read physically possible.
The allowed odd parts of \(Q\) are a finite set. For each normalized
profile, the second read determines the padding \(k\) uniquely,
and thus fixes its third-read candidate.

This gives an exact divisor-graph compiler. For every divisor \(d>1\)
of \(N\), build the finite catalogue of profiles of strong tournaments
with \(H=d\) by labelled enumeration up to the strong-order bound.
Convolve these catalogues on the divisor DAG for \(N\), multiplying
the reads entrywise. A catalogue element at \(d=N\) is a genuine
primitive strong piece, not a step to a smaller height. Every other
factorization uses strictly smaller ordinary counts. For a prime
target there can be no decomposition into smaller nontrivial
strong components; a useful recursive method for prime targets
must split a strong component across its internal interfaces,
not merely identify its ordinal factors.

## Literal small fibres, all orders

For the following heights the exact finite normalized base profiles
are as follows. Every displayed pair \((Q_0,R_0)\) gives the full
padding ray above, and none are omitted.
\[
\begin{array}{c|l}
N & (Q_0,R_0)\\\hline
3&(12,33)\\
5&(32,117)\\
9&(84,411),\ (144,1089)\\
11&(100,465)\\
13&(104,471)
\end{array}
\]
The height-nine second pair is the product of two triangles.
For these heights no strong component has order six or more, since
\(5^{5/3}>13\). The attached check enumerates all orientations through
order five and counts path covers by independent subset endpoint DP.
