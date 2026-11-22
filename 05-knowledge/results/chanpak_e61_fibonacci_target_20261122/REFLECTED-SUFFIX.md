# Exact inversion along reflected halves of a Fibonacci suffix

**PROVED.** Fix a receiver with lower chart \(L(y)=K_{y\bmod2}+cy+h f(y)\), \(c>h>0\). Split a positive address into a leading \(b\)-bit word \(s\), a central word of width \(g\) with value \(m\), and a trailing \(b\)-bit word \(r\). Put \(B=2^b,\ A=2^g,\ X=As+m,\ v=\operatorname{rev}_b(r)\). Specifying the signed integer difference \(v-s\) or \(v-(B-1-s)\) leaves only a bounded number of candidates to test for any prescribed absolute height. At fixed chart this bound is independent of \(b,g,m\). Thus fixed reflected or reverse-complemented halves, and arbitrary carries between them described by one signed integer, have an exact inverse without enumerating the halves. [Original receiver and split identity](README.md)

Write \(v=\epsilon s+d\), where \(\epsilon\in\{-1,1\}\). The splitting identity yields
\[
 L(y)=K_{r\bmod2}+G(s)+c r+h\bigl(f(X)+f(r)+[r<B/2]\,[X\text{ odd}]\bigr),
 \qquad G(s)=X\big[(c+h)B-h(1+v)\big].
\]
The constraints \(B/2\le s<B\) and \(0\le v<B\) give a single integer interval. Using \(0\le f(u)\le u\), a solution at \(N\) must have
\[
 0\le N-\min_j K_j-G(s)
 \le (c+h)(B-1)+hAB+|K_1-K_0|.
\]
For neighboring feasible values of \(s\), the increase in \(G\) is at least \(A(c-h)B\). Indeed the second factor is at least \(cB\), and its per-step decrease is at most \(h\) with \(X<AB\). The display therefore selects at most one plus the floor of the window width divided by \(A(c-h)B\) candidates. Intersect with the feasible interval and substitute their actual words in the original formula. The implementation estimates a quadratic root, corrects it with exact integer comparisons, and works backwards only through the certified window. For a fixed offset it evaluates all solutions. Across different offsets or signs distinct presentations can overlap.

**FINITE-EXACT at exponent 61.** For each of the two charts below the search covered both signs, all middle bits (one for the singleton and none for the four-vertex receiver), and all reflected half differences whose signed nonadjacent binary expansion has at most five nonzero digits. Both searches were empty.

| Receiver | Address width | Distinct offsets per sign | Hits |
|---|---:|---:|---:|
| singleton | 115 | 101,995,775 | 0 |
| rows \((2,12,9,1)\) | 112 | 92,581,247 | 0 |

Enumerate all sums of at most five signed powers of two on nonadjacent positions \(0,\dots,b\), retaining \(|\delta|<B\). Use \(d=\delta\) for direct reflection and \(d=B-1+\delta\) for complemented reflection. Every integer with absolute value below \(B\) has its unique nonadjacent form with highest bit at most \(b\), by the usual carry recursion at odd residues modulo four. Two reflected strings disagreeing in at most five positions always fall in this test: their numerical difference has a signed binary expansion with at most five digits, and conversion to nonadjacent form cannot increase this minimum weight. The actual test also includes differences with long carry chains and arbitrarily many runs.

Run `c++ -O3 -std=c++17 folded_suffix.cpp -o /tmp/e61-folded && /tmp/e61-folded 5` from this directory. The small controls compare complete fixed-difference inverses to direct counts for several middle widths and every pair of halves through five bits. A full unrestricted inverse could enumerate every offset in one sign, but there are exponentially many more offsets than in the reported computation.

## Why a distance multiset also gives an exact inverse

For the marked positions of the ordinary binary address (numbered from the low bit), let \(m_d\) count unordered pairs of marks at distance \(d>0\). Then
\[
 f(y)=y-[y\text{ odd}]-D(y),\qquad
 D(y)=\sum_{d>0}\bigl(2^{d-1}+[d=1]\bigr)m_d .
\]
In the recurrence defining \(f\), a one in position \(i\) contributes weight \(2^{i-j-1}+[i=j+1]\) for each zero at \(j<i\). If all lower positions were zero, the sum would be \(2^i\) for \(i>0\) and zero for \(i=0\). Subtracting pairs of actual marks gives the formula. It holds for every word, with or without leading zeros.

In particular all addresses with the same multiset of pairwise distances have the same \(D\). At prescribed \(N\) such a multiset gives **at most one address of each parity**, computed by exact division
\[
 y=(N-K_j+h j+hD)/(c+h).
\]
Check the parity, distance multiset and width on that one integer. This bypasses reconstruction of every point configuration with those distances, the usual turnpike ambiguity. With padded reversal \(v\) of the whole word the same identity gives \(f(y)-f(v)=y-v-[y\text{ odd}]+[v\text{ odd}]\). For the singleton \(K_1-K_0=h\), so the height is exactly \(43+106y-25D(y)\); two equal-distance words differ in height by exactly 106 times their numerical displacement. This selector solves each prescribed distance class. Which class can give \(E_{61}-1\) remains open. [Earlier turnpike and receiver connection](../chanpak_reader_retention_20261115/COMPENSATING_READER.md)
