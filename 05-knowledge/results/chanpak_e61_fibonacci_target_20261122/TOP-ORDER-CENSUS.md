# A complete inverse for the two top-order receivers

**FINITE-EXACT.** For every odd \(11\le p\le37\), no normalized construction from the inherited six-vertex core of order \(2p+1\) has the single-arc counts \(E_p-1,E_p+1\), where \(E_p=2^{p-1}(2^p-1)\). For the rooted four-vertex receiver the exclusion extends through \(p=45\), again for every odd exponent. The preceding [receiver bound and mod-three obstruction](SUFFIX-INVERSE.md) prove that the singleton and this four-vertex receiver exhaust the possibilities at that order. At exponent 31 there are exact centers at [orders 60 and 62](../chanpak_e31_physical_20261122/README.md). The unrestricted exponent-61 problem remains **OPEN**.

## Splitting at a fixed width

The singleton address has \(2p-7\) bits. The four-vertex address has \(2p-10\) bits and ends in five ones. Write it as \(y=32u+31\). Then \(f(y)=f(u)\), the new slope is \(c'=32\cdot1021\), and both intercepts become \(940+31\cdot1021\). For the singleton put \(u=y\) and \(c'=81\). Set \(u=Bz+r\), with \(B=2^b\), and apply the two-ended identity of [SUFFIX-INVERSE.md](SUFFIX-INVERSE.md) to the reduced address \(u\).

Each of the \(B\) low words gives one explicit interval for \(z\) and a unique residue for \(z\bmod h\), since \(\gcd(c'B,h)=1\) for both charts. Intersect the interval with the required high-word width, then evaluate every surviving quotient. If \(N=E_p-1\), the interval length is at most \(Nh/(c'B)^2\). Thus at most \(1+N/(c'B)^2\) quotients per low word need evaluation after the congruence. For exact inversion take the least power of two \(B\) with \((c'B)^2\ge N\), limited by the available word width. A table of 16-bit words supplies both \(f\) and bit reversal. Each solution is replayed against the original uncompressed recurrence.

The complete runs below use every low word. The entries are exact counts of quotients evaluated. Every entry returned zero addresses.

| \(p\) | Singleton | Four vertices | \(p\) | Singleton | Four vertices |
|---:|---:|---:|---:|---:|---:|
| 11 | 6 | 0 | 29 | 2,000,913 | 6,347 |
| 13 | 32 | 0 | 31 | 8,003,884 | 24,981 |
| 15 | 124 | 0 | 33 | 32,021,270 | 101,186 |
| 17 | 487 | 0 | 35 | 128,055,399 | 404,103 |
| 19 | 1,976 | 4 | 37 | 512,222,183 | 1,615,522 |
| 21 | 7,859 | 31 | 39 | — | 6,460,726 |
| 23 | 31,212 | 99 | 41 | — | 25,850,446 |
| 25 | 124,882 | 404 | 43 | — | 103,388,715 |
| 27 | 499,941 | 1,539 | 45 | — | 413,530,483 |

At exponent 61 the exact inverse calls for \(2^{55}\) singleton low words or \(2^{46}\) reduced four-vertex low words. Shards partition the work by the low-word residue. The small exponent runs establish how often the top order alone can miss the target despite allowing every suffix.

## How close an actual edge comes

For every odd exponent from 11 to 31 we also enumerated all addresses with \(\lvert H_{\rm lo}-(E_p-1)\rvert\le512\) in these **two** charts. To run an interval query, expand both quotient bounds and test every integer between them. Requiring the five terminal ones in the four-vertex chart remains valid throughout this window because \(E_p-1+512<2^{2p-1}\); the proof bounding \(f(y)\) applies verbatim. Other receiver types can become admissible away from the exact center, so this table describes these two charts only.

For an edge with counts \(E_p+\delta-1,E_p+\delta+1\), the table gives the closest attained \(\delta\) on each side. A dash means no address within the tested distance.

| \(p\) | Below | Above | \(p\) | Below | Above |
|---:|---:|---:|---:|---:|---:|
| 11 | −50 | +46 | 23 | −50 | +168 |
| 13 | — | +2 | 25 | −410 | +160 |
| 15 | −174 | +4 | 27 | −18 | +114 |
| 17 | −36 | +152 | 29 | −76 | +62 |
| 19 | −78 | +8 | 31 | −42 | +60 |
| 21 | −10 | +88 | | | |

At \(p=13\), singleton address 393899 gives the strong 27-vertex pair **33,550,337 and 33,550,339**. Its lower endpoint is \(E_{13}+1\). At \(p=29\), four-vertex address 141016812830847 gives the edge centered at \(E_{29}+62\). At \(p=31\), singleton addresses 26085436658688557 and 25984452563041549 give the two closest centers shown. In each example the actual two tournaments constructed by [check.py](check.py) are strong and differ at the marked arc. Order-63 failure at \(p=31\) shows why searching other orders matters even when the desired edge is attainable.

From this directory reproduce the exact searches or the nearest edges with the [standalone inverse](top_order_inverse.cpp).

~~~sh
c++ -O3 -std=c++17 top_order_inverse.cpp -o /tmp/top_order_inverse
for p in $(seq 11 2 37); do /tmp/top_order_inverse "$p"; done
for p in $(seq 39 2 45); do /tmp/top_order_inverse "$p" --receiver 4; done
for p in $(seq 11 2 31); do /tmp/top_order_inverse "$p" --radius 512; done
~~~

The optional flag `--shard I M` visits \(r\equiv I\bmod M\). Every invocation checks exact and interval queries against a literal full-width search at 10, 13 and 16 bits, including planted solutions.
