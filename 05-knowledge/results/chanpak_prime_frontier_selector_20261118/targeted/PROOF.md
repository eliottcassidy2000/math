# Explicit strong witnesses immediately beyond the two audited/reported caps

Prior scope is important. The [Session203 corrected independent replay](https://github.com/eliottcassidy2000/math/blob/main/session203/coverage-recertification/README.md), blob adc288964a14c472f3e9ec05bff7deb633850971, covers every odd through (2^{32}-1=4,294,967,295) except 7 and 21. The [THM-4885 addendum](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-4885-one-apex-join-representation-mass-pi2-over-12.md) reports a separate single bucketed run through (2^{35}-1=34,359,738,367); its full height was not independently recomputed. I did not find a higher contiguous old certificate. Neither cap identifies a counterexample.

This packet checks explicit strong witnesses for the **next 1,000 odd targets** beyond each cap:

* 4,294,967,297…4,294,969,295, all by F1, extending the independently replayed lower cap by just this finite band;
* 34,359,738,369…34,359,740,367, 999 by F1 and one by F2. The full prefix to this new high endpoint inherits the single-run evidence below it; the new band itself is a literal per-value certificate.

The first prime and 7-times-prime axes are:

| reference cap | target | factor type | F1 b,e | order |
|---|---:|---|---|---:|
| 2^32−1 | 4,294,967,299 | 7·613,566,757 | 2,1502296101 | 34 |
| 2^32−1 | 4,294,967,311 | prime | 2,1847524150 | 34 |
| 2^35−1 | 34,359,738,371 | 7·4,908,534,053 | 9,2371737543 | 38 |
| 2^35−1 | 34,359,738,421 | prime | 71,297178177 | 38 |

Minimality of the displayed prime *targets beyond each old cap*, and primality of the indicated quotients, are checked by trial division. No claim of minimum witness order is made.

## Literal orientation and exact count

Write C(u) for the one-apex tournament whose old vertices are (0,1,...,\operatorname{bitlen}(u)-1), increasing by all arcs. The root loses precisely to the old vertex i when bit i of u is 1. If u=0 the canonical block has just its root. Its Hamilton count F(u) obeys F(0)=1, F(2z+1)=F(z), F(2z)=F(z)+z+[z odd] for positive even argument. The usual explicit one-apex word sum counts each inverted tag pair at distance d with weight (2^{d-1}+[d=1]); this independently verifies scalar F for every witness.

For F1, left block is C(b−1), right block is the converse of C(e−1), and **all interblock arcs point from right to left except the one arc from left root x to right root y**. A path avoiding x→y is ordinal right→left, and the only omitted possibility switches directly from y to x. A root-start path in C(u) exists uniquely iff u is even, since the remaining transitive path must start at vertex 0, or be empty. Thus paths avoiding the bridge number (F(b−1)F(e−1)-[b,e\text{ odd}]).

A path using x→y has forced order of runs (right optional, left ending x, right beginning y, left optional). The independent ordered rooted two-cover counts are b and e: in a C(u) cover with first slot ending at the root, either its old prefix is empty or its maximal old index i has bit1, with any choice of the i smaller indices. So the count is (b e+F(b−1)F(e−1)-[b,e\text{ odd}]). The replay materializes row orientations and independently tests strong connectivity. No (2^{34}) or (2^{38}) subset enumeration is claimed; the separator argument computes the exact count.

The single F1 miss in the second band was 34,359,739,333 (prime). Its F2 witness uses left C(v) for v=1,012,954,553, right rows (2,4,9,3,15) rooted at vertex4, and the same right→left orientation except the left-root→right-root bridge. The right profile is H=5, root-end count 0, and rooted two-cover count with second slot beginning at root (R=32), checked by ordinary path subset DP on the five-vertex right bank. Since F(v)=389,038,721, the exact result is (5F(v)+32(v+1)=34,359,739,333) at order36. This repairs the **failed template**, not a gap in the tournament spectrum.

## Targeted inverse and replay

For F1, exact N implies (b e≤N), so order the choices with (b≤e) and search (b≤\sqrt N). For a fixed b put (A=F(b−1)), (e−1=Lq+r) with L=2^k, k≥1 and 0≤r<L, (Z_r) the k-bit reverse-complement of r, and ε=[r<L/2]. The padded suffix identity gives (F(e−1)=F(q)+qZ_r+F(r)-1+ε[q\text{ odd}]). Hence exact N is equivalent to

    C = Bq + A(F(q)+ε[q odd]),
    B = bL + A Z_r,
    C = N − b(r+1) − A(F(r)−1) + (b mod2)[r even].

The bound (1≤F(q)+ε[q\text{ odd}]≤q+2) forces
(\max(0,\lceil(C−2A)/(B+A)\rceil)≤q≤\lfloor(C−A)/B\rfloor).
Choose L as the first power of two at least 2 with (L^2b≥N), iterate r, and verify every surviving candidate exactly. For F2 the same bounds apply with coefficients replaced by its fixed rooted right profile. All arithmetic at the two tested bands fits 64-bit signed integers. The point of the inverse is that the full billions-range census is not needed to land at a specified target.

Run `python replay.py` for all 2,000 physical controls and prime checks. To regenerate F1 candidates run `g++ -O3 search.cpp -o /tmp/frontier_search` and pass the desired explicit odd targets. The F2 inverse `search_f2.cpp` takes the [inherited rooted bank](https://github.com/eliottcassidy2000/math/blob/main/session203/coverage-recertification/session203-independent-audit/aud_bank7_inherited.txt) then the target. These certificates establish only the stated finite bands, not universal selection.
