# Twenty-six physical ways to read p and 7p from the same source

This widens the [first rotational affine pair](DIHEDRAL_HEIGHT.md). For every row in the accompanying table and every actual incoming tournament X with h=H(X), B=f_X(2), there are literal outputs
```
p = alpha B + beta h,       7p.
```
The upper output has an actual seven-rotation, not merely a count divisible by seven.

For the upper, insert the seven-carousel at the marked place c and X at u. The other quotient vertices remain singletons. The ordered nonempty path-cover counts of the carousel for one through five runs are 175,1638,6678,15120,20160, all divisible by seven via the free rotation. The two quotient marks alone are transitive, so a path has at most n−1 runs of either mark. Retain only quotients allowing at most two X-runs. An exact run tally gives normalized weights r1,r2 for one and two X-runs, hence slope r2 and intercept r1−2r2. For the lower use the familiar single-apex signed wrapper and the literal contact set in the table. Match intercepts and slopes up to powers of two, which are aligned by ordinal padding of X on either side.

The finite discovery scope is all ordered pairs of marks in quotients of orders 3–6 and all nonempty contact sets in lower tournaments through order 8. The upper has 87 coprime (odd slope part, intercept) keys; 26 match the lower bank. Independent regeneration used induced restrictions of the archived order-7 representatives for small quotients, and all one-vertex extensions at order 8. Every table row is replayed from the literal masks; the census is only a means of obtaining the table, not a prerequisite for using the constructions.

## Exact nontransitive input selection

For P_m, reverse the longest arc of a transitive tournament of order m+2. Its source reads are h=1+2^m, B=4·2^m+2·3^m−2, as in [Prime-safe phase](PRIME_SAFE_PHASE.md). A prescribed target p on any row can be solved by testing whether (p−beta h)/(alpha B) is a positive power of two, whose exponent k is the shared ordinal padding. This is a **complete inverse within the two-base family**: only alpha·3^m≤p need be tested. Indeed beta>−alpha for every lower wrapper: its one-run template count before subtracting 2G is h(D)+G+Q, so beta=h(D)+Q−G>−G. Also B≥2h. Thus p>alpha·2^k·B/2>alpha·2^k·3^m. The search uses O(log p) checks per row and includes valuation collisions. On the valuation-separated branch k+v_2(alpha B)>m, the valuation of p−beta already determines m for odd beta.

New absolute placements include these primes (complete trial division in replay):

| (alpha,beta) | (m,k) | p | 7p |
|---|---:|---:|---:|
| (284,25) | (10,7) | 4,441,949,209 | 31,093,644,463 |
| (852,−25) | (10,8) | 26,651,515,879 | 186,560,611,153 |
| (1111,593) | (8,3) | 125,864,273 | 881,049,911 |
| (1729,−877) | (8,3) | 195,414,419 | 1,367,900,933 |
| (1554,−577) | (10,6) | 12,152,136,127 | 85,064,952,889 |

For these cases, the lower automorphism group is trivial and the upper is exactly C7: both group orders divide the associated path count and p is greater than either vertex order. A reflection on the entire upper graph is not asserted.

## Arbitrarily deep 3-adic and 5-adic phase with genuine seven-rotation

For **any** i,j,L≥1 and unit target t modulo 2^L·3^i·5^j, the bank's row (5400,−877) gives actual same-source witnesses p' and 7p' with p'≡t modulo that modulus. This strengthens the earlier phase construction by preserving an actual seven-rotation and treating arbitrary depths at all three primes.

First reserve r=max(0,ceil((j−2)/2)) copies of the count-11 core, whose height is 100. This will multiply h by 11^r, but makes 5400B vanish modulo 5^j because 5400 is a multiple of 25. Append L ordinal singletons to make the height term vanish modulo 2^L. Use the literal [Prime-safe phase](PRIME_SAFE_PHASE.md) program to match the binary target and the sign modulo 3 for h, after accounting for multiplication by −877·11^r.

To set the residue modulo 5^j, use a factor h=1+3·2^c with c≥L and c≡5 mod 20. It is 22 modulo 25; its reduction modulo 5 has order 4 and its fourth power is not 1 modulo 25. Thus it generates all units modulo 5^j. Join the required number of these literal blocks before adding the reserved count-11 cores and singletons. This does not change the binary phase or sign modulo 3.

Finally choose a≥L with a≡1 mod 4 and a≥2, and use the block with h=1+3·2^a and height B given in the phase note. Its order is a+3, which is even, so B≡h≡1 modulo 3. Join e=4·5^(j−1) copies to make one repeatable correction block (h_*,B_*). It satisfies h_*≡1 modulo 5^j and 2^L, and v_3(h_*−1)=1 since 3 does not divide e. Appending n such correction blocks leaves both previously programmed coordinates alone. Its effect at 3 is
```
p_n = 5400 B_0 B_*^n − 877 h_0 h_*^n.
```
For any n,r≥0, the difference p_(n+3^r)−p_n has 3-adic valuation exactly r+1. The second term has that valuation by the lifting-the-exponent identity; the first has valuation at least r+4 because B_*≡1 modulo 3 and 27 divides 5400. So the three choices n,n+3^r,n+2·3^r have distinct next triadic digits, while preserving all earlier ones. This programs the requested exponent i in place, even though the height term no longer vanishes modulo 3^i in general.

This is constructive in all three coordinates at once, but the number of copies can grow quickly with the requested precision. It does not claim absolute coverage, prime outputs in every congruence class, or a small graph for a very deep phase request. The sharper gap is now explicit: even paired witnesses with genuine group action can freely hit any of these phases, while placing a general prime requires coordinating h and B at a prescribed absolute height.

Run `python3 paired_bank_replay.py` with the table and earlier packet scripts.
