# The seven odd slots after 9 as an independently liftable prime wheel

The odd interior of \((9,25)\) reads \(11,13,15,17,19,21,23\), giving **P,P,C,P,P,C,P**. At display coordinate \(x\) the shifted tournament question asks for height \(x-2\). Thus height 7 lies under the lower square 9, and height 21 under the last slot 23. The square 9 is the first odd square greater than one. The [inherited bracket theorem](../procgen_numerology_20260923_odd_square_brackets.md) concerns **primes** sharing an interval between consecutive odd squares with a larger positive integer multiple. It suffices to test their doubles; the list is \(2,3,11\). With the same strict interval and no primality assumption, the complete list is 2, 3, 4, 10, 11 and 12. Indeed the first two brackets give those six, and for each odd lower root at least 5 twice its square exceeds the next odd square. [Integer multiplier classification](../session74_eleven_mechanisms_20261027/session74-inheritance/INHERITANCE.md) · [THM-1370, the actual height exclusions](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md)

## Three repeated square windows give actual integer realizations

For an odd root \(q\ge3\) write \(x=q^2\) and examine the first seven odd interior slots of \((q^2,(q+2)^2)\). Let
\[
 A=(2,4,8,10,14),\qquad S=9+A=(11,13,17,19,23).
\]
Having all five \(x+a\), \(a\in A\), prime forces \(3\mid q\): otherwise \(q^2+2\) is composite. Consequently \(q^2+6,q^2+12\) are both composite and \(q^2\equiv9\pmod{24}\). Modulo 5 the five prime slots exclude squares 0 and 1, leaving \(q\equiv\pm2\pmod5\). Thus \(q\equiv\pm3\pmod{30}\). When \(q>3\) the new first slot must also avoid the old prime 11, so \(q\not\equiv\pm3\pmod{11}\). Copying the original square phase to arbitrary precision at 11 already kills this first prime.

**FINITE-EXACT.** The first four roots with the actual prime pattern are 3, 32,397, 213,237 and 254,577. At each of the last three roots **every one of the nine odd heights from \(q^2-2\) through \(q^2+14\)** has an explicit two-apex tournament. This includes the two translated hole positions and every in-between integer of the correct parity.

| Root | Lowest height | Highest height |
|---:|---:|---:|
|32,397|1,049,565,607|1,049,565,623|
|213,237|45,470,018,167|45,470,018,183|
|254,577|64,809,448,927|64,809,448,943|

The last two intervals lie beyond \(2^{35}\). [Replay and all 27 address pairs](brackets.py). The finite proof sieves the possible roots, trial-divides the five candidates, and uses the existing literal one-apex/one-apex join at each successful height. For example the middle row's two translated hole positions have addresses \((2,19162820870)\) and \((2,21249748028)\). A separate search found **12 actual one-arc flips** inside these windows. In root order as in the table, the certified lower offsets for \(q^2+d\to q^2+d+2\) are respectively \(d=-2,0,4,6,8\), \(d=-2,0,4,8\), and \(d=0,6,12\). Each is a pair of strong tournaments on one vertex set and a single specified reversed arc. The twelve other pairs in the three windows have ordinary endpoint witnesses and their adjacency remains undecided here. The script constructs the respective boards from the tagged cores and their two flanks and checks the single edit by comparing the restricted outgoing masks. [Star-pair identity and verification method](../chanpak_gap_star_ports_20261122/PERFECT-BAND.md) [Construction and its arc-load recurrence](../chanpak_prime_frontier_selector_20261118/README.md) · [Earlier high prefix and finite bands](../chanpak_prime_frontier_selector_20261118/targeted/PROOF.md)

For unrestricted integer roots the five quadratics are irreducible and jointly admissible at every prime. The root 3 avoids a zero in their product modulo 2, 3, 5 and 7, and each prime at least 11 has more residue classes than the sum of the five degrees. Thus the polynomial prime-value conjecture predicts infinitely many such roots. This application concerns all integers as roots; its hypotheses give no prediction by themselves for roots restricted to a multiplicative tournament family. [Bateman and Horn, *A heuristic asymptotic formula concerning the distribution of prime numbers*](https://doi.org/10.1090/S0025-5718-1962-0148632-7)

## Exactly 24 of 32 patterns with the six tournament heights

Here the root \(q\) itself is a Hamilton-path count. Take finite ordered joins of the six strong tournaments of heights \(5,17,37,77,129,751\). Each has transitive odd-precedence tournament, the tournament directing a pair according to odd relative order in all Hamilton words. Ordinal join preserves this condition and multiplies heights. Hence \(q\) and \(q^2\) have rigid witnesses; \(q^2\) is the height under the first displayed slot. [Six-block construction](../chanpak_tiling_viewpoint_20261123/PRECEDENCE-FIBER.md)

Define a five-bit observation of the seven-slot interval. Require all five positions \(q^2+a\) to survive trial division by 2, 3, 5 and 7. For each \(p\in S\), the bit is **on** when \(p\) divides its own slot \(q^2+(p-9)\); when the bit is off, require \(p\) to divide none of the five slots. When on, it automatically misses the four other slots since their nonzero differences are even and less than \(2p\). A zero bit asserts only survival of these divisibility tests. At the original root 3, all five bits are on through equality with their primes. At every larger root an on-bit certifies a composite.

**PROVED.** The six-block family produces exactly **24** of the 32 possible masks. Its missing masks are precisely the eight containing both 13 and 17. Adding one further block type of height 41 gives all 32; at most one copy of 41 is needed. Seven is the smallest possible order for such an enlargement using transitive-precedence tournaments, even if arbitrarily many new types of smaller order are offered.

For the obstruction, small-prime survival forces \(q\) nonsquare modulo 5, and excludes the factors 5 and 77 (they kill slots 10 and 14 through 5 and 7, respectively). Marking 17 excludes the factor 17 since it requires \(q\equiv\pm3\pmod{17}\). On the three remaining generators \(37,129,751\), the quadratic characters modulo 5 and 13 agree. The 13 mark would require \(q\equiv\pm3\pmod{13}\), both of which are squares. This contradicts the character at 5. The full transitive-precedence census through order six has heights exactly \(1,5,9,17,37\); the additional square 9 preserves the character equality. At seven the canonical outgoing masks \((118,120,58,113,64,80,4)\) give height 41 with transitive precedence. Its characters at 5 and 13 are opposite.

Here is an explicit exponent decoder proving sufficiency. For a unit modulo \(p\) take its discrete logarithm modulo \((p-1)/2\); this identifies a residue with its negative and retains precisely the information in its square. Using primitive roots \((2,3,2,2,3,2,5)\), respectively, gives these small coordinates.

|Multiplier|5|7|11|13|17|19|23|
|---:|---:|---:|---:|---:|---:|---:|---:|
|129|0|1|3|0|3|2|10|
|37|1|2|2|1|1|0|10|
|751|0|2|3|4|1|8|6|
|41|0|0|3|1|3|4|1|
|17|1|1|4|2|zero residue|1|7|

Let the desired bits be \((b_{11},b_{13},b_{17},b_{19},b_{23})\). Use
\[
 q=129u\,37^i751^j.
\]
For the enlarged family take \(u=41^{b_{13}}\). Aim for coordinates at the five original primes, in order,
\[
 (1+2b_{11},\ 1+3b_{13},\ b_{17},\ 4b_{19},\ 5b_{23}).
\]
The on-coordinates give square 9 and the indicated original divisor. The off-coordinates avoid all five negative offsets. Write \(r_p\) for the coordinate of \(129u\). Choose \(j\equiv r_{19}-4b_{19}\pmod9\) and \(j\equiv b_{17}-r_{17}-1\pmod2\). Now determine \(i\) modulo \(8,3,5,11\) using the target coordinates at \(17,13,11,23\), respectively. At 13 its mod-two part already holds because \(i\) is odd. The coefficients of \(i\) in the four congruences are \(1,1,2,10\), all invertible. CRT supplies nonnegative exponents. At 5 the resulting coordinate is 1; at 7 it is \(b_{13}\), and both values there are safe. At 2 and 3 oddness and divisibility of \(q\) by 3 suffice.

For the original six-block family the same construction with \(u=1\) gives every mask with \(b_{13}=0\). If \(b_{13}=1,b_{17}=0\), choose \(u=17\). Its zero residue at 17 automatically misses all five slots. Retain the target residues at 11, 13, 19 and 23, take \(j\equiv r_{19}-4b_{19}\pmod9\), and replace the modulo-eight condition on \(i\) by evenness. The modulo-three condition now gives \(i+j\equiv2\pmod3\), and the resulting coordinate at 7 is 0. This covers the remaining eight masks and proves sharpness.

## One integer direction lifts all five local depths independently

**PROVED.** Within each attainable mask prescribe an arbitrary positive valuation \(t_p\) at each on-prime. All those valuations hold simultaneously at infinitely many actual rigid root heights, using the same block families as above. The other original primes and the small-prime conditions retain their prescribed readings.

Take a root \(q\) from the decoder and the fixed multiplier
\[
 L=(37^6\,751)^{720}.
\]
It is a product of available tournament heights and is 1 modulo each of \(5,7,11,13,17,19,23\). Its first lifting digits \((L-1)/p\bmod p\) at \(11,13,17,19,23\) are respectively \(10,10,15,12,12\), all nonzero. Thus at every on-prime the map \(k\mapsto L^{2k}\) generates the full principal-unit group modulo \(p^{t+1}\) for any \(t\ge1\) (binomial lifting). Since \(-a/q^2\equiv1\pmod p\), the condition \(p^t\mid(q L^k)^2+a\) occupies a single exponent class modulo \(p^{t-1}\), and exactly one of its \(p\) lifts is excluded for exact valuation \(t\). These independent exponent periods are powers of distinct primes. CRT constructs \(k\), with arbitrarily many positive continuations. One can choose one common order for a finite set of root witnesses by padding with singletons, and likewise for their squares.

This is one precise sense in which the first 25 give an unbounded construction. One small window encodes five locally independent integer depths. The local tests recognize forced composite slots; the occurrence of five *simultaneous primes* infinitely often at square roots of this form is open. The 11-vertex threshold for the 12 residues concerns the original odd-precedence fiber at a single order. Here the extra seven-vertex block solves a distinct obstruction to simultaneous directions, using two characters outside \(2,3,11\). [Eleven-vertex threshold and local generation](../chanpak_tiling_viewpoint_20261123/PRECEDENCE-FIBER.md)

Run `python3 -B brackets.py` from this directory for the modular decoders, the minimum-order input, the lifting digits, the three genuine repetitions and the 27 tournament realizations.
