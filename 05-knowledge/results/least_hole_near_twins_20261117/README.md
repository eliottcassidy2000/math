# Least-hole constraints and physical near twins

**No new tournament-spectrum counterexample was found.** A least hole cannot be forced from its predecessor by adding vertices. There is a stronger port-rewiring mechanism, but its exact inverse becomes a Pell-type conic, and it cannot land at three prescribed decimal-centred prime pairs. Different rewiring families do give repeatable pairs on both prime axes.

## 1. What a third hole would have to look like

The [corrected finite certificate](../../../session203/coverage-recertification/README.md) stops at 4,294,967,295. Any least additional hole N exceeds that bound. Closure under ordinal products and the seeds 49,63,343 imply N is either an odd prime other than 7 or 7p with odd prime p≠3,7. Such targets are indecomposable as products of two admissible odd numbers greater than one. If a witness exists it therefore has one nontrivial SCC, apart from ordinal singleton padding. [Busch's bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) gives a finite order limit 1+3 log₅N for that strong SCC. So an exclusion at an actual fixed target is finitely certifiable in principle. The classifications below exclude only specific grammars, not every strong tournament under that bound.

Multiplicativity by itself cannot make an additional hole propagate. The abstract set of all positive odds except 7,21 and a single prime p>21 is already product-closed and still has all parity-compatible residue classes and all positive even differences. The same is true if the extra hole is 7p instead of p. These are logical countermodels to a scalar argument from closure, not claims about the tournament spectrum. An actual fixed-congruence exclusion is already ruled out by the [one-apex residue placement theorem](../seven_gap_atlas_20261117/README.md): in any fixed normalized interior slice of the dyadic gaps it can realize every prescribed compatible residue at sufficiently large scales. A third hole would have to be an absolute landing failure, not a forbidden odd residue class.

There is a general small-arm bound for the two-core route. In H=ad+be−cf, both terminal ports on nontrivial cores obey c<a and f<d; hence be≤N−1. At least one weighted cover count is ≤√(N−1). For an ordinal decomposition of a pointed core, every earlier SCC contributes its two-slot count multiplicatively to the weighted end count b, since each slot concatenates its pieces in SCC order. The last SCC contributes its own marked b, at least its H. Every preceding layer of r vertices contributes at least 2^r (every subset has a Hamiltonian path in both slots). Busch then gives b≥5^((order−1)/3). Thus at least one core in any nontrivial two-core seam for N has order ≤1+(3/2) log₅(N−1). This reduces one side of the inverse, without proving a viable matching side.

### Why a predecessor induction needs real rewiring

By the [weighted odd-cycle packing identity](https://arxiv.org/abs/2307.05569), induced additions never decrease H. If adding a single vertex v gains exactly 2, there is precisely one odd pack using it: a single triangle v→u→w→v. Write U=N⁺(v), W=N⁻(v). There is only one forward arc U→W, namely u→w. The graphs U−u and W−w are transitive; any triangle in either would give a disjoint two-cycle pack together with the new triangle. There is no three-vertex path in U ending at u or in W starting at w, because either would produce an odd five-cycle through v. Thus u has at most one predecessor z inside U (if present, z is a source), and w at most one successor q inside W (if present, q is a sink). They cannot both be present, as v→z→u→w→q→v would be a five-cycle. Consequently U,W are transitive and the old graph's only possible cycle is a single triangle using u→w and the sole extra vertex z or q. Its H was 1 or 3.

For an induced extension by arbitrarily many vertices, temporarily add the vertices one at a time. All gains are even and nonnegative. A total gain 2 would have a single +2 step, whose old count must be 1 or 3. Therefore if N>5 is a least extra hole, the existence of an N−2 witness cannot deliver N by any induced extension. An exact successor proof must be allowed to re-encode or reverse arcs.

## 2. Odd-square boundary and a real two-arc near-twin switch

For a prime p, the phrase “has a multiple below the next odd square” is exact when the multiple is the least proper ordinary multiple 2p. In the band (2r−1)²≤p<(2r+1)², such a p must lie in the initial interval ending at floor(((2r+1)²−1)/2). Only the first two bands can have members, because 2(2r−1)²≥(2r+1)² for r≥3. Filtering those intervals for primes leaves precisely 2,3,11. The claim does not hold for all integers (4 is a counterexample to that reading), nor if the multiple is required to be odd. More generally replace 2 by a multiplier k>1; the candidate bands satisfy r<(√k+1)/(2(√k−1)).

A different, literal use of the square boundary: take the extreme-reversal cores R_j and R_(j−1)^op from the [two-core construction](../metagraph_two_core_20261117/README.md), and switch only the terminal port of R_j from 0 to 1. This reverses two cross arcs, and the two counts are

    L_j=2^(2j)+3·2^j+1,  L_j+2,       j≥2.

The switch is exactly +2 because the right cover/start pair is (2^(j−1)+1,1), while the left cover increases by 1 and its start count increases by 2^(j−1)−1. The pair lies strictly between consecutive odd squares (2^j+1)² and (2^j+3)². The lower never vanishes mod 7, and the upper does precisely when j≡0 mod 3. Exact pairs include (29,31), (89,7·13), (4289,7·613). Both members avoid divisibility by 3. This is repeatable neighboring realization, not a proof that primes keep occurring on the ray.

## 3. A larger +2 grammar meets an exact arithmetic obstruction

Let X and E be any nonempty tournaments, with (u,R)=(H(X),P₂(X)) and (v,S)=(H(E),P₂(E)); the two slots may be empty. Make a left cyclic block tournament C₃[0, TT₁→X, z], with 0→middle→z→0. Use port 0 or the initial singleton of the middle. The profiles are

    (a,a,u),  (a,a+u,R),     a=2R+u.

The first coordinate comes from the singleton-tail read on (middle→z) or the two-slot triangle law. At port 0 all two-part covers concatenate: the only bad possible starting vertex for the second run would be z, but a first run ending at 0 without z is singleton and z cannot start the remaining path. From 0, a Hamiltonian path traverses the whole middle and ends at z, with u choices. At the middle source, its only in-neighbor is 0, so the sole failed concatenations have first run singleton and second run starting 0, giving u failed joins. The start count is R, since the first middle source path traverses two optional X portions before z→0 and after it.

The right profile, read on the converse before gluing, comes from C₃[0,E,z] at port 0 and is (S+v,S+v,v). The lower count and change of port are

    N=2(2R+u)(S+v)−uv,       Δ=uS−vR+2uv.

In particular Δ=2 is the determinant-one condition u(S/2)−v((R−2u)/2)=1. This is the Farey-neighbor relation on **physical** profile coordinates. It forces coprime u,v. Once X and N are fixed, the other two counts are forced:

    v=(uN−4a)/(4R²−2Ru−3u²),      S=(v(R−2u)+2)/u.

There is a stronger inverse test without searching any core. Put K=uv and W=4vR−K+4. Eliminating S gives

    W²=13K²+(4N−24)K+16,
    1≤K≤N/29, K odd, W+K≡0 (mod 4).

The bound uses R≥2u and S≥2v, the contributions with one slot empty. Equivalently (13K+2N−12)²−13W²=4((N−6)²−52). Therefore the square or Pell test is a necessary preliminary to looking for a physical factor pair K=uv. It is not sufficient by itself.

There are **no** K,W passing that test for the lower endpoint 10,007, 30,011, or 1,000,000,007. The attached complete integer-square replay checks every odd K up to N/29 at each prescribed value, independently of tournament order, with multiplication checked at the square-root boundary. Controls 29,89,1253 do pass; 1253→1255 has the actual typed pair (u,R)=(1,8) and (v,S)=(5,32). The inherited [binary-apex selector](../target_core_selection_20261117/README.md) gives actual witnesses for both endpoints of all three excluded decimal-near pairs. Thus even an all-order exact +2 grammar can miss a whole pair of known witnesses. The failure is in the proposed transportation of marked counts, not an exclusion from the spectrum. In particular the billion twins also have order-20 strong witnesses via cyclic gluing. Their two-cover pairs are (1991,67712),(707,15880) and (1889,62576),(787,17216), respectively. The literal nine-vertex inputs (plus one ordinal singleton in the left arm) are in the replay; the cyclic read is (R−u)(S−v)+2uv.

## 4. Change the whole reader without changing its scale

A left R_k (with k+2 vertices) at the last natural port has (a,b,c)=(2^k+1,2^(k+1)−1,1). For right profile (d,e,f) the seam reads

    (d+2e)·2^k + d−e−f.

Thus a change of reader with equal d+2e moves only the intercept, exactly at every scale. Literal small right cores give:
- Profiles (9,21,1) and (11,20,2) share coefficient 51 and give 51·2^k−13 and 51·2^k−11. At k=5 these are the actual twins 1619 and 1621.
- Profiles (33,64,4) and (45,58,8) share coefficient 161 and give 7(23·2^k−5) and 7(23·2^k−3). At k=3 the prime quotients are 179 and 181. Every odd k avoids divisibility by 3 in both quotients, but primality across scales is unproved.

The prime-axis and sevenfold-axis examples use the same growing left family, not a separate fit per height. Yet their slopes and intercepts are fixed. The remaining obstruction is to select or compose readers so that the slope and intercept are forced by an arbitrary prescribed target, not just available for a previously selected ray.

A genuine decimal-scale comparison uses the five-vertex rows (24,1,3,6,14), with profile (H,P₂)=(11,100). Ordinally take j copies and use the singleton cyclic reader. It gives exactly 10^(2j)+11^j. This matches the decimal scale but not a constant offset. More generally, on a single stock (u,R) with 1<u<R, the difference of any two fixed cyclic linear readers is A R^j+B u^j. A nonzero such expression cannot stay equal to 2 at unbounded j, by dominance of the largest surviving exponential. Thus simply reusing the same ordinal stock cannot explain fixed decimal offsets; different cores or scale-dependent readers really are necessary in that model.

Replay: `python replay.py` and optionally `python -O replay.py`; conic: `g++ -O3 conic.cpp -o conic && ./conic`. The earlier [Hadamard lower-only switch cube](../hadamard_lower_cube_20261117/README.md) is another independent way to change a reader, with a demonstrated fixed-reader gap rather than arbitrary target placement.
