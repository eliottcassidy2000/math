# Two apex contacts turn any even displacement into a literal equal-order pair

The underlying one-apex count F is inherited. The earlier all-even-gap results (Sessions199/221, recalled in [Session243](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session243_prime_gap_bridges_20261110/README.md)) are not claims of novelty here. The new restricted realization uses only a transitive deletion, **two changed contacts**, positive equal-order starts, an explicit sparse absolute inverse for the unit gap, and an exact (4^j±1) dimension line. The [older repunit degree law THM-4749](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-4749-the-repunit-degree-law-phi-2-of-2-is-the-triangle-phi-3-of-2-is-the-gap-and-fermat-products-reverse-the-order.md) already studies (2^k±1) and specifically the values 63/65, but the statistic there is multiplicity and component factorization, not this literal two-contact line. The archived [odd-square triangle dynamics](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/odd_square_20260921_triangles.md) uses 63 and 65 as conductors of modular six-cycles, another exact but different observable.

## Physical geometry and count

Take a transitive spine and add apex v. For positive integers u,w, the left entries of its out-tag word are the standard binary digits of u (most significant first), and the right are the reversed bitwise complement of the standard digits of w. Bit 1 means v beats that vertex. For zero omit that flank. The first vertex on a positive left flank is beaten by v and the last vertex on a positive right flank beats v, forming a directed spanning cycle through the original increasing spine whenever both flanks are positive. In particular these tournaments are strong.

Write (F(0)=1, F(2q+1)=F(q), F(2q)=F(q)+q+[q odd]) for positive argument. Then
    H(u,w) = uw + F(u) + F(w) - 1 + [u,w both odd].
For a direct proof, all directed odd cycles share the apex. The odd cycle packing expansion therefore is 1 plus twice their number (with no higher packing terms). Each cycle picks an earlier tag-1 and a later tag-0, then an even-sized subset between them. The weight after multiplying by 2 is (2^d+[d=0]) where d is the number of intervening vertices. On one individual flank these total (F(u)-1) or (F(w)-1); the binary recurrence follows by adding a least-significant tag, whose extra sum at a zero tag is the prefix value plus its final bit. The right flank gives the same law by reverse complement. Across the two flanks the power separates into the product of the digit place values; summing yields uw. The adjacent-pair bonus survives just when both least significant bits are 1.

## Independent scales versus fixed displacement

For arbitrary integers a,b≥1 consider the pair of addresses

    lower-label=(2a,2b+1),    upper-label=(2a+1,2b).

“Lower-label” names the construction rather than presupposing which output is smaller. Both widths stay unchanged and only their least significant bits move: exactly two ordinary apex arcs change, one in each flank. Both tournaments have the common framed Hamilton path v followed by the spine in increasing order; the flipped arcs are chords in that frame. Thus the two framed tilings differ in two entries. This is **not** asserted to be an edge of the actual tiling/converse metagraph.

Directly using the recursion,

    L(a,b)=4ab+3a+[a odd]+F(a)+F(b)-1,
    U(a,b)=4ab+3b+[b odd]+F(a)+F(b)-1,
    U−L=3(b−a)+[b odd]−[a odd].

The common substrate is a symmetric product, and the imbalance is read only by linear local tags. This decouples absolute scale from displacement partially: a,b can grow by the same amount with their difference and parity fixed. For any signed even δ, if δ≡0 mod6 take d=δ/3 even and b=a+d. If δ≡2 mod6 take d=(δ+1)/3 odd and a odd. If δ≡4 mod6 take d=(δ−1)/3 odd and a even. In every case choose a=max(1,1-d) or a larger value with the appropriate parity: both output tournaments are strong, have equal order, and their difference is precisely δ. This works at infinitely many increasing orders, but it does not choose either absolute endpoint at will.

The two edits interact: applying just the u edit first changes H by (2b+2-a-[a odd]); applying just the w edit first changes H by (b+[b\text{ odd}]-2a). The mixed increment is always −2 from the product and parity term. Thus the fixed small gap is cancellation between two potentially large single-arc changes, not a hidden single arc of that small strength.

## The exact unit inverse and dimension line

For δ=+2 the parameter is uniquely b=a+1 with a odd; L is

    L(a,a+1)=4a²+7a+F(a)+F(a+1).

For positive odd a, (1≤F(a)≤(a+1)/2), and (1≤F(a+1)≤a+2), by elementary induction on the recursion. Consequently

    (2a+1)² < L(a,a+1) < U(a,a+1) < (2a+3)².

The intervals for different odd a are disjoint. Inverting either endpoint at a prescribed N therefore takes a single candidate (a = floor((isqrt(N)-1)/2)): require positive odd a, compute L and U and test equality. The lower endpoints through X number (sqrt X/4+O(1)), since (L=4a²+O(a)) and a is odd; the union of both endpoints has (sqrt X/2+O(1)) values. The order of each witness is (3+bitlen(a)+bitlen(a+1)=log₂ H+O(1)). This sparse family cannot by itself give odd-count coverage.

In particular choose (a=2^k-1,b=2^k), (k≥1). Then (F(a)=1), (F(b)=b+1), and

    (H_lower,H_upper)=(4^(k+1)−1,4^(k+1)+1), at physical order 2k+4.

Increasing k by one increases each of the two dimensions by one vertex and multiplies the common center by four. The two affine terminal corrections remain −1 and +1 exactly. At (k=2), (63↔65) has literal rows

    (126,124,248,240,96,64,128,51),
    (126,124,120,112,96,64,128,63).

The first identifies by isomorphism with the second order-8 count-63 birth orphan from the [previous atlas](https://github.com/eliottcassidy2000/math/tree/68836af2f93be6da276c2af03fe690939297c1dd/05-knowledge/results/chanpak_birth_atlas_metagraph_20261117), rows (8,1,3,198,15,31,55,119). It does not identify with the other count-63 class. This puts a particular exceptional **host**, not merely the number 63, on the infinite dimension line. This observation gives no counterpart for 7 or 21 and no proof of Chan–Pak.

Run `python gap.py` for literal DP controls; `python gap.py --gap -14 --step 4` returns a signed-gap pair and literal rows, and `python gap.py --target 63` runs the exact +2-family inverse.

## Equivalent adjacent-center experiment

On the power line the two entire tag words simplify to `1^k 00 1^k 0` and `1^(2k+2)0`. More generally choose independently positive dimensions p,q, and compare `1^p 00 1^q 0` against `1^(p+q+2)0`; zero positions are adjacent, and the order is p+q+4. Odd cycles with the final zero lose weights (2^{q+1}+2^q=3·2^q) when the two center flags change from 1 to 0. Odd cycles ending at the new first zero gain total (2^p) from the p prefix starts, and those ending at the second gain (2(2^p−1)). Thus the unswitched count is (2^{p+q+2}+1), and the switched count is (2^{p+q+2}+3(2^p-2^q)-1). Equality of the two physical dimensions cancels the two exponential contributions exactly, so the two-contact loss is 2 precisely on p=q. This is the literal physical power-balancing mechanism; nearby unequal dimensions do not give that unit repair.
