# Extra gaps, multiplicative ideals, and the actual 7·3^k ray

Let S be the set of Hamilton path counts of nonempty tournaments. Ordinal joins multiply counts. We use the established exclusions 7 and 21, and the realized values 3, 49, 63, 343. A direct replay checks the literal witnesses for 49 and 343 below, without requiring any large range census:

| H | outgoing masks |
|---:|---|
| 49 | (36,109,8,81,103,76,5) |
| 63 | (238,220,216,16,225,206,8,72) |
| 343 | (250,136,99,148,70,90,138,52) |

Each uses bit j for arc i→j in row i. These are witnesses, not order minima claims for 49 and 343.

## Exact reduction to two atomic families

Chan–Pak is equivalent to both of these tests: every odd prime other than 7 occurs, and every 7p with p≥5 prime, p≠7, occurs. The necessity is immediate. For sufficiency, all odd integers coprime to 7 are products of the available primes. All powers 7^e for e≥2 follow from 49 and 343 (every e≥2 is 2a+3b). If the 7-valuation is one and the other part has some prime p≥5, take 7p as a factor and fill the remaining primes. If the other part is a power of 3, the only excluded exponents are 0 and 1, and 63 fills the exponent-2 seed followed by multiplying by 3.

The same argument is target-by-target: any extra gap has an extra gap divisor from one of the two atomic families. In particular the least extra gap, or a sole extra gap if one exists, lies in one of those families. This specializes the earlier [finite-hole atom viewpoint](../../../01-canon/theorems/THM-2500-finite-hole-relative-atoms-prime-rays-and-cover-bounds.md); it does not assert the existence of an extra gap.

## Conditional classification of finite versus infinite extra gaps

Assume all odd primes except 7 are attained. Then every extra hole is 7m with gcd(m,7)=1 and v₃(m)≤1. Let P consist of those primes p≥5, p≠7 for which 7p is absent. Every prime factor of m other than 3 belongs to P, because any attained 7p would fill 7m. If P is finite, there are at most O((log X)^|P|) holes up to X, regardless of whether those holes eventually stop.

Under these same assumptions and finite P, there are infinitely many extra holes **if and only if** at least one p in P has no hit anywhere on 7p^j, j≥1. If every such ray hits at e_p, an omitted m cannot have v_p(m)≥e_p, by ordinal multiplication, so only finitely many exponent vectors remain. Conversely a wholly empty ray already gives infinitely many holes. On each ray the omitted exponents form an initial segment, possibly infinite. The first hit is necessarily realized by a strong component of that exact count: any proper SCC factor containing 7 would have a smaller index on the same ray, impossible. This locates where genuinely strong new constructions are indispensable.

Closure alone cannot improve this classification. For any divisor-downset D of the integers coprime to 14, containing 1 and 3 but not 9, the positive odd integers with the numbers 7m for m∈D deleted are multiplicatively closed. If a product were a deleted 7m, exactly one factor would have 7-valuation one and its 7-free part would divide m, making that factor deleted as well. This formal countermodel shows that neither a singleton extra 7p hole nor an infinite ray of holes is ruled out by ordinal closure alone; it is not claimed to be a tournament spectrum.

## A three-step strong tripling cube at the genuine 63 seed

In the displayed order-eight seed T with H=63, three disjoint transitive pairs are uniform modules:
* {0,5} has sole outside in-neighbor 4;
* {6,7} has sole outside out-neighbor 3;
* for {1,2}, every arc from an outside out-neighbor back to an outside in-neighbor starts at vertex 4.

All these statements survive replacing any of the other listed pairs by a directed triangle. A Hamilton path can therefore visit each listed module at most twice: every separate revisit needs the named singleton contact, and that singleton is used only once. The ordered one-/two-nonempty-cover counts of a transitive pair are (1,2), and those of a directed triangle are (3,6). Thus replacing any chosen pair by a directed triangle triples H exactly, independently of the other choices. The replacement can be done by adding a single vertex, retaining every old arc. It is uniform to the exterior and closes the original pair into a directed triangle. The original quotient is strong, so all inflations remain strong.

This proves one literal nested strong chain

    63 → 189 → 567 → 1701.

The graph family forms a cube of three independent switches. This is a positive statement at the actual exceptional chain, not just a count divisible by 7. The replay tests all cut choices in this specific one-apex grammar. From the supplied 63 graph the only tripling extensions are precisely the three pair-closures. Every first extension has exactly the other two choices; every second extension exactly the remaining choice; every third extension has none. The 3! histories are just the switch orders. The partition-by-seam quadratic form used in the replay checks every nontrivial cut at each step, so this is an exact exhaustion of successive single-vertex triplings starting from this graph.

By the inherited order bound and classification, the exhaustion applies to **every** strong H=63 seed: the [strong order bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) H≥5^((n−1)/3) rules out n≥9, and the inherited [H=63 census](../../../01-canon/theorems/THM-344-n8-h63-complete-omega.md) gives no order≤7 hit and exactly two order-eight classes. Their score sequences are converse, and the seed displayed here is one of them. Converse transfers the entire one-apex extension tree. Thus no fourth consecutive *one-vertex tripling* from a strong 63 seed is possible. This does not exclude a strong H=5103 tournament or a multi-vertex/changing-cut route.

## Arbitrarily long strong tripling chains with a changing 7-free core

In contrast to exhaustion at 63, there is no universal bound on the length of strong one-vertex ×3 chains at 7-divisible counts. In the [first selective chart](../session321_selective_returns_and_fano_wiring_20261113/SELECTIVE_RAYS.md) use a singleton U, so A=15 and h=1. Let the side reservoir V be an ordinal join of r marked transitive pairs. Closing j of them into directed triangles gives H(V)=3^j and Q_V(2)=4^r·3^j by ordinal multiplication. Hence the exact strong counts are

    3^j(15·4^r−1),   0≤j≤r.

Each closure adds just one vertex. When r≡3 modulo 21 the base has exact 7-valuation one: 4^21≡1 modulo 49 and 15·4^3−1=959=7·137 is not divisible by 49. The base is always a 3-adic unit, so the 3-valuation rises by exactly one each time. Thus there are arbitrarily long strong chains on 7-divisible targets with controlled valuations, but the cofactor changes with the desired depth. This construction does not prove an infinite chain for any fixed 7-free cofactor.

## A finite-cage criterion when more primes are missing

More generally let P be the finite set of missing odd primes. If some p in P has no attained positive power, its powers already give infinitely many holes. Otherwise choose an attained p^(a_p) for each p in P and form the finite cage M of products of p^(r_p) with 0≤r_p<a_p. This need not be the cage of all missing divisors from THM-2500: it is a residual coordinate box after removing attained powers.

The following two conditions are together equivalent to finitely many holes.

1. For each m in M, mq is attained for all but finitely many primes q outside P. Let R be the finite union of exceptional outside primes across all m.
2. For every m in M and every base b in {p^(a_p): p in P} ∪ R, some m b^k with k≥0 is attained.

Necessity follows from cofiniteness by considering each unbounded axis and prime tail. For sufficiency, factor any candidate n as m times the stripped attained powers and all its prime factors outside P. If some such prime q lies outside R, the product mq is attained by condition 1, and all remaining factors are attained. Thus holes are supported on P∪R. Within each residual m, their coordinates are exponents of the finitely many attained bases b. The attained exponent vectors are coordinatewise upward closed; condition 2 bounds the exponent of each coordinate in any hole, proving finiteness. This criterion exposes exactly which prime tails and fixed-base axes have to be supplied by additional structure; it does not assert those physical hypotheses.

## Isolation forces strong witnesses at neighboring targets

Suppose the only extra gap is a prime p. Then 3p, p² and p³ are attained and every one of those counts needs a strong core: each nontrivial factorization uses p (the p³ splits are p·p² and p·p·p). More generally pq needs a strong core for every prime q distinct from 7 and p, whether or not we assume p is the sole extra gap, provided pq itself is attained. Hence any finite-gap scenario with a missing prime implies a whole infinite family of forced strong realizations, not an infinite family of forced omissions.

If the only extra gap is 7p with p≥5, p≠7, the compulsory strong counts are 7p² and 21p. Their only proper two-factor splits use 7, 21 or 7p, so no product of smaller SCCs can supply them. Conversely, assume all primes except 7 and all 7q other than 21,49,7p are covered, and that 7p² and 21p are covered. Then every odd count other than 7,21,7p is covered: when v7=1 and a new prime q divides the cofactor, use 7q; for a cofactor supported on {3,p}, use 63, 21p or 7p² when its exponent vector is big enough. Higher 7-valuations use 49 and 343. This is an exact constructive two-target test for isolation conditional on the other frontier.


The independent [ordered Busch seam](BUSCH_INVERSE.md#7-an-unbounded-strong-ray-starting-exactly-at-63) does give an infinite two-vertex return starting at the literal strong 63 seed; it is not a fourth one-vertex step.
