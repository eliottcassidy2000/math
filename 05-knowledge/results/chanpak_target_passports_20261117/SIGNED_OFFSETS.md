# A varying-offset sevenfold atlas and the cost of using the wrong sign

Write (u,R)=(H(B),P₂(B)) for a nonempty physical block. The [negative reader](../chanpak_changed_contacts_20261117/MERSENNE_MODULES.md) is 7(64R−u); the [positive reader](../chanpak_changed_contacts_20261117/PLUS_RAY.md) is 7(128R+u). Both give strong tournaments. Feeding the negative reader B followed by one ordinal singleton gives the equal-center pair `7(128R±u)`. This section combines the two exact readers with older small dyadic two-cover profiles, rechecked in the replay.

## Choosing the source from the absolute target

These tiles have the following literal profiles:

| H | P₂ | outgoing masks |
|---:|---:|---|
| 5 | 2^5 | (2,4,9,3) |
| 19 | 2^8 | (40,1,3,6,15,30) |
| 223 | 2^12 | (2,40,11,49,71,213,143,31) |
| 599 | 2^13 | (96,21,25,51,225,198,142,15) |

Ordinal sums multiply both coordinates (each path slot concatenates independently). Thus for arbitrary nonnegative multiplicities a,b,c,d, not all zero, build

    u=5^a·19^b·223^c·599^d,   R=2^r,   r=5a+8b+12c+13d.

For every k≥r+7, pad by k−r−7 singletons and apply the two equal-center readers. This realizes **both** `7(2^k−u)` and `7(2^k+u)`, with physical strong outputs. Only the minus sign can use the marginal height k=r+6 with the unpadded block and the original coefficient 64. For u=1 use transitive B (the separate signed-ray proofs cover the smaller admissible heights).

This is a deterministic exact selector on its domain: for a proposed quotient p, try the dyadic powers immediately below and above p, factor the difference over the four listed primes, and check the corresponding height budget. No residue approximation is used. The restriction to neighboring powers suffices for this bank because the admissible offset is at most `2^k/128` even at the smaller minus height (splitting gives 2u≤R). The allowed offsets genuinely grow and are not a finite list of rays, but many p lie outside the domain.

Literal prime quotients include `2^12−5=4091`, `2^15−19=32749`, `2^28−599=268434857`, `2^24−95=16777121` and on the other sign `2^20+223=1048799`, `2^27+599=134218327`, `2^25+95=33554527`. The last offset 95 uses two different tiles. At a shared power-of-two center exactly one sign is divisible by 3, because u is 3-free. These observations are not prime-infinitude claims.

## A general noncancellation bound in both signs

Suppose a fixed reader `p=λR+εu` with (λ,ε)=(64,−1) or (128,+1) hits a prescribed quotient `p=2^k+δ` with fixed δ. For n=|B|, splitting each full path in n+1 positions gives R≥(n+1)u. In particular p≥127u in either case. At order n≥2, R≡0 (mod 4) by complement pairing and Rédei parity, so `u≡εp (mod 4λ)`.

Remove ordinal singleton components from B. The Hamiltonian count u is unchanged and only powers of two are removed from R. The remaining order is at most `3 log₃ u` by [Busch's strong bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3), using equality for triangles and the bound for components of order at least four. The [odd-pack expansion](https://arxiv.org/abs/2307.05569) gives P₂≤2^n u on an n-vertex block, since the product read of each bipartition is at most the whole count. It follows that `oddpart(R)≤u^(1+3 log₃ 2)`.

Now `λR=2^k+δ−εu`. For all sufficiently large k depending on δ, `|δ−εu|≤2^(k−1)`. If δ−εu is nonzero, the valuation of the right hand side is exactly the valuation of δ−εu. Consequently

    oddpart(R) ≥ 2^(k−1) / ((1+|δ|)u),
    u ≥ (2^(k−1)/(1+|δ|))^[1/(2+3 log₃ 2)].

Thus a fixed-offset reader with bounded source count can persist at unbounded heights **only through exact cancellation δ=εu**, and even then it needs a physical dyadic two-cover profile. Using the opposite sign for a nonzero δ forces polynomial growth in the source count. Also u≤n! gives n=Ω_δ(log p/loglog p) on this noncancelling branch, and the split inequality then gives u=O_δ(p loglog p/log p). These are necessary bounds, not a theorem of physical placement. Compare the [contact ledger](CONTACT_PASSPORT.md), where the same exponent comes from bounding a deletion's odd slope capacity rather than a valuation.

Replay: `python replay_signed_offsets.py`.
