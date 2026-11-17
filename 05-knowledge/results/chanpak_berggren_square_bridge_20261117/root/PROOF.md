# An exact Berggren address in strong tournament pairs

The useful transplant from primitive triples is not “three recipes.” It is **three disjoint image regions plus an inverse that strictly decreases height**. Here is a physical tournament family with that same address property, and a final scalar read that exposes it without needing to know the rooted tournament in advance.

## The arithmetic decoder

Use Euclid parameters m>n>0, coprime and of opposite parity, taking the odd leg m²−n² and even leg 2mn. In this convention the three [Berggren moves](../../session256_chan_pak_triples_20261111/spinor/PROOF.md) on the parameter pair are

    A: (m,n) → (2m−n,m),
    B: (m,n) → (2m+n,m),
    C: (m,n) → (m+2n,n).

They preserve the primitive parity lattice and are unimodular, with determinant ±1. The slope m/n after the move lies respectively in (1,2), (2,3), or (3,∞), so it identifies the last letter. Inverting inside its interval decreases m. Slope 2 gives only the root (2,1), and slope 3 would give the excluded non-opposite-parity pair (3,1). Thus every primitive triple has one exactly decodable word. The root gives (3,4,5). The methods of ordinary binary Stern–Brocot decoding and this ternary oriented parity decoding are related but use different chambers.

## A different, strongly connected tournament address

This is a rooted migration construction, **not** the earlier [one-apex substitution](../../berggren_target_return_20261111/berggren/PROOF.md), whose scalar H has collisions and whose chart can change under isomorphism. Use the [rooted migration law](../../chanpak_bank_fusion_20261117/migration/PROOF.md) with a width-5 transitive spine and a new root. The state (h,R,S) records H, the covered path-end weight at the marked root, and the full path-start count; it also retains one full path ending at the root. Start from the order-two transitive graph at its sink, so (h,R,S)=(1,2,0).

All our addresses have leading and trailing bit 1. The leading bit ensures the new root reaches the entire old graph through the spine. Every old vertex reaches the previous root and then the new one along the assured ending path and the return arc. Thus each child is strong. The trailing bit removes S from the h update, so the two observed coordinates alone update as actual matrices:

| letter | address | matrix on (h,R) | range of R'/h' for R/h≥1 |
|---|---:|---|---|
| A | 17 | ((9,18),(15,32)) | [47/27,16/9) |
| B | 19 | ((5,20),(7,32)) | [39/25,8/5) |
| C | 21 | ((9,22),(11,32)) | [43/31,16/11) |

Each letter adds six vertices, i.e. three pairs, preserving even order. The matrices are just [[F(u),u+1],[32−rev₅(u),32]]. Their positive determinants (18,20,46) make the projective maps increasing. The displayed bands are disjoint, lie above 1, and avoid the seed ratio 2. Decode by the unique ratio band, invert the corresponding matrix in exact integers, check the positive cone again, and continue. h decreases at each valid inverse. This recognizes the entire orbit, rejecting not just wrong ratios but failed integer divisions and wrong seeds. The same word that names a primitive triple therefore names a distinct physical rooted tournament profile, at order 2+6d for depth d. Neither the slopes nor the determinants are identified with those of actual triangles.

## Write both coordinates into one Hamilton count

At depth d, all ratios are ≤2 including the seed, so h≤53^d<2^(6d+1). Append once more by migration, now with width 6d+1 and all its bits 1. Its F value is 1, so the resulting literal strong tournament on **12d+4 vertices** has

    N = h + 2^(6d+1) R.

For the known d, Euclidean division of N recovers R as quotient and h as remainder, then the ratio decoder above recovers all letters. This gives **3^d distinct ordinary Hamilton counts at each such even order**, not just 3^d rooted states. It is a direct code from primitive triples of a given depth to actual strong tournament counts. At d=1 the three counts are 10157, 9133 and 9653; 9133 is prime by trial division. Before the final read, A and B both have h=45, so dropping the side coordinate would visibly erase the address.

The last all-one step writes the state into one count, rather than carrying the same projective state forward. This tradeoff is important. It does not make the count a minimal-order witness, or choose a word for an arbitrary target. In fact each inner step multiplies h by at least 25, and the final read adds a factor at least 2·64^d; this entire scalar family has at most O(X^(log(3)/log(1600))) outputs up to X. It is a certified addressable family, not a universal selector.

## Every odd dyadic phase nevertheless occurs

Even the sealed family meets every odd residue modulo every power of two. This is a local coverage theorem, not an absolute selector. Use the A matrix, whose characteristic polynomial is x²−41x+18. Modulo 2^a it has one odd root λ and one even root μ: modulo 2 the roots are simple, so lift one bit at a time. We have λ≡3 mod8, μ≡6 mod8, and the difference is invertible. By the order-two recurrence, after t A letters the first coordinate is

    c λ^t + e μ^t mod 2^a,
    c = (h(A v)−μ h(v))/(λ−μ),     v=(h,R).

The μ term vanishes for t≥a. The unit c satisfies c≡7h+2R mod8: the raw seed has c≡3, while the prefix AC gives (h,R)=(2143,3023), hence c≡7. For a≥3 the cyclic subgroup generated by λ consists exactly of the units congruent to 1 or 3 mod8. Indeed v₂(λ²−1)=3 and each repeated squaring raises it by one, so the subgroup has order 2^(a−2), the cardinality of those two residue classes. The two choices of prefix therefore give both unit cosets. Taking t arbitrarily large in the selected period gives every odd residue; the sealing term is zero mod2^a once 6d+1≥a. The cases a≤2 follow by reduction.

Thus this concrete ternary address tree is globally sparse but has no dyadic hole. Recovering a word from a valid state is genuinely different from proving that an arbitrary prescribed absolute target is such a state.

This is why the geometry of the 3,4,5 example matters: two square-sum reflections compose to a rotation on a precisely specified residue space, and Berggren inverts a direction on a precisely specified parity lattice. Here likewise a marked two-coordinate interface, not the scalar total, is the domain with a clean inverse. The main remaining tournament question is to obtain such a division of domains for the core of an arbitrary prescribed target.

Run `python replay.py`; it checks the primitive and tournament inverse on 3,280 words, exact mixed-radix distinctness at fixed depth, three literal strong order-eight cores and their three literal order-sixteen reads, and the prime certificate.
