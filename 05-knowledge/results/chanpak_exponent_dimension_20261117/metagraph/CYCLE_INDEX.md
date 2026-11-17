# Exponent coordinates for the self-converse vertex count

Let k=floor(n/2), E=binom(n,2), r the standard reversal, and U the
unframed orientations anti-preserved by r. The preceding proof
identifies self-converse classes with centralizer orbits on U. This
rewrites the classical self-converse count (Harary–Palmer,
[Enumeration of self-converse digraphs](https://doi.org/10.1112/S0025579300003910);
see also [Self-Converse Tournaments](https://www.cambridge.org/core/journals/canadian-mathematical-bulletin/article/selfconverse-tournaments/ED7BC2A1F6FF2D7F31547D771CE69B80))
in coordinates suited for comparing with the two-reflection law. No
novelty claim is made for the classical count or leading equivalent.

An odd-order centralizer element g permutes the k r-pairs with odd
cycle lengths. A cycle of l pairs has 2^(l-1) lifts: the product of
the swap signs around the cycle must be positive. Otherwise g would
have even order. Thus g has two l-cycles on vertices for each l-cycle
of pairs, and one extra fixed vertex at odd n. Conversely these
parameters give every odd-order g in the centralizer.

Use e(g) as in PROOF.md (g-edge orbits) and s(g) for the number of
r-pair cycles of g. On each g-edge orbit fixed by r, r fixes every
unordered edge (commuting with odd-order rotation). These are exactly the
orbits of matched edges; r swaps their endpoints. So every fixed orbit
permits an orientation; the other orbits occur in pairs. Hence
|Fix_U(g)| = 2^((e(g)+s(g))/2). For even-order g this fixed set is
empty. Burnside now gives, with m_l multiplicities of odd pair-cycle
lengths and s=sum m_l,

    SC(n) = 2^(-k) sum_{sum l m_l=k}
              2^((e+s)/2) prod_l ( (2^(l-1)/l)^m_l / m_l! ).

In the merged metagraph the node count is (A000568(n)+SC(n))/2. Here
a contribution is 2 raised to a dimension only *after* checking
fixed-edge signs, not because a cycle length or support size itself
becomes the dimension.

## First three corrections are explicit

Divide SC(n) by B=2^((E+k)/2)/(2^k k!). The identity contribution is 1.
For support t, the relative term of each g is 2^-d with
d=(E-e+k-s)/2. The first three classes of moved pair partitions are:

| moved pair cycles | vertex support | number of centralizer elements | d |
|---|---:|---:|---:|
| 3 | 6 | 8 binom(k,3) | 2n-6 |
| 5 | 10 | 384 binom(k,5) | 4n-20 |
| 3,3 | 12 | 640 binom(k,6) | 4n-24 |

To calculate d, all touched edges in these rows have orbit size p
(p=3 or p=5), with no fixed touched edges. Thus
E-e=(1-1/p)(t(n-t)+binom(t,2)); add k-s (respectively 2,4,4) and
halve. The lift counts are respectively 2!*2^2, 4!*2^4, and
(1/2)binom(6,3)(2!*2^2)^2 times the choice of support pairs.
After these terms any support is at least 14. The touched-support
bound and summation in PROOF.md, now with t0=14, give relative
remainder O(2^(-14n/3+O(log n))), nonnegative. The formulas concern
the self-converse sector, not the total A000568 error or the
Hamiltonian count spectrum.

This exposes a distinction lost in the leading equivalent: at odd n a
second reflection can appear with loss n-2, but correction to the
unlabelled SC count only begins with symmetries commuting with r, at
loss 2n-6, the same leading loss as at even n.
