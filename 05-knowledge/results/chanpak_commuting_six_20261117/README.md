# A minimal +6 repair and a physical commuting prime square

The one-apex +12 limit is not a limit for unrestricted tournaments. Reversing one arc between vertices 4 and 7 in the two strong order-eight outgoing-mask tournaments

    B−=(22,132,136,147,166,79,31,97)
    B+=(22,132,136,147,38,79,31,113)

changes H from 383 to 389 and leaves P₂ exactly 5960. This is the smallest possible **nonzero** H change at fixed order with P₂ preserved, even for unrestricted edits: oddness forces it even and P₂≡(−1)^order H mod3 forces it to be divisible by 3.

The optional path-slot polynomial difference is 2q²(4−q²). In fact independent subset DP gives nonempty slot differences (6,−12,−72,−48,0,…), and the optional polynomial is their sum in the binomial basis. In particular higher slot reads need not be preserved.

## Growth and commuting vertical scaling

Put a transitive block T of order k+1 before B and a singleton r after B, with T→B→r. Let r beat precisely T's first vertex t. Excluding r→t yields the ordinal-chain count H(B). A path using r→t puts the remaining k transitive vertices and B into two labeled slots, giving 2^k P₂(B), via the ordinary single-reversal transitive-exterior construction. Any full path of B together with the ordered T and r is a directed Hamilton cycle, so these tournaments are strong. The same B arc always changes their exact absolute H by 6. At k=2 both reads, 24223 and 24229, are prime and checked by literal subset DP. At k=10 the prime reads are 6103423 and 6103429.

There is more than an arithmetic rectangle. Pad B on the right by three ordinal singletons and put it in the [five-arc sevenfold scaffold](../chanpak_changed_contacts_20261117/PLUS_RAY.md), which reads first 128P₂+H and then 7(128P₂+H). The first is 2^10 P₂(B)+H(B). The internal arc edit and the five exterior edits use disjoint arcs on the same labeled set, giving an exact commuting square of strong tournaments at counts

    6103423      6103429
    42723961     42724003

with a single horizontal edit at each level and the same five vertical edits in each column. This is a physical positive result, not a claim of infinitely many prime squares or of placing a square at a prescribed prime.

## Commutation does not preserve order at every intermediate exterior

For any subset of the five vertical edits, the reader still takes at most two runs through its variable block, hence is a linear read αP₂+βH. Its count difference along the horizontal edit is 6β. The original vertical order has β-values 1,5,−9,109,13,7, so the horizontal ordering reverses at an intermediate step even though both vertical tracks are strictly increasing.

The supplied exact skeleton replay computes all 32 linear reads by the two independent transitive substitutions. Across 120 orders of the same five vertical edits: 36 keep β>0 throughout; 29 are vertically increasing for all sources; **none** keep both properties even on this particular prime rectangle. This is a scoped route obstruction, not a ban on larger reroutings. It exposes why disjoint edits alone do not give a monotone count grid.

Replays: `python replay.py` checks the literal small witnesses, the full slot polynomial and the four prime/read values; `python ladder.py` checks the exterior ordering obstruction.
