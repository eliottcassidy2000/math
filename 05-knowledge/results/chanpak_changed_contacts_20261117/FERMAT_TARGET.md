# The Fermat target lands when the right core changes

The [transitive negative reader bank](OPPOSITE_SIGNS.md) cannot hit the prime quotient 65537. Here are **two distinct strong tournaments** that do hit the count 7·65537=458759, of orders 17 and 18, using the positive scaffold and nontransitive right readers. The count 458759 was already inside the established finite range; what is new here is the exact profile-based construction and an obstruction explaining why simply growing the earlier reader cannot generate the pure Fermat ray.

## Keep the scaffold, change the right core

Continue to use left masks (8,1,3,6), replace vertex 0 by B, and mark vertex 1 as terminal. Write P=P₂(B), h=H(B). For a pointed tournament, define its marked triple as full path count, ordered two-cover count with the first path ending at the port and second optional, and full path count starting at the port. The full left triple is

    (a,b,c) = (3P−h, 4P+h, P−h).

The first two readings were established in [the positive seam](OPPOSITE_SIGNS.md). The third has the same at-most-two-module-runs bound by the single exterior exit. On TT₁ and TT₂ its respective values are 1 and 3, giving `h+(P−2h)=P−h`. The literal replay checks all three readings on the nontransitive modules used below.

Take a pointed right reader of profile (d,e,f), but place its **converse** on the actual right of the gluing. All right-to-left arcs are present except for the single seam from the left port to the right port. The count is

    ad+be−cf = (3d+4e−f)P + (−d+e+f)h.

For paths avoiding the seam, each full right path and full left path concatenate except exactly the tuples where the right ends at its port and the left starts at its port; this gives ad−cf. A path using the seam gives a left cover and a right cover, the latter read by reversing its visits, giving be. Both the left tournament and each right reader here have a Hamiltonian path ending at the port. The reversed right terminal path starts at its port and completes a spanning directed cycle via the seam, so both outputs are strong.

| Right reader masks | Port | (d,e,f) | Admissible block (P,h) | quotient |
|---|---:|---|---|---|
| (126,84,40,18,36,74,28) | 5 | (43,106,0) | (824,49) | 79·824+9·49=65537 |
| (254,160,18,166,170,4,190,36) | 7 | (19,208,0) | (512,19) | 127·512+27·19=65537 |

The first block has rows (2,20,9,19,37,79,31). The second is (40,1,3,6,15,30) followed by one ordinal sink. The independent subset-DP replay constructs and checks both full tournaments, including strength. All masks are zero-indexed outgoing arcs.

These are physical all-order profile readers: the first fixed core returns 7(79P+9h) for arbitrary nonempty B, and the second returns 7(127P+27h). More generally the two coefficients are each divisible by 7 exactly when f≡0 and e≡d (mod 7). Writing f=7k and e=d+7(j−k), the quotient is (d+4j−5k)P+jh, with 0≤k≤j because e≥d.

## The second witness has a literal carry

The second reader hits `128P+1` whenever the block satisfies `27h=P+1`. Here 27·19=513=512+1 and 128·512+1=65537. To use this mechanism on a larger pure Fermat form, one could seek a block with P=2^s and h=(2^s+1)/27, so s≡9 (mod 18). That arithmetic condition does not produce a block. In particular ordinal-appending any nonempty tournament C cannot preserve the carry: if P=27h−1 and the two profile coordinates multiply, preservation would require

    27h(P₂(C)−H(C)) = P₂(C)−1,

impossible since P₂(C)≥2H(C). Thus this particular carry cannot renew by ordinal padding; it needs an actual change of block structure or reader.

## Why no fixed right reader gives the full pure Fermat ray

Suppose the block is TT_m and a fixed arbitrary right reader produces `7(2^(m+s)+1)` at two distinct m for a fixed integer shift s. Since h=1 and P=2^m, linear comparison forces −d+e+f=7 and slope 3d+4e−f=7·2^s. The right profile always has e≥d: cut a full path immediately after the port, obtaining a distinct ordered terminal cover. Hence f≤7. The slope is an integer, so s≥0; then substituting e=d+7−f shows the slope is 7d+28−5f and thus f≡0 (mod 7). There are only two cases:

- If f=0, the normalized slope is d+4, an odd integer ≥5, not a power of two.
- If f=7, then e=d. Cutting is now a bijection, so **every** cover concatenates. In particular the cover consisting of the singleton port and any Hamiltonian path of the deletion concatenates. Therefore f is exactly the Hamiltonian path count of the deletion, impossible at 7.

The same argument excludes any fixed right reader with a fixed ordinal core padded by TT_m: the constant coefficient must be −7 or +7 if one matches the corresponding signed ray; here it is nonnegative (`e+f−d≥0`). For the +7 case the core's odd count must divide 7, hence be 1 because 7 is a hole. This reduces to the transitive case. The obstruction uses the first forbidden count, not a finite failed search.

Replay with `python replay_fermat_target.py`, also verified under `-O`. Universal selection of a block for a prescribed prime is still open.
