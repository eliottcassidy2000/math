# Opposite axes, the Catalan unit gap, and choosing the terminal core

The proposed statement is

    (every odd prime p≠7 is a count) OR (every 7q, q odd prime ≠3, is a count).

It is not enough to prove that for every prime p, at least one of p and 7p is a count. The requested disjunction forbids a pair of missing targets p and 7q for independently chosen primes p,q. Nothing here proves the full disjunction. It gives structural reasons why two substantial approaches do not suffice and an exact obstruction at a prime where both unrestricted witnesses exist.

- [All-depth branching barrier on both axes](BRANCHING.md) and [rational weight checks](check_weights.py): arbitrary ternary cyclic nesting, arbitrary transitive leaves, and a bounded number of arbitrary one-apex ordinal sources still give only O(X^0.99) candidate counts on the stated axes.
- [Complementary unary wrapping barrier](UNARY_WRAPS.md): the narrower, changing-width nested reader with one arbitrary ordinal one-apex base has a better exponent.
- [Exact terminal-core obstruction at 257](TARGET257.md) and [replay](target257.cpp): the normalized reader from the previous session misses a prescribed prime and its sevenfold companion, even though both exist.

- [Apex and one spine tag classified at every order](APEX_PAIR.md): the only internal seed family is the known 11·2^k+5 ray.
- [Physical Bézout repair](BEZOUT.md): changes of +776 and −778 preserve the surrounding interface and give an exact −2, with literal twin-prime outputs.

- [A third-slot reader](THIRD_SLOT.md) and [literal replay](replay_third_slot.py): two marked spine vertices give a physical prime/sevenfold family that genuinely reads three cover coordinates.

## Why cross-axis obstruction is a strong-component question

Ordinal composition multiplies counts. A prime count p, if realized, is carried by one strong component; the remaining components are singletons. The same holds for 7q with q≠3,7, since the component product cannot split off a factor 7. For q=7 the same argument applies to 49.

If p and 7q are missing, a tournament counted by 7pq with p,q distinct and outside {3,7} can only have a single nontrivial component of count 7pq, or two of counts 7p and q in either order. The other factorizations split off a 7, a p or a 7q. When p=q and both p,7p are missing, a 7p² witness would have to carry the full product in one strong component. Thus seeing the large composite in the spectrum does not necessarily choose the intended branch. Any contradiction based on factoring 7pq must control which strong factors occur, not merely prove existence of a composite witness. In fact the arithmetic monoid obtained by deleting {7,21,p,7q} from the odd integers is multiplicatively closed for distinct primes p,q outside {3,7}. This shows why compositional closure alone cannot prove even the requested either-or.

## The 8/9 relation, with the coordinate that must be retained

The arithmetic picture is 5+3=1+7=8=2³ and 7+2=9=3². The tournament recontact menu from the previous session realizes the odd multiples 1,3,5,7,9 of the same source count, not the even 8p or 2p as individual Hamilton counts. There is also a literal 9p=8p+p split in one graph: take the source followed by two ordinal directed triangles. Each pair of triangle starts contributes p, yielding a 3-by-3 array of p-path classes. Mark one pair and the remaining eight cells can be labelled by three bits. Without a mark the two independent triangle rotations act transitively on the nine cells, so there is no invariant choice of a distinguished cell. This is a concrete reason a unit correction needs a boundary marker, not just a scalar equality. The additive decompositions of 8p in the recontact menu instead concern different graphs; their parts are not canonically subsets of this cube.

The archived [Fermat–Catalan clock result](../../../01-canon/theorems/THM-4516-perfect-power-clocks-are-fermat-catalan-identities.md) tracks a raw gap 2^K−q^X against a separately chosen carry. Here 9p−8p=p; normalizing by p produces the unit gap but does not supply the marker in a physical path cover. Also 8n is a cube and 9n a square simultaneously iff n is a sixth power (divide the putative cube root by 2 and square root by 3); actual consecutiveness then forces n=1. This is the precise bounded analogy to Mihăilescu's theorem, not a reduction to it. See [Mihăilescu's proof announcement](https://www.math.mcgill.ca/goren/cnta7/invited/Mihailescu.pdf). In the target problem, the comparable retained coordinate is the defect U+E−H(F) in the terminal core. The 257 obstruction shows that solving the scalar read alone can leave this coordinate nonzero.
