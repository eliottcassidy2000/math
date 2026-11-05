# Independent review: rank-one rooted orbits and strong generating functions

**PASS. No correction or blocker.** The rooted/unmarked correspondence, exact stabilizers and centralizers, fixed-word mass, and ordinary/weighted strong-output generating functions are valid with the stated rank-one and anti-fixed-root scope.

## Frozen source

Reviewed in `session197-orbit-wildcard/rank-one-orbits/`:

- `PROOF.md`:4eee0617d36e819f3e3c5bed13e8a4cb20370e8cf36ffaf3f8771092b25fe245
- `verify.py`:ae2abdf96fc6f55a24dfc9a37c59fcf69672fae295c21a419a387427c95047ab
- `receipt.json`:376e0603699c8ace29595392331d8b3c24c5420230fb7e1e483ed2911f19d037

Copied source normal and Python-optimized outputs reproduce the frozen receipt exactly. The named binary-rigidity proof and the already reviewed anti-involution mass theorem were also read for their specific dependencies. Original source files remain unchanged.

## Bijection and group-action fibres

For a marked rank-one output, the out-neighborhood and the unique out-to-in crossing arc recover the rooted side(S,a) intrinsically. Thus the rooted inverse does not depend on an arbitrary choice of anti-involution. The unmarking step correctly uses conjugacy in the full automorphism group: transport the source anti-involution through an unmarked isomorphism, then conjugate it to the target anti-involution. Its unique fixed vertex must follow, recovering a rooted isomorphism. This proves that forgetting the output root introduces no additional identifications.

With G=Aut(S,a), an automorphism fixing the output root preserves its two neighborhoods and both unique crossing endpoints. The two side automorphisms are otherwise independent, giving G×G. Conjugation by the canonical anti-involution swaps the factors, so its centralizer is the diagonal G. Every full-group automorphism commuting with the anti-involution fixes its sole fixed vertex; therefore this is also the centralizer in the full group, not merely in the stabilizer. The count of root-fixing anti-involutions is|G|, as stated.

The fixed-word proof handles all possible half-words. A fixed first half contains exactly one member of each paired vertex set. It must end in the in-neighborhood. Any visit to the out-neighborhood would require the unique crossing arc to return, but that would use both members of one pair in the first half. Hence the half is entirely the in-neighborhood, giving exactly H(S) fixed words. The argument includes t=1 and excludes t=0 correctly.

G acts freely on Hamilton words and on the complete labelled two-path endpoint covers used for q. Thus both H(S)/|G| and q/|G| are integers, and the former is odd. The root-stabilizer orbit formula follows by dividing the rank-one square identity by|G|². The source correctly avoids substituting the stabilizer for the full unmarked automorphism group; C3 is a decisive boundary example.

## Mass, strongness and series

Normalizing a side Hamilton word forces exactly t-1 arcs and leaves(t-1)(t-2)/2 free bits, with t choices of marked position. Dividing by the free rooted automorphism action gives the exact rank-one c-mass. Combining this with the established full self-converse normalization count gives the stated sampling fraction. It is a fraction of normalized objects, not a uniform unlabelled-class fraction.

Strongness is equivalent to every side vertex reaching a, because the distinguished crossing is the only exit from the out-neighborhood. Converse symmetry supplies the corresponding access to every in-neighbor. In a tournament this is exactly membership of a in the terminal SCC.

Ordered SCC decomposition then gives R=A²R_str and C=A R_str=R/A for ordinary isomorphism classes. Rooted class counting is by vertex orbits, so R is not xA'. The weighted quantity H/|Aut| multiplies over SCCs, and orbit-stabilizers instead give the derivative xM' for the weighted rooted series. Thus C_mass=xM'/M is correct. Empty prefixes/suffixes are accounted for by A(0)=M(0)=1; no empty rank-one output is introduced.

The inherited 7p restriction is used at its correct scope. The binary square formula h²+2^m q² implies that7 divides both h andq whenever it divides the output, so49 divides the output; positivity then excludes the exceptional potential49. No unrestricted7p obstruction is asserted.

## Independent controls beyond the source range

`check_orbits.py` independently regenerates the complete unlabelled side bank through order5 by one-vertex extension and degree-class canonicalization. It builds outputs with the root last, rather than first, and checks67 distinct rooted side classes.

The independent checker enumerates full output automorphism groups and full anti-involution sets by adjacency-constrained backtracking, not just port-preserving permutations. It verifies:

- the root stabilizer G×G and the full centralizer G;
- the full anti-involution conjugacy-orbit size|Aut(T)|/|G|;
- equality of the root orbit with the fixed vertices of all anti-involutions;
- complete unmarked-output canonicalization injectivity through side order5;
- the q and H orbit divisibilities and actual full-graph Hamilton counts;
-189138 symmetric half-word candidates;
- terminal-SCC strongness and both generating-function convolutions.

The inventories extend the source by one order:

- rooted side classes:1,2,4,12,48;
- strong rank-one classes:1,1,2,7,33;
- full rank-one c-mass:1,2,6,32,320;
- strong rank-one c-mass:1,1,4,25,281.

There were213 full automorphisms and117 full anti-involutions across the67 outputs. The C3 full-group/stabilizer distinction is explicitly reproduced. Normal and Python-optimized independent outputs agree.

Independent checker SHA256:25e3c1134f5722d7a7d4e8bf5802cc0d2273c9e2e03994b5f64f00f3a1f2678b.
Independent receipt SHA256:863fb8c6ad25b97b6d14e5950e29d5e4a96b47aada459f171eedf54fa8dc345a.

The optional parity consequence suggested by the author is separately written in `PARITY.md`; it is not an alteration to the reviewed source packet.
