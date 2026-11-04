# Nine-name structural families and their singleton-sensitive anchor cycles

**Status: bounded finite-exact reconnaissance; independent review pending.** This packet is not a census of prime tournaments of order nine and does not advance the global highest-changed-node frontier beyond nine. Its displayed identities apply to arbitrary child orders by the previously reviewed all-order anchor theorem, not by a bound on physical child size.

## 1. Physical conventions and the exact finite family

Old quotient names are0 through8; row bit j in row i means i->j. The two actual ports are x->y, and a_i,c_i mean all physical vertices in child i point to x,y respectively. Source and target share the same physical child sets. Corresponding source and target child profiles match completely; different child positions may have unrelated profiles and orders. Every port query lifts to a fixed ambient by having both ports inherit the core root's exterior incidences.

For every subset D choosing one element from each pair {1,8},{2,7},{3,6},{4,5}, let C(D) have i->i+d modulo9 for d in D. This is exactly sixteen labeled circulant orientations for this cyclic action. In addition we retain the two explicitly constructed tournaments

    U9=(16,49,115,247,480,449,387,263,15),
    W9=(254,508,248,496,224,448,128,256,85).

U9 is obtained from the carousel with D={1,2,3,4} by reversing its first four consecutive names. W9 has transitive names0,...,7 and name8 dominates precisely the even names. No classification of all critical tournaments is assumed from these descriptions.

The mask naming a circulant chooses d rather than9-d when its bit d-1 is one. Multiplication of names by any unit u in{1,2,4,5,7,8} is an explicit isomorphism from C(D) to C(uD). Its inverse transports each entire physical child set and both masks. There is never an independent relabeling of the two responses.

## 2. The two prime circulant families and their structural anchor cycles

Exact testing of every nontrivial proper subset for the module property finds twelve prime labeled circulants. They form two unit-multiplier orbits:

    {0,3,5,10,12,15}, represented by C0=(480,449,387,263,15,30,60,120,240);
    {1,4,7,8,11,14}, represented by C1=(226,452,393,275,39,78,156,312,113).

The Hamilton counts of the representatives are3267 and3357, respectively, so these two families cannot be isomorphic. The other four choices, masks2,6,9,13, are composite. Their only proper nonsingleton modules have physical masks73,146,292, the three residue classes modulo3. This establishes the scope directly on the sixteen displayed orientations; no unverified complete tournament enumeration is used.

All sixteen circulants have endpoint parity p=0. Vertex translations make every start count equal, every end count equal, and their sums are the same Hamilton count; hence their integer start and end vectors coincide. Their second-endpoint residues also vanish, as checked by the complete repeated-run primitives. Thus the uniform nontrivial-child endpoint test alone sees none of them.

Put v_j=1 when child j is nonsingleton and v_j=0 for a singleton. For C0, the complete exact anchor primitives give, for every i modulo9,

    Delta(0,e_i)/2 = v_(i+4)+v_(i+5) over F2.             (1)

The contrast in(1) is understood modulo four. Consecutive pairs(i+4,i+5) form a connected9-cycle. Therefore all nine first anchors vanish if and only if every v_j has the same value.

For C1, the corresponding second anchors give

    Delta(e_i,0)/4 = v_(i-1)+v_(i+1) over F2.             (2)

These contrasts are known to be divisible by four and are understood modulo eight. Their pairs form the connected step-two cycle modulo9. Again all nine displayed anchors vanish if and only if every v_j has the same value.

The first-anchor matrix for C1 alone has rank six and does not suffice; its second-anchor matrix has rank eight. For C0 the first-anchor matrix has rank eight. The four composite circulants instead have rank-six anchors and eight blind masks, exactly the assignments constant on each of the three displayed modules. This distinction is retained rather than folding composite examples into the prime claim.

It follows that every prime circulant in this finite family is distinguished through arbitrary matching children whenever some, but not all, children are nonsingletons. The exact physical query may depend on that singleton-position pattern. The two uniform patterns V=0 and V=511 remain blind to this particular anchor mechanism. The complete named banks close V=0 separately: at actual query(a,c)=(1,2), the singleton contrasts for C0 and C1 are218 and-298, both two modulo four. The universal modulo-four higher-cover residues therefore lift these gaps through every fixed ambient. Consequently every singleton-containing assignment is covered; only V=511 is left by this packet. This remaining scope does not assert a full-response collision or exclude other separating methods.

## 3. Two other displayed prime-nine quotients

Direct module tests show that U9 and W9 are prime. For U9, endpoint parity is one at names4 and8 and zero elsewhere. Thus the actual query(a,c)=(16,0) has complete augmented-profile contrast(2,0,...) modulo four for every allowed child assignment. This follows from the general endpoint theorem and survives every fixed ambient.

For W9, endpoint parity is zero at every name. Its complete affine anchor identities include

    Delta(0,256)/2 = v_0+v_7 over F2,
    Delta(4,0)/4   = 1+v_0+v_7 over F2.                 (3)

These two bits sum to one. Hence the fixed two-query battery separates every singleton-position pattern and all actual child profiles compatible with it. The first line has complete profile(2b,0,...) modulo four; the second has complete profile(4b,0,...) modulo eight. Both therefore survive every fixed ambient.

## 4. Why the finite identities have unbounded-child and ambient scope

We use Section2 of the separately reviewed `prime8/PROOF.md` only as an all-order lemma; its quotient-order-eight classification is not being extrapolated. With exact endpoint counts S_i,E_i, start repeated-child skew R_ij and one-port repeated-child skew K_ij, it gives

    t_i(V)=S_i-E_i+2 sum_j v_j R_ij modulo4,
    d_i(V)=D_i+2 sum_j v_j K_ij modulo4,
    Delta(0,e_i)=d_i(V) modulo4,
    Delta(e_i,0)=2t_i(V)+2d_i(V) modulo8 when p_i=0.

This reduction is valid for arbitrary nonempty child tournaments: normalized q2 is2v_j modulo four, the q3 contributions cancel by return-segment interchange, and all higher contributions vanish. The common odd Hamilton normalization fixes the resulting residues. For a nonzero second anchor, the odd-total cycle correction vanishes because y has only x entering it and the cycle-parity contrast is p_i=0. Its complete cover profile is consequently pure in the Hamilton coordinate, not merely a scalar witness. Fixed ambient substitution multiplies that coordinate by an odd Hamilton count.

The retained quotient-word primitives determine every affine coefficient. Evaluating all512 vectors v verifies the stated two uniform blind patterns and the W9 contradiction; the affine identities, not those512 choices interpreted as bounded child examples, supply the unrestricted child-profile scope.

## 5. Source-relative original-module consequence and boundaries

For any covered source core and fixed ambient, take its displayed baseline battery and pair each query with J(a,c)=(511-c,511-a) on core bits only, keeping exterior incidences unchanged. For original-module XOR families supported inside the core and selecting the root oddly, every proper original module lies within a prime child. Target child profiles dominate their corresponding source profiles. Complete converse-complement duality supplies a baseline-negative member, and nonnegative run-substitution coefficients prevent these profile increases from canceling it.

Thus the paired source-relative batteries have at most two queries for U9, four for W9, and eighteen for a covered mixed-position prime circulant source. For the all-singleton circulant source, pair its displayed(1,2) query with J, giving at most two queries. The battery does not depend on the selected family, although a baseline-negative choice may depend on the fixed ambient. This is not pairwise injectivity among arbitrary family outputs, does not permit operations on exterior names or newly created modules, and does not prove the global original-module reconstruction conjecture.

## 6. Reproducible finite evidence

`families.json` retains all18 constructed tournaments, every exact9-by9 joint endpoint matrix, every proper nonsingleton module, complete anchor primitives, and every one of512 singleton-position evaluations per model. The source and target quotient orientations and actual physical names remain explicit. The C++ primitive uses exact integer multiset dynamic programming, including repeated occurrences as indistinguishable visits.

`named_banks.json` retains complete named response basis values and Boolean coefficients for both orientations of C0,C1,U9,W9. Every named cubic and quartic coefficient, including every zero, is retained. These are exact banks for singleton children on the nine named vertices; they are not full banks for every large substitution. Coarse-pattern conclusions elsewhere are never substituted for these complete named data.

The finite evidence is deliberately bounded to the stated family. Its strongest structural feature is that the affine residual equations are the edges of a connected cycle, leaving precisely the two uniform singleton-position patterns for the prime circulants. The separately retained singleton nonanchor probe then closes the all-singleton pattern. This does not state that all prime-nine quotients share that structure.

Run `python3 verify.py` and `python3 -O verify.py`; both use only Python standard-library tools and an installed C++17 compiler. An independent9! permutation implementation checks all joint endpoint coordinates. The complete named Boolean polynomials reconstruct2,097,152 response values; reproducible hashes of their little-endian32-bit arrays are retained without discarding any coefficient. Physical controls retain actual augmented rows, child sets, port masks and exact Hamilton counts, including cyclic fixed-ambient lifts, nontransitive children and unequal child orders. All counted physical graphs have at most20 vertices, so their Hamilton counts fit unsigned64-bit integers.
