# Where the exterior stops being a product: an exact colour-interface test

**Scope.** A graph-theoretic reduction and sharp three-step nonzero-gain example for the odd-degree Hadwiger–Nelson carrier programme. A five-cycle is not a six-colour obstruction. No improvement to the chromatic lower bound is claimed.

## 1. The right state is a translation, not a contact count

Retain the setting of THM-4661: a conjugation-stable subfield F of C avoids sqrt(-3), and M/F is finite of odd degree. An exterior point has at most one unit neighbour in F. Indeed, for a,b in F, simultaneous unit distances force (z-a)/(b-a) to satisfy a quadratic with coefficients in F; its degree divides the odd degree of M/F. Consequently a nonzero additive coset of F in M has *at most one unit vector*. If distinct u,v of norm one lay in the same nonzero coset, u would be a common unit neighbour of 0 and u-v in F.

Take a finite unit patch and collapse the points of each additive F-coset to a vertex. Do not include loops; retain all coset-pairs with an actual cross edge. For each oriented edge i→j, the unique unit vector u_ij in the corresponding nonzero difference coset is the displacement of *every* crossing edge. Each cross contact is therefore part of a single matching. Choose reference points r_i in their cosets, r_0=0 for the old coset F, and write v_ij=u_ij+r_i-r_j in F. Changing r_i changes the gains by a coboundary. The sum along a closed walk is intrinsic; call it its holonomy.

The following are equivalent, also modulo any additive subgroup J≤F:

1. The holonomy on every coset walk is zero modulo J.
2. The holonomy on every fundamental cycle of one spanning forest is zero modulo J.
3. Translating the references by elements of F makes every cross-edge gain zero modulo J.

For the nontrivial implication, take the sum along tree paths as the potential, then subtract it from the gains (or lift each potential in F/J to F and translate). In particular, exact zero holonomy is exactly the test for embedding these labelled patch edges into a Cartesian product of the base unit graph and the coset graph by translations in the fibres. With arbitrary finite patches this is a *subgraph* of that product, not necessarily the whole product.

**Colouring filter.** Suppose the old assignment extends to a q-colouring c of the entire base unit graph that is invariant under translation by J. If the coset graph is q-colourable and holonomy is zero modulo J, the patch preserves that assignment. After the above translation, colour x in coset i by c(x-r_i)+p_i mod q, with p a proper q-colouring of the coset graph, p_0=0 and r_0=0 in the old component. The two endpoints of every cross edge have equal c-values before the different shifts p_i,p_j. Edges within a coset obey the old colouring rule. This can be tested without any search through the interior colourings.

For the six-colour equality-tester search with nonadjacent old ports, if F has a five-colouring, first recolour one port with the spare sixth colour, giving a global separating six-colouring. Thus no *exactly* flat patch with a six-colourable quotient is a tester. More finely, a conjugation-fixed prime with a five-colourable residue unit graph gives a base colouring invariant under the additive maximal ideal J of its valuation ring. All base unit edges are integral and have nonzero residue; colour each additive coset of the valuation ring by the residue colouring after choosing a representative. If the ports differ modulo J, recolouring the whole J-coset of one port is allowed (it is an independent set), and preserves J-periodicity. Hence the six-colour filter also applies to holonomy tested modulo J instead of exactly. If the input is an arbitrary prescribed partial colouring without a global extension, or the ports coincide modulo J, this refined criterion does not supply one.

## 2. Earliest nontrivial return: quotient length three, two exterior vertices

Let f(t)=30t³−5t²−10t−1. It is primitive, cubic, irreducible mod 11, and has discriminant 70,700=100·707>0. Thus it has three distinct real roots t_i. Put F=Q(i,√707) and take M=F(t_1). Since [F:Q]=4, no cubic subfield of Q-algebraic numbers lies in F, so f also has no root over F and is irreducible over F. Its discriminant becomes a square in F, making M the cyclic cubic splitting extension. The quadratic subfields of F are Q(i), Q(√707), Q(√−707), so sqrt(-3) is absent from F and M.

For u_i=(1+i t_i)/(1−i t_i), each u_i is a unit and lies outside F by the inverse Möbius transformation. Substitution in the cubic X³−2X²+2DX−D, with D=(10+i)/(10−i), gives
(10−i)(1−it)³ times the substituted cubic equal to −2if(t). Hence the u_i sum to 2. The five successive vertices
0, u_1, u_1+u_2, 2, 1
are distinct (the middle two are exterior and differ by a unit), and all five consecutive edges are unit. No chord is possible, because any pentagon chord creates a triangle and M lacks sqrt(-3). The cross-coset edges follow a triangle of distinct cosets F, u_1+F, (u_1+u_2)+F, F, with holonomy exactly 2. The two middle cosets differ because u_2 is exterior. Three is minimal: a two-step return traverses the same translation matching backwards and cancels. Two exterior vertices are also minimal for an induced mixed pentagon: one exterior vertex would have two neighbours in F.

This realizes a geometric, nontrivial return but not a coupling core: THM-4661 already rules out a six-colour boundary deletion with this many exterior points. It supplies a short positive control missed by exactly flat probes, not a new candidate for the chromatic threshold.

## 3. A nonzero gain already in the old cyclotomic carrier

Take F=Q(ζ5,√−11), M=Q(ζ55). Write ζ=ζ11, and let R={1,3,4,5,9} be the quadratic residues mod 11. The ordered sum of the five exterior units ζ^r around the associated path returns to the original additive F-coset with displacement
a=Σ_R ζ^r=(-1+√−11)/2,
with the square-root choice fixed by ζ. This is the quadratic Gaussian period: the conjugate sum b comes from nonresidues, a+b=−1 and ab=3 (among the 25 sums r+n, each nonzero residue occurs twice and 0 occurs five times). Every proper nonempty subset sum is outside F. Indeed a subset whose sum is F-fixed must be unchanged by multiplication by R: compare the two coefficient vectors, whose difference has zero constant coefficient and whose associated degree-at-most-10 polynomial is divisible by Φ11. The transitive R-action on the five exponents then forces empty or full. This proves the five cosets on the quotient loop are distinct. The ordinary lifted path need not itself close in five edges; its endpoints differ by a. The gain has rational norm 3 in Q(√−11), so it is nonzero at the earlier conjugation-fixed residue prime over 19 as well. It makes the modular-flatness criterion fail, not the colours.

**Analogy, not identification.** In a tournament seam two scalar budgets can agree while the marked endpoint data make their allowed repeats different. In the unit-distance patch the discarded datum is the signed additive gain on a quotient return. In the square-sum construction an endpoint-preserving return is the datum that justifies iterating a block arrangement, not a certificate for one large size. Here the simple quotient picture suffices only while the translated contacts have zero gain, and only subject to the quotient colourability condition. The first active quotient triangle and the much later possible list-critical obstruction are different thresholds.

**Provenance.** THM-4661 / six_colour_odd_carrier_twenty_vertices_20261016 gives the odd-degree gate, types of list obstructions, and §5's product diagnosis for specific sampled patches. hn_elliptic_bridge_20261016c proves the earlier full boundary bound. The present argument retains the exact additive transport of every cross-matching, supplies its relative-mod-prime colouring test, and constructs two short nonzero-gain controls. These do not claim priority over standard voltage-graph switching.
