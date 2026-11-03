# Torus knots separate three kinds of tournament memory

**Status:** Proved interfaces, with exhaustive finite controls. The standard knot-group, torus-cover and CRT facts are not claimed as new mathematics. The contribution is a precise interface: actual flip tori can retain winding while remaining completely H-blind; a different actual flip torus measures arithmetic synchronization through H; and a torus-knot complement cannot carry the full nonabelian signed-unit transport group. These are three different spaces and must not be silently identified.

## 1. The three spaces

1. A polygon-product torus whose vertices are tournaments and whose edges are single arc flips. Mixed squares specify commuting path homotopies. Its fundamental group is free abelian.
2. A simple primitive-slope curve on a **marked embedded torus in S³**, hence a torus knot. Knot type uses the ambient embedding and meridian/longitude marking, not merely an abstract cyclic path.
3. The complement in S³ of that knot. Its fundamental group is generally nonabelian and differs from both the torus group and the knot-circle group.

The polygon interiors are formal cells; no tournament is assigned to their interior points. A winding curve in the continuous torus is not automatically a literal flip path. An edge-path discretization is required, and its image may revisit vertices, so it is not automatically a knot embedding.

## 2. Exact slope lifts in the inherited physical cover

S170's literal six-vertex seed hexagon gives the physical cover

    R²/(6Z)² -> R²/(3Z)²,

with four pinned physical sheets, deck group (Z/2)², and H identically 29²=841. The six masks are 153,8345,8337,8369,177,185; their ordinal sums give its vertices. This realization and its full constant path-cover profile are inherited, not rediscovered here.

For a general finite-index lattice Λ≤Z², consider the torus cover R²/Λ -> R²/Z². Let w=(p,q) be primitive, and let K_w be a translate of its simple slope circle in the base. Put

    d_w = min{d>0 : d w ∈ Λ}.

**Theorem 1 (slope restriction).** The inverse image of K_w has [Z²:Λ]/d_w circle components. Each maps to K_w with degree d_w and represents the primitive lattice vector d_w w ∈ Λ. Its components are indexed by

    Z² / (Λ + Z w).

**Proof.** Monodromy of the base circle acts on Z²/Λ by translation by w. Every orbit has length d_w, so the orbit count and component indexing follow. A component closes at displacement d_w w. If d_w w=k v for v∈Λ and k>1, primitivity of w forces v=e w for an integer e=d_w/k<d_w, contradicting minimality. Thus the lifted vector is primitive in Λ. ∎

For Λ=mZ², primitivity of w implies d_w=m: if m divides both dp and dq, Bezout gives m|d. Hence there are exactly m lifted circles, each of degree m and, in the scaled marking of the source torus, each of slope (p,q).

For the S170 four-sheet cover, every primitive slope has exactly **two components and two-lap return**. The specific return is w mod 2, so the three nonzero parity classes identify three different pairings of the four physical sheets. Component count alone forgets even this pairing.

For an N-torus isotropic m-cover, the identical proof gives m^(N−1) components of a primitive slope circle, all of degree m. In S170's m=2 physical model, this is 2^(N−1) components and a nontrivial E_n deck translation.

## 3. What finite sheet memory provably cannot see

Fix the standard unknotted torus in S³ and its usual marking. For each k≥1, use the primitive positive slope

    w_k=(1+2mk,2).

All these slopes have the same mod-m deck return (1,2), the same degree-m lift, and the same m lifted-component count. They are nevertheless distinct torus knots: their classical genera are

    g(T(1+2mk,2))=mk.

Thus a fixed finite deck signature has fibers of unbounded knot genus. In S170's actual tournament torus with m=2, all these slope descriptions also have the same constant H=841 at every visited vertex. The same is true of its entire common path-cover profile.

The standard genus formula is an external classical input, not something established by the verifier. The blindness statement follows exactly from the displayed family and Theorem 1. It does not say that an unmarked flip word already defines a knot. The marking and ambient embedding must first be supplied.

A sharper organizational distinction: winding w is complete for based homotopy classes of loops in the torus; residue w mod m is complete only for deck return; neither is a substitute for an embedding of the curve into S³. The torus-knot genus comparison uses simple representatives of primitive slopes, not arbitrary self-intersecting paths with those windings.

## 4. A torus-knot complement cannot generate nonabelian nilpotent transport

Let Γ_(p,q)=<a,b | a^p=b^q>, with p,q positive and coprime, be the classical torus-knot group. Let G be a finite nilpotent group.

**Theorem 2 (exact representation collapse).** Every homomorphism Γ_(p,q)->G has cyclic image. More precisely, the map

    G -> Hom(Γ_(p,q),G),     g -> (a↦g^q, b↦g^p)

is a bijection. Consequently there are exactly |G| based representations, and simultaneous-conjugacy classes of representations correspond exactly to conjugacy classes of elements of G.

**Proof for a finite ℓ-group.** Let A,B be the generator images and z=A^p=B^q. Since z is a power of each generator, it is central in <A,B>. At least one of p,q is prime to ℓ; suppose p is. The order of A is an ℓ-power, so choose u with pu≡1 modulo ord(A). Then A=z^u, and A commutes with B. The case q prime to ℓ is symmetric. Choose integers r,s with rq+sp=1 and put g=A^r B^s. Using commutativity and A^p=B^q gives

    g^q=A^(rq) B^(sq)=A^(1−sp) B^(sq)=A,
    g^p=A^(rp) B^(sp)=B^(rq+sp)=B.

Every pair is therefore uniquely of the claimed form: uniqueness follows from g=(g^q)^r(g^p)^s. For finite nilpotent G, project to its direct-product Sylow factors. The ℓ-group argument in each factor yields the same claim in their product. Conjugation equivariance is immediate. ∎

**Corollary 3 (signed-unit obstruction).** For every inherited finite 2-group G_n generated by signed Cayley–Dickson left transports, any torus-knot-complement representation into G_n has cyclic image. Therefore it cannot realize full nonabelian G_n, a nontrivial image of a knot-group commutator, or the original noncommuting generating transport system used in S167–S170. Any holonomy obtained this way factors through the knot group's abelianization Z. An isolated associator-defect permutation R is not forbidden: like every g∈G_n, it can be the image of an abelianization generator and generate cyclic holonomy. The obstruction concerns the simultaneous noncommuting system, not the eligibility of individual permutations.

This is stronger than saying that a torus itself permits only commuting holonomy: the complement is nonabelian, but its finite nilpotent targets still collapse. The theorem is a restriction on the target, not a claim that knot complements have abelian fundamental groups.

**Hostile boundary.** The trefoil Γ_(2,3) maps onto S₃ by sending a to a transposition and b to a 3-cycle. Both sides of a²=b³ are identity. There are twelve such representations into S₃, six noncommuting, rather than |S₃|=6. Thus dropping nilpotence breaks the result.

**Constructive alternative.** The complement of an r-component unlink has free fundamental group F_r. Any r-generated G_n is therefore a quotient of that complement group. This supplies a topological source for nonabelian finite transport, but it is a link-complement construction, not yet a tournament gadget or a torus-knot construction. The Hopf-link complement, with group Z², would still force commuting images.

Do not use <a,b|a^p=b^q> with noncoprime p,q as a torus-link-group shortcut: its abelianization has torsion, whereas a link complement has free abelian first homology. Actual link-group presentations require the appropriate additional structure.

## 5. A new H-sensitive literal tournament torus

The inherited six-vertex gadget intentionally hides all structural-count variation. To make winding synchronization visible to H, use a different three-vertex block.

On vertices 0,1,2 fix 0->1. Let x indicate 0->2 and y indicate 1->2. Order the four choices as

    (x,y)=(0,1),(1,1),(1,0),(0,0),

and call the corresponding tournaments C_0,C_1,C_2,C_3. Consecutive indices modulo four differ by exactly one arc. C_0 is cyclic; the other three are transitive. Their Hamilton-path counts are

    h=(3,1,1,1),     h(j)=1+2[j=0].

For b∈(Z/4)^k take the ordinal sum

    U_b=C_(b1) ⊕ ... ⊕ C_(bk).

**Theorem 4 (counting torus).** The displayed labelled tournaments induce the flip graph C₄^□k=Q_(2k). With the polygon-product cell structure (do not fill the individual polygon circles), this is the one-skeleton of a k-torus. Every edge is a literal single-arc reversal, and

    H(U_b)=∏_i h(b_i)=3^(#{i:b_i=0}).

**Proof.** Two displayed states differ only in their two free bits per block, so their Hamming-distance-one adjacency is exactly Q_(2k). Every Hamilton path in an ordinal sum must traverse all vertices of its first block, then its second, and so on: it cannot return to an earlier block. Restriction and concatenation are inverse bijections between global Hamilton paths and tuples of block Hamilton paths. Multiply the block counts. ∎

This torus is **labelled**, and its class quotient loses much of its geometry: C_1,C_2,C_3 are isomorphic. It is not a new nontrivial cover of a tournament-isomorphism graph. Q_(2k) also admits a filled cube complex that is contractible; that is a different choice of two-cells. The torus construction retains coordinate polygon circles and fills only mixed-coordinate squares.

## 6. Hamilton counts detect CRT synchronization exactly

Sample the above torus by the finite affine orbit

    b(t)=b+t a mod 4,     t=0,...,L−1,
    L=lcm_i (4/gcd(a_i,4)).

Here a_i,b_i are residues represented by integers 0,1,2,3. A sample step may change several arcs; it is a **stroboscopic observation**, not a claim that the simultaneous step is a single flip. Any chosen sequential implementation gives literal flips, but adds intermediate samples not counted below.

For each coordinate set g_i=gcd(a_i,4), m_i=4/g_i. The zero condition a_i t+b_i≡0 mod4 is impossible if g_i does not divide b_i. Otherwise it is t≡r_i mod m_i, where

    r_i=−(b_i/g_i)(a_i/g_i)^(-1) mod m_i,

with the unique residue 0 when m_i=1. Call S compatible if every chosen equation is solvable and r_i≡r_j modulo gcd(m_i,m_j) for all i,j∈S. Put lcm(empty)=1.

**Theorem 5 (exact orbit census).**

    Σ_(t=0)^(L−1) H(U_(b+ta))
      = Σ_(S compatible) 2^|S| L/lcm_(i∈S)(m_i).

**Proof.** Expand ∏_i(1+2[b_i(t)=0]) as a sum over subsets S. The generalized Chinese remainder theorem says the simultaneous equations have either zero solutions, precisely when incompatible, or one residue class modulo lcm_(i∈S)m_i. Since that modulus divides L, there are exactly L/lcm_(i∈S)m_i solutions in the full orbit. ∎

For k=2 and a=(1,1), aligned peaks b=(0,0) yield counts 9,1,1,1 and orbit sum 12. Offset peaks b=(0,1) yield 3,1,1,3 and sum 8. Thus the same direction and period have different H-census because phase alignment differs. This supplies an actual arithmetic synchronization observable, rather than merely calling a formal torus a counting mechanism.

For k coordinates with a_i odd and every b_i=0, each m_i=4 and all peak times coincide, yielding orbit sum 3^k+3. At the opposite extreme, when the peak residues are distributed among the four times with occupancies n_0,...,n_3, the sum is Σ_j 3^(n_j). Therefore the exact counting problem is a convex occupancy problem: fixed total k is minimized when occupancies differ by at most one and maximized when all are equal in time. The minimal value for k=4u+r, 0≤r<4, is (4−r)3^u+r3^(u+1), and the maximal value is 3^k+3. This follows by moving one coordinate from an occupancy x to y with x≥y+2, which strictly decreases 3^x+3^y.

## 7. What this opens and what remains open

- **Established:** a marked slope's primitive lattice lift, exact component count, finite-memory blindness, and a target-sensitive knot-complement obstruction.
- **Established in real tournaments:** a 3k-vertex literal torus with H-dependent CRT synchronization; the 6N-vertex constant-profile torus remains the inherited contrasting example.
- **Open application:** combine constant-profile transport blocks with H-sensitive readout blocks through an explicitly authorized scheduler, then classify what a readout can learn about a transported sheet. A mere ordinal join gives independent product factors, not coupling, so no computation claim is made yet.
- **Open application:** build a physically intrinsic nonnilpotent tournament transport target admitting a nonabelian torus-knot representation. The S₃ trefoil example specifies the algebraic requirement, not a physical realization.
- **Open application:** use an actual link-complement cell presentation to organize noncommuting tournament circuits. Free-group realization alone does not identify gates with diagram crossings.

## 8. Verification and sources

Run `python3 check.py > certificate.json` and `python3 -O check.py > certificate.optimized.json`; compare the files. Checks use explicit exceptions, not removable assertions.

The verifier independently constructs the 128-element octonion transport group, checks every generator pair for six coprime knot presentations against the exact 128-representation parametrization, and includes the nonnilpotent S₃ hostile control. It computes Hamilton counts by subset dynamic programming on all displayed tournaments through nine vertices, verifies every displayed positive-coordinate step is exactly one arc reversal, checks every phase and direction through three counting blocks against the CRT formula, and directly enumerates finite-sheet slope orbits through modulus twelve. These are controls; the general statements rest on the proofs.

Sources:

- S170, `monodromy/NOTE.md`, §§3–5: exact physical C₆^N -> C₃^N torus, constant profiles and E_n deck charts. Published session index: https://github.com/eliottcassidy2000/math/blob/2d137c00eaf36a6fbe9f06044ac56d2c3bf7c7d6/00-navigation/S170.md
- David S. Bradley, *Torus Knots* (1971): primitive slopes and torus-knot group presentation, https://digitalcommons.usu.edu/gradreports/1130/
- Stanislav Jabuka and Cornelia A. Van Cott, *Comparing nonorientable three genus and nonorientable four genus of torus knots*: introduction recalls the classical orientable genus formula, https://arxiv.org/abs/1907.12970
- Allen Hatcher, *Algebraic Topology*, §1.3: standard monodromy classification of covers, https://pi.math.cornell.edu/~hatcher/AT/ATch1.pdf
