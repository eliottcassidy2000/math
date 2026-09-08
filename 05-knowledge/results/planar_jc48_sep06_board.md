# Planar Jacobian conjecture: 48-hour research session

**Status: ACTIVE RESEARCH / JC(2) OPEN.** Started September 6, 2026,
20:40 UTC, from 2d3c53942d51. Continue through September 8, 2026,
20:40 UTC (14:40 America/Denver). The user redirected this task from the
previous broad portfolio to the complex planar Jacobian conjecture.
Worktree: `/tmp/math-wt-planar-jacobian-sep06`; root owns Git.

The target is: every polynomial map (P,Q): C^2 -> C^2 with nonzero
constant Jacobian has a polynomial inverse. Restricted chart closures,
finite jets, rational primitives and formal solutions do not settle it.

## Inheritance and current concept board

The closest mechanisms are THM-4438's selected row-15 lift, THM-4426's
restoration of omitted source coordinates, THM-4411/4412's collision and
seminormal suspension, THM-3770's principal-part equalizer, and THM-3412's
linear-variable response arms. The linked proof notes give exact slugs.
Hostiles are fixed-prefix rank claims applied to moving prefixes,
rational exactness without polynomial descent, scalar collision periods
after a dimension change, and abstract deleted-sheet counts without an
actual affine source. MISTAKE-237 forbids an NC2-to-JC implication;
MISTAKE-374 separates generic divergence vanishing from polynomial extension.
The least-used sidecars tested here were earlier response kernels,
quadratic tangent evaluation, and transverse derivatives of component jets.

| Lane | Preserved object and result | Next decisive question |
|---|---|---|
| Anchor: moving source | The finite five-dimensional response remains valid, but every nonconstant polynomial Hamiltonian in K+Delta K[p,y] has no nonzero rational scalar time | Change the invariant or carrier; different polynomial lifts inside the cusp ideal are now closed, and compositions still need their actual invariant/regularity sidecars |
| Niche: collision | All scalar two-form periods miss exactly the quadrics through tangent directions in a three-dimensional target | Carry the quadratic data through an actual earlier source deformation and preserve descent |
| Anchor: infinity | The whole one-odd-cusp, two-ordinary-cusp and three-ordinary-cusp (4,6) classes are excluded as irreducible whole Keller support; the last uses actual marked D4 and arbitrary retained subsets | The adjacent(5,3,3) family has an actual marked H4 quotient; involutions, single cycles and mixed(3)(2),(3)(3),(4)(2) meridians are excluded; moved support is at least7. The audited mapping-degree floor is16; the necessary scalar rows at16 still need actual cycle/incidence structure |
| Wildcard: response connection | The fixed smooth polynomial family F=arb+c+h(b), h in b³C[b], has unit class of exact primary order2; the canonical derivative generates its full component-labelled arm | Map an actual source-normal response to this paid operator without changing its polynomial ring |
| Wildcard: actual boundary surface | Shared-root jets are necessary; all single-point octic boundary divisors are excluded, including infinity, with unrestricted polynomial-mate degree; the constant-D finite6/infinity2 class is closed at every finite location | General leading exactness closes all(4,2,2) positions and classifies every degree<=8 leading differential. The5+1+1+1 class is polynomially closed with sharp rational hostiles. Exact pencils classify every global quadratic rational mate; spend those structures on the surviving quartic lower sections |

## Audited checkpoint

1. [Compensated weight-22 transport](planar_jc48_sep06_weight14.md)
   is **PROVED**, with an independent analytic audit and a separate
   Fraction implementation. The exact packet
   `p^5 y^4-(7/10)p y^6-(508/135)p^2 y^6`
   supplies the response formerly paid by weight24 after changing unknown
   rows12–15. A finite additive translation, with coefficient
   `235202/27945`, transports every point of THM-4438's boundary G_m and its
   A10 terminal fiber. It changes the prefix; the old least-weight theorem
   in the frozen-prefix packet remains correct. No weight22 minimality or
   later-row solution is proved. Producer169 and independent referee93
   exact gates pass.
2. [Collision quadrics](planar_jc48_sep06_collision.md) is **PROVED**:
   for spanning tangents in a three-dimensional target, the kernel of all
   scalar-relation wedge periods modulo common motion is exactly the space
   of quadrics vanishing on those directions. Arbitrarily many conic
   directions remain blind; six suitable directions give a sharp complete
   test. The actual seminormal suspension has the all-order graph criterion
   `H(-1,s)=H(0,s)=H(1,s)`. Producer257 gates and an independent analytic/
   source audit pass. This is an ambient collision theorem, not JC(2).
3. [Vertical torsion and its connection](planar_jc48_sep06_torsion.md)
   is **PROVED**, with classical relative-exactness overlap explicitly
   credited. For smooth P with rational Hamiltonian constants k(P), its
   response torsion has one full principal-part arm per fiber component
   modulo the diagonal. The intrinsic operator
   `nabla[a]=[V(a)+div(V)a]`, V(P)=1, differentiates those parts exactly.
   Nonzero primary order n becomes n+1. This makes the repo's special
   divergence ladders presentation-free. Generic cohomology and inverse
   effectivity are still separate. Producer278 gates and independent audit
   pass; no claim of literature priority is made.
4. [Nodal nonproperness](planar_jc48_sep06_infinity.md) is a **RECOVERED
   COROLLARY with an independently audited shorter proof**. The nodal
   ledger `d-1=sum delta_i-sum overlap_p` and the irreducible wholly-nodal
   exclusion were already consequences of public Paper II's Corollaries
   5.2 and4.2. Keeping actual local subsets gives a short degree-two/purity
   finish. Intrinsic pole pair(2,3) cannot be the sole component, so that
   sole-component reduced-ratio sector has rho>=2. Width>=6 was already
   known and is not new progress. Producer4223 exact controls pass.
5. [Hamiltonian source carriers](planar_jc48_sep06_hamiltonian.md) is
   **PROVED + INDEPENDENTLY AUDITED**. A terminal depth correction realizes
   the weight22 transport by a five-term Hamiltonian in k[p,y]. Put
   D=p^3-y^2. For fixed H, its infinitesimal source form is preserved
   exactly when S=c+D R and p divides J_(p,y)(H,S). Preservation for every
   H is equivalent to S in k+p^2 D k[p,y]. Every nonconstant generator in
   even the fixed-H class is not locally nilpotent: the invariant source
   factor D=t^3(1+x^2 t)^2 and factorial closure would otherwise force
   both x and t invariant. The displayed finite generator instead has a
   nonzero cusp-divisor residue and fails source preservation at later
   orders. These are precise stopping reasons and repaired formal spaces,
   not an obstruction to every possible polynomial repair. Producer241
   gates pass; the referee also used independent Fraction arithmetic.
6. [Explicit (4,6) family](planar_jc48_sep06_curve_probe.md) is
   **PROVED + INDEPENDENTLY AUDITED**. For
   `(U,V)=(t^4+t,t^6+t^2+lambda*t)`, the exact collision algebra has length12
   for every lambda. The disjoint exceptional sets consist of three cusp,
   four ordinary-triple and six tacnode parameters; the generic curve has
   six nodes. Retaining one exceptional point's actual fibre gives
   `1=(r-1)delta+n_z+sum_node_overlaps`. This excludes the triple and
   tangency cases as sole Keller nonproperness supports, without assuming
   a specialization identity there. The whole family is therefore excluded
   except at the three roots of `128lambda^3-288lambda-283` by this first
   argument. These have one cusp and five nodes; their necessary passport is
   `n_cusp+sum_overlaps=1` and `d<=2a`. Producer214 gates and two independent
   analytic reads pass. Item7 below closes this intermediate boundary.
   This does not classify all intrinsic (4,6) curves.
7. [Actual cusp passport](planar_jc48_sep06_cusp_passport.md) is
   **PROVED + INDEPENDENTLY AUDITED**, with classical low-degree theorems
   **CITED**. Re-access one smooth cusp point after a chosen loop and
   transport its actual retained subset and meridian together. For the
   marked ordinary-cusp generators, `A'=sigma*tau*A` and the actual cusp
   count is exactly `|A intersect A'|`. The whole-support Euler budget
   forces at most one label outside their union. A complete permutation
   argument then forces mapping degree2,3 or4, excluded by Orevkov and
   Domrina. Thus a whole irreducible A1-normalized support with exactly
   one ordinary cusp and any positive number of ordinary nodes is
   excluded. In particular **every lambda** in item6's explicit family
   is excluded as sole support. The genuine four-sheet local model is
   retained as a hostile: the classical global theorem is essential.
   Producer16815 exact gates and the independent text/source audit pass.
8. [Odd-cusp divisor spectrum](planar_jc48_sep06_odd_cusp.md) is
   **PROVED + INDEPENDENTLY AUDITED**, with primary topology and low-degree
   exclusions **CITED**. For a sole irreducible A1-normalized support with
   one analytic cusp v^2=u^m, m odd, and N>=1 ordinary nodes, the necessary
   local spectrum is d=q+n, where q>=3 divides m and n is the actual cusp
   count, zero or one. The cited degree2–4 exclusions remove q=3. With
   N>=2, n=1 and d=q+1: **d-1 divides m**, the generic actual count and
   missing length both equal d/2, and every node has empty actual fibre.
   The proof reduces the two meridians to cycles sharing one outside label;
   their product is a q-cycle and actual re-access forces the divisibility.
   A genuine abstract q=5, d=6 passport shows why the argument stops at
   higher cusps. Proper divisors and the formal q=1 purity boundary are
   retained. Producer4072 gates and the independent audit pass. This is a
   necessary spectrum, not a realization or general JC(2) exclusion.

9. [Discrete carrier rigidity](planar_jc48_sep06_discrete_carrier.md) is
   **PROVED + INDEPENDENTLY AUDITED**. Even one-way polynomial automorphism
   descent into k[p,y] forces (x,t)->(lambda^-1 x,lambda² t). Jacobian one
   leaves only the identity. The proof retains the collapsed A1 and G_m
   components; rational and nontrivial formal controls mark its boundary.
10. [Five-node global curve](planar_jc48_sep06_global_curve.md) is
    **PROVED + INDEPENDENTLY AUDITED**: Nori's actual affine-complement
    criterion gives square12>2N=10 and excludes the earlier target.
11. [Alternating monodromy and one boundary](planar_jc48_sep06_alternating.md)
    proves the all-size three-cycle support lemma and audits the actual
    degree-six A6 / one-(e,f)=(3,1)-boundary profile. Public antecedents
    for that profile are credited. The class group and canonical Weil
    class are Z[D] and 2[D]; no smoothness of an actual envelope is assumed.
12. [Linear resolution criterion](planar_jc48_sep06_resolution_budget.md)
    proves D²-2N=3e-2+2g-sum m_j on the declared actual resolution.
    Positive margin excludes whole Keller support by classical Nori.
    The explicit four-node sextic has margin zero; the ordinary cusp is
    a nonabelian equality hostile. This test alone gives no equality result.
13. [Marked infinity plumbing](planar_jc48_sep06_boundary_plumbing.md)
    is **PROVED + INDEPENDENTLY AUDITED**. For that four-node sextic,
    two tetrahedral pieces share a three-cycle and can move at most five
    labels. The actual infinity epimorphism excludes its required A6
    covering. The whole same-marked-tree/same-passport class is also
    excluded. A marked A5 image fixing one label prevents any abelianity
    overclaim. All4,000 raw assignments are checked by20,419 exact gates.
14. [Independent rational braid route](planar_jc48_sep06_boundary_braid.md)
    corroborates the same exclusion using21,503 exactly certified path
    segments and172,168 gates. Four actual loop constraints leave only
    the constant tuple among all40³ normalized three-cycle assignments.
    Rational Rouché tubes, marked transport and the monic section were
    independently audited. Numerical proposals supply no proof authority.
15. [Classical boundary surface and shears](planar_jc48_sep06_dg_surface.md)
    gives an explicit V4 with A2 complement, boundary A1, class group Z[D],
    canonical class2[D], Euler2 and a globally exact source form. The
    shear (x+f(t),t) extends exactly when f(0)=0. Every nonzero shear
    changes the old carrier, but its entire orbit algebra still collapses D.
    These surface invariants alone do not exclude a finite envelope.
16. [Actual finite quartic and relative primitive](planar_jc48_sep06_dg_finite_map.md)
    gives a finite flat degree-four map (t,b) on that surface, free global
    basis and complete discriminant16TB(1-4TB)². Interior ramification
    remains. Universally, no global g solves dt wedge dg=c omega, c!=0;
    every finite (F(t),g) map ramifies inside A2 and has even boundary index.
17. [Next actual infinity survivor](planar_jc48_sep06_next_infinity.md)
    has one(2,5) cusp, three nodes and infinity(2,11), Nori margin-2.
    Its actual marked infinity group surjects to A6 with a single-three-cycle
    meridian. Two order60 pieces can together move all six labels: this
    is a precise hostile to extending item13. Its required extension to the
    full affine complement is now **EXCLUDED** by item18; the actual boundary
    representation and geometric computation remain correct.

18. [Two actual meridians](planar_jc48_sep07_next_braid.md) are
    **PROVED + INDEPENDENTLY AUDITED** for the three-node curve. Two smooth
    critical loops, with14887 rationally certified segments, give
    `m2=m4` and `m1=m4 m3 m4 m3^-1 m4^-1`. Thus the affine group is generated
    by two positive meridians. More generally, r such meridians for a
    transitive degree-d action moving at most delta labels each require
    `r(delta-1)>=d-1`, d>1. Two generators therefore exclude every
    irreducible whole-support passport with retained count at least d/2.
    Disjoint three-cycles defeat the stronger common-fixed-label assertion
    and are retained. The exact producer passes119191 gates.
19. [Completed source image](planar_jc48_sep07_completed_response.md) is
    **PROVED FINITE-ROW + INDEPENDENTLY AUDITED**. Preserving the coordinate
    prefix through row9 reduces every polynomial carrier to exactly30
    relevant monomials of valuations6–11. Its complete image inside the
    incoming response family has five source parameters and six coordinate
    kernel directions. A six-term Rstar supplies a completed replacement
    for the high payer, while the old exact packet is outside the image.
    Every nonzero high replacement necessarily changes row10. No later
    Keller equations or polynomial termination follow. Producer4025 and
    independent Fraction reconstruction20 gates pass.
20. [Arbitrary polynomial y-linear coefficients](planar_jc48_sep07_carrier_frontier.md)
    give **PROVED + INDEPENDENTLY AUDITED** genus at least6 for every
    `I=p²Delta(A(p)+B(p)y)`, B!=0, including common repeated factors and
    zeros of B at zero. The exact formula retains their lost local
    ramification; the infinity index2 also proves geometric integrality.
    Every nonzero scalar time for a nonconstant outer Hamiltonian f(I)
    cannot have both source images rational. Producer1018 gates pass.
    Rstar in item19 is not y-linear; its separate genus analysis is item22.

21. [Whole one-(2,5)-cusp (4,6) class](planar_jc48_sep07_twofive_sextics.md)
    is **PROVED + TWO INDEPENDENT AUDITS PASS**: an irreducible curve with
    birational polynomial normalization of degrees4,6, one finite(2,5)
    cusp and at least two ordinary nodes cannot be the whole Keller
    nonproperness support. Four exhaustive infinity types7,9,11,13 use
    Nori, old plumbing, connected equisingular meridian transport, and a
    new marked infinity obstruction. The actual infinity13 census has1120
    assignments, all fixing a label. The beta=1/4 good-family hostile pays
    a lost elimination denominator. Producer3385 gates pass. This does
    not classify all rational sextics or higher finite cusps.
22. [Actual supplier genus](planar_jc48_sep07_supplier_genus.md) is
    **PROVED + INDEPENDENTLY AUDITED**: item19's six-term Rstar has generic
    invariant curve of genus27. More generally
    `I=p²y³Delta(A(p)+b y²)`, deg A=m>=4, b!=0, has genus
    `(9m+21+(m mod2)-gcd(3,m+5))/2`. Every nonzero scalar time of a
    nonconstant outer Hamiltonian f(I) is nonrational in at least one
    source image. The exact finite response remains valid. The degree3
    coalescence hostile has genus12, invalidating extrapolated24.
    Producer130 gates pass. All five finite-image parameters and
    different-invariant compositions remain outside this classification.

23. [Universal carrier genus](planar_jc48_sep07_universal_genus.md) is
    **PROVED + INDEPENDENTLY AUDITED**: every nonzero polynomial divisible
    by p Delta has every geometric generic component of genus at least two.
    One actual exceptional covering surface supplies the lower bound for
    arbitrary coefficients, degrees and repeated factors. Finite relative
    constants handle outer powers. Thus every polynomial lift in the full
    universal carrier has no nonzero rational time. This closes the precise
    local-cusp question left open by the incoming genus work;105404 gates pass.
24. [Higher odd cusp exhaustion](planar_jc48_sep07_higher_odd.md) is
    **PROVED + INDEPENDENTLY AUDITED**. For finite odd cusp order m>=7,
    the only (finite cusp, infinity cusp, node count) triples are
    (7,7,4),(7,9,3),(9,7,3). The nonzero ninth coefficient excludes m>=11.
    Connected good coefficient families retain the exceptional U=t4 chart
    and transport the actual positive meridians. Producer65 gates pass.
25. [Three actual higher-cusp braid certificates](planar_jc48_sep07_higher_braid.md)
    are **PROVED + INDEPENDENTLY AUDITED**: six rational paths with31807
    segments pass286548 gates and force two positive meridians in all
    three representatives. The finite-nine basis change is explicitly a
    basis change, not an invented projection loop. Together with item24
    and the earlier ordinary/(2,5) cases, this excludes the whole (4,6)
    one-finite-odd-cusp, at-least-two-node class as whole Keller support.
26. [Full cusp-ideal scalar-time obstruction](planar_jc48_sep08_cusp_ideal.md)
    is **PROVED + INDEPENDENTLY AUDITED** for every nonconstant
    S in K+Delta K[p,y], over every characteristic-zero field. Generic
    components have genus at least one; possible elliptic components have
    genuine vector-field poles at p=0. A commuting selfmap preserves the
    nonempty polar divisor, forcing finite order and contradicting the
    actual infinite-order scalar flow. This covers every infinitesimal
    fixed-H preserver without claiming its higher source/depth preservation.
    Literal/log comparison and finite constants are paid;6956 gates pass.
27. [DG filtration and recovered cubic exclusion](planar_jc48_sep08_dg_quadratic.md)
    is **PROVED + RECOVERED COROLLARY + INDEPENDENTLY AUDITED**. The basis
    x^a t^(n-j)(1+x²t)^j, a<=2n,j<=n, is exhaustive in every degree.
    THM-2063/2071/2118's full one-coordinate theorems imply that a global
    Keller pair could have no nonconstant pencil member of source degree
    at most three: its polynomial inverse would extend x=1/r across D.
    This credits recovered canon and the incoming source-linear theorem;
    it is not a new low-degree planar theorem. Producer122 gates pass.
28. [Two-cusp sextic](planar_jc48_sep08_two_cusps.md) has **PROVED +
    INDEPENDENTLY AUDITED** geometry: U=t4-2t²,
    V=t6-3t4/2+t3/3-t, two ordinary cusps, four nodes and infinity(2,9).
    The unique smooth vertical fold supplies at most three positive
    meridians. Item31 now proves its full group is Z and excludes it as
    whole Keller support; the original geometry bundle is retained unchanged.
    The correct conditional ledger is1=-a+n1+n2+sum omega, so the
    one-cusp half-retained bound cannot simply transfer. Producer42 gates pass.

29. [Base-only quartic Jacobians and global approximate root](planar_jc48_sep08_quartic_boundary.md)
    is **PROVED + INDEPENDENTLY AUDITED**. The inherited finite-pole proofs
    need only a nonzero base-dependent Jacobian, not its invertibility at
    the tested place. Applied to both full DG charts, including lambda r²,
    they force F=H²+L with global H in L2 and nonconstant L in L1 for any
    hypothetical quartic pair. Neither auxiliary is asserted to inherit a
    Keller mate or to lie in the original pencil. Producer126 gates pass.
30. [Uniform two-cusp actual-sheet passport](planar_jc48_sep08_two_cusp_passport.md)
    is **PROVED + INDEPENDENTLY AUDITED**. Two ordinary cusps and N>=2
    nodes force d in {2a-1,2a} and every deleted node overlap zero. Equality
    in the cusp injection identifies all fixed labels as actual and forces
    even nontrivial cycles; commuting-node parity pays the otherwise open
    N=2 case. With three positive meridians, cited degrees2–4 and the new
    degree5 edge obstruction leave d>=6. The d6 A6 passport is an exact
    abstract survivor, not a representation of the actual curve group.
    Producer1051 gates and an independent unfiltered node-parity check pass.
31. [Actual two-cusp affine braid group](planar_jc48_sep08_two_cusp_braid.md)
    is **PROVED + INDEPENDENTLY AUDITED**. Four common-base loops with
    77174 rational segments force all four positive meridians equal in an
    arbitrary group; defining-polynomial winding proves the group is Z.
    A transitive cyclic cover of degree>1 cannot retain a smooth sheet, so
    this literal sextic is excluded as whole Keller nonproperness support.
    Both co-projected cusps and all critical projection values are retained;
    the unused two complex-node scout words remain HEURISTIC. Producer
    and independent referee pass694669 gates in each normal/optimized run.
32. [Entire boundary-linear carrier](planar_jc48_sep08_boundary_linear.md)
    is **PROVED + INDEPENDENTLY AUDITED**. Every global F of degree<=1
    in r has form rb C(b)+B(b), and no polynomial source mate exists in
    any degree. A twice-F derivative and rational-map degree force one
    denominator root; source criticality reduces its multiplicity to one.
    In the rationally integrable case, two reduced components of one fibre
    make the required double-pole repair incompatible. The quartic
    F=b4-v has a global root b², nonconstant boundary restriction, no
    affine-source critical point, and a rational mate, yet no polynomial
    mate. It is critical at one point of the added boundary. Producer149
    gates and independent analytic/source audit pass.

33. [Quartic relative differential](planar_jc48_sep08_quartic_differential.md)
    is **PROVED + INDEPENDENTLY AUDITED**. For the global square prefix,
    octic boundary roots of multiplicity at most two, disjoint from the
    linear remainder's boundary roots, make the normalized generic relative
    form holomorphic and exclude a global mate. Fully global multiplicity3
    and6 controls exhibit respectively a logarithmic residue and residue-free
    double poles. Thus disjoint boundary roots alone do not imply local
    regularity. Producer78 exact gates and independent audit pass.
34. [Actual unit torsion and connection](planar_jc48_sep08_unit_torsion.md)
    is **PROVED + INDEPENDENTLY AUDITED**. For a!=0 and every h in b³C[b],
    F=arb+c+h(b) is smooth on the affine source. The actual unit class is
    -a/2 times the component-labelled (F-c)^-2 principal part: its exact
    primary order is2, and its j-th connection derivative has order j+2.
    The polynomial primitive of (F-c)² and all special-fibre components
    are retained. The h=0 antecedent is credited. This pays a fixed
    polynomial for the torsion lane, not its map from moving source
    responses. Producer86 gates and independent audit pass.
35. [Componentwise primitive pole degree](planar_jc48_sep08_pole_degree.md)
    is **PROVED + INDEPENDENTLY AUDITED** for rational mates. If F|D is
    nonconstant, a generic transverse D point forces a primitive to have
    local degree3. Its degree equals its total boundary pole order on
    that same compact component. The full multiplicity6 hostile has
    normalized differential orders (-2,-2,0,3,3), so pole degree at most2
    excludes even a rational mate despite vanishing residues. Its affine
    critical points already exclude polynomial mates; that weaker fact is
    not new. The same restriction argument strengthens item33 to rational
    mates. Producer68 gates and independent audit pass.
36. [Complete two-cusp parameter strata](planar_jc48_sep08_two_cusp_family.md)
    is **PROVED + INDEPENDENTLY AUDITED**. Every birational polynomial
    (4,6) curve with two ordinary cusps and otherwise ordinary nodes is
    represented in an explicit three-parameter family; infinity types
    (2,7),(2,9),(2,11) are exhaustive, with5,4,3 nodes respectively.
    The complete good strata are connected, retaining projection-degenerate
    parameters s=+/-1. Strict finite braid certificates and relative marked
    resolutions give actual affine group Z for every curve in the first
    two strata, excluding each as whole Keller support. A good infinity11
    representative is exact, but its actual group remains outside this
    checkpoint. Producer77 gates and independent audit pass.

37. [Exact moving-centre tubes](planar_jc48_sep08_moving_tubes.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. Expanding along
    affine root centres retains cancellation of common motion before taking
    coefficient bounds. Exact pair separation, concentric endpoint gluing
    and a common-base configuration homotopy pay the actual braid. One full
    inherited9553-segment path is recertified in92 moving segments with the
    identical word. Long translation, swapped labels and a curved root with
    vanishing first residual are controls. Producer10324 gates and an
    independent source/method audit pass; no old witness is changed.
38. [Quartic common-boundary-root necessity](planar_jc48_sep08_quartic_common_root.md)
    is **PROVED + INDEPENDENTLY AUDITED** for arbitrary rational mates.
    If the actual binary octic N|S and quartic M|S have no common zero,
    every boundary branch is regular except a controlled double-root
    Newton face. Its generic Morse split bounds the total primitive pole
    degree by2, across all octic multiplicities. Disjointness also forces
    F|D nonconstant, supplying local degree3 on the same compact component.
    Thus any rational mate requires a shared root. The actual h4+h rational
    mate with nonzero sections and nonconstant F|D attains degree3 and marks
    this condition's boundary. The degree principle's canonical antecedents
    are credited. Producer568 gates, a complete root audit and a second
    independent analytic read pass. This does not identify W with an
    actual Keller envelope.
39. [Infinity-eleven actual braid and whole two-cusp class](planar_jc48_sep08_infinity11_braid.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. Five exact
    common-base moving witnesses,736 segments in total, force all four
    positive meridians equal; winding gives actual affine group Z. The
    good infinity11 supplier and every node/cusp/infinity branch are
    independently checked. Item36's connected marked-stratum transport
    now closes infinity11 as well. Therefore the entire birational
    polynomial (4,6) class with exactly two ordinary cusps and otherwise
    ordinary nodes has group Z and is excluded as whole irreducible Keller
    nonproperness support. Producer and independent referee pass6793 gates
    per mode; an independent crossing extraction recovers all five words.

40. [Shared-root first jets](planar_jc48_sep08_shared_roots.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. A simple octic
    root or a unit first normal derivative makes every shared-root branch
    regular. At a finite singular shared point, any nonzero quadratic
    tangential coefficient of N, first tangential coefficient of its normal
    derivative, or first tangential coefficient of M produces a genuine
    logarithmic residue. Thus a rational mate requires m>=3,j>=2,n>=2 at
    every remaining such point. If all octic multiplicities are at most2,
    no rational mate exists for arbitrary M, even when F|D is constant.
    The full infinity factor and the nonzero-M degree3 rational hostile
    are retained. Producer454 gates and independent analytic/source audit
    pass; this is a necessary jet condition, not a shared-root closure.

41. [Nonabelian six-word controls](planar_jc48_sep08_three_cusp_group.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED** as group algebra.
    The literal scout presentation surjects onto B3; an S4 control shows
    that neither cyclicity nor generation by any one cusp pair follows.
    Its three cusp images are proper S3 subgroups of the full image.
    All geometric words remain HEURISTIC in this bundle. Producer115
    gates and an independent full group/source audit pass.
42. [Constant-D finite-six quartics](planar_jc48_sep08_const_d_quartic.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. If N|S=a x6,
    M(0,0)!=0 and F|D is constant, no rational mate exists. The complete
    six-parameter numerator and five-parameter remainder include all lower
    normal jets. Four actual infinity zeros replace the unavailable D
    supplier; the last residue-free layer is an explicit elliptic curve.
    Its entire two-pole function space is generated by1 and1/w, and the
    only primitive candidate fails the original fibre-value coefficient
    comparison. Producer97 gates and an independent full audit pass.
43. [Every finite location of the sixfold point](planar_jc48_sep08_const_d_translation.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. The preceding
    exclusion holds for N|S=a(x-p)^6 at every finite complex p. Exact
    coefficient boxes and constant unit minors show that the full family,
    after source-plane translation, belongs to the auxiliary fixed-zero
    global family. The translation preserves the rational Jacobian while
    changing boundary values; its explicit -2p/r pole forbids asserting
    a surface automorphism. Producer51 gates and an independent audit pass.
44. [Exact D4 and conditional degree-eight floor](planar_jc48_sep08_three_cusp_passport.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**, with geometric
    application **CONDITIONAL**. Reversible conjugations identify the full
    six-word presentation with Artin D4, retaining the original local
    pairs for sheet counts. Every transposition image moves at most four
    labels; the three-cusp Euler ledger, saturated node parity and
    three-cycle support geometry exclude mapping degrees5,6,7 whenever
    the declared common-access geometry is supplied. Cited degrees2--4
    then give degree>=8 for N>=3 nodes. The S4 actual-ledger hostile
    explains the degree-four consumer. Producer789 gates and full audit pass.
45. [Finite D4 fixed-sheet Euler classifier](planar_jc48_sep08_d4_involution.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED** on the literal
    signed-permutation group of order192. A central involution and one
    retained fixed letter reduce every possible Euler-one transitive
    action to26 stabilizers; only degrees1 and4 survive. A separate
    signed-kernel/lift and character calculation reproduces all26 profiles.
    Its Keller consumer needs the actual six pairs, a finite-group entry,
    and equality of retained sets with full fixed sets. Saturation does
    not force involutivity. Producer5631 gates and independent audit pass.

46. [Actual three-cusp braid and local pairs](planar_jc48_sep08_three_cusp_braid.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. Six common-base
    moving witnesses have807 segments. Eighteen whole complex-disk root
    clusters pay the actual cusp and node pairs; two inside-pair Hurwitz
    changes preserve the required actual counts. The literal curve has
    an actual marked D4 quotient and mapping-degree floor8 as whole
    Keller support. Equality with D4 is not asserted. Producer9534 gates
    and independent original-monomial/word/marking reconstruction pass.
47. [Complete fixed-zero octuple class](planar_jc48_sep08_finite_eight.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. For global
    F=H²+L with N|S=alpha*x8, alpha!=0, no polynomial mate of any degree
    exists. Full section-space entry, all normalized j2 branches, actual
    D pole capacity, source criticality, weighted quadratic field trace
    and the boundary-linear theorem exhaust the coefficients. Genuine
    rational mates remain inside this polynomial-only theorem. Producer92
    gates and an independent22-check section/trace reconstruction pass.
48. [Whole ordinary three-cusp family](planar_jc48_sep08_three_cusp_family.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. The complete
    birational (4,6) class with three ordinary cusps and otherwise nodes
    has only infinity7/9 and four/three nodes. Proper intrinsic good loci,
    strict local-certificate persistence and marked resolution transport
    pay the D4 quotient and degree floor8 throughout both good strata.
    Harmless cusp co-projections and a vanishing intermediate infinity
    coefficient are retained. Producer102 gates and an independent full
    audit with23 additional original-integral/series checks pass.
49. [D4 with arbitrary retained subsets](planar_jc48_sep08_d4_retention.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. Cusp and node
    deficit bounds replace the full-fixed-set hypothesis. All53 literal
    reflection-containing stabilizers and211 retained-size trials force
    degree4 whenever the Euler lower bound is at most one; extra nodes
    contribute nonnegative overlap. A separate full signed-kernel/lift
    audit reproduces every subgroup and count. The cited Coxeter supplier
    and classical degree-four exclusion therefore rule out involutive
    meridians for every actual marked three-cusp curve in item48, without
    identifying retained and fixed sheets. Producer12509 gates pass.

50. [Every finite octuple location](planar_jc48_sep08_octuple_transport.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. Corrected actual
    global functions h=(x-p)^4t+(x-p)^2-2p(x-p) and v=h/(x-p) preserve
    the exact rational field and volume coefficient. Weighted quadratic
    trace and the complete three-dimensional genus-two primitive space
    close the new rational cases. A same-fibre pole comparison on the
    original source closes the remaining polynomial case. Every finite
    p is covered, while rational mates remain genuine hostiles to a
    stronger conclusion. Producer85 gates and independent full audit pass.
51. [No single-cycle meridian and degree floor10](planar_jc48_sep08_cycle_support.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. Braided single
    cycles overlap in at least half their support; commuting leaf cycles
    have identical or disjoint supports. Actual retention deficits then
    force an Euler contradiction at every length, except the classical
    degree-four control. The independent involution exclusion also closes
    degrees8,9. A transitive eight-letter order-four control has Euler2,
    so sharp support and transitivity alone still fail. Producer373 gates
    and a separate integer-ledger/cycle-type audit pass.

52. [All-degree marked three-cusp closure](planar_jc48_sep08_three_cusp_closure.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. The entire
    ordinary three-cusp (4,6) class is excluded as whole irreducible
    Keller support. For arbitrary cycle type, actual retained deficits,
    central support incidences, leaf unions and the complete node Euler
    ledger reduce moved counts to2,4,6. The proved small-cycle suppliers
    close2,4; the last formal D15/k7 row requires three singleton node
    overlaps, incompatible with commuting supports. Both good family
    strata are paid by item48. Producer43138 gates and an independent
    full analytic/source audit with separate universe counts pass.
53. [Infinity octuple with its actual weight](planar_jc48_sep08_infinity_octuple.md)
    is **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. The actual
    surface inversion changes volume to u² du wedge dT. Complete local
    degrees2/3/4 and the last balanced Morse layer are retained; the
    remaining quadratic field gives a genus-two holomorphic form or a
    nonzero logarithmic residue. The final composite case is excluded
    in the original polynomial ring. With item50 this closes every
    single-point boundary octic, but rational hostiles remain. Producer141
    gates and ten independent old-basis/field checks pass.

54. [Distinguished two-fourfold roots](planar_jc48_sep08_four_four.md) is
   **PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED**. For the complete
   global class `F=H²+L`, `N|S=alpha*x⁴`, no polynomial mate of any degree
   exists. Every branch at infinity has a regular differential after the
   actual volume weight. The finite branches force the complete residue
   identity `g=C f′`; all survivors have actual source critical points.
   The nonconstant-L example `F=w⁴+w`, `G=1/[x(4w³+1)]`, `w=x²t`,
   preserves the rational-mate boundary. All96 normal/optimized gates and
   the independent analytic/source audit pass. Item56 separately pays every finite location without an assumed surface
   symmetry.

55. [The mixed-cusp boundary](planar_jc48_sep08_three_cusp_boundary.md)
   is **PROVED + INDEPENDENTLY AUDITED**. The three bad-jet lines have
   exactly three quadratic-cover exceptions; every birational good member
   has finite cusp types(5,3,3), infinity7 and three nodes. Invertible
   affine relabelings reduce them to one connected intrinsic good family,
   retaining the harmless s=0 projection collision. All122 exact gates,
   both literal controls and the minimal braid-five hostile pass.
56. [Every finite/infinity four-four pair](planar_jc48_sep08_four_four_transport.md)
   is **PROVED + INDEPENDENTLY AUDITED** to have no polynomial mate of
   any degree. Full shifted global sections preserve the actual residue;
   its survivors give a critical quintic. Two leading coefficients exclude
   all21 forbidden-root factorizations, and parity excludes all four cubic
   cases. All151 gates and an extra original-source critical control pass.
57. [H4 single cycles](planar_jc48_sep08_h4_single_cycle.md) is **PROVED**:
   four equal-length single nontrivial cycles obeying the Artin H4
   relations are all equal. Odd-braid cycle runs and saturated ordinary
   alternation force three disjoint half-size images into one support,
   a contradiction unless every support coincides. A centralizer argument
   then gives equality. All15,111 controls and independent audits pass;
   transitivity therefore precludes any retained fixed sheet.
58. [Actual mixed-cusp H4](planar_jc48_sep08_mixed_cusp_braid.md) is
   **PROVED + INDEPENDENTLY AUDITED**. Six new rational moving witnesses
   with1,517 segments and17,100 exact gates pay all complex-disk clusters,
   positive meridian pairs and the marked Artin H4 quotient. The three
   cusp edges have labels3,3,5; nodes give the complementary commuting
   pairs. Item55 transports this quotient through the whole intrinsic
   good family. This is a quotient conclusion, not a free-lasso or
   fundamental-group equality claim.
59. [H4 reflection retention](planar_jc48_sep08_h4_involution.md) is
   **PROVED + CITED classical group order + INDEPENDENTLY AUDITED**.
   Every nontrivial transitive Coxeter H4 action has reflection fixed ratio
   at most4/15. An exact120-root image of order14400 and all221 closed
   reflection sets recover the bound16/60; an independent published-table
   fingerprint agrees. The actual three-node Euler identity then gives
   `1>=13D/15>1`. All67,850 gates pass. Together with item57, this excludes
   every involutive or single-cycle meridian type for the actual mixed
   family, in every degree. Other mixed cycle types remain open.

60. [Every two-fourfold-root octic](planar_jc48_sep08_two_finite_four_four.md)
   is now excluded for **polynomial mates in every degree and at every
   boundary location** in the stated full global quartic class. The new
   two-distinct-finite-root theorem excludes even rational mates, including
   constant L: an active root gives a degree-four first coefficient that
   cannot be a constant multiple of the cubic fibre derivative. Its
   individual residue numerator has leading32alpha⁴(p-q)¹⁵. All60 exact
   gates and a full independent audit pass. The finite/infinity supplier
   pays the remaining locations and retains its rational-mate hostile.
   Other multiplicity partitions and unrestricted JC(2) remain open.

61. [Odd-cusp arbitrary retention](planar_jc48_sep08_odd_retention.md)
   is **PROVED + INDEPENDENTLY AUDITED**. An odd braid of length2r+1
   gives `(r+1)n>=(2r+1)k-rD` and the stronger one-deficit bound
   `n>=k-floor(rt/(r+1))`, for actual re-accessed subsets. Full support
   incidence and the Euler budget exclude every D<=11 in the marked
   mixed-cusp H4 family. Its mapping degree is therefore at least12.
   All34,983 gates pass; independent conjugacy-weighted enumeration
   recovers the complete finite universe. An unbounded cardinality-only
   survivor records why cycle order is the next missing information.
62. [Leading-coefficient exactness](planar_jc48_sep08_leading_exactness.md)
   is **PROVED + INDEPENDENTLY AUDITED** in every degree: if a rational
   mate exists for F of t-degree n and leading A, then dx/a has a
   primitive in C(x)(a), a^n=A. One formal Laurent coefficient proves
   this, without convergence or a mate-degree bound. It is an iff for
   the pure-monomial model only. All241 gates and the necessary-only
   hostile pass. In the actual global quartic class, every finite double
   root of the boundary octic is forbidden. Hence all(4,2,2) placements
   exclude even rational mates; a lone double point at infinity is not
   covered by that finite-root statement.
63. [Complete radical differential classification](planar_jc48_sep08_boundary_exactness.md)
   is **PROVED + INDEPENDENTLY AUDITED** for every nonzero polynomial
   N of degree at most eight. All67 multiplicity partitions reduce to
   the displayed17 possible exact types, with precise position equations
   for four exceptional types. The5+1+1+1 elliptic case is exact iff
   its transformed cubic has no linear or quadratic coefficient; the
   complete primitive space is paid. All195 gates and independent
   partition/Fraction checks pass. This supplies the complete leading
   necessary test, not a full two-variable mate classification.

64. [Every5+1+1+1 boundary placement](planar_jc48_sep08_five_one_one_one.md)
   excludes polynomial mates of every degree in the full global quartic
   class. The next formal coefficient forces the complete degree-four M
   row to have five counted zeros, hence L is constant and the polynomial
   Jacobian has a nonconstant factor. All62 gates and a full independent
   audit pass. The rational boundary is sharp inside this exact partition:
   an explicit global H has rational mate1/(3delta*h), and H²+c inherits
   a rational mate. Its degree-three coordinate subfield is retained.
65. [H4 mixed(3)(2) equality](planar_jc48_sep08_h4_mixed32.md)
   is **PROVED + INDEPENDENTLY AUDITED** in every ambient size. Each pair
   reduces to at most ten moved labels; all5040 partners give complete
   cycle-intersection tables. Ordinary triples and complementary
   commutations force a common transposition, then the proved single-cycle
   theorem forces all four generators equal. All51,345 gates and a separate
   literal-action census pass. The minimal five-letter hostile records why
   taking powers does not generally preserve the fifth braid relation.
66. [Exact pencils and full quadratic rational mates](planar_jc48_sep08_exact_pencils.md)
   are **PROVED + INDEPENDENTLY AUDITED**. Exactly six two-dimensional
   spaces of degree-at-most-eight polynomials have generic exact radical
   differentials; every nonzero member is exact. Universal primitives
   retain the pencil parameter rationally. The full quadratic discriminant
   then gives an iff criterion for rational mates of every global DG
   quadratic, including its proportional boundary. Generic geometric
   genus is at most one, and elliptic examples occur in the actual global
   ring. All64 gates pass. This also classifies rational mates for constant-L
   global quartics; no nonconstant-L or globally regular mate is inferred.

67. [Sharp genus capacity in every degree](planar_jc48_sep08_genus_capacity.md)
   is **PROVED + INDEPENDENTLY AUDITED**. Exact dx/sqrt(N) has genus
   at most floor((n-4)/4) for even degree n>=4; odd degree forces a
   pure root. The bound is attained in every even degree. In degree4m,
   maximum genus m-1 forces exactly
   `N=u^(2m+1)(a u^(2m-1)+b)`, a,b nonzero, with its explicit primitive.
   The complete odd pole space pays every position equation. All5081
   gates and independent additive partition DP pass. Exact genus-two
   degree-twelve examples retain the degree boundary of item66.
68. [H4 mixed(3)(3) equality](planar_jc48_sep08_h4_mixed33.md)
   is **PROVED + INDEPENDENTLY AUDITED** in every ambient size. The
   complete ordinary-pair joint orbits and semiregular centralizers force
   the first two generators equal; complementary commuting/odd relations
   then give equality of all four. All448,965 gates and an independent
   census recover all36,960 pairs and61 centralizer orbits. Equal-length
   cycle mixing is explicitly retained, with six-letter hostiles to a
   falsely assumed common two-block decomposition. Thus this type is
   excluded for the actual mixed-cusp covering in every degree.

69. [All-finite quartic(4,3,1)](planar_jc48_sep08_four_three_one.md)
   is **PROVED + INDEPENDENTLY AUDITED**. A rational mate forces L
   constant; all polynomial mates are excluded in the complete global
   square-prefix class, for every finite root placement and any mate
   degree. Moving residues, an order-four zero, geometric integrality
   and the entire genus-one primitive space leave an exact coefficient
   contradiction. The anchored degeneration retains its extra parameter.
   All107 gates and independent original-source reductions pass. Actual
   rational constant-L mates in the same partition exist for every p.
70. [Every graph-complement surface W_m](planar_jc48_sep08_dg_genus.md)
   has its complete filtration, sharp quadratic rational-mate genus bound
   m-1, and an exact extremal-pencil iff **PROVED + INDEPENDENTLY AUDITED**.
   Explicit rational mates attain equality for every m. The fixed high
   root and all coefficient descent are paid, as is the field degree
   2m-1 of the sharp-family map. The stronger polynomial-mate exclusion
   is explicitly a recovered THM-2071 corollary. All1550 gates pass.
71. [All-exponent binomial exactness](planar_jc48_sep08_binomial_exactness.md)
   is **PROVED + INDEPENDENTLY AUDITED**. For nonzero a,b and integers
   A>=0,r>=1, du/sqrt(u^A(a u^r+b)) is algebraically exact precisely
   at(0,1) or A=r+2+2r ell, ell>=0. The complete affine primitive ring
   and exponent classes prove the iff; an explicit recurrence pays
   coefficient descent. It gives whole sparse pencils and a sharp
   same-curve/different-differential hostile. All1654 gates and1170
   independent original-coordinate Laurent systems pass. Classical
   binomial integration is credited; no external priority is asserted.

72. [Fourfold infinity in(4,3,1)](planar_jc48_sep08_four_three_one_infinity.md)
   is **PROVED + INDEPENDENTLY AUDITED**. Complete weighted infinity
   regularity and primitive capacity force constant D; the reduced
   M/N has a compulsory residue unless L is constant. Together with
   item69 and the odd-degree leading table, this closes all boundary
   locations for polynomial mates. Literal quadratic trace and all-p
   rational constant-L controls provide separate sharp checks. All87
   gates and the independent full local/global audit pass.
73. [H4 mixed(4)(2) equality](planar_jc48_sep08_h4_mixed42.md)
   is **PROVED + INDEPENDENTLY AUDITED** in every ambient degree. The
   unique cycle blocks, exact pair inventories and the absence of an
   allowed semiregular restriction to six labels force a common
   transposition; the paid single-cycle theorem then gives equality.
   All343032 gates and an independent83160-pair census pass. The two
   minimal six-label odd-pair hostiles preserve the need to prove the
   common factor before removing it. All meridian types moving at most
   six sheets are now excluded for the actual mixed-cusp covering.

74. [Complete shared constant-D finite-six](planar_jc48_sep08_constant_d_shared_six.md)
   is **PROVED + INDEPENDENTLY AUDITED**, with an exact five-case
   rational-mate existence table. Together with the earlier M-unit
   supplier it closes the entire constant-D finite6/infinity2 polynomial
   class for every finite location. Conic logarithms, a residue forcing
   an actual critical line, and same-fibre pole repair pay all cases.
   Actual nonconstant-L rational submersions remain as sharp hostiles.
   All68 gates and independent original-source identities pass.
75. [Actual mixed-cusp degree at least sixteen](planar_jc48_sep08_h4_degree16.md)
   is **PROVED + INDEPENDENTLY AUDITED**. The small-type closures force
   moved support at least7. The complete scalar head leaves only three
   rows below16; actual complementary supports and the two one-sheet
   deficits in the partial-retention row contradict Euler value1.
   All561 gates and an independent scalar-existence census pass.
   The necessary scalar rows at16 are retained without claiming a
   realized permutation passport or Keller map.

76. [All-finite5+3 quartic boundary](planar_jc48_sep08_five_three.md)
   is **PROVED + INDEPENDENTLY AUDITED**: a rational mate forces L
   constant, and every polynomial mate is excluded. Moving-root and
   higher inverse residues reduce three coefficient strata to complete
   genus1,2,3 primitive spaces. The genus3 repair is regular at the
   ordinary affine points as well as the boundary. All117 gates and
   independently reconstructed matrices pass. Actual all-p rational
   constant-L examples preserve the precise boundary; either infinity
   placement remains separate open work.

77. [Every6+1+1 boundary placement](planar_jc48_sep08_six_one_one.md)
   is **PROVED + INDEPENDENTLY AUDITED** for polynomial exclusion.
   In the all-finite class, rational mates force L constant through
   a moving-centre logarithm, an actual quadratic trace, higher
   inverse residues and an explicit Euclidean contradiction. Both
   infinity placements fail leading exactness. All63 gates and an
   independent81-check formal-residue/Bezout reconstruction pass;
   all-p rational constant-L examples retain the precise boundary.
78. [Every3+3+2 boundary placement](planar_jc48_sep08_three_three_two.md)
   is **PROVED + INDEPENDENTLY AUDITED**. Only the double-infinity
   location can have rational mates; one active triple forces
   constant D, genus one and a complete two-dimensional primitive
   space. The formal generic bracket excludes nonconstant L.
   All56 gates pass, including a special level where the bracket
   becomes constant and cannot replace the generic fibre.
79. [Every4+2+1+1 boundary placement](planar_jc48_sep08_four_two_one_one.md)
   is **PROVED + INDEPENDENTLY AUDITED**. The leading midpoint
   condition leaves a full family with genus two in both distinct
   infinity regimes. Its complete basis1,1/u,H/u has two successive
   bracket obstructions. All41 gates and independent original-chart
   identities pass. Every polynomial mate is excluded; rational
   constant-L examples realize the remaining position conditions.

The [347-artifact manifest](planar_jc48_sep06_manifest.json) pins79 proof notes,
80 audit notes,81 sources,81 frozen outputs and26 compressed witnesses.
The81 programs report2524164 exact gates per complete normal or optimized
run. Each frozen output agrees in both modes. All335 inherited pins remain
unchanged. These are scoped proof controls, not a census of Keller maps;
no external priority claim or full JC(2) result is asserted.


## Literature and correction recovery

The [focused source sidecar](../reference/CORE-PAPERS-PLANAR-JC.md) records
Bonnet's relative-exactness precedent, the exact public Paper II/III
antecedents, and the Jelonek version guard. The current arXiv record of
2011.03472 is withdrawn; correct version3 and the published weaker
statement are the inputs. The withdrawn stronger connectedness claim is
not used. The nodal argument is credited as recovered work, after the
independent referee found the direct antecedent beyond the first overview.

Methods used: Search the statement before the method; Compute the repair
quotient before testing the residual defect; Separate descent, ambient
scale and regularity debt. No new META-PATTERNS card is promoted on this
single checkpoint. Recovered results are retained with their scope and
provenance, rather than counted as new exclusions.

## Immediate continuation

The moving-prefix calculation has first unretained background terms at
bracket row15 and defect row16. The exact Hamiltonian mechanism is now
known, and its two distinct failures are proved. Do not search for a
locally nilpotent Hamiltonian inside the fixed-H carrier: that class is
excluded. The new discrete theorem also excludes every nonidentity
polynomial symplectic automorphism with even one-way preservation of the
whole carrier. A productive rational-time repair must escape the entire
polynomial cusp ideal or use a composition requiring a new invariant
argument. Preserving only one fixed source does not escape item26's
obstruction. A successful finite table is not polynomial termination.

At infinity the exact family and ordinary-cusp-plus-nodes class are closed.
The higher-cusp spectrum is proved. Both the five-node target and its four-node equality successor are now
closed. The subsequent three-node target is now closed as well: its necessary
d=6, a=delta=3 passport fails the two-meridian graph bound.
Keep the actual subset, conjugating loop and fibre counts visible. Its
abstract passport needs a full complement representation, an actual affine
source and polynomial completion before it becomes a Keller candidate.
With several nonnodal points the one-point positivity argument changes;
retain their number and all actual special fibres before generalizing.

Incoming commit efdf4fe524 was read and integrated before checkpoint
03d16f65c7. Its [joint interval compiler](continuing5_20260906_synthesis.md)
retains a common translate where independent extrema lose it. The parallel
method here is to retain `A'=gA` when comparing two actual sheet subsets.
This is a connection between proof operations, not a map between LRC rows
and Keller maps: both require the common choice before evaluating a joint
predicate. The four-sheet cusp hostile is the decisive test here.

The subsequent incoming b6fddf843d synthesis was also read. It identifies
an exact shared algebraic operation: the norm of p^3+2 in our rank-six
pair algebra is the cusp polynomial `128lambda^3-288lambda-283`.
This unit test detects diagonal collisions, but its lambda=-1 and2
hostiles retain a triple point and a tacnode. The complete incidence and
actual-fibre data remain necessary. That incoming norm interpretation
agrees with our certified resultant identity and does not change the
scope of the later cusp exclusion.

### Incoming response and completion results

Incoming `1100abe5ba` was read and independently compared by two agents
before the integration checkpoint. Its [earlier-memory theorem](planar_jc_long_20260906_memory_earlier.md)
classifies the complete valuation-at-least-eleven response through row15:
weight22 gives G_m x A^9 with separate A^10 coordinate fibres; weight21
still fails in that declared horizon. Seven equations retain the delayed
odd channel. A polynomial inverse with square-zero correction proves a
constant response quotient at every background specialization while
keeping that background in the actual lifts. This strengthens the finite
classification around our transport without asserting later-row completion.

The [maximal depth carrier](planar_jc_long_20260906_depth_carrier.md) is
exactly k+p²k[p,y]. Its intersection with preservation of even one affine
source is precisely our universal carrier k+p²Delta k[p,y], where
Delta=p³-y² is a function, not the boundary divisor D. The
[completed Hamiltonian theorem](planar_jc_long_20260906_hamiltonian.md)
integrates this smaller carrier at scalar times on the actual Delta-adic
completion. It does not identify that whole completion with a Laurent
chart or make completed coefficients polynomial.

The [nonrational-time theorem](planar_jc_long_20260906_nonrational.md)
excludes every nonzero rational scalar time for S=f(p^a Delta), a>=2,
f nonconstant, using a generic invariant curve of genus ceil(a/2)+1.
This includes our named formal S=p²Delta control and is stronger than
non-local-nilpotence. Its rational u-flow hostile is exactly our discrete
carrier control. Our polynomial DG shears are compatible: their
Hamiltonians F(t) are outside k[p,y] and their invariant has genus zero.
Arbitrary R, compositions and changed carriers remain outside that theorem.

That concrete source image test is now paid by item19:
`R -> pi15{-u/2+H_xi,p²Delta R}` retains every induced coordinate correction
and original bracket/depth/prefix equation. Its five-dimensional image
contains a replacement of the high payer; polynomial specialization and
later equations remain unpaid. Our torsion connection requires an
actual fixed smooth polynomial and labelled principal parts. Neither
Delta, the background xi, nor the square-zero lift matrix is automatically
that connection.

Incoming `7c557fb04c` and its [continuing8 synthesis](continuing8_20260906_synthesis.md)
were also read. Its exact source-chart transfer boundaries preserve the
operator, characteristic, lattice and zero tail separately; they give no
characteristic-zero Keller obstruction. This agrees with the distinct
regularity and carrier hypotheses above. No frozen artifact was changed
by these incoming scope comparisons.

### Closed targets and the next geometric scale

The curve `(t^4+t^2,t^6+t^5+t^2)` is **EXCLUDED** as whole Keller
nonproperness support by Nori: its resolved square12 exceeds twice its
five-node count10. The curve `(t^4+t^3+t^2,2t^6+3t^5+2t^3+2t^2)` is
also **EXCLUDED**, despite Nori equality: its marked infinity group
cannot supply the actual six-sheet representation. Exact global braid
constraints independently corroborate this second closure.

The [third actual target](planar_jc48_sep06_next_infinity.md), now closed, is

```text
U=t^4+t^3+t^2,
V=16t^6+24t^5-19t^3-19t^2.
```

Its **audited geometry** is one affine(2,5) cusp, three ordinary nodes,
and infinity(2,11). The resolution cost is18 and D²=4, giving Nori
margin-2. Its actual marked infinity tree has an explicit surjective A6
representation with a single-three-cycle meridian. This boundary witness
**does not extend** to the required representation of the full affine
complement: the two exactly certified smooth-critical loops force two
positive meridian generators. The boundary/affine distinction supplied
the successful operation rather than an untracked numerical heuristic.

The [whole one-(2,5)-cusp (4,6) class](planar_jc48_sep07_twofive_sextics.md)
is now **EXCLUDED** as whole Keller support.
Its normal forms have exactly four infinity types7,9,11,13. The first two
inherit established exclusions; the fourth has a new marked A6 obstruction;
the third transports the actual two-meridian certificate through a connected
equisingular family. Independent geometric and full-proof audits both pass.
Items24-25 now close every higher finite odd cusp in the same degree pair
and node range. Multiple nonnodal points, fewer nodes and higher degree
normalizations remain separate targets.

In parallel, the explicit DG surface's boundary separator supplies a
real finite map but forces retained-source ramification in the proven
first-coordinate class. A new pair must escape both the orbit carrier
and C[t] as its first-coordinate restriction. Neither condition alone
pays the bracket-one equation, source entry or the actual finite envelope.

Incoming fe3ac58164's continuing7 synthesis was read before this work.
Its full singular-fibre and norm controls reinforce the need to preserve
actual fibres before taking a scalar quotient. Our shared pair equation,
with a changing residual polynomial, is a concrete internal transport;
no LRC-to-JC implication is claimed. The board remains the current truth
source; earlier scouting targets are superseded by these closures.

The heartbeat now runs this board every30 minutes through the stated
48-hour cutoff. It stays quiet unless there is substantive progress, a
correction, completion, failure or required user action. At the cutoff,
push the final coherent checkpoint and pause the heartbeat.


### Incoming genus work and cross-lane reassessment

Incoming `a1f02fa5ef50` was read before this checkpoint. Its
[carrier genus theorem](continuing9_20260907_flow_carrier_genus.md) treats
univariate multipliers, all affine multipliers, and coprime monomial outer
invariants; its [constant-slope extension](continuing9_20260907_ylinear_genus.md)
allows arbitrary A(p). Item20 extends the latter to variable B(p), retaining
every common-factor valuation. The index-two infinity branch replaces a
failed zero-valuation primitivity test. Item19's six-term supplier has
`Rstar=y³(Q4(p)+b y²)` and belongs to neither y-linear class; a distinct
degree-seven generic-curve calculation in item22 now excludes its nonzero
rational times while preserving the proved finite-row action.

Reassessment of every live lane: the source image now pays finite
compatibility but exposes earlier motion and rationality as separate debts;
the infinity lane uses the full affine group after its boundary quotient
lost decisive relations; collision quadrics still need an actual transported
deformation; the torsion connection still needs a fixed smooth polynomial;
and the V4 boundary surface still needs a coordinate pair outside the proven
collapsed/prescribed-coordinate classes. Shared procedure is not a theorem
identifying those five objects.


The former finite(2,7) scouts are now closed by items24-25, together with
the finite(2,9) stratum. Their numerical words have been replaced by full
rational path certificates and independently audited actual group maps.
The next geometric object was item28's two-ordinary-cusp sextic; item31
now pays its simultaneous access paths and proves its group is Z. Its
single smooth projection fold alone only gave three meridians. Item30's
abstract A6 survivor explains the additional global information required.
These curves and finite passports are not Keller maps.

Incoming `89aca08819` and `4441760e1a` were read before synchronization.
Their [continuing9 synthesis](continuing9_20260907_synthesis.md) confirms the
primitive-invariant and variable-coefficient boundaries used above; the new
[no-three-in-line repair theorem](continuing10_20260907_no3_repair.md) retains
original points in a capacity-two graph before bounding an adaptive repair.
Our meridian graph retains actual moved sheets before testing transitivity.
These are related proof operations, with different vertices and predicates;
neither graph supplies a mathematical map to the other problem. The new
fixed-moment corrections do not change any dependency in this checkpoint.

### September 8 recovery and next decisive operations

The account usage limit interrupted the agents after the September 7
02:37 UTC work. Missed heartbeat invocations are not counted as research
or verification. Work resumed September 8 at01:35 UTC by reading the
saved proof drafts, completing fresh independent replays, and auditing
the unfinished dependencies. The original September 8 20:40 UTC cutoff
remains in force.

Incoming `e291fcc945` was read before this checkpoint. Its
[source-linear global-carrier theorem](continuing10_20260907_dg_linear_carrier.md)
has an actual map into our W lane: the same two-chart intersection ring
and the same bracket. The recovered quadratic/cubic plane theorems then
strengthen its first-coordinate exclusion via the missing global inverse
coordinate. Its critical-free, boundary-separating function with a rational
conjugate remains the hostile: regularity on the source chart and generic
A1 geometry do not pay the two special-fibre pole conditions.
Its [joint D-moment octant exclusion](continuing10_20260907_nonpositive_circuits.md)
and [native wedge topology](continuing10_20260907_lrc_wedge_topology.md)
also restore simultaneous information before a quotient. That procedural
comparison supplies no map from D moments or runner labels to the present
curve or source-ring objects, and no JC consequence is inferred.

Next operations compare all five live lanes. The quartic common-root
and shared-root first-jet gates remain the entry constraints. Items42--43
close the constant-D finite6/infinity2 stratum with M a unit
at its finite point, at every finite location, including the genuine residue-free elliptic layer. Item47
closes the complete fixed-zero octuple polynomial-mate class; item50
pays every finite point using changed global carriers, the complete
three-dimensional genus-two primitive space, and original-source pole
repair. Item53 pays the infinity-point octuple with the actual surface
inversion's changed volume form, completing all single-point divisors.
Items54,56 and60 close every two-fourfold-root octic at all boundary
locations. The degree-four moving-residue obstruction closes the two-finite
case even for rational mates. Items62--63 now close every4,2,2 position, every6,2 case with its double
point finite, and every leading differential outside the complete table.
Item64 closes5,1,1,1 by its full later-coefficient five-zero argument,
with rational sharpness realized inside the same partition. Item66 pays
all exact degree-eight pencils and every global quadratic rational mate.
Item67 pays the all-degree genus bound and complete extremal sparse-curve
locus. Item70 pays its actual-surface extension to z=x^m, the full
extremal-pencil iff and sharp rational families. Item71 pays all-exponent
binomial exactness. Item69 closes all-finite4,3,1 with both complete
primitive systems, geometric integrality and full boundary checks.
Item72 closes fourfold infinity as well, completing all locations.
The all-finite5,3 and6,1,1 strata are under new coefficient and primitive
analysis. Item74 closes the complete M-shared constant-D finite6/infinity2
complement, with genuine rational submersions preserved as hostiles.
Together with items42–43, the entire constant-D polynomial class is closed.
The nonconstant-D sixfold class remains a separate entry. The other surviving leading strata retain their
lower-coefficient and global descent obligations. The degree
of F|D is distinct from a primitive's local degree3, which comes from
the order-two original differential; no degree-four primitive is assumed.
Item39 closes the entire two-ordinary-cusp (4,6) class. Items46--48 pay
the actual three-cusp moving paths, local pairs and full family transport;
item52 now closes the entire ordinary three-cusp class in every degree.
The finite involution and single-cycle results remain audited antecedents.
Items55 and58 pay the complete adjacent(5,3,3)/infinity7 good family
and its actual marked H4 quotient. Items57 and59 exclude single-cycle
and involutive meridians. The remaining target is genuinely mixed cycle
structure with arbitrary actual retained subsets. Item61 proves the sharp odd-braid retention bounds and mapping-degree
floor12. Item65 closes the first mixed(3)(2) cycle type by actual
cycle-intersection data. Item68 closes the unordered(3)(3) type by preserving local block mixing
and spending surrounding centralizers. Item73 closes(4)(2) as well.
All moved-support sizes at most six are eliminated; item75 proves
the stronger actual mapping-degree floor16 with all retention cases.
Other mixed cycle structures remain open. The unbounded cardinality relaxation still leaves actual
cycle order as the decisive sidecar. The exact pair(123),(345) still refutes
transfer of the ordinary half-support bound. The S4 hostile
prevents cyclicity or local cusp generation
shortcuts. For moving sources, every polynomial
Hamiltonian inside the cusp ideal is closed for nonzero rational time.
Different-invariant compositions and changed carriers retain their
invariant and regularity debts. Collision quadrics still need a
transported deformation. Item34 pays a fixed smooth polynomial and all
component principal parts for the torsion connection, while its map
from moving source responses remains open. No curve or finite-group
result supplies that map. The META-PATTERNS cards used remain Search
the statement before the method, Inventory retained power classes,
and Separate descent from regularity. No new card is promoted here.


### Superseded experiment retained for provenance

The [earlier symmetric4+2+2 critical-polynomial draft](../../04-computation/planar_jc48_sep08_four_two_two.py)
was overtaken by item62, which excludes the entire partition even
rationally at all locations. No companion proof or independent audit was
completed for that narrower route. Its77 exact algebra controls and
forbidden-root example remain a research draft, explicitly outside the
audited manifest and outside the proved dependency graph.


## September8 pre-promotion correction: complete pole-space membership

The all-finite5,3 scout rejected a genus-three candidate that repaired a
boundary principal part by adding P(1)/(u-1): it introduced poles at the
two ordinary affine points of that same vertical fibre. No such result
was promoted. The [mistakes entry](../../01-canon/MISTAKES.md) records
root's initial acceptance and the producer's correction. The replacement
uses only powers of u in its denominator and cancels the boundary jet
with a polynomial/global section combination; item76 and its independent
audit now prove full pole-space membership and the complete exclusion.


## September8 complete remaining DG quartic entry

Items77–79 close three full binary boundary types while retaining their
actual rational exceptions. The final audit of the two5+3 infinity
placements is underway; those placements have nonconstant-L rational
mates, so only the polynomial conclusion can be combined across all
locations. Once accepted, the complete leading differential table will
leave only binary7+1,6+2 and5+2+1. The constant-D6+2 class is already
closed by item74 and its M-unit companion; nonconstant D is separate.
The current portfolio has5+2+1 residue/primitive geometry, the full
remaining6+2 coefficient class, and an orthogonal H4 mixed5,2 group lane.
