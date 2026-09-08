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
| Anchor: infinity | The whole birational (4,6) class with one finite odd (2,m) cusp and at least two nodes is excluded; the new two-cusp sextic has actual affine group Z | Extend the finite certified relations under coefficient perturbation and identify the complete good two-cusp parameter strata; the local A6 passport alone loses these global relations |
| Wildcard: response connection | Full torsion is component-labelled principal parts; the canonical derivative raises every nonzero primary height | Map the actual source-normal response to this operator without changing its polynomial ring |
| Wildcard: actual boundary surface | Every hypothetical quartic pair has F=H²+L with global H in L2 and nonconstant L in L1; all coordinates linear in r in the boundary chart are excluded, with unrestricted mate | Test the actual relative differential on compact generic fibres; distinguish regular boundary branches from nonzero logarithmic residues and retain both components of a special fibre |

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

The [142-artifact manifest](planar_jc48_sep06_manifest.json) pins32 proof notes,
33 independent audits,34 sources,34 outputs and nine compressed certificates.
The34 programs report1444297 exact gates per complete normal or optimized run. Each frozen output agrees in normal and
optimized modes; the braid verifiers also check their pinned witnesses.
These are scoped proof controls, not a census of Keller maps. No external
priority claim is made.

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

Next operations compare all five live lanes. Item29 pays globality of the
quartic approximate root using the base-only Jacobian extension; F global
alone remains insufficient, as its explicit hostile shows. Retain global
H in L2, L in L1, and the unrestricted mate. On the compact surface,
the generic fibre has equation N²+M S³-c S4=0. Simple zeros of N on S
with M nonzero are a positive signal for regularity of the relative
differential; multiple zeros can create nonzero logarithmic residues,
so a claim of automatic regularity must not be extrapolated. This local
direction is under investigation, not part of the proved manifest.
For the two-cusp curve, the four certified finite relations are stable
under sufficiently small coefficient perturbations. Determine the complete
good parameter strata before promoting that openness to a whole-family
statement; exceptional repeated roots of U' and infinity strata remain
explicit obligations. For moving sources, all
polynomial lifts inside the cusp ideal are closed for rational time;
different-invariant compositions and changed global carriers need new
arguments. Collision quadrics still need an actual transported deformation;
the torsion connection still needs a fixed smooth polynomial and all
component principal parts. Neither new genus nor braid result pays those
missing maps. The META-PATTERNS cards used remain Search the statement
before the method, Inventory retained power classes, and Separate descent
from regularity. No new card is promoted from these related curve arguments.
