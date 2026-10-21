# Session23 extension: actual maps, flow generating functions, and information that cannot be discarded

The session began on 2026-10-21 at19:23 UTC and was extended at20:08 UTC at the user's request. This report records the extension's mathematical results, not an execution-status or session-END notice. The [first-phase synthesis](../session23_generating_function_synthesis_20261021/README.md) records the earlier denominator, marked-path, renewal and symmetry-slack results.

The strongest extension result is a new explicit carrier for the original operation: bounded integral flows with their selected Hamiltonian path, capacity, and labelled demand retained. The tetrahedral/Fano ideas also yielded a genuine duality, a rigorous geometric obstruction and replacement, and exact operation-relative generating functions. Classical ingredients are credited; a targeted repo search is not a literature-priority claim.

## Counting variables must stay typed

- n: tournament order; z in the SCC packets marks vertices.
- M: a uniform integer edge capacity; w marks this capacity in Ehrhart series.
- j: total arc use of an integer flow; q marks this use.
- k: required colors per vertex in an acyclic multicoloring; u marks k.
- r: either the half-order in n=2r+1, or a substitution depth when explicitly stated. These are separate uses, not a numerical identification.

A rational function in one of these variables need not be the same enumerator, or even have the same rationality behavior, in another.

## 1. An exact map on the original fixed-path operation

The [capacity-flow packet](../capacity_flow_metagraph_transport_20261021/README.md) fixes the incidence convention Bf=inflow minus outflow. If a set D of arcs is reversed, complement their flow coordinates inside a retained capacity M:

    f'_e=M-f_e on D;   f'_e=f_e elsewhere.

Then exactly

    B' f' = Bf - M B_D 1.

This is an integral affine bijection, involutive when both the orientation and demand are updated. It is not a claim that arbitrary partial reversal preserves circulations.

For a regular tournament T and an actual directed Hamiltonian path P, take D to be all nonpath arcs. The path itself stays unchanged, and

    bounded circulations on T
        <-> bounded flows on its nonpath-complement partner
            with demand M(e_last-e_first).

The partner has outdegree r+1 at the first endpoint, r internally, and r-1 at the last endpoint. It is not regular. Holding the demand at zero would be a genuine type error: at order5,capacity1 the source has24 circulations and the partner has24 flows of the correct endpoint demand, but only14 circulations.

The score/demand invariant is M s(T)+b, where s is the labelled outdegree vector. The score-complement identity itself is inherited from THM-646; the flow map is the integral lift. In centered coordinates j_e=f_e-M/2, reversal simply sends j_e to -j_e on D. Centered divergence and every fixed-edge quadratic energy sum r_e j_e² are preserved. This is orientation-covariant current bookkeeping, not a new Ohm's-law or voltage realization theorem. The uncentered endpoint demand changes with the chosen baseline.

### Universal capacity counts, and their exact limitation

Any two tournaments with the same labelled score vector have integrally affine-isomorphic bounded-circulation polytopes. For regular tournaments, the entire unweighted capacity count therefore depends only on n, not the isomorphism class. At capacity1 it is exactly the number of labelled regular tournaments:24 at order5 and2640 at order7.

This is a strong information-loss theorem. The cyclic regular7 tournament has H=175 and Paley7 has H=189, yet their capacity counts agree for EVERY M. Their capacity-one total-use coefficient at q⁴ differs,28 versus21. Scalar capacity counts cannot recover Hamilton-path counts, much less complete metagraph wiring. The edge coordinates and their changed weights carry real information.

For odd n>=3, let d=(n-1)(n-2)/2. The full capacity polytope is integral of dimension d, so its count N_n(M) is an Ehrhart polynomial of exact degree d. Its series has reduced denominator (1-w)^(d+1), and its numerator is palindromic of exact degree d-1. The proof uses all-ones translation of interior points at dilation M to closed points at dilation M-2. The edgeless order1 case is explicitly separate.

At order5,

    N_5(M)=(M+1)²[3(M+1)⁴+5(M+1)²+4]/12,

    sum N_5(M)w^M
       =(1+17w+72w²+72w³+17w⁴+w⁵)/(1-w)^7.

The degree bound is proved first; exact determining values and an independent count establish the coefficients. This is not interpolation from an assumed recurrence.

### Keeping the SAME Hamiltonian path: a different numerator law

Impose f_e=0 on P. At capacity1 these flows biject to regular tournaments containing that same ordered path. The count equals the sum of H(T)/|Aut(T)| over regular isomorphism classes:3 at order5 and91 at order7.

For every odd n>=5, T minus P is strongly connected. A cut with s vertices has s(n-s)/2 outgoing arcs in a regular tournament, and the path can remove at most min(s,n-s), leaving a positive number. Hence the restricted polytope has dimension

    delta=(n-1)(n-4)/2.

Its Ehrhart series has reduced denominator (1-w)^(delta+1). Its exact codegree is3: dilation2 would require the unbalanced all-ones nonpath vector; at dilation3, add the indicator of a first-to-last path in T minus P to repair that imbalance. Thus the numerator has exact degree delta-2.

This numerator is NOT generally palindromic. At order7 it has degree7, constant coefficient1 and leading coefficient34. The leading coefficient counts the interior points at dilation3, equivalently the path-containing tournaments with score vector (r-1,r,...,r,r+1). The coefficient of w is91-10=81; other coefficients are not claimed computed here.

At order5 the restricted count is binomial(M+2,2), with series1/(1-w)^3. The cyclic base has a triangle and pentagon sharing a path; their nonnegative multiplicities a,b satisfy a+b<=M. Its capacity/use series is

    1/[(1-w)(1-wq³)(1-wq⁵)].

The three actual Hamilton-path presentations have capacity-one use polynomials1+q³+q⁵,1+q⁴+q⁵,1+q³+q⁴. They agree when q=1 but not before that specialization. This is a small exact example of why the selected path cannot be discarded.

## 2. All-order regular-flow reciprocity and symmetry rarity

The [reciprocity packet](../regular_tournament_flow_reciprocity_20261021/README.md) proves a different theorem, now with q marking total arc use. Let T be regular of order n=2r+1, with m=nr arcs and circulation dimension d=r(2r-1). For every ordinary automorphism or anti-automorphism a, its fixed-flow series satisfies

    R_a(q^-1)=(-1)^d q^m R_a(q).

Classical lattice-cone reciprocity supplies the interior relation. Regularity supplies the all-ones translation. The tournament-specific step proves that every fixed space has the same dimension parity: ordinary automorphisms have odd order, while an anti-involution has fixed circulation dimension r², and the general anti-map reduces to an involution and a commuting odd-order part. Therefore ordinary and anti-coset Burnside averages share the same sign and shift.

A matching denominator-factor parity forces a palindromic numerator of a predictable degree. This does not permit an arbitrary denominator presentation to be called minimal. Strong nonregular T4 is a hostile: its raw and anti-fixed flow series have the same shift10 but opposite reciprocity signs, so their average has no common sign.

For each fixed self-converse regular T of order2r+1>=5, the fraction of ordinary flow-orbits with total use AT MOST M that are converse-fixed obeys

    C_T M^[-r(r-1)] (1+O(1/M)),   C_T>0.

The same exponent holds for capacity boxes, generally with another constant. The proof uses fixed-space dimensions and lattice-volume asymptotics. It does not need pointwise coefficient aperiodicity. At order7 the cumulative exponent is -6. For R5 the cumulative total-use constant is3150/17, whereas the pointwise coefficient-ratio constant below is2100/17. These two statements must not be conflated.

Three explicit Paley7 fixed-cone series accompany the theorem:

    order7 ordinary: 1/(1-q⁷)³;
    order3 ordinary: (1+q³+q⁶+q⁹+q¹²)/[(1-q³)(1-q⁶)²(1-q⁹)²];
    order6 anti: (1+q¹⁵)/[(1-q⁶)(1-q¹²)(1-q¹⁸)].

All have sign minus and shift21. These are fixed-action corrections, not the still-uncomputed complete raw Paley7 numerator.

## 3. Explicit R5 formulas, with certified geometry rather than a fitted pattern

The [pentagon packet](../pentagon_flow_generating_functions_20261021/README.md) derives the complete raw circulation series R(q) on cyclic R5:

    R(q)=[1+2q³+3q⁴+q⁵+3q⁶+q⁷+3q⁸+2q⁹+q¹²]
           /[(1-q³)³(1-q⁴)²(1-q⁵)].

The ordinary flow-orbit series and converse-fixed ordinary-orbit series are

    f(q)=[R(q)+4/(1-q⁵)²]/5,

    k(q)=[1+q⁴+q⁵+q⁶+q⁷+q⁸+q¹²]
           /[(1-q³)(1-q⁵)(1-q⁶)(1-q⁸)].

The raw cone is an integral lift of the fractional stable-set polytope of a pentagon. Its half-integral vertex is essential. An exact25-cone unimodular triangulation and half-open boundary rule prove the numerator. Independent pairwise separation and facet certificates verify coverage and nonoverlap. A separate five-case lattice decomposition derives k without using the triangulation.

For R,f,k,(f+k)/2, the minimal constant-coefficient recurrence orders are22,26,22,32, and the minimal quasiperiods are60,60,120,120. In particular averaging can introduce a larger reduced denominator than either initially displayed ingredient.

The coefficients satisfy

    [q^j]R ~ 17 j⁵/259200,
    [q^j]f ~ 17 j⁵/1296000,
    [q^j]k ~ 7 j³/4320,

so [q^j]k/[q^j]f ~ (2100/17)j^-2. These are proved pointwise for this explicit R5 series, unlike the more general cumulative theorem.

## 4. Rational flow decorations require three different half-word carriers

The [flow-decorated SCC packet](../flow_decorated_tournament_series_20261021/README.md) begins with the unique strong order4 tournament. Its nonnegative circulation semigroup is free on two triangles and one four-cycle, with weights3,3,4. The anti-map swaps the two triangle coordinates. Thus its raw and converse-fixed series are respectively

    1/[(1-q³)²(1-q⁴)],   1/[(1-q⁶)(1-q⁴)].

Every nonnegative circulation vanishes on arcs between distinct strong components. Flow decorations therefore factor over the actual ordered SCC word. Fixed-cone rationality and Burnside, combined with the inherited finite strong-atom/center bounds, prove rational decorated series at every fixed UNDERLYING h and every fixed UNDERLYING fundamental black degree, including degree0. With the R5 formulas above, the degree0 series becomes entirely explicit.

The crucial distinction is between independent opposite-half decorations and a genuinely mirrored decorated half. Their neutral denominators are

    1-z²-z⁶/(1-q³)²,
    1-z²-z⁶/(1-q⁶),

respectively. A self-converse underlying tournament may carry a non-self-converse flow. Decorating C3 followed by C3 with unequal cycle multiplicities is the minimal hostile. Substituting one denominator for the other loses actual objects.

No formula here redefines h using the smaller stabilizer of a decorated object. Nor does a circulation count Euler tours: support may be disconnected, and traversal order is forgotten.

## 5. The tetrahedron produces two genuinely dual Fano models

The [duality bridge](../tetrahedral_fano_duality_20261021/README.md) connects two previously separate constructions exactly.

The edge model has six tetrahedral edges and a formal seventh point. Its lines are four face triples and three opposite-edge-pair triples through the formal point. The cycle model has the seven nonzero binary cycle-space elements: four face triangles and three Hamiltonian quadrilaterals, with XOR defining its lines.

The canonical map sends a face cycle to its face line and a quadrilateral cycle to the line through its two OMITTED opposite edges and the formal point. Cycle-plane lines become edge-plane concurrence pencils. This is S4-equivariant projective duality, not an equivariant point bijection. Such a point bijection is impossible: the natural point-action orbit sizes are1+6 in one model and4+3 in the other.

This supplies both an actual map and a sharp obstruction to a too-strong identification. It does not turn the seven metric feature points of the proposed sphere construction into either model automatically.

### The geometry claim itself and its positive replacement

The [geometry packet](../tetrahedral_inversion_spheres_20261021/README.md) disproves the universal seven-label sphere claim, even with all constructions finite. A rational unit-sphere tetrahedron has B,H′,F distinct and collinear, so no sphere can contain them. This counterexample is independent of the ambiguity in which original-sphere intersection is called A.

A complete exceptional-locus theorem is proved for the second-intersection interpretation. Its valid solutions have D=E and B=F, hence exactly five distinct points; those points are noncoplanar and therefore not concyclic. Degenerate centroid/tangent-plane cases and the nondegenerate rank-loss endpoint are explicitly handled.

There is also a positive generic six-point theorem. Normalize the original sphere to center0 and radius1, write v=V0, s=V1+V2+V3 and w=O′. When v,s,w are independent, B,C,D,E,F,4C lie on an explicit sphere invariant under radius2 inversion. The inverse pairs are(D,E),(B,F),(C,4C). The omitted A and H′ lie strictly outside this sphere. This is a separate construction, not a silent alteration of the original assertion.

## 6. Fano coding gives exact new move graphs, with a retained return obstruction

The [Fano move packet](../fano_tournament_move_quotients_20261021/README.md) uses tetrahedral edge incidence rather than merely seven labels.

Coherent boundary orientation of the four tetrahedral faces forces all eight cyclic completions to have H=171, not Paley's189. Reversing one face changes the relevant parity and permits exactly two Paley completions. The old16/112 Fano-orientation census and Hamming-code locus are inherited and explicitly credited.

Flipping one whole cyclic Fano line gives a Q7 graph, isometrically embedded in the full labelled regular7 triangle-reversal graph. Its syndrome quotient is K8. A triangle returning to the same syndrome lifts to a nonzero weight-three codeword, so it does not return to the original presentation. Adding parity restores bipartiteness but still does not retain the full codeword. This is a concrete return-state obstruction, not an analogy.

The full regular7 triangle-move quotient, in H189,H171,H175 order, is

    [0,14,0; 2,9,3; 0,7,7],

with class sizes240,1680,720 and eigenvalues14,1+sqrt15,1-sqrt15. The H-weighted move-word series starting at Paley is

    (189-630z-1470z²)/[(1-14z)(1-2z-14z²)].

Paley hitting means are7 in the restricted Fano cube, versus10 and12 from the other two classes in the full move graph. These are different operations with different transition laws. Triangle reversal is NOT the original fixed-path nonpath-complement operation; the capacity packet is the direct bridge to that latter operation.

Every anti-automorphism of a seven-vertex tournament is an odd permutation, whereas the Fano collineation group is even. Hence no anti-action preserves one Fano plane. Paley has exactly two contained Fano designs, and its anti-actions exchange them. A chosen design is real additional data.

## 7. Three kinds of color, and a rationality boundary

The [signed-rigidity packet](../signed_tournament_palette_rigidity_20261021/README.md) reads the actual uploaded Gluck paper. It is a finite-permutation-group paper, not the source of the stated tetrahedron geometry. Its classical odd-order regular-subset theorem implies ordinary binary rigidity of tournaments and faithful residual marked-path actions. A third color kills the remaining possible anti-involution for n>=2. Kernel exceptions for empty/singleton carriers are retained.

For the Paley ordinary group, signed Paley group and full Fano group, the numbers of regular named c-colorings are respectively

    c⁷-7c³+6c,
    c⁷-7c⁴-7c³+7c²+6c,
    c⁷-21c⁵+28c⁴+56c³-112c²+48c.

Their least distinguishing palette sizes are2,3,4. A single equivariantly chosen subset cannot break a nontrivial automorphism group: the equivariance requirement would make it invariant. The full marked family or its orbit must be retained. Existing ordinary two-color results are credited, not billed as new.

The [acyclic-palette packet](../paley_multicoloring_and_carry_20261021/README.md) concerns a different predicate. Paley7 has exact k-fold acyclic coloring number ceil(7k/3), with explicit small certificates. Its palette-budget generating function is

    (3u+2u²+2u³)/[(1-u)(1-u³)],

of minimal recurrence order4. Substitution obeys chi_d(P7[U])=ceil(7 chi_d(U)/3), so chi_d(P7[P7])=7 rather than9.

Iterating substitution gives a_0=1, a_(t+1)=ceil(7a_t/3), starting1,3,7,17,40,94,... . Its depth generating function is NOT rational. More generally the ceiling recurrence for any reduced rational multiplier A/B>1 with B>1 and positive integer initial value is not C-finite. A hypothetical recurrence would make its bounded carry eventually periodic; a periodic block then forces an integer with arbitrarily large powers of B in its denominator, a contradiction.

A finite carry alphabet is therefore not an autonomous finite-state model of the integer process. This does not prove a Collatz theorem. Distinguishing colors, acyclic colors and proper unit-distance colors remain three separate conditions; none of these numerical palette values is a Hadwiger–Nelson bound.

## 8. Current repo work, remaining obligations, and precise next targets

During publication, incoming commit1e25cb70aec303cb6e8515bad610c9e99ddb5f7f updated the base-p Collatz/coalescence audit. It was fetched and preserved, not overwritten. Its narrowed Haar/almost-sure scope and known sign-conjugation barrier reinforce an existing guardrail: measure-theoretic or quotient statements do not automatically cover every positive integer. That incoming theorem was not re-proved in this session, and the rounded-carry example is not identified with its maps.

What remains open in this extension:

1. Complete original metagraph wiring, all bridges, and general realization criteria are not classified. Scalar counts alone are demonstrably insufficient.
2. The uniform153N source-budget obligation and long-horizon R′ bound are not advanced or certified by these flow analogies. Native mod15/mod35 guards and actual clocks have not been replaced.
3. No general A000568 speedup has been proved or benchmarked. Capacity1 counts labelled score fibres, not unlabelled tournament classes.
4. No Lonely Runner, Hadwiger–Nelson, Collatz, or universal Hamiltonian-path-spectrum conjecture is settled here.
5. The new order7 fixed-path Ehrhart numerator has known degree7, first coefficients1,81 and last coefficient34; its remaining coefficients are a concrete finite next target.
6. The raw Paley7 flow numerator and anti-involution fixed-cone numerator remain uncomputed. Known dimensions and reciprocity sharply constrain them; the smaller fixed-action formulas are already explicit.
7. The higher-order triangle-move quotient must be tested for equitability before any order7 recurrence is generalized. Same H alone is not a valid state partition without proof.

The reusable method is to specify the objects, the exact operation, the transported observables, and what is lost. This session supplies successful maps, as well as counterexamples that prevent an attractive but invalid transfer.
