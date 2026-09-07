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
| Anchor: moving source | Incoming work classifies the valuation-eleven weight22 response and integrates the universal carrier in the actual completion; polynomial carrier automorphisms are scalings | Test whether the required finite response lies in the completed carrier image, retaining all coordinate corrections; seek genuinely different carriers for polynomial operations |
| Niche: collision | All scalar two-form periods miss exactly the quadrics through tangent directions in a three-dimensional target | Carry the quadratic data through an actual earlier source deformation and preserve descent |
| Anchor: infinity | Nori excludes the five-node target; marked infinity and exact braids exclude the four-node equality target | Test the full affine complement of the new three-node sextic, whose actual infinity group admits A6 |
| Wildcard: response connection | Full torsion is component-labelled principal parts; the canonical derivative raises every nonzero primary height | Map the actual source-normal response to this operator without changing its polynomial ring |
| Wildcard: actual boundary surface | The classical V4 surface realizes the bare constraints and has a finite quartic map; every shear-orbit carrier still collapses its boundary | Both coordinates must escape the collapsed algebra; a first coordinate in C[t] forces interior ramification and even boundary index |

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
    is a precise hostile to extending item13. Full affine-complement
    realization and a polynomial Keller source remain **OPEN**.

The [71-artifact manifest](planar_jc48_sep06_manifest.json) pins seventeen
proof notes, seventeen independent audits, eighteen sources and eighteen
outputs, plus the compressed rational-path certificate. The eighteen
programs report **221,396** exact gates per complete normal or optimized
replay. Each output is byte-identical in both modes; the braid verifier
also checks the pinned witness. These are scoped proof controls, not a
census of Keller maps. No external priority claim is made.

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
whole carrier. A productive repair must change the carrier, preserve only
one genuine fixed source, or remain explicitly formal. A successful finite
table is not polynomial termination.

At infinity the exact family and ordinary-cusp-plus-nodes class are closed.
The higher-cusp spectrum is proved. Both the five-node target and its four-node equality successor are now
closed. The new three-node target below forces d=6, a=delta=3, one actual
cusp point and three omitted node points.
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

A concrete new source consumer is the image test
`R -> pi15{-u/2+H_xi,p²Delta R}` against the required seven-equation
response, retaining all induced coordinate corrections and original
bracket/depth/prefix equations. A positive answer would connect the
finite response to a completed operation; polynomial specialization and
later equations would remain unpaid. Our torsion connection requires an
actual fixed smooth polynomial and labelled principal parts. Neither
Delta, the background xi, nor the square-zero lift matrix is automatically
that connection.

Incoming `7c557fb04c` and its [continuing8 synthesis](continuing8_20260906_synthesis.md)
were also read. Its exact source-chart transfer boundaries preserve the
operator, characteristic, lattice and zero tail separately; they give no
characteristic-zero Keller obstruction. This agrees with the distinct
regularity and carrier hypotheses above. No frozen artifact was changed
by these incoming scope comparisons.

### Closed targets and the current concrete survivor

The curve `(t^4+t^2,t^6+t^5+t^2)` is **EXCLUDED** as whole Keller
nonproperness support by Nori: its resolved square12 exceeds twice its
five-node count10. The curve `(t^4+t^3+t^2,2t^6+3t^5+2t^3+2t^2)` is
also **EXCLUDED**, despite Nori equality: its marked infinity group
cannot supply the actual six-sheet representation. Exact global braid
constraints independently corroborate this second closure.

The [current actual target](planar_jc48_sep06_next_infinity.md) is

```text
U=t^4+t^3+t^2,
V=16t^6+24t^5-19t^3-19t^2.
```

Its **audited geometry** is one affine(2,5) cusp, three ordinary nodes,
and infinity(2,11). The resolution cost is18 and D²=4, giving Nori
margin-2. Its actual marked infinity tree has an explicit surjective A6
representation with a single-three-cycle meridian. This boundary witness
is **not** a representation of the full affine complement until it
satisfies the additional actual relations. The next cheap decisive test
is a rationally certified global braid constraint, adapting the validated
method rather than treating numerical tracking as proof.

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
