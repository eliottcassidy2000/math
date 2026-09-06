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
| Anchor: moving source | Weight22 transport has a Hamiltonian realization; exact source-preserving generators all fail local nilpotence | Find a finite or formal carrier that pays later rows while retaining actual depth; additive polynomial flows are excluded in the proved class |
| Niche: collision | All scalar two-form periods miss exactly the quadrics through tangent directions in a three-dimensional target | Carry the quadratic data through an actual earlier source deformation and preserve descent |
| Niche: infinity | An explicit intrinsic (4,6) family is excluded as sole support except at three cusp parameters | Test the actual cusp-plus-five-node passport, retaining local monodromy and affine fibres |
| Wildcard: response connection | Full torsion is component-labelled principal parts; the canonical derivative raises every nonzero primary height | Map the actual source-normal response to this operator without changing its polynomial ring |
| Bridge: finite versus global | Each quotient now has an explicit omitted coordinate | Test the first later row or global extension where that coordinate matters |

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
   except at the three roots of `128lambda^3-288lambda-283`. These have
   one cusp and five nodes and remain OPEN; their necessary passport is
   `n_cusp+sum_overlaps=1` and `d<=2a`. Producer214 gates and two independent
   analytic reads pass. This does not classify all intrinsic (4,6) curves.

The [26-artifact manifest](planar_jc48_sep06_manifest.json) pins six proof
notes, six independent audits, seven sources and seven outputs. The seven
programs report **5,475** exact gates in each combined normal/optimized
run. Their outputs and the immutable input certificate agree byte for byte.
Computations support the stated proofs; they are not a census of Keller maps.

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
excluded. Test source-preserving formal generators against the next row
and projected depth, or change the carrier while retaining its actual
source equation. A successful finite table is not polynomial termination.

At infinity the exact family has only three cusp parameters left. Keep
the five individual node overlaps and the actual cusp fibre visible;
meridian fixed points alone do not identify actual retained sheets. The
next cheap probe is a faithful cusp monodromy passport consistent with
the global Euler equation. An abstract passport needs an actual affine
source and polynomial completion before it becomes a Keller candidate.
With several nonnodal points the one-point positivity argument changes;
retain their number and all actual special fibres before generalizing.

The heartbeat now runs this board every30 minutes through the stated
48-hour cutoff. It stays quiet unless there is substantive progress, a
correction, completion, failure or required user action. At the cutoff,
push the final coherent checkpoint and pause the heartbeat.
