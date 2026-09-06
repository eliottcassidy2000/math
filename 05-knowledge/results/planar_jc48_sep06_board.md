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
| Niche: infinity | Ordinary-cusp-plus-nodes supports are excluded; a single (2,m) cusp with at least two nodes forces d-1 to divide m | Investigate the actual degree-six source for the explicit (2,5)-cusp/five-node target |
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

The [34-artifact manifest](planar_jc48_sep06_manifest.json) pins eight proof
notes, eight independent audits, nine sources and nine outputs. The nine
programs report **26,362** exact gates in each combined normal/optimized
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

At infinity the exact family and ordinary-cusp-plus-nodes class are closed.
The higher-cusp spectrum is now proved; the concrete next target below
forces d=6, a=delta=3, one actual cusp point and five omitted node points.
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

### A concrete higher-cusp target

**VERIFIED target-curve control / Keller realization OPEN.** The curve
`(U,V)=(t^4+t^2,t^6+t^5+t^2)` has one (2,5) cusp and **five** ordinary
nodes. The initial guess of four nodes missed the pair-sum root p=1,
whose parameters satisfy t^2-t+1=0 and have image (-1,1). The geometry
agent derived the complete pair quotient; root independently checked the
ten identities below. This is a concrete target for the higher-cusp
passport, not a constructed Keller map. These scouting checks are separate
from the frozen-program gate total above.

The critical gcd is t and the triangular target shear is
`V-U+U^2=t^5+3t^6+t^8`, giving the sole (2,5) cusp. The nonzero pair sums
have five simple roots and two distinct parameters each. The displayed
Groebner ideal excludes off-diagonal tangencies. A possible triple fibre
would force both `A^2+A-1=0` and `(A-2)(A^2+1)=0`, where A=a+2 and a is
the common target's first coordinate. Their Bezout identity below excludes
it. Thus the five pairs have distinct node images.
Monicity of U gives a finite normalization, and the zero-dimensional
collision scheme gives birationality and intrinsic pole pair (4,6).

Reproduce from any checkout with SymPy:

```sh
python3 - <<'PY'
import sympy as S
s,t,p,A=S.symbols('s t p A')
U=t**4+t**2; V=t**6+t**5+t**2
N=S.cancel((U.subs(t,s)-U)/(s-t)); M=S.cancel((V.subs(t,s)-V)/(s-t))
D=S.diff(U,t).subs(t,s)*S.diff(V,t)-S.diff(U,t)*S.diff(V,t).subs(t,s)
H=-(p-1)*(p**4+2*p**3+4*p**2+8*p+1)
checks=[S.gcd(S.diff(U,t),S.diff(V,t))==t,
 S.expand(V-U+U**2)==t**5+3*t**6+t**8, S.gcd(U,V)==t**2,
 S.expand(N-(s+t)*(s*s+t*t+1))==0, S.expand(M.subs(t,-s))==s**4,
 S.rem(S.expand(4*M.subs(t,p-s)-H),s*s-p*s+(p*p+1)/2,s)==0,
 S.discriminant(H,p)==-11776000, S.resultant(H,p*p+2,p)==123,
 S.groebner([N,M,D],s,t,domain=S.QQ)==S.groebner([s+t,t**4],s,t,domain=S.QQ),
 S.expand((A*A-A-1)*(A*A+A-1)-(A+2)*(A-2)*(A*A+1))==5]
if not all(checks): raise RuntimeError('higher-cusp target control failed')
print('PASS',len(checks),'exact target-curve identities; no Keller realization')
PY
```

The heartbeat now runs this board every30 minutes through the stated
48-hour cutoff. It stays quiet unless there is substantive progress, a
correction, completion, failure or required user action. At the cutoff,
push the final coherent checkpoint and pause the heartbeat.
