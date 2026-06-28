# LRC14 Shadow-Charge Conservation Atlas

The useful shift today was to stop asking which frame is "the" frame.  The
better question is which proof charge a frame carries, and what happens to that
charge when the frame is projected.

Post-fetch, the incoming index-theorem and contact-holonomy work sharpened
this.  The atlas was renumbered to HYP-3400 after upstream claimed HYP-3260
for the unit-nullspace/blind-residue layer and HYP-3310 for the C6
residue-magnitude factorization inside the index/Chebyshev/contact-holonomy,
floor-census, tight-locus-manifold, contact-graph, and observability/Morse
sequence.
That was not just bookkeeping.  Analytic Cech/Euler index, topological
Borsuk-Ulam degree, and Gauss-sum index became the top reservoir,
`index_theorem_degree_charge`; mac-mini's finite tight-locus plus margin and
the HYP-3310/HYP-3300/HYP-3266/HYP-3265/HYP-3260/HYP-3254/HYP-3256/HYP-3258/HYP-3259
C6 residue/magnitude, observability/Morse, proof-obligation, contact-graph,
unit-nullspace, and manifold split became
`uniform_margin_floor_charge`; and the HYP-3253 shell-lag/contact result became
`contact_holonomy_curvature_charge`.  The debt is now sharper: nonzero degree
must still be connected to an actual lonely point, and any lag/residue quotient
with curvature must be killed by exact `zeta_7` holonomy, lifted to endpoint
cells, or named.

This makes the topology language less ornamental.  HYP-3242's Cech/Euler hole,
HYP-3243's topes and graph carriers, HYP-3244's lift/descent span, HYP-3245's
equioscillation/autocorrelation pair, HYP-3238's even-positive/odd-negative
duality, and HYP-3236's conductance graph are not rival stories.  They are
different boundary measurements on one packet.

The synthesis I like most:

```text
Cech boundary
Green current
autocorrelation lag transport
normal-fan slack
root-motion collision
contact-holonomy curvature
state-lift obstruction
unit-contact graph case split
observability/Morse descent
```

should be read as a conservation law around a finite obstruction.  If a non-AP
trap moves autocorrelation mass outward, it should also have a Green bottleneck,
a normal-fan slack coordinate, a topological owner change, a root-motion debt,
contact-holonomy curvature discharge, or a state-lift label.  The proof search
should try to make that "or" into a finite typed disjunction with no unnamed
remainder.

The information-theory angle also became sharper.  HYP-3201 already showed
that failed commutativity, associativity, idempotence, and distributivity are
not philosophical failures; they are conditional-entropy residues on quotient
fibers.  The same rule applies to LRC:

```text
positive graphing forgets negative covariance
even moments forget odd sign
Euler characteristic can cancel Betti data
root radius forgets endpoint owner
autocorrelation forgets witness address
lag/residue quotient forgets contact holonomy
raw A000568 count forgets tiling fiber
```

Each of those failures is useful only if it emits a sidecar.  That is the
"no naked quotient" lemma in embryo.

The strongest proof-facing packet now looks like this:

```text
index theorem / degree / Gauss-sum charge
+ finite tight-locus / uniform-margin floor
+ zeta7 contact-holonomy curvature sidecar
+ Phi_14/Phi_{14d} witness address
+ tiling/half-tiling descent certificate
+ Cech open-hole status
+ D7 odd sign payload
+ Lee-Yang/Joukowski/HB root sidecar
+ Green/Dirichlet current sidecar
+ discrepancy/Hensel bulk floor
+ finite-address/state-lift exit
```

This is bulky, but that is exactly the point.  The scalar routes failed because
the proof object is not scalar.  The next compression should be a theorem about
why this packet can be reduced, not an unpriced wish that it already has been.

Risk: the charge language could become just another metaphor.  The antidote is
the executable ledger: for each reservoir, name preserved fields, destroyed
fields, transfers, and proof exits.  The next version should stop being a
synthetic atlas and become a trap-by-trap table.

Concrete next computation:

1. For each HYP-3202 non-AP trap, compute a row with autocorrelation transport,
   Green resistance excess, Toeplitz slack, normal-fan first failed coordinate,
   D7 sign payload, root-motion class, and contact-holonomy curvature status.
2. Add a `charge_balance` column: which finite proof exit does the trap route
   to first?
3. Mark any row whose charge is not discharged as `named_residual_debt`.
4. In parallel, test index equality as descriptor plus the S-dependent floor
   as proof.
5. If the residual set is empty in the bounded bank, try to express the table
   as a symbolic chamber theorem.

The cleanest possible proof route from here is not "find a magic scalar."  It is
closer to:

```text
build the full packet
prove charge conservation under legal maps
show every primitive packet exits by index nonvanishing, witness, bulk, trap discharge, holonomy discharge, descent, or H=7 contradiction
prove no named residual debt remains
```

That is still a long road.  But the road has better signs now.
