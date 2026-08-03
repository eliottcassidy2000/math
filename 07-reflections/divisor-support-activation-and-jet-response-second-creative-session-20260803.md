# Divisor support, activation graphs, and jet response

**Status:** CURRENT SYNTHESIS over
[THM-3320](../01-canon/theorems/THM-3320-projected-k3-z216-fourth-ruler-prefix-and-affine-multicover-closure.md),
[THM-3325](../01-canon/theorems/THM-3325-composite-even-flat-boundary-blindness-and-interior-collapse.md),
and
[THM-3326](../01-canon/theorems/THM-3326-linear-in-z-unit-response-trichotomy-and-jet-torsion.md).
The first theorem remains projected LRC bookkeeping; the latter two are new
uniform mechanisms.  None closes LRC(14) or JC(2).

## Portfolio outcome

### Anchor: the next `z1=216` work prefix

The fourth deterministic prefix closes the singleton row `141` and the full
three-row `L=30576` family `133,219,359`:

```text
230 states = 172 crude + 58 status + 0 residual,
ledger/wall/families: 373157/353/31 -> 373153/349/29.     (1)
```

The important structural correction is not the count.  Two certificates that
need support three in the integral zero-constant template have exact
half-affine support-two duals.  Circuit arity therefore belongs to the chosen
dual coordinate system, not just to the packet.  Row `138` remains a genuine
support-three control after the affine coordinate is restored.

### Niche: a composite boundary false-positive flat collapses uniformly

For every odd prime `p`, the normalized even-coordinate flat modulo `2p` has

```text
(2p)^(p-1)-(2p-1)^(p-1)                                 (2)
```

right-boundary blockers, but zero is its only full-cell blocker.  At `p=7`,
this is `2,702,726` nonscalar boundary false positives removed at once.

The uniform proof came from replacing the six-row `n=14` table by its native
operation.  A half-turn at a chosen protector cell activates another
half-turn on odd shifts only along

```text
j -> j/2, 3j/2.                                         (3)
```

Every legal arrow lowers `nu_2(j)`.  Choosing an active minimum destroys all
other half-turn pressure on the odd shifts, and ordinary residues cannot
cover the remaining parity class.

### Wildcard: the torsion ladder is a marked formal-coordinate object

For `P=f(x)+lambda*x^d z` at a simple zero of `f-b`, the annihilator lengths
from THM-3318 survive unchanged:

```text
Ann(theta)=((P-b)^(d-1)),       Ann(mu)=((P-b)^d).       (4)
```

The scalar bridge does not survive.  It becomes

```text
(P-b)mu=-F(P-b)theta,
F(T)=(d-1)T psi'(T)/psi(T) mod T^(d-1),                 (5)
```

where `psi=(f-b)^(-1)` compositionally.  The pair `(f'(0),F)` reconstructs
the full truncated source jet.  Abstractly the module is unchanged; with its
canonical generators it remembers exactly what the abstract isomorphism
forgets.

Deforming the vertical coefficient instead produces a sharper bifurcation.
For every gradient-unimodular `P=f+gz`, constant `g` gives an exact unit
response, a single repeated geometric root gives `(4)`, and two or more roots
make `K[P]theta` free.  The last branch is the algebraic de Rham residue of
`dx/g`: no nonzero polynomial in `P` can remove its logarithmic part.

## Live concept board

| Concept | Native representation | Preserved invariant | Operation that revealed it | Information lost by the coarse view |
|---|---|---|---|---|
| projected `z216` packet | common-status table | exact necessary-screen emptiness | add an affine dual constant | template support is not intrinsic dual support |
| composite boundary flat | endpoint-owner hypergraph | fixed-boundary blocking | pass to interior protector cells | boundary forgets bin parity |
| half-turn pressure | directed activation graph | odd-shift danger | choose minimum `nu_2` | gcd cost alone sees only the number `p` |
| one-root Hamiltonian block | localized cokernel | annihilator lengths | mark canonical generators | abstract cyclic type forgets the source jet |
| multi-root Hamiltonian block | rational differential `dx/g` | residue vector | split geometric root support | pole order alone forgets logarithmic residues |

The board suggests one disciplined meta-move: before compressing an
obstruction to a scalar cost, record both its **support components** and the
operation by which one component activates another.  This is a research
heuristic, not an identification of LRC residues with Jacobian divisors.

## Three precise cross-connections

### 1. Support count comes before multiplicity

**Source:** the roots of `g` in THM-3326.

**Target:** endpoint-owner blocks in THM-3325.

**Map:** only the abstract incidence pattern is transferred: a support
component is a geometric root on the Jacobian side and a CRT owner block on
the LRC side.

**Preserved predicate:** whether a coarse local order can be killed by a
global response.

**Destroyed information:** algebraic residue values on one side and physical
time, width, and endpoint origin on the other.

**Needed sidecars:** partial-fraction residues for `dx/g`; bin parity and the
activation graph for micro-staircases.

**Cheapest decisive tests:** two repeated roots for `g`; one zero even
coordinate at the `a=p` boundary.

The common lesson is exact but the objects are not: one component permits a
finite killing ladder, whereas multiple components can turn the same-looking
local multiplicity into a qualitatively different obstruction.

### 2. Certificate complexity depends on the response coordinates

**Source:** THM-3320's integral versus affine common-table duals.

**Target:** THM-3326's abstract torsion block versus its marked bridge.

**Map:** restore one omitted coordinate and recompute minimal support.

**Preserved predicate:** infeasibility in the status table; exactness in the
Hamiltonian cokernel.

**Destroyed information:** an affine constant in the first case and the
formal source jet in the second.

**Needed sidecars:** the rational affine coefficient; the response polynomial
`F`.

**Cheapest decisive tests:** THM-3320's two apparent support-three packets;
`d=3,f=x+a x^2`, where `F=2(1-aT)`.

This explains why a minimal support number should never be quoted without
naming the ambient response space.

### 3. Well-founded selectors can replace large union bounds

**Source:** THM-3325's half-turn activation graph.

**Target:** future large-family status compilation after THM-3320.

**Map:** orient proof obligations by an operation that strictly decreases a
well-founded invariant, then cache certificates by orbit rather than by row.

**Preserved predicate:** existence of one safe shift or one exact status
contradiction.

**Destroyed information:** none if the arrows retain the full packet labels;
an unlabelled DAG would be insufficient.

**Needed sidecars:** `nu_2` and coordinate labels in the first setting;
capacity signatures and `S_4` orbit labels in the second.

**Cheapest decisive test:** the next nineteen-row `gcd72/L7056` family.

This third connection is a proposed scheduling mechanism, not yet a theorem
about the nineteen-row family.

## Repaired near misses

1. **Boundary sufficiency at composites is false by millions, not by one
   curiosity.**  The vector `(0,1,0,...)` was only the smallest witness to the
   full family `(2)`.
2. **The THM-3318 bridge is not deformation-invariant.**  Its torsion lengths
   survive, but the bridge records the nonlinear source jet.
3. **Gradient smoothness does not imply finite response torsion.**  The smooth
   two-root family `g=x^2(1+x)^2` has a free unit-response submodule.
4. **Integral template arity is not full-dual arity.**  An affine constant
   lowers two fourth-prefix circuits from apparent support three to two, while
   the row-138 hostile proves that support three still genuinely occurs.

## Next proof-facing experiments

1. **Physical LRC lift:** construct a speed-to-residue map landing in, or
   transverse to, the even flat while retaining cell width and endpoint
   origin.  Without this sidecar THM-3325 cannot touch LRC(14).
2. **Beyond `2p`:** replace the single half-turn graph by the divisor-lattice
   activation complex for `n=pq` and prime powers.  Test whether a lexicographic
   valuation gives a well-founded selector or whether a genuine cycle is the
   first obstruction.
3. **Next projected prefix:** compile the nineteen `L7056` rows by capacity
   signature and `S_4` orbit, retaining the affine dual coordinate from the
   start.
4. **Multi-root divergence:** THM-3326 classifies `theta`, not the canonical
   divergence class `mu`, in the multiple-root branch.  Compute its residue
   vector and determine whether `K[P]mu` has a free part plus root-supported
   torsion.
5. **Quadratic vertical perturbation:** for `P=f+g z+h z^2`, identify the
   algebraic differential replacing `dx/g` and ask whether branch residues or
   trace residues give the next exact response trichotomy.

The session did not produce a claimed route from one conjecture to the other.
It produced two reusable mechanisms—support-sensitive response and
well-founded activation—and one exact projected advance, each with the lost
coordinates and next decisive test exposed.
