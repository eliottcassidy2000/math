# The exceptional critical deck is an algebraic etale two-clutch germ

**Status:** structural synthesis for
[THM-3319](../01-canon/theorems/THM-3319-exceptional-quadratic-two-clutch-etale-persistence.md).
The exact evidence is the
[two-clutch tangent scout](../04-computation/jc_exceptional_quadratic_two_clutch_formal_persistence_scout_20260803.py)
and its
[frozen output](../05-knowledge/results/jc_exceptional_quadratic_two_clutch_formal_persistence_scout_20260803.out).
Nothing here proves `JC(2)` or `DC(2)`.

## Inheritance and the change of question

THM-3306 found a transverse coefficient-pair base locus on the fixed
`C=c+x,E'=1` slice.  THM-3309 pulled it to a nonsplit quadratic common-root
deck and found the true obstruction:

```text
P_x=P_z=0, so gradient unimodularity fails before mu(P) exists.       (1)
```

The canonical hostile was that a resolving cover can expose rather than cure
a defect.  The corrected near miss was the temptation to call the inverse
different a Keller cofactor.  The least-used sidecar was the *physical clutch
plane* already present in THM-3289:

```text
C=c+d x,                    E'=k.                         (2)
```

The previous priority said to release one of `d,k` and ask whether the deck
disappears.  Transversality suggests a stronger question.  The base ideal has
two equations in the two internal coordinates `(x,c)`.  Why freeze either
external slope?  Release both and let `(x,c)` move.

## The exact mechanism

Let `a,b` be the two coefficients of the linear subresultant row.  At the
degree-36 residue field `A_i`, the fixed point satisfies

```text
a=b=0,                     det partial_(x,c)(a,b)=a_x b_c in A_i^*.   (3)
```

For `p=d,k`, the physical derivatives include the chain rule
`partial_d C=x`.  The infinitesimal equations are

```text
a_x dot x_p+a_p=0,
b_x dot x_p+b_c dot c_p+b_p=0.                            (4)
```

Both exact accessory-field computations give four nonzero velocities.  More
importantly,

```text
dot x_d dot c_k-dot x_k dot c_d !=0.                      (5)
```

Thus the clutch plane is not absorbed through one hidden scalar direction:
it moves the critical base point with full tangent rank.  More is already
available before completion.  The relative Jacobian criterion makes the
algebraic zero scheme `V(a,b)` etale over the `(d,k)`-plane at the degree-36
point.  Completing that algebraic germ gives the unique series
`x(d,k),c(d,k)`.

On that algebraic germ the linear row vanishes, the quadratic row stays
unit-leading, and the cubic degrees stay three.  The quadratic row is
therefore the exact gcd over the germ's function field.  Its discriminant is
a unit with nonsquare residue.  Idempotents in the local finite-etale algebra
are detected on the residue field, so the nonsplit special fibre lifts as a
connected rank-two algebraic cover rather than two chosen branches.

This changes the geometric verdict:

```text
fixed isolated critical pair                REFUTED as the right picture;
algebraic etale two-parameter critical germ PROVED in the declared fields. (6)
```

The fixed `(x,c)` hostile remains useful.  If either `d` or `k` moves while
the internal point is frozen, at least one linear-row coefficient immediately
becomes nonzero.  The deck persists only because the exact implicit tangent
repairs that failure.  This distinguishes genuine deformation from a copied
fixed certificate.

## Connection contract

```text
source:      fixed transverse base point plus physical clutch plane;
target:      connected rank-two algebraic critical cover;
map:         project V(a,b) to (d,k), then complete at the degree-36 point;
preserved:   quadratic gcd, deck exchange, true gradient vanishing;
destroyed:   global component identity, rational section, distant walls;
sidecar:     moving internal base point (x(d,k),c(d,k));
test:        exact two-clutch tangent determinant and residue discriminant. (7)
```

The determinant in `(5)` is the new coordinate.  Trace, norm, or one branch
alone would not show that both physical directions survive independently.

## Two distinct first-obstruction gates

The concurrent
[THM-3318](../01-canon/theorems/THM-3318-hamiltonian-divergence-torsion-ladder-for-x-plus-xr-z.md)
provides a useful orthogonal control.  Its family has a polynomial gradient
Bezout row, so `mu(P)` is defined; the obstruction is a nonzero special-fibre
torsion class with exact annihilator `(P^r)`.  THM-3319 instead has a connected
critical deck on which the gradient itself vanishes.

```text
THM-3319: gradient ideal proper  -> mu(P) undefined;
THM-3318: gradient ideal unit    -> mu(P) defined but nonzero torsion.      (8)
```

This is a genuine taxonomy, not a reduction between the families.  The two
objects live on different carriers and THM-3318 explicitly supplies no map to
the exceptional deck.  The common research procedure is to test the first
gate before computing a later invariant.

## Cheapest next tests

1. **Global-component test.**  Saturate and eliminate `(a,b)` over
   `K_i[d,k]` to identify the closure containing the etale germ, its degree,
   and the first locus where the projection stops being etale.
2. **Ramification/owner test.**  Compute the quadratic discriminant and the
   finite-clutch/owner resultants along one exact algebraic curve through the
   etale germ.  The present theorem gives only local units.
3. **Gate stratification.**  For new planar source-fibre families, first split
   the parameter locus into gradient-proper and gradient-unimodular strata;
   compute the Hamiltonian divergence class only on the second.
4. **Hostile against global optimism.**  Search for the first parameter value
   where the component closure meets a wall, ramifies, or changes gcd degree.
   Local algebraic persistence is not global continuation.

The reusable move is a version of “audit the next native operation”: when an
exceptional point is transverse in its internal coordinates, release all
available external directions before declaring it isolated.  The correct
object may be the algebraic etale graph of the base ideal, with the deck
carried over it, rather than one slice or one selected root.

## Reproduction

Run

```text
python3 04-computation/jc_exceptional_quadratic_two_clutch_formal_persistence_scout_20260803.py
python3 -O 04-computation/jc_exceptional_quadratic_two_clutch_formal_persistence_scout_20260803.py
```

and compare LF-normalized bytes with the frozen output.  The companion uses
both symbolic chain-rule differentiation and exact central differences for
the tangent packet, has no Python assertion nodes, and contains no floating
literals.
