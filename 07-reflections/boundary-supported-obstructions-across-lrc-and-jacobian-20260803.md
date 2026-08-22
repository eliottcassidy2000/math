# Boundary-supported obstructions across the runner and Jacobian lanes

**Status: SYNTHESIS / NEXT-OPERATION NOTE.**  This note records a shared
research move, not a reduction between conjectures.  The prime
micro-staircase theorem, the `z1=216` gcd prefix, and the Hamiltonian cokernel
family all improved when the forgotten boundary data was retained as an
explicit small object.  There is no map from runner residues to polynomial
Hamiltonian classes, and no consequence transfers between the conjectures.

## The three new boundary objects

### Prime micro-staircase

For a normalized residue vector `w`, right-boundary coverage at
`a=-j^(-1)` forces

```text
{l/w_l,(l-j)/w_l:l in supp(w)}=F_p.                     (1)
```

The collision set

```text
X={l/w_l:l in supp(w)}                                  (2)
```

is the sidecar that makes the quotient faithful.  Off `X`, the maps
`l -> l-yw_l` must permute the support.  Their product polynomial is then a
nonzero constant off `X` and zero on `X`, contradicting its low-degree field
sum.  THM-3316 therefore closes every prime modulus in the finite cell model.

At composite `14`, inverses and the field sum fail together: the one-defect
vector `(0,1,0,...,0)` passes all right-boundary tests.  The missing coordinate
is not merely support size; it is the CRT/valuation type of `(l,w_l)` and the
collision set on each unit stratum.

### Projected `z1=216` wall

Canonical [THM-3313](../01-canon/theorems/THM-3313-projected-k3-z216-third-ruler-cost-prefix-multicover-closure.md),
with a concurrent independent audit, closes the next four gcd-family rows and
changes

```text
ledger 373161 -> 373157,
wall 357 -> 353,
families 33 -> 31.                                      (3)
```

The one formerly weighted status closes only after coupling two nested load
thresholds.  A high-tail point outside a coordinate union must contribute an
extra incidence, producing the exact gap `510`.  Here the sidecar is the
location of high-tail mass relative to the coordinate union; separate scalar
tail counts forget the overlap that proves emptiness.

### Hamiltonian cokernel

For

```text
P=x+lambda x^r z,
```

THM-3318 identifies two classes in `C_P=K[x,z]/D_P(K[x,z])`:

```text
Ann(mu)=(P^r),       Ann(theta)=(P^(r-1)),
P mu=-(r-1)theta.                                       (4)
```

Both become exact off `P=0`, yet their Laurent primitives have unavoidable
`x=0` poles.  The sidecar is the exact `P`-adic torsion order, not just the
fact that generic localization solves the equation.  The generic chart
forgets precisely the extension obstruction carried by the special fibre.

## Shared method, not shared object

The common operation is:

```text
localize or quotient
    -> solve the generic/interior problem
    -> retain the smallest exceptional-set sidecar
    -> test whether its mass/torsion can vanish globally.                (5)
```

| lane | generic simplification | exceptional sidecar | decisive test | lost information |
|---|---|---|---|---|
| prime micro | scale shifts in `F_p` | collision set `X` | low-degree field sum | width, endpoints, physical speed |
| composite wall | threshold marginals | high-tail overlap + gcd family | two-threshold incidence | chronology, owners, physical cover |
| Jacobian | invert `x` or `P` | `P`-primary cokernel order | pole valuation | polynomial extension across `x=0` |

This pattern is useful only at the level of proof design.  The source rings,
targets, predicates, and failure modes are different.  In particular,
`P`-torsion is not a runner gcd class, and the runner collision set is not a
cohomology class.

## Three next operations

1. **Composite micro Bockstein.**  Split support pairs `(l,w_l)` by
   `gcd(l,n),gcd(w_l,n)`.  On each unit CRT factor form the collision set
   `X_d={l/w_l}`.  The prime proof becomes the zeroth-moment equation
   `|X_d|=0` in the residue field.  When that equation is characteristic-
   torsion, add one adjacent Farey interior cell and ask whether it supplies a
   first-moment correction incompatible across CRT factors.  The `n=14`
   one-defect hostile is the mandatory negative control.

2. **Cokernel torsion spectroscopy.**  For a new gradient-unimodular `P`, do
   not begin with a large mate ansatz.  Compute the `K[P]`-orders of
   `theta=[1]` and the canonical divergence class, then locate the divisor
   that kills them.  A finite torsion block gives a bounded extension problem;
   a class surviving generic localization requires a different mechanism.
   The `r=2` member of THM-3318 is the hostile showing that a nontrivial deck
   is unnecessary.

3. **Nested-tail closure before more atlas expansion.**  Enumerate the
   four-bit two-threshold incidence inequalities suggested by the new
   `z1=216` state and replay them against the eleven weighted cores in the
   previous prefix.  A positive result would replace solver certificates by a
   short submodular taxonomy; a negative result should name the first
   three-threshold or owner sidecar before more rows are opened.

The strongest immediate target remains LRC(14), but the prime theorem changes
its diagnosis: scalar rigidity is not intrinsically a large-cell phenomenon.
The obstruction at `14` is specifically composite leakage plus the missing
speed/residue lift.  On the Jacobian side, the torsion ladder similarly says
that generic solvability is cheap; extension across the deleted divisor is
the load-bearing coordinate.
