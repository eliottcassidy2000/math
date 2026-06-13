# Convolution Factor Capture And Hidden Coefficient Tilings

The useful move today was to stop treating the coefficient tiling as the whole
object.  A coefficient vector is only a boundary-total shadow.  If

```text
a_k = sum_{i+j=k} b_i c_j,
```

then the polynomial has a hidden two-dimensional lift: the multiplication grid
with diagonal sums `a_k`.  Reducibility is exactly the existence of a
nontrivial integral lift.  Irreducibility is the absence of that lift.

That makes the user's triangular picture sharper.  HYP-2450 projected fixed
Hamiltonian-path tournament tilings to coefficient profiles and magnetization
slices.  HYP-2451 adds the residue/valuation split-survivor carrier.  HYP-2452
asks whether those profiles have an integral factor tiling underneath them and
uses value-factor capture as a pruning score.  The pilot is intentionally
modest but solid: for primitive degree `<=5`, linear/quadratic factor search
is complete, and the stored scan has zero mismatches against Sympy across
`3856` degree-4 rows and `2016` degree-5 rows.

The Singh-style value witness fits as a pruning lens.  At an integer `m`, the
prime factorization of `f(m)` supplies tokens.  Any hypothetical factorization
must allocate those tokens among factor values.  Low `Omega(f(m))` therefore
means few possible factor slots and few allocations.  It is not a proof
without the paper's hypotheses, but it is a cheap way to decide which hidden
lifts deserve attention.

The residue and sign-cube tournaments are useful mostly as warnings.  The
simple residue tournament `r -> s` when `v_p(f(r)) < v_p(f(s))` detects fixed
divisor obstruction, but the score is too scalar and mostly transitive after
tie-breaking.  The sign-cube tournament by `f(B)` is also transitive.  The
data still matter, but the proof-bearing cycles must come from richer
obligations: convolution compatibility, Newton boundary layers, recombination
subsets, or support realization.

This is why the same idea points back to LRC14 and the `[72,36,16]` code.  A
raw denominator ledger or a weight enumerator is a boundary-total object.  The
hard question is whether a hidden incidence lift exists.  For LRC14 that lift
records which runners, owners, carries, and divisor fibers block which unit
twists.  For the length-72 code it records whether the scalar Type II
enumerator can be realized by binary supports and designs.  The scalar can be
beautiful and still fail because no lift exists.

The next good technical step is an ILP/SAT formulation of the convolution
equations.  It should start as a bounded degree-6 extension of the exact
degree-4/5 oracle, then branch to Newton-slope boundary lifting for sparse
multivariate polynomials.  If that works, the code-support problem starts to
look less like a mystical missing object and more like a family of hidden-lift
feasibility problems with increasingly rich coefficient rings.
