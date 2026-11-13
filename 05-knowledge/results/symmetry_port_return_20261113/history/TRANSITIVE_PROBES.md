# Finite physical probes and simultaneous separation of port responses

**Status: PROVED.** The path-run substitution law is inherited from
[THM-3121](../../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).
The new interface here is a complete finite test using actual transitive
modules, and a simultaneous physical separator for finitely many nonidentity
guards. The historical prompt is the bounded simultaneous coefficient-grid
construction in
[THM-2270](../../../../01-canon/theorems/THM-2270-simultaneous-balanced-cut-relation-and-six-uniform-orientation.md).
This supplies an exact identity/separation tool, not a prescribed-target or
primality theorem.

## 1. The run-bound sidecar is essential

Let a fixed uniform context have variable nonempty module sites
`M_1,...,M_d`. A context may include a fixed named root outside these sites
and a fixed finite root-child program at fixed or fresh exterior ports.
Suppose its response has certified run caps `R_i`: in every counted object,
module `i` is split into at most `R_i` nonempty directed path runs. The
same fixed caps must work independently of module sizes.

This condition holds for one variable module with a fixed finite exterior.
For Hamilton paths it also holds when the quotient induced by all variable
sites is transitive: between two runs in one site there must then be a
fixed exterior vertex, so `R_i <= f+1` if there are `f` fixed vertices.
For a fixed number of ordered path slots, add the number of slots to the
corresponding separation budget. One may instead supply an exact finite
quotient-word certificate, as in the current clock and polarized frames.

An arbitrary cyclic quotient of unbounded variable modules need not have
such caps: `C3[M_1,M_2,M_3]` permits arbitrarily many balanced runs. The
theorem below must not be applied there with a guessed finite cutoff.

Write `c_r(M)` for ordered covers by exactly `r` nonempty directed paths.
Splitting paths into maximal module runs gives

```
F(M_1,...,M_d)
  = sum_(1<=r_i<=R_i) a_(r_1,...,r_d) product_i c_(r_i)(M_i).   (1)
```

For an actual count, the skeleton coefficients are nonnegative integers.
Differences or integer linear combinations of count observations have
integer coefficients of either sign. Fixed exterior endpoint or deletion
conditions merely restrict the finite skeletons, so use the same form.
Marks inside a variable module require extra endpoint-cover coordinates
and are not covered by (1) without that further argument.

## 2. Complete transitive test for a universal identity

For the transitive tournament `TT_n`,

```
c_r(TT_n)=r! S(n,r),                                  (2)
```

where `S(n,r)` is a Stirling number of the second kind. Indeed, partition
the vertices into `r` nonempty labelled sets; each set has exactly one
directed path in its inherited transitive order.

The matrix with rows `n=1,...,R` and columns `r=1,...,R` in (2) is lower
triangular with diagonal `1!,2!,...,R!`. It is invertible over the rationals.
Taking tensor products proves:

**Finite physical identity theorem.** Two responses with the run caps
above agree for every tuple of nonempty tournament modules if and only if
they agree on the actual transitive grid

```
(TT_(n_1),...,TT_(n_d)),       1<=n_i<=R_i.             (3)
```

Thus a proof of the run caps plus these finitely many exact physical tests
is a proof of the universal module identity, not an empirical fit. The grid
recovers the entire signed skeleton coefficient tensor by triangular
back-substitution. Conversely, a nonzero tensor has a distinguishing
physical substitution on (3).

Every probe module has Hamiltonian count one. Nonetheless the family of
such probes detects every bounded-run response coefficient. This is a
precise reason not to confuse small scalar `H` with a weak observer.

## 3. One physical substitution can satisfy finitely many nonzero guards

Let `F_1,...,F_s` be nonzero responses of form (1), with common caps `R_i`
after zero-padding their tensors. The problem is to make all `F_j` nonzero
at the **same** physical substitution.

By inclusion-exclusion for onto maps,

```
r!S(n,r)=sum_(k=1)^r (-1)^(r-k) binom(r,k) k^n,
                                                    n>=1. (4)
```

Therefore each transitive response is a finite exponential polynomial with
bases `1,...,R_i` in coordinate `n_i`. This change of basis is invertible.
Their product is not the zero exponential polynomial. For an elementary
algebraic justification, express every positive integer base by its prime
exponents. Finite exponential polynomials become ordinary polynomials in
one formal variable per prime and coordinate. Distinct bases are distinct
monomials, and this polynomial ring has no zero divisors.

Let

```
B_i={k_1...k_s : 1<=k_j<=R_i},       b_i=|B_i|.        (5)
```

The product has bases in `B_i`. Evaluating in coordinate `i` at
`n_i=1,...,b_i` gives an invertible Vandermonde matrix: the bases are
distinct positive integers, and multiplying the usual exponent-zero
Vandermonde columns by their bases preserves invertibility. Tensor these
matrices. A nonzero product cannot vanish on the whole grid. Hence:

**Simultaneous physical separator.** There exist actual module orders

```
1<=n_i<=b_i<=binom(R_i+s-1,s)                         (6)
```

such that every `F_j(TT_(n_1),...,TT_(n_d))` is nonzero. Exhausting this
finite grid is a constructive search. The last bound counts multisets of
`s` bases; multiplication can identify several multisets, so it is only an
upper bound.

This transports the **simultaneous** part of THM-2270, not its LRC equation
or six-uniform sign pattern. There the parameter grid chooses one integral
relation avoiding all balanced-cut hyperplanes. Here the grid chooses one
tuple of actual modules avoiding all listed response equalities. The
coefficient algebra changes from ordinary polynomials in signed relation
coordinates to exponential polynomials at positive physical sizes.

## 4. A small sharp simultaneity boundary

Use the actual contexts

```
W(M)=C3[M,1,1],        L_3(M)=C3 -> M,
L_5(M)=Q4 -> M,        H(Q4)=5.
```

Their responses are

```
H(W(M))=3c_1(M)+c_2(M),
H(L_3(M))=3c_1(M),     H(L_5(M))=5c_1(M).             (7)
```

The two nonzero differences on `TT_n` are

```
F_1(n)=2^n-2,      F_2(n)=2^n-4.
```

Both have run cap two. The identity-test grid `n=1,2` distinguishes each
one separately, but no member distinguishes both simultaneously: the first
vanishes at `1` and the second at `2`. Their product has bases `{1,2,4}`;
`n=3` is the first common separator, producing actual counts `9,3,5`.
Thus substituting the single-response grid for the simultaneous grid would
be a genuine logical error. The bound (6) is sharp in this small example.

### A six-probe certificate for the current physical clock frame

The inherited frame `Pi[M,N]` has ordered sites
`a,v,c,u,d0,d1,x0,f0,f1`, all arcs forward except `d0,d1->a` and
`f0,f1->c`; put `M` at `u` and `N` at `v`. Every return to `N` uses the
unique fixed gateway `a`. Every return to `M` uses an as-yet unused `a`
or `c`, since re-entering through `N` first requires `a`. Thus the actual
run caps are `(3,2)`.

The checker constructs all six graphs with `M=TT_n`, `N=TT_m`,
`1<=n<=3`, `1<=m<=2` and counts their Hamilton paths by vertex subset DP:

```
361   721
1081  2161
3241  6481.
```

Triangular inversion gives the complete nonempty-run coefficient tensor

```
361  180
360  180
120   60.
```

Using `Q_M(3)=3c_1(M)+3c_2(M)+c_3(M)` and
`Q_N(2)=2c_1(N)+c_2(N)`, it is exactly

```
H(Pi[M,N])=60 Q_M(3)Q_N(2)+H(M)H(N).                (8)
```

Formula (8) is inherited from the
[clock boundary](../../target_collatz_transitions_20261113/collatz/CLOCK_BOUNDARY.md),
not claimed anew. The point is that the certified run bounds make these
six literal probes a complete all-module proof certificate. No growing
module census or informal extrapolation is needed.

## 5. How this changes the current search

Use [the marked hypergraph congruence](MARKED_HYPERGRAPH.md) to discard
neutral changes which cannot affect a protected port. When an operation
releases a port star or a module incidence, keep its finite run tensor.
Equations (3) and (6) then provide two positive tools:

* certify an all-module neutral identity from a proved cap and a finite
  basis of physical tests;
* construct one physical probe where finitely many genuinely different
  port/context responses all separate, rather than testing each on a
  different graph.

Nonzero integer differences are not automatically `+2`, `-2`, or units
modulo a specified prime. Modulo a prime, distinct positive bases can
coincide or vanish, so the Vandermonde proof must be redone on the actual
residue carrier. No exact-target renewal, arbitrary gcd-one guarantee, or
Chan–Pak coverage follows from (6) alone.

`transitive_probes.py` checks the complete basis identity through cap six,
literal counts for (7), the sharp simultaneity control, a signed two-module
tensor reconstruction, and all six literal graphs for (8). Run it from the repository root; its exact receipt
is `transitive_receipt.json`.
