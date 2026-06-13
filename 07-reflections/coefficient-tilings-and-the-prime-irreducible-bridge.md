# Coefficient Tilings And The Prime-Irreducible Bridge

The useful turn in this session was to stop asking whether a polynomial
"is" a tournament.  That direction collapses too quickly.  The better object
is a quotient:

```text
fixed-path tournament tiling -> coefficient profile -> polynomial certificate.
```

The user's triangular model lands cleanly in the fixed Hamiltonian-path gauge.
For `N` vertices, the gap diagonals have sizes `N-1, N-2, ..., 1`; for degree
5, this is exactly the `5,4,3,2,1` coefficient tower.  The base path can be
treated as the constant term, while higher diagonals become mutable
coefficients.

Two maps matter immediately.  The count map `c_d = #forward arcs on gap d`
is a Cohn digit map.  In base `b >= N`, prime values of the digit profile
certify irreducibility of the digit polynomial.  The centered map
`A_d = 2c_d-(N-d)` is the sign/magnitude map: coefficient signs are layer
majorities, and fixed coefficient magnitudes are slices through the tiling
hypercube.

The finite atlas is small but already clarifying.  Full grids through `N=7`
show zero Cohn mismatches, as they should.  Fixed-path `N=6` has 120 profiles
over 1024 tilings; 57 positive-degree Cohn-prime profiles; 96 centered
irreducible profiles; and 859 centered-irreducible tilings by weight.  More
importantly, 91 of the 120 fixed-path profiles contain multiple `H` values,
with a maximum spread of 34.  This is the scalar/support split in miniature:
the coefficient polynomial is a shadow, and the tournament fiber is where the
lost support data lives.

This reframes HYP-2448.  Singh's value-factor criterion, Cohn's digit
criterion, and Iravanian's real-factor recombination are not competing
stories.  They are different retained channels around one quotient:

```text
integer polynomial -> integer values -> factor/prime certificates
real/finite-field atoms -> recombination side channel -> integer factors
tournament tiling -> coefficient quotient -> irreducibility shadow
```

Church's product-quotient warning fits the same shape: scalar supersingularity
can pass while the retained diagonal Frobenius channel blocks unirationality.
Here, scalar coefficient data can pass irreducibility tests while the tiling
fiber still carries the Hamiltonian-path and support information needed for
LRC14 or `[72,36,16]`.

The most actionable proof route is now three-stage:

1. Freeze the exact Cohn-diagonal lemma.  This gives a clean theorem: prime
   base-value plus diagonal address implies irreducibility of the corresponding
   digit polynomial.
2. Study centered magnitude slices.  The parity-minimum slice at `N=6` has
   only 8 distinct polynomials and all are irreducible in the pilot; that is a
   natural target for a small structural classification.
3. Audit fibers, not just profiles.  For LRC14, fibers should remember which
   resource blocked which denominator.  For `[72,36,16]`, fibers should
   remember support/design/matroid realizability beyond the scalar weight
   enumerator.

The "turn it on its head" version is also worth keeping: vertices need not be
polynomials or coefficient profiles.  They can be sign assignments at fixed
magnitude, roots ordered by argument, finite primes ordered by factorization
type, Newton edges ordered by slope, derivative obligations, or code support
moves.  The procedural tournament in the script ranked these by a majority of
novelty, testability, bridge strength, support transfer, and risk.  It was
transitive, which is a useful first spine rather than a disappointment:
diagonal Cohn map first, centered magnetization second, fiber entropy third.

The next good session should try to make the centered-slice phenomenon less
anecdotal.  If some magnitude vectors force irreducibility because every
possible factor split violates a layer parity or layer-majority constraint,
then the user's tiling picture would become more than a metaphor: it would be
a genuine finite obstruction calculus for factorization.
