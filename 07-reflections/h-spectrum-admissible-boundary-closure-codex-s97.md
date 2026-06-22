# H-spectrum admissible boundary closure

The boundary-function analogy sharpened the productive-closure question.  The
issue is not whether a quotient is beautiful, equinumerous, or graph-theoretic;
the issue is whether the approach to the boundary preserves the invariant we
want to read.

For bounded analytic/harmonic functions on the disk, Fatou-type theorems give
stable almost-everywhere boundary limits along admissible approaches.  Lindelof
principles add that if a bounded holomorphic function has a curve limit at a
point, the non-tangential limit agrees there.  In the opposite direction,
Bagemihl-Seidel and Kaczynski-style boundary-function theory studies how
flexible curvilinear approach families can realize wild boundary data,
including for functions with values in metric or Riemann-sphere-like targets.

That is the right moral for the `H` spectrum:

- **Admissible approach:** strong-component condensation plus OCF packets.
  This has an actual factor ledger, `H(T)=prod_i H(C_i)`, so a target value
  decomposes into a finite list of factor routes.
- **Boundary regularity:** once a single strong core is forced, THM-115's Moon
  pancyclic bound raises `H` above `21` for `n>=9`, while THM-079 handles
  the finite base.
- **Wild curvilinear approach:** fixed-path even-graph smoothing/contraction.
  HYP-2872 showed these moves preserve GF(2) parity data but can collapse many
  `H` values to one core and can even move `H` upward under contraction.

So productive closure is an admissible-boundary statement, not a minor-closed
shadow statement.  The fixed-path/even-graph quotient remains useful as an
address quotient, but it cannot be promoted to a proof carrier without a
preservation theorem, and S96 shows that theorem is false for the obvious
degree-2 and contraction moves.

The sharper S33 minor-order prompt fits after correction: `contract` should be
read as strong-component condensation, where `H` multiplies and the factor
ledger is explicit.  It should not be read as arbitrary GF(2) contraction in
the even-graph shadow.

The S97 arithmetic ledger is small but clean.  For `H=7`, the prime target
forces a single forbidden strong value, closed by THM-200/THM-343.  For
`H=21`, the product route `3*7` is closed by the same `H=7` theorem, and the
single-core route is closed by THM-079 for `n<=8` plus THM-115 for `n>=9`.
This is exactly the productive closure needed for the known permanent gaps.

The important correction from concurrent S33 is that this is not a `7Z`
ideal.  Strong `n=7` tournaments already realize `H` values divisible by `7`,
including `35,49,133,147,175`.  Therefore the hidden structure is finite
low-boundary obstruction followed by direct strong-core re-entry, not a global
multiplicative prime obstruction.

## Assumption challenge / Tournament Analysis

Candidate vertices considered during this session:

- tournament vertices;
- fixed Hamiltonian path free arcs;
- even-graph triangle-basis bits;
- strong components;
- OCF odd-cycle conflict components;
- factor routes for a target `H`;
- large-core Moon-boundary obligations;
- wild quotient operations;
- boundary approach classes;
- proof obligations.

The chosen proof-carrier vertices are

```text
factor route -> single strong core -> finite base -> Moon boundary -> wild quotient guardrail.
```

The proof-carrier tournament is transitive.  Its leading vertex is the
strong-component/OCF factor route because it preserves the `H` predicate.  The
fixed-path/even-graph quotient is demoted to an address layer because it
destroys `H` while preserving only parity/cycle-space information.

Challenged assumption: tournament vertices and even-graph edges need not be the
right vertices for the proof.  For this closure problem, the correct vertices
are routes and approach classes.

## Sources checked for the boundary analogy

- Bagemihl-Seidel boundary properties: https://eudml.org/doc/169458
- Kaczynski boundary functions: https://www.jstor.org/stable/1995093
- Fatou/Lindelof background in Fulkerson dissertation: https://oaktrust.library.tamu.edu/bitstream/handle/1969.1/86034/Fulkerson.pdf?sequence=1
- Radial limits / Fatou theorem summary: https://www.numdam.org/item/10.5802/crmath.609.pdf
