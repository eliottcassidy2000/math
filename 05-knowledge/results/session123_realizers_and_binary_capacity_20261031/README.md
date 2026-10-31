# Session123: realizer algorithms and a complete binary capacity law

2026-10-31. **Scoped proofs, exact algorithms, and independent audits.**
Baseline: `683d9e32c0071ee028a010759a1fec40141fdbf8`.

## Main results

1. **[One exact formula for every capacity budget](monoidal_ramsey/ALL_BUDGET_KERNEL.md).**
   For the specified singleton/join/union/full-C5 class,

       F(A,B)=sum_(0<=i<A,0<=j<B) 2^popcount(i & j).

   Actual optimal graph constructors, origin translation, and a fast digit
   algorithm accompany the uniform parity-lifting proof. The explicitly
   enlarged [hereditary class](monoidal_ramsey/HEREDITARY_COMPLETION.md) has
   the same capacity, though it has more objects.
2. **[Complete minimal generators and binary tensor](tensor/TENSOR_AND_GENERATORS.md).**
   The irreducible budgets are exactly1 through the least power of two at
   least A. Capacity matrices are determinant-one integral Gram matrices.
   [Independent capacity/tensor audit](audit_ramsey/AUDIT.md).
3. **[Polynomial-bit common-time realizer](lrc/CONSTRUCTIVE_CONDUCTOR.md).**
   Exact integral bisection handles strict budgets and equality without
   enumerating huge periodic combs. A [full family selector](lrc/FAMILY_SELECTOR.md)
   constructs witnesses for every previously certified eight-core/five-tail
   family, including arbitrary large binary inputs. [Independent audit](audit_lrc/AUDIT.md).
4. **[Protected shared-gate repair](circuits/RESOURCE_COST_AND_BOUNDARY_KERNEL.md).**
   Retaining externally consumed products P in a forest cone, the rewrite
   saves |E\P|-nu(E\P) gates. Its exact obstruction is a matching; simultaneous
   output kernels preserve protected boundary variables. [Independent audit](audit_interfaces/AUDIT.md).
5. **[Realizability and computable sections](realizability/WITNESS_REALIZABILITY.md).**
   External surjectivity is distinguished from internal/tracked realization.
   An exact integer-floor compiler supplies a separate fixed-dimension
   selection theorem with classical Lenstra input.

## The categorical mechanisms and their limits

[Closed profile convolution](monoidal_ramsey/CLOSED_PROFILE_STRUCTURE.md) has
an actual internal-hom adjunction, but it also exhibits residual collapse
and a7-versus8 failure of careless saturation. [Same-clock predicates](lrc/SAME_CLOCK_CATEGORIES.md)
have a Boolean residual only when endpoint information is retained; they
also give an explicit continuous phase-section obstruction. Rational
computability and continuous choice are different requirements.

No new realizability topos is asserted. None of these results settles
unrestricted Ramsey numbers, LRC, or general NP circuit lower bounds.
[Primary sources and prior-art boundary](SOURCES_AND_BOUNDARIES.md).
[Research synthesis and remaining scope](SYNTHESIS.md).

## Reproduce

    python3 -B verify.py

The master command runs all15 declared finite checkers normally and with
`-O`, compares every stored receipt, and validates the sealed files.
`--checkpoint-directory PATH` additionally retains each successful output.
The packet uses the Python standard library. Classical realizability and
fixed-dimension integer-programming theorems are cited inputs, not implemented
external solvers. The inherited Session122 universal family certificate is
not silently counted as a new replay.

## Inheritance and concept board

* Anchor: common clock predicates -> actual witness, retaining endpoints.
  Hostile: separate successful attachments may have no common time.
* Niche: capacity profiles -> exact binary kernel and constructors, retaining
  the mixed-parity cover. Raw internal hom adds no missing witness data.
* Wildcard: independent bits tensor exactly; a graph-specific inequality is
  still needed to turn a numerical kernel into exact realizable capacity.
* Realizability: external existence, a tracked section, and a continuous
  section have different quantifiers and different counterexamples.
* Resources: protected outputs and sharing survive a rewrite; coordinate
  compression alone does not pay original-basis gate costs.

Method cards: retain observer/recurrence/finite head; correct the object before
sharpening the technique; search the statement before claiming priority.
