# Lucas-prime symmetries, full marking holonomy, and joint-port moments

2026-10-28. **PROVED scoped / CITED classification and prime-divisor density /
FINITE-EXACT controls; same-agent alternate audit.** No external-priority claim.
Initial source4961bb2f354e; integrated concurrent source59cb22d0e9, including
the recovered ordered-run witnesses. No general Collatz, LRC(14), planar
Jacobian, or six-colour theorem is claimed.

## Direct progress

1. [All Lucas-prime reflection groups](REFLECTION_GROUPS.md). If p=L_n is
   prime, n is odd and n>=7, then **any p-cycle plus word reflection** generates
   A_p, except that n==11 mod12 gives S_p. The fixed-word count F_((n+3)/2)
   lies strictly between1 and sqrt(p); the proof excludes every proper case
   in Jones's primitive-group classification. The arithmetic eleven-state
   PSL(2,11) case is exceptional. Even power-of-two indices admit a parallel
   edge-reflection theorem with a different, explicitly chosen cycle carrier.
2. [Prime-specific projective clocks](PRIME_ATLAS.md). The six-mirror's
   fixed-point type and PSL/PGL membership are classified modulo24. Lucas
   divisibility is exactly a midpoint condition in the cyclic group generated
   by the Fibonacci matrix in PGL2(F_r). The order, initial valuation, and
   quadratic splitting character give exact divisor towers. Moree's cited
   density2/3 concerns prime **divisors**, not prime-valued Lucas terms.
3. [Full reversal-marking holonomy](FULL_HOLONOMY.md). With K the inverted
   twisted subgroup and C the marked centralizer, the loop group is
   C intersect <K>. For W_k=C3 iterated k times it is the entire C_k,
   including the nonabelian order243 group at27 vertices. Every correction
   has a constructive loop of at most3[1+(k-1)3^k] square-root transport
   edges. This counts the stated re-marking moves, not tournament complements.
4. [All-order joint-port recurrence](PORT_MOMENTS.md). At fixed r, the
   number of centralizer orbits on ordered r-tuples is a degree-r polynomial
   in k with leading coefficient2^r r!. One/two-port counts are2k+1 and
   8k^2+1. The excess over independent marginals is exactly4k(k-1), the
   variance of a specified fixed-vertex distribution. A critical ternary
   Galton-Watson process realizes the uniform-automorphism law exactly.
5. [The76-state composite boundary](COMPOSITE_76.md). The period-nine ring
   is F4 times F19. Reflection is Frobenius on F4; at each nonzero F19 label,
   all four hidden tags give different reflected F19 outputs. Eight chart
   types with nine phases reconstruct the72 nonzero states, with four
   additional zero-fibre states. Composite degree also has an explicit
   imprimitive hostile to the low-fixed-point group test.

## Inheritance and concept board

**Anchor:** robust Lucas-prime growth mechanisms rather than extrapolating
eleven. **Niche:** the entire reversal holonomy group and its repair cost.
**Wildcard:** the composite76 chart and fixed-port branching moments.

Closest proved mechanism: the all-period reader and source-preserving ideal
lifts in [Session74](../session74_eleven_mechanisms_20261027/session74-clocks/RESULT.md),
and the square-root transport of [Session90](../session90_even_graphs_odd_descent_20261029/carriers/TRANSPORT_AND_HOLONOMY.md).
Canonical hostile: the [A29 word certificate](../lucas_mirror_gap_spaces_20261028/LUCAS_MIRRORS.md).
Corrected near miss: the March Fibonacci/Lucas-prime independence heuristic,
repaired by [the October11 audit](../abc_difference_families_20261011c.md).
Least-used sidecars: the actual centralizer permutation and the discarded
F4 coordinate, beyond their cardinalities.

| Live object | Invariant / operation | Failure tested | Exact repair or outcome |
|---|---|---|---|
| Lucas words | reflection fixed count; addition cycle | eleven closure extrapolation | all-prime A/S classification |
| Projective mirror | determinant squareclass; fixed lines | confusing two involutions at11 | retain action and twist |
| Prime divisors | projective order; valuation lift | support versus multiplicity | (d_r,e_r,chi_r), plus lifted state when queried |
| Reversal marks | square-root transports | endpoint-only port alignment | full centralizer correction, charged loop |
| Shared ports | diagonal group orbits | products of marginal quotients | moment recurrence and joint charges |
| Period-nine CRT | Frobenius and chart phase | F19 alone forgets mirror output | four-valued tag; eight-chart table |

Every connection specifies its actual carrier and allowed operations in the
linked note. The full group does not preserve all original metric, toggle,
source-weight, or integer-clock guards. The probability law in the branching
model is uniform group elements, not a law on actual positive Collatz sources.

## Prior and incoming work used, with scope retained

- Historical March `fib_lucas_primes_s90bo.out` and
  `fibonacci_lucas_primes_s116l.out` supplied Pell/split-prime provenance.
  Their heuristics were not promoted; current corrections take precedence.
- June [HYP-2187](../../hypotheses/HYP-2187-equinumerosity-equidecomposability-fiber-bridge.md)
  motivated a finite piecewise-group criterion. The allowed group and pieces
  are now specified; no Euclidean scissors equivalence is imported.
- Session90's order81,27-mark hostile is strengthened to the complete
  holonomy group at every depth; its existing nonflatness and group formulas
  are credited rather than rediscovered.
- Session91's golden-fibre controls prevent treating the even golden quotient
  as a field of cardinality L_n. Its characteristic and source tags matter.
- Concurrent Session92's [joint-closure algebra](../session92_joint_order_carriers_20261029/algebra/JOINT_CLOSURE_ALGEBRA.md)
  explains why shared observations need product closure. Its newly recovered
  [six-order witness](../session92_joint_order_carriers_20261029/runs/RECOVERED_RUN_WITNESS.md)
  has the same C3 orientation quotient as three distinct ports in one child.
  The surrounding Hamilton-count problem remains a separate operation.
- Session92's divisor extension-count repair preserves a declared source law
  under reassociation. Our based holonomy repair preserves an actual marked
  permutation under transport. Both retain missing composition data; no
  category equivalence or equality of their state spaces is claimed.

This is targeted recovery, not a claim to have read every historical mention.
META-PATTERNS used: search the statement before the method; retain local
incidence through counting; keep operation-compatible sidecars. No new
universal method card is promoted from this session alone.

## Reproduction and frontier

Run [reproduce.py](reproduce.py) from any directory. It checks producer and
alternate receipts in normal and optimized Python modes. The producer uses
the standard library; the alternate path uses SymPy for the small11/29
group orders, explicit Jordan certificates at47/199, and a separate recursive
construction of the centralizers. See [AUDIT.md](AUDIT.md) for exact universes,
dependency limits, hashes and proof-boundary checks.

Next meaningful targets: canonical tuple representatives with weighted native
operations; optimal marking-loop lengths; compatible CRT charts at general
composite periods; and interaction between valuation/order records and prime
value sieves. The port **counting** recurrence is closed here. Lucas-prime
infinitude, an A000568 speedup for arbitrary tournaments, and universal
Collatz entry remain open in this work.
