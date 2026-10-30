# Odd descent, dyadic obstructions, and local relaxation certificates

**PROVED scoped synthesis / FINITE-EXACT controls / CITED algorithmic background /
OPEN major applications. Same-agent proof review; no formal proof build or
literature-priority claim.**

Odd cycles, odd functions and odd integers meet at a precise interface:
**sign reversal introduces a factor of two, while an odd number of compatible
copies can be rounded without losing a fixed constraint.** The action, the
constraint, and the marked values must be specified. Oddness alone is insufficient.

The strongest assembled result is a complete finite feasibility test for
constraints preserved by odd coordinatewise medians: first find the maximal
support face of their local marginal polytope, then test the powers of two in
the Smith normal form of its equality matrix. An odd-denominator solution
rounds to an actual assignment. The proof is self-contained in
[PROOFS.md, P3–P4](PROOFS.md); symmetric-polymorphism and LP/affine-relaxation
methods are established background, explicitly credited in [SOURCES.md](SOURCES.md).

The bridge also covers binary affine constraints: use XOR of an odd number
of copies in place of median. The repo's seven-coordinate Hamming syndrome
fibres preserve all6528 tested triples under XOR, but median fails already
on the three valid masks0,7,25. A constraint-preserving symmetric operation,
rather than oddness of values alone, is the shared ingredient.

Other exact conclusions:

- **Sharp cover-denominator exponent (P4b).** For any finite local constraint
  model with feasible marginal polytope, the maximal-support Smith data gives
  the smallest power of two dividing a satisfying cover degree. This does
  not determine the smallest total degree or its connectedness.

- **Odd-cover descent and two-cover compression (P1).** For fixed strictly
  monotone edge functions, every odd-sheet matching-cover solution gives a base
  solution by fibre medians. A solution on any finite cover gives a solution on
  the canonical two-sheet cover by retaining minima and maxima. This is a
  feasibility certificate for a fixed network, not a carrier closed under
  arbitrary later attachments.
- **Signed-affine classification (P2).** For integral equations
  `x_v = s_e x_u + b_e`, `s_e in {+1,-1}`, integer feasibility is equivalent
  to real feasibility plus feasibility modulo two. The minimum possible finite
  cover degree is absent, one, or two. A tree-coordinate content formula counts
  solutions over every `Z/q` and classifies their odd-prime support.
- **Weighted boundary (P5).** For `d` odd and any nonnegative convex unary cost,
  median cost is at most `2d/(d+1)` times average cost, sharply. For a
  marked-point indicator the ratio can instead be `d`, unbounded.
- **Integral odd/even decomposition (P6).** Every exchanged pair contributes
  exactly one parity compatibility condition. Odd-degree transfer preserves
  two-primary torsion, subject to an actual restriction/transfer identity.

[RECOVERY.md](RECOVERY.md) connects June's odd-function dictionary, the repaired
odd-cycle digit tower, September's odd-residue floor formula, October's
odd-cover and shared-clock results, and current Collatz/LRC certificates.
[PROPOSALS.md](PROPOSALS.md) gives concrete next proof obligations for
clocked rounding, collision-safe labels and marked tournament responses.

The exhaustive signed-affine census contains 15,756 systems: 4,504 integral,
2,128 real but half-integral only, and 9,124 real-infeasible. All 15,625
three-sheet presentations on the declared three-vertex universe obey descent;
64 of 729 two-sheet presentations repair an integral obstruction. Further
checks cover every modulus 2 through 12 on 131 systems, 625 local support-face
problems, convex costs and the unchanged 81/27 tournament control.
Among the425 real-feasible support-face examples, the minimum denominator
exponent is0/1/2 in377/36/12 cases. The checker performs196,773 local guard
evaluations, with identical normal and optimized receipts.
The precise counts and universes are in [certificate.json](certificate.json)
and [AUDIT.md](AUDIT.md).

The failure examples are useful: unclocked Collatz, square-sum relations with
changing square labels, all-different constraints, wrapped runner phases and
unlabelled tournament majority each break an explicitly identified hypothesis.
No proof of Collatz, LRC(14), the graceful tree conjecture, Sha finiteness, or
a new unrestricted A000568 recurrence follows from this packet.
