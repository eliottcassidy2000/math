# Independent review: prime-field rigidity and global periods

Date: 2026-10-19. Verdict: **PASS** for the final addendum, including its new
wheel-span section. No new chromatic bound or minimum-defect value is established.

## Mathematical checks

- The full norm-one graph over F_(p²), for prime p>=5, has rigidity rank
  2p²-3. The proof legitimately needs only three nonparallel direction pairs.
  Its affine-function step uses the prime-field assumption; it makes no claim
  for nonprime q. The exceptional p=3 rank is 12.
- The normalization is z=c+Rh, with c=(N(X_w-X_v)-1)/p and R the trace
  derivative. Reversing an edge changes both derivative factors' signs and
  changes neither its row nor its norm error. The straight-cycle error sum is
  -2, so equality forces exactly one error, of value -2, on each cycle.
  The two-direction elimination leaves 2p offsets and rank 2p-3.
- The seam's carry and residual signs are correct. In particular its w-carry
  is -1_(x=0)-1_(y=p-1), and the residual has x-period -1 and y-period +1.
  These are opposite manifestations of one remaining obstruction, not two
  independent missing conditions after all three cycle families are imposed.
- The final wheel argument is exact. A vanishing weighted wheel combination
  first forces separated center weights f(x)+g(y), then affine weights, giving
  a three-dimensional kernel and wheel rank p²-3. After two-family integration,
  wheel equations force T=F(x)+G(y), diagonal cycles impose sum F+sum G=0,
  and full solvability is equivalent to the single remaining period sum F=0.
  The seam proves this last functional is nonzero. Thus cycles plus wheels
  have rank p²+2 in the triangular stress space of dimension p²+3.
- The triangular-subgraph and full-graph counts remain distinct. At 29 the
  former has stress dimension 844, wheel rank 838, and cycle-plus-wheel rank
  843. All five complete family stress spaces have disjoint edge supports,
  so their direct sum has dimension 4,220. Its codimension in the full stress
  space of dimension 10,936 is 6,716. Cycles plus wheels alone span 4,215.
- The full 435-edge construction meets every full-graph straight cycle once
  and passes all 4,205 family wheel checks, yet fails the rigidity equations.
  This is a counterexample to sufficiency of those checks. It does not rule
  out another 435-edge support, including a possible optimum at 29.

## Independent finite-exact controls

Run `python3 verify_independent.py`. It uses only the standard library, reads
no producer artifacts, and writes its own JSON report beside itself.

1. Full-graph rigidity ranks at p=3,5,7,11,13,17,19,23,29 are checked through
   exact two-direction elimination; separate full matrices are reduced at
   p=3,5,7. At 29 the output is rank 1,679, stress dimension 10,936, and
   10,501 independent conditions beyond disjoint cycles.
2. The earlier q=3 and q=5 positive section witnesses are checked directly
   from their lifted coordinates. Every cycle has exactly one error -2;
   the affine carry equation and orientation invariance are also checked.
   The q=5 correction coordinates are copied explicitly from the preceding
   checkpoint, rather than imported from its verifier.
3. Seam controls at p=5,11,17,23,29 check all cycles, integer wheel support
   identities, actual carry-wheel syndromes, the residual identity, and both
   period sums. Exact reduced coefficient and augmented ranks differ by one.
4. The full p=29 construction is independently regenerated with the same
   deterministic family representatives. It has 12,615 edges, 435 selected
   edges, 435 valid cycle checks, and 4,205 valid wheel support/carry checks.
   Its entire reduced system has rank 55 and augmented rank 56.
5. Wheel and cycle-span ranks are directly computed at p=5,7,11,17,29,
   including exact wheel equilibrium and all three affine dependencies.
   At 29 these give 838 and 843. The p=7 check is an extra algebraic
   triangular-framework control, outside the addendum's anisotropic
   p=2 modulo 3 scope; no field-norm conclusion is inferred from it.

The producer's final standard-library verifier was also rerun successfully.
Its output is recorded separately as `producer-verified.out`.

## Scope and publication

The final section explicitly retains prime p>=5, p=2 modulo 3 for its
anisotropic triangular setting. The local wheel coefficient was written
explicitly to eliminate a difference-operator sign ambiguity; c and R were
defined explicitly at their first use. There is no remaining mathematical
blocker. The universal arguments are proofs; these finite controls are
regression checks, not replacements for those proofs.

Suggested independent publication subset: this report,
`verify_independent.py`, `independent-results.json`, and
`verify_independent.out`. No Git files were changed by this review.
