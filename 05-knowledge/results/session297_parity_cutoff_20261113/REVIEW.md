# Independent structural review

**Pass, with scope notes.** Reviewed CYCLIC, CUTOFF, and BRIDGES, and
replayed the direct chord-mask census at orders 3–8. No numerical result
here is used to infer an all-order rule.

- **Cyclic origin.** An odd-order rotation subgroup acts freely on
  unordered chords because a stabilizing nonidentity translation would
  interchange their endpoints. A reflection-fixed orbit contains
  exactly one fixed chord (solve the involution's affine equation over
  the odd rotation subgroup). With reflection offset -1, there is one
  fixed cycle edge at odd order and two at even order; subtracting
  from the swapped endpoint pairs leaves exactly floor((n-3)/2)
  fixed chords. An invariant chord's orientation is genuinely free:
  endpoint swap and arc reversal cancel. The cyclic marked object
  must not be confused with the unmarked converse orbit.
- **Blue incidence.** Independent adjacency/permutation code checked
  all labelled framed tilings through order 5 and 32 seeded samples at
  each of orders 6, 7, 8. It distinguished literal presentations from
  vertex-name permutations. In all even cases the half-rotation of each
  blue closed frame was blue and gave a different presentation. In all
  tested cases the blue-closed count modulo 2 was 0 at even order and
  was the Hamilton-cycle count modulo 2 on self-converse nodes at odd
  order. The structural proof, not these samples, yields the parity.
- **Cutoff.** Modulus-sized approximation must be strict: equality in
  the error bound still permits an adjacent lift. The odd composite
  control is multiplicatively closed, congruence-dense on the odd coset,
  and has every even difference by CRT, yet has all odd primes missing.
  The physical C2 case is inherited from Session199, and the one-clock
  regularity objection must remain scoped to the fixed-prefix strategy
  of Session197. The rejecting-language construction uses canonical
  positive digit words, so an accepting lasso cannot create spurious
  leading-zero words or finitely many integer values. Cofiniteness of
  a physical subfamily suffices for global coverage; a lasso of holes in
  that subfamily is not a lasso of global tournament omissions.
- **Prime wheels.** At any odd prime power the forbidden unit classes
  for the second summand coincide precisely when the prime divides the
  target. At 2^a the x-variable is determined by the y-variable, giving
  one x for each odd y exactly on the asserted parity. To justify
  Dirichlet density even if a target is nonunit, choose two admissible
  summand classes in each local factor first, combine by CRT, then
  independently choose arbitrarily large primes in the two unit classes.
  Equality of wheel counts never asserts that 2q is prime.
- **Gap and curve provenance.** Recovered the exact {4,11} inverse
  Collatz cycle, the two LRC CRT pairs, the normalized tournament
  {3,10} gaps, and the square-sum 23/25 counts in the cited source
  paths. Inherited ball-volume recurrence and MRB interval root are
  correctly scoped. The moving-radius proof uses strict Gamma
  log-convexity for monotonicity and
  Gamma(x+1/2)/Gamma(x)=sqrt(x)(1-1/(8x)+O(x^-2)) with x=(n+2)/2;
  squaring gives the stated n+3/2 expansion. The only ties are
  consecutive dimensions and radius exactly at a crossover.

No scalar gap cutoff or Chan–Pak proof follows from this review.
