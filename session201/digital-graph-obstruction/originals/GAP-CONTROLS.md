# Two existing gaps are not explained by the triple-prime screen

These are exact checks of existing targets, not a new occupancy census.
They use the previously audited c=3/c=5 factorization from
`../session201-twelve-factorization/PROOF.md`.

- 2647 is prime, but 3n-2=7939=17*467 and 5n-4=13231=101*131.
  Since n is prime, only the odd/odd factorization can apply. For c=3,
  the shifted divisors17 and467 request source sums6 and156. The source-sum
  fiber at6 is {3}; the fiber at156 is empty. For c=5, neither101 nor131
  is congruent to9 modulo10, so no admissible nontrivial factor pair exists.

- 2134117=643*3319, 3n-2=6402349=23*278363, and 5n-4=10670581 is prime.
  For c=3, the odd/odd pair already fails at source sum(23+1)/3=8, whose
  fiber is empty; no mixed factor has the required residue5 modulo6.
  For c=5, the prime shifted target excludes odd/odd indices. The only
  nontrivial mixed candidate is643*3319 with source sums643 and664;
  both fibers are empty. Neither c divides n, excluding the both-even case.

All source-fiber claims follow from the complete finite intervals proved in
the factorization packet, with literal g evaluation. The factorizations can
be verified by multiplication and trial division. These observations do not
change the separately established all-height gap status of either target.
