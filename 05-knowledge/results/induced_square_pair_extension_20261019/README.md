# Square-scaled induced extensions

2026-10-19. Proofs independently reviewed; exact controls pass. No external priority claim.

## Main result

Read [SCALED_PAIR_EXTENSION_AND_TWO_DEGENERATE.md](SCALED_PAIR_EXTENSION_AND_TWO_DEGENERATE.md) first. It proves universal induced two-neighbor extension after a common square scaling of any finite seed, induced realizability of every finite two-degenerate graph, and an explicit height below 4096n^4 for the exact one-subdivision of every n-vertex graph. The particular seed {1,2,3} shows why a universal three-neighbor version fails even after square scaling, using an inherited congruent-number obstruction.

The finite two-degenerate construction has an exponential bound on label bit length; it has no claimed polynomial height or countable limit.

## Proof and review inventory

- [SCALED_PAIR_EXTENSION_AND_TWO_DEGENERATE.md](SCALED_PAIR_EXTENSION_AND_TWO_DEGENERATE.md): principal proof
- [INDEPENDENT_REVIEW.md](INDEPENDENT_REVIEW.md) and [review-root/PROOF_REVIEW.md](review-root/PROOF_REVIEW.md): independent audits
- [ODD_PRIME_ONE_SUBDIVISION.md](ODD_PRIME_ONE_SUBDIVISION.md): separate prime-power construction
- [MOD7_MINIMAL_RESIDUE_ADDENDUM.md](MOD7_MINIMAL_RESIDUE_ADDENDUM.md): optional precisely scoped minimality fact
- [provenance/ONE_SUBDIVISION_DYADIC_GUARD.md](provenance/ONE_SUBDIVISION_DYADIC_GUARD.md): independently derived earlier construction; its weaker constant is superseded

## Reproduce from this directory

    python3 -B verify_one_subdivision.py > /tmp/one-subdivision.json
    python3 -B -O verify_one_subdivision.py > /tmp/one-subdivision-optimized.json
    cmp /tmp/one-subdivision.json /tmp/one-subdivision-optimized.json
    cmp /tmp/one-subdivision.json one_subdivision_results.json
    python3 -B verify_two_degenerate.py > /tmp/two-degenerate.json
    python3 -B -O verify_two_degenerate.py > /tmp/two-degenerate-optimized.json
    cmp /tmp/two-degenerate.json /tmp/two-degenerate-optimized.json
    cmp /tmp/two-degenerate.json two_degenerate_results.json
    python3 -B verify_prime_family.py > /tmp/prime-family.json
    python3 -B -O verify_prime_family.py > /tmp/prime-family-optimized.json
    cmp /tmp/prime-family.json /tmp/prime-family-optimized.json
    cmp /tmp/prime-family.json prime_family_results.json

The three programs record 7,991,806 explicit checks in total. SOURCE_PROVENANCE.json records the exact reviewed source hashes and publication-only presentation changes.
