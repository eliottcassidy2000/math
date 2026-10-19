# Sharp induced-book height

The square-sum graph on positive integer labels contains an induced book with `m` pages at minimal maximum label

`H(m) = exp((1/log(2) + o(1)) log(m) log log(m))`.

Equivalently, the maximum number of induced-book pages available among labels at most `N` has logarithm

`log b(N) ~ log(2) log(N) / log log(N)`.

The proof is elementary and includes an explicit construction for every `m`. The common-neighbor factorization is classical; the matching construction combines a divisor concentration argument with a modulus-3 independence condition. A consequence is that all finite 2-degenerate graphs cannot have a uniform polynomial label bound, even though every such graph is induced-realizable.

## Contents

- `BOOK_HEIGHT_ASYMPTOTIC.md`: full proof, exact finite parameters, scope, and primary background sources.
- `verify_book_height.py`: primary standard-library exact checker.
- `exact_results.json`, `exact_results_optimized.json`: identical normal and optimized receipts, **4,081,492** checks each.
- `review/independent_controls.py`: separately written exact checker.
- `review/independent_results.json`, `review/independent_results_optimized.json`: identical independent receipts, **105,750** checks each.
- `review/ARITHMETIC_REVIEW.md` and `review/ROOT_PROOF_REVIEW.md`: independent proof audits.
- `PUBLICATION_FILES.txt`, `FROZEN_SHA256.txt`: publication contents and integrity manifest.

## Reproduce

Run from this directory:

```sh
python3 verify_book_height.py > /tmp/book-primary.json
python3 -O verify_book_height.py > /tmp/book-primary-optimized.json
cmp exact_results.json /tmp/book-primary.json
cmp /tmp/book-primary.json /tmp/book-primary-optimized.json
python3 review/independent_controls.py > /tmp/book-independent.json
python3 -O review/independent_controls.py > /tmp/book-independent-optimized.json
cmp review/independent_results.json /tmp/book-independent.json
cmp /tmp/book-independent.json /tmp/book-independent-optimized.json
sha256sum -c FROZEN_SHA256.txt
```

The combined suites execute **4,187,242** exact checks per normal or optimized run. Both scripts use explicit checks rather than disabled-by-optimization assertions. The primary suite completely enumerates divisors for `k=4,...,18`, verifies common-neighbor formulas on 16,340 anchor pairs, checks all unrestricted parity capacities for differences through 5,000, samples larger parameters through `k=256`, and fully checks induced adjacency in 40 books with up to 1,024 pages. It includes a concrete eight-page certificate. These are finite implementation checks; the written arguments prove the unbounded conclusions.

No previous frozen packet is modified by this addendum. There is no claim that the sharp book height is the optimal universal height for the entire 2-degenerate graph class.
