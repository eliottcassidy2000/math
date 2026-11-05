# Independent audit: square top amplitudes for marked self-converse rank two

Status: **PASS, no correction required.** No frozen source files were changed.

## Frozen sources

Files in `session197-prime-atoms/`:

- `SYMMETRIC_SQUARE.md`: SHA256 `933f74d20c1fa3aee6654274d0517bc849aede53242277ae4044fe132c26d62e`.
- `verify_symmetric_square.py`: `1453fddaed72057c8cf62195088b4d416ab584dd999bb9137c16ea97da1d2ae7`.
- `verify_symmetric_square.cpp`: `a24a0666a6b6755a5666f86864ea46a0e80168633ed338c720ce0663967589c0`.
- `receipt-symmetric-square.json`: `3ce8e75314c2c5814363c9b0aae6b1007a3e40f45c5b603d146ad5bc2e03ce66`.

The named root-run, symmetric-bound and even-self-converse dependencies were previously reviewed; their relevant statements were checked against the proof here.

## All-order proof review

The root-fixed anti-involution pairs O with I and makes the O-to-I matrix symmetric, including its diagonal. The classification is exhaustive: two disjoint off-diagonal edges force matching four, a triangle forces matching three, and a pairwise-intersecting triangle-free edge family is a star. Additional loops are restricted exactly as stated. The edge-only, one-loop, full-block and two-isolated-loop cases are all covered, with no missing self-paired crossing arcs.

In a star, the two internal root fragments must separate the center vertices u and tau(u). Neither fragment can use both, because then the other has no crossing. After deleting those two vertices all residual cross arcs point I' to O'. This forces the stated prefix-center-tail forms, an all-I' outer prefix and an all-O' outer suffix. The two three-path half-covers are related by reversal followed by the anti-involution. Their choices are independent, and either order of the center-bearing fragments is allowed. This proves exactly 2q^2, including all empty-prefix and empty-suffix cases.

For two isolated loops and the full block, each internal fragment must contain precisely one active O vertex and one active I vertex. They are adjacent at the crossing; any intervening residual vertex would require an unavailable extra crossing. The same half-cover factorization applies. There is one matching in the diagonal case and two in the full block, yielding respectively 2q^2 and 4q^2 after the two fragment orders. The q>0 constructions are valid even when the residual set is empty.

For a negative-unit ternary clock, the even-clone congruence forces A=2 modulo4. Thus only 2q^2 with q odd remains. The prior lower bound 78 leaves only A98 and A162 below242. Removing outside singleton SCCs is justified by the gcd of the first two clone counts. The invariant strong root component has odd order and H<=485. The exact comparison 5^12>485^3 excludes every order at least13, so the complete order11 check closes the remaining gap. No attainment at242 is claimed.

The root-unfixed hostile and the separate positive-unit example correctly distinguish the theorem's hypotheses and sign. They do not contradict the square factorization or imply a broader unit-sign obstruction.

## Source replays

Normal and Python-optimized source runs reproduced the frozen receipt byte for byte. This rechecks all symmetric binary matrices through half-order five, the 4,146 direct root-run versus half-cover controls through order nine, 12,438 exact clone counts, the root-unfixed hostile and the positive-unit scope example.

## Independent order-eleven certificate

`check_symmetric_square.cpp` is a separate implementation. It enumerates every symmetric binary 5-by-5 matrix and obtains its matching number by subset-state matching, rather than generating only the classified templates or using augmenting paths. Exactly 160 have matching number two. It then enumerates every one of the 1,024 O-side tournaments, producing the complete 163,840 physical bank.

Half-cover counts are computed by literal directed-path enumeration on O and direct assignment of the remaining vertices to three labelled slots. This uses neither the source subset endpoint table nor its nested subset-cover summation. The independent candidate counts are exactly 2,640 for q7 and 5,600 for q9.

For all 8,240 candidates, the complete physical tournament is constructed with the root first, and its anti-involution is checked explicitly. H is counted with odd-cycle collections rather than either Hamilton endpoint recurrence. Counts greater than485 are safely saturated using positive cycle-packing contributions; the target values293 and485 remain exact. There are no first-equation survivors. Optimized and unoptimized library-checking builds produce identical receipts.

Independent checker SHA256: `fa963849027e47c32102006458896eda70bfd4783d386bc8d58aef9ed2b6dae3`.
Independent receipt `square-independent.txt` SHA256: `6cc2cb2142ad51d4ca14743ab9bc538cbed0ed4d692d7bf43c140ee51ef393ff`.

Both the structural square theorem and the all-order A>=242 consequence pass at the frozen hashes above.
