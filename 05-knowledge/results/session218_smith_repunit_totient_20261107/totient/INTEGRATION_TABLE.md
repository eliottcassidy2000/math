# Even nontotients through 76: exact integration

Status: FINITE-EXACT. The universe is every even integer 2 through 76, with complete inverse-totient enumeration filtering empty fibers. F counts prime factors with multiplicity. Base digits are integer arrays, not decimal strings. Repair a means phi(x)=2^a n; every smaller a was checked, so a found within 1..8 is the exact first positive repair. A failure through 8 would only be bounded-not-found. No unbounded existence claim is made.

| n | factorization | F(n) | Smith base n−F+1 | digits | factor-digit weight | prime half? | phi(n) | first repair a | preimage x |
|---|---|---|---|---|---|---|---|---|---|
| 14 | 2 * 7 | 9 | 6 | [2, 2] | 4 | True | 6 | 1 | 29 |
| 26 | 2 * 13 | 15 | 12 | [2, 2] | 4 | True | 12 | 1 | 53 |
| 34 | 2 * 17 | 19 | 16 | [2, 2] | 4 | True | 16 | 2 | 137 |
| 38 | 2 * 19 | 21 | 18 | [2, 2] | 4 | True | 18 | 5 | 1217 |
| 50 | 2 * 5^2 | 12 | 39 | [1, 11] | 12 | False | 20 | 1 | 101 |
| 62 | 2 * 31 | 33 | 30 | [2, 2] | 4 | True | 30 | 7 | 7937 |
| 68 | 2^2 * 17 | 21 | 48 | [1, 20] | 21 | False | 32 | 1 | 137 |
| 74 | 2 * 37 | 39 | 36 | [2, 2] | 4 | True | 36 | 1 | 149 |
| 76 | 2^2 * 19 | 23 | 54 | [1, 22] | 23 | False | 36 | 4 | 1217 |

The moving-base Smith theorem of session218-smith/BASE_COMPLETENESS.md applies to every even nontotient: each is composite (2 is a totient and other even integers are composite), and the even excluded Smith numbers 6,8,12 are themselves totients. This does not characterize nontotients or imply a common-base Smith property. The phi(n) column concerns n as input; the nontotient filter concerns n as output.

Reproduce: python3 session218-scaling/integration_table.py. Uses exact prime-support inverse enumeration from verify.py, independently checked there against a complete totient sieve for targets <=200. Larger repair targets retain the same proven-complete source enumeration, with explicit preimages checked independently by Euler’s totient formula.
