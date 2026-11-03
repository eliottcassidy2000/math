# Independent audit

Status: PASS, 2026-11-03.

The proof that equal traces force x=y or xy=1 is correct over F_7417. The inversion fixed points are 1 and -1, so the trace image has 1237 elements. The first invisible pair occurs at distance 618 and thus odd length 1237. Even lengths cannot have an invisible pair, since their exponent gaps are odd whereas every exponent giving -1 is 1236 modulo 2472. This is distinct from collisions between full words, which already occur in the displayed 7-digit census.

An independent standard-library program enumerates every decimal half-word directly through 13 digits, constructs its mirrored integer and tallies residues. All six full census summaries agree with the dynamic-programming certificate. This is an independent computational path rather than a second DP implementation. The program independently enumerates each unit's multiplicative order for moduli 240, 1001 and 2222 and confirms group orders and exponents. The abstract product decompositions in README follow from CRT and the standard unit groups of the stated prime powers.

The decimal subgroup has order 2472, index 3 in the cyclic field unit group of order 7416. The associated dihedral action has 4944 elements. AGL(1,7417) has 55,004,472 elements. These are correctly distinguished from the trace quotient set.

The automaton's leading-digit restriction and center convolution are correct for m>=1 and arbitrary modulus. Its displayed census does not cover one-digit numbers; that would require a nonzero central digit. The README only computes m=1 through 6.

Banks–Hart–Sakata's actual primary source, https://arxiv.org/html/math/0405056v1 , introduction and Corollary 4.4, confirms the stated p>10 and order>=3 sqrt(p) criterion and the scope of asymptotic equidistribution of all palindromes. No sparse-family or prime-density conclusion follows.

Reproduction: python3 check.py and python3 -O check.py; python3 independent_check.py and python3 -O independent_check.py. Each normal/optimized output pair agrees byte-for-byte. Saved receipts are audit_normal.json, audit_optimized.json, independent_results.json, and independent_results_optimized.json.
