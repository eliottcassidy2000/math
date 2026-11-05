# Independent audit: physical retained-root digit blocks

Verdict: PASS.

- DIGIT_BLOCK.md SHA256: aec3252026cf755687fe20767830cfebc2a888cddd292b21f8c3a236d8d99fc4
- check_digit_block.py SHA256: bd6d1076ea80c5e8586979a0d691a82dc31411c06fc978904047af3c4c92c21b
- Both supplied receipts SHA256: 50b861fdcd769c6d64789d53dc599e58b6195c9fbad11c9515b404cb9678dfde

The chronological convention is consistent: the earliest letter is the least
significant bit of A, whereas the accumulated U subtraction is the fixed-length
bit reversal R. The U/E/d formulas follow directly from the physical recurrences.
For D, the earlier-letter interactions and adjacent E correction give exactly
the displayed C(A). The recursion for f(A)=g(A+1) expands into the complementary
bit-pair sum. Adding those sums gives A-p, which proves the D formula at all
lengths. Leading-zero cases and k=1 are included. The zero-length block is
correctly handled separately.

The complement identity has the correct sign. Pair terms sum to A-R and adjacent
terms telescope to t-p. The final wording correctly subtracts the complemented
expansion from the original. This arithmetic cancellation alone does not assert
physical closure under subtraction.

Both fixed-left C2 readouts are correct. Substituting the full source formula
produces the displayed constant terms -5d and -7d and all other coefficients.
The reviewed complete O/Z/J grammar justifies splitting every nontrivial source
word at J letters, allowing zero-length blocks between consecutive J moves.
Thus finite tuples of the stated integer parameters are unconditionally physical,
and they cover every nontrivial source. This is a parameterization of the fixed-
left subclass, not a proof that its scalar range is cofinite. The singleton is
properly excluded from physical Z application.

Isolated normal and optimized checker runs reproduce both frozen receipts and
all four release-file hashes. They test 8,190 complements, 32,760 linear basis
blocks and 17,892 physical-prefix blocks. The algebraic proof, not this finite
range, supports the all-length conclusion.

An additional no-import check independently constructs physical source graphs,
counts their complete four-coordinate profiles, and constructs actual N/M grafts.
It uses the bit-pair expansion rather than the recursive g evaluator. All182
retained-root blocks over13 prefix graphs agree, and all364 fixed-left grafts
have exactly the predicted Hamilton count and C=2, up to order11. See
../digit-independent.py and ../DIGIT_INDEPENDENT.json.

No conductor, cofinite coverage theorem, or unrestricted source-membership
assumption has been added.
