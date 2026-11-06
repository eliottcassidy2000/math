# Independent continuation audit: PASS

Verified all nine order-nine witnesses using a separately written memoized endpoint recurrence (`independent_audit.py`). Checked full strong connectivity and checked each 36-bit adjacency code restricts to the stated order-eight base and has the stated incident pattern. Witnesses realize exactly611,615,617,619,623,625,635,647,655.

Reran the author's all-9!-permutations checker; output matches frozen check.out byte-for-byte. Recompiled and reran fill.cpp on inputs.txt; output matches fill.out byte-for-byte, finding all nine values after2917 extensions. This targeted witness search is not a full order-nine support or class census and must not be described as one.

The inference m(h)=9 for these nine counts is rigorous: the independently audited exhaustive order-eight support excludes them, and source padding implies all smaller supports are subsets of it. Combined with first occurrences in the independently audited orders1..8 census, this yields the complete329-entry minimum-order table for attainable h in[1,661]. The table's output was compared entrywise against the independently recomputed mapping. Even h and7,21 are impossible by inherited theorems; h1 has minimum positive order1. Every listed inverse order-set is its upward ray by source padding.

The cyclic-fan formula is correct for k>=1: three contiguous-chain placements plus2^k−2 proper nonempty splits give2^k+1 paths. This is the classical nearly-transitive tournament NT_(k+2), ordered as y,chain,x with sole backward arc x→y, and should be credited to the inherited Moon/Busch formula recorded in THM-4744. It supplies an upper bound, not an exact general minimum; k5 is already strict because m33=6<7.

No new global coverage or closed-form claim is supported or made. No publication performed by the reviewer.

## Expanded proof addendum

Reviewed final PROOF.md's additional deductions. PASS. The prime-power recurrence for positive exponent is the exact strong-atom recurrence specialized to powers of one prime. Independently checked the seed minima m(3^j)=3,5,6,7,8 for j1..5 and equality with strong minima. Ordered joins give the8q+d_r upper bound. For any attained h>1 all powers are attained, so Fekete applies to a finite nonnegative subadditive sequence; the existence and infimum characterization of its limit are valid. Positivity or attainment is correctly not asserted.

Independently compared all329 full first occurrences with strong first occurrences from the complete n<=8 census and the nine verified strong n9 witnesses; equality holds entrywise. This finite evidence does not establish the global strong-first-occurrence question, correctly left open. The singleton is strong under the packet's convention; the global question is explicitly h>1.
