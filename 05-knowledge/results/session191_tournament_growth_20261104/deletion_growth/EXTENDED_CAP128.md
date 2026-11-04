# Extend the strong-deletion increment atlas through order twelve

**FINITE-EXACT; independent replay passed.** This is a continuation of
the separately audited cap44 packet, not a dependency of that packet.

Use the exact same marked-OCF monotonicity and joint strong-deletion
completeness theorem in PROOF.md, changing only the even cap
from 44 to 128 and continuing to order 12. No H cutoff is imposed.

## Restricted strong-deletion pairs

The complete rooted bank sizes are

| order | rooted strong-deletion pairs with Delta<=128 | minimum |
|---:|---:|---:|
| 4 | 2 | 2 |
| 5 | 14 | 4 |
| 6 | 140 | 6 |
| 7 | 1959 | 10 |
| 8 | 5142 | 16 |
| 9 | 5121 | 26 |
| 10 | 2996 | 46 |
| 11 | 705 | 78 |
| 12 | 0 | none |

The step to order 11 tests 2,996*1,022=3,061,912 nonconstant ears, of which
3,055,920 preserve strongness after deleting the root. Exactly 3,548 retained
presentations quotient to 705 rooted classes. The minimum is 78, attained
with H=225 and H(T-0)=147.

The step to order 12 tests705*2,046=1,442,430 nonconstant ears, of which
1,441,020 preserve strongness after deleting the root. None has Delta<=128.
The minimum among examined ears is 130, attained with H=375 and H(T-0)=245.

Thus, the exact minimum over strong deletion
pairs at orders11 and12 is respectively78 and130, and every such pair at
all orders>=12 satisfies Delta>=130. The all-order exclusion follows from
the proved descent, not an extrapolation of the observed minimum.

Witness outgoing rows, with marked root 0, are:

    n11: (1038,4,8,1042,1159,1055,1087,1391,1663,1279,6),
    n12: (2062,4,8,2066,2183,2079,2111,2671,2303,3455,2559,6).

The corresponding interval-chain witnesses are IC(2,3,2,3) and
IC(3,2,3,3), deleting the middle vertex of the internal triangle. Their
counts are independently determined by the all-order chain formula.

## All roots, including weak root-deletions

The independent-from-the-restricted-bank all-root pipeline finishes
through order 11, with class counts at orders3,...,11

    1,4,24,194,2392,6242,5962,3168,705.

Its 705-class order 11 bank is byte-identical to the restricted order 11 bank.
At order 12 the restricted run already counts 1,441,020 ears. An independent complementary count checks precisely the 1,410 remaining ears, those whose
root-deletion is weak. None has Delta<=128; their minimum is 250, with
H=375 and H(T-0)=125. These two disjoint cases cover all 1,442,430 ears.

Therefore the completed computation, independently validated,
strengthens the result to **every vertex of every strong tournament**:

* exact minima for orders3,...,12 are
  2,2,4,6,10,16,26,46,78,130;
* at every order>=12, Delta(T,v)>=130 for every v.

The invariant-preserving quotient is still rooted isomorphism. No H cutoff,
ratio conjecture or claim that every old vertex stays strongly deletable is
used. The all-order implication is rooted strong descent from the allroots
addendum.

To reproduce both complete enumerations using the scripts accompanying this note:

    g++ -O3 -std=c++17 rooted_ears.cpp -o rooted_ears
    python3 run_rooted.py 128 12
    python3 run_allroots.py 128 12

The programs generate their own banks and logs. The independent replay used a separate implementation and separately checked the weak-root-deletion ears.
