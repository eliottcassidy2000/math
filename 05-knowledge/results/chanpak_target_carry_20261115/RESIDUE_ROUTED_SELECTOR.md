# Ask for one terminal residue, not a menu of all phases

The [reversed carry selector](REVERSED_CARRY_SELECTOR.md) minimizes
the model distance from a given target on an address interval.
This upgrade optimizes **only over physical pairs with
F(x,s−x)≡N mod 2M**, independently of an all-syndrome certificate
or complementary alignment. It can reject, approximate, or
certify equality on the interval.

## Read the checksum from the same end as the carries

Scan the digits of x from the low end. After j bits, put
R_j=sum_{i<j}(1−b_i)2^(j−1−i), the reversed low zero mask,
and let ell be the most recent low bit. When adding the new
higher bit b_j, the old lower part keeps its original checksum.
Its new contribution from bit b_j equals
```
p_new−p_old = b_j [R_j+1−ell]  (j>=1);
R_(j+1) = 2R_j+1−b_j.
```
For j=0 the increment is zero, equivalently initialize ell=1
and R_0=0. This follows from the two-site formula for p:
the old zeros at distance d from the new higher 1 have
weight 2^(d−1) except that distance 1 has weight 2 instead of
1. This is why the last observed digit is needed.

Keep the four channels modulo 2M: x's already accumulated address,
the two reversed zero masks, and the combined unhalved checksum.
Add a bit b to the address with weight 2^j; keep each tape's
most recent bit. There is no need to store y's address as another
coordinate because y=s−x at the terminal check. At the endpoint
compute the polynomial part B+Qx(s−x)+Tx+S(s−x) and add 45 times
the accumulated checksum, all modulo 2M. Accept **only the
requested residue** N. No other residue gets a target copy.

## Exact functional graph and useful outcomes

Add the preceding channels to the [54-state carry and comparison
graph](REVERSED_CARRY_SELECTOR.md), together with the two last
bits, and do the same predecessor/successor optimization.
There are at most 216·(2M)^4 states per layer. The input length
is O(log s), and all transition arithmetic is modular plus the
bounded-size integer carry comparisons. For fixed M this is
linear in the number of columns in a unit-cost state count.
For growing M it is a finite exact procedure, not a near-linear
bound in log M.

Choosing the first finishing edge produces one successor on the
viable graph and a physical (x,y) at the end. If neither search
has a path, the entire interval has no exact witness. Otherwise
let R be the best model discrepancy and E its certified envelope
from the reversed-carry note. Then:

- If R>E, the entire interval has no exact witness, even though
  one or more congruent pairs may be present.
- If R+E<2M, the returned physical pair is exact: its integer
  difference from N is a multiple of 2M of absolute value less
  than 2M.
- Otherwise only the bound R+E on the achieved error is proved.

This does not prove that the target has a favorable interval
or modulus. In particular, choosing M merely larger than the
envelope need not create any finishing path. The benefit is to
keep the one-edge choice, exact congruence, and integer envelope
inside a single inverse calculation rather than separately
programming all outgoing phases.

## Replay

`replay_residue_selector.py` compares the residue-constrained
nearest carry with brute force on small prescribed address windows,
checks empty-basin rejection, and verifies each exact-equality
certificate by evaluating the actual count formula.
