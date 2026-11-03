# One-apex spectrum computation and reproducible count

## Status and prior work

The global classification remains open: the conjecture is that every positive odd integer other than 7 and 21 is a tournament Hamiltonian-path count. The omission of 7 and 21 is recorded as proved in `math/01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`. THM-4749 reports a finite strong-spectrum certification below 80,265 and attributes an earlier 2016 inventory of all odd values below 80,557 except the two gaps to bof and Royle. Neither finite record proves asymptotic density or completeness. See also THM-4748 for one-apex and other low-feedback constructions. The calculations below are a restricted-family experiment, not a new global completeness bound.

## A direct one-apex identity

Let vertices 0,...,m−1 form a transitive chain, and let an apex x beat precisely those chain vertices whose bit in mask M is one. Assume bit 0 is one and bit m−1 zero. Then the tournament is strongly connected and its number H of Hamiltonian paths is

    H = 2 + sum_i [i+1 < m, bit(i)=0, bit(i+1)=1]
          + sum_{0 <= j < i < m} [bit(j)=1, bit(i)=0] 2^(i-j-1).

Proof: Remove the apex from a path. The chain vertices partition into an increasing run A before x and an increasing run B after x. A run may be empty. If A is empty there is one path since x beats vertex 0; if B is empty there is one path since vertex m−1 beats x. When both are nonempty put i=max A and j=min B. If i<j then j=i+1 and the partition is unique, subject to i→x→j. If j<i, vertices below j must lie in A, vertices above i in B, and each of the i−j−1 intervening vertices independently chooses a run. The same edge conditions apply. This accounts for every path exactly once. Strong connectivity follows because x reaches every chain vertex via vertex 0 and every chain vertex reaches x via vertex m−1.

**Finite-decision corollary.** The single pair (j,i)=(0,m−1) contributes 2^(m−2), so H ≥ 2+2^(m−2). Thus a strong one-apex witness for target h has m ≤ 2+floor(log2(h−2)); to decide membership in this subfamily it suffices to check finitely many masks. In particular h=13 forces m≤5, and exhaustive checks show no such mask. This obstruction is only to this particular family, not to tournaments generally: 13 is a known tournament count. The bound also permits a direct finite search for any fixed h, with O(h log² h) elementary arithmetic operations using the displayed formula.

`check.py` independently tests the identity against Held–Karp subset dynamic programming for all 127 strong one-apex signatures with 2≤m≤8. `apex.cpp` enumerates all such signatures through m=22, retaining values ≤2,000,000, and writes the first witness of each value to `apex-record.csv`. Run `g++ -O3 apex.cpp -o /tmp/apex && /tmp/apex`. Counts of newly encountered values for m=2,...,22 are 1,1,2,4,9,12,28,52,106,198,411,780,1580,3060,6163,12273,24509,48990,98035,193961,107629. The final count is cut off by the bound. These signatures do not realize 13 at any order in the searched range, so this family by itself is far from a classification.

## Explicit two-apex example

Let chain vertices be 0,...,11 with smaller beating larger, let apex x beat precisely the vertices in bitmask 2663 and y beat precisely those in bitmask 1740, and let x→y. A subset-DP check yields H=80,269; the tournament is strongly connected. It is an independently checkable example, not a new range record: 80,269 is already included in the previously reported finite inventory. The original 16-vertex random-search witness had masks 14951 and 14028 for a 14-chain; its last two chain vertices are singleton strong components and may be removed without changing H, giving the 14-vertex witness above.

To verify the example execute `python3 check.py`; it includes an explicit second subset-DP calculation and strong-connectivity check. `twoapex.cpp` is the seeded random search that produced the longer witness; its exact seed, inputs, and iteration are recorded in the source. No claim of exhaustion is made for the two-apex search.
