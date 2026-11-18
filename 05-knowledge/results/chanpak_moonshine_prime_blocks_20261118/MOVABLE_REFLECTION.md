# Let the reflected arm be a block

This is a companion to the [genuine dihedral family](DIHEDRAL_HEIGHT.md). Here reflection exchanges two frames, not necessarily vertices of one fixed host. The new freedom is to let the entire reflected arm be a nonempty tournament instead of a singleton.

Use the four-vertex crossed quotient A=(10,4,9,2) and switch vertex 3 to obtain A'=(2,12,1,5). Insert L,R,W at 0,2,3. Vertex 1 stays a singleton. In each frame the three block positions are transitive after deleting vertex 1. The full switch reverses every arc between W and its exterior and no internal arc.

Let F_U(q) count ordered covers by q optional paths of U. For q≥1 the change A minus A' is exactly
```
q sum_{s=±1} (F_W(q+s)−F_W(q))
   [F_L(q)F_R(q+s)−F_R(q)F_L(q+s)].
```
To see it, use the odd-cycle packing formula. Every directed cycle using multiple blocks or meeting vertex 1 uses the singleton 1, unless it is wholly inside a block. The cyclic color words after 1 are RL, RW, RLW in A and RL, WL, WRL in A'. The contributions from RL and the three-block word agree since each participating block uses one run. For the two remaining words, pick nonempty internal supports of total even size, to make the cycle odd after including 1. The unrestricted and signed nonempty sums on a block are precisely F_U(q+1)−F_U(q) and F_U(q−1)−F_U(q), by slot convolution and F_S(−1)=(−1)^|S|H(S). The even selector contributes one half, canceled by the 2q weight of the crossing odd cycle. This proves the formula without assuming the blocks themselves are transitive.

Write h,B,C for reads at 1,2,3, u=B−h, x=C/B and y=h/B. Exact neutrality at q=2 is the slope match
```
rho = (B_W−h_W)/(C_W−B_W) = (x_R−x_L)/(y_R−y_L)
```
when the latter difference has nonzero denominator. The direct cross-multiplied formula applies in every case. At q=1 the count on the first side is
```
H_A = 2 h_L h_R h_W + B_W u_L u_R + u_W h_L u_R,
```
and the second swaps L and R. For example the expression follows by enumerating compressed words in which each block has one or two runs, then setting the two-run read to B−2h. The difference A' minus A is u_W(B_L h_R−h_L B_R).

## Choose the arm from the absolute target

Fix the two flanks with positive slope rho=t/s in lowest terms. A neutral arm must have three linked readings
```
(h,B,C)_W = (h,h+t z,h+(t+s)z), z≥1.
```
Consequently a prescribed target N on the first side imposes the positive linear equation
```
N = (u_L u_R + 2 h_L h_R) h + t B_L u_R z.
```
Its finite integer solution progression gives an exact list of profile requests, rather than an interval around N. The issue is to realize one of those linked triples by an actual W. No fixed flanks are claimed to do so for all integers.

**Physical checks.** A two-vertex transitive W has profile (1,4,9), so t/s=3/5. Two nine-vertex flanks in the replay yield one switch with exact ordinary counts 1,729,645,521 and 1,734,383,121, and common height 630,250,594,716. These counts are composite. A different pair of even-order flanks with the same arm gives H=11,010,786,697,063 and 11,015,945,813,863 at height 112,463,595,873,847,652. Neither is divisible by 3 or has a common divisor with the height. Apply one ordinal singleton to the first host and the known paired slope-6 readers: their output is the prime 1,349,574,161,272,868,887 and its exact sevenfold partner 9,447,019,128,910,082,209. The primality certificate uses the complete factorization of p−1 and base 3.

The arm can genuinely be nontransitive. The strong four-vertex W=(6,12,8,1) reads (5,32,117) and selects t/s=27/85. The replay verifies one compatible flank pair and exact counts 96,833,175,110,435 and 96,875,688,057,635 at common height 1,578,628,339,040,509,540. This is a literal larger reflected arm, not a prime example: these endpoints are composite. The replay checks 129 smaller literal slot counts against the formula and the compressed transfer, plus the larger symbolic values.

Replay: `python3 movable_arm.py` alongside the earlier packet files.
