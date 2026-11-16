# A second fresh ternary digit by changing the continuation

**Positive result.** The 163-path seed of the [first splice](../ternary_deposit/README.md) has a 493-path $+4$ child. Re-encode that scalar count with a different rooted tournament. It can now go to **1479, 1481, or 1483**, with one and the same new continuation $(E,O,g)=(1107,1122,15)$. The original child's continuation was $(252,267,15)$; the re-encoded parent has $(369,374,5)$; the three new outputs triple this **changed** continuation. They do not triple the original continuation. This is a second exact digit at a matched target, not a general second-round rule.

## Exact families and graphs

The terminal rows below are outgoing-neighbor masks, vertex 0 in the low bit. The root is the last vertex and the marked slot is given in the table. A multiplier $m>1$ means prepend an ordinal strong component of count $m$ (a triangle for 3, the unique strong four-vertex type for 5). That prefix is forced into the complementary path for $E,O$, so it multiplies all four coordinates. The remaining rows are strong without a prefix.

| Base output | $m$ | Mark | Terminal rows |
|---:|---:|---:|---|
| 1479 | 3 | 8 | 510,1004,1016,496,482,448,896,256,512,185 |
| 1481 | 1 | 10 | 2014,4092,4056,4048,4064,4045,1920,1792,3584,1024,2048,705 |
| 1483 | 1 | 10 | 2042,4092,4089,4016,4064,1984,1928,1792,3584,1024,2048,737 |
| 1469 | 1 | 9 | 2014,1020,984,976,2016,973,896,1792,512,1024,366 |
| 715 | 5 | 7 | 254,252,504,240,224,448,128,256,91 |
| 717 | 1 | 9 | 1022,2044,2040,2032,1888,960,1920,784,512,1024,417 |
| 719 | 1 | 9 | 2038,1020,1016,2033,2016,960,896,1792,512,1024,358 |

In each, replace the mark with any nonempty tournament $M$. Let $u=H(M)$ and $R=P_M(2)$, the ordered optional two-path cover count. The three upper outputs have profiles
$$
(1107R-735u+du,\ 1107R/2,\ 1107R/2+15u,\ 15u),
\qquad d=0,2,4.
$$
At the same continuation the extra 1469 row gives offset $-10u$. The terminal component of the first row in this top band is itself an exact **new source**:
$$
P=(369R-245u,\ 369R/2,\ 369R/2+5u,\ 5u).
$$
Thus the three profiles are $3P+(du,0,0,0)$, with no loss of the absolute anchor in this new family. For a transitive block of length $s\ge1$, $u=1,R=2^s$, so the units stay absolute 2 and 4 at unbounded sizes.

The lower band similarly has count $545R-375u+du$, $d=0,2,4$, with continuation $(545R/2,545R/2+5u,5u)$. At base size, 715,717,719 cross the old digit-block boundary. Since $717=3\cdot239$, this is a signed $-2,0,+2$ window at the first splice's 239 branch, not its requested $0,+2,+4$ window. It is not a precise tripling of a predecessor continuation, since the transitive band's $g$ is 5.

**Proof.** In every terminal template the mark's unique outward neighbor is the root. Each nonterminal module run consumes this root, so a Hamilton path can use at most two runs. In a cover whose first path starts at the root, any module run in either path has to be terminal. The only possible patterns are one run in one path or a split, one in each. After deletion of the root there can only be one run. The checker exhausts the finite ordinary-vertex/run skeletons, then uses one-run weight $u$ and two-run weight $R-2u$. For split runs the signed first-path weight is $-(1+(-1)^{|M|})u$ by the [alternating path-slot identity](../triple_splice/README.md#arbitrary-modules-exact-digit-cost). A run wholly in the first path has weight $(-1)^{|M|}u$ and one wholly in the second has weight $u$. The coefficient identities are checked symbolically for each row. A separate direct subset DP verifies literal transitive and nontransitive samples. Run **node check_bands.js**.

## What rewires and what does not align

For the literal 493 source, keeping its root matched to the output root, each of the three displayed terminal constructions can be reached with exactly four old-arc flips for an optimal old-to-new vertex assignment, then adding the remaining vertices. For the 1479 child add one terminal vertex and the three prefix vertices; for 1481 and 1483 add three vertices. The prefix arcs all go forward and do not affect the old-arc cost. This scope is the three selected terminal witnesses with fixed root matching, not a minimum over all tournaments at each count. The exact old-to-new assignments and arc reversals are in the output of **c++ -O3 -std=c++17 edit.cpp -o /tmp/recontact; /tmp/recontact**.

The new ready source $369\cdot2^s-245$ and the old $+4$ child $252\cdot2^t-11$ coincide only at $s=t=1$. Equality would require
$$
369\cdot2^s-252\cdot2^t=234.
$$
The right side has valuation one and the second term on the left is divisible by eight. Hence $s=1$, which forces $t=1$. Thus the new band is infinite **and** the 493 handoff is exact, but the handoff is not infinite along the old row. Neither the first +4 choice nor a second dyadic inflation can be silently treated as a scalar-only move.

See [DYADIC.md](DYADIC.md) for the relation to signed Collatz, shifted doubling rows, and the archive's fixed-root response obstructions. See [SMALL_PROFILES.md](SMALL_PROFILES.md) for the all-order boundary in six low deletion layers of the two literal nine-count windows.

**Next-height follow-up.** The original source of count 331 has a 997-path +4 child. An independent changed continuation reopens it at 2991,2993,2995, with common (2897,2912,15). See [HEIGHT_TWO.md](HEIGHT_TWO.md) for the exact graphs and a third infinite replacement band. The two old-row matches do not prove alignment at all heights.
