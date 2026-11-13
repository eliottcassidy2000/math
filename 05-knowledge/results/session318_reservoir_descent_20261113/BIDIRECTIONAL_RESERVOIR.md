# A reservoir with two read directions

This calculation uses the same physical cyclic buffer Z as the pure-end
switch: w -> a -> TT_m -> w, with the chain b_1 -> ... -> b_m,
m>=1. Let the apex cut omit exactly the chain sites marked by binary
digits d_1,...,d_m (both w and a are in the out-cut). Write r=h_Z-s_Z
for the missing-start count, and E for the two-slot end-marked response
of the terminal atlas. Then exactly
```
r = sum_j d_j 2^(m-j)
E = sum_j d_j [2^(m-j+1) 3^(j-1) + 2^j - [j=m]].
```
Thus r is ordinary binary read from the start, while the dominant part
of E is a redundant base-3/2 read from the end. This is not the same
as a binary and ternary independent digit bank.

**Proof.** A full path starting at b_j has all earlier chain vertices
after the forced pair w,a, and any subset of later ones before that
pair. Its count is 2^(m-j). Starts at a and w have one each, so h_Z
is 2^m+1. To count a distinguished nonempty slot ending at b_j, let
l=j−1 and q=m−j. If that slot uses neither w nor a, its chosen chain
vertices other than b_j can be any subset of the l earlier vertices.
The other slot contains the two special vertices and has 2^k+1 paths
when its chain size k is positive, but 1 path when k=0. Its summed
contribution is 2^q 3^l + 2^l − [q=0]. If the distinguished slot
uses only a, it contributes 2^l; if only w, it contributes zero.
If it uses both, w must be followed by a. Each selected earlier
chain vertex can go before w or after a, and each later chain vertex
can go before w or to the other slot. This case contributes 2^q 3^l.
Add the four cases and sum over selected ends. This derivation includes
the last-chain boundary correction.

## Joint approximate control, exact finite decoder, and limitation

For a prescribed normalized pair (u,v) in [0,1] times [0,2] and any
m>=2L, there is a literal cut satisfying
```
|r/2^m - u| <= 2^(-L)
|E/3^m - v| <= 2(2/3)^L + 2(2/3)^m.
```
Choose the first L bits of the binary expansion of u (at 1 use all
ones). Independently choose the last L bits to match the first L
digits of an expansion v=sum_(i>=1) f_i(2/3)^i. Such expansions
fill [0,2] because the two first-digit images [0,4/3] and [2/3,2]
cover the interval; the residual after L choices has diameter
2(2/3)^L. Choose the middle bits freely. The main part of E/3^m
uses exactly those reversed digits; the normalized remaining terms
are bounded by 2(2/3)^m. This proves joint control in the normalized
limit, not just two separate marginal limits.

**Exact ending:** at fixed m and exact r in [0,2^m−1] there is just
one cut of this specified kind, recovered by the binary expansion of
r. Its E is forced by the second sum. The apparent freedom of the
limiting rectangle disappears at exact precision. This is a concrete
obstruction and a design clue: leave r in an interval during coarse
choice and introduce a new physical buffer before demanding its final
bits, or change m or the cut family. Neither the rectangle nor the
existence of every gap 2r at sufficiently large m proves that the
corresponding absolute targets fill an interval.
