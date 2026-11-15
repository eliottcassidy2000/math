# Recover every physical arm from its full prefix response

The [all-prefix response theorem](EQUAL_RESPONSE.md) shows that two arms can differ by a constant on every prefix only if they have the same slope and mixed coefficient (L,t). The following inverse works backward from that desired response and removes the apparent need to enumerate large cores when L grows but t does not. It is a complement to [the leading-coefficient inverse](ABSOLUTE_CARRY.md), which instead fixes an absolute target and a gap.

Use d=(Q_U−h_U)Q_W, t=h_U h_W and a=h_W(Q_U−h_U). For a suffix v of width k, write r for its reversed value. The response relation is
```
L=d·2^k+t(2^k−1−r),   0≤r<2^k.
```
It implies d·2^k≤L, so k has O(log L) choices. For each k there is **at most one** candidate d: it must be in the half-open interval
```
L/2^k−t < d ≤ L/2^k
```
and simultaneously d≡2^(−k)L mod t. Here t is odd and the inverse is modulo t (when t=1 all integers meet the congruence). Pick the sole residue representative in that interval, reject unless positive even and consistent with the exact response relation, and recover r and the literal suffix by reversal, keeping width. This includes k=0.

Next factor t=h_U h_W. For each h, [fixed-count tightness](../chanpak_crossing_scales_20261115/FIXED_COUNT_TIGHTNESS.md) gives the full Q bank as a finite union q_0·2^s with q_0 from nonsingleton strong components of total order at most floor(3 log_3 h). Store a witness for each reduced q_0. For h=1 use the empty reduced tournament and require s≥1 after padding to obtain a nonempty core.

For each pair of reduced q_0 values, the parity of Q_U−h_U is odd. Thus
```
s_W = v₂(d)−v₂(q_W0).
```
Reject if it is negative or gives an empty W. Form Q_W=q_W0·2^s_W, require Q_W divides d, then set the only possible Q_U=d/Q_W+h_U. Test whether Q_U/q_U0 is an allowed power of two. These checks produce literal U,W by adding the decoded ordinal singletons. The constant c of each output arm is then evaluated by subsets or the exact binary recurrence. In particular, for fixed t this is a finite-fiber precomputation followed by O_t(log L) candidate profile tests. The precomputation need not be cheap as t grows. The enumerator gives the complete constant list at a specified response, including any exact +2, not a theorem that such a list contains +2.

**Example.** The known exact fibers H=1,3,9 use reduced bases 1,12,{84,144}, respectively (for 9 use a literal H=9 strong core or two triangles). The inverse gives exactly two constants 207 and 243 at (L,t)=(252,9). These are the physical +36 arms; no hidden high-order padding yields another constant at that response. For changing t the factored profile banks also change, so this does not rule out an all-prefix +2 at a different response.

The replay exhausts low d,t,k congruence controls and compares the decoded constants for several small fixed-t responses to a separate direct enumeration over widths and padded profiles.
