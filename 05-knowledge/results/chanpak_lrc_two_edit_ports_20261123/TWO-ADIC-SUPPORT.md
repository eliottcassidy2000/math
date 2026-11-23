# A two-adic support cap for a literal prime-to-sevenfold square

**PROVED**, using the [cross-component singleton budget](../chanpak_kaprekar_decimal_20261123/FIRST-JUMP-SQUARES.md#cross-component-gain-and-singleton-budget). The cap concerns a given two-edit witness or the corresponding decision problem, not existence for every prime.

Let \(p\ge5\) be prime, \(p\ne7\), and let the opposite endpoints of a literal two-flip square have heights \(p,7p\). Each has exactly one nonsingleton SCC by [atomicity](../chanpak_two_holes_atoms_20261123/README.md); call their supports \(A,B\). Set \(b(t)=\lfloor1+3\log_5 t\rfloor\). The reduced host \(A\cup B\) has order at most
\[
 \min\{b(p)+b(7p)-3,\ b(7p)+2\}.
\]
The two marks survive this restriction whenever a one-edit witness is impossible.

If neither mark is internal to \(A\), then \(A\subseteq B\) since the strong induced copy survives. If both are internal, then \(B\subseteq A\) by ordinal separation. Otherwise do the internal reversal first. If it preserves the strength of \(A\), again \(A\subseteq B\). In the remaining case let \(T\) be that split intermediate, with height \(q\), and let \(S\) be the union of its nonsingleton SCCs. The case \(q=1\) is impossible: a single flip from transitive has count 1 or \(1+2^j\), never zero mod seven. Each of the two edits from \(T\) back to \(p\) or forward to \(7p\) crosses SCCs of \(T\), merging a full ordinal interval; the new unique nonsingleton SCC contains \(S\). Both gains are positive by the rooted gain formula (an adjacent singleton swap would leave the height and all SCC orders unchanged).

The two gains differ by \(6p\), so exactly one has 2-adic valuation one. Neither endpoint is divisible by three, so the ratio of that arm is not three. The valuation-one arm absorbs at most two singleton vertices. In fact, if \(T\) has two or more nonsingleton SCCs it absorbs at most **one**: absorbing two without two singleton endpoints would require a singleton endpoint plus a singleton middle, exhausting the budget; there would be no nontrivial middle and only one available nontrivial endpoint. Two singleton endpoints cannot enclose a nontrivial component at valuation one.

If \(q\equiv p\pmod4\), the short arm is the one to \(7p\), and \(|A\cup B|\le |A|+2\le b(p)+2\). If \(q\equiv7p\pmod4\), the short arm is the one to \(p\), and the analogous cap is \(b(7p)+2\). Both +2s are +1 on the multicore branch. The older cap follows from \(|A\cap B|\ge3\): in the split case \(S\) contains a directed triangle, and the nested cases are immediate. Inducing on \(A\cup B\) preserves both endpoint heights since discarded vertices are ordinal singleton pads in each endpoint.

**Search form.** It suffices to start from a strong core at one of the two endpoint heights, add 0, 1, or 2 singleton pads at its ends, and test pairs of arcs of that host. On the nonnested branch, use the endpoint opposite the valuation-one arm as the starting core; on the nested branch use the larger support and no pads. This is an all-order, per-prime finite decision procedure for distance at most two, not a uniform construction.

**Sharp support example.** On 12 vertices use row bitmasks \((4094,4092,4088,4080,4064,3392,3200,3360,3648,3296,2048,0)\) at the split corner of height 13. The flips \((7,11)\), \((0,7)\) individually give 43 and 301; doing both gives 1099. The supports at 43 and 301 are respectively vertices 5 through 11 and 0 through 9, so the union exceeds the larger support by two. The short gain is \(43-13=30\). This is not a lower bound on the smallest host for 43; the [order-eight certificate](two-port-witnesses.js) already works.
