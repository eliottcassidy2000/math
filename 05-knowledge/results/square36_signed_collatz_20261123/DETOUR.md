# Shortest square-sum detours across an even perfect number

**PROVED.** For every odd exponent \(p\ge3\), the two odd numbers \(E_p-1,E_p+1\), where \(E_p=2^{p-1}(2^p-1)\), have distance exactly three in the square-sum graph on positive odd integers. They lie opposite each other on an induced six-cycle. The formulas work for composite \(2^p-1\) as well.

Write \(r=2^{(p-1)/2}\), so \(E=r^2(2r^2-1)\). For each of the two signs put \(C_\pm=2r^2(r\pm1)^2\). The path
\[
                       E-1,\ C_\pm+1,\ C_\pm-1,\ E+1
\]
has outer edge sums \(E+C_\pm=r^2(2r\pm1)^2\) and middle sum \(2C_\pm=[2r(r\pm1)]^2\). Since \(r\ge2\), the two paths have distinct interior labels, all positive. A direct edge would require \(2E\) to be square, although its odd part \(2r^2-1\) is 3 modulo 4. A two-edge path would give squares differing by 2, impossible modulo 4. This proves minimality. The only remaining possible chords between opposite mod-four types have sum \(C_++C_-=4r^2(r^2+1)\), also nonsquare. Same-type sums are 2 modulo 4. These are ordinary undirected square-sum paths. Requiring each step to be a signed Collatz move would forbid even two consecutive steps. [Signed square steps](SIGNED-SQUARES.md)

## At exponent 61 the four intermediate counts have explicit tournament witnesses

**FINITE-EXACT.** Let \(a=E_{61}-1\), \(b=E_{61}+1\). The following two paths use different interior labels from the uniform formulas. All four interior labels are attained Hamilton-path counts, as is \(a\). Write the path in either row as \(a,u^2,s^2-u^2,b\).

| \(u\) | Middle square root \(s\) | First edge root | Last edge root | Orders realizing the two interior counts |
|---:|---:|---:|---:|---:|
|416971683598778710015|1000314801110273882448|416974871400173936320|909267273607204199216|156, 160|
|18198865913699551345185|31770071780973054887002|18198865986738607702920|26041097194277759566134|180, 173|

All six integers are distinct. Exactly the six displayed cycle edges are squares among their fifteen pairwise sums. Thus **membership of this specified induced six-cycle in the square-sum graph on attained counts is equivalent to ordinary tournament realizability of the single integer \(b\)**. Padding with initial sources would put any six realizers on a common order. Individual realization of the central pair would still leave their one-arc adjacency open.

A strong tournament on 93 vertices already realizes \(a\). For \(b=2432582681\cdot1092853292237112554142488617\), the smaller factor has a strong 34-vertex realization. A realization of the larger factor would suffice by ordinal sum. The existing prime-fiber results leave both that route and direct strong realization of \(b\) open. [Exact lower construction](../chanpak_e61_compact_volchkov_20261122/STRONG-93.md) · [Upper-prime fibers and explicit small-factor witness](../chanpak_fermat_solids_seams_20261122/UPPER-PRIME-PROFILES.md)

## The divisor search

Choosing a square interior height turns the other two edges into a factorization problem. For \(d\mid a\) with \(0<d<\sqrt a\), put
\[
                    u=(a/d-d)/2.
\]
Then \(a+u^2=(u+d)^2\). Since \(a\equiv3\pmod4\), the two complementary divisors have opposite odd residues modulo 4 and \(u\) is odd. Two more edges with roots \(s,t\) require precisely \(t^2-s^2=b-u^2\) and \(s>u\). Enumerate factorizations of \((b-u^2)/4\) as \(mn\) with \(m,n\) of equal parity, take \(t=|m+n|\), \(s=|m-n|\), and require the stated inequality. This covers all possible positive even roots, up to swapping or jointly negating the factors. Both intermediate heights then factor predictably as \(u^2\) and \((s-u)(s+u)\).

The two displayed records use \(d=3187801395226305\) and \(73039056357735\). Each of \(u,s-u,s+u\) factors into the following tournament heights. The factor 49 is used whole.

| Path | \(u\) | \(s-u\) | \(s+u\) |
|---|---|---|---|
|First|5 · 83 · 7243 · 2129389 · 65145583|49 · 47 · 48119 · 10663217 · 493657|109 · 197573 · 1391029 · 47311571|
|Second|3 · 5 · 37 · 5051 · 27127 · 77551 · 3085921|31² · 379 · 17107 · 411883 · 5288203|47 · 1597 · 84121 · 732541 · 10803413|

Directed ordinal sums multiply Hamilton counts because a spanning path must traverse the blocks in order. Each listed factor has its own recorded tournament address, with no appeal to a conjectural coverage interval. For 25 distinct factors the certificate uses the two-root construction. For the four remaining heights 83, 109, 1597 and 84121 it uses one four-vertex core with apex masks 10, 1, 11 and 4 and right addresses 5, 10, 83 and 7593 respectively. The replay builds the actual rows, uses the inherited block path counts, constructs each ordinal sum and checks the complete six-label graph. [Two-root constructor](../chanpak_prime_frontier_selector_20261118/join_loads.py) · [Direct apex path decomposition](../chanpak_decimal_kernel_20261122/e15_pair.py) · [Certificate](detour.py)

Run `python3 -B detour.py` from this directory.
