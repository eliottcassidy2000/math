# Two nonbinary layers and an exact even-perfect-midpoint screen

2026-11-26. Extends the [commuting unit forks](COMMUTING-UNIT-FORKS.md). These are explicit constructions and scoped negative inverses, not a proof of the Chan–Pak conjecture.

## 1. Nest the seam without disturbing the endpoint channel

Fix one of the three outer left cores, with profile $(H,E,R)$ and isolated edit $\Delta E=2$. For the outer right block, take the seam of a rooted core $C$ with terminal count one and a converse binary block at address $c$. Denote the inner core's path count by $h$, root starts by $e_+$, ending-first cover by $r$, starting-second cover by $v$, and all ordered two-slot covers (allowing empty slots) by $Q$. By the seam and cover formulas, the outer count for odd $c$ is exactly

    N(c) = A F(c−1)+B c−D,
    A=H h,   B=H r+R Q,   D=H e_+ + R v+E.

The outer right root still has one terminal path. Thus changing $E$ shifts $N$ by two independently of the same inner spoke inverse. Given $\ell=2^i$, $i\ge1$, the same unit equation with $(A,B)$ works: $B\ell\equiv\sigma\pmod A$, and $\beta=(B\ell-\sigma)/A+Z+\epsilon$. Orient the outer edit against the sign $\sigma$. Every vertex of the inner core can reach its root by the unique ending path, so the bridge and highest right bit make the nested graph strongly connected whenever the constructed address has positive prefix. The disjointness of the two physical edges is unchanged.

For example, use outer $(H,R,E)=(25,28,8)$ and the seven-vertex inner masks $(126,60,56,48,96,64,14)$ rooted at 6. Its profile $(h,e_+,r,v,Q)$ is $(23,0,50,29,592)$, and $(A,B,D)=(575,17826,820)$. With $c=253$ and $i=1$, the exact four counts are $(4618983,4618981,4618985,4618983)$. The replay checks 110 separate literal nested orientations against the profile identity, as well as this profile-level fork and exact-inverse positive controls.

## 2. A finite bank with no cutoff on an ordered prefix

Prepending an ordered chain of $k\ge1$ vertices to the inner core keeps $h$ and $v$ unchanged, makes the root-start count zero, and multiplies both $r$ and $Q$ by $2^k$: the prefix vertices can be assigned independently to the two slots. So the effective slope becomes $2^kB$ and the offset becomes $Rv+E$. This covers arbitrarily long such prefixes, not only the enumerated base orders.

For target midpoint $M$, apply the [square-root window](README.md) to $M+D$ and the affine count $A F(c-1)+Bc$. In the prefix case, first group by $u=B2^{k+i}$, since $w=u+AZ$ and $u\le w<3u/2$ when $k\ge1$. This gives a short interval of powers for $u$. Reject using the congruence before examining the finite splits $k+i$; each split then only reads the small square-root window of suffixes.

I enumerated every uniquely terminal core through order 12 by the [Fibonacci grammar](UNIQUE-TERMINAL-GRAMMAR.md). Combined with the three outer cores this gave 53,103 distinct unpadded parameter triples and 53,097 padded ones. **None crosses the even-perfect midpoint at exponents 19, 31, 61, 89, 107 or 127**, on either right-edit copy of the fork. At exponent 61 there were 22 padded leading congruence lanes, but not a single compatible split of the power between prefix and spoke. This is a restriction on the stated construction only.

## 3. The basic three cores: an all-exponent scale restriction

For the nonnested cores, let $U=2^{p-1}$ and $M=U(2U-1)$. The square-root center differs from $U\sqrt H$ by less than $\sqrt H/2$. Together with the radius $12H$ and $R\ell\le w<(R+H)\ell$, this forces $i=p-4$ for every $p\ge12$: for the three profiles the lower candidate $\ell\le U/16$ has $(R+H)/16<\sqrt H$, and the upper candidate $\ell\ge U/4$ has $R/4>\sqrt H$, with strict gaps exceeding the radius divided by $U$ already at $p=12$.

The $(25,28)$ and $(29,31)$ congruence classes then both require **even** $p$. The $(23,24)$ profile can survive at odd $p$ only for $p\equiv5\pmod{22}$. At all 228 such candidate exponents through 5,000 the exact inverse found no crossings after 21,860 suffix reads on the two left orientations. Direct checks cover the six smaller odd exponents below 12. This search is only for the signed binary spoke in these three forks. The modular progressions give one fast way to reject targets, not a global obstruction to other families.

Replay: `node nested_fork_inverse.js 12 --perfect`.
