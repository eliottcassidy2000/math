# Reopening the next point of the old dyadic row

The [first ternary construction](../ternary_deposit/README.md) has source
$q=84\cdot2^s-5$. The [493 handoff](README.md) is at $s=1$. At $s=2$,
the original source is 331 and its $+4$ child has 997 paths. A second
exact fresh digit is possible here too, with a different continuation:
$$
997\ \rightsquigarrow\ 2991,2993,2995,\qquad(E,O,g)=(2897,2912,15).
$$
The arrow means literal re-encoding with a new host, not retaining the
493 reader or assuming an induced extension of the original 997 graph.
All three outputs here are strong. The new continuation is not
three times an integer continuation: neither 2897 nor 2912 is
divisible by three.

Use outgoing masks, vertex 0 the low bit, last vertex the root and
penultimate vertex the marked slot.

| Count | Row masks |
|---:|---|
| 2991 | 16350,16380,16344,16336,8160,8141,16256,7936,7680,7168,6144,4096,8192,4016 |
| 2993 | 8174,4084,4088,8178,8161,8128,3968,3840,3584,3072,2048,4096,1990 |
| 2995 | 2014,2044,2040,2000,4064,4041,1920,1792,1536,1024,2048,975 |

The same unique-exit run argument now certifies a third infinite
replacement band. Replace the marked slot by any nonempty
tournament of count $u$ and two-optional-slot count $R$. Then
$$
(H,E,O,g)=(2897R-2803u+du,\ 2897R/2,\ 2897R/2+15u,\ 15u),
\qquad d=0,2,4.
$$
The exact Hamilton two-run coefficient is 2897 in each witness.
Root deletion has one-run coefficient 15. Rooted cover patterns
$(0,1),(1,0),(1,1)$ have unsigned counts $(2912,2897,2897)$
and signed counts $(0,15,15)$; the signed filling formula gives
$E-O=-15u$. **node check_height_two.js** exhausts these finite
skeletons and verifies several literal substitutions by independent
subset DP; it does not require replaying the discovery search.

For transitive modules these are three consecutive odds at
$2897\cdot2^t-2803$, $t\ge1$. They meet the requested old branch
$756\cdot2^s-33+\{0,2,4\}$ exactly at $s=2,t=1$ when the
three left endpoints are equated. Indeed
$2897\cdot2^t-756\cdot2^s=2770$, whose right side has valuation
one; the second term is always divisible by four, so $t=1$, then
$s=2$. Together with the first-height match, this shows that changing
the reader can reopen *two* successive points of the old dyadic
row. A single common reader family that aligns for all heights is
still unproved.
