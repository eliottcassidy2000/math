# Three module amplifiers for the 67/69 pair

## Exact scope

For every nonempty tournament B let h=H(B) and
Q=sum_{A⊆V(B)} H(B[A]) H(B−A), with H(empty)=1.
There are literal tournaments with the six counts

    sQ+67h, sQ+69h,  s=72,73,74.

At a transitive module on k≥1 vertices, h=1 and Q=2^k, so each slope
gives an exact +2 pair at every such height. At the singleton source
instead of the wrapped module, the same six readers realize 67 and 69.
In particular, 67 has an explicit eight-vertex witness. The literal
certificate is port_switch.py; each code is an adjacency matrix, not a
proposed arithmetic coordinate.

These are exact positive families, not a range theorem. A target N
still requires the appropriate physical (h,Q), not an arbitrary pair
of integers satisfying an affine equality. At slope 72 the 69 branch
is exactly 3(23h+24Q), whereas the 67 branch does not have that forced
factor. With ordinal modules h and Q multiply componentwise, so a
common divisor in (h,Q) remains a permanent divisor of either readout.

## Source interface

Use the rooted cyclic substitution T=C3[r,w,B], with r→w→B→r.
For the base-2 port-chain transfer in THM-5070, the state is
(H,E,O): E (respectively O) sums H(T−P) over directed rooted paths P
of even (respectively odd) number of vertices. The full state of T is

    (h+Q, Q/2, h+Q/2).

To see why the even and odd halves of Q agree, use the odd-cycle packing
formula H(U)=sum_C 2^|C|, where C ranges over vertex-disjoint directed
odd cycles of U; all remaining vertices are singleton blocks. In
the signed convolution sum_A (-1)^|A| H(A)H(B−A), fix the resulting
packing on B. Every one of its odd-sized blocks independently
chooses side; summing its two signs gives zero. Since B is nonempty
every term cancels. The same unsigned expansion gives Q≥2^|B|.

A source path starting at r either stops immediately (odd contribution
h), or next uses w and then follows a Hamilton path on A⊆B; its parity
is |A| and its deleted-tournament weight gives exactly one convolution
summand. For H(T) itself, packings avoiding r live completely inside
B; if the packing contains a cycle through r, that cycle must use w
plus a Hamilton path on A⊆B of odd positive size. Each such cycle gets
weight 2 and can be accompanied by any packing in B−A. Thus H(T)=h+Q.

## Actual reader cores

Use the inherited HYP-9640 / THM-5070 port-chain operation: take ordinal
sum with the next block but reverse the edge from the current port
to the indicated entry of that next block. Continue from its indicated
exit port. In each case first attach an order-five core and then attach
the same order-two block with adjacency 0010 and entry/exit (1,0).
Vertices are zero-based; adjacency strings are row-major. The effective
first rows of the two-block products are:

| Count term | Core code | Entry, exit | Effective first row |
|---|---|---|---|
| 72Q+67h | 0000110010110001010001110 | 0,1 | (31,46,36) |
| 72Q+69h | 0000110010110001010001110 | 4,4 | (33,42,36) |
| 73Q+67h | 0001110001110000110000110 | 0,1 | (33,46,34) |
| 73Q+69h | 0010110000010001110001110 | 0,4 | (29,48,40) |
| 74Q+67h | 0000110010110001010001110 | 1,4 | (33,48,34) |
| 74Q+69h | 0001010000110000110011110 | 4,4 | (37,42,32) |

The rows were calculated using the original exact
kps_s18_ternary_chains_digraph_weights_20261109/chain_blocks.py.
For instance the two slope-72 readers have differential row
(2,−4,0), so the matching reservoir direction is not a coincidence
of integer factorizations: it is annihilated by the change of reader.
On the literal source it reads 2(H−2E)=2h. The finite target inverse
loop can restrict |B|≤floor(log_2(N/s)), but this still leaves the
attainability of residual module coordinates as the central issue.

## Exact obstruction in the rooted O/Z/J language

The pair distinguishes a limitation and a repair. The eleven-block
O/Z/J bank contains 69, for example OO·JOZZZ, with terminal state
(H,E,L,h_root)=(69,1,73,1). It misses 67 even if ALL run caps and
block restrictions are removed. This claim follows from exact finite
closure, not from merely inspecting a list of short words. At bound
N=69, keep every state with H≤N and L−E≤N−H and record (without
expanding) each terminal violating the latter inequality. On these
states E≥0, h_root≤H and E≤H, so expanding cases have L≤N and the
closure is finite. O leaves H unchanged and cannot lower L−E, while
the first later Z or J increases H by that gap. Thus pruning cannot
lose any new reachable count. The replay closes 116 expandable states
and 176 terminal states and never records 67. This is a finite exact
separation of the grammars; it says nothing by itself about cofiniteness
of any full tournament spectrum.

## Literal mixed-base modules

Take B=C3[1,1,TT_m] with m>=1, itself the cyclic wrap of an
m-vertex transitive module. Then h=2^m+1 and
Q=2^(m+2)+2(3^m−1). For Q, assign every odd-cycle packing block
to either side of the convolution. A directed cycle through the
two singletons and an odd a-subset of the spine has weight
2^(m+2−a) after the side assignment; there are binom(m,a) such
cycles, and all of them intersect. The sum of the odd binomial
terms is 2(3^m−1). Consequently the six readers also realize the
mixed-base families

    (4s+c)·2^m + 2s·3^m + c−2s
    for s=72,73,74 and c=67,69.

For the slope-72 pair these are 355·2^m+144·3^m−77 and
357·2^m+144·3^m−75. The natural factors 2 and 3 here measure,
respectively, a free two-way side choice and the binomial sum
that allows an odd spine subset to enter the unique overlapping
cycle. This explains the powers structurally but not target
coverage. In particular, to get a prime from this lane one
necessary test is gcd(2^m+1,3^m−3)=1, since Q≡2(3^m−3) mod h.
