# Two large contacts cancel to an exact +2 in one reader

Use the same seed throughout: outgoing rows (2,4,25,19,35,15),
root 3, and the usual initial synchronization zero. In the
[reader coordinates](../variable_reader_coverage_20261116/PROOF.md)
this gives (I,Z,g)=(-11,82,9). The n-bit address b, followed by
the usual read zero, gives the literal count

    H=I+Z(2D-rev_n(b))+g(2b-J(b)),  D=2^n,
    J(b)=sum_{j>i} [2^(j-i-1)+1_{j=i+1}] b_i b_j.

Take u≡5 mod 6, u≥5, put V=2^u, i=u+3, W=2^i=8V, and let the
address consist of u high bits U, a three-bit patch, and i low
bits L. Write X=rev_i(L). Change the patch from **010 to 111**.
The two edited arcs both meet the reader root; no seed or exterior
arc changes. The changes in b and its reversal are 5W and 5V.
The new pair interactions add 6 internally, 5(X+U) from the
exponential cross kernel, and (X mod 2)+(U mod 2) at the two
adjacent contacts. The exact count difference is therefore

    5(18W-82V) - 9[5(X+U)+6+(X mod 2)+(U mod 2)].

The phase assumption is precisely 2^u≡5 mod9. It makes
K=(62V−13)/9 an odd integer. Choose **any** integer U in [0,V),
set X=K−U and reverse X to obtain L. Indeed V<K<W, and the
two adjacency bits now sum to one. The displayed difference
is exactly 2. The final read makes both outputs strong.
At u=5 there are 32 choices and the first supplied order is 24:
U=0 gives the two addresses 731 and 2011, with actual
Hamilton counts 6,152,939 and 6,152,941.

## The choices here give distinct absolute anchors

At patch 010 the only mixed interaction terms are 2(X+U) and
8XU. In H_0, the smooth U-varying part after substituting
X=K−U is

    (1808V−72K)U+72U².

The remainder is 9(2L−J(L))−82 rev_u(U)−9J(U), which ranges
over an interval of length at most 235V since 0≤J(w)≤w.
The smooth polynomial has derivative at least 1232V (K<8V),
so consecutive U give strictly increasing counts, in fact
with spacing Θ(V). Thus the construction gives V unambiguous
physical +2 anchors at scale Θ(V²), not merely V programs
whose outputs might all collide. The lower-order spacing
and large gaps between scales prevent claiming exact placement.

## What the switch does not preserve

A second synchronized reading by another zero has
H(w0)−H(w)=2ZD+2g b(w). The two repaired counts at the first
reading therefore differ by 2+90W after the next read, not by 2.
The fixed exterior is protected against the local two-bit edit,
but this is not an exact return with arbitrary later continuation.

There is a general reason for this failure within a fixed-width
binary address. If synchronized sources have the same odd g,
width n, and final Z, then Z=2^n Z_0+g b. Equal Z for two such
words forces 2^n|Δb because gcd(g,2)=1. Both addresses are in
[0,2^n), so their addresses and initial Z_0 are equal. A
nontrivial edit with a fixed width and fixed g cannot be
continued simply by pretending the final slope is unchanged.
A change of reader or width, or another interface, is necessary.

Replay: `python check_seam.py` tests the count identity across
many phases and checks both order-24 counts directly from the
literal adjacency by splitting Hamilton paths at the root.
