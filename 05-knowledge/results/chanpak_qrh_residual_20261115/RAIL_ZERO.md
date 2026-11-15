# A second low-bit rail: target height ending in zeros

The [adaptive complementary selector](ADAPTIVE.md) keys on the residue
−1 modulo g. There is a second exact cancellation, for totals S=x+y
congruent to zero modulo an even power g of two, on the odd-address grid.

Write x=gu+b with b odd, hence 1<=b<g, and write y=gv+c where
c=g−b is also odd and u+v=S/g−1. Using the definitions of A in the
[low transfer](COMPLEMENT.md), we have

    c=bar(b−1), A_c=(g−1)−A_(b−1),
    A_b−A_(b−1)=−g/2, A_b+A_c=g/2−1.

The second equation changes only the least significant input bit, which
is the most significant weight in A. The sum of slope terms therefore
becomes

    A_b u+A_c v=(g/2−1)(u+v)/2+(A_b−A_c)(u−v)/2.

The first term is constant as x moves by two and y moves by minus two,
even across high-part block boundaries. As before the remaining
checksum has size O(d+g+S/g); the quadratic changes O(d+1) on the
step-two grid. Use the first odd x>=S/2, changing the initial d by at
most 1, and bisect signs between the above and below endpoints.

With s0=ceil sqrt(N/108), choose S=s0+delta and delta=(−s0) mod g.
The end of the sign search remains d=O(sqrt(L(delta+1))) where
L=sqrt(N/432), giving the same error formula

    O(sqrt(L(delta+1))+g+L/g).

Minimizing over both rails is a literal target-dependent selector; it
does not interpolate between them. For example s0=2^(2r), g=2^r has
delta=0, so the conclusion is O(N^(1/4)) on the entire interval
108(s0−1)^2 < N <= 108s0². The earlier −1 rail covers the corresponding
all-ones blocks. These pockets have the same three exact count choices.

Replay: `python zero_rail.py` checks all odd b through g=1024, 30,000
consecutive targets, and 256 aligned zero-ending target heights.
