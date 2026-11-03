# Complete response ranks at every attachment depth

**PROVED, scoped all-size; same-agent proof review with independent finite
controls.** HB5 below extends [HB1--HB4](HALVING_BOUNDARY.md).

## HB5. A binomial classification

For N>=1 use T_x=C3[TT_(2N) with adjacent pair bits x,1,1], retaining the
two singleton anchor names. Let W_(N,r) be the rational span of the constant
baseline and all Hamilton-count functions of x obtained by at most r named
births with fixed cuts. One may restrict to contexts whose every intermediate
graph is strong. In either menu,

    W_(N,r) = span_Q {chi_D : D subset [N], |D|<=2r},
    R_(N,r) = dim W_(N,r) = sum_(j=0)^min(2r,N) binom(N,j).

Here chi_D(x)=(-1)^(sum_(i in D) x_i). HB1 proves containment. For equality
it suffices to produce one context with a nonzero D coefficient for each
allowed D. The context bank is closed under translating pair bits: permute
the corresponding old names in every cut. Fourier averaging over these
translations then isolates chi_D. This is a rational span argument, not a
single physical context implementing a signed average.

List D in ordinal order. Group consecutive selected pairs two at a time,
leaving one single pair if |D| is odd. Each two-pair group occupies the
interval between its two selected pairs; intervening pairs are unselected.
Add one reader beaten only by the first named vertex of each selected pair
in that interval, beating every other interval vertex. Put the reader in
the ordinal sum relative to outside vertices. For a one-pair group use the
same first-vertex cut. Keep all readers in the large C3 module. This uses
ceil(|D|/2)<=r births and all intermediate graphs are strong. Unselected
pair reversals are invisible, since the cuts are uniform on those pairs.

For a single-pair interval the (H,F2) difference is (-2,-4): its three-vertex
responses are (1,8),(3,12). For a two-pair interval with ell intervening
vertices, apply HB1 cancellation to terms in which a selected pair block
does not touch the reader. Every remaining term has both selected bits
equal to 1, hence positive sign. For H the two possible five-vertex blocks
are first-pair/reader/second-pair and its reverse. The forward order allows
no intervening vertex before or after the block, so contributes [ell=0].
In the reverse order each intervening vertex can go before or after the
block, each portion in its unique increasing order. Therefore

    mixed H = 2^ell + [ell=0] > 0.

For F2, the two ordered possibly-empty path bins, the remaining terms again
all have positive sign. Placing a Hamilton path in one bin and leaving the
other empty gives mixed F2 >= 2*mixed H > 0. No fitted F2 formula is needed.
Uncovered ordinal gaps have H=1 and positive F2.

Across the intervals H and F2 multiply, so their D coefficients are
products of the local differences. Both have the same nonzero sign:
positive for even |D|, negative for odd |D|. The wrapper has H=H(S)+F2(S),
so its coefficient cannot cancel. This supplies each chi_D. QED.

`rank_depth.py` exhausts every one-birth cut for N=1,...,5, giving ranks
2,4,7,11,16, and checks the spaced-reader formula for ell=0,...,6. HB4
independently checks depth-two full-parity readers. These computations
control the proof rather than provide an extrapolation premise.

## HB6. Recurrences, completion and the even-subset grading

Each new depth layer has exactly

    R_(N,r)-R_(N,r-1)=binom(N+1,2r), r>=1

directions, by Pascal's identity on degrees 2r-1 and2r. An explicit index
map leaves an even-cardinality Walsh subset D unchanged and appends a new
symbol * to an odd-cardinality D. The result is an even subset of N+1
symbols, with size at most 2r. Thus these spaces have the grading of the
even exterior-vector basis on N+1 symbols. This preserves indices and
grading, not multiplication: Walsh characters multiply by symmetric
difference, while overlapping exterior monomials have zero product.
No algebra isomorphism or tournament/even-graph class bijection is asserted.

At fixed r, R_(N,r) is a polynomial in N of degree 2r, hence its exact
minimal constant-coefficient rational recurrence is (E-1)^(2r+1), where
E increments N. With the formal convention R_(0,r)=1, its ordinary GF is

    sum_(N>=0) R_(N,r) z^N = sum_(j=0)^(2r) z^j/(1-z)^(j+1).

At fixed N, the depth GF is

    sum_(r>=0) R_(N,r) y^r
      = (sum_s binom(N+1,2s)y^s)/(1-y).

All rational joint weights become observable exactly when r>=ceil(N/2).
With total mass retained as baseline, full weight recovery then needs
exactly 2^N-1 additional scalar contexts: dimension gives necessity and a
basis of the spanning context bank gives sufficiency. Below that depth,
full weight recovery is impossible. Positive laws already witness failure,
by adding a sufficiently small hidden direction to the uniform law.

This is linear-law completion, distinct from the often smaller context
count for identifying individual states in OBSERVER_COMPLETION. Integral
Smith layers and extra actual complements remain OPEN. The rank recurrence
counts response directions on this supplied strong family; it is not an
A000568 class-count recurrence or a full Hamilton spectrum theorem.
