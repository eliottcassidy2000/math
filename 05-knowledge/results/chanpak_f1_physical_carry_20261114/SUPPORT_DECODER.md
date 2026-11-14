# The pointed support decoder and a rigid triangle-only carry family

A finite-field point count at q=1 forgets precisely the positive phase
ranks. In the repo's [pointed-monoid model](../session59_f1_geometry_and_information_20261025/SYNTHESIS.md),
a one-unit phase instead leaves one point **per actual support**. This
suggests making the arithmetic digits be concrete contact choices,
with the directed sense of each contact retained. Here are two exact
instances: a complete inverse for any fixed profile, and an all-order
rigidity result in which signed positional digits recover the contact
supports that matter. Neither uses a literal field with one element.

## Exact physical decoder at fixed H, F2, and exposed positive response

For a puncture tournament B and an added apex v, let U be the set of old
vertices beaten by v. Write h=H(B), f=F2(B), where F2 is the sum of
H(B[S])H(B−S) over all subsets S, and write P(B,U) for the sum over
u in U of all positively weighted subset paths starting at u.
The weight of a path on S is H(B−S). The dual Q(B,U) is the same
response in the converse with contacts replaced by their complement.
The ordinary augmented count K is H(B+v). The six-coordinate
switch interfaces additionally preserve signed responses S,T; the
inverse below can filter by those as well.

Every puncture has an ordered strong-component presentation
\[
 TT_{l_0}\to B_1\to TT_{l_1}\to\cdots\to B_m\to TT_{l_m},
\]
where B_i are precisely its nonsingleton strong components. For any
fixed h their ordered list has a finite exact library: ∏H(B_i)=h,
m≤floor(log_3 h), and each strong component with count h_i has order
at most 3(h_i−1)/2. The last bound needs no minimal-order census:
every vertex of a nonsingleton strong tournament belongs to a directed
triangle, and the [odd-cycle expansion](https://arxiv.org/abs/2307.05569)
gives h_i≥1+2 times its number of directed triangles. One may use
[Busch's exponential bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
to build a much smaller library. Count and store all the contact masks
on all listed cores. The singletons are not library items; they are
the arbitrarily variable runs.

For the selected list, the sum of the run lengths has to be
\[
 D=l_0+\cdots+l_m=\log_2\bigl(f/\prod_i F2(B_i)\bigr)
\]
and must be a nonnegative integer. Enumerate its weak compositions.
For each and for each choice of **core** contact masks, there is at most
one possible choice of **all singleton** contact bits. It is found as
follows:

* For an ordinal join A→C, the positive start sum is
  \(P(A\to C)=F2(C)P(A)+H(A)P(C)\). Paths starting in A
  split C over their suffix and complementary path; paths starting in C
  leave A unselected. This factorization includes empty sides.
* When peeling a terminal singleton run of length l, its local P is the
  ordinary binary address of its l contact bits. With current count h',
  \(P_{\rm in}=2^lP_{\rm prefix}+h' z\).
  Recover \(z\equiv P_{\rm in}(h')^{-1}\pmod {2^l}\) with
  0≤z<2^l. Since tournament counts are odd, this is unique. Peel by
  exact nonnegative division.
* When peeling a nonsingleton, its response P(C) is already tabulated.
  Subtract \(H(A)P(C)\), then divide by \(F2(C)\), rejecting if this is
  negative or nonintegral. The process ends at the empty prefix with
  P=0.

The output is the actual directed puncture and its apex contact mask.
For a prescribed switch block digit, recompute K and check all unused
coordinates exactly. If all possibilities reject, the **entire** fixed
fiber fails, not just graphs below an arbitrary order bound. This is an
all-order source compiler for a fixed interface, not a proof that the
right interface exists at a specified Chan–Pak target.

For fixed h, the total number of candidate masks in this compiler is
\[
 O_h\big((1+\log f)^{\lfloor\log_3 h\rfloor}\big).
\]
The bound applies before filtering Q or signed endpoints and counts
presentations, so it also bounds distinct K digits. It follows from
the number of weak compositions of D into m+1 parts and the finite
core-mask library. If the center lies between nonempty positive
transitive wings in the two-wing chart, f≤N follows from the fxy term.
Thus arbitrary neutral padding at constant h does not create an
exponentially growing supply of independent physical choices at a
single exposed start response.

## Triangle-only punctures have no nontrivial independent carry digit

Suppose every nonsingleton component of the puncture is a directed
triangle. This allows arbitrarily many triangles and arbitrary
singleton padding. Then **(h,f,P,Q) determines K** without any of the
signed endpoint readouts. In particular no fixed-profile K switch is
possible among triangle-only punctures.

Encode a singleton by its contact bit, and encode a triangle by the
two binary digits of its number k of contacts (0≤k≤3). Let A be the
concatenated word of declared length L. For a triangle, (h,f)=(3,12)
and each individual unsigned rooted start response is 3, so the
normalized forward digit is k and its radix is f/h=4. A singleton's
radix is 2. Therefore f/h=2^L, and A is exactly the L-bit expansion of
P/h by the ordinal factorization.

Now form \(\operatorname{rc}(A)\), the word in opposite order with
complemented bits, and let \(Z\) be its integer value. Physically the
converse with complementary contacts uses triangle digit 3−k, not
the reverse-complement of the two original bits. These agree for
uniform contact k=0,3, but differ by +1 at k=1 and by −1 at k=2.
Thus if j is the starting position of a two-bit block measured from
the **left** starting at zero,
\[
               Q/h-Z=\sum_{j\text{ mixed triangle}}\epsilon_j\,2^j,
       \qquad \epsilon_j=+1\ (k=1),\quad -1\ (k=2).
\]
The exponent is j because the reverse ordering puts the low digit
of this reflected slot at binary place j from the right. The
nonzero signed positions are never adjacent: two different triangle
windows are disjoint. Their signed expansion is unique. In a nonzero
integer, its 2-adic valuation locates the first nonzero digit; after
dividing by that power, its residue mod 4 gives the sign (\(+1\) or \(-1\)).
Remove that signed digit and iterate. For zero only the all-zero
expansion is possible by the lowest nonzero digit test. This recovers
every **mixed** triangle position, including its contact number.

Uniform triangles can be erased into two equally contacted
singletons after normalizing the whole augmented count by h. To see
this, take the [common-apex gluing formula](../chanpak_movable_reader_returns_20261114/INDEPENDENT_RETURNS.md).
After division by the product h, it depends only on each block's
normalized \(f/h,P/h,Q/h,S/h,T/h,K/h\).
A directed triangle with all contacts zero or all one agrees in
these six data with the transitive pair carrying the same uniform
contact: both have normalized radix 4, K/h=1, unsigned digits
0 and 3 in the appropriate order, and signed digits 0 and −1 in
the appropriate order. Replace the uniform triangle by that pair.
For a mixed triangle the position and k have been recovered;
any two masks of that size are cyclic rotations and therefore
give the same physical block. All remaining places are singletons.
The normalized ordered graph thereby becomes determined by the
decoded marked word. Multiplying by the original h recovers K.

The support-based lesson is constructive, but negative for this
family: equal phase mass is not an arbitrary addition operation.
Here the directed conjugate response remembers exactly which
bit windows are mixed, so the apparently free ternary modules
supply **no** physical carry choice. The two new coprime nonlinear
switch centers lie outside this rigid family.
