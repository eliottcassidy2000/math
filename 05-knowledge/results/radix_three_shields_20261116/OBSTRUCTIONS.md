# Where multiplication by 3 stops acting like a radix

## A literal 3^k amplifier, without a fresh low digit

Let the rooted tournament state be (I,Z,g) as in the
[variable-reader block law](../variable_reader_coverage_20261116/BLOCKS.md).
A uniform near-sink block B sends it to

    (H(B) I, P_B(2) Z, H(B) g),

and the physical read immediately after the append is
H(B) I + P_B(2) Z. A directed triangle C3 has H=3 and P(2)=12.
Starting with the synchronized four choices I=−75+2i, Z=334,
g=15, i=0,1,2,3 from the
[five-count bank](../variable_reader_coverage_20261116/FIVE.md),
append k near-sink triangles. For k≥1 the actual four counts are

    334·12^k + 3^k(−75+2i),         i=0,1,2,3.

So factors 3^k are exact at every height: at k=1 these are
3783, 3789, 3795, 3801. The new triangle is attached through
the retained root; this is not merely a formal ordinal multiplier.
But it scales the entire old window. It does not create an
independent +2 dial beside the amplified one.

In fact every uniform block of Hamilton count divisible by 3
makes **all three** resulting coordinates zero mod 3. The
reason is P_B(2)≡P_B(−1)=(-1)^|B| H(B) mod 3: in the
positive odd-cycle packing expression for the slot polynomial,
each exponent has parity |B|. The ordinary-sink block law
also has every coefficient divisible by 3 when H(B) is.
Subsequent ordinary or near-sink integral transfers at the
same root cannot leave this sublattice. No additional binary
letters at this root refill the unit digit. To make a radix
the next rule has to read a new coordinate (for instance a
different root), or be genuinely nonuniform.

A symmetric array of exact synchronized blocks is governed
by weights w_j=∏_(i≠j) h_i, where h_i are fixed local Hamilton
counts. It cannot give 1,3,9 even abstractly: the occurrence
of w_j=1 forces h_i=1 for every other i, making all their
weights equal to h_j. The physical same-profile obstruction
is worse: the transitive h=1 block has fixed contact mask when
its full boundary address is fixed. This narrow obstruction
does not apply to the [two live dials](TWO_DIALS.md), whose
upper tape uses a different response, not this product law.

## Even a cyclic triangle shield cannot multiply a four-count window twice

The [triangle wrapper](TRIANGLE_WINDOW.md)
W(U)=C3[U,1,1] acts on its ordered-slot polynomial by the exact rule

    (M f)(z) = z[f(z+1)+z f(z)−f(z−1)].

Take the difference f of the polynomials of any two choices in a
palette. Nonempty modules have f(0)=0. Exact gap multiplication
by 3 on the first application means f(2)=2f(1), since
(M f)(1)=f(2)+f(1). A second exact factor 3 requires f(3)=0:
directly (M² f)(1)=2f(3)+5f(2)−f(1).

If the two initial witnesses have the same order parity,
f(3)≡f(−2)=±f(2) mod 5. Thus f(3)=0 and f(2)=2f(1)
force f(1) divisible by 5. Any four choices with initial
counts A,A+2,A+4,A+6 contain two of the same order parity,
whose nonzero difference belongs to ±{2,4,6}. Hence **no
four-count window of consecutive odds at arbitrary orders
can scale by 3 twice via this same cyclic wrapper**.
For a third exact step, the next triangular relation gives
f(4)=2f(1); comparing 4 and −3 modulo 7 imposes a further
factor 7 on the original same-parity gap. The obstruction
comes from curve parity and a finite congruence, not just
eventual rate of growth.

The concrete palette with f(z)=2z has steps 2,6,30,246,... .
Even if those three first steps could be made independent,
the two lower width-four dials fill shifts 0 through 24,
whereas the third starts at 30. The two values 26 and 28 are
missing at the onset. One must reset the profile rather than
declare the second scale to be 18.
