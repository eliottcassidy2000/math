# Pair programming without losing the absolute bracket

Separate modular programming of two addresses changes their sum. In the
balanced two-tape construction, that moves the baseline at first order
and loses the tangent gain. The redesign programs a pair using
**complementary** suffixes. It still places both counts near the prescribed
integer, with error O_m(N^(2/5)) for each fixed odd modulus m.

## A controllable complementary suffix

Write p for the binary surcharge and q0(x)=p(x)/2. All arithmetic in
this subsection is modulo the fixed odd m (the m=1 case is trivial).
Start with any two nonnegative prefix addresses x,y with x+y≡−1 mod m.
When appending a bit b to x, append 1−b to y. After ell steps the
integer sum is exactly 2^ell(x+y+1)−1, no matter which choices were
made. In particular the sum remains −1 mod m.

Normalize x to −1 mod m and odd using ceil(log2(2m)) complementary
digits: the x suffix can have any value up to 2^ell−1, so choose
the desired residual address modulo 2m. Now work with x odd and
y even, and keep z=q0(x)+q0(y)+(x+1)/2 mod m. The two words on x
A=1 and B=01 (with their complements on y) act by
```
A(x,z)=(2x+1,z),    B(x,z)=(4x+1,z+x+1).
```
For A the checksum grows by y/2, and for B it grows by
(x+1)/2+y+1; inserting y=−1−x and the change of (x+1)/2
gives the displayed actions. The two words end in 1.

Let D=B A^(−2), with rightmost operation first.
Then D shifts x by −2 and adds (x+1)/4 to z.
Both A and D are invertible by positive repetitions: A has order
dividing ord_m(2), and D^m is identity (sum the affine vertical gains).
Both E=A D A^(−1) and D² shift x by −4; their respective vertical gains
are (x+1)/8 and x/2. Hence D^(−2)E fixes x and translates z by
(1−3x)/8, equal to the unit 1/2 at x=−1. Repetitions prescribe any z.
All cycles can be padded to the same physical length by A at x=−1,
where A fixes both modular coordinates. Use a final fixed-length
word to choose x mod m (precompensate for its deterministic surcharge).
This constructs complementary suffixes of one fixed length ell=ell(m)
giving any pair of prescribed residues for x and q0(x)+q0(y), from
any pair of prefixes of sum −1 mod m.

Equivalently, put D0=2^ell. If a desired exact sum K is
`D0(z0+1)−1` with z0≡−1 mod m, each block
```
U=D0 j+b, V=K−U,      0<=b<D0, 0<=j<=z0,
```
contains a programmed pair. The starting prefixes are j,z0−j and
the trailing addresses are b,D0−1−b.

## Physical bracket at a prescribed odd target

Use the constants B,Q,T,S,h from [two-tape steering](TWO_TAPE_STEERING.md).
For sufficiently large odd N and a fixed odd m, take a power of two q
of order N^(1/10), and K of the above form with z0≡−1 mod m. Choose K
just above 2 ceil(sqrt((N−B)/Q))/q + 4D0 + 10: the rounding adds at
most mD0. Let the two physical addresses be
```
x=qU+v,  y=qV+v,  U+V=K.
```
Here v, common to both addresses, is a bounded integer chosen next.
Put d=gcd(45,m). Choose 0<=v<d so
`B+(T+S)v−Sq≡N mod d`. This is possible because gcd(T+S,45)=1.
Then v<q for large N. Prescribe U=0 mod m using the suffix compiler,
so V=−1 mod m. The non-checksum part of F(x,y) is now fixed mod m;
its reduction mod d is B+(T+S)v−Sq because 45 divides Q.
By the binary split with K odd,
```
p(x)+p(y)=2p(v)+[q−1−rev(v)]K+delta + p(U)+p(V),
```
where delta is one minus the top bit of the r-bit low word, q=2^r.
The first three terms are fixed, and the variable term gives 90 times
q0(U)+q0(V) in F. Since gcd(90,m)=d for odd m, choose its residue to
make F≡N mod m. This uses actual words. Counts are odd, so the
congruence is also mod 2m.

To place the integer, program one pair at every block j from floor(z0/2)
through z0. The first programmed x and y lie within O_m(q) of their
common midpoint a=qK/2+v. Its product gives a count strictly above N:
the choice of K gives a positive margin of order q sqrt(N), whereas
the product loss off the midpoint is only O_m(q²) and all remaining
terms in F are nonnegative. At the last block, y=O_m(q) and
x=O_m(sqrt(N)); F=O_m(q sqrt(N))<N. Bisect on the **actual** counts,
preserving opposite sides, without assuming monotonicity, until the
two block indices are adjacent.

Near their crossing write x=a+t,y=a−t. The common center term differs
from N by O_m(q sqrt(N)); the smooth dependence on t is
−Qt²+(T−S)t. The unmatched surcharge 45[p(U)+p(V)] lies in [0,45K].
Since one bracketing value is at least N, its |t| is
O_m(q^(1/2)N^(1/4)); the neighboring block changes t by O_m(q).
Thus adjacent counts at the crossing differ by at most
```
O_m(q^(3/2)N^(1/4) + q² + sqrt(N)/q) = O_m(N^(2/5)).
```
This proves witnesses on both sides, each congruent to N mod 2m.
For fixed m the algorithm uses only logarithmically many block probes
and bounded-length modular word programming, with no enumeration of
tournaments. It improves the congruent one-tape square-root estimate,
though rounding the high sum to a complementary cylinder gives a
weaker rate than the unconstrained 3/8 bracket.

The theorem is not uniform in the modulus. It does not allow choosing
m larger than the error while holding the construction cost fixed,
and it does not produce exact counts or tighter bounds at all dyadic
moduli. The next problem is an adaptive choice of the carry cell:
program just the needed scalar congruence with variable depth, instead
of forcing two free residues in every cylinder.

## Arbitrary modulus and a growing dyadic precision

The [paired parity extension](PAIR_PROGRAMMING.md#the-parity-coupling-pays-its-own-cost-all-moduli-for-placement)
lets the target retain the same O_M(N^(2/5)) bound modulo 2M for
**every fixed M**, not just odd M. If M=2^a m, require q divisible
by 2^(a+1); this holds automatically at sufficiently large N. The
non-checksum part then does not see U and V modulo 2^(a+1),
while the odd part is controlled exactly as above. Fix the common
low v using d=gcd(45,M)=gcd(45,m). After division by two of the
even congruence discrepancy, the adjustable q0(U)+q0(V) has
coefficient 45 modulo M, and the d test is sufficient to solve
for it. The sharper two-tape error and the full fixed-modulus
control therefore coexist, with real witnesses on both sides.

For pure dyadic M=2^a, the proof gives length ell=O((a+1)^2):
there is no odd loop and each reset word and the number of
available reset bits are O(a+1). This quantifies how far the
precision can grow before it overwhelms the real bracket.
Write D0=2^ell and take a=a(N)=o(sqrt(log N)). Then D0=N^{o(1)}.
With q of order N^(1/10), the center overshoot is O(qD0 sqrt N)
and an adjacent block changes t by O(qD0); the surcharge remains
O(sqrt(N)/q). Thus the proof gives **O(N^(2/5+o(1))) physical
error while the congruence holds modulo 2^(a+1)**. The initial
product margin and final below-target check still hold because
qD0=N^(1/10+o(1)) is o(sqrt N). This is a genuine increasing
precision result; it still does not reach a modulus larger than
the final error, which would force equality.

This also identifies the redesign needed to reach equality. For
pure dyadic moduli the all-syndrome compiler needs D0 at least M,
just to offer M choices in an aligned block. Exposing D0 in the
smooth-step budget gives the two competing terms q^(3/2)D0^(3/2)N^(1/4)
and sqrt(N)/q. Their optimum is of order N^(2/5)D0^(3/5), even if
q could be chosen freely. Making this **estimate** smaller than M
would require M of order at least N when D0 is at least M, already
outside this near-balanced proof. This is not a lower bound on
actual errors; it isolates the cost of programming *all* syndromes
in each block. To beat it, select only the needed syndrome along a
target-dependent carry cell or let a separate physical operation
pay the resulting error.
