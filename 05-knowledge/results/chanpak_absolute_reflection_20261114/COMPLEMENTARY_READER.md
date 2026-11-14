# Complementary right readers with the exact same ordinary count

There are literal strong one-apex tournaments with two complementary
right readers, one of each minimal width, and the same ordinary count
but different positive rooted responses. This holds in a two-parameter
family and the equality within its reader rail has a complete decoder.
It is a positive exact return, not coverage of an arbitrary odd target.

Use the binary-spine function g(0)=1 and
g(2x+1)=g(x), g(2x)=g(x)+x+(x mod 2). For an n-bit word w,
allowing leading zeroes, write rev_n(w) for reversal, Hbit and Lbit
for its two boundary bits, and bar(w)=2^n-1-w. Then

    g(w)-g(bar(w)) = w-rev_n(w)+Hbit-Lbit.

Indeed g-1 counts binary pairs of the form 10 with weight 2 at
distance one and 2^(d-1) at distance d≥2, by unwinding the recurrence.
The complement difference replaces each 10/01 indicator difference
with the difference of the two bits themselves. The running sums of
the pair weights through distance d are 2^d for d≥1 and zero for d=0.
This gives the stated coefficients, including the two boundary bits.

Choose one central singleton which beats the apex, a left transitive
wing of positive address x and a right transitive wing of positive
reverse-complement address y in minimal width. The whole spine is
transitive. Its contact word from the apex is x in binary, then 0,
then the reverse-complement of y in binary. Its count and rooted
response are (from the bilateral formulas, or the one-apex recurrence)

    H(x,y)=g(x)+g(y)-1+(2y+1)x+(x mod 2),
    R(x,y)=(2x+1)2^bitlen(y)-rev_bitlen(y)(y).

For fixed k≥2 and n≥k+2 consider complementary right readers
y=2^(n-1)+2^k-1 and z=2^(n-1)-2^k (the latter has minimal width n-1).
When computing the complement difference use width n for both, which
does not affect g. Thus y-z=2^(k+1)-1 and

    g(y)-g(z)=2^k-2-2^(n-1)+2^(n-k).

The two ordinary counts agree exactly iff

    x=(2^(k-1)-1)(2^(n-k-1)-1)/(2^(k+1)-1).

This is a positive integer iff n=(m+1)(k+1) for some m≥1.
Indeed gcd(2^(k-1)-1,2^(k+1)-1) equals 1 or 3. The multiplicative
order of 2 modulo the reduced denominator is k+1: it divides k+1,
and a proper divisor is at most (k+1)/2, whereas the reduced
denominator exceeds 2^((k+1)/2)-1 (for parity even k the denominator
is intact; the smallest divided case is k=3). Hence precisely multiples
of k+1 work. In that case x is odd, and the absolute results are

    H(x,y)=H(x,z)=g(x)+2^n x+2^(n-2)+2-2^(k-1),
    R(x,y)=2^(n+1)x+2^(n-k)-1,
    R(x,z)=2^n x+2^(n-1)-2^(n-1-k)+1.

Both contact words begin with 1 and end with 0, making both literal
tournaments strong. The reader widths differ by one, so the tournaments are not
isomorphic even without the root. At fixed k and increasing n,
H is of order 4^n while g(x) and the additive terms are O(2^n).
The two exact rooted responses are therefore 2H+O(sqrt(H)) and
H+O(sqrt(H)), with zero change to the ordinary count. For k=2 the initial examples have H=81,4741,300065;
for k=3,m=1 there is H=831 with responses 1567 and 881; for k=4,m=1,
H=7419. To decide whether prescribed H belongs to this construction,
test only n≤floor(log2 H), each k+1≥3 dividing n with quotient≥2,
and compare the absolute readout. No congruence substitutes for
this test, and the decoder does not guarantee acceptance.

Run replay_complement.py for exhaustive complement identities through
width 14, 88 exact rail pairs, and four literal subset-DP pairs.
