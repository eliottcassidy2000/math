# Joint absolute bracketing and exact residue for every fixed odd modulus

For every fixed odd modulus m>1 and all sufficiently large N,
the same six-vertex repair gives a full three-odd window on each
side of N, both within O_m(N^(1/3)). The two address pairs can have
the **same actual sum**, and the same designated center choice is
congruent to N modulo m on both sides. All choices of center still
give physical counts at steps +2. The result does not take m to
scale with N.

Use the [paired finite-state transfer](QRH_TRANSFER.md): at a
prescribed address sum S, one can replace a fixed-length suffix
with a complementary suffix that prescribes x and
w=p(x)+p(y) modulo m, exactly. The resulting x is odd and y is
even; the high prefixes and actual S do not change. The transfer
works for all odd m. Its earlier corollary only used m coprime to
30 because it inverted 45 directly. Here is how to handle the
remaining primes while keeping the bracketing geometry.

Choose j=2N mod 3 so that k_j=135+2j equals N mod 3.
For x odd, y even, at fixed S the readout is

    H_j = Q_S(x) + k_j + 45w,
    Q_S(x) = 432x(S−x) −129x +258S +15.

To prescribe w modulo m, the only necessary constraint on x and S
is Q_S(x)+k_j=N mod d, where d=gcd(45,m).
* Modulo 3 this is automatic.
* Modulo 9 the polynomial reduces to `−3x+6S+6`, so the required
  multiple of 3 can be chosen by setting x mod 3.
* Modulo 5 choose x=0 mod 5; then the polynomial is 3S. If 5
  divides m, impose `3S=N−k_j mod 5`.

The constraints on x can be combined by ordinary CRT. Having
chosen x, divide through by d and invert 45/d modulo m/d to set
w; the exact paired transfer realizes these choices at the literal
sum. Start with S≥ceil(sqrt(N/108)) and S≡−1 mod g with
g≈N^(1/6) a power of two as in [complementary bracketing](COMPLEMENT.md).
The optional mod-5 constraint moves S by at most another 4g.
In particular it retains S−sqrt(N/108)=O(g).

## Representatives on one fixed-sum grid

For fixed m, let the universal steering suffix consume B bits,
including one complementary bit that forces boundary parity.
B depends only on m and fits in the log₂(g)-bit complementary
band for sufficiently large N. Partition x into consecutive
blocks of length 2^B. Above the diagonal x≈S/2, choose one
steered representative from each block. The x-prefix is held
fixed within the block; y's suffix is complementary and its high
prefix fixed by S. Every representative has the same prescribed
residue in the same center j.

The first full block above the diagonal has imbalance bounded
in terms of m, so its whole triplet is above N for sufficiently
large N: the positive linear term is of order sqrt(N).
An imbalance O(sqrt(sqrt(N) g)) suffices to force the upper
envelope `432xy+129x+258y+45S+45+139` below N. Take a whole
block past the cutoff. Along this block sequence, the quadratic
changes between representatives by O_m(N^(1/3)): their address
gap is at most 2^(B+1). The uniform checksum remainder after
subtracting the common complementary main term is O(N^(1/3))
independently of their suffix choices. Thus consecutive counts
also differ by at most O_m(N^(1/3)). Bisect by sign on the lowest
and highest center branches to get a whole window on each side.
All these graphs are strong and have O(log N) vertices.
No monotonicity of the checksum and no claimed zeta zero-free
region are used.

## Short certified guard compiler

For small m, a direct finite-state dynamic program is much faster
than positive words obtained from group inverses. At current sum
s mod m, use state (x,w,e) with e the parity of x and w the
checksum-sum change. Appending b=0 or 1 to x (and 1−b to y) gives

    b=0: (2x,w+x+e,0)
    b=1: (2x+1,w+s−x+1−e,1)
    new sum = 2s+1.

Enumerate attainable states after exactly a proposed number of steps
from each initial residue pair and complementary parity. Keeping one
physical word per state is sufficient. For m=7,11,13 the least
universal lengths are 13,12,13 bits respectively (each final x-odd
state is reachable). For other fixed moduli the algebraic compiler
already proves some finite guard exists. These finite constants do
not claim uniformity in m.

Replay:

    python verify_short.py   # DP certification and 156 same-sum brackets
    python full_odd.py       # arithmetic and literal checks, including 3 and 5

The [capacity bound](CAPACITY.md) explains the separation from exact
placement: even every residue of every fixed modulus does not make
this prescribed geometric lane dense.
