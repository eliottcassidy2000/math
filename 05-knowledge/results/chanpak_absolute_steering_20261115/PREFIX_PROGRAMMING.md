# Arbitrary-high-prefix programming and congruent square-root placement

Use the known physical height-zero marked 45-core from the unit-slide
packet. Prepending a transitive tape of address x produces exactly

    F(x) = 2855 + 4034 x + 45 p(x),

where p(0)=0, p(2u+1)=p(u), p(2u)=p(u)+u+(u mod 2).
Use the nonzero address's actual binary expansion so the tournament is
strong. Write q=p/2, an integer. The exact recurrence for appending one
binary digit is

    1: (x,q) -> (2x+1,q)
    0: (x,q) -> (2x,q+(x+(x mod 2))/2).

## Two independent modular outputs in any prefix cylinder

**Theorem.** For each M>=1 there is a fixed suffix length L=L(M)
such that, starting with *any* high binary prefix, one can append
exactly L bits and obtain any pair (x mod 2M, q mod M).
In particular each aligned block of 2^L addresses contains a witness
for each pair. This uses actual digit appending, not independent
choice of two unattached arithmetic parameters.

Write M=2^a m with m odd. First use a fixed-length binary suffix
of length ceil(log2(2M)) to normalize x to -1 mod 2M:
choose its binary value as the least nonnegative representative
of -1 minus 2^length times the initial x. This necessarily ends in 1.

For the odd part, act on states with odd last bit. Blocks A=1 and
B=01 act modulo m by
    A(x,q)=(2x+1,q), B(x,q)=(4x+1,q+(x+1)/2).
Let D=B∘A^-2; thus D(x,q)=(x-2,q+(x+1)/8).
All inverses are physically positive repetitions: A has order ord_m(2),
and D^m is identity, since its q increment is m(x-m+2)/8 mod m.
For E=A∘D∘A^-1 the x displacement is -4 and q increment (x+1)/16;
D² also shifts x by -4 and has q increment x/4. Thus D^-2∘E
fixes x and feeds q by (1-3x)/16, a unit 1/4 at x=-1.
At most m-1 repeats prescribe q mod m with x=-1 mod m.
Pad with enough ones to restore x=-1 mod 2M; this leaves q unchanged.
If m=1 the odd stage is empty.

The dyadic repair needs to change q without changing either x residue
or q mod m. For j=0,...,a-1 construct a reset word of ones except for
two initial zero digits whose positions differ by d=j+2, followed
by additional isolated zeros. Positions are numbered from the right
at 0. Let K be a multiple of ord_m(2) at least a+2, and L0 a multiple
of ord_m(2) at least a+1 (use period 1 if m=1).
Use top zero positions L0+mK+d and L0+mK. Use r=(-1-2^d) mod m
further zero positions L0+iK, 0<=i<r. The complement address D0
is divisible by m (its contributions are 2^d,1,r mod m) and by 2^(a+1).
Appending this word therefore preserves x=-1 mod 2M, by
    x_new = 2^length (x+1) - 1 - D0.
At a zero digit, if R is the complement address of the word *before*
that digit, the contribution to q is -(R+(R mod 2))/2 mod 2^a.
This follows since x at that point is -1-R mod 2^(a+1). The first
zero contributes 0 and the second contributes -2^j (the prefix
defect there is 2^(d-1)); all additional zeros contribute 0 because
the least significant nonzero prefix defect has index at least a+1.
Repeat the reset word m times and its q change also vanishes mod m.
The resulting dyadic steps -m2^j form a binary digit palette modulo
2^a, because -m is a unit. Choose each at most once. Both coordinates
are now correct for x=-1 and any desired q mod M.

Finally, using a fixed-length suffix of length ceil(log2(2M)), append
the binary value congruent to the desired final x plus 2^length mod
2M. Its q drift is a deterministic function of x=-1 mod 2M and its
bits, so program q for that drift in the preceding stages. All variable
programming lengths can be made the same by padding with ones just
before this final suffix, where x=-1 mod 2M, so padding changes neither
coordinate. This proves the fixed-length claim including M=1. For
powers of two the sum of all reset lengths is O(a²). In general the
construction gives L=O(m² + (a+1)m²(a+m) + log M).

## Correct congruence and absolute square-root error at once

**Theorem.** For each fixed M, every sufficiently large odd N has
a physical strong tournament in this single fixed-core family, with
count congruent to N mod 2M and within O_M(sqrt(N)) of N. Moreover
there are witnesses on both sides of N at that distance.

To see this, append r>=1 ones to the programmable address u. Then
    x=2^r(u+1)-1, p(x)=p(u),
    F(x)=2855+4034(2^r(u+1)-1)+90q(u).
The congruence condition mod 2M, after division by two, is
an ordinary linear equation on (u mod 2M,q(u) mod M): the two
coefficients are 2017·2^r and 45, coprime since gcd(2017,45)=1.
It has a prescribed pair by Bezout, then a physical u anywhere in
any aligned interval of length D=2^L by the theorem above.

Let z be the real value of u giving F=N if q were zero. Choose the
first aligned block of length D starting at or above z for the
upper witness. For the lower witness choose the last aligned block
ending at or below z-45(z+1)/(4034·2^r). Within that block the omitted
surcharge 45p(u) is at most 45u<=45(z+1), since 0<=p(u)<=u. For the upper bound by induction, an odd v
obeys p(v)<= (v-1)/2, so the extra +1 at a step x=2v
with odd v is paid for; even v needs no such extra. Both witnesses are within O(4034·2^r D + 45(z+1)) of N.
Taking 2^r of order sqrt(N/D) proves the claim; for fixed M, N can be
taken large enough that all relevant u are positive.

At fixed M this combines congruence, relative scale, and a real physical
cut, keeping the remaining debt O(sqrt(N)), not O(N). It is not an exact
placement theorem: the h45 consecutive slide gives only k half-steps at
height k and changes its placement profile. Nor may we automatically make
M exceed the error: the all-pair programming certificate requires at
least 2M² distinct suffixes, and the proof's error is proportional to
sqrt(N·2^L). Its generic modular lift cannot catch its own error.
The next goal is an independently variable physical repair of this
scale, or a target-directed programmer that retains only the relevant
one-dimensional discrepancy, not two freely programmable coordinates.
