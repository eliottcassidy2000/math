# Matching the target phase without losing the one-third absolute bound

Combine the [one-third bracket](../arithmetic/SELECTOR.md) with the
[complementary-word controller](../bank/PROOF.md), on the original
balanced reflected seam. We prove:

**For every fixed positive modulus M and every sufficiently large odd
target N, there is a physical free-seam pair (H,H+2) with
H=N (mod M) and |H-N|=O_M(N^(1/3)).** The choices retain the exact
+2 switch. This is a simultaneous assertion, stronger than finding
unrelated witnesses for the bracket and the congruence. It does not
make the actual error zero.

The lower read for positive opposite-parity x,y is the symmetric
expression H=81+142(x+y)+256xy+5(F(x)+F(y)), physically putting the
odd address first. Write M=2^s m with m odd. We may increase s to 2
if necessary. Modulo m choose y=0 and use the Bezout identity
142(-2)+5(57)=1 to choose x=-2(N-81), Fsum=57(N-81). In particular
the desired sum x+y modulo m is now known.

## A constant-width complementary controller along the right lines

Let Q=x+y+1, Z=x+1, f=F(x)+F(y)-Q+Z. Before insertion start
from positive prefixes X odd, Y even with X+Y+1 divisible by 4.
Append a word to X and its bitwise complement to Y. One bit
always doubles Q, and f transformations for regulators A=1,
B=01,T=001 are the ones established in the companion proof.
The odd-modulus regulator words control Z and f independently
with Q fixed, and their lengths are multiples of d=ord_m 2
(take d=1 if m=1). They preserve f modulo 4 when Q is a multiple
of 4 and X is odd.

In fact **every word ending in 1** preserves f modulo 4 from
this starting state. A 0-bit changes f by -1 modulo 4 when
crossing from odd x to even x; another 0 contributes zero.
A 1-bit contributes +1 when returning to odd x, and otherwise
zero. This follows from the one-bit formulas
f'=f+2Z-Q+[x odd]-2 for a 0 and f'=f+[y odd] for a 1.
The net change to a parity-preserving word is zero.

Therefore we do not try to steer f mod 4: we steer the final
odd address's low two bits. Once the exact total S=x+y with
S=-1 (mod 2^s) has been fixed, for a proposed x_2=x (mod 2^s)
the required Fsum modulo 2^s is

    (N-81-142*S-256*x_2*(S-x_2))/5.

It is even, and its mod-4 value does not depend on x_2. Since
Q=0 mod 4, exactly one of x_2=1,3 makes the required
f=Fsum+Z agree with the initial f mod 4. Use that x_2 as the
fixed s-bit terminal word ending in 1. (There is no need to
keep x_2 independent of N or the prefixes.)

Take a long guard A^t with length divisible by d to zero
both Q,Z modulo 2^s without changing the odd state or f.
Compute the known influence of the terminal word from that
reset state; modulo m the terminal word is invertible.
Now steer the odd (Z,f) states to their needed preimages
with the odd regulator. At the dyadic part the preterminal
desired f agrees mod 4 as above. For r>=s and d|(r+3),
W=1^r001 changes f by -4 modulo 2^s; W^m is an odd-state
identity, and changes f by -4m on the dyadic state. These
moves select the higher dyadic bits. Finally use another
long odd-neutral A guard to reset Q,Z at the dyadic part.
All allowed operations ending in 1 have the expected
invariant mod 4. The argument also permits padding with A^d
inside the reset state to make the total length exactly
K=K(M), the same for all prefix/target data. All moving
parts take lengths 0 mod d except the fixed s-bit suffix,
so this padding can be arranged uniformly; the target odd Q
only depends on K mod d and is known once this residue is set.
This is the same finite-clock construction as the companion
proof with no initial toggle, and the replay has an explicit
upper bound for K. The word leaves the new total exactly

    S = D(X+Y)+(D-1),   D=2^K.

## Insert the controller after absolute bracketing

Choose dyadic L comparable to N^(1/6) with 4D and 2^s
dividing L (for fixed M and large N). Choose S as the first
integer at least ceil(sqrt(N/64)) with S=-1 (mod L) and
S=x_odd+y_odd (mod m) as prescribed above. The steps between
eligible lines are Lm. Now perform the same sign-crossing
argument as the one-third theorem on x+y=S. The center lies
above N, and moving left O_M(sqrt(SL)) in address units
brackets N. This follows directly from the quadratic product:
the line overshoot from its ideal center costs O_M(SL), and
the Fsum and linear terms cost O(S). Adjacent points in this
bracket differ by O_M(sqrt(SL)+L+S/L) by the complementary-carry
lemma. Binary sign bracketing picks an address pair within this
absolute error, without assuming monotonicity.

Pick a positive odd prefix X within O(1) of the approximate
x/D, and set Y=(S+1)/D-1-X; it is positive even, and
X+Y+1=0 mod4. Apply the controller. Its output is within O_M(1)
of the original address x, and y moves in the opposite direction
because the total remains *exactly* S. Thus only O_M(1) further
adjacent steps are added, all lying in the same carry-controlled
central band (up to a constant). Since S=Theta(sqrt(N)) and
L=Theta(N^(1/6)), the final physical lower member simultaneously
matches the full modulus and stays within O_M(N^(1/3)).

**Remaining obstruction.** The controller length and implicit
constant depend on M. Feeding it a modulus larger than the
observed error is circular: the allowed scale loss can grow
far faster than that modulus. A new quantitative or physical
repair would be needed to turn this into exact targeting.
