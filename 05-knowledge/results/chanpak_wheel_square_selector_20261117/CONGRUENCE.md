# The protected-fringe congruence theorem

The earlier [signed two-flank law](../chanpak_decimal_core_selection_20261117/FLANKS.md) has
a fixed physical central tournament with marked profile
`(h,P,L,R,Lminus,Rminus,kappa)`. For positive left/right addresses
`u,w` encoded by their binary digits, the actual count is

    H = kappa + h(F(u)+F(w)-2) + R u + L w + P u w
              - Rminus [u odd] - Lminus [w odd].

Here F(0)=1, F(2x+1)=F(x), and F(2x)=F(x)+x+(x mod 2) for x>0.
Every F is odd by induction. The suffix of the right *address* is
taken before the physical right word's reverse-complement.

**Theorem.** Prescribe any binary prefix starting with 1 and any
nonempty binary suffix, independently for each flank. Their suffixes
fix parities e,f. Given any positive modulus M, put
`g=gcd(M,h,P,L,R)`. The counts in this protected-fringe family
attain exactly the odd residues (more precisely, residues of odd
integers modulo M) satisfying

    H = kappa - 2h - Rminus e - Lminus f  (mod g).

Realization is in both flanks, using actual binary words and the
unchanged physical core. Each allowed residue has arbitrarily long
witnesses, with two freely interleavable, distinct, equal-length
middle blocks available on either flank. This is a local theorem,
not a claim of exact placement or density in the integers.

## One binary address as a controllable state

For any modulus M, any prefix and suffix as above, the possible
pairs `(x mod M, F(x) mod M)` are precisely the pairs consistent
with (i) the suffix's forced low bits of x modulo the two-part
of M, and (ii) odd F. If the suffix is shorter than the two-part's
bit width, extend it on the left by the missing desired bits;
the proof now assumes the suffix has at least that width.

First work modulo the odd part m, with internal state odd x;
append a 1 to the prescribed prefix if it ends in 0. Set y=x+1,
z=F(x), and use maps in chronological digit order:
   
    A = 1   : (y,z) -> (2y,z)
    B = 01  : (y,z) -> (4y-2,z+y)
    T = 001 : (y,z) -> (8y-6,z+3y-2).

They are permutations of the finite state space over Z/mZ. An
inverse used algebraically below means a positive power of the
same physical word modulo m; it is not an inverse graph operation.
With composition written right to left,

    C = A^(-2) B   : (y,z) -> (y-1/2,z+y)
    Q = A^(-1) C A : (y,z) -> (y-1/4,z+2y)
    E = A^(-3) T   : (y,z) -> (y-3/4,z+3y-2)
    V = Q^(-3) E   : (y,z) -> (y,z-3y-1/2).

Indeed Q^3 has shift -3/4 and z-increment 6y-3/2; subtract its
readout at the same input as E. Go to y=0 with powers of C,
steer z there with V (unit increment -1/2), then go to the
desired y by C with its known z drift compensated beforehand.
No odd prime, including 3, obstructs the address state.

Write M=2^a m. A fixed suffix is invertible over the odd ring
on the pair (x,z): the previous bit is known at each step,
and the linear part for each appended bit has determinant 2.
So first choose the odd state that will give the desired pair
after the suffix. Let d=ord_m(2) (d=1 if m=1). Now append a
run of t ones with t>=a and d|t. This resets y to 0 modulo 2^a
and acts identically modulo m, without changing z. Choose
r>=a with d|(r+3), and put W=1^r001. Modulo m,

    W : (y,z) -> (y-6,z+3*2^r*y-2).

Consequently W^m is the identity modulo m: its z drift sums a
constant and an arithmetic progression, and m divides m(m-1)/2
because m is odd. Every copy of W adds exactly -2 to z modulo
2^a, since its initial ones reset y before 001. Thus W^m adds
-2m, generating every even dyadic increment. After the
appropriate number of repeats, append the t ones once more to
reset y without touching z. Now append the fixed suffix. Its
low bits fix the desired x, and the drift it imposes on z is
known from the reset y. This proves the simultaneous address claim
even for 2^a times an arbitrary odd-prime wheel.

**Repeatable liberty.** Before the final suffix, regard the
synced state y=0 modulo 2^a as protected. The blocks
`W^(m*2^max(a-1,0)) 1^t` and the all-ones word of the same
length both act identically on the odd state and on the protected
dyadic state. They are distinct, and inserting arbitrary words
in these two blocks makes an exact binary self-similarity in
*addresses*. It is not asserted that all readouts are distinct
at a fixed height, or that the sizes of the exact H gaps are
controlled by the modulus. Across growing lengths reads are
unbounded for the literal cores here (P>0 and both addresses positive).

## Projecting the protected state to the count

Necessity of the congruence is immediate. At each odd prime
power dividing M, the variable part is
`h(z1+z2)+R u+L w+Puw`, with all four coordinates independently
assignable. If h has minimum valuation, vary z1. Otherwise use
R with w=0 or L with u=0 when one of those has minimum valuation.
If P is strictly lower than h,R,L, set u=1 and vary w in the
coefficient L+P. If all vanish modulo the prime power the
image is just the zero ideal. This is exactly the ideal
generated by h,R,L,P. At the dyadic part, h is odd and z1
can range over the odd classes with both addresses otherwise
fixed. It follows that the actual count can hit every odd
class there. Combine by CRT; the address theorem makes this
a simultaneous physical construction, not independent local
witnesses that might conflict.

For example the triangle core with profile (3,12,3,6,-1,-2,5)
has g=gcd(M,3). Its four parity settings give 2,1,0,2 modulo
3 respectively in order (0,0),(1,0),(0,1),(1,1). Therefore
switching only those boundary choices gives *all odd residues
at every modulus*, already with the same three-vertex central
tournament. In contrast the earlier six-vertex core with
(45,432,129,258,-15,-30,137 or 139) has g=gcd(M,3) and neither
mask hits residue 0 modulo 3 for any suffixes.

This is a sharp warning for the universal selector: the triangle
does not hit equality H=N for all N just because it has no
local holes. In particular 7 and 21, known impossible as exact
tournament counts, both have protected-fringe congruence
witnesses modulo every M. Congruence filtering may eliminate
a core/parity, but no finite first-prime wheel can select an
exact witness when the local projection has already saturated.
What is missing is a size-controlled lift carrying a signed
absolute interval, not more residue states.
