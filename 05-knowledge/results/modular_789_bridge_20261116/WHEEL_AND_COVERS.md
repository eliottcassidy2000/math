# Two sharper obstructions: no map between the levels, no fixed finite cycle wheel

## Genus doubling does not even give a branched holomorphic map

There is no nonconstant holomorphic map X(9)→X(8). The genera
10 and 5 would force degree two, by Riemann–Hurwitz, and exactly
two simple branch points. The deck involution would therefore
have exactly two fixed points on X(9).

The full automorphism group of X(9) is the visible modular group
G_9, of order 324. Here is a direct argument in this case: Hurwitz's
84(g−1) bound is 756, so the visible subgroup has index at most
two. If index two, it is normal and the extra automorphism
descends to X(1), necessarily fixing the three branch values
whose inertia orders are distinct (2,3,9). It descends to the
identity and hence is already a deck transformation. Compare
[Bars, Kontogeorgis and Xarles, automorphisms of X(N)](https://mat.uab.cat/~francesc/xenemodular2.pdf).

Any involution with a fixed point is conjugate to S, the generator
of the order-two inertia: the other two branch orders are odd.
The projective centralizer of S has order 12. Indeed AS=SA and
AS=−SA give respectively matrices with determinant condition
a²+b²=1 and a²+b²=−1 modulo nine. Each equation has 12 solutions:
there are four solutions mod three with nonzero gradient, each
having three lifts mod nine. Quotient by ±I. Fixed cosets in the
fiber over the order-two branch value number 12/2=6.
So an involution has either zero or six fixed points, never two.
This rules out the proposed geometric doubling even allowing
branching. The different unramified level-72 cover from MODULAR.md
is not subject to this argument.

## Any fixed congruence wheel accepts arbitrarily long noninteger cycles

Here work with the positive shortcut on rationals of odd denominator.
Parity is the numerator parity in reduced form, agreeing with the
2-adic parity. Given any modulus M and real lower bound K≥1, there
are arbitrarily long primitive positive rational cycles entirely
above K with no integral vertex, whose denominators are coprime
to M. In particular their cyclic word equations and all states
reduce legitimately on the fixed wheel. This is an obstruction to
congruence-only rejection, not a claim that every condition imposed
on a wheel must accept them.

Let L be the unit exponent of the part of M coprime to six,
using L=1 if that part is one. Choose a≡1 and b≡0 modulo L and
take the word consisting of a ones followed by b−a zeroes.
Its fixed point for the positive sign is
```
D = 2^b−3^a, B = 3^a−2^a,
x = B/D = (1−(2/3)^a)/(2^(b−a log₂3)−1).
```
D is odd, prime to three, and equals −2 modulo the odd part of
M prime to six. Thus it is invertible on the entire wheel. The
affine inverse letters 2y and (2y−1)/3 land in the prescribed
even or odd class and each contracts 2-adically. Their composite
has a unique fixed point, namely x. In particular the prescribed
parity word is actually followed in the rational dynamics.

Choose arbitrarily large a=1+Lj with b a multiple of L so that
2^(b−a log₂3)−1 lies between 1/(4K) and 1/(2K).
This is possible by density of irrational rotation, since
log₂3 is irrational (unique prime factorization).
Then for large a, K<x<4K; the odd block rises, the even block
falls back to x, and x is the minimum. There is one cyclic run
of ones, so the word is primitive. For an integer with a
initial ones the equation T^a(x)+1=(3/2)^a(x+1) forces
x≡−1 mod2^a, contrary to 0<x<4K when a is large. No other
state can be an integer, since the forward orbit of an integer
would remain integral and return to x.

For the joint wheels 63 and 72 and K=20, take a=1465,b=2322:
floor(x)=47, D mod63=19, D mod72=37. A cycle proof that hopes
to reject all such candidates must retain something beyond
a predetermined finite modulus, for example a guard growing
with the word plus an exact integer-height estimate.

The fixed-wheel issue is complementary to the prime-seven
carry test. For a binary word of length b with a≥1 ones,
D=2^b−3^a is automatically a unit at 2 and 3, so its linear
period equation has one solution modulo every 2^u3^v, regardless
of the ordered carry. Mod7 may still exclude a word: 1100 has
D=7 and B=5. The parity-lift interpretation at the even prime
must still be kept; this is a periodic equation, not a
deterministic forward map on residues.

Run `python3 check_wheel_cover.py` for the exact matrix
and rational-cycle controls.
