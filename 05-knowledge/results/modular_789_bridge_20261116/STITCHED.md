# Even multiplication tables and the marked ±1 affine clock

Scope: an exact finite-ring comparison, not a return theorem for the
integer Collatz map or a construction of Hamilton path counts.

## An even table is an odd-style sign triangle with a marked parity edge

The [archived stitched-sheets note](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_spectrum_equality_and_stitched_sheets_20261005.md)
gives the odd center of the multiplication table, with rows and columns
indexed 0 through M−1. When M=2h+1, the indices h and h+1 are negatives;
their products read ±1/4 mod M. The two unreduced affine letters
A₋(x)=3x−1 and A₊(x)=3x+1 exchange the marks 1/4 and −1/4 when applied
in that order. The sign conjugacy x↦−x exchanges A₊ and A₋.

At M=2h the midpoint index h is fixed by negation, not paired with
a second index. The independent sign changes and the swap of row and
column always have the fundamental domain 0≤i≤j≤h. The new h+1
orbits on j=h, relative to the triangle for M=2h−1, supply the exact
channel
```
i h = 0 if i is even,    i h = h if i is odd    (mod 2h).
```
Indeed 2h=0 in this ring. When h is odd, h²=h mod 2h, and CRT
identifies h as the idempotent (1 mod 2, 0 mod h): the new edge
reads parity. When h is even, h²=0 mod 2h, so this nonzero midpoint
is square-zero, rather than a projection onto the whole dyadic
factor. This extends the archive's off-cross sign-triangle count
from the case of odd h to every even modulus.

One must distinguish ordinary negation from a parity-corrected
reflection. The alternating two-step map S=A₊∘A₋=9x−2 has a fixed
point exactly when 8x=2 mod M, equivalently gcd(8,M) divides 2.
For M odd its sole marked start is 1/4. For M=2h with h odd there
are exactly two starts, the two CRT lifts of 1/4 mod h; the first
letter sends each to h−x, not −x, since 4x−1 is an odd multiple
of h. Thus the even central cross marks precisely the carry for
the stitched reflection. When 4 divides M there is no two-step
stitched fixed point; longer cycles take its place.

## A marked affine coordinate gives the entire dyadic clock

For an odd integer a>1 let A₊,A₋ be the raw letters ax±1 and let
S_a=A₊∘A₋. Then
```
S_a(x)=a²x−(a−1),
(a+1)S_a(x)−1=a²((a+1)x−1).
```
The map x mod M ↦ z=(a+1)x−1 mod (a+1)M is bijective onto the
marked slice z≡−1 mod a+1. It intertwines S_a with multiplication
by a² on that slice, even when a+1 is not a unit mod M. In
particular z is always odd. For pure dyadic modulus M=2^s, S_a^j(x)
can equal x mod M exactly when
```
v₂(a^(2j)−1) ≥ s + v₂(a+1).
```
Since a²≡1 mod 8, iteration of the factorization of squares gives
v₂(a^(2j)−1)=v₂(a−1)+v₂(a+1)+v₂(j). Thus every point has the
same exact dyadic period 2^max(0,s−v₂(a−1)). This is a clock
in the marked state, not only a bound on the order of the multiplier.

For a=3 this also gives a simple full classification at
M=2^s 3^b u, gcd(u,6)=1. Write z=4x−1:
* Dyadic period 2^max(0,s−1).
* At 3^b, the coordinate enters zero after
  ceil(max(0,b−v₃(z))/2) steps (interpret v₃(0)=∞).
* The u-period is ord_[u/gcd(u,z)](9), with ord_1=1,
  since 4 is invertible mod u and multiplication by 9 is a unit.
* The complete preperiod is the middle number and the complete
  eventual period is the lcm of the two displayed periods.

Modulo 72 = 8·9 this means every x enters one of two period-four
cycles within one S-step: all periodic x satisfy x≡7 mod 9
and the mod-8 component has two 4-cycles. Under the sign-marked
one-step alternation each cycle doubles in length. This is an
exact use of the same 8·9 factorization as the level-six group:
the two cycles of S_3 are (7,61,43,25) and (16,70,52,34).
Reversing the sign order negates both cycles. This is not an
identification of the group action with either clock.

There is a literal distinction between 7=8−1 and 9=8+1.
Their square multipliers are both 1 mod 8, but the marked
translation yields periods 4 and 1 respectively. Modulo 72,
S_7 has six 4-cycles and four 12-cycles, all periodic: modulo
9 there are three fixed points and two 3-cycles, which combine
with the two 4-cycles modulo 8. S_9 has eight fixed points,
entered after at most one step. Meanwhile 3 is the only odd a>1
for which a²−1 is a power of two: a−1 and a+1 would both
be powers of two differing by two, hence 2 and 4. This is the
elementary square/dyadic corner of the Catalan theme, giving
exactly 9−8=1 without importing Mihăilescu's theorem.

The replay's all-modulus test compares this arithmetic with
literal orbits. No modular period is a real return:
S_a^j(x)=x as rationals would force z=0, or x=1/(a+1),
not an integer.

## What most quickly breaks a transfer to Collatz

For a fixed sign, actual shortcut Collatz sends even x to x/2
and odd x to (3x±1)/2. It does not descend to a deterministic
map modulo any M>1. For even M choose x and x+M even: their
outputs differ by M/2, despite the same branch and residue.
For odd M, an even x and x+M use different branches. For their
successors to agree mod M one would need 2x±1=0 mod M for every
even x, impossible as even x covers all residues. The missing
input is a dyadic lift, not only a residue or a sign marker:
even M already fixes the sign of the parity branch and still
requires one extra input bit per shortcut division to read an
output mod M. The odd affine letters and the actual divisions
must not be interchanged.

The [archived Carmichael carry synthesis](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/carmichael_collatz_carry_20261019/INHERITANCE_AND_BRIDGES.md)
already distinguishes full modular-affine return from integer
descent; its 7→11→17 mod-5 return grows on the integers.
The finite periods above are therefore a usable residue screen for
the marked alternating word, not a proof that every admissible
absolute target has a tournament witness. The odd unit-period
factor is bounded by the exponent of the unit group; that
Carmichael-function fact is not Carmichael's totient conjecture.
In fact λ(72)=6 and ord_72(49)=3, yet S_7 has period-twelve
orbits: S_7³(x)=x−18 mod 72 is an order-four translation.
The marked coordinate has conductor (7+1)·72=576, and at that
conductor ord_576(49)=12. This is an exact example of why a
multiplier's order at the original modulus need not control the
affine composition's order.
