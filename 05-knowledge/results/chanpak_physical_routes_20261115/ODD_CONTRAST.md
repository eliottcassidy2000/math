# Short marked contrasts at the moving odd modulus 2^s−3

The [equal-Q route](EQUAL_Q.md) uses the step 6(2^s−3). The previous
modular compiler for fixed odd factors did not give a linear guard
for these moving factors. Here is a precise partial tool: a pair
of literal one-apex controls with unchanged **ordinary and reversed
addresses modulo m**, but whose Hamilton-count **difference** can
be programmed to any residue modulo m, where m=2^s−3. The whole
guard has O(s) sites. This is a programmable contrast, not an
absolute residue or small-step guarantee for an individual count.

Use bit positions from the least significant end. For s≥9 set

    A(z)=1−z^(2s−3)−z^(2s),
    C(z)=1+z^4+z^6,
    δ(z)=A(z)(z^(4s)−C(z)),
    δ_low(z)=−A(z) C(z).

The low block has degree 2s+6 and the high starts at 4s.
All its coefficients are 0,±1: the only close pair of offsets in
A differs by 3, whereas the differences of offsets in C are 2,4,6.
Thus δ is a literal change of bits, with +1 starting at 0
and −1 starting at 1. It is neutral to both address readings
modulo m: A(1/2)=0 mod m, because 2^(2s)≡9 and
2^(2s−3)≡9/8, and (2^(4s)−C(2))≡81−81=0 mod m.
The two evaluations detect both readings, since reversal only
multiplies the latter by a power of 2. In particular arbitrary
exterior padding in the two-cut seam changes this contrast by
zero modulo m (use fixed endpoint guards for parity).

The gap of positions j from 2s+8 through 4s−2 contains at least
s consecutive sites, none changed by δ. As a function of an
untoggled bit x_j there, the contrast H(after)−H(before) is
affine, with coefficient

    −2 Σ_i c(|j−i|) δ_i
      = −2^(j−1) δ_low(1/2) − 2^(−j−1) δ_high(2)
      ≡ 2^(−j−1) δ_low(2) mod m.

This uses the [one-apex count formula](../chanpak_graceful_ternary_20261115/GRACEFUL.md)
and the distance-at-least-two guard. Quadratic terms involving
only unchanged bits cancel from the contrast. The prefactor
δ_low(2)≡81·73/8 is a unit: 3 does not divide m, and the powers
of 2 modulo 73 are 1,2,4,8,16,32,64,55,37, never 3.
If the s consecutive bits are read from highest index down,
their coefficients form a unit times 1,2,…,2^(s−1) (successive
powers of two). Set all other bits, evaluate the base contrast,
divide the target difference modulo m by the unit, and write
its residue in binary across those s sites. This gives every
prescribed contrast residue. Adding a fixed first 1 and final
0 makes both witnesses strong.

The distinction is the same as the signed-count caution in the
[Hasse/Artin and Tunnell comparisons](NUMBER_THEORY.md): the
difference cancels unary terms the individual readout retains.
In particular two counts with equal residues can still be
separated by many moduli. This guard supplies a short exact
test for a new moving modulus, not an integer crossing on its
own.
