# When even independent powers of 2 and 3 cannot make a return

One tempting repair is to decouple the two exponents of the
[synchronized port](../chanpak_physical_routes_20261115/SYNCHRONIZED_PORT.md):
suppose a finite bank of marked readers really did have independent
addresses a,b, with scalar

    F_i(a,b)=A_i·2^a+B_i·3^b+C_i,      A_i,B_i positive integers.

This is an idealized carrier model, **not** an assertion that
the independent addresses physically exist. Even granting them, a
finite bank of this form cannot have bounded-step transport deep in
both clocks unless only its finite offsets are changing.

## Bounded-step rigidity theorem

Fix the finite bank and a maximum step K. There is a threshold L
such that, if all four exponents in a move i,a,b→j,a',b' are
at least L and |F_j−F_i|≤K, then

    A_i 2^a = A_j 2^a',       B_i 3^b = B_j 3^b',
    F_j−F_i = C_j−C_i.

Thus any physical path within this model, staying deep in both
addresses and using bounded changes, preserves the two exponential
magnitudes. On a return to the same reader its net displacement
is exactly zero. Even independent control of the two exponents does
not by itself make an integer voltage.

**Proof.** Fix i,j and the integer step, from a finite list.
Subtract the offset difference to obtain a fixed integer d. The
equation on the four terms is

    A_j 2^a' − A_i 2^a + B_j 3^b' − B_i 3^b = d.

If d≠0 and no nonempty subsum on the left vanishes, divide by d
and use the finite-rank multiplicative-group theorem of
[Evertse–Schlickewei–Schmidt](https://annals.math.princeton.edu/2002/155-3/p04).
There are only finitely many such solutions. If a proper subsum
vanishes, it has size two or three. A zero triple leaves one
term equal to d and bounds its exponent. A mixed-base zero pair
fixes its two exponents: the ratio equation 2^a/3^b=fixed has at
most one solution. A same-base zero pair leaves the nonzero d
as a difference of two huge powers of the other base with fixed
integer coefficients. This is impossible once their common
power exceeds |d| and thus cannot divide it.

If d=0 and there is no proper vanishing subsum, divide by one
term and use the same theorem with three ratios. It gives
finitely many normalized lists. Each list determines at most
one absolute list: two solutions with the same ratios would
have a common scale factor equal both to a power of 2 and to
a power of 3, hence one. In the degenerate case the terms must
come in zero pairs. Mixed-base pairs are bounded as above.
Only cancellation by matching each same-base pair persists
deep in both exponents. Taking a threshold beyond the finite
exceptions, across the finite list of readers and possible
steps, proves the result.

## Role of this result

This follows the October [S-unit repair algebra at the −1 anchor](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-4640-the-repair-algebra-at-minus-one-collisions-are-increasing-two-three-representations-ladders-and-gersonides-carries.md):
distinguish *degenerate exact cancellations* from a large
congruence population. The new use is to classify bounded
transport in the ideal independent-clock replacement, not to
bound the Collatz collision spectrum. The finite list of offsets
can carry a small local strip (as the real tournament examples
do), but cannot carry a positive return at the same state.

The result does not address a modulus M growing with the height
unless a uniform bound on the exceptional solutions is separately
proved. It does not address changing coefficients, nonlinear
readouts involving cross terms, or moving down to the boundary
where one exponent stays small. Those are concrete escape routes
for Chan–Pak. No exact quantitative height cutoff is asserted
from the S-unit finiteness theorem alone.
