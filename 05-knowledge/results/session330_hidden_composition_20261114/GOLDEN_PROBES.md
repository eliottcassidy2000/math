# Golden depth locates an exact same-order tournament split

The eleven-square packet (eleven_squares_routing_bridges_20261007, §2)
retains a marked lift because equal unmarked tori can give unequal
clearances. Its Fibonacci matrix is literally useful as an arithmetic
coordinate test here, but no claim that a planar packing represents a
tournament is made.

## Two exact probes recover the hidden profile

An order-n tournament has H Hamiltonian paths and Q ordered optional
two-slot path covers. For a tournament W with profile (h,q), take the
literal cyclic substitution C3[U,W,1]. Its output is

    Y_W(U)=(q-h)(Q-H)+2hH=(q-h)Q+(3h-q)H.

This is the S329 two-reservoir law including both reservoirs split.
Two tests with probes W,V of the **same** order give outputs at the
same final order. Their 2-by-2 Jacobian in the abstract variables H,Q
has determinant 2(h_W q_V-h_V q_W). This is ordinary linear inversion,
not an invocation of the planar Jacobian conjecture.

For orders n≥3, every probe obeys 4|q and q ≡ (-1)^n h mod 3.
The latter follows from parity of the slot-polynomial exponents. So
any equal-order pair's probe determinant is divisible by 12. In fact
for n even and n≥4, 8|q: the empty packing contributes 2^n, and every
nonempty packing of c disjoint odd directed cycles contributes
2^c multiplied by 2 to a positive even power, since an individual
cycle of length l reduces the original exponent by l-1. Thus even
orders force determinant a multiple of 24.

Set z=(q-4(-1)^n h)/12. The normalisation is integral and the probe
determinant divided by 12 is the ordinary determinant in (h,z).

There are literal strong order-seven probes attaining determinant 12:

    outgoing rows                    (h,q)       (h,z)
    (2,12,9,49,39,71,31)             (61,944)     (61,99)
    (32,33,11,19,7,92,31)            (69,1068)    (69,112).

Indeed 61·112−69·99=1. Responses are Y1=883Q−761H and
Y2=999Q−861H; the exact inverse is

    H=(999Y1−883Y2)/24,   Q=(861Y1−761Y2)/24.

This has the smallest nonzero Jacobian determinant allowed for any
equal-order pair at orders ≥3. Simple transitive probes of smaller
orders are not counted in that claim. Order seven is the first strong
layer where the determinant 12 is attained: there is one strong
probe profile at orders three and four, all five at order five have
nonzero paired determinants of absolute value at least 24, and order six is excluded by the
even-order property. The replay obtains minima 24 at orders five and
six. The inverse gives arithmetic tests and a complete recovery of
(H,Q), not of U's isomorphism class or its full slot polynomial.

A one-channel alias really hides this state: both strong order-seven
sources with rows (2,4,81,7,43,15,59) and (2,4,65,23,39,15,59) give
the count 777 when W is a singleton, because H+Q=777 at profiles
(37,740) and (45,732). Thus equal size and equal scalar measurement
do not locate the hidden profile. The exact same-order two-probe map
always distinguishes these profiles. This corresponds narrowly to the
older double-torus Jacobian work (THM-3942): two detectors can have
independent coordinates yet a separate admissibility condition must
still be established. An inverse output passing integer congruences
is not a physically supplied graph.

## The profile that keeps both physical phase guards

The deep pair is not appropriate for every number-theoretic target:
the second probe has 3 dividing its Hamiltonian count, and the
inherited parity rule makes that entire channel divisible by 3.
In the seven-vertex census there is exactly one determinant-optimal
pair with neither count divisible by 3 (in fact both are coprime to 99):

    (10,36,33,22,7,89,31)            (89,1240)    (89,133)
    (16,17,11,35,108,71,15)          (91,1268)    (91,136).

Its normalised vectors are G²(45,44), G²(46,45), again adjacent at
their common golden inverse stage. In this case the physical outputs
are Z1=1151Q−973H and Z2=1177Q−995H, whose matrix has determinant 24
and exact inverse

    H=(1177Z1−1151Z2)/24,   Q=(995Z1−973Z2)/24.

For an even-order source with 3∤H, both outputs avoid 3, since
both probe orders are odd. The second output additionally has
Z2 ≡ 6H mod 11 because 1177=107·11. Thus even-order sources with
H coprime to 99 have a deterministic 99-unit second output on this
optimal two-coordinate chart, without acquiring another slot or
changing the output size relative to the first experiment. This
retains a usable phase, not an absolute prescribed target.

## Absolute-target constraint: the difference must be paid

Fix the first deep reading Y1 and ask whether changing a source of
the same order can alter the second reading without changing Y1.
Coprimality of 883 and 761 gives the exact rational-lattice direction
(ΔH,ΔQ)=(883t,761t), t in Z. When the source order is odd and at
least three, the parity and mod-3 conditions are exactly t ≡ 0
mod 4. When the source order is even and at least four, they force
t ≡ 0 mod 24, and these divisibilities are sufficient for the
linear congruences. The companion reading changes by −24t; its
smallest nonzero arithmetical jumps are 96 or 576, respectively.
No assertion that a physical source realizes those jumps is made.

There is no distinct source profile at the same first reading
near 30,000: for any source of order at least three the path-cut
bound gives Y1≥2771H. Two possible odd H values would differ by
at least 3532 just from 4|ΔQ, and the larger H would be at least
3533. In particular all such first readings below 2771·3533 admit
at most one profile, even allowing different source orders. A
phase overlap at this scale therefore cannot be implemented as a
nontrivial sliding source along this particular measurement line.
This does not exclude changing to a different cut or retaining a
higher slot coordinate concealed in the same profile.

## A precise Fibonacci coincidence, with a cutoff

Let G(h,z)=(z,h+z), the same Fibonacci matrix used in the marked
eleven-square packet. Call the depth of a positive profile vector the
number of successive inverse steps G^(-1)(h,z)=(z-h,h) for which both
coordinates remain strictly positive. These inverses test cone
membership; they are not physical tournament moves.

Of the 131 distinct strong order-seven profiles, exactly two have
depth seven and none have higher depth. They are exactly the two
determinant-optimal probes above. In fact their seventh inverse
images are the adjacent lattice vectors (6,1),(7,1), whose determinant
is -1. The odd power G^7 flips its sign to +1. Thus the continued
fraction chart gives a selection of a primitive experiment pair, not
merely two vaguely similar ratios. At orders five and six the maximum
depths are respectively two and zero.

Unimodularity supplies a literal integer-coordinate chart: every lattice
point in the closed cone on (61,99),(69,112) has a unique expansion in
these two rays with nonnegative integer coefficients. This is a useful
screen for proposed same-order moves. It does **not** imply that adding
the coefficient vectors corresponds to gluing tournaments; that is
exactly a missing physical realization gate. The two ratios are on the
same side of the golden ratio: their signed golden norms
z²−hz−h² are 41 and 55, respectively. The first and second seeds above
have norms −41 and −55, and each golden move reverses the norm sign.
Thus 55 tracks a literal quadratic-norm coordinate here; its identity
with the eleven-square packet's 55 is a coincidence in meaning.

There is also an all-order cutoff **for this golden lens**. Depth
at least seven implies (h,z)=G^7(a,b) with a,b>0. Its two column
slopes are 13/8 and 21/13, so

    21/13 < z/h < 13/8.

The path-cut injection gives Q≥(n+1)H. For odd order n the slope
z/h=(Q/H+4)/12 is therefore at least (n+5)/12; for even order it
is at least (n-3)/12. Comparing to the strict upper endpoint shows
that depth seven is impossible at odd orders beyond 13 or even
orders beyond 22. These necessary order cutoffs are not claims of
existence anywhere else.

More generally, a fixed Fibonacci recurrence, or any of its
rational affine conjugates on positive raw (H,Q) profiles, cannot
grow forever. The stable component contracts; if the expanding
component is zero the integer sequence eventually stabilizes. The
expanding eigenline is irrational, so a rational conjugacy cannot
make either of its coordinates zero. A genuinely growing orbit
would have H unbounded and Q/H tending to a finite constant. But
unbounded H forces unbounded minimum vertex order, and Q/H≥n+1
contradicts the finite limiting slope. To learn from the packing
transfer at arbitrary scale, use signed differences or a moving
frame, not this fixed cone on raw profiles.

**Scope.** The same-order split and readout are exact, as are the
odd/even lattice conditions and the depth cutoff. Recovering a marked
interface is not the Chan–Pak coverage or regeneration theorem.
