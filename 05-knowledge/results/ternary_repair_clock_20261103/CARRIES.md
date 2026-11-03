# The smallest closed integer correction module, and its observation quotients

Status: **PROVED** for integer linear combinations of corrections under
single-bridge source contexts. These minima are not classifications of
individual tournaments, joint attachment laws, or arbitrary re-markings.

## C1. Actual repairs force a three-dimensional integer module

The four realized repairs in THEORY R3 have delta=7kappa with

    kappa=(0,12,2),(0,10,2),(0,13,8),(0,14,5).         (C1)

They already span the entire plane 7{(0,x,y):x,y in Z} under integer linear
combinations. For their last two coordinates K_i, explicitly

    (1,0)=14K_1-18K_2+K_3,
    (0,1)=K_4-2K_1+10(1,0).

Closing this plane under q2 alone gives

    L7=7{(2u,x,y):u,x,y in Z}.                        (C2)

Indeed the first coordinates of (0,1,0)P_q2 and (0,0,1)P_q2 are 8 and 6.
Subtract their x,y parts, then take their difference, to obtain (2,0,0).
Conversely P[1,0] and P[2,0] are even, so this index-two lattice is invariant.

This invariance extends to **every actual supplied source block**, not only
the enumerated alphabet: the universal block compiler and (R3) give

    P[2,0]=h(c a)-h(c),
    P[1,0]=h(c aa)-h(c a),

both even by Hamilton-path parity. The same three calibrations give integer
entries in the full matrix. Thus arbitrary such blocks preserve L7. This
does not assert closure under an arbitrary new path complement or root
switch, which may require a different response carrier.

L7 is the smallest closed **integer module** containing these four realized
corrections. An arbitrary signed combination in L7 need not be a difference
of two individual legal source states. Realization guards are still needed.

## C2. Exact finite-modulus minima

Use the basis (14,0,0),(0,7,0),(0,0,7) of L7. Under repeated q2, the first
three H observations have matrix

    [[14,70,1022],[0,56,1008],[0,42,840]],
    Smith factors (14,14,336).                       (C3)

All later observations are integer combinations of these three, by
(X-1)(X^2-18X+5). Therefore, modulo any positive integer m, the exact
smallest additive observation quotient for q2 alone is

    (Z/(m/gcd(m,14)))^2 x Z/(m/gcd(m,336)).           (C4)

If all source contexts are permitted, empty/a/aa already give

    [[14,14,14],[0,14,14],[0,0,14]],
    Smith factors (14,14,14).                        (C5)

Every subsequent H response is 14 times an integer linear form, while these
three forms are unimodular after dividing by 14. Hence the exact quotient is

    (Z/(m/gcd(m,14)))^3.                             (C6)

For m=3, q2 alone has 9 observer values versus 27 with all source contexts.
For m=16 the comparison is 64 versus 512. For m=7 both see no carry at all.
The residue clock can therefore be exact and still miss later integer
responses. The minimum concerns arbitrary integer combinations in L7; the
reachable set of a particular nonlinear source grammar can be smaller.

## C3. A nonzero repair carry cannot remain in a finite error set

Let kappa=(0,r,s) be any nonzero integer vector. Then the sequence
kappa P_q2^k is unbounded; indeed its H observations have a nonzero
expanding component with root 9+2sqrt(19).

The eigenvalues are 1 and 9+/-2sqrt(19). The left 1-eigenspace is Q tau,
tau=(2,-1,0). If the expanding component of a rational vector vanished,
Galois conjugation would force its contracting component to vanish too.
The vector would then be a multiple of tau, impossible for a nonzero
vector of first coordinate 0. The H observation sees both quadratic
eigenspaces, as det O=96 in THEORY R1 proves. This establishes the claim.

For the four repairs, the first q2 step observes excesses over the ideal
affine rule of 756,644,1064,994, respectively. These are 7 times
108,92,152,142. Further steps obey the inherited cubic recurrence. Thus
a finite list of fixed error vectors is not closed under arbitrary q2
runs. This does not rule out compact formulas storing an initial carry
and a clock exponent, nor does it prove an information lower bound of k
bits for that special trajectory.

## C4. One additional exact offset transporter, with a singular residue phase

The complete chart audit through seven finds exactly three individual
symbols satisfying tau P=+/-tau:

    b: -,       q2: +,       s78: +.

Here s78 is an omitted six-vertex strong source chart, with

    P_s78=[[41,21,28],[80,43,56],[70,38,49]].           (C7)

Therefore the inherited two offset constructions extend to **every** prefix
over {b,q2,s78}, with sign (-1)^(number of b symbols). This is an infinite
family of actual strong sources, not a finite numerical coincidence. It
still fails to close at the 95/97 states by THEORY R1. The classification
is of individual letters, not all longer products returning tau to itself.

In integer coordinates (Y,X2,g) from THEORY R5, the new action is

    [[77,152,34],[28,55,12],[0,0,1]].                 (C8)

Modulo 7 its bulk has rank one: Y'=0, X2'=5Y+6X2. The full linear map has
a one-dimensional kernel generated in these coordinates by (1,5,4), or
by (5,4,1) in the original v coordinates. Reversible residue transport
therefore needs one F7 kernel coordinate; assigning a nonzero multiplicative
phase to every state would be false.

On the actual H31 seed, s78 gives

    v=(4141,2203,2877),  alpha=0,  g=4 mod 7.           (C9)

This graph is independently counted. It shows why the full clock must
include the zero field element. Subsequent q2 steps keep alpha zero,
although the integer source counts continue to grow. Frobenius, the
multiplication clock, and the integer carry are distinct operations/data.

## C5. What an all-size proof still needs

The carry formulas are exact evaluators of supplied legal replacements.
They do not produce a replacement when the target fibre is empty. To close
an induction, supply a source family and an algorithm that finds a legal
u=3v+j tau+7kappa for both digits from every family state, with u again in
the family. The four-level certificate verifies this only on its finite
tree. For the full odd spectrum one must additionally supply the 3H branch
and the seed range required by the inherited ternary induction target.

The new module identifies what an integer-linear continuation must retain;
the exact inverse and native commands identify what counts as realization.
Neither substitutes for the missing all-size source-existence theorem.
