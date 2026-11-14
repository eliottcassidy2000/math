# Exact inverse-Collatz rays in the tournament binary chart

This is an equality of physical word recursions, not a reduction
between the Collatz conjecture and Chan–Pak.

## Parity chooses the inverse phase, then the word updates

A normalized strong one-apex/transitive-spine word starts in 1
and ends in 0. Its binary address P has unique form 2^s x for
s≥1 and x odd positive. The exact outputs are

    H_s(x)=g(x)+(2^s−1)x+1,       R_s(x)=2^s x+1.

This follows by iterating g(2z)=g(z)+z+(z mod 2) and
g(2z+1)=g(z). Now write 3x+1=2^k y with y odd. Necessarily
3∤y. Conversely an odd target y prime to 3 has a positive odd
inverse x=(2^k y−1)/3 exactly for positive even k if y≡1 mod 3,
or positive odd k if y≡−1. The parity requirement comes from
2^k≡(−1)^k modulo 3. Raising k by two gives x→4x+1, i.e. append
bits 01 to x. For odd x, g(4x+1)=g(x)+x+1. The literal
tournament transport is therefore

    H_s(4x+1)−H_s(x) = (3·2^s−2)x+2^s,
    R_s(4x+1) = 4R_s(x)+2^s−3.

In particular s=1 gives (H,R)→(H+2R,4R−1). Each normalized
word belongs to one such inverse-Collatz ray indexed by (y,s).
No forward convergence or randomness assumption is used.

For fixed base x₀ and s, let t count steps down the inverse ray.
Put A=(3·2^s−2)(3x₀+1) and C=9H₀−A. Then

    9H_t = A·4^t + 6t + C.

Indeed x_t=((3x₀+1)4^t−1)/3 and one sums the exact increment.
Here A≡1 modulo 3, A>0 and C odd. The sequence increases
strictly, so an absolute target on a specified ray has an exact
monotone inverse: binary search t and insist on equality.

## The two local modes lift jointly at every precision

For odd N modulo 2^b, take t≥ceil(b/2). Then H_t≡N iff

    3t ≡ (9N−C)/2  (mod 2^(b−1)).

Every odd residue appears. For the ternary test put α=(A+2)/3.
By expanding (1+3)^t modulo 27, the three output residues
relative to H₀ at t=0,1,2 are 0, α, 2α+1 modulo 3.
For any α these form exactly two residues.

Both residues lift to all their 3-adic digits. For j≥1 write
u_j=(4^(3^j)−1)/3^(j+1). The base u₁=7 and cubing gives
u_(j+1)=u_j+3^(j+1)u_j²+3^(2j+1)u_j³, so u_j≡7 mod 9.
Consequently the exact increment gives

    H_(t+z·3^j)−H_t ≡ z·3^j(α−1+t)   (mod 3^(j+1)).

For clarity, before division by 9 the difference is congruent
to 3^(j+1)z(A·4^t u_j+2) modulo 3^(j+3); the bracket is
divisible by 3 and its quotient gives the displayed coefficient.
The singular class t≡1−α mod 3 is in the repeated-output pair,
not the unique-output class t≡2−α. Thus each of the two output
classes has a nonsingular preimage, and its remaining digits
lift uniquely from that preimage. The CRT combines this one
t mod 3^a with the dyadic constraint; choosing a large enough
representative makes the witness deep as required.

In fact just **two fixed physical rays** are jointly complete
on odd residues modulo every 2^b·3^a. Use x₀=1 (odd y=1)
and s=1 or s=2. The first has H₀=3 and (A,C)=(16,11), with
ternary image {0,1}. The second has H₀=5 and (A,C)=(40,5),
with ternary image {1,2}. The CRT choice is made on one of
these two actual rays, not by combining different graphs
that separately pass the two congruences.

## A fixed admissible target with exponentially expensive matches

For an odd fixed N not on a given ray, matches at all large
dyadic precisions are deep. Otherwise 4^t divides
6t+C−9N, impossible for unbounded t since the latter grows
only linearly (unless it vanishes at a fixed t). Bounded t
can match at arbitrarily large precision only at equality.
In the deep case, with M=2^b,

    3t = (9N−C)/2 + j·M/2.

For sufficiently large b, j≥1 because t≥ceil(b/2), so
t≥M/6−O(1) and H≥exp(Ω(M)). This is a scalar obstruction
for a fixed ray, not the pinned-three-response bound
from the earlier Sun-inspired construction.

The same separation is explicit at **N=63**, known to be
attained in other tournaments. The two given rays pass through
(3,9,31,117,459) and (5,19,73,287,1141) at t=0,…,4,
so neither contains 63. The listed terms never match 63
modulo 2^b for b≥6. When t≥5, |6t+C−567|<4^t and the
difference is nonzero, excluding a shallow match. Thus the
first matching t in the first ray, for b≥6, is

    t = (278+j·2^(b−1))/3,
    j=1 when b is odd, and j=2 when b is even.

The second ray's first t is exactly one larger (replace 278
by 281). Hence the first ray gives the least count across
the two, exactly (16·4^t+6t+11)/9. The least error is genuinely
exponential in the dyadic modulus, even for this one prescribed
admissible N and even though the union is jointly residue-complete
for 2 and 3. Higher precision cannot turn these two physical
rays into a universal exact selector.

This mirrors the arithmetic role of parity in inverse Collatz:
it selects a ternary mode, while an independent ordinary-size
argument would still be needed to find a small witness.
For constructive height windows of a different kind see
[TIGHT_COMPLEMENT.md](TIGHT_COMPLEMENT.md).

Run **python3 replay_collatz_runs.py** for recurrence, Hensel,
modular CRT and literal small tournament controls.
