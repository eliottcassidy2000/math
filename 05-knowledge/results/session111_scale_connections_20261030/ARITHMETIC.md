# Gaussian doubling, aliquot sources and an Eisenstein obstruction

**PROVED elementary identities; FINITE-EXACT controls. Classical ingredients are not claimed as new.**

## Inheritance

The earlier source snapshot is [5ab37f308003e8634d13ea4c0dfe8f1558c11054](https://github.com/eliottcassidy2000/math/tree/5ab37f308003e8634d13ea4c0dfe8f1558c11054).
Relevant paths:
- 05-knowledge/results/collatz_sqdbl_20260925_squares_doubles_foundry.md.
- 01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md.
- 01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md.
- 01-canon/theorems/THM-4805-aliquot-sources-of-2-5-7-21-twin-prime-collisions-the-powers-and-doublings-chain-structure-and-abc-hits-on-the-row-of-one.md.
- 01-canon/theorems/THM-4826-which-grid-carries-a-norm-shell-one-to-two-split-cubic-character-of-two-centred-units-complete-the-quartic-designs-four-eisenstein-fourth-powers-miss-five-mod-eight.md.

MISTAKE-675 and MISTAKE-690 qualify prior-art and interpretation around the latter sources. We use only the elementary identities proved below. For classical Gaussian factorization background, see the [formalized Gaussian-integer development](https://isa-afp.org/entries/Gaussian_Integers.html).

## The exact transport

Let R(n)={(a,b):0≤a≤b, a²+b²=n}. The map
D(a,b)=(b−a,b+a)
is a bijection R(n)→R(2n), with D²=2I. The inverse on (u,v) is ((v−u)/2,(v+u)/2); equal parity follows from u²+v² being even. This is Gaussian multiplication by 1+i followed by a lattice symmetry.

The complete shells are:
- 25: (0,5),(3,4).
- 50: (5,5),(1,7).
- 100: (0,10),(6,8).
- 200: (10,10),(2,14).

For all j≥0,
R(25·2^(2j))={(0,5·2^j),(3·2^j,4·2^j)},
R(25·2^(2j+1))={(5·2^j,5·2^j),(2^j,7·2^j)}.

Axis and diagonal solutions interchange. Positive unequal solutions are preserved, but positivity alone is not. The hostile R(9)={(0,3)}, R(18)={(3,3)} rules out an unequal representation for every doubled square. For odd n, the content of D^k(a,b) is 2^floor(k/2) gcd(a,b), so primitivity also needs its scale retained.

A classical exact count sharpens the boundary. Unique factorization in Z[i] gives
r_2(2^k m²)=4 Π_(p≡1 mod4)(2v_p(m)+1)
for ordered signed representations. The ramified prime2 has no splitting choice; each split prime contributes its exponent choices. There is one special four-element axis or diagonal orbit and all other sign/order orbits have size 8. Consequently the number of positive unequal representations is
[Π_(p≡1 mod4)(2v_p(m)+1)−1]/2,
independent of k. Such a representation exists exactly when m has a prime divisor congruent to 1 modulo 4.

## Euclid and Berggren: preserve shape, retain scale

For Φ(a,b)=(b²−a²,2ab,a²+b²)=(A,B,C), direct expansion gives
Φ(D(a,b))=(2B,2A,2C).
Thus Gaussian doubling preserves primitive right-triangle shape, with a leg swap. It is not a new step down Berggren's tree.

For example, (3,4) lifts to (7,24,25), while (1,7) lifts to (48,14,50), which reduces to the same triangle with swapped legs. The separate inherited Berggren examples (3,4,5),(5,12,13),(15,8,17) explain the hypotenuse squares25,169,289. Their Gaussian transports are
169=5²+12² →338=7²+17²,
289=8²+15² →578=7²+23².
Taking a Berggren child and doubling a square-sum representation are distinct maps.

## Aliquot sums require an odd-part sidecar

For m odd,
σ(2^k m)=(2^(k+1)−1)σ(m),
s(2^k m)=2^k(2σ(m)−m)−σ(m).
This follows by writing each divisor uniquely as2^j d with d|m.

At m=25, σ(m)=31, so s(25·2^k)=37·2^k−31. The values are6,43,117,265. The inherited collision s(25)=s(6)=6 does not propagate by plain doubling. The exact correction is
s(2n)−2s(n)=σ(oddpart(n)).
The missing coordinate is the divisor sum of the odd part, not merely dyadic depth.

## Eisenstein hostile

For Q(a,b)=a²−ab+b², parity gives Q even iff a,b are both even. Therefore E(4n)=2E(n), but E(n) is empty whenever v_2(n) is odd.

Modulo5, Q=0 forces a=b=0 because the discriminant −3=2 is not a square. Hence E(25) consists of the six multiples by 5 of the units: (±5,0),(0,±5),(5,5),(−5,−5). The shell sizes along25,50,100,200 are6,0,6,0.

Gaussian doubling works because2 ramifies; the Eisenstein obstruction occurs because2 is inert. Coordinate dimension alone does not transfer the theorem.

## Verification

[arithmetic.js](arithmetic.js) independently enumerates Gaussian shells through n=10000, Eisenstein scaling through n=500, the aliquot law for odd m≤49 and k≤6, and 500 instances of the classical square-shell count. The universal conclusions rest on the proofs above, not finite extrapolation. No tournament enumeration or Collatz convergence follows from the norm identities.
