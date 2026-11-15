# Exact rational twists on the cheapest bottleneck-free floor chains

Wave 28 (THM-5202) identifies the alternating span-two ribbon floor units
A=0.2 and B=1.2. Wave 29 (THM-5212(A)) proves that these are the only
rank-at-least-two words in the {A,B} unit language, and identifies their
two nonzero characteristic roots. We use these results here. This note
concerns averages **over the depth of a chain**, not a probability law
with branching at each depth.

Let h_j=H(0(AB)^j) and g_j=H(0(BA)^j). At j≥1 both obey
`h_(j+2)=9h_(j+1)+8h_j`; the two initial pairs are (9,89) and
(5,53). They follow from the integer count coordinates of THM-5212(A);
the zero eigenvalue in the three-coordinate model is discarded after
the first pair, not at j=0.

For X=h or g, let M_X(a/q) be the limiting complete mean of
exp(2πia X_j/q), averaging over j. It exists since the sequence is
eventually periodic mod q. For gcd(a,q)=1 the following hold.

| Conductor | h | g |
|---|---|---|
| 3 | 1/4 | −1/2 |
| 5 | −1/24 | −1/24 |
| 8 | exp(2πia/8) | exp(2πi5a/8) |
| 3^r | zero for r≥3 | zero for r≥4 |
| 2^r | zero for r≥5 | zero for r≥5 |
| 5^r, r≥2 | magnitude at most (1/6)·5^(−(r−2)/2) | same |

In particular, a balanced mixture of the two families cancels the
primitive 1/8 twist exactly, while each individual family has modulus
one there. Neither family branches at all. Thus finite or even
high-conductor modular cancellation along size cannot be identified
with archimedean spectral density, and finite-denominator twists may
survive in a floor language even once all its chain bottlenecks are
removed. The number of chain values at most X is O(log X).

## Proof

Let T=[[0,1],[8,9]] act on the column (X_j,X_(j+1)).
Modulo 3 its period on both initial pairs is 4, with outputs
(0,2,0,1) and (2,2,1,1). Modulo 5 both state orbits have period 24 and
cover all nonzero vectors; the first-coordinate histogram is (4,5,5,5,5).
Modulo 8 the recurrence is constant. These three finite statements can
all be checked by multiplication, as in the verifier.

**Dyadic conductor.** The two 2-adic roots α,β of
f(x)=x²−9x−8 lift its two simple roots mod 2. The unit root α satisfies
v₂(α−1)=4 because f(1)=−16 and f(α)−f(1)=(α−1)(α−8) with α−8 odd.
The other root satisfies v₂(β)=3 by αβ=−8. Therefore
X_j=Cα^(j−1)+Dβ^(j−1), where C,D∈Z₂ since α−β is odd, and C is
odd since X_1 is odd. For j≥1+ceil(r/3) only the unit term
survives mod 2^r. The elementary identity
v₂(α^(2^t)−1)=4+t (successive squaring) shows that the unit term
is a bijection over a complete period of length 2^(r−4) on one
coset mod 16. Summing a primitive additive character over that coset
gives zero for r≥5.

**Ternary conductor.** T^4=I+9D with D=[[79,97],[776,952]].
Split j into the four orbit phases, then write the readout as
F(t)=e₁ᵀ(I+9D)^t v for t in Z₃ and the corresponding starting
phase vector v. Set s=v₃(e₁ᵀDv). For h all phases have s=0;
for g the successive phases have s=(0,1,0,1). These are exactly
verified modulo 9 with the four starting vectors.
For s≤1, (F(t)−F(0))/3^(2+s) is a 3-adic isometry:
its leading linear coefficient is a unit, while every remaining
term is a strict 3-adic contraction. To see this directly use
L=log(I+9D), with L=9D mod 3^4, and expand exp(tL); the higher
coefficients L^i/i! have valuation ≥4 for i≥2, and integral
monomials are 1-Lipschitz. Lifting digits gives a bijection modulo
each 3^k. The phase output is therefore uniform on a coset modulo
3^(2+s); this proves the stated cancellations. More precisely the
h cosets modulo 9 are 0,8,0,1; the g cosets are 5 mod 9, 26 mod 27,
4 mod 9, and 1 mod 27.

**Quinary conductor.** Write T^24=I+5E, where
E≡[[3,2],[1,1]] mod 5. Again split into the 24 orbit phases with
vectors v. Exactly four satisfy e₁ᵀEv=0 mod 5: the nonzero
vectors on a line. On any other phase, t↦e₁ᵀT^(24t)v is
uniform on a coset mod 5, by the same isometry argument, so
contributes zero at r≥2. On a critical phase define
F(t)=(e₁ᵀT^(24t)v−e₁ᵀv)/25. Write L=log(I+5E), so L≡5E mod 25.
The linear coefficient of F is integral by the critical condition;
the quadratic coefficient is e₁ᵀE²v/2 mod 5 and is a unit:
the first rows (3,2) and (1,3) of E and E² are independent mod 5.
Every higher coefficient is divisible by 5, since k−v₅(k!)≥3
for k≥3. In particular F' is an isometry, by reduction to its unit
linear term modulo 5 and the contraction of all other coefficients.

The unique zero t₀ of F' exists by its digit-lifting isometry.
With a=F''(t₀)/2 a unit, write
F(t₀+z)=F(t₀)+az²U(z), U≡1 mod 5 coefficientwise. The square root
V=√U is 1 modulo 5 coefficientwise, by the convergent binomial
series. The change z↦zV(z) is an isometry and gives a bijection
mod every 5^R, so the complete primitive sum of this phase has
normalized modulus exactly 5^(−R/2), R=r−2. Indeed the square of
the quadratic Gauss sum is evaluated by the invertible change
(u,v)↦(u−v,u+v) and additive orthogonality. Weighting the critical
phases by 4/24 yields the bound. The log/exp series used above
converge entrywise on these respective prime-adic matrix balls.


None of these cancellations gives a polynomial bound on collisions
at a prescribed absolute count. In the rank-two chains the contrast
is exact: size-indexed congruence equidistribution can coexist with
only logarithmically many integer values. A density argument still
needs many genuinely distinct branches at comparable absolute height.

## Mixed rational phases

The conclusions persist with explicitly bounded losses when the
conductor is 2^a·3^b·5^c (no other primes) and the numerator is coprime
to it. Sufficient exact-vanishing conditions are:
- a≥8, for either chain;
- b≥4 for h or b≥5 for g.

If c≥2, the bound (1/6)·5^(−(c−2)/2) holds independently of a,b.

For the proof let j run over a common multiple of all modular periods,
after the dyadic transient. Periods can be taken to divide respectively
2^max(a−4,0), 4·3^max(b−2,0), and 24·5^max(c−1,0). Freeze j mod 24
(or include a harmless factor 24 in the averaging period). The three
remaining primary components of the depth index are independent by the
Chinese remainder theorem. The dyadic coordinate sees the congruence
j mod 8, which replaces its coset mod 16 by one mod 128, giving a≥8.
The ternary coordinate sees j mod 3, multiplying a phase's slope scale
by 3, hence the stated thresholds. The quinary coordinate has no extra
5 restriction from the other two periods, so its phase estimate is
unchanged. Exponential phases split by CRT across the three prime powers;
their numerators remain units. Finally, for a uniformly random phase
mod 24 the critical quinary weight remains 4/24. These bounds explain
why the small 3-and-5 major arcs can be present while deep smooth-modulus
cancellation survives; they do not constrain denominators involving
other primes.

## Exact divisibility reappears at the depth index

Let U_0=0,U_1=1 and U_(j+2)=9U_(j+1)+8U_j. Then
h_j=U_(j+1), whereas g_j=U_(j+1)−4U_j. In particular, any prime
count on the AB chain must have j+1 prime. For a proper divisor d>1
of j+1=dm, put α,β equal to the two (real) characteristic roots.
The quotient U_(dm)/U_d is the symmetric polynomial
Σ_(i=0)^(m−1)(α^d)^(m−1−i)(β^d)^i, hence an integer
(its elementary symmetric inputs α^d+β^d and (−8)^d are integral).
Both U_d>1 and U_(dm)>U_d follow from the positive recurrence.
Thus chain-level rank-two irreducibility does not prevent an
arithmetic factorization along periodically spaced depth indices.
By the standard upper bound for the prime-counting function, the
AB chain contributes at most O(log X / log log X) prime counts up
to X. No corresponding primality criterion is asserted for BA.
