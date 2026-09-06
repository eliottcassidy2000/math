# A prime divisor supplies complementary norm units at unbounded heights

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The result below supplies a new sufficient condition for nonvanishing of
the regular complementary norm. It therefore proves exact carried-boundary
multiplicities in an unbounded-height family. It does not prove the sign of
that norm, every odd-parameter unit, or general Laurent noncancellation.

## Inheritance and the coordinate retained

The closest proved mechanisms are the reflected regular factorial rows in
[continuing4 regular duality](continuing4_20260906_regular_duality.md) and the
all-height boundary-jet formula in
[continuing5 complementary norm jets](continuing5_20260906_complementary_norm_jets.md).
The latter makes the missing object precise: only the unit
`N_H(2r+1) != 0` is needed for the boundary order. The full positive
characteristic-coefficient bank at H<=6 is stronger than that consumer needs.

The canonical hostile is the lost inverse carry in the original row;
the boundary-jet theorem retains it before passing to the regular complement.
The corrected near miss is separate real-rootedness or a positive trace
being mistaken for cross-row coprimality. The least-used sidecar here is the
prime residue of the common upper factorial endpoint, rather than its real
size or a new symbolic height bank.

Source: the two rational factorial rows. Target: their actual resultant.
Map: reduction in the local ring Z_(p), followed by the resultant polynomial.
Preserved predicate: a nonzero resultant modulo p implies a nonzero rational
resultant. Destroyed information: real sign, root ordering, and the converse.
Required sidecar: the prime p, its valuation at the shared endpoint, all
factorial denominators, and monicity. The cheapest decisive hostile deletes
the valuation threshold at H=z=1,p=5; the doubled row then does not become a
monomial, although the rational resultant remains nonzero.

[Independent proof and exact audit](continuing6_20260906_norm_units_audit.md) passes.

## 1. Prime-divisor separation theorem

Let H>=1 and z>=1 be integers. Set S=z+2H and write the monic rows as

    Phi_(H,z)(t) = sum_(j=0)^H H! S! t^j /
                    [j! (z+2j)! (3H-3j)!],
    Psi_(H,z)(t) = Phi_(2H,2z)(t).

These are exactly the regular rows already proved in the reflected-family
dictionary, including the constant terms. Let

    N_H(z) = det(-multiplication by Psi in Q[t]/(Phi)).

Suppose p is a prime divisor of 2S-1, p>4H, and

    e = v_p(2S-1) > floor(6H/p).                       (1)

Then both rows belong to Z_(p)[t], every coefficient of Phi is a p-adic
unit, and every nonleading coefficient of Psi is in p Z_(p). More exactly,

    v_p([t^j]Phi) = 0                    (0<=j<=H),
    v_p([t^j]Psi) = e-floor((6H-3j)/p)   (0<=j<2H).

Consequently

    Psi mod p = t^(2H),
    N_H(z) mod p = (-1)^H Phi_(H,z)(0)^(2H) != 0.     (2)

In particular Phi and Psi are coprime over Q, hence have no common complex
root. The conclusion uses no real-root theorem and no primitive-charge
condition. The sign of N_H(z) is not asserted.

Because p>4H, condition (1) has the particularly simple sufficient forms

* p>6H and p divides 2z+4H-1; or
* p>4H and p^2 divides 2z+4H-1.

The second form is useful when 4H<p<=6H. No primality assumption on
2z+4H-1 itself is needed.

### Proof

The prime is odd. Since 2S is congruent to1 modulo p, the standard residue
of S is (p+1)/2. In the coefficient of t^j in Phi, the upper factorial
ratio is the product of the last 2(H-j) integers ending at S. It is a
nonempty product of residues strictly between0 and p, or an empty product:
its length is at most2H<(p+1)/2. No factor is divisible by p. The factorials
H!,j!, and (3H-3j)! are also units because 3H<p. This proves the first
valuation formula, including the constant term.

For j<2H in Psi, its upper factorial ratio is the product of the last
4H-2j integers ending at 2S. The interval has length at most4H<p, contains
2S-1, and contains no other multiple of p. Its valuation is therefore e.
The numerator factorial (2H)! and denominator j! are units. Finally,
6H-3j<3p/2<p^2, so Legendre's elementary factorial count gives exactly
floor((6H-3j)/p) for the last denominator valuation. This proves the second
formula. Hypothesis (1) makes it positive for every nonleading index.

The resultant is a polynomial with integer coefficients in the coefficients
of the two monic polynomials, so it commutes with reduction in Z_(p).
For a monic Phi of degree H,

    Res(Phi,t^(2H)) = Phi(0)^(2H).

The characteristic norm N is (-1)^H times that resultant. Its reduction
is nonzero because Phi(0) is a unit. This proves (2) and coprimality.

## 2. Exact carried-boundary multiplicities

For the carried family and norm c_(h,h)(x) in the preceding boundary-jet
theorem, take integers h>=1 and 1<=r<=h. Its complement has

    H=h-r, z=2r+1, S=2h+1, 2S-1=4h+1.

For H=0 the complementary norm is1. For H>=1 the same theorem proves that
the order of c_(h,h) at x=-r is exactly r-1 if and only if N_H(2r+1) is
nonzero. Thus the following sufficient arithmetic condition proves that
exact order:

    p | 4h+1, p>4(h-r),
    v_p(4h+1)>floor(6(h-r)/p)                          (3)

for at least one prime p. In particular, it suffices that the largest prime
factor of 4h+1 is greater than6(h-r).

If 4h+1 is prime, (3) applies at every

    ceil(h/3) <= r <= h,

or equivalently every complementary height 0<=H<=floor(2h/3).
This includes arbitrarily large complementary heights. For example h=13,
p=53 gives H=8,r=5,z=11, beyond the former H<=6 bank.

There are unboundedly many primes congruent to1 modulo4: if a finite list
contained all of them, every odd prime divisor of (2 times their product)^2+1
would again be1 modulo4 and outside the list. Here the residue of the
squared factor has order4, so4 divides that prime minus1. Hence the
prime-index consequence is an unbounded family, not a finite prime bank.

The prior theorem's full leading-jet formula, all factorial normalizations,
and inverse-carry treatment remain necessary dependencies. Only nonvanishing
has been supplied here. A residue being nonzero does not select the real
sign of the leading jet. The previously proved seven diagonals retain
their stronger real-sign certificates.

## 3. Scope, failure, and reproducibility

The producer checks every eligible prime in the explicit universe H=1..24,
odd z=1..101, and primes through301. This includes every prime divisor of
2z+4H-1 in that universe because that positive integer is at most297.
There are492 eligible parameter-prime triples,292 with H>=7, including four
cases where the square-prime extension is needed. Sixteen periodic controls
at H=7,12,20,31 show that z need not lie below p. The checked finite universe
supports the algebraic proof; it is not its unbounded justification.

Every coefficient valuation is checked exactly using rational factorial
coefficients. The producer then constructs the actual multiplication matrix
of the reduced doubled response and compares its determinant with (2).
Normal and optimized runs retain the same always-active gates.

The minimal threshold hostile is H=z=1,p=5,e=1. It has

    Phi=t+1, Psi=t^2+10t+1,
    Psi mod5=t^2+1, Psi(-1)=-8.

Deleting e>floor(6H/p) loses the monomial implication. The strongest survivor
is a nonzero rational response in this example; failure of this certificate
does not imply cancellation. A second hostile to a larger modular-unit
claim occurs at H=26,z=3,p=109: the two reductions share a linear factor.
The independent audit gives the exact reduced factor t-5 modulo109.
The same rows have N_26(3)=65 modulo157, proving rational coprimality.
As a separate FINITE-EXACT consequence, the carried boundary(h,r)=(27,1)
has order0. This shows the prime-divisor condition is sufficient rather
than necessary; it is not a common-root assertion over Q. The general unit for every odd z>=3 and H>=7 remains OPEN.

Reproduction:

    python -B 04-computation/continuing6_20260906_norm_units.py
    python -B -O 04-computation/continuing6_20260906_norm_units.py

The exact bank and norm residues are retained in
`continuing6_20260906_norm_units_certificate.json`. Raw source, output, and
certificate hashes belong in the continuing6 manifest; no byte-normalized
artifact is substituted for frozen evidence.
