# A nearby factorial endpoint reduces the norm to a smaller quotient

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**

[Independent proof and exact audit](continuing7_20260906_norm_rank_drop_audit.md) passes.

## Inheritance and the retained consumer

The nearest proved mechanism is
[continuing6 prime-divisor units](continuing6_20260906_norm_units.md).
Its monomial reduction supplies actual complementary norm units, consumed by
[continuing5 boundary jets](continuing5_20260906_complementary_norm_jets.md).
The canonical hostile is a zero modular resultant that does not imply a zero
rational resultant. The corrected near miss is treating individual row root
geometry as cross-row coprimality. The underused operation is a nearby odd
factorial endpoint: it gives a small residual quotient instead of a monomial.

Source: the same two rational factorial rows. Target: their characteristic
norm. Map: reduce at an odd upper-endpoint residue and factor out the zero-root
block of the doubled row. Preserved predicate: the exact resultant modulo p.
Destroyed information: real sign and the converse from a bad prime.
Required sidecars: H,m,p, the constant unit, and the residual determinant.
The cheap hostile at H=2,p=9601 shows why the last sidecar cannot be dropped.

## 1. Exact rank reduction at every short odd endpoint

Let H,z be positive integers, S=z+2H, and

    Phi(t)=sum_(k=0)^H (H)_k (S)_(2k) t^(H-k)/(3k)!,
    Psi(t)=Phi_(2H,2z)(t),
    N_H(z)=(-1)^H Res(Phi,Psi).

All falling products and the norm convention agree with the inherited rows.
Let p>6H be prime and 0<=m<2H an integer such that

    2S == 2m+1 (mod p).

Define the monic degree-m polynomial over F_p

    R_(H,m)(t)=sum_(k=0)^m
                (2H)_k (2m+1)_(2k) t^(m-k)/(3k)!.

Then every coefficient of Phi is a p-adic unit, both rows are p-integral,
and the entire doubled row satisfies

    Psi mod p = t^(2H-m) R_(H,m)(t).                 (1)

Consequently, with bars denoting reduction,

    N_H(z) == (-1)^H Phi_bar(0)^(2H-m)
                        Res(R_(H,m),Phi_bar) (mod p). (2)

The resultant on the right is the determinant of multiplication by Phi_bar
in the degree-m quotient F_p[t]/(R). For m=0 the empty determinant is1.
Thus a degree-H response determinant can be tested in a degree-m quotient.
The quotient is smaller when m<H; the identity also holds for H<=m<2H.
The modular norm is nonzero exactly when this residual determinant is nonzero,
which sufficiently proves rational nonvanishing. Rank reduction alone is not
a nonvanishing theorem.

### Proof

The standard residue of S is (p+2m+1)/2: the alternative half-integer is not
an integer, and 2m+1<p. This residue exceeds3H and is less than p. The last
2H factors below S therefore avoid p. All denominator factorials are below p.
Every first-row coefficient is a unit. In the doubled row, reduction of its
falling upper product equals (2m+1)_(2k). It vanishes exactly when k>=m+1;
for k<=m it gives the displayed coefficient. All other denominator and
height factors remain units, proving (1).

Resultants commute with integral reduction and multiply across factors.
Since Res(Phi,t)=(-1)^H Phi(0) and
Res(Phi,R)=(-1)^(Hm) Res(R,Phi), their two signs combine with the
characteristic sign (-1)^H to give exactly (2). This argument keeps the
actual consumer, rather than just the vanishing pattern of coefficients.

## 2. Rank one supplies a second unbounded family

For m=1, R=t+2H. Put

    A_H=sum_(k=0)^H (H)_k (3/2)_(2k) (-2H)^(H-k)/(3k)!,

a rational number in lowest terms. If p>6H, p divides 2z+4H-3, and p does
not divide the numerator of A_H, then N_H(z) is nonzero. Its exact residue
is given by (2), using the single evaluation Phi_bar(-2H).

The rational exceptional value A_H is nonzero for EVERY H>=1. Indeed set

    a_k=(H)_k (3/2)_(2k)/[(3k)!(2H)^k].

Then a_0=1, a_1=1/16, and a_k>0 for every1<=k<=H. For1<=k<H,

    a_(k+1)/a_k = (H-k)/(2H)
       * (2k-3/2)(2k-1/2)/[(3k+1)(3k+2)(3k+3)] <1/2.

The finite alternating sum is therefore between15/16 and1, with the upper
inequality strict. Hence

    15/16 <= A_H/(-2H)^H <1.                         (3)

For each H there are primes greater than6H outside the finite prime divisors
of this nonzero numerator. For any such p, z=(p+3)/2-2H is positive, and all
positive translates z+lp satisfy the same certificate. Odd z occur infinitely
often, by choosing the parity of l. This gives an unbounded-height family
of exact boundary orders through the inherited jet theorem. No statement
about the real sign of the positive-integer norm follows from (3): its
half-integer endpoint evaluation serves only as a modular residual.

One explicit carried example is (h,r)=(11,4), hence H=7,z=9. The new prime
p=43 divides4h-1 and supplies N_7(9)=16 modulo43. The former endpoint4h+1=45
has no prime factor satisfying the previous unit criterion. Thus this is
an actual additional consumer, not only a different proof of an eligible
old-endpoint case. Its norm order at x=-4 is exactly3.

## 3. A minimal-height obstruction to automatic rank-one units

At H=2,

    Phi_(H,S=3/2)(t)=t^2+t/4+1/640,
    A_2=Phi(-4)=9601/640.

Take z=4798,p=9601. The prime exceeds6H and divides2z+4H-3, but both
reductions have the factor t+4. The modular norm is zero. These same rational
rows have norm20135 modulo28807, so they are coprime over Q. The first failed
implication is from a sparse doubled row to a unit residual. The strongest
survivor is (2), with its residual determinant retained. The missing
coordinate is precisely the exceptional numerator; the open question is
whether useful structure controls its prime factors at higher residual rank.
This hostile has the smallest possible height: A_1=-15/8 has no numerator
prime greater than6, whereas the prime9601 in A_2 exceeds12.

## 4. Exact scope and reproduction

The producer checks all172 triples obtained from H=1..18,
m=0..min(4,2H-1), and the first two primes in (6H,6H+44]. It reconstructs
both rational factorial rows, verifies (1), and compares the full degree-H
multiplication determinant with the degree-m determinant in (2). Rank-one
periodic lifts are checked from the actual rational rows. Eighty exact
half-endpoint evaluations check (3), alongside the exceptional-prime hostile
and its independent good-prime control. These finite checks support the
algebraic proof; they do not replace its unbounded quantifiers.

    python -B 04-computation/continuing7_20260906_norm_rank_drop.py
    python -B -O 04-computation/continuing7_20260906_norm_rank_drop.py

All gates remain active under optimization. The certificate is adjacent to
the source. General odd-z nonvanishing, the positive real norm sign at all
heights, and full Laurent noncancellation remain OPEN.
