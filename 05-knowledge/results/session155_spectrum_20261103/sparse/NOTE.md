# Congruence coverage with zero-density exact output

Let F be the union of the two literal strong-tournament companion families from the S154 mixed-prefix construction, using all finite q2/triangle prefix words with **even** triangle parity. These have column seed (Y,X)=(25,34), actions

    q: (Y,X) -> (11Y+4X,18Y+7X),
    b: (Y,X) -> (5Y+2X,8Y+3X),

and counts H=3(Y+6)+2 or H=3(Y+6)+4. The source construction and the existence of the actual companion tournaments are inherited from `golden_inverse_defects_20261103/TRANSPORT.md` at f6d5304849e50b29df6e368239e92fffa92e57e3 and the S154 mixed-prefix packet at fb8bc7775a0678ecf6ef4169faff452220573572. The assertions below concern this named family, not the full tournament spectrum.

## Sparse-bank theorem

Let δ be the unique real number with

    17^(−δ)+8^(−δ)=1,       δ≈0.284524642884874.

For every σ>δ and X>0, the number of **distinct** counts H∈F with H≤X is at most

    2 (8X/579)^σ / (1−17^(−σ)−8^(−σ)).

Consequently F has zero natural density even among positive odd integers. In particular, the full modular coverage of the named family at powers of 7,13,19,31 does not by itself provide an exact representation of every odd integer or every 7p.

**Proof.** Every state has positive Y,X. After any first action, the new state satisfies 2X≥3Y: under q, 2(18Y+7X)−3(11Y+4X)=3Y+2X>0; under b, 2(8Y+3X)−3(5Y+2X)=Y>0. Hence every subsequent q action multiplies Y by at least 17, and every subsequent b action multiplies Y by at least 8.

Let r,s be the total numbers of q,b letters. The empty word has Y=25≥193/8. The first q step gives Y=411≥(193/8)17, and the first b step gives Y=193=(193/8)8. The preceding paragraph and induction give the exact bound

    Y ≥ (193/8) 17^r 8^s

for every word, regardless of parity. For even parity the two companion counts are 3Y+20 and 3Y+22, so both are at least (579/8)17^r8^s.

The number of distinct counts is no larger than the number of word/offset presentations. For σ>δ, apply 1_{H≤X}≤(X/H)^σ, discard the even-parity restriction, and sum over all finite words. The number-weighted series is

    sum_w 17^(−σ r(w))8^(−σ s(w))
      = sum_{n≥0}(17^(−σ)+8^(−σ))^n
      = 1/(1−17^(−σ)−8^(−σ)).

There are two offsets, giving the displayed bound. The defining left-hand side of the δ equation is strictly decreasing from 2 at δ=0 to 0, so δ exists uniquely and lies below 1; hence one may choose δ<σ<1. The conclusion is independent of possible collisions between distinct word presentations. The exact inequalities imply the theorem; numerical checking is not needed for the proof.

## Exact 7p consequence in each fixed reduced class

Fix nonnegative exponents b,c,d and M=7*13^b*19^c*31^d. Let j be any residue with gcd(j,M)=1. The count of primes p≤x, p≡j mod M, for which 7p belongs to F is Oσ(x^σ) for every σ>δ, by the theorem with X=7x. The prime number theorem in arithmetic progressions says the total number of primes p≤x in that class is asymptotic to x/(φ(M) log x). Choosing δ<σ<1 therefore shows that **a proportion tending to 1** of primes in every such fixed reduced class have 7p outside F.

The class must be reduced for the cited prime asymptotic; no assertion about an unreduced class is implied. A prime p in a reduced class modulo M is automatically different from 7, so 7p has exact 7-adic valuation 1. The modular result nevertheless supplies infinitely many family counts with exact valuation 1 and H/7 in that class; this does not say any such quotient is prime. The result does not rule out 7p in other tournament families, even if 7p is missing here.

## Check

    python3 check.py

The script enumerates all 65,535 action words through length 15 using exact integers and verifies the cone and growth bounds, and numerically locates δ at 60-digit Decimal precision. Enumeration is a diagnostic and does not replace the proof. The exact prime-class consequence uses the classical prime number theorem for arithmetic progressions, not this computation.
