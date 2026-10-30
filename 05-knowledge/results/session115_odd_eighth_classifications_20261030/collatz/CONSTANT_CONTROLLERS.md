# Constant-letter controllers: all stalls, all entries, and a typed eighth

**PROVED special-class classification.** Every controller consisting of a repeated fixed positive valuation a has an explicit maximal all-layer mean over all entries. Its sharp original-clock bound is5/8, attained only by valuation1 with an all-one entry. Arbitrary controllers and Collatz convergence remain open.

This uses the exact marked-periodic theorem from [Session114](../../session114_continuation_symmetry_20261030/collatz/PERIODIC_CONTENT_RATES.md). The new family-level conclusions are the all-stall bound and equality classification, plus a precise parameter-space interpretation of3/8. Time mean, clock-normalized rate, and parameter-space mean are different objects.

## 1. Full maximal-rate formula

For the one-letter controller (a), F_a(X)=(3X+1)/2^a. Its multiplier has residue order

 t_a=2 if a=1 mod4;  t_a=4 if a is even;  t_a=1 if a=3 mod4.

The affine lifting depth is

 r_a=v5(9−4^a) if a=1 mod4;
 r_a=v5(81−16^a) if a is even;
 r_a=0 if a=3 mod4.

For every finite positive entry word u,

 **R(u,(a)) <= R(empty,(a))=(r_a+1/4)/t_a.**           (1)

Equality holds exactly when u is a power of the one-letter word a, including the empty word. For a controller a^d the period mean is d times the displayed rate and its clock is da, so the normalized conclusion is unchanged.

For a nonmaximal entry, the exact value is the Session114 geometric cutoff f(t_a,r_a,D_a(u)). Its depth is obtained without subgroup enumeration from the following entry data. If (P_u,Q_u,B_u) are its affine coefficients, set

 E_a(u)=(3−2^a)B_u−(P_u−Q_u),
 (R,S,C)=(3P_u,2^a Q_u,P_u+2^a B_u),
 H=F_(a u).

The commutator valuation is v5(E_a(u)), since the commutator numerator for F_a and H is2^a E_a(u). Apply the explicit depth rule in the cited theorem to this H. This is a finite exact classification of every entry's rate, not a finite search over entries.

**Proof of (1) and equality.** There is one phase. Empty entry puts that phase in the cyclic subgroup at every depth, giving cutoff infinity; the geometric sum is strictly increasing in cutoff. Equality therefore requires the entry's phase to commute with F_a as a rational affine map. Faithfulness of positive valuation words implies that the word a u commutes with a. The elementary word centralizer is exactly the powers of a, so u must be such a power. Conversely powers plainly attain the infinite cutoff. The residue order and r_a formulas follow directly from3/2^a modulo5. In the residue1 case the translation1/2^a is a unit, so the affine depth is0 even when the multiplier itself has a deep congruence.

## 2. The sharp all-stall clock bound

At a=1, t=2,r=1 and the maximal rate is5/8. At a=3 mod4 the rate per clock is1/(4a)<1/2. For all remaining a>=2, t is2 or4. Only2^{at} dominates3^t, and

 5^r <= 2^{at}−3^t <2^{at}<5^(at/2).

Thus r<at/2. Since at is even and r integral, r<=at/2−1. Therefore

 **R(u,(a))/a <= (r+1/4)/(ta)
                 <=1/2−3/(4ta)<1/2**                 (2)

for these a. Combining cases proves the sharp5/8 bound and its unique equality family. No bound on the depth of a stall beyond the elementary integer-size inequality was assumed.

Examples (maximal entries): a=1 gives5/8 per clock; a=2 gives9/32; a=3 gives1/12; a=5 gives1/8. At a=17, the mean per block is9/8 but the clock rate is only9/136. Erasing the clock would confuse distinct normalization problems.

## 3. Exact stall addresses and a genuine Haar3/8

In each fixed residue class a modulo4 other than3, the map from the five-adic exponent coordinate to the relevant principal unit is an analytic bijection with unit derivative after dividing by5. Equivalently, the elementary LTE identity shows that2^4 generates1+5 Z_5. Hence for every k>=1, exactly one residue class modulo5^(k−1) within that fixed modulo4 class has r_a>=k. The conditional Haar law is

 Pr(r_a>=k)=5^(1−k),  Pr(r_a=k)=4/5^k.

This statement can be read purely at finite precision: uniform residues modulo4·5^E give exactly those proportions for k<=E+1. It is not an unproved statement about averaging an unbounded function over an initial integer interval.

On the profinite exponent space Z/4 × Z_5, extend the numerical expression for the maximal **block mean** in (1) by these valuation formulas. This extended function is integrable; each possible infinite-valuation parameter is a singleton of measure zero. A noninteger profinite parameter is not being called a legal Collatz valuation or an actual controller. Conditional on a being even, t=4 and

 E[r_a]=sum_(k>=1)5^(1−k)=5/4,
 **E[R(empty,(a)) | a even]=(5/4+1/4)/4=3/8.**        (3)

Conditional on a=1 mod4 the corresponding mean is3/4; on a=3 mod4 it is1/4. The total parameter-Haar mean is7/16. These follow from monotone sums of the finite residue proportions. They are neither orbit means for one Collatz source nor ordinary averages of clock-normalized rates in a.

The distinguished all-one controller has a different exact decomposition:

 5/8 = sum_(e>=1) 1/(2·5^(e−1)) =1/2+1/8.

Its high-layer tail is1/8. Formula(3)'s3/8 is a parameter-Haar mean over a different family of return clocks. These are genuine arithmetic constructions of eighths, but their measures and normalization differ. The formal complements3/8=1−5/8 and7/8=1−1/8 do not automatically describe return events or transfer a theorem to another field.

## Reproduction

`check_constant.py` checks the residue classification, all-entry maximum/equality on finite universes, arbitrary-depth formula controls, exact finite residue proportions, the elementary uniform bound, and the distinct eighth examples. Written proofs supply all-parameter coverage. Normal and optimized checks use explicit exceptions.
