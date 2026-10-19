# Sharp two-thirds scaling of transient amplification

19 October 2026. **PROVED, with independent audit PASS.** Extension of the global dyadic cofactor classification.

Let F(n)=n-p+2c-1 for positive even n=cp, c the exact power of two and p odd. Let M(n)=max_{j>=0} F^j(n), which is finite by the global classification. This is the cofactor map; it agrees with the aliquot map only while p is an odd prime.

**Theorem.**

    limsup_{n even, n->infinity} (M(n)-3n)/n^(2/3) = 3.

In particular the coefficient 3 and exponent 2/3 arise from an exact resonant family, rather than a fitted growth law. The formula does not assert this behavior for every start, and it does not assert any analogous general aliquot bound.

## 1. Bound the excursion before the smaller checkpoint

Use the cases of COFACTOR_MAP_GLOBAL_CLASSIFICATION.md. The descent algorithm selects a smaller checkpoint in at most six steps, except at fixed points and the explicitly handled finite exceptions.

In every nonresonant case the values up to the selected smaller checkpoint are at most 3n-2:

- If the first step descends, this is immediate.
- If q=p+1 is not a power of two, the first image is at most 3n-2 and the second is below n.
- If p=1 and n>=32, the three displayed images are bounded by 3n-2.
- If q=d>=4 is a power of two and d<c<d², the first image is n+2c-d<=5n/3, the second satisfies n_2-n=(c/d)(d-1)+d<=n, and the third is below n. Thus all are below 3n-2.
- If c>d², the same first and second bounds apply. For h=c/d²>=4 the third is below n; for h=2 the third is n+d<2n and the fourth is below n.

Fixed points have no excursion. The only larger envelopes occur at the two resonances:

    c=d:    n=d(d-1),       block maximum <=3d²-2;
    c=d²:   n=d²(d-1),      block maximum <=3d³-2.

For d>=8 these follow by reaching the power of two d² or d³ and using its exact three-step formula. For d=4 the displayed exceptional chains give the same maxima: M(12)=46 and the six-step block from 48 has maximum 190. The remaining exceptional starts have maximum at most 46, and specifically those below 12 have maximum at most 22.

## 2. Propagate through decreasing checkpoints

All subsequent checkpoints m are smaller than the original n. The nonresonant envelope is at most 3n-2. At a first resonance m=d(d-1)<=n,

    3d²-2=3m+3d-2 <=3n+3sqrt(n)+1,

because d<=sqrt(n)+1. At a second resonance m=d²(d-1)<=n, the inequality (d-1)^3<=m gives d<=n^(1/3)+1, hence

    3d³-2 = 3m+3d²-2
           <=3n+3n^(2/3)+6n^(1/3)+1.

The latter bound dominates the first-resonance bound for n>=1. Including the finite exceptional trajectories therefore yields the universal envelope

    M(n)<=max(46, 3n+3n^(2/3)+6n^(1/3)+1).

It follows that the stated limsup is at most 3.

## 3. Exact resonant lower bound

For every dyadic d>=8 take n_d=d²(d-1). The genuine consecutive orbit begins

    n_d -> d[d(d+1)-1] -> d³ -> 3d³-2.

Consequently

    (M(n_d)-3n_d)/n_d^(2/3)
      >= (3d²-2)/(d²(d-1))^(2/3) -> 3.

This proves the matching lower bound and the theorem.

There is also an exact maximum on the resonant family:

    M(d²(d-1))=3d³-2 for every dyadic d>=8.

After five steps the smaller checkpoint is m=3d³/4+4<d³. The universal envelope from Section 2 gives

    M(m) < (9/4)d³+3d²+6d+13 <3d³-2.

The constant alternative 46 is smaller as well. For the final strict inequality the gap is (3/4)d³-3d²-6d-15, positive at d=8 and positive thereafter: for d>=8 it is at least 3d²-6d-15>0. All earlier values of the five-step block are at most 3d³-2. Thus no later iterate exceeds the displayed peak.

The prime guard is not asserted along this family: its second checkpoint d³ has odd part 1. These resonances certify the global auxiliary-map envelope, not a native aliquot orbit of the same length.

## 4. Connection to the broader program

The global six-step descent proves eventual behavior; retaining intermediate values reveals a sharper scale that endpoint-only contraction would discard. The exponent 2/3 comes from the exact cubic scale n=d³-d² and the quadratic correction d². This offers a concrete template for analyzing transient height in another arithmetic map: classify resonant valuation transitions before promoting a return multiplier to a uniform height estimate. No transfer of this exponent to Collatz is justified by the present proof.
