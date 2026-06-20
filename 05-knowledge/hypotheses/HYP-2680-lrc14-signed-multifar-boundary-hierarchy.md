---
id: HYP-2680
title: LRC(14) signed multi-far boundary hierarchy
status: OPEN; exact Newton/Stirling boundary formula, analytic resonance bound missing
source: codex-2026-06-20-S51
depends_on:
  - THM-548
  - HYP-2679
  - HYP-2678
  - HYP-2677
  - HYP-2676
  - HYP-2675
  - THM-547
  - THM-546
  - THM-531
  - HYP-2648
  - HYP-2639
related:
  - HYP-2638
  - HYP-2637
  - HYP-2636
  - HYP-2633
  - HYP-2614
  - HYP-2608
  - OPEN-Q-108
---

# HYP-2680 - Signed Multi-Far Boundary Hierarchy

## Claim Being Tested

The signed two-far bound of THM-548 and HYP-2679 is the second Newton
difference of one general boundary hierarchy.  If `B` is a bounded core and
`F` is a far set, define

```text
Delta_S(B) = sum_{T subset S} (-1)^{|S|-|T|} p0(B union T).
```

The exact finite expansion is

```text
p0(B union F) = sum_{S subset F, |S| <= 6} Delta_S(B),
```

because there are only six inner sectors to miss.

Let `p_t(B)` be the exact measure of points where `B*x` misses exactly `t` of
the six inner sectors.  In the fully decorrelated far-runner limit, the signed
`s`-far Newton term should be

```text
Phi_s(B) = 7^-s * sum_{t=1}^s (-1)^(s-t) * t! * S(s,t) * p_t(B),
```

where `S(s,t)` is a Stirling number of the second kind.  This recovers the
known cases

```text
Phi_1(B) = p1(B)/7
Phi_2(B) = (2*p2(B) - p1(B))/49
Phi_3(B) = (p1(B) - 6*p2(B) + 6*p3(B))/343.
```

The corresponding fully decorrelated boundary value for `r` far runners is
the Newton sum

```text
P_r(B) = p0(B) + sum_{s=1}^{min(6,r)} binom(r,s) * Phi_s(B),
```

equivalent to THM-548's formula

```text
sum_t p_t(B) * sum_{i=0}^t (-1)^i binom(t,i) (1 - i/7)^r.
```

## Proof Obligation

The remaining analytic content is not the formula above; it is the signed
resonance discrepancy

```text
R_S(B; f_1,...,f_s) = Delta_S(B) - Phi_s(B).
```

For `s=2`, THM-548 identifies frequencies `m*f_1+n*f_2` and a sharp
two-far constant `C_2=13/(4*7^3)`.  The next target is the three-far analogue:
bound

```text
Delta_{u,v,w}(B) - (p1 - 6*p2 + 6*p3)/343
```

by a signed Abel packet whose denominators are controlled by all small linear
forms

```text
m*u+n*v+l*w.
```

This is stricter than pairwise resonance control.  A triple can have no zero
pair relation and still carry a three-body relation such as `u - 2v + w = 0`.
That is the exact place where the user's "two things with opposite bounded
signs" becomes a relation-lattice statement: rank one gives a finite
scale-invariant ledger; higher rank should pay a signed dimension penalty.

## Expected Split

```text
No small relation among far speeds
  -> signed Abel / Koksma decay around Phi_s(B).

Rank-one far relation lattice
  -> Freiman/scale reduction to a finite boundary atlas.

Rank >= 2 relation lattice
  -> multi-form discrepancy bound with dimension penalty.
```

The rank split preserves the LRC predicate `p0(B union F) <= cap_k` through the
Newton expansion and destroys only the order in which far speeds are added.
That loss is acceptable because the mixed differences are symmetric after
summing over subsets.

## Challenged Assumption

Do not assume tournament vertices must be runners or far speeds.  In this
route the useful vertices can be:

- Newton orders `s=1..6`;
- relation-lattice ranks;
- proof obligations `Phi_s`, `R_s`, and direct cap margin;
- small linear forms among far speeds;
- state-word packets of the bounded core.

The quotient preserves signed boundary corrections.  It destroys individual
runner identity except through the retained relation forms.

## Next Computation

Reserve `codex-2026-06-20-S51` for an exact three-far boundary-hierarchy
script that will:

1. verify the Stirling coefficients through `s=6`;
2. scan exact three-far triples against `Phi_3(B)` on the known true-wide
   cores;
3. record pair and triple resonance distances separately;
4. run Tournament Analysis on proof-obligation vertices rather than runners.

No proof of LRC(14) is claimed here.
