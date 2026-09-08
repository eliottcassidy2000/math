---
id: THM-4458
title: "Sharp constant-gauge one-sided adverse-leak budget"
status: >
  PROVED ANALYTIC + FINITE-EXACT + INDEPENDENTLY AUDITED.
  Sharp real-profile variance inequality; weighted-median and dual compilers.
  Strict improvement over the old norm gate on an auxiliary density control.
  Canonical owner/word realization and LRC(14) remain OPEN.
source: euler-bridge-20260908
depends_on: []
scoped_consumers:
  - THM-3682-lrc-pure-role-to-lawful-target-leak-tariff
related:
  - THM-3674-sharp-successor-variance-drift-and-target-energy-tariff
script: 04-computation/lrc14_euler_adverse_leak_20260908.py
output: 05-knowledge/results/lrc14_euler_adverse_leak_20260908.out
script_sha256: 9082548c166e4f3235f4a3781b4c49737b830060258d5a90670b1b7d1bb42111
output_sha256: e0c4b570b64f8b89a6ec5041d1b8c98eddb1ac8ea149f8ee783f3fce8853aa52
hash_basis: raw LF bytes
audit: >
  Independent root review of the hinge/median identity, duality,
  sublinearity, all-budget estimate, sharpness, and inherited normalizations.
  7020 profile pairs checked by three exact rational paths; auxiliary circle
  intervals reconstructed independently. Explicit checks stay active under
  python -O, and outputs agree. No Lean claim and no PDE dependency.
---

# THM-4458 -- one-sided adverse-leak budget

**Status: PROVED ELEMENTARY + FINITE-EXACT within the statements below;
independently reviewed and promoted on 2026-09-08. LRC(14) remains
OPEN.** No Euler blowup theorem is used as a dependency. The external paper
supplies a research operation; all transferred inequalities are proved here.

## Inheritance, portfolio, and connection contract

- **Anchor:** sharpen the live distinction between a frozen role profile and
  an actual covariant target profile. Closest proved mechanism:
  [THM-3682](../../01-canon/theorems/THM-3682-lrc-pure-role-to-lawful-target-leak-tariff.md),
  including its corrected present-Q versus delayed-Q(Rx) typing.
- **Niche:** integral rank-one shears in the stabilizer of an LRC orbit.
  Closest relation carrier:
  [THM-4009](../../01-canon/theorems/THM-4009-euclidean-covering-transference-short-relation-compression.md).
- **Wildcard:** replace the full leak norm by a one-sided cost after removing
  its arbitrary constant gauge. The minimization becomes an exact weighted
  median, not a Fourier calculation or a numerical optimizer.
- **Canonical hostile:** THM-3682's translated danger density has a positive
  frozen-role energy but a completely constant physical marginal.
- **Corrected near miss:** THM-3682 cannot move a current-scale danger factor
  inside Q(Rx) when 13 divides R. Also retain MISTAKE-547: changing an
  observer changes its cancellation budget; this note fixes one R throughout.
- **Least-used relevant sidecar:** the signed product of the fixed centered
  role with the physical leak, retaining the source/target/word action.

The five-concept board after the probe is:

| Concept | Result of comparison |
|---|---|
| Full leak energy | Recovered as an upper bound on the new adverse cost |
| Positive/negative part | Favorable leakage is free, with an exact sharp bound |
| Constant gauge | A weighted median is the exact minimizing gauge |
| Physical phase restriction | Indicator-density control realizes the gain, but canonical owner realization is not supplied |
| Rank-one shear | Arbitrarily large ambient amplification can fix every physical phase |

The [Euler paper](https://cdn.openai.com/pdf/315b36cd-ec98-4023-8342-93345194ece1/euler.pdf),
Section 2.2, in its leading pressure-increment discussion, keeps
only the adverse positive pressure contribution in the coercivity budget.
The source operation is a one-sided estimate of a quadratic response. Here
the target is the projection of the actual profile C onto a fixed centered
role r. The map sends a leak L to the signed products -r_s(L_s-b), charges
their positive parts, and minimizes over the irrelevant constant b.
It preserves a certified lower bound for physical profile variance. It
forgets variance orthogonal to r, and supplies neither the actual LRC
owner/word intertwiner nor a bound on a canonical row's adverse budget.
The cheapest decisive test is the translated-density cancellation hostile;
the favorable-leak control tests whether this is stronger than the norm gate.

The source's moving phase and incompressible PDE are not transported by this
map. In particular, no claim about the external paper's complete proof is
needed or made. The phase-average/physical-phase distinction in its
independent-angle construction and subsequent graph restriction matches a **typing warning** here,
not an identification of the two systems.

## 1. Exact adverse budget and weighted-median compiler

Let p>=2 and R,L be arbitrary real p-vectors. Define

```
C = R+L,       r = R-mean(R),       E = mean(r^2)>0,
E(X) = mean((X-mean(X))^2),
B_R(L) = min_(b in R) mean_s (-r_s(L_s-b))_+ .                 (1)
```

The subscript R is load-bearing: changing the observer changes the budget.
No primality, nonnegativity, physical realization, or LRC hypothesis is
needed for the following lemma. The restriction to real profiles is explicit;
the weighted-median statement is not asserted for arbitrary complex profiles.

**Proposition 1 (PROVED).** The minimum in (1) is attained at every weighted
median of the values L_s with weights |r_s|. More explicitly,

```
B_R(L) = (min_b mean_s |r_s| |L_s-b| - mean_s r_s L_s)/2.       (2)
```

The function B_R is nonnegative, positively homogeneous, subadditive, and
unchanged by adding a constant to L. It satisfies

```
B_R(R)=0,       B_R(-R)=E,
B_R(L)^2 <= E E(L).                                          (3)
```

**Proof.** Use x_+=(|x|+x)/2 and sum r_s=0 to obtain (2). The first term
is weighted absolute deviation, so its minimum is at a weighted median.
For completeness, its left/right derivatives straddle zero exactly when
the weight strictly to either side of b is at most half the total weight.
Since r is nonzero and centered, the positive and negative weights both
have total sum |r_s|/2, and the objective is coercive as |b| tends to infinity.

Constant invariance and positive homogeneity follow by changing b.
The inequality (x+y)_+<=x_++y_+ and gauges b+c prove subadditivity.
For L=R, choosing b=mean(R) makes every summand zero. For L=-R, every
gauge costs at least mean(r_s R_s)=E, while b=-mean(R) attains E.
Finally choose b=mean(L), bound the positive part by |r_s||L_s-mean(L)|,
and apply Cauchy--Schwarz. This proves (3). QED.

An independent dual description, useful for exact auditing, is

```
B_R(L) = max { -mean_s z_s r_s L_s :
                 0<=z_s<=1, sum_s z_s r_s=0 }.                (4)
```

Weak duality follows term by term from the hinge maximum. At a minimizing
median choose z_s=1 where -r_s(L_s-b)>0 and z_s=0 where it is negative;
on tied coordinates choose intermediate z_s to make sum z_s r_s=0.
The median subgradient condition permits exactly this choice, proving
strong duality without importing a general minimax theorem. Vertices of
this box intersected with one hyperplane have at most one fractional
coordinate; coordinates with r_s=0 can be ignored.

## 2. Sharp improvement over the quadratic leak tariff

**Proposition 2 (PROVED).** For every pair R,L above,

```
                 E(C) >= (B_R(L)-E)^2 / E.                    (5)
```

Consequently B_R(L)!=E is a sufficient nonconstancy certificate. The
practically useful one-sided gate is B_R(L)<E: adverse leakage has not paid
the complete cancellation invoice. The exact equality B_R(L)=E is necessary
for constant C, but not sufficient.

**Proof.** Since L=C-R and -R=L-C, subadditivity and (3) give

```
B_R(L) <= E+B_R(C),
E <= B_R(L)+B_R(-C).
```

Both B_R(C) and B_R(-C) are at most sqrt(E E(C)), by (3).
Thus |B_R(L)-E|<=sqrt(E E(C)); squaring proves (5). QED.

This lower bound dominates the old THM-3682 reverse-triangle bound

```
E(C) >= (sqrt(E)-sqrt(E(L)))_+^2.                             (6)
```

If E(L)<E, then 0<=B_R(L)<=sqrt(E E(L))<E; substituting the upper
bound into (5) proves the comparison. If E(L)>=E, (6) is zero and the
comparison is automatic. Thus the improvement has a precise logical
direction, not merely favorable examples.

**Exact value versus bound.** Equation (5) uses the exact B_R(L). If only
an upper bound B_R(L)<=b_0 is certified, the universally valid lower
bound is (E-b_0)_+^2/E. One cannot square b_0-E when b_0>E: exact
cancellation has B_R(L)=E, also obeys the loose bound B_R(L)<=2E,
and yet E(C)=0. More generally an enclosure B_R(L) in [b_-,b_+]
gives dist(E,[b_-,b_+])^2/E. The lower endpoint matters for a
certificate on the B_R(L)>E side.

When the complete R,L profiles are already explicit, E(C) is directly
computable and stronger than a lower bound. The operational gain is the
possibility of certifying a one-sided adverse budget from partial or
structural information, without reconstructing the full target profile.

**Sharpness at every adverse budget.** Fix any nonconstant R and theta>=0.
Take L=-theta R plus any constant. Then B_R(L)=theta E,
E(C)=(1-theta)^2 E, and (5) is an equality. No larger bound depending
only on E and B_R(L) is possible.

**Failure boundary.** R=(-1,0,1), L=(1,1,-1) gives C=(0,1,0),
B_R(L)=E=2/3, but E(C)=2/9. The budget is a sharp sufficient detector;
it is not a complete classifier of nonconstant profiles. At E=0 all r
and all budgets vanish, and (5) is not defined or asserted.

## 3. Exact physical-density control that the old gate misses

Use p=13 and normalized Haar measure on R/Z. Put

```
D={y: ||y||<1/14},       A_s={y: ||y-s/13||<1/14},
w_0=1_D,                R_s=int 1_D 1_(A_s),
S={2,3,4,5,6},
w_s=1_(A_s) if s in S, otherwise 1_D,
C_s=int w_s 1_(A_s),     L_s=C_s-R_s.                        (7)
```

All densities are indicators, all have mass 1/7, and the distinguished
density at s=0 is unchanged. The exact overlaps are

```
R_0=1/7,       R_1=R_12=6/91,       R_s=0 otherwise;
L_s=1/7 for s in S,                 L_s=0 otherwise.
```

Direct rational interval arithmetic gives

```
E(R) = 2508/1399489,
E(L) = 40/8281 > E(R),
E(C) = 6018/1399489,
B_R(L) = 125/107653 = (1625/2508) E(R),
minimizing median b=0,
(E(R)-B_R(L))^2/E(R) = 779689/3509918412 > 0.                (8)
```

Thus the old norm gate supplies zero while (5) certifies nonconstancy.
This is an actual shifted-danger **auxiliary density realization**. It is
not asserted to arise from the canonical scalar-cover owner action.

The inherited hostile sets w_s=1_(A_s) for every s. Then C_s=1/7 is
constant, and B_R(L)=E(L)=E(R). The new inequality respects the exact
cancellation boundary. Favorable leakage L=2R instead has B_R(L)=0,
E(L)=4E(R), and E(C)=9E(R): a large leak norm need not be adverse.

## 4. Correctly typed LRC consequence and remaining obligation

Apply (5) to the **same present-Q auxiliary row** in THM-3682, with its
actual r and its actual leak formula. Write

```
Gamma = (B_R(L)-E(R))^2/E(R).
```

The already proved THM-3682/THM-3674 implications then give

```
Var(S_Q)>=Gamma/13,
D_Q>=Gamma/2028,
E_dt,Q>=Gamma/26364.                                         (9)
```

No new field, Fourier normalization, or physical evaluation is introduced.
In the useful one-sided regime, a certified budget bound B_R(L)<=b_0
and the inherited role floor E(R)>=e_0=121 rho^2/2028>b_0 imply
Gamma>=(e_0-b_0)^2/e_0. The function (e-b_0)^2/e is increasing for
e>=b_0, which justifies substituting the role floor.

**OPEN:** control this budget for an actual covering-row profile and supply
the delayed-Q(Rx) target/word intertwiner. Equation (9) does not replace
these tasks. The mechanism is a reduction of the debt coordinate, not an
exclusion of any of the 165 valuation profiles or a proof of LRC(14).

## 5. Integral shear hostile: amplification can fix the full orbit

**Proposition 3 (PROVED).** Let s be an integer speed vector and let
nonzero integer vectors a,m satisfy m.a=m.s=0. For any integer k put
S_k=I+k a m^T. Then det(S_k)=1, S_k^{-1}=I-k a m^T, and S_k s=s.
The induced torus automorphism fixes every point t s of the physical LRC
orbit. Consequently, for every set B in the torus,

```
{t: ts in S_k B} = {t: ts in B}.                             (10)
```

**Proof.** (a m^T)^2=0 gives the inverse; the rank-one determinant
identity gives det(S_k)=1+m.(k a)=1. The equality S_k s=s follows
from m.s=0. Apply S_k^{-1} to the membership ts in S_k B. QED.

For the AP13 hostile s=(1,...,13), take m=(1,-2,1,0,...,0) and a=e_4.
Then ||S_k-I||_op=|k| sqrt(6) is unbounded, while the complete physical
phase function remains min_i ||it|| and has maximum exactly 1/14.
For the safe box B=[1/14,13/14]^13, the exact transported inequalities
retain normals S_k^{-T}e_i. In particular the fourth normal is e_4-km;
its pairing with s is still 4. Erasing those normal/owner labels changes
the target set and invalidates the transfer.

This is the algebraic rank-one, orthogonal-pair analogue of the paper's
moving-frame syntax. It is **not** a map from an Euler solution to an LRC
row. Its useful conclusion is a precise stopping obstruction: relation-
stabilizing shear growth alone can create no new physical arrival. The
labelled inequality sidecar is essential.

## Reproduction, search, and scope ledger

Run:

```
python 04-computation/lrc14_euler_adverse_leak_20260908.py
```

Source:
[lrc14_euler_adverse_leak_20260908.py](../../04-computation/lrc14_euler_adverse_leak_20260908.py).
Frozen output:
[lrc14_euler_adverse_leak_20260908.out](../../05-knowledge/results/lrc14_euler_adverse_leak_20260908.out).
The exact finite universe is all R,L in {-1,0,1}^p for p=2,3,4,
excluding constant R, plus the displayed rational sharpness, density,
cancellation and favorable-leak controls. Hinge minimization, a weighted
median, and an independently enumerated dual polytope agree. Circle-interval
intersections independently reconstruct the p=13 density profile.
The deliberate loose-upper-bound hostile checks the distinction after (6).
All checks use explicit exceptions and remain active under Python -O. The
normal and optimized complete audit streams agree, ending in PASS.
Raw LF SHA-256 pins are:

```
script 9082548c166e4f3235f4a3781b4c49737b830060258d5a90670b1b7d1bb42111
output e0c4b570b64f8b89a6ec5041d1b8c98eddb1ac8ea149f8ee783f3fce8853aa52
```

The startup, relevant current frontier, guardrails, research protocol, core
papers, and targeted message peek were read. Searches for one-sided budgets,
adverse/anti-alignment, shear-orbit stabilizers, phase graphs, signed endpoint
budgets, and leak-energy formulas recovered THM-4002, THM-3682, THM-4449,
THM-4451, and the pertinent corrections. Their existing endpoint, positive
overlap, and norm results are credited rather than renamed. No literature
priority claim is made for this elementary weighted-median inequality.
META-PATTERNS used: search the statement before the method; correct the
object before sharpening the technique; type every analogy and implication.
No new method card is promoted from this single transfer.
