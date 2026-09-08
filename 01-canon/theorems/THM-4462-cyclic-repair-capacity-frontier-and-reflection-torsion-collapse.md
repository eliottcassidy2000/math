---
id: THM-4462
title: "Cyclic repair capacity frontier and reflection torsion collapse"
status: >
  PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  A mean-zero integral cycle profile is Laplacian-repairable exactly when
  its first moment vanishes modulo the cycle order; reflection annihilates
  that obstruction at odd order. Bounded real divergence repair has an exact
  cyclic arc criterion. The odd-order symmetric danger-overlap role has a
  sharp two-regime minimum-variance law for every uniform flux capacity,
  with cancellation first at 125/1183 for the p=13 control. A physical
  owner/word flux construction, the delayed-word intertwiner, and LRC(14)
  remain OPEN.
source: cross-concepts-20260908
depends_on: []
related:
  - THM-3990-componentwise-harmonic-obstruction-and-repair-quotient
  - THM-3682-lrc-pure-role-to-lawful-target-leak-tariff
  - THM-4458-lrc-one-sided-adverse-leak-budget
  - THM-3392-bipartite-sign-lift-and-synchronization-loss
  - THM-2195-transitive-quotients-exactly-control-universal-substitution-products
  - THM-2348-prime-type-rectangularity-and-target-token-conditioning
script: 04-computation/cross_concepts_tournament_20260908.py
output: 05-knowledge/results/cross_concepts_tournament_20260908.out
script_sha256: 4423e2065406547eb731d87e5942e15bd1e33cdf9f8717b934165ab20afce747
output_sha256: f764419858d8dfa25cac814dd3b3bce8a0f1141c601b55531973971cad3d6662
semantic_sha256: 683412cf9fb449d5ab35090b13ce3c82c10d1a65e137de7b316e543aaf947d48
hash_basis: raw LF bytes
audit: >
  Independent donor-memory review accepted the overlap scaling, both phase
  transitions, explicit feasible fluxes, first-variation signs, odd-order
  torsion collapse, p=5 boundary, and p=13 integral potential without a
  mathematical correction. The 12965 explicit gates remain active under -O;
  normal, optimized, and stored LF transcripts agree. Independent paths
  include reduced-Laplacian elimination, cyclic-arc enumeration, exact
  box optimality certificates, and exhaustive small integer flux searches.
---

# THM-4462 -- cyclic repair capacity and reflection torsion collapse

**PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
The application is a relaxation of profile transport: an owner/word
construction supplying the bounded local flux is **OPEN**. LRC(14) remains
**OPEN**. No global novelty claim is made.

The useful change of question is from whether a centered role can be repaired
to the exact capacity needed to repair it. An odd-cycle reflection symmetry
kills the integral repair obstruction entirely. Nevertheless the familiar
13-shift danger profile has a sharp nonzero cancellation threshold, and its
best residual variance can be computed for every allowed capacity.

## Inheritance, portfolio, and search boundary

The closest proved mechanism is [THM-3990 / componentwise-harmonic-obstruction-and-repair-quotient](THM-3990-componentwise-harmonic-obstruction-and-repair-quotient.md):
real Laplacian repair leaves component means, while integral repair can retain
torsion. The least-used relevant sidecar here is the **size of the repair**
after that quotient vanishes. The canonical hostile is [THM-3682 / pure-role-to-lawful-target-leak-tariff](THM-3682-lrc-pure-role-to-lawful-target-leak-tariff.md),
whose translated danger densities make the physical profile constant. The
corrected near miss is its present-Q versus delayed-Q(Rx) distinction; no
cyclic graph construction here supplies that missing intertwiner.

Other recovered mechanisms were [THM-3392 / bipartite-sign-lift-and-synchronization-loss](THM-3392-bipartite-sign-lift-and-synchronization-loss.md),
whose ternary support partition records exactly what a scalar relaxation
forgets; [THM-2195 / transitive-quotients-exactly-control-universal-substitution-products](THM-2195-transitive-quotients-exactly-control-universal-substitution-products.md),
whose quotient operation can make expensive marked changes cheap; and
[THM-2348 / prime-type-rectangularity-and-target-token-conditioning](THM-2348-prime-type-rectangularity-and-target-token-conditioning.md),
which requires factorization under every perturbation. Their lesson is to
retain the allowable operation, not to impose a tournament on phase labels.
There is no intrinsic tournament orientation in this problem.

The anchor is LRC profile transport, the niche is cyclic integral repair, and
the wildcard is a sharp constrained diffusion problem on the same labels.
The concept board after the calculation is:

| Concept | Consequence of the calculation |
|---|---|
| Integral repair quotient | Reflection forces its class to zero at odd order |
| Source/target cancellation | Can occur integrally, with an explicit potential |
| Edge flux capacity | Exact obstruction is a finite set of cyclic arc cuts |
| Profile energy | Canonical profile has a sharp two-regime variance law |
| Physical owner/word action | A bounded-flux intertwiner still has to be supplied |

Targeted canon searches for cycle/cyclic Laplacian, critical-group reflection,
bounded flux, cumulative leak and prefix oscillation found the general repair
mechanism above, but no copy of the formulas below. THM-2953's Reynolds
Laplacian is an observation-kernel problem; THM-4110's reciprocal phase sheets
are a different weighted lattice. THM-3273/3277's internal critical quotient
does not identify these physical phase labels. Matching correction searches
retained the warnings in the opening Euler-bridge entry of MISTAKES.md and
the THM-3682 scale repair. This is a scoped search, not a literature claim.

The META-PATTERNS cards used were “Compute the repair quotient before testing
the residual defect” and “Type every analogy and every implication.” No new
card is proposed on the evidence of this lane alone.

## 1. Exact integral cycle quotient and its reflection boundary

Let n>=3, label the cycle by s in Z/n, and set

```
(Delta x)_s = 2x_s-x_(s-1)-x_(s+1),
(div f)_s = f_s-f_(s-1).
```

**Proposition 1 (PROVED, iff).** For a in Z^n with sum a_s=0, the equation
Delta x=-a has an integral solution exactly when

```
q(a) = sum_(s=0)^(n-1) s*a_s = 0 mod n.                 (1)
```

If a_s=a_(-s), then 2q(a)=0 mod n. In particular **every such symmetric
profile has an integral repair when n is odd**; primality is not required.

**Proof.** Put A_s=sum_(j=0)^s a_j, so A_(n-1)=0. Every solution of
div f=-a is f_s=h-A_s for one scalar h. Such a flux is of the form
f_s=x_s-x_(s+1) on a cyclic potential exactly when sum f_s=0, forcing
h=(sum A_s)/n. Since

```
sum_s A_s = sum_j (n-j)a_j = -sum_j j*a_j,
```

h is integral exactly under (1). If h is integral, integrate f starting
with x_0=0; sum f=0 makes x periodic and integral. Necessity follows from
the same equations. Reflection s->-s negates q while fixing a, hence
q=-q. At odd n, multiplication by 2 is invertible modulo n. QED.

The missing-symmetry hostile a=(1,-1,0) at n=3 has q=2. The even-order,
reflection-symmetric hostile a=(1,0,-1,0) at n=4 has q=2. Both have real
repairs but no integral repair. The strongest survivor without oddness is
the two-torsion restriction 2q=0, not universal vanishing.

This is the familiar cycle critical-group computation proved directly;
the contribution is its consequence for the symmetric role and its repair
cost below. It is not an obstruction to general LRC rows or to real fluxes.

## 2. Bounded flux has an exact arc-cut criterion

Let a in R^n have zero sum and give edge s (joining s to s+1) a finite
nonnegative capacity k_s. The relaxed corrected profile is c=a+div f with
|f_s|<=k_s. Its mean remains zero.

**Proposition 2 (PROVED, iff).** Complete cancellation c=0 is feasible iff

```
intersection_s [A_s-k_s, A_s+k_s] is nonempty,            (2)
```

equivalently, iff every proper nonempty cyclic arc I has

```
|sum_(s in I) a_s| <= k_(left boundary)+k_(right boundary). (3)
```

For uniform capacity k, the exact cancellation threshold is

```
k_* = (max_s A_s-min_s A_s)/2.                           (4)
```

**Proof.** The cancellation flux f_s=h-A_s respects capacities exactly
when h belongs to (2). A finite family of closed real intervals intersects
iff every pair intersects. Pairwise intersection says
|A_i-A_j|<=k_i+k_j, and A_i-A_j is the sum over one cyclic arc, up to sign.
This gives (3); (4) is the smallest radius covering the prefix values. QED.

In particular every feasible correction satisfies the energy certificate

```
mean(c^2) >=
 (|sum_I a_s|-k_left-k_right)_+^2 / (|I|*(n-|I|)).       (5)
```

To prove it, sum div f over I and use the two boundary capacities. Apply
Cauchy--Schwarz to the zero-mean vector c and the centered indicator of I,
whose squared norm is |I|*(n-|I|)/n. The positive part is essential when
only an upper capacity bound is known. This is not THM-4458's invalid
substitution of an upper bound into an unrestricted squared expression.

General network cut duality is not being advertised as new. On a cycle,
the elementary prefix formula retains the exact circulation coordinate,
and the following canonical family permits the entire optimum to be solved.

## 3. Sharp variance law for the symmetric danger-overlap family

Let p>=5 be odd and use the circle with normalized Haar measure. For

```
D={y: ||y||<1/(p+1)},
R_s=integral 1_D(y)*1_D(y-s/p) dy,
```

only s=0,+1,-1 give nonzero overlap:

```
R_0=2/(p+1),       R_(+1)=R_(-1)=(p-1)/(p*(p+1)).       (6)
```

Put d=p^2*(p+1)/2 and a=d*(R-mean R). Then

```
a_0=(p-1)^2=:u,
a_(+1)=a_(-1)=(p^2-5p+2)/2=:v,
a_s=-(2p-1)=:-b at the remaining m=p-3 positions.
A=u+2v=m*b,
k_0=p*(p+1)/4=(u-v)/2,
k_*=A/2.                                                (7)
```

For real k>=0 define the complete relaxed minimum

```
V_p(k)=min_(|f_s|<=k) mean((a+div f)^2)/d^2.             (8)
```

Thus actual normalized leakage is div f/d, with edge capacity k/d.

**Proposition 3 (PROVED, sharp for every k).**

```
V_p(k) = [(u-2k)^2+2v^2+m*(-b+2k/m)^2]/(p*d^2)
                                             for 0<=k<=k_0;
         (A-2k)^2/(3m*d^2)                    for k_0<=k<=k_*;
         0                                   for k>=k_* . (9)
```

The minimizing residual profile is unique and equals

```
c_0=u-2k, c_(+1)=c_(-1)=v, c_other=-b+2k/m    (k<=k_0);
c_0=c_(+1)=c_(-1)=(A-2k)/3,
c_other=-(A-2k)/m                            (k_0<=k<=k_*);
c=0                                         (k>=k_*).   (10)
```

**Explicit proof of feasibility and optimality.** For 0<=k<=k_* let

```
f_1=-k,
f_s=-k+2k*(s-1)/m          for 2<=s<=p-2,
f_0=-k, f_(p-1)=k          for k<=k_0,
f_0=-(2k_0+k)/3,
f_(p-1)=(2k_0+k)/3         for k>=k_0.                   (11)
```

The middle progression reaches +k at s=p-2. The second choice of f_0
has magnitude at most k precisely when k>=k_0. Direct differences give
(10), including at both joining endpoints. For k>=k_*, retain the
cancellation flux at k_*.

For an arbitrary competing feasible flux g, write h=g-f. Summation by
parts gives

```
sum (a+div g)^2 - sum c^2
 = 2 sum_s (c_s-c_(s+1))*(g_s-f_s) + sum (div h)^2.      (12)
```

In (11), whenever c_s>c_(s+1), the flux is f_s=-k; when the inequality
reverses, f_s=+k. At equal residual values the term vanishes. Every term
in the first sum of (12) is therefore nonnegative. This proves global
optimality directly, with no optimizer or general duality theorem. Strict
convexity of the squared norm in the residual gives its uniqueness.
Squaring (10) gives (9). The hypotheses p>=5 ensure v>0 and k_0<k_*.
QED.

In the second regime (9) is exactly the arc certificate (5) for the m
consecutive negative sites, so that particular cut certificate is sharp.
In the first regime there is also an unsaturated interior equalization
obligation; one cut alone need not recover the full optimum.

## 4. The exact p=13 comparison and what torsion forgets

At p=13, d=1183 and

```
a=(144,53,-25,-25,-25,-25,-25,-25,-25,-25,-25,-25,53).
k_0=91/2,       k_*=125.
```

Thus (9) reads

```
[(144-2k)^2+2*53^2+10*(-25+k/5)^2]/(13*1183^2),  k<=91/2;
(250-2k)^2/(30*1183^2),                         91/2<=k<=125;
0,                                             k>=125. (13)
```

At k=0 this is 2508/1399489, the energy in
[THM-4458 / one-sided-adverse-leak-budget](THM-4458-lrc-one-sided-adverse-leak-budget.md).
At k=100 it is 250/4198467. Cancellation first occurs at actual flux
capacity 125/1183.

The integral repair quotient is already zero. An explicit integral
potential solving Delta x=-a is

```
x=(0,72,197,297,372,422,447,447,422,372,297,197,72),
f_s=x_s-x_(s+1)
 =(-72,-125,-100,-75,-50,-25,0,25,50,75,100,125,72).      (14)
```

All potentials differ by a constant, so their range is always 447. Their
minimum real sup norm modulo constants is 447/2, and the minimum integral
sup norm is 224. Uniform flux capacity 125 is attained by (14) and required
by the ten-site negative arc. The corresponding general potential range is

```
range(x)=(p-1)*(2p^2-3p-1)/8.                            (15)
```

Indeed by symmetry the potential increases to the two middle sites;
its height is u/2+b*(p-3)*(p-1)/8, which simplifies to (15).

The THM-3682 hostile, w_s(y)=1_D(y-s/p), has constant C_s=2/(p+1).
Its centered leak is exactly -a/d, hence its minimum-flux representation
is (14)/d. The capacity threshold respects a real density cancellation
example; it does not forbid a known feasible target by mistyping it.

## 5. Connection contract and next decisive obligation

**Source -> target:** THM-3990's integral graph repair problem -> the
cyclic phase profile in THM-3682/4458.

**Map:** clear the explicitly fixed denominator d, center the role, and
write centered leakage as div f/d on the **fixed cyclic ordering** of the
phase labels. The potential subclass is f_s=x_s-x_(s+1).

**Preserved predicate:** cancellation is exactly div f=-a. With a certified
capacity bound, (9) gives the sharp residual variance for this role.

**Destroyed information:** owner identity, physical word action, nonnegative
density realization, and the current-versus-delayed scale distinction are
not encoded by a free flux. The unrestricted integral quotient also loses
the size of its repairing potential. No tournament is asserted.

**Needed sidecar:** a canonical construction of f from the actual moving
densities/owner action, together with an a priori edge-capacity estimate.
Every zero-mean leak has some cyclic flux, so existence alone adds nothing.
The cycle ordering and edge norm are not invariant under arbitrary
permutations; a proposed action must preserve that order or transport the
metric/capacities too. An integral potential is an extra hypothesis and is
not automatically supplied by the density model.

**Cheapest decisive test:** construct that flux on the present-Q auxiliary
table and test whether its gauge-minimized maximum magnitude is strictly
below 125/1183 in the canonical overlap control. Failure identifies which
arc transports the missing mass. A successful bound on actual canonical
rows, followed by the still-missing delayed-word intertwiner, would be a
real next result. Neither is claimed here.

**Stopping reason for the torsion lane:** the reflection lemma annihilates
the entire candidate obstruction, not merely a finite sample. The stronger
surviving question is quantitative admissibility of the repair. The exact
bounded-flux law supplies the replacement object and its equality boundary.

## 6. Reproduction and scope

Run:

```
python3 04-computation/cross_concepts_tournament_20260908.py
python3 -O 04-computation/cross_concepts_tournament_20260908.py
```

The transcript is
[cross_concepts_tournament_20260908.out](../../05-knowledge/results/cross_concepts_tournament_20260908.out).
The dependency-free verifier uses explicit exceptions, so -O retains every
gate. Its declared finite universes include 1089 mean-zero integer vectors,
2430 heterogeneous capacity instances, 140 exact optimality certificates
across odd p=5..31, and 3369 enumerated small integer-flux vectors.
The integral image is independently checked by reduced-Laplacian Gaussian
elimination through n=5. The capacity criterion is checked by separately
enumerating cyclic arcs. The profile optimum is checked both by exact box
first variation and by the small direct flux search. These finite audits
support, but do not replace, the symbolic proofs above.

The donor-memory agent independently audited the overlap normalization,
both transition values, explicit feasible fluxes, first-variation signs,
p=5 boundary, and p=13 integral potential. Verdict: **PASS with no
mathematical correction**. Normal and optimized outputs agree with the
stored LF transcript and semantic hash
`683412cf9fb449d5ab35090b13ce3c82c10d1a65e137de7b316e543aaf947d48`.
