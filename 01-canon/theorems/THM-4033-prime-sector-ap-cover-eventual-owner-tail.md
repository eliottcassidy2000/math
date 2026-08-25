---
id: THM-4033
title: "Prime-sector AP-cover eventual finite-owner tail"
status: >
  PROVED + VERIFIED-EXACT. For every fixed prime P, the AP sector-cover
  deficit has an eventual finite-owner decomposition over exactly the
  reduced rationals of denominator below P, is eventually phase-rational
  with period dividing lcm(1,...,P-1), and has an explicit positive
  1/m constant and two-sided rational bounds. The onset is existential and
  pointwise in P. The observed odd-prime sharp onset (P^2+3)/4 is
  FINITE-EXACT only; the P=2 edge is m_0=2.
source: root + prime_sector_theorem / generated sequence task, 2026-08-24
audit: >
  PASS. Exact owner coefficients and the totient constant agree through
  prime 31; the denominator-P two-sided robustness trace passes; and a
  separate direct rational-wall engine agrees on the sampled tails for
  P=2,3,5,7,11,13,17,19,23 while retaining every tested pre-onset
  mismatch. Exact guard containment and owner nonoverlap pass through prime
  31, and modulus four is a composite hostile. Normal and optimized outputs
  are byte-identical.
depends_on: []
related:
  - THM-4029-lrc14-ap-cover-twelve-owner-rational-tail
  - THM-4031-lrc14-endpoint-owner-rational-tail
script: 04-computation/prime_sector_ap_cover_eventual_owner_tail_thm4033.py
output: 05-knowledge/results/prime_sector_ap_cover_eventual_owner_tail_thm4033.out
script_sha256: 77c3de5f38e553839dbabb44cd71bf405fba3dd0389a63a8e633e0b4397fa33b
output_sha256: 39c330e077a0a66c74eb6358679459b372035888484fa9e5627f8d3c5fa88025
hash_basis: raw LF bytes
---

# THM-4033 -- prime-sector AP-cover eventual finite-owner tail

**PROVED + VERIFIED-EXACT.** The persistent-owner classification, compactness
argument, finite-phase selector lemma, and closed constant below are proved.
The sharper numerical onset in Section 8 remains computational evidence only.

## 1. Object and statement

Fix a prime `P`.  For `m>=1`, let

```text
a_P(m)=meas{x in R/Z:
  {floor(P e x) mod P:0<=e<m}=Z/PZ},
D_P(m)=1-a_P(m),
n=m-1.
```

The persistent owners are

```text
B_P={p/q in R/Z:gcd(p,q)=1, 1<=q<P}.
```

For `p/q in B_P`, use the lift `theta_0=Pp/q` on the circle of length `P`.
For `0<=s<q`, put

```text
A_s=Pps/q,  b_s=floor(A_s) mod P,  f_s=frac(A_s),
E_s(n)=n-((n-s) mod q).
```

For a missing sector `r` (that is, `r` is not among the `b_s`), let

```text
d^+_(s,r) in {1,...,P-1} represent r-b_s mod P,
d^-_(s,r) in {1,...,P-1} represent b_s-r mod P,
c^+_(s,r)=d^+_(s,r)-f_s,
c^-_(s,r)=d^-_(s,r)-1+f_s,

rho^+_(p/q)(n)=max_r min_s c^+_(s,r)/E_s(n),
rho^-_(p/q)(n)=max_r min_s c^-_(s,r)/E_s(n).       (1)
```

All denominators in (1) are positive once `n>=P-1`.

### Theorem

For each fixed prime `P`, there is a finite geometric onset `M_P^geom` such
that, for every `m>=M_P^geom`, the noncover set, modulo finitely many wall
points, is the disjoint union of the owner pieces

```text
[theta_0-rho^-_(p/q)(n), theta_0+rho^+_(p/q)(n))
```

on the length-`P` theta circle.  In particular,

```text
D_P(m)=(1/P) sum_(p/q in B_P)
             (rho^-_(p/q)(m-1)+rho^+_(p/q)(m-1)).   (2)
```

Let

```text
L_P=lcm(1,2,...,P-1).
```

For each phase `r mod L_P`, there is a (phase-dependent) finite onset
`M_(P,r)` such that on that phase every signed owner radius in (1) is one
rational term

```text
C_(p/q,sign,r)/(n-c_(p/q,sign,r)),
C_(...)>=0, 0<=c_(...)<=q-1.                       (3)
```

Thus `D_P` is eventually `L_P`-phase rational.  The assertion is pointwise in
`P`: no uniform numerical onset in `P` is claimed.

Finally,

```text
D_P(m)=C_P/m+O_P(m^-2),

C_P=(1/P) sum_(q=1)^(P-1) phi(q)[2(P-q)-1]/q.       (4)
```

After the phase onsets, the nonnegative representation also gives

```text
C_P/(m-1) <= D_P(m) <= C_P/(m-P+1).                 (5)
```

## 2. Persistent-owner classification, including `q=P`

Write the orbit in the original `x` circle.  An irrational orbit is dense,
so it enters the interior of every one of the `P` sectors.  A rational of
reduced denominator `q>P` has the grid of `q` equally spaced points.  Every
open sector has length `1/P>1/q`, hence contains a grid point strictly in its
interior.  In either case, select one strict interior hit per sector.  Those
finitely many hits persist under a small perturbation, so coverage is locally
robust at a finite horizon.

The equality case `q=P` needs a separate argument because every grid point is
on a sector wall.  Let `x=p/P+epsilon`, or equivalently
`theta=p+delta` with `delta=P epsilon`.  For `0<delta<1/P`, times
`e=0,...,P-1` give residues

```text
ep mod P,
```

which cover all sectors.  For `-1/P<delta<0`, times `e=1,...,P-1` give

```text
ep-1 mod P,
```

which give every residue except `P-1`, while time `e=P` gives `P-1`.
At `delta=0`, times through `P-1` cover directly.  Therefore one common
horizon, `e<=P`, covers an entire two-sided neighbourhood.  The negative-side
time `P` is essential; omitting it is the boundary bug in the informal
argument.

If `q<P`, the orbit has only `q` points and cannot occupy `P` sectors.  These
and only these rationals are persistent noncovers, proving the classification.

## 3. Local track lemma and endpoint conventions

Because `P` is prime and `q<P`, `gcd(P,q)=1`.  Hence `f_s>0` for every
`s>0`, and the `q` values `b_s` are distinct.  Choose a small owner arc on
which, for both drift signs,

```text
q|delta|<1,
s delta < 1-f_s       (delta>0, s>0),
s |delta| < f_s       (delta<0, s>0).               (6)
```

Such an arc exists; the owner arcs can also be chosen pairwise disjoint.

For `e=s+kq`,

```text
floor(e(theta_0+delta)) mod P
 = b_s+floor(f_s+e delta) mod P.                    (7)
```

The inequalities (6) preserve `b_s` at the starting time `e=s`.  Since
`q|delta|<1`, consecutive values on a congruence track change by at most one
integer, so no sector can be skipped.

On the positive side, track `s` has reached missing sector `r` by time `n`
exactly when

```text
delta >= c^+_(s,r)/E_s(n).                          (8)
```

On the negative side, with `eta=-delta>0`, it has reached `r` exactly when

```text
eta > c^-_(s,r)/E_s(n).                             (9)
```

The strict sign in (9) is forced by `floor(z)<=-d iff z<-d+1`.  Taking the
inner minimum over tracks and outer maximum over missing sectors proves (1).
It also fixes the half-open convention: the negative-radius endpoint remains
noncover, while the positive-radius endpoint is covering.  Changing either
endpoint changes no measure.

The coefficient `c^-` can be zero.  This is not a division or topology
failure: it says every negative perturbation, however small, immediately
supplies that missing sector.  The owner itself remains a measure-zero
noncover point.

## 4. Compactness closes the global gap

Shrink the disjoint guard arcs from Section 3 once more, and call the smaller
arcs `V_b`.  Their complement `K` is compact and contains no persistent
owner.  By Section 2, every point of `K` has a neighbourhood covered at some
finite horizon.  A finite subcover of `K`, followed by the maximum of those
horizons, gives one finite horizon covering all of `K`.

Meanwhile every numerator in (1) is at most `P`, and
`E_s(n)>=n-q+1`.  Hence every radius tends to zero.  At a second finite
horizon all owner radii lie inside their corresponding `V_b`.  Beyond the
maximum of the two horizons, Section 3 describes every local noncover and
`K` contains none.  This proves the disjoint union and (2).

This compactness step proves existence for each fixed `P`; it does not hide a
uniform cutoff.  Making it effective would require quantitative cover
witnesses on `K`, which the theorem does not claim.

## 5. Finite-phase selector lemma

Fix `n mod L_P`.  For every `q<P` and `s<q`,

```text
E_s(n)=n-c_s,  0<=c_s<=q-1,                         (10)
```

with `c_s` constant on the phase.  Each candidate in (1) is therefore
`C/(n-c)`.  For two candidates, cross-multiplication gives an affine
comparison:

```text
C/(n-c) <= D/(n-d)
iff (C-D)n <= Cd-Dc.                                (11)
```

It can switch at most once unless it is an identity.  There are only finitely
many owner, side, missing-sector, and track candidates.  Past the last
pairwise crossing on a fixed phase, their entire order is fixed; hence the
nested max-min in (1) selects one fixed candidate.  This proves (3), including
the claimed finite but phase- and `P`-dependent onset.

The sum of the selected coefficients, divided by `P`, is `C_P` from Section
6.  Since every offset lies between `0` and `P-2`, termwise comparison gives
(5).

## 6. Closed prime constant

For a fixed denominator `q<P`, multiplication by `Pp` permutes the residues
modulo `q`.  Thus every reduced numerator has the same unordered theta-grid

```text
y_j=Pj/q,  0<=j<q.                                  (12)
```

The grid spacing is `h=P/q>1`.  The positive leading radius is the greatest
clockwise distance from a grid point to the left boundary of a missing
sector.  In a gap ending at `y_(j+1)`, the last possible missing-sector
boundary is `floor(y_(j+1))-1`, at distance

```text
h-1-frac(y_(j+1)).                                  (13)
```

The fractional parts are `0,1/q,...,(q-1)/q`; the cyclic gap ending at `0`
attains the maximum.  Therefore

```text
K_q^+=(P-q)/q.                                      (14)
```

For negative drift, the first missing sector after `y_j` has right boundary
`floor(y_j)+2`.  Its counterclockwise distance from `y_(j+1)` is

```text
h-2+frac(y_j).                                      (15)
```

The maximum fractional part `(q-1)/q` attains the maximum (zero is allowed at
`q=P-1`), so

```text
K_q^-=(P-q-1)/q.                                    (16)
```

As `n` tends to infinity, `n/E_s(n)->1`; finite max and min commute with this
limit.  Thus one owner of denominator `q` contributes

```text
(K_q^++K_q^-)/P=[2(P-q)-1]/(Pq)                    (17)
```

to `lim mD_P(m)`.  There are `phi(q)` such owners, proving (4).

## 7. Small primes and zero-width sides

- `P=2`: `B_2={0}`, `rho^+=1/n`, and `rho^-=0`.  In fact
  `D_2(m)=1/[2(m-1)]` for every `m>=2`, and `C_2=1/2`.
- `P=3`: the owners are `0` and `1/2`.  If `E_1(n)` is the largest odd
  integer at most `n`, the exact specialization from `m>=3` is

```text
D_3(m)=1/n+1/[6E_1(n)],
C_3=7/6.                                            (18)
```

  The denominator-two owner has `rho^-=0` exactly.
- In general, every denominator `q=P-1` owner has a zero negative radius.
  Its unique missing sector is hit at time `q` by every negative perturbation,
  while the rational owner itself stays noncover.  Formula (2) is a measure
  formula and handles this singleton correctly.

Thus neither `P=2` nor `P=3` requires an exception to the theorem; they only
expose why zero coefficients and half-open endpoints must be allowed.

## 8. Proved local cutoff half and open attraction lemma

The canonical audit script uses explicit `require`, so `python`
and `python -O` execute the same checks.  It independently:

1. enumerates all owner max-min coefficients and the totient formula through
   prime `31`;
2. checks the two-sided `q=P` residue trace through prime `31`;
3. compares the owner formula with a direct rational-wall cover engine for
   `P=2,3,5,7,11,13,17,19,23`;
4. checks cutoff-guard containment and positive owner nonoverlap through
   prime `31`; and
5. preserves `P=7,m=12` and composite modulus four as hostiles.

There is now a proved local half of the sharp-cutoff pattern. For odd prime
`P`, put

```text
n_0=max_(1<=q<P) q(P-q)=(P^2-1)/4,
m_0=n_0+1=(P^2+3)/4.                                (19)
```

For an owner of denominator `q>=2`, the positive and negative start/no-skip
guards satisfy

```text
g_q^+,g_q^- >=1/[q(q-1)].                            (20)
```

Indeed every nonzero fractional part `f_s` is in
`{1/q,...,(q-1)/q}`, while `1<=s<=q-1`. Section 6 and
`E_s(n)>=n-q+1` give

```text
rho^+ <=(P-q)/[q(n-q+1)],
rho^- <=(P-q-1)/[q(n-q+1)].                          (21)
```

At `n=n_0`, the positive radius is at most its guard because

```text
(q-1)(P-q+1)<=floor(P^2/4).                          (22)
```

The negative radius is strictly inside its guard because
`(q-1)(P-q)<q(P-q)<=n_0`. The case `q=1` is direct. Thus the included
negative endpoint remains in the exact track region, while equality on the
positive side is harmless because that endpoint is excluded. Every proposed
half-open owner piece is therefore guarded at `n_0`, and remains so at all
later horizons. Middle denominators make `(22)` sharp and explain the
quadratic scale.

The only missing all-prime implication is the following **OPEN finite
Sturmian attraction lemma**. Let `G_(p/q)` be the one-sided guarded cylinder
around `Pp/q` on which all `q` starting representatives are unchanged. Then

```text
noncover at n_0 and theta off every floor wall
  => theta lies in a unique G_(p/q) with q<P.          (23)
```

A repeated-sector pigeonhole gives some `q<P` with `|q theta-Pp|<1`, but
does not preserve the starting representatives: the intercept can absorb
the return error. Exact hostile prefixes already occur for `P=5` on
`theta in (4/3,3/2)`. A proof of `(23)` must retain the grid alignment in
the explicit three-gap/Farey parameterization. If `(23)` holds, the local
result above proves the owner formula from `m_0` onward. It does not by
itself prove that every earlier row fails.

Through its audited horizons, the direct engine finds sampled onset `m_0` for
the tested pairs

```text
(3,3), (5,7), (7,13), (11,31), (13,43), (17,73),
(19,91), (23,133),
```

with `P=2,m_0=2` as the even edge.  Every tested `m` from `P` through
`m_0-1` is a genuine mismatch, while all tested values from `m_0` onward
agree.

Primality is load-bearing. At modulus `H=4`, the denominator-two owner has a
nonzero track starting on a wall; negative drift loses an occupied base
sector. The naive prime formula undercounts by

```text
1/(4E_even)                                            (24)
```

through the audited tail and is not an eventual composite theorem.

Run:

```powershell
python 04-computation/prime_sector_ap_cover_eventual_owner_tail_thm4033.py
python -O 04-computation/prime_sector_ap_cover_eventual_owner_tail_thm4033.py
```

The frozen normal and optimized outputs are byte-identical.

```text
script sha256 = 77c3de5f38e553839dbabb44cd71bf405fba3dd0389a63a8e633e0b4397fa33b
output sha256 = 39c330e077a0a66c74eb6358679459b372035888484fa9e5627f8d3c5fa88025
```

The symbolic argument proves the stated eventual theorem; the computation is
regression evidence for its formulas and boundaries. **QED.**
