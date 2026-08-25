---
id: THM-4033
title: "Prime-sector AP-cover sharp finite-owner tail"
status: >
  PROVED + VERIFIED-EXACT. For every fixed prime P, the AP sector-cover
  deficit has a finite-owner decomposition over exactly the
  reduced rationals of denominator below P, is eventually phase-rational
  with period dividing lcm(1,...,P-1), and has an explicit positive
  1/m constant and two-sided rational bounds. For every odd prime the
  decomposition holds from m_0=(P^2+3)/4. For P>=5 an explicit positive
  open interval, uniformly defined in P and outside every owner piece,
  proves geometric failure at m_0-1. For P=3, m_0=3 is the lower-domain
  edge because the formula is undefined at m_0-1. The P=2 edge is m_0=2.
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
  - THM-4038-ap-deficit-holonomic-sixty-phase-law
script: 04-computation/prime_sector_ap_cover_eventual_owner_tail_thm4033.py
output: 05-knowledge/results/prime_sector_ap_cover_eventual_owner_tail_thm4033.out
script_sha256: 8c73998c615d3709819898494ffede3a773c2efaf415b1cc3f4d6e3e5c923c61
output_sha256: 3670f67e6c3b4bfb1621f6b7f77c944ca3693b22563e6d95e0fdddac552e04af
hash_basis: raw LF bytes
---

# THM-4033 -- prime-sector AP-cover sharp finite-owner tail

**PROVED + VERIFIED-EXACT.** The persistent-owner classification, sharp
three-gap attraction lemma, finite-phase selector lemma, and closed constant
below are proved.  The owner formula holds from the stated onset for every odd
prime.  Its geometric failure one row earlier is proved for odd `P>=5`;
`P=3` is the lower-domain edge.  The stronger claim that the scalar formula
mismatches at every earlier horizon remains only finite-exact.

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

For an odd prime `P`, put

```text
n_0=(P^2-1)/4,  m_0=n_0+1=(P^2+3)/4.
```

For every `m>=m_0`, the noncover set, modulo finitely many wall points, is
the disjoint union of the owner pieces

```text
[theta_0-rho^-_(p/q)(n), theta_0+rho^+_(p/q)(n))
```

on the length-`P` theta circle.  In particular,

```text
D_P(m)=(1/P) sum_(p/q in B_P)
             (rho^-_(p/q)(m-1)+rho^+_(p/q)(m-1)).   (2)
```

For `P>=5`, the geometric owner decomposition fails at `m=m_0-1`: an
explicit positive open interval, uniformly defined in `P`, is still noncover
and lies outside every owner piece.  For `P=3`, `m_0=3` is instead the
lower-domain edge, since (1) has a zero denominator at `m_0-1`.  This does not
assert that the scalar sum in (2) fails at every `P<=m<m_0`; that stronger
statement is finite-exact through the audited primes.  For `P=2`, (2) holds
for every `m>=2`.

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

Thus `D_P` is eventually `L_P`-phase rational.  The geometric part has the
explicit uniform-in-`P` formula above; the selector stabilization within each
phase is still only asserted with its finite phase-dependent onset.

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

## 4. Sharp three-gap attraction closes the global gap

### 4.1 Maximal guards and their disjointness

Let `q>=2`, write `A=Pp`, and let

```text
u/s < A/q < v/t,  (u,s)+(v,t)=(A,q),  vs-ut=1
```

be the two Stern--Brocot parents of `A/q` (with the determinant orientation
chosen to match the displayed order).  The maximal strict arc on which all
starting representatives in (6) are unchanged is exactly

```text
G_(p/q)=(u/s,v/t).                                  (G1)
```

Indeed the endpoint distances from `A/q` are `1/(qs)` and `1/(qt)`, and
the Farey-parent property says that no floor wall of denominator below `q`
lies between an endpoint and the centre.  The `q=1` guard is the unit arc
`(A-1,A+1)`.  Thus (G1) is precisely the maximal version of (6), not an
extra neighbourhood chosen by compactness.

These owner guards have disjoint interiors.  To see the arithmetic reason,
write their primitive parent vectors as `U=(u,s)`, `V=(v,t)`.  Every reduced
rational in the parent cylinder has a unique expression

```text
W=xU+yV,  x,y positive integers.
```

Since `u+v=A=0 (mod P)` and unimodularity forces `u!=0 (mod P)`, the numerator
of `W` is divisible by `P` exactly when `x=y (mod P)`.  If `x=y`,
primitivity gives `x=y=1` and recovers the centre.  Otherwise
`|x-y|>=P`, so the denominator `xs+yt` is greater than `P`.  Hence no
distinct owner of denominator below `P` lies in an owner cylinder.
Stern--Brocot parent cylinders are nested or disjoint, so the owners form an
antichain and their guard interiors are disjoint.  On the length-`P` circle,
the `q=1` guard is also separate: the lower Stern--Brocot parent of every
other positive owner is at least `1` (the nearest guard begins exactly at
`1`), with the reflected statement at the other seam.

### 4.2 Exact large-gap attraction lemma

Put `P=2h+1` and

```text
Q=h(h+1)+1=(P^2+3)/4.                              (G2)
```

Normalize `alpha=theta/P` to the unit circle and consider the `Q` points
`0,alpha,...,(Q-1)alpha`.  We use the following exact form of the three-gap
lemma.  When these points are distinct, the nearest clockwise and
counterclockwise returns determine Farey neighbours

```text
a/q < alpha < b/r,
A=q alpha-a,  B=b-r alpha,
rA+qB=1,  q+r>=Q.                                  (G3)
```

The gap lengths and their multiplicities are

```text
A : Q-q,   B : Q-r,   A+B : q+r-Q.                 (G4)
```

This is the standard subtractive proof of the three-gap theorem: translate a
gap until one endpoint reaches the origin; the two origin gaps have lengths
`A,B`, every other gap is one of them or their sum, and counting their
allowable translates gives (G4).  The Farey determinant gives the identity in
(G3).

We claim that every gap strictly longer than `1/P` forces `alpha` into one
of the guards (G1), after dividing its theta endpoints by `P`.

First suppose an `A+B` gap occurs and `A+B>1/P`.  Then `q+r>Q`.
Both denominators cannot be at least `P`, since (G3) would give

```text
1=rA+qB >= P(A+B)>1.
```

Except for the immediate `P=3,5` check below, both cannot be below `P`
either, because `q+r<=2P-2<Q`.  Assume `q<P<=r`; the other orientation is
the reflection.  For `q>1`, let `t in {1,...,q-1}` be the right-parent
denominator of the theta owner `Pa/q`, so

```text
Pa t == -1 (mod q).
```

The right guard condition is `A<1/(Pt)`.  If it failed, then (G3) and
`A+B>1/P` would give

```text
q(A+B)=1-(r-q)A>q/P,
r-q<t(P-q).                                        (G5)
```

The Farey determinant in (G3) also gives `ar==-1 (mod q)`, hence
`r==Pt (mod q)`.  Write `r=Pt-kq`.  Inequality (G5) forces `k>=t`, so

```text
q+r <= q+t(P-q) <= q+(q-1)(P-q) <= Q.              (G6)
```

The last inequality is exact: for integral `q`,

```text
Q-[q+(q-1)(P-q)]=(q-h-1)(q-h-2)>=0.
```

This contradicts `q+r>Q`.  For `q=1`, failure of its guard means
`A>=1/P`, already incompatible with `r>=P` in (G3).  Thus the small
Farey endpoint attracts `alpha` into its owner guard.

It remains to treat `q+r=Q`, when the `A+B` multiplicity in (G4) is zero.
Suppose `A>1/P`; the `B` case is symmetric.  Equation (G3) gives `r<P`.
For `r>1`, let `s` be the left-parent denominator of the upper theta owner
`Pb/r`, so `Pb s==1 (mod r)`.  Outside its guard, `B>=1/(Ps)`.  But

```text
qB<(P-r)/P,  q==Ps (mod r).
```

Writing `q=Ps-kr` forces `k>=s+1`, and therefore

```text
Q=q+r <= s(P-r) <=(r-1)(P-r)<=h^2<Q,               (G7)
```

a contradiction.  If `r=1`, then `q=Q-1` and

```text
B<(P-1)/[P(Q-1)]<=1/P,
```

so the denominator-one guard contains `alpha`.  For `P=3,5`, the only way
both denominators below `P` could have `q+r>Q` forces equal denominators,
contrary to the Farey determinant; the zero-sum case is already covered by
(G7).  This proves the attraction lemma for every odd prime.

If `alpha` is rational with reduced denominator `d<Q`, the first `Q` orbit
points contain the full `d`-grid.  For `d<P` this is an owner; for `d>=P`
its gaps are at most `1/P`.  The boundary `d=Q` is likewise the full
`Q`-grid, and `Q>=P`.  Thus rational repetitions and the equality boundary
introduce no omitted case; the distinct-point argument applies to irrational
`alpha` and to rational denominator greater than `Q`.

### 4.3 From attraction to the exact onset

A missing half-open sector `[j/P,(j+1)/P)` forces a gap strictly longer than
`1/P`: its left endpoint cannot be occupied, while the successor may only
equal the excluded right endpoint.  Section 4.2 therefore puts every
noncover at `m=Q` into a unique maximal owner guard.  Section 3 gives its
exact half-open piece there, and Section 4.1 lets the lengths sum.  For
`m>Q`, a noncover was already a noncover at `Q`, so the same guarded local
description applies.  This proves (2) for every `m>=m_0=Q`.

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

There is also a qualitative holonomic corollary.  For each fixed prime `P`,
the sequence `D_P(m)` is P-recursive over the rationals, equivalently its
ordinary generating function is D-finite.  Indeed, on each eventual residue
class modulo `L_P`, (3) makes `D_P` a finite sum of rational functions of the
phase index, hence P-recursive.  Finite interlacing of the residue classes is
D-finite, and changing the finite prefix preserves D-finiteness.  In
particular `P=7` has a proved qualitative holonomic law from a phase modulus
dividing `60`.  This corollary by itself supplies neither an explicit
recurrence nor a minimal phase.  THM-4038 independently supplies exactly that
`P=7` refinement: from `n=12` (the present sharp onset) it proves minimal
phase `60`, an explicit polynomial-coefficient recurrence, and a D-finite but
nonalgebraic ordinary generating function.  No corresponding minimal-phase
claim for general prime `P` is made here.

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

## 8. Sharpness one row before onset

### 8.1 Every proposed piece is still guarded

Let `P=2h+1>=5` and inspect the row immediately before onset:

```text
n_-=m_0-2=h(h+1)-1.
```

For an owner of denominator `q>=2`, its two guard radii satisfy

```text
g_q^+,g_q^- >=1/[q(q-1)].                            (S1)
```

Indeed the nonzero `f_s` are multiples of `1/q`.  Section 6 and
`E_s(n)>=n-q+1` give

```text
rho^+ <=(P-q)/[q(n-q+1)],
rho^- <=(P-q-1)/[q(n-q+1)].                          (S2)
```

Every negative radius is strictly guarded at `n_-`, because

```text
q+(q-1)(P-q-1)<=h^2+1<h(h+1)=n_-+1.                 (S3)
```

The same argument guards the positive radius unless `q=h+1` or `q=h+2`:

```text
h(h+1)-[q+(q-1)(P-q)]
  =(q-h-1)(q-h-2)-1,                                (S4)
```

and away from those two integral values the product on the right before
subtracting `1` is at least `2`.

The two middle denominators are exact rather than exceptional.  If
`q=h+1`, then `P=2q-1`, `n_-=q^2-q-1`, and

```text
E_s(n_-)=q(q-2)+s.
```

For each missing sector, choose a track minimizing its leading positive
coefficient `c`.  Section 6 gives `c<=(q-1)/q<1`; track `s=0` has a positive
integral coefficient and cannot be chosen.  Hence `s>=1`,
`E_s>=(q-1)^2`, and `c/E_s<=1/[q(q-1)]`.

If `q=h+2`, then `P=2q-3`, `n_-=q^2-3q+1`, and

```text
E_0=q(q-3),  E_1=n_-,
E_s=q^2-4q+s  (2<=s<q).                             (S5)
```

Now a leading-minimizing track has `c<=(q-3)/q<1`, so again `s!=0`.
Tracks `s=1` and `s>=3` satisfy the desired guard bound directly.  For
`s=2`, if `c<(q-3)/q`, the `1/q` lattice gives `c<=(q-4)/q`, which also
suffices.  In the equality case, write `A=Pp (mod q)`.  Then
`c=1-f_2=(q-3)/q`, so `2A==3 (mod q)`.  If `t` is the right-parent
denominator, `At==-1 (mod q)` and the exact positive guard is `1/(qt)`.
The value `t=q-1` would force `A==1`, contradicting `2A==3`; hence
`t<=q-2`.  Finally

```text
(q-3)(q-2)<=q^2-4q+2=E_2,
```

so `c/E_2<=1/[q(q-2)]<=1/(qt)`.  For the denominator-one owner the radii
are `(P-1)/n_-` and `(P-2)/n_-`, both at most its unit guard radius (the
negative one strictly).  Thus every half-open owner piece at `n_-` is inside
its open guard: the negative endpoint is strictly inside by (S3), while
equality at a positive endpoint is harmless because that endpoint is
excluded.

### 8.2 An explicit uniformly defined owner-piece-free noncover interval

Put

```text
a=2h^2-h+1,  d=h^2-1,  c=2h^3-h^2-h+1,
L=a/h,  R=c/d,  G=(L,R).                            (S6)
```

The exact identities are

```text
hc-ad=1,  a==2 (mod P),  c==1 (mod P).              (S7)
```

We first show that sector `1` is absent throughout `G` for every
`1<=e<=n_-`.  If an interval `[j/e,(j+1)/e)` with `j==1 (mod P)` met `G`,
then the integers

```text
U=ce-jd,  V=hj-ae
```

would satisfy `U>=1`, `V>=1-h`, and unimodularity would give

```text
(j,e)=U(a,h)+V(c,d),
e=Uh+Vd,  2U+V==1 (mod P).                          (S8)
```

If `V>=1`, then `e>=h+d=n_-`; equality forces `U=V=1`, whose residue in
(S8) is `3`, not `1`.  If `V=0`, then `2U==1`, so `U>=h+1` and
`e>=h(h+1)>n_-`.  Finally write `V=-w` with `1<=w<h`.  Positivity forces
`U=wh+ell` with `ell>=0`, and then `e=ell h+w<=n_-` forces `ell<=h`.
The congruence in (S8) becomes

```text
2(ell-w)==1 (mod 2h+1),
```

which has no solution for `-(h-1)<=ell-w<=h-1`.  This proves that all of
`G` is noncover at `m_0-1`.  The next time `n_0=h(h+1)` starts the sector-one
interval exactly at `L`: `n_0L=(h+1)a==1 (mod P)`, and
`R-L=1/[h(h^2-1)]<=1/n_0`.

It remains to exclude owner pieces, not just owner centres.  The fraction
`L=a/h` is the right parent of the owner

```text
W_0=Ph/(h+1),  a(h+1)-Ph^2=1,
```

so the guard of `W_0` ends at `L`.  Also `G` is itself a Stern--Brocot
parent cylinder by (S7).  If another owner guard met `G`, laminarity would
force one cylinder to contain the other.  A guard contained in `G` would
put its owner centre there, but every reduced fraction strictly between
the Farey neighbours `L,R` has denominator at least

```text
h+d=h(h+1)-1=n_->=P,
```

so it cannot be an owner.  Conversely, a guard containing `G` would overlap
the guard of `W_0` unless its left endpoint were exactly `L`.  In that last
case, if its owner were `Pp/q`, the parent determinant would give

```text
Pph-aq=1,
```

hence `-2q==1 (mod P)` and therefore `q=h`; but a proper child of the parent
`a/h` has denominator greater than `h`.  This is impossible.  Thus `G`
meets no owner guard, and Section 8.1 shows that it meets no proposed owner
piece.  Since `G` is a positive open set of genuine noncovers, the explicit
geometric union fails at `m_0-1`.

For `P=3`, this is a lower-domain edge rather than the same interval
mechanism: at `m_0-1=2` only the two times `e=0,1` are present, so every
point is noncover, while the denominator-two track formula has `E_0(1)=0`
and the stated guarded decomposition is not valid.  This completes the
geometric sharpness proof for `P>=5` and identifies the `P=3` domain edge.

### 8.3 Exact audit and the composite hostile

The canonical audit script uses explicit `require`, so `python` and
`python -O` execute the same checks.  It independently:

1. enumerates all owner max-min coefficients and the totient formula through
   prime `31`;
2. checks the two-sided `q=P` residue trace through prime `31`;
3. compares the owner formula with a direct rational-wall cover engine for
   `P=2,3,5,7,11,13,17,19,23`;
4. checks onset and pre-onset guard containment, owner nonoverlap, the exact
   inequalities used above, and the interval (S6); and
5. preserves `P=7,m=12` and composite modulus four as hostiles.

The direct engine finds sampled onset `m_0` for

```text
(3,3), (5,7), (7,13), (11,31), (13,43), (17,73),
(19,91), (23,133),
```

with `P=2,m_0=2` as the even edge.  Every tested `m` from `P` through
`m_0-1` is a genuine scalar mismatch, while all tested values from `m_0`
onward agree.  Only this finite-exact range, not an all-prime all-earlier
scalar mismatch theorem, is claimed.

Primality is load-bearing. At modulus `H=4`, the denominator-two owner has a
nonzero track starting on a wall; negative drift loses an occupied base
sector. The naive prime formula undercounts by

```text
1/(4E_even)                                            (S9)
```

through the audited tail and is not an eventual composite theorem.

Run:

```powershell
python 04-computation/prime_sector_ap_cover_eventual_owner_tail_thm4033.py
python -O 04-computation/prime_sector_ap_cover_eventual_owner_tail_thm4033.py
```

The frozen normal and optimized outputs are byte-identical.

```text
script sha256 = 8c73998c615d3709819898494ffede3a773c2efaf415b1cc3f4d6e3e5c923c61
output sha256 = 3670f67e6c3b4bfb1621f6b7f77c944ca3693b22563e6d95e0fdddac552e04af
```

The symbolic argument proves the stated sharp-onset and eventual-selector
theorem; the computation is regression evidence for its formulas and
boundaries. **QED.**
