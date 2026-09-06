# An admissible singular wall in the original shared-root charts

**Status: PROVED exact chart identities + FINITE-EXACT algebraic witness + INDEPENDENTLY AUDITED.** Solving the original phase equation for x introduces genuine singular walls. One such wall contains a point with positive simple beta roots, both required interlacers, and z>0. It cannot be discarded by dividing by the coefficient of x. At its selected original phase, the full Laurent response is strictly negative. General sign on the shared-positive-root stratum remains OPEN.

This is an admissible point of the anchored coefficient model. No integer-support Laurent realization or LRC instance is claimed.

[Independent proof and exact audit](continuing6_20260906_shared_wall_audit.md) passes.

## 1. Inheritance and the precise question

The model is the one in [continuing3_20260906_residue_floor.md](continuing3_20260906_residue_floor.md) and [continuing5_20260906_moments_rectangle.md](continuing5_20260906_moments_rectangle.md):

```
B(v)=v^5-13v^4+55v^3-xv^2+yv-z,
C(v)=v^4-12v^3+45v^2-(2x/3)v+3y/7,
D(v)=v^4-11v^3+36v^2-(5x/12)v+y/7.
```

Beta roots are nonnegative with sum13 and square sum59, and C,D weakly interlace B. The original positive phase s satisfies

`g(s)=z s^4-(12/7)y s^3+x s^2-10s+1/11=0`.

The current parent session has separately handled z=0 and repeated beta roots; their new proofs are not used here. The immediate remaining coordinate is a shared positive root r of B with C or D. The closest proved mechanism is exact same-zero elimination in the rectangle proof; its corrected near miss is switching to a transformed first-zero condition. The least-used sidecar here is the zero coefficient of the affine x equation, together with its constant term.

The board is: the original phase; a shared positive beta root; the two linear chart equations; the coefficient denominator wall; full root/interlacing admissibility; and the original response with its inverse carry. The source-to-target map eliminates y,z using B(r)=C(r)=0 or B(r)=D(r)=0. It preserves those exact equalities. Solving the remaining phase for x can destroy an entire singular fibre. The necessary sidecar is the simultaneous coefficient/constant system before any division. The cheapest decisive test is one exact admissible point on that fibre.

## 2. Both original charts and their singular denominators

For r>0, the C-shared chart is exactly

```
y=(14/9)xr-(7/3)r^2(r^2-12r+45),
z=r^2[5x/9-r(4r^2-45r+150)/3].
```

The D-shared chart is exactly

```
y=(35/12)xr-7r^2(r^2-11r+36),
z=r^2[23x/12-r(6r^2-64r+197)].
```

These solve the shared-root equations in both directions. Put u=rs. On the C chart,

`g(s)=x u^2 A_C(u)/(9r^2)-H_C(r,u)/(33r)`,

where

```
A_C=5u^2-24u+9,
H_C=44r^2u^4-132r^2u^3-495ru^4+1584ru^3-3r
    +1650u^4-5940u^3+330u.
```

On the D chart,

`g(s)=x u^2 A_D(u)/(12r^2)-H_D(r,u)/(11r)`,

where

```
A_D=23u^2-60u+12,
H_D=66r^2u^4-132r^2u^3-704ru^4+1452ru^3-r
    +2167u^4-4752u^3+110u.
```

Thus the proposed coefficients s^2(5(rs)^2-24rs+9)/9 and s^2(23(rs)^2-60rs+12)/12 are correct. Away from A_C=0 or A_D=0, one may solve respectively

`x=3r H_C/(11u^2 A_C)` or `x=12r H_D/(11u^2 A_D)`.

On the denominator wall, this division is invalid. If H is nonzero there is no original phase. If H also vanishes, the phase equation imposes no additional condition on x along that shared-root chart. The latter possibility really occurs.

## 3. The positive singular pairs are finite and exactly isolated

Eliminating u from A_C=H_C=0 gives, up to a nonzero rational factor,

`P_C(r)=705672r^4-15079284r^3+115450055r^2-375696750r+441317250`.

Its u coordinate is recovered by

`u=(11286r^2-109085r+245025)/(27126r^2-261360r+586025)`.

For D the corresponding formulas are

```
P_D(r)=120434688r^4-2310536448r^3+15142183583r^2
       -39453584808r+35022385200,
u=(446688r^2-4023865r+7744176)
  /(1978416r^2-17739216r+34339426).
```

Both reconstruction denominators are coprime to their quartic. Substitution into A and H gives zero modulo that quartic. Exact Sturm counts and rational interval evaluation prove that there are exactly four positive pairs for each chart. Their r coordinates lie in the following rational intervals; each endpoint has denominator10^9.

| Chart | Lower numerator | Upper numerator |
|---|---:|---:|
|C|3528453826|3528453827|
|C|3550824804|3550824805|
|C|6081020695|6081020696|
|C|8208387541|8208387542|
|D|2049574992|2049574993|
|D|2831681376|2831681377|
|D|6130326664|6130326665|
|D|8173391713|8173391714|

Every recovered u is positive, and s=u/r is therefore positive. The two largest r values exceed the beta-root bound71/10 and are inadmissible. Indeed, Cauchy on the other four roots gives M^2+(13-M)^2/4<=59, hence M<71/10. No admissibility claim is made for all the other six. The next section certifies one of them in full.

## 4. A positive, simple, interlacing witness on the C wall

Let r be the unique root of P_C in

`(3528453826758,3528453826759)/10^12`.

Let u be the C rational reconstruction above, set s=u/r, and choose

```
x=(3/5)r(4r^2-45r+150)+1/10,
y=(14/9)xr-(7/3)r^2(r^2-12r+45),
z=r^2/18.
```

The last formula follows exactly from the C chart. Thus z>0, both shared-root equations hold, A_C(u)=H_C(r,u)=0, and the **original** phase equation g(s)=0. The approximate value s=0.1162053256 is descriptive only; every gate uses rational intervals or exact polynomial identities.

After removing the known factor(v-r), the beta quotient has one root in each of

`(18,20)/1000, (669,671)/1000, (2450,2452)/1000, (6331,6333)/1000`.

The C quotient has one root in each of

`(388,390)/1000, (1951,1953)/1000, (6130,6132)/1000`.

D has one root in each of

`(181,183)/1000, (1493,1495)/1000, (3383,3385)/1000, (5939,5941)/1000`.

For each bracket the source evaluates the endpoint polynomial as a rational interval over the entire isolated r interval. The signs are strictly opposite. The number of disjoint brackets equals the degree, so these are all the roots and each is simple. The known r is distinct from every beta quotient root. Their ordering proves that B has five distinct positive roots, C weakly interlaces B with its third root equal to the fourth beta root r, and D strictly interlaces B.

The beta coefficients automatically give the anchors sum13 and square sum59. The same exact intervals also show

`86<x<87, 0<y<40, 0<z<1`.

There is also a nontrivial admissible x-interval around the displayed value.
Keep r,u fixed and vary x along the C-shared chart: both shared-root equations
and the original phase persist identically. The strict rational root-bracket
signs persist by continuity, so B stays positive/simple, C retains its shared
root and other strict brackets, and D stays strictly interlacing. The negative
response likewise persists on a possibly smaller interval. No explicit radius
or sign on the entire singular fibre is claimed.

Thus even strong basic coefficient caps and full root admissibility do not remove this wall. No converse realization as an integer-support source is inferred.

## 5. The original response is negative at this wall point

The verifier independently reconstructs the inherited full Laurent model. With

```
beta=t^-1+13+55t+xt^2+yt^3+zt^4,
c=t^-1+12+45t+(2x/3)t^2+(3y/7)t^3,
d=t^-1+11+36t+(5x/12)t^2+(y/7)t^3,
```

the raw second response is beta^2+2tcd. Its coefficientwise carrier is the sum of the odd-binomial square and the shifted even-binomial square for length14. The resulting Q has the full Laurent support-1 through8, including coefficient28 at t^-1. The first polynomial is checked to satisfy P(-s)=2002g(s), so no transformed first-root condition is substituted.

Direct rational interval evaluation of this complete Q at t=-s proves

`-5634725<Q(-s)<-5634723`.

The negative value is a positive control for the desired sign theorem, not a sign obstruction. The obstruction is to removing the singular fibre by dividing the phase equation. The response is checked at the one selected original phase; no unproved claim about every phase of this witness is needed.

## 6. Reproduction and the exact next question

```
python continuing6_20260906_shared_wall.py
python -O continuing6_20260906_shared_wall.py
```

The script uses SymPy for exact polynomial identities and Sturm counts, and standard-library rational intervals for admissibility and the full response. Normal and optimized runs pass65 always-active gates with identical actual LF output. No producer implementation is imported. The certificate records the defining equations, every rational isolation and sign bracket, and exact rational response bounds.

- Source SHA256:`ef604a3c9276171b844f6e260bc788a17160ad6e85a70f3fbd131ee4f3c362c4`.
- Output SHA256:`4803b7577514138cbe0fca4ef339b0e2cb690ba421e2276eb823846ef0d33572`.
- Certificate SHA256:`6adab7b43f30296cdac428596ce7c046d163751b6d48fea4db02fe17a100dce4`.

The precise next object is the remaining x-dependent response on these finite singular(r,u) fibres, restricted to root/interlacing admissibility. The generic charts and the singular fibres must be treated separately. This note establishes one admissible fibre point with the desired sign, and stops there; the full shared-positive-root sign theorem remains OPEN. No clock scan, maintained-file edit, or theorem-ID claim was made.
