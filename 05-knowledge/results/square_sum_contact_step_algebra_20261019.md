# Small-translate formulas and the obstruction to a uniform bound

Status: proved symbolically. This concerns the specified anchor point and the
specified translating point, not all possible direct rational constructions.
No external theorem beyond the elliptic group law is needed.

Put

- f(u)=4u(u²−1)/(u²+1)²;
- k=(t²+2t−1)/(t²+1), K=k²;
- E: Y²=(X−K)(X²−1), P=(0,k);
- Q=((t²+4t+1)/(t²+1),8t(t+1)/(t²+1)²).

The point Q is the image of the anchor (u,z)=(t,t²+1). Let u_q(t),z_q(t)
be the inverse quartic coordinates of Q+qP. All formulas below are rational
identities on their ordinary charts and extend by continuity where stated.

## 1. Direct translation law on the quartic

On z²=K(u²+1)²−4u(u²−1), put a=K(u²−1)−u. Then

u(R+P)=(a+kz)/(u(u−2K))=(2Ku+1)/(a−kz),

u(R−P)=(a−kz)/(u(u−2K))=(2Ku+1)/(a+kz).

The equality of each pair follows from

a²−k²z²=u(u−2K)(2Ku+1).

If v=u(R+P), its z-coordinate is

z(R+P)=[K(v²−1)−v−uv(v−2K)]/k.

If v=u(R−P), its z-coordinate is the negative of this expression.
These formulas follow either from the group law and birational map or by
solving the inverse translation equation for z'. The companion checker
verifies them against the original elliptic formulas, not just against the
quartic equation.

## 2. Exact q=1 and q=2 formulas

u_1=(2t³+t²−2t+1)/(t(t²−4t+1)),

z_1=−(t⁶+24t⁵−13t⁴−16t³+19t²−8t+1)/(t²(t²−4t+1)²),

u_−1=(t²+4t+1)/(t³+2t²+t−2),

z_−1=−(t⁶+8t⁵+19t⁴+16t³−13t²−24t+1)/(t³+2t²+t−2)².

For q=2, define

A=2t³+t²−2t+1,
B=t²−4t+1,
C=5t⁷+14t⁶+15t⁵−34t⁴−t³+18t²−11t+2,
D=t⁶+28t⁵+19t⁴−40t³+19t²−4t+1.

Then u_2=BC/(AD), and

u_−2=−(t³+2t²+t−2)(t⁶+4t⁵+19t⁴+40t³+19t²−28t+1)
       /((t²+4t+1)(2t⁷+11t⁶+18t⁵+t⁴−34t³−15t²+14t−5)).

## 3. Universal failure at q=1

Every admissible t lies in (1,4/3) or (7,∞). Indeed f is increasing up to
1+sqrt(2), decreasing thereafter, and f(4/3)=f(7)=336/625>1/2.

For 1<t<4/3:

- u_1<−1, because u_1+1=(t−1)(3t²−1)/(t(t²−4t+1))<0;
- 3/2<u_−1<3. The upper bound follows from
  3(t³+2t²+t−2)−(t²+4t+1)=(t−1)(3t²+8t+7)>0.
  For the lower bound, 2(t²+4t+1)−3(t³+2t²+t−2)
  =−3t³−4t²+5t+8 is decreasing and at t=4/3 equals 4/9>0.

For t>7:

- 0<u_−1<1, because denominator minus numerator=(t+1)(t²−3)>0;
- 2<u_1<5. The first difference has numerator 9t²−4t+1>0.
  The second has numerator 3t³−21t²+7t−1, positive at 7 and increasing
  thereafter.

On (−∞,−1) and (0,1), f is negative. On [3/2,3] and [2,5], f is at
least 120/169>1/2, by the same unimodality. Hence neither Q+P nor Q−P
can yield an admissible neighbor, even if a negative u is normalized using
u↦−1/u, which preserves f.

## 4. A rigorous q=2 success patch

Every rational anchor 24≤t≤32 has an admissible distinct neighbor at Q+2P.
Indeed t itself is admissible because 0<f(t)<4/t≤1/6. The denominator AD is
positive for t≥24. The derivative of u_2 has denominator A²D² and numerator

297t¹⁶+784t¹⁵+408t¹⁴+9072t¹³−13268t¹²−19632t¹¹+35208t¹⁰
−10192t⁹−11354t⁸+2352t⁷+15336t⁶−18352t⁵+10124t⁴−2960t³
+376t²+16t−7.

In powers of t−2 its coefficients, descending, are

297,10288,166488,1680336,11873516,62316528,250947144,788978960,
1951478374,3798930064,5786178600,6808009584,6055240652,
3929641296,1752327096,479014832,60367929.

Thus u_2 is strictly increasing for t≥2 wherever defined. Direct exact
endpoint evaluations give

u_2(24)=12369373894202/11829214068337>1,
u_2(32)=168006958644642/135115086262081<5/4.

Consequently 1<u_2<5/4 and 0<f(u_2)<f(5/4)=720/1681<1/2.
The contact equation is r²=k²−f(u_2)=1+f(t)−f(u_2).

To verify distinctness, the numerator of f(u_2)−f(t), up to a nonzero
positive factor and square denominator, is

(7t⁴−16t³+6t²−1)
(t⁴+8t³+2t²−8t+1)
(3t⁴+8t³+6t²−8t+3)
(t⁶−26t⁵+17t⁴+12t³−17t²+6t−1)
(t⁶+24t⁵−13t⁴−16t³+19t²−8t+1).

The first three factors are positive for t≥24. Each of the last two is
monic with constant ±1; a rational zero must be ±1. Neither vanishes at a
rational t≥24. Hence the neighbor's spread is distinct throughout this
rational patch. The patch is only a sufficient interval, not claimed maximal.

## 5. No finite uniform translate bound

Let Pell numbers be p₀=0, p₁=1, pₙ₊₁=2pₙ+pₙ₋₁. For every fixed q≥1,

lim[t→∞] (u_q,z_q)=(p_(q+1)/p_q,h(p_(q+1)/p_q)),
lim[t→∞] (u_−q,z_−q)=(−p_(q−1)/p_q,h(−p_(q−1)/p_q)),

where h(u)=u²−2u−1.

Proof: k→1, and the limiting quartic factors as z²=h(u)². The exact q=1
formulas give the bases (2,−1) and (0,−1), both on z=h(u). On that component,
the translation laws reduce to

F(u)=2+1/u for +P, and G(u)=1/(u−2) for −P.

The apparent cancellation at u=2 in the original +P chart causes no issue:
the conjugate chart denominator a−kz tends to u, which is nonzero at every
positive Pell ratio. For −P, the original chart is regular at every negative
Pell ratio. At the initial u=0, its conjugate denominator tends to −2 and
is regular. The polynomial formula for z' gives h(F(u)) and h(G(u)) in the
limits, so the branch remains z=h. Thus induction proves all stated limits.
Every fixed finite list of translations is defined on these charts for all
sufficiently large t.

The Pell identity

p_(q+1)²−2p_(q+1)p_q−p_q²=(−1)^q

and the identity 1−f(u)=h(u)²/(u²+1)² give

lim f(u_q)=1−1/(p_(q+1)²+p_q²)²≥24/25, for q≥1;
lim f(u_−q)=1−1/(p_(q−1)²+p_q²)²≥24/25, for q≥2.

For the remaining translate q=−1, its exact formula gives 0<u_−1<1 for
t>sqrt(3), so f(u_−1)<0, although its limit is zero.

Therefore, for every finite B, all 0<|q|≤B fail admissibility for every
sufficiently large rational t. Such anchors exist, e.g. arbitrarily large
integers t. The q=0 point has the anchor's own spread and is excluded.

Conclusion: the least successful |q| in this prescribed elliptic search is
unbounded over admissible rational anchors. This obstruction persists after
all spread-preserving reparametrizations and primitive normalization. It does
not rule out bounded translate counts depending on height or spread, an
alternative starting point, a different translating point, or a separate
rational construction. It does not contradict eventual success for each
fixed anchor from real-group recurrence.
