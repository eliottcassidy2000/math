---
id: THM-1226
title: THE GCD-PERIOD PROJECTIVE-CHARGE OBSTRUCTION — exact vertex loads, disconnected rescue, primitive phase fibers, and relation-cycle holonomy
status: PROVED (exact arithmetic obstruction and protected-needle embedding; conditional disconnected-G_gt transfer; sharp forest crown; AP-free strict-high/heavy-circuit separation; THM-605 phase-fiber refinement; primitive-coordinate and relation-cycle identities; THM-1237 scale/gcd dichotomy; centered blocker-cycle positive holonomy). Does not prove the remaining address-compression/alternate-gap descent or LRC(14)
source: codex-2026-07-19-S82
depends_on: [THM-605, THM-1166, THM-1218, THM-1221, THM-1233, THM-1236, THM-1237, THM-1240, THM-1241]
related: [THM-599, THM-773, THM-1234, THM-1238, THM-1242, MISTAKE-184, MISTAKE-185, HYP-7678, HYP-7870]
script: 04-computation/lrc14_gcd_period_projective_charge_obstruction_referee_codex_S82.py
output: 05-knowledge/results/lrc14_gcd_period_projective_charge_obstruction_referee_codex_S82.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCGCDPeriodProjectiveCharge.lean
script_sha256: 9ae71334d58a3d5613cb4197724e4af5d29c95f2e0d87e5a8ce6f9b0341f55a7
output_sha256: cfc88e17e0177a4288f3521ff3437a063a9ba68f6646d251adf41c3b39609503
formalization_sha256: f23ffb3b6e331f73e1be5549d092431d6e4fa93f7a0b3c7ea9c0d59b007031aa
---

# THM-1226 — the gcd-period projective-charge obstruction

## Statement

For an edge `s_i=g_ij a_ij`, `s_j=g_ij b_ij`, with
`g_ij=gcd(s_i,s_j)` and `(a_ij,b_ij)=1`, put

```text
eta_ij=rho_ij(1-rho_ij),
kappa_ij=eta_ij ab/(a+b).
```

Then the periodic positioning error factors exactly as

```text
eta_ij/g=kappa_ij(1/s_i+1/s_j).                         (1)
```

For a tree `T`, define the projective load at vertex `i` by

```text
Lambda_i(T)=sum_(ij in T) kappa_ij.                      (2)
```

Writing `H=sum_i 1/s_i`, the total THM-1166 period error is therefore

```text
E_T=sum_(ij in T) eta_ij/g_ij
   =sum_i Lambda_i(T)/s_i.                               (3)
```

There is no absolute finite `C` such that every seven-speed packet admits a
tree satisfying both

```text
sum_(ij in T) rho_ij >=15/154,
E_T<=C H.                                                (4)
```

This remains false even when **every** spanning tree clears the first
inequality and the strict-high graph `G_gt={rho>1/63}` is complete.

The exact finite projective branches of THM-1221 do give a complementary
positive theorem.  If `G_gt` is disconnected, some THM-1221 floor tree obeys

```text
E_T<=C_disc H,        C_disc=448916/194775.              (5)
```

Consequently every interval `I` of length `L` satisfies the conditional
positioned-tree bound

```text
sum_(ij in T) mu(I intersect D_si intersect D_sj)
 >=(15/154)L-(448916/194775)H.                           (6)
```

If the seven combs cover `I`, THM-1166's forest-Hunter cap improves the
reverse tree estimate from `H/7` to `6H/49`.  Therefore

```text
H/L>=59625/1485836.                                      (6a)
```

For a protected needle with `L>=1/(7m)`, this gives

```text
mH>=59625/10400852,
min_i(s_i)/m<72805964/59625.                             (6b)
```

On the separated slow-gap branch `a=min(S)>=13m`, where the other six
deleted combs cover a complete `a`-gap, THM-1233 further gives

```text
max(S)/m<34145997116/11925,                              (6c)
```

and THM-1241/1236 force one of five macroscopic cut indices together with the
exact continuous-pivot median invoice in that packet.

More generally, after shrinking the protected interval to the canonical
length `|I|=1/(7m)`, the branch `a=min(S)<13m` makes `(F7)` automatic for
every `C>=195/1078`.  Thus every genuinely nontrivial connected-`G_gt` row
also lies in the separated complete-gap regime.

THM-1237 adds a second exact dichotomy for any covered protected needle.  If
`T0` is its THM-1221 floor tree, then

```text
min(S)<616m/5
or some ij in T0 has gcd(s_i,s_j)<=286m/5.               (6d)
```

Finally, the exact phase quotient is not a bare gcd period.  For every
primitive channel `(p,q)`, a unimodular change of coordinates writes a speed
pair as a skew flow `(kt,ht)`, where `h=qa-pb`.  Its static transverse fiber
has zero minimum exactly on THM-605's nine channels `p+q<=7`.  A relation
cycle on those channels has an exact multiplicative-holonomy identity: off
product-one holonomy it gives an explicit speed crown, while product-one
holonomy forces a signed determinant cocycle.  Thus the connected branch is
no longer merely “unbounded projective height.”  Its exact remaining object
is an addressed local system carrying the full gcd/exact-relation sheet,
primitive low-height detuning, tooth phase, and deletion/blocker circuit.

## 1. Exact charge and oriented forms

Since `s_i=g a` and `s_j=g b`,

```text
1/s_i+1/s_j=(a+b)/(gab).
```

Multiplication by `eta ab/(a+b)` proves (1).  Summing it over a tree proves
(3).  Equivalently, after rooting `T` and charging every edge to its child,

```text
E_T=sum_(i != root)
 [eta_(i,p(i)) * (s_i/gcd(s_i,s_p(i)))]/s_i.             (7)
```

This gives a precise sufficient replacement for the false universal claim.
If a rho-heavy rooted tree has every tail height

```text
s_i/gcd(s_i,s_p(i))<=K,
```

then `rho<=1/14` gives `eta<=13/196`, and hence

```text
E_T<=(13K/196)H.                                        (8)
```

The relevant carrier is therefore not merely a low/high edge color.  It is a
reduced ratio channel together with its projective height and an oriented
vertex-load assignment.

## 2. The translated-LCM counterfamily

Let

```text
A=27720 N,       N>=1,
R={1,2,3,5,7,11,13},
S_A={A+r:r in R}.                                       (9)
```

Every difference of two elements of `R` divides `27720`, while the elements
of `R` are pairwise coprime.  If `r<s`, then

```text
gcd(A+r,A+s)
 =gcd(A+r,s-r)
 =gcd(r,s-r)
 =gcd(r,s)=1.                                           (10)
```

Thus all twenty-one pair periods are the whole circle.  THM-1166's folded
formula and defect bound give, for every edge,

```text
rho(A+r,A+s)
 =1/49+[F(r+s)-F(s-r)]/[196(A+r)(A+s)]
 >x_A:=1/49-1/(4A^2).                                  (11)
```

Here `A=0 mod 14`.  The exact twenty-one folded corrections, in lexicographic
pair order on `R`, are

```text
(20,16,8,0,-16,-24,32,16,0,-32,-20,
 24,0,-48,-16,0,-24,-8,0,0,16).                        (12)
```

Since

```text
1/49-5/308=9/2156,
```

the inequality `9A^2>539` implies `x_A>5/308`.  It already holds at
`A=27720`.  Hence `G_gt=K_7`, and every spanning tree satisfies

```text
sum_(ij in T)rho_ij>6*(5/308)=15/154.                   (13)
```

On the other hand every `g_ij=1`.  Because `rho<1/2`, the map
`x -> x(1-x)` is increasing on the relevant range, so

```text
E_T>6 x_A(1-x_A)>6(5/308)(303/308)=4545/47432.          (14)
```

But

```text
H_A=sum_(r in R)1/(A+r)<7/A.                            (15)
```

Consequently, uniformly over **all** spanning trees,

```text
E_T/H_A >(6A/7)x_A(1-x_A)
         >(4545/332024)A -> infinity.                   (16)
```

Given any proposed finite `C`, choosing
`A>332024 C/4545` contradicts the second inequality in (4).  More precisely,
the fixed correction bank (12) gives the uniform asymptotics

```text
rho_ij=1/49+O(A^-2),
E_T=288/2401+O(A^-2),
H_A=7/A-42/A^2+O(A^-3),
E_T/H_A=(288/16807)A+O(1).                              (17)
```

No choice of maximum tree, Pareto tree, root, or charge orientation repairs
this family: both reduced endpoint heights on every edge are order `A`.

## 3. A protected covered needle where the error loses the answer

The obstruction above is to the unsigned period-error proof, not to HYP-7678's
phase-aware F7 statement.  This can be seen inside the exact protected-needle
geometry.  Put

```text
q=A+1,
S_A=q+{0,1,2,4,6,10,12},
W={ (3q+1)/2+k : 0<=k<=5 },
m=(3q+11)/2,
I=[1/q-1/(14m), 1/q+1/(14m)].                           (18)
```

The number `q` is odd, so the six core speeds are integers, and
`|I|=1/(7m)`.  The core's uniform clearance beyond `1/14` is bounded below by

```text
(5q-77)/(14q)>0                                         (19)
```

for `q>=17`.  For every deleted speed `q+d`, `d<=12`, its distance from the
integer center throughout `I` is at most

```text
12/q+(q+12)/[7(3q+11)].                                 (20)
```

The clearance of `1/14` above (20) is exactly

```text
(q^2-517q-1848)/[14q(3q+11)].                           (21)
```

The numerator is `-288` at `q=520` and `236` at `q=521`, and increases
thereafter.  Our first value is `q=27721`.  Hence `W` is safe throughout `I`
while **all seven** deleted danger combs contain `I`.  Every one of the
twenty-one restricted pair intersections equals `I`, so every local tree has
weight `6|I|`.

Yet the separately summed period-error envelope remains order one and gives
no positive edge credit.  Indeed, when `g=1`,

```text
L rho-rho(1-rho)=rho(L-1+rho).
```

For `L<=139/154` and `rho<=1/14`, this is at most

```text
rho(139/154-13/14)=-2rho/77<0.                          (22)
```

Thus the gcd-period abstraction permits incompatible worst endpoint phases
on six tree edges simultaneously.  The actual wall events in (18) are instead
perfectly aligned.  This proves that the negative result does not refute F7;
it identifies information destroyed before F7 is asked.

## 4. The disconnected strict-spectrum rescue

THM-1221 classifies every disconnected `G_gt` packet into finite projective
branches.  Identity (1) turns those ratio banks directly into positioning
constants.

If the weak-high graph `G_ge={rho>=1/63}` is disconnected, normalization at
its singleton component places all vertices in

```text
V_<={1} union R_<,
```

where `R_<` is the fourteen oriented strict-low ratios.  Exact evaluation of
all pairs in this fifteen-vertex bank gives

```text
max kappa=85975/342804,                                 (23)
```

attained on reduced channel `20:33`.  Since a seven-vertex tree has degree at
most six, (3) gives

```text
E_T<=(85975/57134)H.                                    (24)
```

Now suppose `G_ge` is connected but `G_gt` is disconnected.  In the `1+6`
cut, all vertices lie in `{1}` plus the twenty-four oriented closed-low
ratios.  The complete twenty-five-vertex bank has

```text
max kappa=224458/584325                                 (25)
```

at the vertex pair `(5/9,9/5)`, whose mutual reduced channel is `25:81` and
whose pair mass is `97/4725`.  Therefore

```text
E_T<=(448916/194775)H.                                  (26)
```

THM-1221's `2+5` branch has only twelve normalized packets.  Their complete
pair bank has

```text
max kappa=43774/276507,
E_T<=(87548/92169)H.                                    (27)
```

The `3+4` branch is empty, and both constants in (24),(27) are smaller than
(26).  This proves (5).  Summing THM-1166's interval discrepancy lower bound
on the THM-1221 floor tree proves (6).

If the seven combs cover `I`, THM-1166's adaptive forest-Hunter inequality
supplies the sharper reverse local tree bound `<=6H/49`.  Thus the
disconnected branch has the explicit harmonic crown

```text
H >= [ (15/154)/(448916/194775+6/49) ] L
  = (59625/1485836)L.                                   (28)
```

Since a protected needle has `L>=1/(7m)` and
`H<7/min_i(s_i)`, (28) also gives

```text
mH>=59625/10400852,
min_i(s_i)/m<72805964/59625.                            (28a)
```

If `a=min_i(s_i)>=13m` and the six faster combs cover the complete `a`-gap,
THM-1233's `max_i(s_i)/a<2345` composes with (28a) to give

```text
max_i(s_i)/m
 <2345*(72805964/59625)
 =34145997116/11925.                                   (28aa)
```

THM-1241/1236 additionally force `max(S)/min(S)>70/69`, some adjacent ratio
`>211/210`, and the continuous-pivot median invoice.  Thus the disconnected separated branch is uniformly bounded
relative to the retained core and has only five macroscopic cut positions.

The earlier fragmentation cap `H/7` gives the weaker diagnostic ratio
`417375/10488302`.  The connected strict-high family (9) shows exactly why
neither crown extends through the last component branch using gcd periods
alone.

### Complete strict-high connectivity does not force a heavy circuit

There is a second translated-LCM packet that separates the THM-1221 graph
from THM-1218's deletion-circuit predicate.  Put

```text
A=55440N,
R'={1,2,3,5,11,13,17},
S'_A={A+r:r in R'}.                                     (28b)
```

The offsets are pairwise coprime, and their difference set is

```text
{1,2,3,4,6,8,9,10,11,12,14,15,16},
```

whose every element divides `55440`.  The gcd argument in (10) again makes
all twenty-one translated pairs coprime.  The offsets contain no four-term
arithmetic progression.  Finally THM-1166's coarse floor gives

```text
rho(A+r,A+s)>=1/49-1/[4(A+r)(A+s)]>5/308,              (28c)
```

so `G_gt=K_7` and every tree clears `15/154`.  Nevertheless every quartet is
non-AP.  The contrapositive of THM-1218 therefore bounds every quartet BAD
mass by `60/637`: the heavy-circuit graph is empty.

Thus complete strict-high connectivity cannot manufacture the deletion
circuit needed to spend the positive global incidence margin.  Their
correlation must come from wall phase, endpoint ownership, or a higher
incidence carrier.

## 5. The nontrivial connected branch is a compact phase kernel

Take the protected interval at its canonical sublength

```text
|I|=1/(7m),                   a=min_i(s_i).              (28d)
```

If `a<13m`, then

```text
H>=1/a>1/(13m)=7|I|/13.                                 (28e)
```

For `eta=15/154`, any `C>=13eta/7=195/1078` makes
`eta|I|-CH<0`.  Since restricted pair intersections are nonnegative, every
tree satisfies `(F7)` trivially.  In particular the disconnected constant
`448916/194775` is far above this threshold.

It remains only `a>=13m`.  An interval of length `1/(7m)` then contains a
complete safe gap of the `a`-comb: the worst start-to-full-gap span is
`13/(7a)<=1/(7m)`.  The other six deleted combs must cover that gap.  THM-1233,
THM-1241, and THM-1236 therefore apply.  Writing their ordered speeds as
`d_1<...<d_6`, they give

```text
1<d_i/a<2345,
y0=max(7a/6,d_4),
sum_i |d_i-y|>y/14 for every y>=7a/6,
some q=d_(j+1)-d_j>d_6/210>a/210.                       (28f)
```

with one of five cut indices (and the sharper `a/180`, `a/132`, `a/108`
cut floors as one, two, or at least three pivots cross `7a/6`).

THM-1237's positioned protected-needle debt composes with this reduction.
For its THM-1221 floor tree `T0`, coverage forces

```text
24mH+13m sum_(ij in T0)1/g_ij >=30/11.                  (28g)
```

If `a>=616m/5`, then `H<7/a<=5/(88m)`, so the harmonic term in
(28g) is at most `15/11`.  The six reciprocal gcd edges must pay the other
half:

```text
m sum_(ij in T0)1/g_ij >=15/143.
```

Consequently some floor-tree edge has

```text
g_ij/m<=286/5.                                         (28h)
```

This proves (6d).  The connected row is therefore either absolutely bounded
at its least wall by `616m/5`, or it exposes a mixed clock of size at most
`286m/5`; merely saying “the ratios are compact” misses this second scale.

Thus THM-1226's unbounded fully coalescing counterfamily cannot lie in the
remaining large-scale branch `min(S)>=616m/5`.  Its covered embedding from
Section 3 is shallow (`min(S)/m` tends to `2/3`) and is already F7-vacuous.
The entire unresolved connected branch reduces to a compact six-ratio packet,
one of five macroscopic cuts, and the arithmetic cyclic phase/offset stalk
inside the two at-most-five-speed subclusters.  A
macroscopic difference alone is not a local overlap estimate; the next
section gives an exact counterexample and identifies the missing coordinates.

## 6. The exact phase fiber, its height-seven threshold, and its dynamic limit

Let

```text
B=(-1/14,1/14) in R/Z,
pi_(p,q)(x,y)=qx-py,
```

where `1<=p<=q` and `gcd(p,q)=1`.  Choose integers `u,v` with
`qu-pv=1`.  The map

```text
(tau,z) -> (p tau+u z, q tau+v z)                       (29)
```

is a Haar-preserving torus automorphism and sends `pi_(p,q)` to `z`.
Therefore the pushforward of `1_(B x B) dx dy` has the exact fiber density

```text
F_(p,q)(z)=sum_(n in Z) f_(p,q)(z+n),                   (30)
```

where

```text
f_(p,q)(w)=
  1/(7q),                              |w|<=(q-p)/14,
  ((p+q)/14-|w|)/(pq),  (q-p)/14<=|w|<=(p+q)/14,
  0,                                  otherwise.        (31)
```

Indeed, `qx` and `-py` push the two real radius-`1/14` intervals to
rectangles of widths `q/7` and `p/7`; their convolution is the trapezoid
(31), and reduction modulo one gives (30).  Its plateau area is
`(q-p)/(49q)` and its two shoulders have total area `p/(49q)`, so

```text
integral F_(p,q)=1/49.                                  (32)
```

If `p+q<=7`, the support half-width `(p+q)/14` leaves a zero phase (only the
antipodal seam when equality holds).  If `p+q>7`, a nearest lift of every
`z` has absolute value at most `1/2` and lies strictly inside the support;
more precisely

```text
F_(p,q)(z)>=min{1/(7q),(p+q-7)/(14pq)}>0.              (33)
```

Thus

```text
min_z F_(p,q)(z)=0 iff p+q<=7.                          (34)
```

This recovers, and (30)--(33) explicitly refine, the already-proved static
THM-605 dangerous-pattern theorem.  Its primitive phase channels are exactly

```text
(1,1), (1,2), (1,3), (1,4), (2,3),
(1,5), (1,6), (2,5), (3,4).                            (35)
```

`CombPatterns.lean` already kernel-checks the avoidance construction and the
converse interior-point existence underlying (34); the new referee replays
THM-605's exact minima table from the single periodized trapezoid.

The useful dynamic normal form is equally exact.  For speeds `a,b`, put

```text
h=qa-pb,                 k=ub-va.                       (36)
```

Then

```text
(a,b)=(pk+uh,qk+vh),     gcd(k,h)=gcd(a,b),             (37)
```

and the concrete pair orbit is the skew flow `(tau,z)=(kt,ht)` through the
fiber (29).  In particular

```text
D_a intersect D_b subset {t: ||ht||<(p+q)/14}.          (38)
```

This also shows why relation coefficients must be primitive.  The old
THM-864 localization theorem allowed a multiplied presentation and thereby
bought false `1/y` decay.  At radius `1/13`, its exact row

```text
A=3744, B=3745, p=q=y=12, E=[1/3,1/2]
```

has error `2/507`, while the claimed clean bound is only
`13129/3359265`, smaller by `531/14556815`.  The primitive relation is just
`B-A=1`; twelve copies create no new orbit samples.  MISTAKE-184 records the
correction.

Even primitive height-seven data are not by themselves a finite-interval
inverse.  For `N>=14` and `N=0 (mod 14)`, the pair

```text
a=N+1,       b=2N+1,       2a-b=1
```

has

```text
rho(a,b)=1/49+6/(49ab)>1/63,
D_a intersect D_b=empty on [3/14,11/14].               (39)
```

So a macroscopic ratio cut can still have a macroscopic empty overlap window.
Conversely, the high-height exact relation is also indispensable.  For every
`K>=1`, the pair `(64K,75K)` has an empty-overlap interval

```text
I_K=[407/(896K),407/(896K)+449/(4928K)].                (40)
```

Every nonzero coefficient vector of `l1` height at most seven has phase drift
at least

```text
11K |I_K|=449/448>1,                                   (41)
```

yet the omitted exact relation `(75,-64)` of height `139` gives common period
`1/K`, longer than the interval.  This refutes THM-598's claimed resolved
local independence and THM-602's truncated-lattice fully-resolved branch;
see MISTAKE-185.  It does not affect THM-605's static theorem.

The lesson is exact: the local carrier needs both the full primitive exact
kernel/gcd period with its torsion sheet and the low-height detuned phase
channels.  Neither projection reconstructs the other.

## 7. Primitive relation trees and the multiplicative-holonomy circuit

Suppose a rooted relation tree has primitive edge labels `(p_e,q_e,h_e)`
from (35), with `|h_e|<=Bm`.  Along a path of length at most six, repeated
use of (37) gives

```text
Q_v s_v-P_v s_root=H_v,
P_v,Q_v<=6^6=46656,
|H_v|<=6*6^5 Bm=46656Bm.                               (42)
```

Hence the leading rational tree shape belongs to a finite bank and its
tangent offsets are bounded.  If every `h_e=0`, each edge lies in the
nine-channel bank (35).  Directly evaluating the oriented charge
`eta(rho(p,q))*max(p,q)` gives

```text
(1,1):6/49,       (1,2):13/98,
(1,3),(2,3):20/147,
(1,4),(3,4):27/196,
(1,5),(2,5):34/245,
(1,6):41/294.                                          (43)
```

Thus the exact height-seven constant is `C=41/294`, attained at `(1,6)`.
This sharpens the valid but coarse substitution `K=6` in (8), which gives
`39/98`.

Cycles add the compatibility that a tree lacks.  Orient a relation cycle of
length `r` by

```text
q_i s_(i+1)-p_i s_i=h_i,       s_(r+1)=s_1.            (44)
```

Iterating once around it gives the exact holonomy identity

```text
(prod_i q_i-prod_i p_i)s_1
 =sum_i h_i (prod_(j<i)q_j)(prod_(j>i)p_j).             (45)
```

If `|h_i|<=Delta`, `p_i,q_i<=Q`, and the product holonomy is nontrivial,
integrality yields

```text
s_1<=r Delta Q^(r-1).                                  (46)
```

For `r<=7`, this is `s_1<=326592 Delta` on the nine height-seven channels,
or `s_1<=33787663 Delta` if coefficients through thirteen are retained.  If
the product holonomy is one, (45)'s right side must vanish.  All its weights
are positive, so a nonzero balanced circuit must contain both positive and
negative `h_i`; a monochromatic nonzero relation cycle is automatically in
the bounded branch.

THM-1240 now proves that every six-cover supplies a blocker cycle, and its
centered-spoke data give a natural exact instance of (45).  Fix the canonical
circle lift of the carrier gap, with gap index `0<=k<c`.  Then its selected
phase lies in `0<t_i<1`; writing `Q_i=c+d_i` and `t_i=P_i/Q_i` gives
`0<P_i<Q_i`.  If `j=b(i)` is a chosen blocker, choose the nearest tooth
address `N_i` and signed residue `r_i`.  Since `d_j>c` and the carrier gap
stays more than `1/(14c)` from either endpoint of the unit interval, danger
forces `0<N_i<d_j`, and

```text
P_i d_j-N_i Q_i=r_i,          |r_i|<Q_i/14.
```

Then

```text
P_i Q_j-N_i Q_i=P_i c+r_i=:K_i>0.                      (46a)
```

Positivity follows because `t_i` lies strictly inside the `c`-safe gap, so
`cP_i/Q_i>1/14`, while `r_i/Q_i>-1/14`.  Applying (45) around every directed
blocker cycle proves the new exact, **canonical-gauge** sign law

```text
prod_cycle P_i > prod_cycle N_i.                       (46b)
```

Thus the proved centered cycle is automatically in the positive,
non-balanced holonomy branch in this fixed lift; a product-one signed cocycle
cannot hide it.  The product comparison is not invariant under arbitrary
integer translations of the addresses, which is why the canonical lift is
part of the statement.  What remains is quantitative **address compression**:
the present argument supplies no uniform bound on the reduced sizes of
`P_i,N_i`, so (46) is not yet a uniform crown.
THM-1239 proves that one resonant blocker can erase every crack on an
arbitrarily chosen gap, while THM-1242's `q=15` sunflower proves that
nontrivial quotient periods alone need not leave a beat residue.  A closing
theorem must therefore turn (46b), the separated cut-clock doublet, and the
positioned mass debt into an alternate-gap/address descent.  No scalar
connectivity or low-height inverse can substitute for that global transport.

## 8. Tournament and alternate-vertex audit

For (9), switching at `1/63` marks all twenty-one edges strict-high.  The
speed-order gauge is transitive: score histogram `(0,1,2,3,4,5,6)`, no
directed cycles, seven singleton SCCs, one Hamiltonian path, and zero low-edge
flips.  These fingerprints are constant while `E_T/H` diverges.  Thus the
runner tournament preserves global connectivity and the `15/154` Hunter
credit but destroys the projective height that controls localization.

Using gcd periods as vertices alone is worse on this family: every edge
becomes the same period-one vertex.  That quotient destroys the common wall
alignment in (18).  The opposite truncation to the nine THM-605 channels
destroys the high exact period in (40).

There is a second, non-runner tournament on the nine primitive phase channels.
Use zero-fiber slack `(7-p-q)_+` as the pairwise obstruction observable and
Farey height followed by lexicographic order as the tie gauge.  Its tie
Hamiltonian path is exactly the order in (35); the fingerprint is transitive,
with score histogram `(0,1,2,3,4,5,6,7,8)`, no directed cycles, nine singleton
SCCs, and one Hamiltonian path.  This quotient remembers which static phase
can vanish but still loses `k,h`, the exact gcd sheet, and interval address.

The meaningful cyclic tournament uses **proof obligations** as vertices.
For a centered carrier spoke or active pair, orient `i -> j` when wall `j`
is chosen as its third blocker, and attach

```text
edge sidecar=(primitive p,q; exact k,h; gcd sheet; tooth address;
              signed endpoint cocycle; deletion owner).                    (47)
```

Every finite blocker functional graph has a directed cycle; (45) is the exact
arithmetic readout of such a cycle once its relation labels are supplied.
THM-773 explains why owner sheet, inverse step, endpoint order, and global
carry must remain in its stalk.  The challenged assumptions are now both
directions of the same forgetting error: a global heavy runner tree cannot be
localized by absolute gcd errors, and a static low-height phase atlas cannot
discard the full exact orbit.  The smallest plausible quotient is their
fibered product with the blocker/deletion circuit.

## 9. Formalization and verification boundary

The dependency-free referee checks (1) in `2,691` integer pairs; freezes all
twenty-one corrections in (12); verifies four exact members of the infinite
family; proves the spectral, error, harmonic, and asymptotic-facing
inequalities; checks the AP-free separation packet (28b) at three scales;
checks every core and deleted endpoint in (18); recomputes the strict,
closed, and twelve-packet projective maxima in (23)--(27); verifies the exact
fiber threshold on `278` primitive pairs and replays THM-605's published
minimum bank; checks `4,725` unimodular relation coordinates; proves the two
localization counterexamples exactly; and replays `2,406` relation cycles.

`LRCGCDPeriodProjectiveCharge.lean` kernel-checks the symmetric and oriented
charge identities, the seven-vertex load consumer, all displayed transfer
constants, the localized cleared-denominator consumer, the sharper forest
crown composition, the shallow-F7 vacuity, the nine-channel enumeration,
the sharp `41/294` oriented-charge table, primitive relation coordinates and
their complete common-divisor/gcd invariance, THM-1237 debt composition and
six-edge pigeonhole, the triangle holonomy/sign consumers, the centered
blocker lift, the strict-high quadratic threshold, and both protected-needle
margin polynomials.  It imports
`CombPatterns.lean`, which already formalizes the static avoidance and
interior-point arithmetic behind THM-605.  The measure-theoretic positivity,
periodized-trapezoid disintegration, arbitrary-speed Haar
formula, THM-1221 branch classification, and finite-bank identification
remain explicit paper/exact-referee providers rather than being overstated as
formalized.

Frozen hashes:

```text
source         9ae71334d58a3d5613cb4197724e4af5d29c95f2e0d87e5a8ce6f9b0341f55a7
output         cfc88e17e0177a4288f3521ff3437a063a9ba68f6646d251adf41c3b39609503
formalization  f23ffb3b6e331f73e1be5549d092431d6e4fa93f7a0b3c7ea9c0d59b007031aa
```
