---
id: THM-2247
title: "Nonsplit terminal quartic degree-fourteen closure"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED. In the polynomial exact-square
  prefix, nonsplit terminal quartic branch, no Keller mate has reduced fibre
  degree fourteen. THM-2245's singular spectral quartic normalizes first to
  a conic and then to a residual T-double-cover. On the smooth conic, the
  residual radicand is u^-3 times an explicit sextic. Its unavoidable odd
  valuations at zero and infinity force that sextic to be a scalar square
  if the cover is rational, but two exact coefficient residuals have
  nonzero algebraic norms and rule out every such square. On the reducible
  conic, each component gives a cubic double-cover. It has genus one unless
  the spectral quartic has a quadruple root; in that cusp the retained
  Keller one-form is a seventh-power rational primitive, and polynomiality
  of the quadratic approximate root gives an unavoidable pole. Together
  with THM-2214, a noninvertible survivor in this precise nonsplit terminal
  branch has reduced degree at least eighteen. Split/even descent, other
  short edges, higher degrees, JC(2), and DC(2) remain open.
source: codex-2026-07-25-degree14-singular-cover
depends_on:
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2217-square-prefix-pole-alternative-and-odd-leading-degree-terminal-wall
  - THM-2245-degree-fourteen-spectral-quartic-discriminant-reduction
related:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
script: 04-computation/jc2_nonsplit_terminal_degree14_closure_thm2247.py
output: 05-knowledge/results/jc2_nonsplit_terminal_degree14_closure_thm2247.out
script_sha256: a32904378fc37e6410f6a1f11bd77d43f2d7150cf9143eaa8ff961e4bb661406
output_sha256: eccf74910cc562cc888997e3576b4a117b881145f321775fc5912b3148dbdf4c
hash_basis: working-tree bytes (LF)
---

# THM-2247 -- the nonsplit terminal degree-fourteen branch is empty

THM-2245 reduces the first open nonsplit terminal degree to one singular
spectral quartic, but deliberately stops before analyzing the second double
cover. This theorem computes that cover and closes every one of its
degeneracies. It concerns only the exact-square-prefix terminal branch of
THM-2214/2245. It is not a proof of planar JC.

## 1. Inherited carrier

Work over

```text
R=C[x],                    K=C(x),
```

with `C` algebraically closed of characteristic zero. Use THM-2245's
nonsplit coordinates

```text
P=w^4+2d w^2+q w+(d^2-s),             q^2=T,
Q=E_14+alpha E_10+beta E_6+gamma E_2,              (1)
```

after the exact target-shear quotient and full Faber gauge. The deck sends
`(w,q)` to `(-w,-q)`, fixes `d,s,T`, and is genuine. In particular,

```text
q!=0,                 T!=0.                          (2)
```

Write `Psi_0 in C` for the second constant flux and set

```text
y=7s-2alpha.                                         (3)
```

THM-2245 proves `y!=0`, eliminates the first two fluxes to

```text
(235298T^2+B(y))^2=(10976y)^2 H(y),                 (4)
```

and shows that every survivor has a multiple root in the depressed quartic
`H`. Thus, for some `e,D in C`,

```text
H(y)=425(y-e)^2((y+e)^2+D).                          (5)
```

Choose `c in C` with

```text
c^2=425.                                             (6)
```

The sign of `c` is immaterial and will be audited below.

Matching all four coefficients in (5) gives

```text
beta
 =5(17D+12alpha^2-34e^2)/168,

gamma
 =5(68Dalpha+17De+16alpha^3-136alpha e^2)/1568,

Psi_0
 =-5(289D^2+136Dalpha^2+68Dalpha e-1224De^2
      +16alpha^4-272alpha^2e^2+1088e^4)/10976.       (7)
```

The large cubic in (4) then collapses to

```text
B(y)=13720(16y^3-17De).                              (8)
```

This disappearance of `alpha` is the first structural simplification: the
residual double-cover depends only on the singularity coordinates `(e,D)`,
not on the chosen center of the Faber train.

The trajectory cannot be identically at the singular point `y=e`. If it
were, then `H(y)=0`, hence the square root `W` in THM-2245 would vanish.
Equation (4) would make `T^2` constant. The constant-field argument would
make `T`, and then the anti-invariant `q`, constant, contradicting the
genuine nonsplit deck. We may therefore divide by `y-e` as a rational
function.

Define

```text
v=W/(c(y-e)).
```

Equations (4)--(6) give the first normalization

```text
v^2=(y+e)^2+D,                                      (9)

235298T^2=10976c y(y-e)v-B(y).                     (10)
```

The rest of the proof separates the smooth conic `D!=0` from the reducible
conic `D=0`.

## 2. Smooth conic: the residual sextic cannot be a square

Assume `D!=0`. Parameterize (9) by the nonzero rational function

```text
u=(y+e)+v.
```

Because

```text
((y+e)+v)(v-(y+e))=D,
```

one has

```text
y=(u-D/u)/2-e,              v=(u+D/u)/2.             (11)
```

Substituting (8) and (11) into (10) gives

```text
T^2=2 S_6(u)/(343u^3),                               (12)
```

where

```text
S_6(u)
 =(c-20)u^6
  -6e(c-20)u^5
  +(-Dc+60D+8ce^2-240e^2)u^4
  +10e(-7D+16e^2)u^3
  +D(-Dc-60D+8ce^2+240e^2)u^2
  +6D^2e(c+20)u
  +D^3(c+20).                                       (13)
```

The companion obtains (13) by direct substitution into (10). It records
the unsimplified but identical scalar

```text
1372/235298=2/343.
```

Since

```text
(c-20)(c+20)=25,                                   (14)
```

both factors `c-20` and `c+20` are nonzero. Thus `S_6` has degree exactly
six, has nonzero constant coefficient, and the radicand in (12) has odd
valuation `-3` at both `u=0` and `u=infinity`.

### 2.1 Branch count forces a scalar square

For a rational function `R(u)`, the normalization of

```text
Z^2=R(u)
```

is branched exactly at the points where `ord(R)` is odd. If the branch
support has `r` points, its genus is

```text
g=(r-2)/2.                                           (15)
```

The pair `(u,T) in C(x)^2` in (12) gives a rational map from `P^1` to this
normalization. If `g>=1`, properness extends the map across its poles and
Riemann--Hurwitz makes it constant. Then `u,T`, hence `y,v`, are constant.
THM-2245's Keller one-form

```text
A(2T F'+F T')=2kappa T,              kappa!=0       (16)
```

would have zero left side and nonzero right side, by (2). Therefore the
cover in (12) must have genus zero.

The points zero and infinity already contribute two branches. Hence every
finite nonzero root of `S_6` must have even multiplicity. Equivalently,
over the algebraically closed constant field,

```text
S_6=(c-20)(u^3+r u^2+b u+k)^2.                       (17)
```

No genericity assumption is being made: (17) is forced precisely in the
multiple-root cases where a squarefree-sextic argument would be weakest.

### 2.2 Two coefficient residuals rule out every square

Matching the coefficients of `u^5,u^4,u^3` in (17) uniquely gives

```text
r=-3e,

b=(8Dc+155D-16ce^2-325e^2)/10,

k=e(10Dc+185D-16ce^2-335e^2)/10.                   (18)
```

The residual coefficients of `u^2` and `u` in
the right side of (17) minus (13) are

```text
R_2
 =[69D^2c+1420D^2-250Dce^2-7360De^2
   +293ce^4+5500e^4]/4,                             (19)

R_1
 =e[71D^2c+1130D^2-282Dce^2-5730De^2
    +259ce^4+5380e^4]/2.                            (20)
```

If `e=0`, then

```text
R_2=D^2(69c+1420)/4.
```

The algebraic norm of the remaining factor is

```text
1420^2-425*69^2=-7025!=0,                           (21)
```

so `R_2` cannot vanish.

If `e!=0`, put `z=D/e^2`. Dividing (19) by `e^4/4` and (20) by `e^5/2`
gives the two quadratics

```text
p_2(z)
 =(69c+1420)z^2+(-250c-7360)z+(293c+5500),

p_1(z)
 =(71c+1130)z^2+(-282c-5730)z+(259c+5380).          (22)
```

Their resultant, reduced with `c^2=425`, is

```text
Res_z(p_2,p_1)=-8670000(110c+2281).                  (23)
```

Again the final factor is nonzero, because its norm is

```text
2281^2-425*110^2=60461!=0.                          (24)
```

Thus (19) and (20) never vanish simultaneously. The forced square (17) is
impossible, closing every `D!=0` trajectory.

## 3. Reducible conic: genus one except at the quadruple cusp

Now assume `D=0`. Equation (9) factors in the function field, so for one
fixed sign `sigma in {+1,-1}`,

```text
v=sigma(y+e).                                       (25)
```

Substitution into (10) gives the exact cubic cover

```text
T^2
 =(16/343)y[(sigma c-20)y^2-sigma c e^2].           (26)
```

Both signs in (25) are retained; no component of the reducible conic has
been discarded.

If `e!=0`, the cubic in (26) has the three distinct roots

```text
0,      +/- sqrt(sigma c e^2/(sigma c-20)).
```

Indeed `sigma c`, `sigma c-20`, and `e` are nonzero. Equivalently, its
discriminant is

```text
4(sigma c-20)(sigma c)^3 e^6!=0.                    (27)
```

The double cover is a smooth genus-one curve. The same
Riemann--Hurwitz/one-form argument following (15) makes the trajectory
constant and contradicts (16). Hence

```text
e=0.                                                 (28)
```

This is exactly the maximally singular spectral quartic

```text
H(y)=425y^4.
```

It remains to show that its rational cusp cannot satisfy the Keller
one-form and polynomiality simultaneously.

## 4. Quadruple cusp: the Keller primitive creates a forbidden pole

Fix the component sign and abbreviate

```text
r=sigma c.
```

At `D=e=0`, coefficient matching (7) gives

```text
beta=5alpha^2/14,             gamma=5alpha^3/98.     (29)
```

Equation (26) becomes

```text
T^2=lambda^2 y^3,
lambda^2=16(r-20)/343!=0.                            (30)
```

Because `y,T` are nonzero rational functions, define

```text
rho=T/(lambda y).
```

Then

```text
y=rho^2,                       T=lambda rho^3.       (31)
```

### 4.1 The third flux is a seventh-power primitive

THM-2245's exact third observable is

```text
F(y,T)
 =-T[-343T^2+640alpha^3-320alpha^2y
      -2688alpha beta+896beta y+6272gamma+160y^3]
   /(8960y).                                         (32)
```

Substitute (29)--(31). All terms involving `alpha` cancel and

```text
F=(r-30)T y^2/560
  =(r-30)lambda rho^7/560.                           (33)
```

Divide the Keller one-form (16) by the nonzero `T`. Since
`T'/T=3rho'/rho`, equation (33) gives

```text
A[2F'+F T'/T]
 =A * [17(r-30)lambda/560] rho^6 rho'
 =2kappa.
```

Therefore

```text
A(rho^7)'=7840kappa/[17(r-30)lambda] in C*.          (34)
```

The denominator is nonzero because `r^2=425!=900`.

Apply the rational-primitive lemma used in THM-2214. If

```text
A in C[x]\{0},        S in C(x),        AS' in C*,
```

then either `A` is constant and `S` is affine-linear, or, after translating
`X=x-xi`,

```text
A=a_0X^m,          m>=2,
S=s_0+s_1X^(1-m),  s_1!=0.                           (35)
```

The constant case is impossible for `S=rho^7`: absence of finite poles
makes `rho` a polynomial, while

```text
deg(rho^7)=7 deg(rho)!=1.
```

In the second case, a binomial with two nonzero terms has simple nonzero
roots and cannot be a seventh power. Hence `s_0=0`, `m-1=7k` for some
`k>=1`, and

```text
rho=Lambda X^(-k),              Lambda!=0.           (36)
```

Thus the Keller one-form forces a pole of `rho`.

### 4.2 Polynomiality of the approximate root forbids that pole

Return to the polynomial exact-square prefix of THM-2214:

```text
L=Az+E in C[x,z],

mathcal H=a(L+s)^2+d in C[x,z],
a=1/T.                                               (37)
```

The first flux formula of THM-2245, under (29)--(31), becomes

```text
dT=(2r-35)y^2/245.                                  (38)
```

Consequently

```text
a=1/(lambda rho^3),
s=(rho^2+2alpha)/7,
d=(2r-35)rho/(245lambda).                           (39)
```

At the pole `X=0` in (36), the terms involving the polynomial `E` in the
constant coefficient of (37) are regular:

```text
aE^2=O(rho^-3),                 2aEs=O(rho^-1).
```

The remaining part is exact:

```text
a s^2+d
 =[2(r-15)/(245lambda)]rho
   +[4alpha/(49lambda)]rho^-1
   +[4alpha^2/(49lambda)]rho^-3.                    (40)
```

The coefficient of the unique polar term is nonzero, since

```text
15^2-425=-200!=0.                                   (41)
```

Nothing regular can cancel it. Thus the constant coefficient of
`mathcal H` has a pole, contradicting `mathcal H in C[x,z]`. This closes
the quadruple cusp and hence the whole `D=0` branch.

## 5. Conclusion and exact frontier

Sections 2--4 exhaust the factorization (5). Therefore the assumptions

```text
polynomial exact-square-prefix quartic
+ genuine nonsplit deck
+ terminal full-Faber normal form
+ reduced mate degree 14
+ nonzero constant Jacobian
```

are contradictory.

THM-2214 already makes degree two tame and excludes reduced degrees six and
ten. The nonsplit parity train has degrees `4R-2`, so the next unclosed
reduced degree in this precise branch is

```text
18.                                                   (42)
```

This does not address the split deck, an even-leading terminal descent,
other short Newton edges, or degree at least eighteen. It proves neither
`JC(2)` nor `DC(2)`.

The hostile boundaries were used rather than suppressed:

* both choices `c=+/-sqrt(425)` obey the nonzero norm checks;
* both reducible-conic components `sigma=+/-1` are replayed;
* `D!=0`, `D=0,e!=0`, and `D=e=0` are separate arguments;
* the smooth-conic proof handles every repeated-root sextic by the forced
  scalar-square comparison; and
* the spectral discriminant is not treated as a continuation state without
  the third flux and Keller one-form.

Run

```bash
python3 04-computation/jc2_nonsplit_terminal_degree14_closure_thm2247.py
python3 -O 04-computation/jc2_nonsplit_terminal_degree14_closure_thm2247.py
```

The exact companion verifies coefficient matching (7)--(8), the conic
substitution (11)--(13), the forced-square residuals and resultant
(18)--(24), both component cubics (26)--(27), the cusp third flux and first
flux, the seventh-power differential, the pole coefficient (40), and all
four `(sign(c),sigma)` hostile replays. Both modes reproduce the frozen
transcript byte for byte. The branch-count, Riemann--Hurwitz,
rational-primitive, and polynomial-pole arguments are the mathematical
proof. QED.
