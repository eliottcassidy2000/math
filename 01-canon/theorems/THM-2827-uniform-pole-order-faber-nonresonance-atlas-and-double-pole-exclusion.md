---
id: THM-2827
title: "Uniform pole-order Faber nonresonance atlas and double-pole exclusion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  normalized reduced Faber degree m=4R-2, R>=7,
  and every local pair (ord V,ord M)=(nu,1), nu>=3, the nonsplit
  polynomial exact-square-prefix chart is impossible when R is 0 or 1
  mod 3.  For R=3k+2 it is impossible off the exact H-resonance rays
  nu=2+(4k+3)delta, a=1+(2k+1)delta.  No existence is asserted on a
  resonance ray.  Applied to THM-2796, every surviving pole part must be
  divisible by 4k+3; in particular double pole parts are excluded in all
  degrees.  MISTAKE-317 records why unique-face support alone did not
  close the resonance.
source: root/uniform-pole-order-faber-nonresonance-atlas-2026-07-28
audit: >
  thm2827-hostile-audit independently rederived the target normalization,
  H extraction, complete valuation trichotomy, exact resonance equivalence
  and local witness, passport divisibility, and primitive Farey boundary;
  normal, optimized, stored, LF hashes, and documentation all pass.
depends_on:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2760-exact-prefix-even-faber-flux-gcd-and-smooth-boundary-exclusion
  - THM-2784-nonsplit-response-square-potential-divisor-and-infinity-classification
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
related:
  - THM-2826-uniform-triple-pole-faber-obstruction-and-simple-pole-chamber-exclusion
  - MISTAKE-317
script: 04-computation/jc_uniform_pole_order_faber_nonresonance_thm2827.py
output: 05-knowledge/results/jc_uniform_pole_order_faber_nonresonance_thm2827.out
script_sha256: 3020a0840ee8f9780a4e21b8cef6a1a06306400ec277cc059a2c23146247d410
output_sha256: 56c94e52ca3acb7b47edb331c569c150f95c26d744be2948c70813573218ad81
hash_basis: LF-normalized bytes
---

# THM-2827 -- the all-order Faber obstruction has one exact resonance ray

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2826 excludes the triple-pole pair `(ord V,ord M)=(3,1)`.  The same
Newton architecture extends to arbitrary `ord V`, but not as an
unconditional closure.  A unique `H` face can have exactly the valuation
required by the response equation.  This theorem gives the complete
nonresonance atlas and isolates that equality ray.

## 1. Uniform Faber formulas

Fix `R>=7`, put

```text
m=4R-2,                         alpha=R-1/2,           (1)
```

and normalize a full reduced representative by the legal constant target
translation:

```text
Q=E_(4R-2)+sum_(j=1)^(R-2)a_jE_(4j-2).               (2)
```

Indeed, before normalization the coefficient of `E_(4R-6)` changes by
`-(R-1/2)c`, so it can always be killed.  The missing `R-1` row makes all
top/lower inequalities below strict.

In the nonsplit polynomial exact-square-prefix chart write

```text
H=Vz^2+B_src z+C_0,              L=A_src z+E_0,
P=H^2+L,                         U^2=V,                (3)

q=A_src/U,                       T=q^2,
d=C_0-B_src^2/(4V),
s=A_src B_src/(2V)-E_0.                                  (4)
```

For the row indexed by `j`, define

```text
A(t)=1+2dt^2+qt^3+(d^2-s)t^4
    =(1+dt^2)^2+t^3(q-st),

A(t)^(j-1/2)=sum_n c_n^(j)t^n,                         (5)

phi_j=4c_(4j-1)^(j)/q,             psi_j=4c_(4j)^(j),
K_j=(4c_(4j+1)^(j)+2d c_(4j-1)^(j))/q.               (6)
```

The inherited flux equations are

```text
phi_Q=0,                          psi_Q in C,
K_Q=R_Q/q,                        A_src K_Q=lambda M,  (7)
```

with `lambda!=0`.  Exact extraction gives

```text
K_j=-(d/2)phi_j+T H_j,                                (8)

H_j =
 4 sum_(h,u,ell>=0; ell+2u=j-2-3h)
   binom(j-1/2,j+1+h) binom(j+1+h,ell)(-1)^ell
   binom(-2-2h,u) T^(2h)(Td)^u s^ell.                 (9)
```

For completeness, `(9)` follows by expanding

```text
(1+dt^2)A(t)^(j-1/2)
 =sum_(v>=0)binom(j-1/2,v)t^(3v)(q-st)^v
             (1+dt^2)^(2j-2v).                       (10)
```

Terms `v<=j` have degree at most `4j`.  Writing `v=j+1+h`, the
`t^(4j+1)` condition is `ell+2u=j-2-3h`, and division by `q^3`
turns the remaining `q` power into `T^(u+2h)`.  Thus

```text
H_j in Q[T,s,Td]                                     (11)
```

for every row.  On `phi_Q=0`,

```text
K_Q=T H_Q,                       H_Q=H_R+sum_(j<=R-2)a_jH_j. (12)
```

## 2. Arbitrary pole order and the exhaustive trichotomy

Let `beta` be a finite response point and set

```text
nu=ord_beta(V)>=3,                 ord_beta(M)=1,
a=ord_beta(A_src),                 b=ord_beta(B_src),  (13)
```

with `b=infinity` if `B_src=0`.  Polynomiality and `(7)` give
`a,b>=0`, with `a` finite.  For `v=ord_beta`,

```text
v(T)=2a-nu,                        v(K_Q)=1-a,
v(d)>=min(0,2b-nu),               v(s)>=min(0,a+b-nu). (14)
```

A negative minimum is exact.  The nonnegative lattice `(a,b)` splits into
three disjoint regions:

```text
polar:       2b<nu and a+b<nu;
pure-q:      not polar and 2a<nu;
regular-H:   not polar and 2a>=nu.                    (15)
```

If the second region is not polar, necessarily `2b>=nu`: otherwise
`a<nu/2` and `b<nu/2` would imply `a+b<nu`.  Hence `(15)` is exhaustive.

## 3. Polar cone: both ends of the prefix polynomial

In the polar cone pass to the local square-root extension and put

```text
omega=B_src/(2U),                  rho=q/omega^3.       (16)
```

The initial source is exactly

```text
in(d)=-omega^2,                    in(s)=q omega.       (17)
```

Write

```text
x=v(s)=a+b-nu<0,                   y=v(rho)=a-3b+nu.   (18)
```

THM-2760 supplies exact-prefix polynomials

```text
in(phi_j)=4 in(s)^(j-1)P_j(rho),
in(psi_j)=4 in(s)^jQ_j(rho),                          (19)

deg P_j=floor((j-1)/3),            deg Q_j=floor(j/3), (20)
P_j(0)Q_j(0)!=0,                   gcd(P_j,Q_j)=1,     (21)
```

and the leading coefficients in `(20)` are nonzero.

If `y>0`, the nonzero constants in `(21)` survive.  If `y<0`, the
nonzero leading coefficients survive.  If `y=0`, coprimality says that
at least one of the two residue values survives.

The top row cannot be cancelled by `(2)`.  For `y<0`, for example,

```text
v(phi_R)-v(phi_j)
 =(R-j)x+(deg P_R-deg P_j)y <0,                       (22)

v(psi_R)-v(psi_j)
 =(R-j)x+(deg Q_R-deg Q_j)y <0                        (23)
```

for every `j<=R-2`: the first term is strictly negative and the second is
nonpositive.  For `y>=0`, the relevant exponent is zero and the strict
term `(R-j)x` remains.  At `y=0`, use whichever top residue in `(21)` is
nonzero.

Thus either `phi_Q` has a nonzero leading pole, contradicting
`phi_Q=0`, or `psi_Q` has one, contradicting `psi_Q in C`.  Every polar
point is excluded.

## 4. Pure-q cone and the unique resonance

In the second cone, `2a<nu` and `2b>=nu`.  Put

```text
delta=nu-2a>0,                    v(q)=-delta/2.        (24)
```

Here `d` is regular.  If `s` has a pole, then

```text
a+nu<3b,                                             (25)
```

because:

- for `nu=2c`, `a<=c-1` and `b>=c`;
- for `nu=2c+1`, `a<=c` and `b>=c+1`.

Inequality `(25)` is exactly

```text
v(q)/3 < v(s)/4.                                     (26)
```

So `qt^3` is strictly more pole-efficient than `st^4`; if `s` is
regular, this is even clearer.  Consequently the unique least face is the
pure-`q` monomial

```text
R=3k+1:  phi_R -> 4binom(R-1/2,4k+1)T^(2k);
R=3k+2:  H_R   -> 4binom(R-1/2,4k+3)T^(2k);
R=3k:    psi_R -> 4binom(R-1/2,4k)T^(2k).             (27)
```

All coefficients are nonzero.  The relevant coefficient index is divisible
by three.  Moving from row `R` to any retained `j<=R-2` lowers that index
by at least eight and lowers the maximal number of `q` factors by at least
three.  Together with strict efficiency `(26)`, this proves strict
top/lower separation for `(27)`.

For `R=3k+1`, `(27)` contradicts `phi_Q=0`.  For `R=3k`, it has
valuation `-2k delta<0` and contradicts `psi_Q in C`.

Let `R=3k+2`.  The actual unique face has

```text
v(H_Q)=-2k delta.                                     (28)
```

But `(12)` and `(14)` prescribe

```text
v(H_Q)=v(K_Q)-v(T)=nu+1-3a=1-a+delta.                 (29)
```

The two values agree exactly when

```text
a=1+(2k+1)delta,
nu=2+(4k+3)delta.                                     (30)
```

Off `(30)`, the pure cone is impossible.  On `(30)`, face uniqueness
only fixes the valuation; it does not satisfy the leading coefficient in
`A_srcK_Q=lambda M`.  The theorem makes no existence claim there.

The resonance is sharp at every scale.  On any ray `(30)`, with local
parameter `tau`, take

```text
V=tau^[2+(4k+3)delta],
A_src=tau^[1+(2k+1)delta],
B_src=C_0=E_0=0,                    M=tau,
Q=E_(4R-2).                                           (31)
```

Then `d=s=0`, `T=tau^(-delta)`, and, with the nonzero constant

```text
c_R=4binom(R-1/2,4k+3),
```

exact extraction gives

```text
phi_Q=psi_Q=0,
H_Q=c_R T^(2k),                    K_Q=c_R T^(2k+1),
A_src K_Q=c_R M.                                      (31a)
```

Thus every local flux and valuation identity used here is compatible on
every resonance.  This is a formal local source witness, not a global
Keller pair.  The first ray in the stated degree range is

```text
R=8, k=2, delta=1, nu=13, a=6,
c_R=-195/131072.                                      (31b)
```

MISTAKE-317 records the failed implication that led to `(31)`: uniqueness
of a face does not by itself make its valuation incompatible with the
required one.

## 5. Regular-H cone

Finally assume `2a>=nu` outside the polar cone.  Then `T` is regular.
Also `s` is regular: if `a+b<nu`, nonpolarity would force `2b>=nu`,
while `2a>=nu` would give `a+b>=nu`.

If `d` is regular, so is `Td`.  If `d` has a pole, nonpolarity gives
`a+b>=nu`, and

```text
v(Td)=(2a-nu)+(2b-nu)=2(a+b-nu)>=0.                  (32)
```

Thus `(11)` makes `H_Q` regular and `(12)` makes `K_Q` regular.  Yet
`a>=ceil(nu/2)>=2`, so `(14)` gives `v(K_Q)=1-a<0`.  This contradiction
excludes the entire third cone.

Sections 3--5 prove the exact atlas:

> For `R=0,1 mod 3`, every `nu>=3` is excluded.  For `R=3k+2`,
> every `nu>=3` is excluded except possibly the rays `(30)`.

## 6. Balanced-passport divisibility and the double-pole corollary

In the THM-2796 balanced normal form,

```text
V=v S D T^2,                      M=S E T,
D=product_j(x-beta_j)^(p_j),      sum_jp_j=N,          (33)
```

and the disjoint squarefree factors give

```text
(ord_beta_j V,ord_beta_j M)=(p_j+2,1).                (34)
```

If `R=0,1 mod 3`, `(34)` supplies an excluded point for every nonconstant
balanced passport, because such a passport has at least one pole.

If `R=3k+2`, equation `(30)` requires

```text
p_j=nu-2=(4k+3)delta_j                             (35)
```

at every pole.  Therefore a surviving passport must satisfy

```text
(4k+3) divides every p_j and N,       N>=4k+3,
h<=N/(4k+3).                                           (36)
```

This is a necessary condition, not a construction.  In particular,
`p_j=1` and `p_j=2` never survive: simple and double pole parts are
excluded in every normalized degree.  Every passport with `N` not divisible
by `4k+3`—in particular every `N<4k+3` passport—is also excluded.  For
`R=8`, the first permitted part is `11`.  This residual is not
passport-vacuous: THM-2796's complete
`e=0` chamber contains the cyclic one-pole carrier with `N=p_1=11`,
whose pole has `(ord V,ord M)=(13,1)`.  What remains open is its entry into
the Faber source chart, not existence of the abstract balanced response.

There is an exact arithmetic explanation for the recurring `2` and `3`.
On a resonance ray,

```text
(a-1)/p=(2k+1)/(4k+3),
(4k+3)-2(2k+1)=1.                                   (37)
```

After primitive reduction, `(2k+1,4k+3)` is a Farey neighbor of the
square-deck wall `(1,2)`.  Thus the `C_3` response class `R=2 mod 3`
selects one parabolic flank of the `C_2` square-root slope.  This is only
an exact interpretation of the resonance arithmetic; it does not close
the ray or identify Keller monodromy with `PSL_2(Z)`.

## 7. Scope and next decisive target

This is a local entry theorem inside the inherited nonsplit polynomial
exact-square-prefix chart.  It does not destroy the abstract rational
response maps of THM-2796, enter arbitrary Keller pairs into the chart,
handle nonpolynomial prefixes, or prove `JC(2)` or `DC(2)`.

The formal family `(31)--(31a)` already shows that the local unit equation
matches, with `lambda=c_R`; there is no remaining one-point coefficient
test.  The next honest target is global entry compatibility: can the cyclic
`N=(4k+3)delta` carrier from THM-2796 admit global polynomial
`A_src,V` with the resonance divisors and the complete Faber/source
one-form simultaneously?  Only a global divisor constraint, a next
response, or an independent chart-entry obstruction can close that lane.

## 8. Exact companion

The companion:

1. reconstructs the relevant row coefficients by exact generalized
   multinomial extraction and independently checks formula `(9)` through
   `R=18`;
2. verifies the exact degree, endpoint-coefficient, and gcd laws
   `(20)--(21)` through `R=18`;
3. exhausts `nu=3,...,25`, a rectangular hostile grid in `(a,b)`, and
   `R=7,...,18`, checking all three regions and every retained lower row;
4. compares the actual `H` face with the prescribed valuation `(29)`,
   printing the equality rays instead of counting them as exclusions;
5. checks the passport divisor `(35)` on every printed ray; and
6. contains neither a Python `assert` node nor a float literal.

The all-order theorem comes from the symbolic inequalities and formulas
above, not the finite ranges.  Run

```text
python 04-computation/jc_uniform_pole_order_faber_nonresonance_thm2827.py
python -O 04-computation/jc_uniform_pole_order_faber_nonresonance_thm2827.py
```

Normal, optimized, and stored transcripts agree exactly.

## 9. Independent hostile audit

An independent audit rederived the target translation, formula `(9)`, the
three-cone exhaustion (including `b=infinity`), every polar endpoint and
unit-residue case, pure-face separation, the iff resonance arithmetic, and
the formal local witness.  It separately checked the THM-2796 local typing,
passport divisibility, cyclic boundary, and primitive Farey interpretation.
Both execution modes byte-match the stored output, both LF hashes match,
and the companion has no truth-bearing assertion or floating literal.

**QED.**
