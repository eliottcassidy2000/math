---
id: THM-3989
title: "Cusp-log Laurent conductor, scalar moment, and nondividing-depth reduction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the exact
  cusp-log chart s=y/p=xt, tau=H/p^2=t, the height-two completion has
  B_2 intersect k[s,tau]=k[p,y], and every tau^-d Laurent coefficient is
  divisible by s^d; the associated depth-d symbol module is exactly
  s^d*k[s]. The entire order-one pole-cancellation syzygy descends back to
  k[p,y]. For a hypothetical Darboux pair, both cusp-pole depths are
  positive, their leading symbols are powers of one polynomial h with
  s^g|h, and the full Laurent equation has exact scalar moment
  sum_i i*a_i*c_-i=-s. Elementary target shears reduce every hypothetical
  pair to positive depths neither of which divides the other, so both are
  at least two and the first reduced arithmetic type is 2:3. At that first
  cell the p=0 normalization seam computes an exact depth-one extension
  class and square-root-lift criterion; the nonliftable simple-base slice
  forces F_p(0,0)=lambda^2/12. The reduced nondividing cells and JC(2)
  remain open.
source: root / post-THM-3985 cusp-log denominator-depth lane, 2026-08-24
audit: >
  PASS (jc-mixed-generator-submersion plus independent every-line hostile
  audit, 2026-08-24). The audits
  rederived the dominant tau=0/H-adic intersection argument, the exact
  negative-coefficient cone and depth-symbol surjectivity, the first-depth
  syzygy, every Laurent convolution sign, the zero integration constant,
  the UFD common-power law and s^g conductor, and terminating target-shear
  descent. A separate generator-degree-six search found no module or
  intersection hostile. Normal, optimized, and frozen outputs byte-match at
  CHECKS=2439. A no-import 13,876-gate path separately checks the bounded
  intersection, syzygy kernels, conductor cone, convolution signs, UFD rows,
  and shear descent. It also repaired the description of the outside-B2
  Laurent hostile. A further 13,007-gate normal/optimized exact companion
  checks the full p=0 seam, the extension class, the first four 2:3 rows,
  the forced p-jet, and the c_2(0) hostile to a discarded F_y claim. Hashes
  agree.
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
  - THM-3985-cusp-plane-rational-time-residue-and-height-two-mixed-submersion
related:
  - THM-3986-cusp-submersion-single-positive-x-monomial-adjacency-criticality
  - THM-3987-gwozdziewicz-every-line-height-two-three-weight-floor
  - THM-3990-componentwise-harmonic-obstruction-and-repair-quotient
script: 04-computation/jc2_cusp_log_laurent_conductor_thm3989.py
output: 05-knowledge/results/jc2_cusp_log_laurent_conductor_thm3989.out
script_sha256: 63215ca6652c620f795af49e2112344df1cfe714044ca753b2a2bed4c22a93ed
output_sha256: e9b27753864beee2934b29a44c061029db740dfad5ceae595d92bf49013fba25
semantic_sha256: f30679954e4f8f42fc9d8c53f17dcfe0f656eb211bd13311cd14b89cdc562984
independent_audit_script: 04-computation/jc2_cusp_log_laurent_conductor_thm3989_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_cusp_log_laurent_conductor_thm3989_independent_audit.out
independent_audit_script_sha256: 8183206aa526cb2b95fd689a0242ba1be150542677fce56245344bae84013b74
independent_audit_output_sha256: 1f5883e70aca0d9efb81095447cb3e4633d8b4b6cf366ca4b871f2c9121e144
independent_audit_semantic_sha256: ef8e587b10aaa8398d69769411ee828985b50311c382cc7e13fe3235beccd612
depth23_supplement_script: 04-computation/jc2_cusp_log_depth23_extension_seam_thm3989_supplement.py
depth23_supplement_output: 05-knowledge/results/jc2_cusp_log_depth23_extension_seam_thm3989_supplement.out
depth23_supplement_script_sha256: ba52ebbb910cbe39a3542f052b317d92c5c22edeb84ffcffec99921cd64d1b67
depth23_supplement_output_sha256: 2a413b04fb083a5efd667877307ebd7602739e46a1eb6b87c44199d5ffecc312
depth23_supplement_semantic_sha256: 2e57c7d6d59027ebb0cb1d48dd1464a2f182fe6066c5225574db7b775eaaa815
hash_basis: raw LF bytes
---

# THM-3989 -- the missing cusp has an exact Laurent conductor

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. In the height-two
completion of THM-3973 put

```text
H=p^3-y^2,
s=y/p=xt,                       tau=H/p^2=t.              (1)
```

Then

```text
p=s^2+tau,          y=s(s^2+tau),
x=s/tau,            u=x^2t=s^2/tau,                     (2)
dx wedge dt=ds wedge dtau/tau.                           (3)
```

Thus `B_2=k[x,u,p,y]` is a Laurent subring of
`k[s,tau,tau^-1]`. For a nonzero `F in B_2`, write its unique finite
expansion

```text
F=sum_i f_i(s) tau^i,                    f_i in k[s],     (4)
```

and define its cusp-pole depth by

```text
d(F)=max(0,-min{i:f_i!=0}).                              (5)
```

The theorem proves six exact statements.

1. The depth-zero intersection and every negative coefficient are

   ```text
   B_2 intersect k[s,tau]=k[p,y],
   [tau^-j]F in s^j k[s]                    for j>=1.     (6)
   ```

   More precisely, if

   ```text
   P_j=B_2 intersect tau^-j k[s,tau],
   ```

   then coefficient extraction gives, as `k`-vector spaces and as
   `k[p,y]`-modules via cusp restriction,

   ```text
   P_j/P_(j-1) isomorphic to s^j k[s]        for j>=1.   (7)
   ```

2. For `G,K in k[p,y]`, an order-one combination loses its cusp pole if
   and only if it already descends to the cusp plane:

   ```text
   xG+uK in k[p,y]
   iff pG+yK in (H)
   iff G=p^2L+yM, K=-yL-pM for some L,M in k[p,y],

   and then xG+uK=yL.                                    (8)
   ```

3. If `A,C in B_2` hypothetically obey `J_(x,t)(A,C)=1`, write

   ```text
   A=sum_i a_i(s)tau^i,             C=sum_j c_j(s)tau^j. (9)
   ```

   Their complete coefficient equations and scalar moment are

   ```text
   sum_(i+j=k)(j a_i'c_j-i a_i c_j')=delta_(k,0),        (10)
   sum_i i a_i c_(-i)=-s.                               (11)
   ```

4. Every such hypothetical pair is equivalent under finitely many
   elementary polynomial target automorphisms to another pair in `B_2`
   whose two depths `d,e` satisfy

   ```text
   d,e>=2,                    d does not divide e,
                              e does not divide d.       (12)
   ```

   If `g=gcd(d,e)` and `a=a_(-d),c=c_(-e)` are the leading
   symbols of a reduced pair, then

   ```text
   a=alpha h^(d/g),          c=beta h^(e/g),
   alpha,beta in k*,         s^g divides h.              (13)
   ```

5. Every Laurent row in `B_2` obeys the normalization seam

   ```text
   sum_i (-1)^i [s^(ell-2i)]f_i=0               for ell>=1.
   ```

   For `P_1 -> s k[s]`, its exact extension class in
   `k[s]/k[s^2,s^3]` is

   ```text
   theta(h)=[s^3]h*[s].
   ```

   Hence a depth-two `A` with `a_-2=h^2` has a square approximate root
   modulo `P_0` exactly when `h|a_-1` and

   ```text
   [s](a_-1/(2h))=[s^3]h.
   ```

6. In the first nondividing cell, the negative coefficient rows begin with
   a formal common square/cube root over `k(s)`, but the class `theta` is
   extra filtered data. In particular, if

   ```text
   A=x^2+lambda*u+F(p,y),                    lambda!=0,
   ```

   has a depth-three `C in B_2` for which all negative Jacobian rows vanish,
   then necessarily

   ```text
   F_p(0,0)=lambda^2/12.
   ```

The smallest unordered depth pair allowed by `(12)` is `(2,3)`. This is a
necessary arithmetic normal form, not a construction of that cell.

## 1. The exact logarithmic chart

THM-3985 gives `y/p=xt` and `p=t+(xt)^2`. Hence `(1)--(2)` follow directly.
Since `x=s/tau` and `t=tau`, differentiation gives

```text
dx=(tau ds-s dtau)/tau^2,
dx wedge dtau=ds wedge dtau/tau,                         (14)
```

which is `(3)`. Consequently

```text
J_(x,t)(A,C)=tau(A_s C_tau-A_tau C_s).                   (15)
```

The missing cusp `H=0` has therefore become the logarithmic divisor
`tau=0`; its residue coordinate is literally `s`.

## 2. The Laurent conductor is exactly s^j k[s]

It is enough to inspect a generator monomial. For `a,b,c,e>=0`, equations
`(2)` give

```text
x^a u^b p^c y^e
 =s^(a+2b+e) tau^(-a-b)(s^2+tau)^(c+e).                 (16)
```

Suppose its `tau^-j` coefficient is nonzero. The selected binomial index is

```text
ell=a+b-j,                   0<=ell<=c+e.                (17)
```

The corresponding power of `s` is

```text
E=a+2b+e+2(c+e-ell)=2c+3e-a+2j.                         (18)
```

Condition `(17)` gives `a<=c+e+j-b`, and therefore

```text
E-j>=c+2e+b>=0.                                         (19)
```

Every negative coefficient in `(16)` is thus divisible by the asserted
power `s^j`; linear combinations prove the inclusion in `(6)`.

It is sharp coefficient by coefficient. At depth `j`, the offsets zero and
one are supplied by

```text
x^j                    and                    x^(j-1)u.  (20)
```

Every integer `n>=2` lies in the cusp semigroup `<2,3>`, say `n=2c+3e`.
Then `x^j p^c y^e` has principal coefficient `s^(j+n)`. These rows span all
of `s^j k[s]`, proving `(7)` as well as sharpness of `(6)`.

### 2.1 Depth zero is exactly the cusp plane

The reverse inclusion in the first row of `(6)` is the delicate direction.
Using THM-3985,

```text
B_2=k[p,y,yp/H,y^2/H],                                  (21)
```

so every `F in B_2` can be written `N/H^D` with `N in k[p,y]`. Along the
dominant component `tau=0`, the generic point has `p=s^2!=0`; because

```text
H=tau p^2,                                               (22)
```

the `tau`-adic valuation there is exactly the `H`-adic valuation on
`k(p,y)`. If `F` has no negative `tau` power, this valuation is nonnegative.
The irreducible cusp polynomial `H` is prime in the UFD `k[p,y]`, so all
denominator powers cancel and `F in k[p,y]`. Conversely `p,y` are visibly
polynomial in `s,tau`. This proves the intersection equality.

The word *dominant* matters in this argument: `(22)` also contains the
contracted locus `p=0`, but the divisor valuation used above is at the
generic point of `H=0`, where `p` is a unit.

## 3. The complete first-depth descent syzygy

For `G,K in k[p,y]`, equations `(21)` give

```text
xG+uK=y(pG+yK)/H.                                       (23)
```

Since `gcd(H,y)=1`, this expression belongs to `k[p,y]` exactly when
`H` divides `pG+yK`. Write `pG+yK=HL`. Then

```text
p(G-p^2L)+y(K+yL)=0.                                    (24)
```

Coprimality of `p,y` gives a polynomial `M` with

```text
G=p^2L+yM,                    K=-yL-pM.                 (25)
```

Conversely `(25)` implies `(24)`, and substitution in `(23)` gives `yL`.
This proves `(8)`. The relation `xy=up` is the `L=0` syzygy; the other
generator `xp^2-uy=y` is the `M=0,L=1` descent. Thus a coordinated first
layer can cancel its pole, but only by returning to the already closed
cusp-plane subring.

## 4. The logarithmic coefficient law and scalar moment

Substitute `(9)` in `(15)`. The coefficient of `tau^k` is

```text
sum_(i+j=k)(j a_i'c_j-i a_i c_j').                      (26)
```

Equating it with one proves `(10)`. At `k=0`, put `j=-i`; then `(26)` is

```text
-d/ds [sum_i i a_i c_(-i)]=1.                           (27)
```

Every summand in brackets vanishes at `s=0`: for `i>0`, row `(6)` applies
to `c_(-i)`, and for `i<0` it applies to `a_i`. The integration constant in
`(27)` is therefore zero, proving `(11)`.

The first-order part of this identity is the exact log-jet determinant

```text
a_1(0)c_(-1)'(0)-a_(-1)'(0)c_1(0)=-1.                  (28)
```

Thus the scalar bracket is carried by opposite shallow Laurent channels
even when both coordinates have much deeper leading poles. Equation `(11)`
retains their coupling; a scalar count of pole depths alone would lose it.

## 5. Common leading powers and target-shear descent

Suppose `(9)` is a Darboux pair and put `d=d(A),e=d(C)`. Neither depth can
be zero. Indeed, Section 2.1 would put that nonconstant coordinate in
`k[p,y]`, while THM-3985 excludes every nonconstant cusp-plane element from
a polynomial Keller pair.

Let

```text
a=a_(-d),                      c=c_(-e).                 (29)
```

The coefficient of `tau^(-d-e)` in `(10)` must vanish, giving

```text
d a c'-e a'c=0.                                        (30)
```

Put `g=gcd(d,e)`. Unique factorization in `k[s]` turns `(30)` into the
common-power law

```text
a=alpha h^(d/g),                 c=beta h^(e/g).        (31)
```

By `(6)`, `s^d|a` and `s^e|c`; either row of `(31)` then gives `s^g|h`.
This proves `(13)` before reduction as well as after it.

Now suppose `d|e` and write `e=qd`. In `(31)` one has

```text
a=alpha h,                       c=beta h^q.             (32)
```

The elementary target shear

```text
C_1=C-(beta/alpha^q)A^q                                 (33)
```

cancels the entire `tau^-e` coefficient, so `d(C_1)<e`. It preserves both
membership in `B_2` and the bracket `J(A,C_1)=1`. The new depth cannot be
zero, by the same THM-3985 argument. If instead `e|d`, apply the symmetric
shear to `A`. Whenever the depths are equal, a linear instance of `(33)`
lowers one of them.

Iterating strictly decreases the positive integer `d+e`, so it terminates.
At termination neither depth divides the other, and in particular neither
equals one. This proves `(12)`. Equivalently, the reduced exponents

```text
d/g, e/g                                                   (34)
```

are coprime integers at least two. The first possible unordered pair is
`2:3`; larger possibilities include `2:5`, `3:4`, and nonprimitive depth
pairs such as `(4,6)`, whose coprime quotient is again `2:3`.

The missing object at this frontier is now exact. The formal common base
`h tau^-g` would continue the Euclidean cancellation, but `(31)` does not
put that approximate root in `B_2`. A proof or counterexample must decide
whether the conductor `s^g|h`, the moment `(11)`, and the higher convolution
rows force such a lift or permit a genuine nondividing pair.

## 6. The p=0 seam and the exact depth-one extension

The divisor chart `tau=0` does not retain the contracted color `p=0`. To
recover it, set `tau=-s^2` in `(2)`. Then

```text
(x,u,p,y)=(-s^-1,-1,0,0).                               (35)
```

Consequently every `F` in `(4)` specializes to a polynomial in `s^-1` and
has no positive `s` power. Extracting `[s^ell]` proves, for every `ell>=1`,

```text
sum_i (-1)^i [s^(ell-2i)]f_i=0,                         (36)
```

where a negative requested exponent contributes zero. The first seam is

```text
[s]f_0=sum_(j>=1)(-1)^(j+1)[s^(2j+1)]f_-j.             (37)
```

Now use the exact sequence supplied by `(7)`:

```text
0 -> P_0 -> P_1 --sigma--> s k[s] -> 0,
sigma(F)=[tau^-1]F.                                     (38)
```

If `sigma(F)=h`, the class of `[tau^0]F` modulo the cusp restriction
`k[s^2,s^3]` is independent of the lift. The following quotient is taken as a
`k`-vector space. Since

```text
k[s]/k[s^2,s^3]=k*[s],
```

equation `(37)` gives the exact extension class

```text
theta(h)=[s^3]h*[s].                                    (39)
```

This is also sufficient, not merely necessary. Monomial symbols have the
explicit lifts

```text
s <- x,       s^2 <- u,       s^(1+2c) <- x p^c,
s^(4+2c) <- x p^c y.                                   (40)
```

Only `xp=s^3 tau^-1+s` has a nonzero missing-`s` residue. Hence, for prescribed
`h in s k[s]` and `q in k[s]`, there is `R in P_1` with negative coefficient
`h` and constant coefficient `q` if and only if

```text
[s]q=[s^3]h.                                             (41)
```

For sufficiency, start with any lift `R_0` of `h`; condition `(41)` puts
`q-[tau^0]R_0` in `k[s^2,s^3]`, so adding the corresponding `G(p,y)` gives
the required lift.

Let `A in P_2` have `a_-2=h^2`. Comparing the two negative coefficients of
`A` and `R^2` now proves the exact square approximate-root criterion:

```text
there is R in P_1 with [tau^-1]R=h and A-R^2 in P_0
iff h divides a_-1 and [s](a_-1/(2h))=[s^3]h.            (42)
```

The smallest hostile is

```text
A=x^2+2u=s^2 tau^-2+2s^2 tau^-1.                        (43)
```

Here `h=s` and `a_-1/(2h)=s`, so `(41)` would read `1=0`. Thus the sharp
symbol module `s k[s]` does not by itself produce a filtered common root.

## 7. The first 2:3 rows and a forced simple-base jet

At depths `(2,3)`, determinant-one target scaling together with a constant
rescaling of the common base normalizes

```text
a_-2=h^2,                         c_-3=h^3.              (44)
```

Indeed, if the original leading coefficients are `alpha*g^2,beta*g^3`, choose
`eta^5=alpha*beta`, put `h=eta*g`, and scale the two targets reciprocally.

The next three negative rows integrate exactly. A linear shear of `C` by `A`
and constant translations of `A,C` kill the three integration constants. With
`q,r,w` initially in the rational field `k(s)`, the result is

```text
a_-1=2hq,                    c_-2=3h^2q,
a_0=q^2+2hr,                 c_-1=3hq^2+3h^2r,
a_1=2hw+2qr,                 c_0=q^3+6hqr+3h^2w.         (45)
```

Before the constants are killed, the three conserved expressions are

```text
(2c_-2/h^2-3a_-1/h)',
(2c_-1/h-3q^2-3a_0)',
(2c_0-3ha_1-6hqr-2q^3)'.                                (46)
```

The weight `-1` row then has general solution

```text
c_1=(3/2)h a_2+(3/2)h r^2+3hqw+3q^2r+kappa/h.          (47)
```

The first four terms are the next coefficient of the formal identity
`C=R^3` when `A=R^2`; the last is the first commuting mismatch
`kappa*R^-1`. Equivalently, adding `delta` to the displayed common-root
coefficient changes the row by `2h(h delta)'`. Crucially, `(45)--(47)` only
construct `q,r,w` in `k(s)`. Polynomiality and membership in `B_2` require
the filtered seams from Section 6 and their higher analogues.

The first seam tax can be computed completely in the nonliftable simple-base
slice

```text
A=x^2+lambda*u+F(p,y),                       lambda!=0. (48)
```

Put

```text
phi=F(s^2,s^3),
psi=F_p(s^2,s^3)+s F_y(s^2,s^3).                         (49)
```

If a depth-three `C in B_2` has `c_-3=beta*s^3`, `beta!=0`, and all negative
Jacobian rows vanish, a linear target shear and the first three integrations
give one scalar `K` and

```text
c_-2=(3/2)beta lambda s^3,
c_-1=(3/2)beta s(phi+lambda^2s^2/4)+Ks,
c_0=beta((3/2)s psi+(3/4)lambda s phi-(1/16)lambda^3s^3)
    +(Klambda/2)s.                                      (50)
```

The `ell=1` instance of `(36)` forces

```text
K=(3/4)beta lambda-(3/2)beta F(0,0).                    (51)
```

The weight `-1` row determines `[s]c_1`. Substitution into the `ell=3` seam
then leaves exactly

```text
beta*lambda*(lambda^2-12F_p(0,0))/32=0.                 (52)
```

Thus `F_p(0,0)=lambda^2/12` because `beta*lambda!=0`. This is necessary only:
it does not construct `C`, solve the scalar moment, or show that the forced
cusp-plane jet sequence terminates.

An exploratory continuation once appeared to force `F_y(0,0)=0` at seam
`ell=4`; that inference is **retracted**. At `ell=4`, the positive Laurent
coefficient contributes the free term `[s^0]c_2`, which absorbs the entire
apparent condition. The supplement has a dedicated regression gate for this
first positive-coefficient boundary.

## 8. Preservation ledger, scope, and reproduction

```text
source          B_2 inside the rational cusp chart
representation  finite Laurent polynomial in (s,tau)
invariant       cusp-pole depth and principal coefficient
operation       polynomial target shear / Euclidean depth descent
preserved       B_2 membership and constant Jacobian
destroyed       individual support labels below the cancelled top pole
sidecar         common base h, conductor s^g, seam theta, and moment -s
cheapest test   first nondividing depth pair 2:3
```

This theorem does not close the `2:3` depth cell, construct a Darboux pair,
or prove `JC(2)`. Its gain is a complete coordinate system for the remaining
denominator problem: the cusp plane is exactly depth zero, every negative
coefficient has a sharp conductor, and all divisible depth pairs are
procedurally reducible.

The companion verifies the chart, generator monomial formula, negative
coefficient divisibility through broad exact grids, symbol surjectivity for
depths `1..9`, depth-zero finite hostiles, the first-depth syzygy, every
Laurent convolution sign, the scalar moment and tame automorphism hostile,
and common-power rows across divisible and nondividing depth pairs. Finite
rows are controls for the symbolic proof, not extrapolation.

Reproduce with

```bash
python3 04-computation/jc2_cusp_log_laurent_conductor_thm3989.py
python3 -O 04-computation/jc2_cusp_log_laurent_conductor_thm3989.py
python3 04-computation/jc2_cusp_log_laurent_conductor_thm3989_independent_audit.py
python3 -O 04-computation/jc2_cusp_log_laurent_conductor_thm3989_independent_audit.py
python3 04-computation/jc2_cusp_log_depth23_extension_seam_thm3989_supplement.py
python3 -O 04-computation/jc2_cusp_log_depth23_extension_seam_thm3989_supplement.py
sha256sum 04-computation/jc2_cusp_log_laurent_conductor_thm3989.py \
  05-knowledge/results/jc2_cusp_log_laurent_conductor_thm3989.out \
  04-computation/jc2_cusp_log_laurent_conductor_thm3989_independent_audit.py \
  05-knowledge/results/jc2_cusp_log_laurent_conductor_thm3989_independent_audit.out \
  04-computation/jc2_cusp_log_depth23_extension_seam_thm3989_supplement.py \
  05-knowledge/results/jc2_cusp_log_depth23_extension_seam_thm3989_supplement.out
python3 agents/check_docs.py
```

**QED.**
