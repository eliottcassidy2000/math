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
  at least two and the first reduced arithmetic type is 2:3. The reduced
  nondividing cells and JC(2) remain open.
source: root / post-THM-3985 cusp-log denominator-depth lane, 2026-08-24
audit: >
  PASS (jc-mixed-generator-submersion, 2026-08-24). The independent audit
  rederived the dominant tau=0/H-adic intersection argument, the exact
  negative-coefficient cone and depth-symbol surjectivity, the first-depth
  syzygy, every Laurent convolution sign, the zero integration constant,
  the UFD common-power law and s^g conductor, and terminating target-shear
  descent. A separate generator-degree-six search found no module or
  intersection hostile. Normal, optimized, and frozen outputs byte-match at
  CHECKS=2439; hashes agree.
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
  - THM-3985-cusp-plane-rational-time-residue-and-height-two-mixed-submersion
related:
  - THM-3986-cusp-submersion-single-positive-x-monomial-adjacency-criticality
  - THM-3987-gwozdziewicz-every-line-height-two-three-weight-floor
script: 04-computation/jc2_cusp_log_laurent_conductor_thm3989.py
output: 05-knowledge/results/jc2_cusp_log_laurent_conductor_thm3989.out
script_sha256: c85b1a267ecaa14f6f892c21310d4fe5f373181391319bcdd3794da072f268e2
output_sha256: e9b27753864beee2934b29a44c061029db740dfad5ceae595d92bf49013fba25
semantic_sha256: f30679954e4f8f42fc9d8c53f17dcfe0f656eb211bd13311cd14b89cdc562984
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

The theorem proves four exact statements.

1. The depth-zero intersection and every negative coefficient are

   ```text
   B_2 intersect k[s,tau]=k[p,y],
   [tau^-j]F in s^j k[s]                    for j>=1.     (6)
   ```

   More precisely, if

   ```text
   P_j=B_2 intersect tau^-j k[s,tau],
   ```

   then coefficient extraction gives

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

## 6. Preservation ledger, scope, and reproduction

```text
source          B_2 inside the rational cusp chart
representation  finite Laurent polynomial in (s,tau)
invariant       cusp-pole depth and principal coefficient
operation       polynomial target shear / Euclidean depth descent
preserved       B_2 membership and constant Jacobian
destroyed       individual support labels below the cancelled top pole
sidecar         common base h, conductor s^g, and scalar moment -s
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
sha256sum 04-computation/jc2_cusp_log_laurent_conductor_thm3989.py \
  05-knowledge/results/jc2_cusp_log_laurent_conductor_thm3989.out
python3 agents/check_docs.py
```

**QED.**
