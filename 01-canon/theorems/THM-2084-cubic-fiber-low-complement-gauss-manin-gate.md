---
id: THM-2084
title: "Cubic-fiber Faber gate through reduced degree eleven"
status: >
  PROVED. Let one output-pencil member P of a complex planar Keller pair
  have degree three along a fixed linear source fiber, and minimize the fiber
  degree of its complement modulo target shears by C[P]. The reduced degree
  is never divisible by three. If it is at most eleven, then it is one or two
  and the pair is tame; the residual degrees four, five, seven, and eight are
  impossible, as are ten and eleven. Hence a hypothetical non-tame planar
  Keller pair with such a cubic member has reduced complementary fiber degree
  at least thirteen. The proof
  identifies the coefficient ladder with the polynomial parts of
  (z^3+pz+q)^(n/3), uses exact boundary resultants to remove every finite pole,
  and then uses the first two omitted Laurent coefficients to rule out the
  polynomial branches at infinity. This is a source-fiber/pencil theorem, not
  a generic cover-degree or full planar-Jacobian theorem.
source: codex-2026-07-22-JC2-cubic-fiber
depends_on:
  - THM-2063
  - THM-2071
related:
  - THM-1345
  - THM-1330
  - THM-2097
  - MISTAKE-229
script: 04-computation/jc2_cubic_fiber_gauss_manin_gate_codex_20260722.py
output: 05-knowledge/results/jc2_cubic_fiber_gauss_manin_gate_codex_20260722.out
script_sha256: 9d42027cb54ddbf2b3a0b649c898f8de9d77ba48a8ec585f5236b9d470026e20
output_sha256: ab12c5b06fbfa1392b115cab5699bf80b1b96bd7f8388e7aa3a083f71d6e487b
hash_basis: repository blobs with LF line endings
---

# THM-2084 -- cubic-fiber Faber gate through reduced degree eleven

## 1. Statement

Let

```text
P=A(x)y^3+B(x)y^2+C(x)y+D(x),        A != 0,
Q in C[x,y],                         J(P,Q)=kappa in C*,
```

and define, in the chosen linear source-fiber direction,

```text
mu = min_H deg_y(Q-H(P)),            H in C[T].             (1)
```

Then:

1. `3` does not divide `mu`;
2. `mu` is not any of `4,5,7,8,10,11`;
3. if `mu<=11`, then `mu` is `1` or `2` and the pair is tame;
4. consequently a non-tame pair in this cubic-fiber stratum has
   `mu>=13`.

The last implication uses THM-2063 at `mu=1` and THM-2071, applied with the
two components interchanged, at `mu=2`.  The new content is the exact cubic
connection and the exclusions of `4,5,7,8,10,11`.

## 2. The cube gate and depressed cubic

Choose `H` attaining (1), replace `Q` by `Q-H(P)`, and write

```text
Q=sum_(j=0)^n q_j(x)y^j,             q_n != 0, n=mu.
```

The coefficient of `y^(n+2)` in the Jacobian equation is

```text
n A' q_n-3A q_n'=0.                                      (2)
```

Thus

```text
(q_n^3/A^n)'=0,                 q_n^3=c A^n, c in C*.     (3)
```

If `n=3m`, then `q_n=d A^m` for a constant `d`, and subtracting `dP^m`
lowers the fiber degree without changing the Jacobian.  This contradicts
minimality.  Hence `3` does not divide `n`.

Since `gcd(n,3)=1`, the valuation of every irreducible factor of `A` in (3)
is divisible by three.  Absorbing a nonzero scalar cube, write

```text
A=U^3,                         U in C[x] nonzero.          (4)
```

Over `K=C(x)` put

```text
h=B/(3U^2),                   z=Uy+h,
p=C/U-3h^2,                  q=D-h^3-ph.                  (5)
```

Then

```text
P=z^3+pz+q.                                                (6)
```

The change `y -> z` has determinant `U`, so the Jacobian equation is

```text
L(Q)=r,
L=(p'z+q') partial_z-(3z^2+p) partial_x,    r=kappa/U,    (7)
```

where `partial_x` holds `z` fixed.  Notice that `L(P)=0`.

## 3. The rank-two connection

Put `t=P`.  The monicity of (6) gives a unique decomposition

```text
Q=AA(x,t)+z BB(x,t)+z^2 CC(x,t),       AA,BB,CC in K[t].  (8)
```

Reducing `L(Q)` modulo `z^3=t-pz-q`, its coefficients in the basis
`1,z,z^2` are

```text
1:   -p AA_x+3(q-t)BB_x+q'BB,
z:    2p BB_x+3(q-t)CC_x+p'BB+2q'CC,
z^2: -3 AA_x+2p CC_x+2p'CC.                              (9)
```

The subscripts in (9) hold `t` fixed.  The last equation in `L(Q)=r` gives

```text
AA=(2/3)p CC+H(t).                                        (10)
```

The term `H(P)` is another target shear and may be discarded.  Equations
(9)--(10) are the promised rank-two Gauss--Manin system.

There is a more revealing all-degree description.  For `3` not dividing
`m`, define the cubic Faber representative

```text
E_m(z;p,q)=Pol_z (z^3+pz+q)^(m/3),                       (11)
```

where `Pol_z` keeps the nonnegative powers in the Laurent expansion at
`z=infinity`.  Explicitly,

```text
E_m=sum_(2i+3j<=m)
      binom(m/3,i+j) binom(i+j,i) p^i q^j z^(m-2i-3j).   (12)
```

The discarded tail starts at `z^-1`.  Since `L(P^(m/3))=0`, applying `L` to
that tail can create at most powers `z^1` and `z^0`.  Therefore

```text
L(E_m)=z Phi_m'+R_m',                                    (13)
Phi_m=3 [z^-1] P^(m/3),        R_m=3 [z^-2] P^(m/3).     (14)
```

This proves the functional form of the coefficient ladder.  Conversely,
if `L(S)` has `z`-degree at most one, compare the top coefficient of `S`.
If `S=s_m(x)z^m+...`, the coefficient of `z^(m+2)` is exactly `-3s_m'`,
so `s_m` is constant.  More generally, the constants of `partial_x` on
`C(x)[t]` are precisely `C[t]`.
For top degree divisible by three subtract a constant power of `P`; otherwise
subtract the corresponding constant multiple of `E_m`.  Repetition gives a
unique linear combination of the `E_m`, modulo `C[P]`.  This is the recursive
normal-form argument behind the finite calculations below.

For the first degrees needed here:

```text
E_1 = z,
E_2 = z^2+(2/3)p,

E_4 = z^4+(4/3)pz^2+(4/3)qz+(2/9)p^2,

E_5 = z^5+(5/3)pz^3+(5/3)qz^2+(5/9)p^2z+(10/9)pq,

E_7 = z^7+(7/3)pz^5+(7/3)qz^4+(14/9)p^2z^3
      +(28/9)pqz^2+((14/81)p^3+(14/9)q^2)z+(14/27)p^2q,

E_8 = z^8+(8/3)pz^6+(8/3)qz^5+(20/9)p^2z^4
      +(40/9)pqz^3+((40/81)p^3+(20/9)q^2)z^2
      +(40/27)p^2qz-(10/243)p^4+(40/27)pq^2.             (15)
```

The two tail coefficients in (14) are

```text
m    Phi_m                                      R_m

1    p                                          q
2    2q                                         -p^2/3
4    (4/3)pq                                    -(4/27)p^3+(2/3)q^2
5    -(5/27)p^3+(5/3)q^2                       -(5/9)p^2q
7    -(7/81)p^4+(14/9)pq^2                     -(28/81)p^3q+(14/27)q^3
8    -(40/81)p^3q+(40/27)q^3                   (8/243)p^5-(20/27)p^2q^2
10   (70/243)pq(-p^3+6q^2)                     (35/2187)(p^6-36p^3q^2+27q^4)
11   (22/2187)(2p^6-90p^3q^2+135q^4)          (44/729)p^2q(2p^3-15q^2).
                                                               (16)
```

Thus a reduced solution of degree `n` in `{4,5,7,8,10,11}` has the form

```text
n=4: Q=cE_4+dE_2+eE_1,
n=5: Q=cE_5+dE_4+eE_2+fE_1,
n=7: Q=cE_7+dE_5+eE_4+fE_2+gE_1,
n=8:  Q=cE_8+dE_7+eE_5+fE_4+gE_2+kE_1,
n=10: Q=cE_10+dE_8+eE_7+fE_5+gE_4+kE_2+lE_1,
n=11: Q=cE_11+dE_10+eE_8+fE_7+gE_5+kE_4+lE_2+mE_1,
                                                        c != 0. (17)
```

after a target shear.  If `Phi` and `R` denote the same linear combinations
of the two columns in (16), equation (7) becomes exactly

```text
Phi'=0,                         R'=kappa/U.              (18)
```

## 4. Polynomiality removes the centering pole

It remains essential that the original `Q` lies in `C[x,y]`, not merely in
`K[z]`.  Let `v` be the valuation at a finite point and suppose first that

```text
v(h)=-H<0.                                                 (19)
```

Because the original constant fiber coefficient is polynomial,

```text
D=P(x,0)=h^3+ph+q                                       (20)
```

has nonnegative valuation.  Write `rho=-v(p)` and compare `rho` with `2H`.

If `rho>2H`, then `q~-ph`.  The unique lowest-valuation term in the top
first integral `Phi_n` is, respectively,

```text
n=4: pq,       n=5: p^3,       n=7: p^4,       n=8: p^3q,
n=10: p^4q,    n=11: p^6.
```

It cannot occur in a constant.  If `rho<2H`, then `q~-h^3`.  For `n=4`,
`E_4(h)` has the unique leading term `-h^4/3`; for `n=5`, `Phi_5` has the
unique leading term `(5/3)q^2`; for `n=8`, `Phi_8` has the unique leading
term `(40/27)q^3`; for `n=11`, `q^4` is unique.  For `n=7`, the term
`p q^2` is unique when `p` has a pole, while for regular `p` the direct
boundary value `E_7(h)~(2/9)h^7` is nonzero.  The same split at `n=10` uses
the unique `p q^3` term for polar `p`, and the nonzero value
`E_10(1;0,-1)` for regular `p`.
Every case contradicts either `Phi'=0` or polynomiality of `Q(x,0)`.

The balanced case `rho=2H` contains the only possible cancellation.  Write

```text
p~a h^2,                         q~-(1+a)h^3.             (21)
```

The top first-integral equation requires
`phi_n(a)=Phi_n(a,-1-a)=0`, while polynomiality of `Q(x,0)` requires
`e_n(a)=E_n(1;a,-1-a)=0`.  Clearing the displayed rational denominators, the
exact numerator pairs used by the referee and their resultants are

```text
n   e_n numerator                         phi_n numerator                    resultant

4   2a^2-3                               -4a(a+1)                            48
5   -(5a^2+10a+6)                        -5(a^3-9(a+1)^2)                   -11025
7   -2(14a^3+21a^2-9)                    -7a(a^3-18(a+1)^2)                 13070456784
8   -5(2a^4-24a^3-72a^2-72a-27)         40(a+1)(a^3-3(a+1)^2)             21902400000000
10  -7(2a^5-45a^4-120a^3-90a^2+18)      70a(a+1)(a^3-6a^2-12a-6)          -285658406085931200000
11  4(22a^5-55a^4-330a^3-495a^2
      -330a-90)                          22(2a^6-90a^5-45a^4+450a^3
                                             +810a^2+540a+135)             -1392286514585108181811200000.
                                                                        (22)
```

All six resultants are nonzero.  Their factorizations are

```text
2^4*3,
-3^2*5^2*7^2,
2^4*3^9*7^3*11^2,
2^12*3^4*5^8*13^2,
-2^9*3^7*5^5*7^10*17^2,
-2^20*3^9*5^5*11^5*13^5*19^2,
```

so (22) is also an exact Euclidean-algorithm certificate, not a numerical
root test.  This eliminates (21).  We conclude

```text
h in C[x].                                                  (23)
```

## 5. Polynomiality removes the depressed-coefficient poles

With `h` polynomial, (5) gives

```text
Up in C[x],                         q+hp in C[x].           (24)
```

Suppose `p` has a finite pole of order `rho>0`.  Then `q=-hp+O(1)` has pole
order at most `rho`.

- At `n=4`, if `q` has a pole then the product `pq` uniquely dominates the
  linear terms in `Phi`; if `q` is regular, `(2/9)p^2` uniquely dominates
  `Q(x,0)`.
- At `n=5`, the `-(5/27)p^3` term uniquely dominates `Phi`.
- At `n=7`, the `-(7/81)p^4` term uniquely dominates `Phi`.
- At `n=8`, the `-(10/243)p^4` term uniquely dominates `E_8(h)`; all other
  terms have at most cubic order in `p` because `q=O(p)`.
- At `n=10`, the top `p^5` term uniquely dominates `E_10(h)`; every other
  monomial has ordinary `(p,q)`-degree at most four.
- At `n=11`, the `p^6` term uniquely dominates `Phi_11`.

Each alternative is impossible.  Hence `p` is polynomial, and (24) then
makes `q` polynomial.  Now `R` in (18) is polynomial, so

```text
kappa/U=R' in C[x].                                       (25)
```

Since `kappa` is a nonzero constant, (25) forces `U` to be constant.  Thus
`z=Uy+h` is an honest triangular polynomial coordinate and the remaining
problem is entirely in `p,q in C[x]`, with `R` affine nonconstant.

## 6. The polynomial branches at infinity are empty

We now combine

```text
Phi(p,q)=constant,                  R(p,q)=linear nonconstant. (26)
```

### Degree four

Here

```text
Phi=(4c/3)pq+e p+2d q.
```

After translating `p,q`, this is a nonzero scalar times
`(p+alpha)(q+beta)` plus a constant.  If its fixed value is nonzero, both
factors are polynomial units and hence constant.  If it is zero, one factor
vanishes identically.  On the remaining branch, `R` has a nonzero quadratic
term in `q` or a nonzero cubic term in `p`, so it cannot be affine.  Thus
degree four is impossible.

### Degree five

Let `a=deg p` and `b=deg q`.  If exactly one is nonconstant, the `p^3` or
`q^2` term in `Phi` is unique.  If both are nonconstant, constancy of `Phi`
forces

```text
3a=2b,                         a=2s, b=3s, s>=1,
q_lead^2=p_lead^3/9.                                      (27)
```

But the leading term of `R` is then the nonzero `-(5c/9)p^2q`, of degree
`7s`, not one.  Degree five is impossible.

### Degree seven

If `q` is constant and `p` is not, the `p^4` term in `Phi` is unique.  If
`p` is constant and `q` is not, `R` has the unique cubic term
`(14c/27)q^3`.  When both are nonconstant, the highest degrees in `Phi` are
`4a` and `a+2b`; equality again forces `a=2s,b=3s`, now with

```text
q_lead^2=p_lead^3/18.                                    (28)
```

Modulo (28), the top part of `R` is

```text
-(77/243)p_lead^3 q_lead != 0,
```

of degree `9s`.  Degree seven is impossible.

### Degree eight

If either `p` or `q` is constant, `Phi` or the unique `(8c/243)p^5` term of
`R` gives the contradiction.  Suppose both are nonconstant.  Let `d` be the
coefficient of `E_7` in (17).

If `d=0`, the upper Newton edge of `Phi` forces

```text
3a=2b,                 q_lead^2=p_lead^3/3.              (29)
```

If `d!=0`, the only slopes at which the maximum among

```text
3a+b, 3b, 4a, a+2b
```

is attained twice are `b/a=1` and `b/a=3/2`.  At slope one the `p^5` term of
`R` is uniquely highest.  At slope `3/2`, relation (29) holds and the top
part of `R` becomes

```text
-(52/243)p_lead^5 != 0.                                  (30)
```

Thus `R` cannot be affine, and degree eight is impossible.

### Degree ten

If one of `p,q` is constant, the unique `q^4` or `p^6` term in `R` gives the
contradiction after specialization.  If both are nonconstant, the two top
degrees in

```text
Phi_10=(70c/243)pq(-p^3+6q^2)+lower terms
```

force `3a=2b`, say `a=2s,b=3s`, and

```text
q_lead^2=p_lead^3/6.                                    (31)
```

On this branch the top part of `R_10` is

```text
-(595/8748)p_lead^6 != 0,
```

of degree `12s`.  Degree ten is impossible.

### Degree eleven

If exactly one of `p,q` is nonconstant, `p^6` or `q^4` uniquely dominates
`Phi_11`.  If both are nonconstant, its three top monomials force
`a=2s,b=3s`.  Put

```text
Y=q_lead^2/p_lead^3.
```

The leading first-integral equation is

```text
2-90Y+135Y^2=0.                                         (32)
```

The leading primitive could vanish only at `Y=2/15`, but substituting this
value into the left side of (32) gives `-38/5`, not zero.  Hence `R` has
nonzero degree `13s`, and degree eleven is impossible.

This proves every new exclusion.

## 7. Tame endpoints and honest frontier

The case `mu=0` is impossible directly: then the reduced mate is `S(x)` and
the coefficient of `y^2` in `J(P,S)=-P_yS'` forces `S'=0`.  Positive multiples
of three were removed in Section 2.  Therefore the only reduced degrees at
most eleven not excluded above are one and two.

At degree one, THM-2063 gives a tame inverse.  At degree two, apply THM-2071
to the Keller pair `(Q,-P)` and then restore the target shear.  This proves the
statement in Section 1.

The proof does **not** close the cubic-fiber stratum.  Its exact next object is
the all-degree Faber ladder (11)--(14), and the first surviving reduced degree
is thirteen.  For example,

```text
Phi_13=(65/6561)p(p^6-63p^3q^2+189q^4),
R_13=(91/6561)q(5p^6-60p^3q^2+27q^4).                 (33)
```

An all-degree closure would follow from two precise noncollision laws:

1. at every finite centering pole, the boundary polynomial
   `E_n(1;a,-1-a)` and the flux polynomial `Phi_n(a,-1-a)` have no common
   root;
2. on every upper Newton branch of a constant linear combination of the
   `Phi_m`, the matching linear combination of the `R_m` cannot be affine.

These are now the irreducible cubic questions.  The quotient to `Phi` alone
forgets the primitive `R`; the pair `(Phi,R)` is the necessary sidecar.  No
generic cover-degree, Jelonek-curve, VC(4), or full JC(2) equivalence is being
claimed.

## 8. Exact referee

The companion script checks over the exact rational polynomial ring:

- the eight identities `L(E_m)=z Phi_m'+R_m'` in (13)--(16);
- the rank-two connection (9);
- all six boundary resultants in (22);
- the nonzero leading primitives in (28)--(32);
- the degree-thirteen frontier covariants in (33);
- and the tame control `P=y^3+x, Q=y` with Jacobian one.

The normal run ends in `RESULT: PASS`. The current referee implements its
checks with Python `assert`, so an optimized `python -O` transcript is vacuous
and must not be cited as independent evidence. The computation checks
identities and hostile boundary controls; the proof above supplies the
all-polynomial quantifiers. QED.
