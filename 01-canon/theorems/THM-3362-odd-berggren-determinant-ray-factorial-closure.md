---
id: THM-3362
title: "Ordered-simplex odd-profile tensor pairs satisfy first-window FC(3)"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the Berggren
  coordinates T=x+y+z, U=x-y-z, V=x+y-z, homogenize any two nonzero real odd
  polynomials h,k along U/T,V/T, multiply them, and allow an arbitrary
  nonnegative radial power T^e.  Their factorial moments factor exactly as a
  radial factorial times the product of the two interval moments.  The first
  three factorial moments then detect every A+B R+C R^2.  Powers P^d of the
  THM-3357 transfer determinant are the specialization h=k=t^d,e=d for
  positive odd d.  Two moments are sharply insufficient; an exact complex
  odd profile defeats first-three detection, while a mixed-parity pair defeats
  the ordered-triangle factorization.
source: root/factorial-jacobian-lrc-threebranch-2026-08-14
depends_on:
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
  - THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit
related:
  - THM-3107-finite-layer-product-gamma-consecutive-width-three-orientation
  - THM-3125-monomial-ray-first-window-factorial-closure-in-three-variables
  - THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go
script: 04-computation/factorial_odd_berggren_determinant_ray_thm3362.py
output: 05-knowledge/results/factorial_odd_berggren_determinant_ray_thm3362.out
script_sha256: 9a2795d1bbf9932ff4aacbd7b19dc21474b549bd34ddfcee2415188765831e0b
output_sha256: ac16812f9083a97fefdb940852322a65e04c83872b154a90e6cde3af68ac231a
hash_basis: LF-normalized bytes
---

# THM-3362 -- ordered-simplex odd-profile tensor pairs satisfy first-window FC(3)

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3357 found a cubic determinant whose odd factorial moments all vanish.
The cubic is not isolated: its ordered-simplex coordinates separate any two
odd one-variable profiles.  A radial factorial gap, together with two
Cauchy--Schwarz inequalities, then orients the first three moments on every
such ray.  This gives an infinite nonlinear, generally nonhomogeneous FC(3)
subclass.  Exact witnesses show that complex coefficients defeat first-three
detection, while mixed parity defeats the factorization mechanism.

## 1. Inheritance and the family

Let

```text
calL:C[x,y,z] -> C,                  calL(x^i y^j z^k)=i!j!k!,      (1)

T=x+y+z,                 U=x-y-z,                 V=x+y-z.          (2)
```

The closest proved mechanisms are THM-3018's exponential/simplex integral,
THM-3125's positive product-Gamma theorem on one **monomial** ray, and
THM-3357's exact law for the transfer determinant

```text
P=TUV=(x-y-z)(x+y-z)(x+y+z).                              (3)
```

The corrected near miss is MISTAKE-350: full `FC(3)` is not homogeneous-only,
and its index is ambient variable count, not monomial-slot count.  The
least-used sidecar taken up here is the ordered-simplex pair `(U/T,V/T)`.
The canonical controls in Section 5 show that an arbitrary complex odd profile
can defeat first-three detection, while a mixed-parity pair loses the
ordered-triangle product formula.

Fix two nonzero odd polynomials

```text
h,k in R[tau],       s=deg h, t=deg k  (so s,t are odd and >=1),    (4)
```

and homogenize them by

```text
H_h(A,T)=T^s h(A/T),       H_k(A,T)=T^t k(A/T) in R[A,T].   (5)
```

For any integer `e>=0`, define

```text
R_(h,k,e)=T^e H_h(U,T)H_k(V,T),             D=s+t+e>=2.     (6)
```

The three linear forms `(U,V,T)` are independent, so `R_(h,k,e)` is a nonzero
homogeneous polynomial of degree `D`.

**Theorem.**  For every family member (6) and all `A,B,C in C`, put

```text
F=A+B R_(h,k,e)+C R_(h,k,e)^2.                              (7)
```

Then

```text
calL(F)=calL(F^2)=calL(F^3)=0       ==>       A=B=C=0.        (8)
```

Consequently each polynomial in (7) satisfies the conclusion of the original
ambient three-variable Factorial Conjecture: if all its positive factorial
moments vanish, then it is zero.  The same family results from replacing
`(h,k)` by `(c_1 h,c_2 k)` for any `c_1,c_2 in C^*`, since
`R_(c_1h,c_2k,e)=c_1c_2 R_(h,k,e)` and `B,C` absorb the resulting scalar
powers.

## 2. Exact ordered-simplex moment compiler

For `r>=0`, the exponential-integral identity gives

```text
calL(R_(h,k,e)^r)
 =integral_([0,infinity)^3) R_(h,k,e)^r exp(-T) dx dy dz.     (9)
```

Write `x=TX`, `y=TY`, `z=TZ` on the standard simplex
`X+Y+Z=1`.  Then `dx dy dz=T^2 dT dX dZ`, while

```text
u=U/T=2X-1,                 v=V/T=1-2Z,                       (10)
X=(u+1)/2,        Z=(1-v)/2,        Y=(v-u)/2.                (11)
```

Thus the simplex becomes the ordered triangle

```text
-1<=u<=v<=1,                    |d(X,Z)/d(u,v)|=1/4,          (12)
```

and

```text
R_(h,k,e)=T^D h(u)k(v).                                      (13)
```

The factors need not be equal.  The reflection

```text
rho(u,v)=(-v,-u)                                             (13a)
```

preserves the ordered triangle.  Because `h,k` are both odd, it sends
`h(u)^r k(v)^r` to `h(v)^r k(u)^r`.  Thus the ordered integrals with the two
profiles exchanged are equal.  Their sum is the integral over the full square,
so each is half the product of the interval integrals.  Separating the radial
Gamma integral proves the exact all-`r` formula

```text
calL(R_(h,k,e)^r)
 =(Dr+2)!/8
   * [integral_(-1)^1 h(tau)^r dtau]
   * [integral_(-1)^1 k(tau)^r dtau].                       (14)
```

At `r=0`, the right side is `2!/8 * 2 * 2=1`, fixing the normalization.
Because both profiles are odd,

```text
calL(R_(h,k,e)^r)=0                       for every odd r.    (15)
```

This is a tensor-product compiler, not a symmetry slogan: it supplies every
even moment as well as the cancellations.

## 3. A three-moment orientation lemma

The remaining algebra is independent of the origin of the moment sequence.
Let a unital complex linear functional on `C[Q]` satisfy

```text
L(Q)=L(Q^3)=L(Q^5)=0,
p=L(Q^2)>0,          q=L(Q^4)>p^2,          w=L(Q^6),          (16)
p w>3q^2.                                                       (17)
```

For `G=A+BQ+CQ^2`, vanishing of the first moment gives `A=-pC`.
After this substitution, direct expansion gives

```text
L(G^2)=pB^2+(q-p^2)C^2,                                     (18)
L(G^3)=3(q-p^2)B^2C+(w-3pq+2p^3)C^3.                       (19)
```

If `C=0`, (18) and the first moment give `B=A=0`.  If `C!=0`, eliminate
`B^2` from (18).  Equation (19) becomes

```text
L(G^3)=C^3/p * J,
J=pw+3p^2q-p^4-3q^2
 =[pw-3q^2]+p^2(3q-p^2)>0.                                  (20)
```

Thus the first three moments have no nonzero common projective zero.  The
argument is over `C`; real positivity is required only for the scalar
orientation in (16)--(17).

## 4. The radial factorial gap proves the hypotheses

Put

```text
I_j(h)=integral_(-1)^1 h(tau)^j dtau,
I_j(k)=integral_(-1)^1 k(tau)^j dtau.                       (21)
```

Since both profiles are real and nonzero, their interval moments of orders
`2,4,6` are positive.  Formula (14) gives

```text
p=(2D+2)! I_2(h)I_2(k)/8,
q=(4D+2)! I_4(h)I_4(k)/8,
w=(6D+2)! I_6(h)I_6(k)/8.                                  (22)
```

Strict variance under the positive exponential measure gives

```text
q-p^2=Var(R_(h,k,e)^2)>0.                                  (23)
```

Cauchy--Schwarz applied separately to `h` and `k` gives

```text
I_4(h)^2<=I_2(h)I_6(h),          I_4(k)^2<=I_2(k)I_6(k),    (23a)
```

while

```text
 p w/q^2
 =[(2D+2)!(6D+2)!/(4D+2)!^2]
   * [I_2(h)I_6(h)/I_4(h)^2] * [I_2(k)I_6(k)/I_4(k)^2]
 >=product_(j=1)^(2D) [(4D+2+j)/(2D+2+j)].                 (24)
```

Each product factor is at least its value at `j=2D`:

```text
(3D+1)/(2D+1)>=7/5                         because D>=2.    (25)
```

Consequently

```text
p w/q^2 >=(7/5)^(2D)>=(7/5)^4=2401/625>3.                 (26)
```

Equations (15), (23), and (26) verify the lemma and prove (8).

## 5. Sharp controls and equality boundaries

### The third moment is necessary

For every real odd-profile-pair generator, set

```text
C=1,                    A=-p,                    B^2=-(q-p^2)/p. (27)
```

A square root exists in `C`, and (18) shows

```text
calL(F)=calL(F^2)=0,                    F!=0.             (28)
```

The third moment is `J/p>0`.  Thus two moments never detect the entire
three-slot subspace; the window length three in (8) is sharp.

### Unequal odd profiles are intrinsic

The factorization is genuinely broader than a tensor square.  Taking

```text
h(tau)=tau,              k(tau)=tau^3,              e=0
```

gives `R_(h,k,0)=UV^3`, of degree four, and (14) yields the exact prefix

```text
[calL(R_(h,k,0)^r)]_(r=0)^4
  =[1,0,86400,0,49249028505600].                           (28a)
```

### The transfer-determinant rays

For every positive odd `d`, take

```text
h(tau)=k(tau)=tau^d,             s=t=d,             e=d.    (29)
```

Then

```text
R_(h,k,e)=T^dU^dV^d=P^d,                                   (30)
```

so (8) specializes to

```text
calL(A+B P^d+C P^(2d))
=calL((A+B P^d+C P^(2d))^2)
=calL((A+B P^d+C P^(2d))^3)=0
                       ==> A=B=C=0.                           (31)
```

This recovers the originally discovered odd determinant-ray theorem and
strictly extends it.  For example, `h=k=tau+tau^3,e=0` gives the degree-six
generator

```text
(UT^2+U^3)(VT^2+V^3),                                      (32)
```

which is not a power of `P`.

### Arbitrary complex odd profiles fail

The real-up-to-phase condition cannot simply be erased.  Let

```text
h_a(tau)=tau+a tau^3,            15a^2+42a+35=0,
a=-7/5 +- (2i sqrt(21))/15.                                    (33)
```

With `h=k=h_a` and `e=0`, put

```text
R_a=(UT^2+aU^3)(VT^2+aV^3).                                  (34)
```

Oddness still kills moments one and three, while

```text
I_2(h_a)=2(15a^2+42a+35)/105=0.                              (35)
```

Hence the nonzero polynomial `R_a` satisfies

```text
calL(R_a)=calL(R_a^2)=calL(R_a^3)=0.                         (36)
```

It is not an FC counterexample.  Modulo the quadratic in (33),

```text
I_4(h_a)=128(417a+560)/1126125!=0,                           (37)
```

because the displayed quadratic and `417a+560` are coprime.  Therefore

```text
calL(R_a^4)=26!/8 * I_4(h_a)^2!=0.                           (38)
```

This is an exact failure witness for the proposed arbitrary-complex
three-moment extension; it does not classify every complex profile that may
still satisfy (8).

### Mixed parity destroys the reflection factorization

Common odd parity, rather than equality of the factors, is load-bearing.  If
`h(tau)=1`, `k(tau)=tau`, and `e=0`, then `R=V` and

```text
calL(V)=1.                                                 (38a)
```

Indeed, the ordered-triangle angular integral is `2/3` before the Jacobian
factor `1/4`, and the radial factor is `3!`.  The false half-square product
would instead contain `integral_(-1)^1 tau dtau=0`.  Thus a mixed-parity pair
does not satisfy (14).

No parity restriction is imposed on `e`.  It changes only the radial degree
`D`; angular oddness supplies (15).

## 6. Connection contract and frontier boundary

The exact transfer is

```text
source:     two real odd interval profiles h,k and radial exponent e
target:     an ambient FC(3) subspace span{1,R_(h,k,e),R_(h,k,e)^2}
map:        ordered-simplex homogenization (5)--(6)
preserves:  the product moment-zero predicate at every r and computes every
            target moment by (14)
sidecar:    h,k, their degrees s,t, common odd parity, and radial exponent e
loses:      the separate interval moments (retaining only their product),
            mixed-parity/nonseparable angular coupling, branch ancestry,
            and off-ray terms
hostiles:   the two-moment point (27), complex profile (33)--(38), and the
            mixed-parity factorization hostile (38a).                  (39)
```

The target moment sequence also forgets the labelled order of `h,k`: exchanging
the profiles generally changes the polynomial but not formula (14).

Although each `R_(h,k,e)` is homogeneous, a general polynomial (7) mixes total
degrees `0,D,2D`; the result is ambient FC(3), not merely HFC(3).  It is also
outside THM-3125's monomial rays.  The interval-moment product in (14), rather
than a product-Gamma factorization, is the new mechanism.

Via THM-3300, substitution

```text
(x,y,z)=(|z_1|^2,|z_2|^2,|z_3|^2)                            (40)
```

reproduces the same first-three detection statement in the
`U(1)^3`-invariant Gaussian subalgebra.  This is a moment identity, not
`FC=GMC`, and it supplies no planar-JC implication.  The determinant pencils
are sharply separate: if `q=det(xL+yM+zR)` is THM-3357's two-dimensional
spinor pencil, then at `(x,y,z)=(1,1,1)`

```text
det Sym^2(L+M+R)=q(1,1,1)^3=1,           P(1,1,1)=-3.       (40a)
```

Thus symmetric square does not commute with weighted summation, and the FC
cubic cannot be identified with the Hessian/Keller quadratic.  Likewise, the
ordered-simplex formula retains no LRC phase, owner, or safe-time exit;
THM-3357's determinant-gate Horn rule is a separate sufficient-gate mechanism.

## 7. Exact verification

The companion performs six exact checks.

1. Direct sparse factorial expansion independently verifies (14) for the
   equal pairs `h=k=tau,tau+tau^3,tau-2tau^3` and the unequal pair
   `h=tau,k=tau^3`, for `e=0,1,2` and `0<=r<=4`.
2. Direct integration and factorial expansion verify the mixed-parity hostile
   `h=1,k=tau` in (38a).
3. Symbolic elimination verifies (18)--(20).
4. Integer/rational arithmetic audits (24)--(26) for every `2<=D<=100` and
   the determinant-ray specialization for every positive odd `d<=99` and
   `0<=r<=6`.
5. Direct arithmetic in `Q[a]/(15a^2+42a+35)` verifies (35)--(38), including
   the nonzero fourth moment.
6. Exact matrix determinants verify the separate-pencil hostile (40a).

The finite ranges audit implementation and indexing; the universal theorem
comes from the displayed identities and inequalities.  Optimized and ordinary
runs must byte-match the stored transcript:

```text
python3 04-computation/factorial_odd_berggren_determinant_ray_thm3362.py
python3 -O 04-computation/factorial_odd_berggren_determinant_ray_thm3362.py
```

**End of proof.**
