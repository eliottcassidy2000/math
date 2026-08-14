---
id: THM-3362
title: "Ordered-simplex odd-profile tensor rays satisfy first-window FC(3)"
status: >
  PROVED + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.  In the Berggren
  coordinates T=x+y+z, U=x-y-z, V=x+y-z, homogenize any nonzero real odd
  polynomial h along U/T and V/T, multiply the two copies, and allow an
  arbitrary nonnegative radial power T^e.  The factorial moments of the
  resulting generator R factor exactly as a radial factorial times the square
  of one interval moment.  The first three factorial moments then detect every
  A+B R+C R^2.  Powers P^d of the THM-3357 transfer determinant are the
  specialization h=t^d,e=d for positive odd d.  Two moments are sharply
  insufficient, and an exact complex odd profile shows that real-up-to-phase
  is load-bearing for the three-moment theorem.
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
script_sha256: b656e386e673a65151d9f41c7e77a3ba06a158d64849816d72da399ed3c7f52d
output_sha256: 2e51a57d1c306c61a50dec9362f2b4036069cded904664e268d8f5e2195bef3a
hash_basis: LF-normalized bytes
---

# THM-3362 -- ordered-simplex odd-profile tensor rays satisfy first-window FC(3)

**PROVED + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-3357 found a cubic determinant whose odd factorial moments all vanish.
The cubic is not isolated: its ordered-simplex coordinates separate two copies
of any odd one-variable profile.  A radial factorial gap, together with one
Cauchy--Schwarz inequality, then orients the first three moments on every such
ray.  This gives an infinite nonlinear, generally nonhomogeneous FC(3)
subclass and an exact boundary where complex angular coefficients defeat the
three-moment strengthening.

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
The canonical hostile in Section 5 shows that arbitrary complex odd profiles
do not inherit the real orientation.

Fix a nonzero odd polynomial

```text
h in R[t],                    s=deg h  (so s is odd and s>=1),       (4)
```

and homogenize it by

```text
H_h(A,T)=T^s h(A/T) in R[A,T].                              (5)
```

For any integer `e>=0`, define

```text
R_(h,e)=T^e H_h(U,T)H_h(V,T),              D=2s+e>=2.       (6)
```

The three linear forms `(U,V,T)` are independent, so `R_(h,e)` is a nonzero
homogeneous polynomial of degree `D`.

**Theorem.**  For every family member (6) and all `A,B,C in C`, put

```text
F=A+B R_(h,e)+C R_(h,e)^2.                                  (7)
```

Then

```text
calL(F)=calL(F^2)=calL(F^3)=0       ==>       A=B=C=0.        (8)
```

Consequently each polynomial in (7) satisfies the conclusion of the original
ambient three-variable Factorial Conjecture: if all its positive factorial
moments vanish, then it is zero.  The same family results from replacing `h`
by `c h` for any `c in C^*`, since `R_(ch,e)=c^2R_(h,e)` and `B,C` absorb the
two scalar powers.

## 2. Exact ordered-simplex moment compiler

For `r>=0`, the exponential-integral identity gives

```text
calL(R_(h,e)^r)
 =integral_([0,infinity)^3) R_(h,e)^r exp(-T) dx dy dz.       (9)
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
R_(h,e)=T^D h(u)h(v).                                         (13)
```

The angular integrand `h(u)^r h(v)^r` is symmetric in `u,v`, so the ordered
triangle contributes half the full square.  Separating the radial Gamma
integral therefore proves the exact all-`r` formula

```text
calL(R_(h,e)^r)
 =(Dr+2)!/8 * [ integral_(-1)^1 h(t)^r dt ]^2.                (14)
```

At `r=0`, the right side is `2!/8 * 2^2=1`, fixing the normalization.
Because `h` is odd,

```text
calL(R_(h,e)^r)=0                         for every odd r.    (15)
```

This is a tensor-square compiler, not a symmetry slogan: it supplies every
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
I_j=integral_(-1)^1 h(t)^j dt.                              (21)
```

Since `h` is real and nonzero, `I_2,I_4,I_6` are positive.  Formula (14)
gives

```text
p=(2D+2)! I_2^2/8,
q=(4D+2)! I_4^2/8,
w=(6D+2)! I_6^2/8.                                         (22)
```

Strict variance under the positive exponential measure gives

```text
q-p^2=Var(R_(h,e)^2)>0.                                    (23)
```

Cauchy--Schwarz gives `I_4^2<=I_2I_6`, while

```text
 p w/q^2
 =[(2D+2)!(6D+2)!/(4D+2)!^2] * (I_2I_6/I_4^2)^2
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

For every real-profile generator, set

```text
C=1,                    A=-p,                    B^2=-(q-p^2)/p. (27)
```

A square root exists in `C`, and (18) shows

```text
calL(F)=calL(F^2)=0,                    F!=0.             (28)
```

The third moment is `J/p>0`.  Thus two moments never detect the entire
three-dimensional ray; the window length three in (8) is sharp.

### The transfer-determinant rays

For every positive odd `d`, take

```text
h(t)=t^d,                    s=d,                    e=d.    (29)
```

Then

```text
R_(h,e)=T^dU^dV^d=P^d,                                      (30)
```

so (8) specializes to

```text
calL(A+B P^d+C P^(2d))
=calL((A+B P^d+C P^(2d))^2)
=calL((A+B P^d+C P^(2d))^3)=0
                       ==> A=B=C=0.                           (31)
```

This recovers the originally discovered odd determinant-ray theorem and
strictly extends it.  For example, `h=t+t^3,e=0` gives the degree-six
generator

```text
(UT^2+U^3)(VT^2+V^3),                                      (32)
```

which is not a power of `P`.

### Arbitrary complex odd profiles fail

The real-up-to-phase condition cannot simply be erased.  Let

```text
h_a(t)=t+a t^3,                  15a^2+42a+35=0,
a=-7/5 +- (2i sqrt(21))/15.                                    (33)
```

With `e=0`, put

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

This is the exact failure boundary for the proposed arbitrary-complex
three-moment extension; it does not classify every complex profile that may
still satisfy (8).

No parity restriction is imposed on `e`.  It changes only the radial degree
`D`; angular oddness supplies (15).

## 6. Connection contract and frontier boundary

The exact transfer is

```text
source:     a real odd interval profile h and radial exponent e
target:     an ambient FC(3) subspace span{1,R_(h,e),R_(h,e)^2}
map:        ordered-simplex homogenization (5)--(6)
preserves:  every factorial moment through the exact tensor square (14)
sidecar:    the real profile h, its degree s, and radial exponent e
loses:      arbitrary angular coupling, branch ancestry, and off-ray terms
hostiles:   the two-moment point (27) and complex profile (33)--(38).      (39)
```

Although each `R_(h,e)` is homogeneous, a general polynomial (7) mixes total
degrees `0,D,2D`; the result is ambient FC(3), not merely HFC(3).  It is also
outside THM-3125's monomial rays.  The product in (14), rather than a
product-Gamma factorization, is the new mechanism.

Via THM-3300, substitution

```text
(x,y,z)=(|z_1|^2,|z_2|^2,|z_3|^2)                            (40)
```

reproduces the same first-three detection statement in the
`U(1)^3`-invariant Gaussian subalgebra.  This is a moment identity, not
`FC=GMC`, and it supplies no planar-JC implication.  Likewise, the ordered
simplex formula retains no LRC phase, owner, or safe-time exit; THM-3357's
determinant-gate Horn rule is a separate sufficient-gate mechanism.

## 7. Exact verification

The companion performs four exact checks.

1. Direct sparse factorial expansion independently verifies (14) for
   `h=t,t+t^3,t-2t^3`, `e=0,1,2`, and `0<=r<=4`.
2. Symbolic elimination verifies (18)--(20).
3. Integer/rational arithmetic audits (24)--(26) for every `2<=D<=100` and
   the determinant-ray specialization for every positive odd `d<=99`.
4. Direct arithmetic in `Q[a]/(15a^2+42a+35)` verifies (35)--(38), including
   the nonzero fourth moment.

The finite ranges audit implementation and indexing; the universal theorem
comes from the displayed identities and inequalities.  Optimized and ordinary
runs must byte-match the stored transcript:

```text
python3 04-computation/factorial_odd_berggren_determinant_ray_thm3362.py
python3 -O 04-computation/factorial_odd_berggren_determinant_ray_thm3362.py
```

**End of proof.**
