---
id: THM-3375
title: "Berggren U-spine: an infinite positive two-distinct-cube Pell ray"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED WITH SCOPE REPAIR
source: kps-s181
depends_on:
  - THM-3370
  - THM-3334
  - THM-463
related:
  - THM-3359
  - THM-3376
companion: 04-computation/berggren_positive_two_cube_pell_ray_kps_s181.py
output: 05-knowledge/results/berggren_positive_two_cube_pell_ray_kps_s181.out
script_sha256: f01b96a8db812da88bbec974f58f17a5d1b794cdd5168617892ae2457a6f81d0
output_sha256: ec04036ad35e463a04210179bed61a54f280e54db4f0642ffc5e315f9f888332
hash_basis: LF-normalized bytes
audit: independent Pell/positivity reconstruction and ambient-versus-orbit local-solubility repair
---

# THM-3375 -- infinitely many positive two-cube points on the U-spine

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED WITH SCOPE REPAIR.**

This closes the positive-infinitude question left open by THM-3370.  The key
move is to let the cube sum `d=x+y` move on a Pell section.  Every fixed-`d`
Pell orbit has only finitely many positive points.

## 1. Explicit family

Put

```text
epsilon=7775+312 sqrt(621),              N(epsilon)=1.  (1)
```

For `j>=0`, define positive integers `W_j,h_j` by

```text
W_j+h_j sqrt(621)
  =(5059+203 sqrt(621)) epsilon^(5j).                   (2)
```

The seed and the fifth-power congruence satisfy

```text
5059^2-621*203^2=2692,
epsilon^5=A+B sqrt(621),
(A,B)=(-1,0) mod 29.                                  (3)
```

Since `29|203`, recurrence under `(A,B)` proves `29|h_j` for every `j`.
Write `u_j=h_j/29` and set

```text
d_j=841u_j^2+2,
e_j=u_j W_j,
a_j=21866u_j^3+81u_j,
q_j=568516u_j^4+2860u_j^2+1,                          (4)
x_j=(d_j+e_j)/2,
y_j=(d_j-e_j)/2,
r_j=(a_j-1)/2.                                        (5)
```

Then all quantities in `(4)--(5)` are integers and

```text
x_j>y_j>0,
x_j^3+y_j^3=a_j^2+2=(2r_j+1)^2+2=Q_(r_j).             (6)
```

The values `r_j` are strictly increasing.  Consequently infinitely many
Berggren U-spine scalars are sums of two distinct positive cubes.

The first member of this particular ray is

```text
2899^3+38312^3
  =7500605^2+2
  =Q_3750302
  =56259075366027.                                     (7)
```

It is not the global first intersection; THM-3370's smaller point remains
`107^3+232^3=Q_1851`.

## 2. Pell and polynomial identities

Taking norms in `(2)` and using `h_j=29u_j` gives

```text
W_j^2-522261u_j^2=2692,
522261=621*29^2.                                      (8)
```

Direct coefficient comparison in `u` proves

```text
a_j^2+2=d_j q_j,
4q_j-d_j^2=3e_j^2.                                   (9)
```

For example, the two nontrivial coefficients in the second identity are

```text
4*568516-841^2=3*522261,
4*2860-4*841=3*2692.                                  (10)
```

Since

```text
x^3+y^3=(x+y)((x+y)^2+3(x-y)^2)/4,                    (11)
```

equations `(9)` imply `(6)`.

The parity is also stable.  The scalar coefficient of `epsilon^5` is odd and
its radical coefficient is even; starting from odd `W_0,h_0`, every `W_j,h_j`
is odd.  Hence `u_j,d_j,e_j,a_j` are odd, so `(5)` is integral.

## 3. The orbit stays strictly inside the positive chamber

Let `rho_j=W_j/h_j`.  If `epsilon^5=A+B sqrt(621)`, multiplication gives

```text
rho_(j+1)=(A rho_j+621B)/(B rho_j+A).                  (12)
```

For `rho>sqrt(621)`,

```text
rho'-sqrt(621)
 =(A-B sqrt(621))(rho-sqrt(621))/(B rho+A)>0,
rho'-rho=B(621-rho^2)/(B rho+A)<0.                    (13)
```

The seed norm is positive, and therefore

```text
sqrt(621)<rho_j<=5059/203<29.                         (14)
```

Since `h_j=29u_j`, equations `(4)` and `(14)` give

```text
0<e_j=h_j W_j/29<h_j^2<d_j.                          (15)
```

Thus `x_j>y_j>0`.  Moreover

```text
e_j/d_j -> sqrt(621)/29 < 1,                          (16)
```

so this is a genuine interior ray, not a family tending to the boundary
`y=0`.  Both coefficients in `(1)` are positive and `epsilon>1`, so `h_j`
grows exponentially.  Since `r_j` is cubic in `u_j`, the family contributes
`Theta(log R)` members with `r_j<=R`.

## 4. Why the previous fixed-sum Pell orbit could not work

For any positive pair put `d=x+y` and `e=x-y`, choosing `e>0` after ordering.
Positivity and distinctness are exactly

```text
0<e<d,                    e=d mod 2.                   (17)
```

For fixed `d`, there are at most `floor((d-1)/2)` possible values of `e`.
Equivalently,

```text
d^3/4 < x^3+y^3 < d^3.                                (18)
```

Thus no fixed-`d` Pell orbit can prove positive infinitude: it eventually
leaves `(17)`.  The moving `d_j=h_j^2+2` in `(4)` is the missing coordinate.

## 5. Discriminant-minus-eight form compiler

Every positive solution in this intersection has odd `d=x+y` and

```text
dq-a^2=2.                                             (19)
```

Hence the positive binary quadratic form `[d,2a,q]` is primitive and has
discriminant `-8`.  A reduced primitive positive form of this discriminant has
`3A^2<=8`, hence `A=1`; reduction then leaves only `[1,0,2]`.  Thus the class
number is one, so there are integers `p,r,s,t` with `pt-rs=+/-1` such that

```text
d=p^2+2r^2,
a=ps+2rt,
q=s^2+2t^2.                                           (20)
```

This is a structural compiler for the norm-factor part of the intersection;
the additional square equation in `(9)` remains load-bearing.

The particular Pell ray comes from a more general reduced-seed ansatz.  Put

```text
d=h^2+2,          a=kd+h,
h=nu,             k=mu,             e=uW.             (21)
```

Then `(19)` is automatic with `q=k^2d+2kh+1`, and the cube-square equation is
equivalent to

```text
3W^2=n^2(4m^2-n^2)u^2+4(2m^2+2mn-n^2).               (22)
```

The successful slope is `(m,n)=(26,29)`, for which `(22)` is exactly `(8)`.
Other slopes with `n/2<m<n`, an integral seed, a divisibility-preserving unit
subsequence, and a ratio below `n` are crisp candidates for independent
positive rays.

## 6. The ambient mod-63 sieve is the complete local obstruction

THM-3370 proved the necessary depth classes

```text
r mod 63 in
{2,4,10,11,16,24,25,31,37,38,46,51,52,58,60}.         (23)
```

They are also sufficient for local solubility of the **ambient U-spine
equation**.  In fact,

```text
X^3+Y^3=Q_r is soluble in every Z_p
iff r mod 63 lies in (23).                             (24)
```

Necessity is exactly the cube-residue obstruction modulo `9` and `7` from
THM-3370.  For sufficiency, fix an integer `r` in `(23)` and a prime `p`.
The equation

```text
X^3+Y^3=Q_r                                             (25)
```

has a solution in `Z_p`.

- At `p=3`, the allowed depths give `Q_r=0` or `2 mod 9`.  If
  `Q_r=2 mod 9`, then `Q_r/2` is a 3-adic cube; take `X=Y`.  If `Q_r=9M`,
  put `X=1,Y=-1+3z`.  Then

  ```text
  X^3+Y^3=9(z-3z^2+3z^3),
  ```

  and the map on the right divided by `9` has derivative `1 mod 3`, so
  Hensel solves it for every `M`.  This explicit mod-`9` argument is needed;
  the original cubic derivatives are divisible by `3`.
- At `p=7`, `(23)` is exactly the cube-sum support, and each supported residue
  has a representative with at least one nonzero coordinate, so Hensel lifts.
- If `p=2`, choose `(1,0)` for an odd target and `(1,1)` for an even target;
  a partial derivative is odd.
- For odd `p=2 mod 3`, cubing is a bijection on `F_p^*`; the zero target uses
  `(1,-1)`.
- For `p=1 mod 3`, `p>=13`, the smooth projective cubic for a nonzero target
  has at least `p+1-2sqrt(p)>3` points by Hasse, while only three lie at
  infinity.  A zero target again uses `(1,-1)`.  Every affine point obtained
  is nonsingular and lifts.

Therefore no congruence sieve based only on ambient cube-sum representability
can improve the exact periodic ceiling

```text
15/63=5/21.                                            (26)
```

The remaining sparsity is global divisor/square geometry, not another local
modulus.

This ambient assertion is not an equidistribution statement about the single
Pell ray `(2)`.  Exact modular dynamics give that ray period `24` modulo `63`
and only

```text
r_j mod 63 in {38,51,60},             r_j=2 mod 5.     (27)
```

Also `epsilon^5=(1,0) mod 11` in scalar/radical coordinates.  These are
orbit-specific scars.  They restrict this distinguished ray even though all
residue classes modulo `5` are locally admissible for the ambient cube-sum
equation.  Quotienting the ambient surface to one Pell orbit loses local
residue coverage.

## 7. Exact audit and scope

The standard-library companion verifies the two norms, the fifth-power
congruence, every coefficient identity, the first eight enormous orbit points,
their positivity, parity, discriminant, primitivity, strict growth, the first
member `(7)`, and the exact mod-`7`, mod-`9`, and CRT supports.  Normal and
optimized outputs are byte-identical.  LF-normalized source/output hashes are

```text
f01b96a8db812da88bbec974f58f17a5d1b794cdd5168617892ae2457a6f81d0
ec04036ad35e463a04210179bed61a54f280e54db4f0642ffc5e315f9f888332.
```

An independent audit reconstructed the Pell family and positivity identities,
then caught and repaired the possible misreading of “local saturation.”  It
proved the ambient iff `(24)`, checked the exceptional-prime lifts explicitly,
and verified hostile controls for every prime below `500` and through modulus
`3^8`.  It separately confirmed that the orbit scars `(27)` are real but do
not narrow the ambient local-solubility set.

This theorem proves positive infinitude but no asymptotic for all positive
intersections.  It does not prove the proposed fixed-`d` uniqueness, nor
connect the Eisenstein factor labels to the coprime Gaussian carrier of
THM-3370.  THM-3376 subsequently proves a second distinct positive ray and a
finite slope atlas through denominator `29`; neither result is a classification
of all intersections.  This theorem supplies no LRC owner/phase, factorial-moment cancellation,
Keller mate, or AMM flow; LRC(14), FC(3), JC(2), and AMM 12592 remain open.
