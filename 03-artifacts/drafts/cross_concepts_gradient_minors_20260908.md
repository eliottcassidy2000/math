# Periodic gradient minors: exact matched-energy shear-mixture control

**Status: PROVED elementary/Fourier identities + FINITE-EXACT independent
audit.** This note audits the root's spatial compatibility addition to
THM-4459. It is a supporting proof/control, not a separate canon namespace.
No Euler trajectory, gradient Young-measure realization, or LRC implication
is asserted. The parent owns integration into the theorem.

## 1. The minor identity survives exactly at W^(1,2)

Let M be a constant real n-by-n matrix, let v belong to W^(1,2)(T^n;R^n),
and write A=M+Dv, with (Dv)_(ij)=partial_j v_i. The normalized torus has
total measure one. Every two-by-two minor satisfies

```
mean(A_ij*A_kl-A_il*A_kj) = M_ij*M_kl-M_il*M_kj.       (1)
```

Indeed, the terms linear in Dv have mean zero. For smooth periodic v,
periodic integration by parts gives

```
integral (partial_j v_i)(partial_l v_k)
 = -integral v_i*(partial_j partial_l v_k)
 = integral (partial_l v_i)(partial_j v_k).
```

This cancels the quadratic difference. Smooth periodic approximation in
W^(1,2) makes the product terms converge in L^1, so (1) holds at the stated
regularity. This proof asserts no integrability or weak continuity of
three-by-three determinants at W^(1,2).

If A has rank at most one almost everywhere, its minors vanish and (1)
forces rank M<=1. If in addition tr A=0 almost everywhere, then tr M=0.
Consequently support in

```
S={a tensor n : a,n in R^3 and a dot n=0}
```

forces M itself to belong to S, including its zero member. This strengthens
the rank-only conclusion when the pure-shear premise is retained.

The field A is periodic. The associated velocity u(x)=Mx+v(x) is
affine-periodic; it is not claimed periodic when M is nonzero. If tr M=0
and div v=0, then div u=0.

## 2. Quantitative compatibility and optimal universal constant

In three dimensions let d_S(A)=inf_(S in S)||A-S||_F. If the singular
values of A are s_1>=s_2>=s_3>=0, then

```
d_S(A)^2 >= s_2^2+s_3^2,
||cof A||_F^2 = s_1^2*s_2^2+s_1^2*s_3^2+s_2^2*s_3^2
 <= (s_1^2+s_2^2+s_3^2)*(s_2^2+s_3^2).
```

The first estimate is the distance to the larger set of all rank-at-most-one
matrices. For completeness, fixing a unit right vector n, the least-squares
choice of a in a tensor n is a=An. Minimizing gives ||A||_F^2-||An||^2;
the maximum of ||An||^2 is s_1^2. The cofactor formula follows by orthogonal
invariance of its Frobenius norm and evaluation on the singular-value
diagonal. Thus

```
||cof A||_F <= ||A||_F*d_S(A).                         (2)
```

Combining (1), the triangle inequality for the integral, and Cauchy--Schwarz
gives

```
||cof M||_F^2
 <= mean(||A||_F^2) * mean(d_S(A)^2).                  (3)
```

The estimate is meaningful for every W^(1,2) field above: 0 belongs to S,
so d_S(A)<=||A||_F and every expression is integrable.

**The universal constant 1 in (3) is optimal even at trace zero.** Take the
constant field M_epsilon=E_12+epsilon E_23, 0<epsilon<=1, and v=0. Its
singular values are 1,epsilon,0. The admissible shear E_12 realizes the
distance epsilon, so the unrestricted rank-one lower estimate is attained.
Exactly,

```
||cof M_epsilon||_F^2=epsilon^2,
||M_epsilon||_F^2=1+epsilon^2,
d_S(M_epsilon)^2=epsilon^2.
```

The ratio of the two sides before inserting the universal constant is
1/(1+epsilon^2), tending to 1. This proves optimality of the coefficient;
it does not assert equality at a nonzero epsilon.

## 3. Explicit divergence-free waves with exact quadratic cancellation

Use T^3=(R/(2*pi Z))^3 and normalized Haar measure. Let

```
M=E_12+E_23+E_31 = [[0,1,0],[0,0,1],[1,0,0]],

v_1= sin(x_1+x_2)+sin(x_1+x_3),
v_2=-sin(x_1+x_2)+sin(x_2+x_3),
v_3=-sin(x_2+x_3)-sin(x_1+x_3).
```

Each wave has amplitude perpendicular to its frequency:

| Frequency k | Amplitude a |
|---|---|
| (1,1,0) | (1,-1,0) |
| (0,1,1) | (0,1,-1) |
| (1,0,1) | (1,0,-1) |

Thus div v=0 pointwise, and A=M+Dv is trace-free. The exact mean data are

```
mean A=M,
mean(cof A)=cof M=M,
||cof M||_F^2=3,
mean ||A||_F^2=3+3*(4/2)=9.                           (4)
```

In particular (3) certifies

```
mean d_S(A)^2 >= 1/3.                                 (5)
```

The cancellation in (1) is nontrivial in the upper-left minor. If
c_ij=cos(x_i+x_j), then

```
A_11=c_12+c_13,       A_22=-c_12+c_23,
A_12=1+c_12,          A_21=-c_12.
mean(A_11*A_22)=-1/2=mean(A_12*A_21).
```

The zero minor is obtained by cancellation of two nonzero quadratic means.
It is not a consequence of treating the derivatives as independent.

## 4. An arbitrary shear mixture has the same mean and energy

Put equal mass 1/3 on

```
3E_12,       3E_23,       3E_31.                       (6)
```

Each is a pure shear; the mixture has mean M, second moment 9, and mean
squared distance to S equal to zero. Its average cofactor is zero, while
cof(mean)=M is nonzero. Thus (6) cannot be the law of a field M+Dv satisfying
the periodic W^(1,2) gradient constraint. Matching the mean and the second
moment does not restore the lost spatial compatibility.

The genuinely compatible positive control is M=E_12 and
v(x)=sin(x_2)e_1. Then A=(1+cos x_2)E_12 is pointwise a pure shear,
cof A is identically zero, and its exact energy is 3/2.

**Connection contract:** arbitrary shear law -> affine-periodic gradient
law, with the same matrix mean and Frobenius second moment. Those two
statistics are preserved in the control, but all spatial incidence and
the exact-gradient constraint were absent in the source law. The necessary
sidecar is the minor identity (1), arising from row-wise gradient
compatibility. The matched-energy pair (4)/(6) is the decisive falsification
of a first-moment-only realization transfer. The pointwise set S alone does
not encode that compatibility.

## 5. Exact evidence

The companion
[cross_concepts_gradient_minors_20260908.py](../../04-computation/cross_concepts_gradient_minors_20260908.py)
uses rational Gaussian Fourier coefficients and multivariate Laurent
convolution. It differentiates sine coefficients, reconstructs A, forms all
nine cofactors, and extracts the exact zero Fourier coefficient. An
independent route evaluates derivatives directly at the 64 points with
x_j in {0,pi/2,pi,3pi/2}; all trigonometric values are integers and the grid
is an exact cubature for the frequency support, not a numerical sample.

Four explicit fields, the three-shear mixture, and eight rational sharpness
parameters pass **444 explicit gates**. The coarse two-grid energy is 15,
not 9, providing an aliasing hostile against an insufficient quadrature.

Run:

```
python3 04-computation/cross_concepts_gradient_minors_20260908.py
python3 -O 04-computation/cross_concepts_gradient_minors_20260908.py
```

Both outputs match the stored LF
[transcript](../../05-knowledge/results/cross_concepts_gradient_minors_20260908.out).
The semantic SHA256 is
`43665a0fca1f35d44b31c0239ce2cfce3720ec7f223de6512fcd033648c5fe81`.
Checks use explicit exceptions and remain active under optimization.
