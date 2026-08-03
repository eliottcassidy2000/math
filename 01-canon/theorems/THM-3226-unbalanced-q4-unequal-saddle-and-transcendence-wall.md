---
id: THM-3226
title: "Unbalanced Q4 unequal saddle and transcendence wall"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The strong unbalanced four-vertex quotient has a unique unequal positive
  diagonal saddle (s,t,t,s), an exact implicit product radius, strict complex
  minimality, and a nondegenerate asymptotic
  C_d*r^(d-5/2)*R^(-r).  Its saddle parameter and scale ratio s/t are both
  transcendental.  Schanuel's conjecture makes the product radius
  transcendental, but unconditional radius transcendence and Q4
  non-P-recursiveness remain OPEN.  The analytic variational saddle theorem
  extends to every strongly connected fixed tournament quotient.
source: root/frontier-sidecars-cont-2026-08-02
audit: >
  A symbolic companion reconstructs the quotient determinant and positive
  pole numerator, the paired boundary, cubic-field norm identities, tangent
  Hessian eigenvalues, and one-block hostile.  Its 100-digit controls isolate
  the unequal saddle, compare the false equal pole, and test the d=1
  asymptotic ratio; those floating checks are NUMERICAL CONTROL ONLY.
  Normal, optimized, and stored transcripts agree byte-for-byte.
  An independent hostile audit rederived the analytic saddle and arithmetic
  proof, strengthened s/t from irrational to transcendental, repaired the raw
  non-C-finite attribution, and supplied the coercive compactness step for
  the quotient-wide theorem.  All corrected scopes were accepted.
depends_on:
  - THM-3213-tournament-normalized-cyclic-diagonal-and-fast-moving-jet-transform
  - THM-3121-path-cover-walk-content-substitution-kernel
related:
  - THM-3202-c3-repeated-join-moving-jet-formula-and-cfinite-obstruction
script: 04-computation/tournament_unbalanced_q4_saddle_thm3226.py
output: 05-knowledge/results/tournament_unbalanced_q4_saddle_thm3226.out
script_sha256: e2f47e09f912fd727955d362ce8f9ca797d8e18167676577e365443747c4a199
output_sha256: 11f660647f32a83695e1400575130544b1e829383163df59cc0610f4d99a1e9b
hash_basis: LF-normalized bytes
---

# THM-3226 -- unbalanced Q4 unequal saddle and transcendence wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Object and exact reduction

Let `Q4` have arcs

```text
1->2->3->1,  3->4,  4->1,  4->2,
```

and let `A` be its adjacency matrix.  Put `u(x)=exp(x)-1` and

```text
W_Q(X)=1^T X(I-A X)^(-1)1.
```

Exact elimination gives

```text
H(X):=det(I-A X)
 =1-X1 X2 X3-X2 X3 X4-X1 X2 X3 X4.                 (1)
```

Moreover `W_Q=N/H`, where `N` has fifteen monomials, all with positive
integer coefficients.  Hence the positive Perron pole is not cancelled, and
`W_Q^d` has pole order exactly `d` for every fixed `d>=1`.

For `V_r=Q4[T_r,T_r,T_r,T_r]`, THM-3213 gives

```text
a_d(r):=F_(V_r)(d)/(r!)^4
 =[x1^r x2^r x3^r x4^r] W_Q(u(x1),...,u(x4))^d.    (2)
```

On the positive singular boundary, (1) becomes

```text
u(x2)u(x3)(exp(x1+x4)-1)=1.                         (3)
```

The diagonal saddle maximizes `x1 x2 x3 x4` subject to (3).  For fixed
`x1+x4`, the constraint is unchanged and AM--GM uniquely maximizes `x1x4`
at `x1=x4`.  For the other pair, set

```text
f(y)=log u(exp(y)).
```

The identity

```text
f''(y)=x[1-(1+x)exp(-x)]/(1-exp(-x))^2 >0,
x=exp(y),                                               (4)
```

shows that, at fixed `u(x2)u(x3)`, Jensen uniquely maximizes `x2x3` at
`x2=x3`.  Thus the saddle is necessarily

```text
c=(s,t,t,s).                                            (5)
```

The exact surface and critical equations are

```text
(exp(2s)-1)(exp(t)-1)^2=1,                              (6)
s/(1-exp(-2s))=t/(1-exp(-t)).                           (7)
```

Put `p=exp(t)-1`.  Then

```text
t=log(1+p),
s=(1/2)log(1+p^(-2)),                                   (8)
```

and `p` is the unique positive root of

```text
E(p):=p(1+p^2)log(1+p^(-2))-2(1+p)log(1+p)=0.           (9)
```

Uniqueness is not numerical.  In logarithmic coordinates `y_i=log x_i`,
the left side of (3) has logarithm

```text
G(y)=log(exp(exp(y1)+exp(y4))-1)
     +f(y2)+f(y3).                                      (10)
```

Its Hessian is positive definite: the last two summands use (4), while the
first is `f(log(exp(y1)+exp(y4)))`, whose Hessian is the positive sum of the
rank-one `f''` term and the transverse log-sum-exp Hessian.  Therefore
`{G<=0}` is strictly convex.  Its support point in direction `(1,1,1,1)` is
unique, and its Lagrange equations are exactly (7).  The product tends to zero
at both ends of the reduced surface, so that support point exists.  This also
proves that (9) has exactly one positive root.

## 2. Exact product radius and numerical isolation

The ordinary generating function `sum_r a_d(r)z^r` has radius

```text
R_Q4=s^2t^2
 =1/4 [log(1+p^(-2)) log(1+p)]^2
 =[(1+p)/(p(1+p^2))]^2 log(1+p)^4,                      (11)
```

where `p` is defined exactly by (9).  A high-precision isolation is

```text
0.3866667229 < p < 0.3866667230,
p = 0.3866667229877961774674905892220983809606625132006...

s = 1.0198605224431466419625042380265599014819888574395...
t = 0.3269028262212605653135267402312502491359268621201...

R_Q4 = 0.1111524174859334724422697977612252181685979362223... . (12)
```

This is an exact implicit analytic value; it is not presently reduced to an
algebraic closed form.

The equal positive pole is a useful hostile.  If `rho` is the positive root of
`rho^4-2rho-1=0`, then the equal pole is
`x_i=log(1+1/rho)=0.5403879692...`; its product is
`0.08527518825...`, strictly below (12), and its logarithmic gradients are
unequal.  It is a singular point, but not the diagonal saddle.

## 3. Strict complex minimality and the nondegenerate saddle

Suppose `|z_i|<=c_i` and `H(u(z_1),...,u(z_4))=0`.  Perron comparison gives

```text
1 <= rho(A diag(|u(z_i)|))
  <= rho(A diag(u(c_i)))=1.                             (13)
```

Strong irreducibility makes equality strict if any coordinate modulus is
smaller.  Hence every `|u(z_i)|=u(c_i)`.  Equality in
`|u(z)|<=u(|z|)<=u(c)` forces `|z|=c`, and equality between the first two
positive Taylor terms of `u` forces `z=c`.  Thus `c` is the only singularity
in its closed polydisc.

For an exact nondegeneracy certificate, define

```text
phi(x)=x/(1-exp(-x)),
chi(x)=x[1-(1+x)exp(-x)]/(1-exp(-x))^2 >0,
kappa=phi(t)=phi(2s)/2.                                 (14)
```

On the orthonormal tangent basis

```text
(1,0,0,-1)/sqrt(2),
(0,1,-1,0)/sqrt(2),
(1,-1,-1,1)/2,
```

the restriction of `Hess G` has eigenvalues

```text
kappa,
chi(t),
(chi(2s)/2+chi(t))/2,                                   (15)
```

all strictly positive.  The negative second variation of
`sum_i log x_i` along `G=0` is (15) divided by `kappa`, so the Gaussian saddle
is nondegenerate.  Numerically (15) is

```text
1.1723410467520150...,
0.1811991183308294...,
0.4980109050638674...,
```

and the Gaussian determinant is `0.06565790039939037...`.

The standard pole-residue calculation therefore gives, for every fixed
`d>=1`,

```text
a_d(r) ~ C_d r^(d-5/2) R_Q4^(-r),   C_d>0.              (16)
```

For `d=1`, an independent exact Stirling-profile expansion gives the
power-corrected ratio `9.0060197654...` at `r=40,41`, already close to
`R_Q4^(-1)=8.9966554270...`.

## 4. New exact arithmetic: the scale ratio is transcendental

Let

```text
q=s/t=(1+p)/(p(1+p^2))>1.                               (17)
```

Then `q` is irrational.  In fact, both `p` and `q` are transcendental.

First, `p<1`: otherwise
`2s=log(1+p^(-2))<=log(2)<=log(1+p)=t`, and strict increase of `phi`
would give `phi(2s)/2<phi(t)`, contradicting (7).  Formula (17) then gives
`q>1`.

To prove irrationality, suppose `q=m/n` in lowest positive terms.  Since
`q>1`, `m>n`.  Equations (6)--(8) give

```text
F(p):=m p^3+(m-n)p-n=0,                                 (18)
p^(2n)(1+p)^(2m)=(1+p^2)^n.                             (19)
```

If `p=a/b` is rational in lowest terms, (19) becomes

```text
a^(2n)(a+b)^(2m)=b^(2m)(a^2+b^2)^n.
```

Coprimality first forces `a=1`, then `b=1`; (19) would then force `n=2m`,
contrary to `m>n`.

If `p` is irrational, `F'(z)=3mz^2+(m-n)>0` on the real line.  Hence (18)
has only one real root.  A reducible rational cubic would have a rational real
root, which would have to be `p`; therefore `F` is irreducible.  In
`K=Q(p)`, exact evaluation of the cubic at `-1,+i,-i` gives

```text
N_K/Q(p)=n/m,
N_K/Q(1+p)=2,
N_K/Q(1+p^2)=2(n/m)^2.                                  (20)
```

Taking norms in (19) again yields `2^(2m)=2^n`, so `n=2m`, the same
contradiction.  Thus `q` is irrational.

If `p` were algebraic, then (17) would make `q` algebraic, while (6) says

```text
1+p^(-2)=(1+p)^(2q).                                    (21)
```

Gelfond--Schneider rules out algebraic irrational `q`, and rational `q` was
just excluded.  Therefore `p` is transcendental.

Finally, if `q` were algebraic, equation `(17)` would make `p` a root of

```text
q p^3+(q-1)p-1=0.
```

It would then be algebraic over the algebraic numbers and hence algebraic,
contradicting the preceding paragraph.  Thus `q=s/t` is transcendental as
well.

## 5. Why the transcendental-radius proof does not yet survive

THM-3213's D-finite singularity lemma says that P-recursiveness of the rational
sequence `(a_d(r))` would force its finite positive radius `R_Q4` to be
algebraic.  In the Perron-balanced case all coordinates equal one logarithm,
so Hermite--Lindemann immediately makes the corresponding power of that
logarithm transcendental.

Here `(11)` is a product of two different logarithmic scales.  The proved
transcendence of `p` and `q` does **not** imply that `R_Q4` is transcendental:
transcendental factors can have an algebraic product.  The tempting
conditional branch “`q` algebraic” is now known to be vacuous, because the
last paragraph of Section 4 proves that `q` is transcendental.

Schanuel's conjecture would nevertheless close the product-radius question.
If `R_Q4` were algebraic, put `L=st`, which is algebraic.  Equations
`s^2=Lq(p)`, `t^2=L/q(p)`, `exp(t)=1+p`, and
`exp(2s)=1+p^(-2)` make `Q(s,t,exp(s),exp(t))` algebraic over
`Qbar(p)`, so its transcendence degree is at most one.  But `s/t=q` is
irrational, hence `s,t` are Q-linearly independent, and Schanuel predicts
transcendence degree at least two.
Unconditionally, however, no theorem used in this workspace excludes that
escape.  Therefore:

```text
Q4 fixed-depth sequence:
  normalized non-C-finite  PROVED by THM-3213's prime shift;
  raw non-C-finite         PROVED by (16), since (r!)^4 makes its growth
                           superexponential;
  non-P-recursive          CONDITIONAL on R_Q4 transcendental
                           (in particular on Schanuel);
  unconditional status    OPEN.                         (23)
```

Fresh denominator primes cannot replace this step, since factorially
hypergeometric P-recursive sequences have the same phenomenon.
Multiplication and division by `(r!)^4` preserve P-recursiveness, so the raw
and factorial-normalized Q4 sequences share the same conditional/open
P-recursive status.

## 6. Strongest quotient-wide generalization

The **analytic saddle theorem** extends from Perron-balanced quotients to
every strongly connected tournament quotient `Q` on `q>=3` vertices.

Let

```text
psi(z)=log rho(A_Q diag(exp(z_i))),
G_Q(y)=psi(f(y_1),...,f(y_q)).                           (24)
```

Perron log-convexity gives `Hess psi>=0`.  If `v,w>0` are Perron vectors and
`pi_i=w_i v_i/(w^T v)`, then

```text
Hess G_Q
 =diag(f') Hess(psi) diag(f') + diag(pi_i f''(y_i))>0.  (25)
```

Camion's Hamilton cycle gives `prod_i u(x_i)<=1` throughout the positive
Neumann sublevel.  Hence

```text
prod_i x_i <= prod_i x_i/u(x_i).
```

If any coordinate tends to infinity, the right side tends to zero because
`x/u(x)->0`; if a coordinate tends to zero while the others stay bounded,
the product itself tends to zero.  A positive comparison point therefore
turns every maximizing sequence into a compact one.  Strict convexity in
`(25)` then gives one and only one maximizer; common scaling places it on the
boundary `rho=1`.  Denote this support point by `c_Q`.  It is characterized by

```text
rho(A_Q diag(u(c_i)))=1,
pi_i c_i u'(c_i)/u(c_i)=kappa  for every i.              (26)
```

The same comparison argument as (13) gives strict complex minimality, and
(25) gives a nondegenerate `(q-1)`-dimensional saddle.  Hence

```text
R_Q=prod_i c_i
   =max {prod_i x_i : rho(A_Q diag(u(x_i)))<=1},          (27)

F_(Q[T_r,...,T_r])(d)/(r!)^q
 ~C_(Q,d) r^(d-1-(q-1)/2) R_Q^(-r), C_(Q,d)>0.          (28)
```

This is an unconditional extension of the exact radius/asymptotic geometry to
**all strong fixed tournament quotients**.  It does not by itself extend the
unconditional non-P-recursiveness theorem: that requires a separate proof
that the variational value (27) is transcendental.  Perron-balanced quotients
remain the broadest combinatorial class certified by the present
transcendental-radius argument because (26) then collapses to one logarithmic
scale.  More generally, any individually
certified one-log algebraic ray saddle has the same transcendence proof.

For `q=1`, there is no Perron pole: `W=X`,

```text
F_(T_r)(d)=d!S(r,d),
sum_r F_(T_r)(d) z^r/r!=(exp(z)-1)^d,
```

so the raw sequence is C-finite and the factorial normalization is
P-recursive.  No strong tournament exists on `q=2`.  Reducible quotients split
by their dominant SCCs and are outside (24)--(28).

## 7. Exact and numerical controls

Run

```text
python3 04-computation/tournament_unbalanced_q4_saddle_thm3226.py
python3 -O 04-computation/tournament_unbalanced_q4_saddle_thm3226.py
```

and compare LF-normalized bytes with the declared output.  The companion
checks the determinant, numerator positivity, paired reduction, exact norm
identities, tangent Hessian signs, and recurrence-class hostiles symbolically
or with integers.  Its 100-digit saddle isolation and coefficient-ratio
comparison are explicitly **NUMERICAL CONTROL ONLY**; the uniqueness,
minimality, asymptotic, scale transcendence, Schanuel implication, and
quotient-wide variational theorem are proved analytically above.  Normal,
optimized, and stored transcripts agree byte-for-byte.

The recurrence taxonomy is deliberately split:

- normalized non-C-finiteness is inherited from THM-3213's prime shift, while
  raw non-C-finiteness follows from the superexponential form of `(16)`;
- non-P-recursiveness follows conditionally from transcendence of `R_Q4`;
- unconditional transcendence of `R_Q4`, and hence unconditional Q4
  non-P-recursiveness, remains **OPEN**.

**QED.**
