# Retaining the energy divisor doubles the global stability constant

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED (analytic proof).**
**Current update:** the [sharp global theorem](open_frontier_sep06_stability_complement.md)
now proves K_best=K3. The K0 inequality below remains valid; its claim that
the optimum is open is superseded by that mutually audited continuation.
This intermediate improvement retains, without editing,
the frozen [quantitative stability theorem](quantitative_stability_empty_core_morning_sep06.md).
No theorem ID or external-priority claim.

## 1. A stronger uniform bound

For any finite real list `r` with `p_1=p_2=1`, put

```text
G(s)=prod_i(1+r_i s),
E=e_2(r_i^2)=(1-p_4)/2>0,
D=-[s^4]G(s)^2=(5-8p_3+3p_4)/6,
J=D/E,       c_*=(13-8sqrt(2))/3.
```

Pad with zeros, and let `a>=b>=0` be the largest two positive entries,
using a zero when there is no second positive entry. As before,

```text
d^2=dist_l2(r,{permutations of (1/sqrt(2),1/sqrt(2),0,...)})^2
   =2-sqrt(2)(a+b).
```

Then

```text
J-c_* >= K_0 d^2,
K_0=(12sqrt(2)-16)/3 =0.3235209161... .                    (1)
```

Thus the proved lower bound for the best constant doubles. With the
inherited three-positive-atom upper obstruction,

```text
K_0 <= K_best <= K_3,
K_3=4sqrt(3)/[3(1+sqrt(2))(1+sqrt(3))]
   =0.3501345012... .                                    (2)
```

The distance power two remains optimal by the inherited actual
finite-degree extremizers. In particular the new estimate gives
`d<=sqrt((J-c_*)/K_0)`, dust square mass at most `(J-c_*)/K_0`, and
signed dust first-moment error at most `(J-c_*)/(sqrt(2)K_0)`.
Zeros, both signs, arbitrary list length and the one-atom energy boundary
are retained. Normalize an original cancellation list by its nonzero
first moment before applying these statements, as in the inherited note.

## 2. Inheritance and the one-sided compression

The recovery pass at HEAD `677bde8eb3ae` found the existing bounds,
the sharp finite-degree carriers, the exact sphere remainder and the
three-positive-atom obstruction, but no improved global constant in the
targeted current results/theorem search. The nearby primary-literature
and moment-image comparison from the inherited note is retained as
context; no external theorem is newly transported here.

The closest mechanism is the previous signed-tail envelope. The canonical
hostile is the one-atom boundary, where the unnormalized moment gap
vanishes while the quotient tends to one. The corrected near miss is the
local unequal-two-atom constant, which the three-atom example disproves
as a global bound. The least-used sidecar is the energy divisor itself:
the previous proof bounded it separately after compressing the tail.

The live board is: square coefficient; positive-root order; energy;
signed dust; integer root multiplicity; uniform-atom boundary. The new
map moves the target energy term into the fourth-moment coefficient
*before* using the tail envelope. It retains a one-sided inequality at
fixed top positive roots and loses tail multiplicities. Section 5 gives
an exact obstruction showing why that remaining loss matters.

## 3. Put the divisor into the objective first

Set

```text
B=sqrt(2)-1, m=B^2=3-2sqrt(2), h=1/sqrt(2),
A=2-sqrt(2), t=a+b, C=2B-mt.
```

Since `t<=sqrt(2)`, we have `C>=A>0`. Literal substitution of (1) and
the formulas for `D,E` shows that its target is equivalent to

```text
F_actual=m(1+t)-p_3+C p_4 >=0.                           (3)
```

Indeed `C=A+(3K_0/8)d^2`, and
`F_actual=(3E/4)(J-c_*-K_0 d^2)`.

For `f_C(u)=u-Cu^2`, the derivative is positive on `[0,b]`. A useful
exact certificate, using `0<=b<=h` and `a>=b`, is

```text
1-2Cb-m
  =2mb(a-b)+4m(h-b)(1+h-b) >=0.                         (4)
```

Thus `f_C'(u)>=m>0` for `0<=u<=b`. Every positive tail entry contributes
`r_i^2 f_C(r_i)<=r_i^2 f_C(b)` to `p_3-Cp_4`; a nonpositive entry
contributes at most zero, hence also at most `r_i^2 f_C(b)`.
Consequently

```text
F_actual >= F(a,b),
F(a,b)=m(1+a+b)-a^3-(1-a^2)b
       +[2B-m(a+b)] [a^4+(1-a^2)b^2].                   (5)
```

It remains to prove a two-variable inequality on
`a>=b>=0`, `a^2+b^2<=1`. No sign replacement or tail deletion has been
made to the actual list.

## 4. An exact monotonicity reduces the envelope to two boundaries

Differentiate the polynomial in (5) at fixed `a`:

```text
partial_b F=(1-a^2) T(a,b),
T=m(1+a^2-2ab-3b^2)+4Bb-1.
```

On the declared domain,

```text
T(a,b)<=U(b)=m(2-6b^2)+4Bb-1,
U-T=m[(1-a^2-b^2)+2b(a-b)]>=0,
U(b)=-6m(b-h)(b-(4+sqrt(2))/6)<=0.                      (6)
```

The final sign follows from `b<=h<(4+sqrt(2))/6`. Therefore `F` is
nonincreasing in `b` at fixed `a`. It suffices to test
`b=min(a,sqrt(1-a^2))`.

If `a<=h`, this is the diagonal. Exact factorization gives

```text
F(a,a)=2m(1-a)(h-a)^2>=0.                               (7)
```

If `a>=h`, set `s=sqrt(1-a^2)`. At `a=1`, `F(1,0)=0` directly.
Otherwise `s>0`, `t=a+s in (1,sqrt(2)]`, and the two-atom energy is
`E_2=a^2s^2`. The inherited sphere identity, or direct elimination of
`as=(t^2-1)/2`, gives

```text
g_2=B-a^3-s^3+A(a^4+s^4)
   =(t-1)^2(sqrt(2)-t)(At+1)/2.
```

When `t<sqrt(2)`, division by the positive factors yields

```text
g_2/[E_2(2-sqrt(2)t)]
   =sqrt(2)(At+1)/(t+1)^2
  >=sqrt(2)/(1+sqrt(2))^2
   =sqrt(2)m=3K_0/4.                                   (8)
```

This proves `F(a,s)>=0` through the exact identity in (3). At
`t=sqrt(2)` both sides of the undivided inequality vanish, so it holds
there as well. Equations (5)--(8) prove (1) for every eligible finite
list. This boundary reduction concerns the envelope, not the still-open
global extremizer classification of the actual root lists.

## 5. The new constant is sharp for this compression, not yet for the roots

The tail envelope becomes an equality at the diagonal `a=b=x` only if
all nonzero tail entries also equal `x` and all entries are nonnegative.
This follows from the strict monotonicity in (4) and the strict loss
at every negative nonzero entry. Thus `N x^2=1` for an integer `N>=2`
is a necessary equality condition for the square-mass relaxation.
It is not sufficient for a finite list with the additional normalization
`p_1=1`: the macroscopic configurations below are realized as dust limits,
with the missing signed first moment kept separately.

The envelope itself does not retain that integer. Its diagonal formally
allows `p_3=x`, `p_4=x^2`, `t=2x` for arbitrary `0<x<h`. For that relaxed
data the quotient satisfies

```text
J=(5-3x)/(3(1+x)),
(J-c_*)/(2-2sqrt(2)x)
  =4/[3sqrt(2)(1+x)(1+h)] -> K_0 as x->h from below.     (9)
```

Therefore no stronger constant can follow from this same two-coordinate,
square-mass envelope alone. The obstructing interval
`1/sqrt(3)<x<1/sqrt(2)` would require a number of equal roots strictly
between two and three. This is an exact failure of the relaxation,
not an actual hostile to a larger constant for finite real-rooted rows.
The missing coordinate is the integer packing of the remaining square
mass into roots no larger than `b`.

For comparison, `N>=3` genuine equal macroscopic atoms of value
`1/sqrt(N)`, with vanishing-square-mass negative dust correcting `p_1`,
give exactly the limiting obstruction

```text
K_N=4sqrt(N)/[3(1+sqrt(2))(1+sqrt(N))].                   (10)
```

It is strictly increasing in integer `N>=3`; hence the smallest such
obstruction is `N=3`, the upper bound in (2). The finite realizing family
has `N` raw roots equal to one and `L` roots equal to `-q`, where

```text
q=(N-1)/[L+sqrt(L(L+N-1)/N)].
```

It solves `N(N-1)-2NLq+L(L-1)q^2=0`, hence has exact `e_2=0`.
Normalize by its positive sum, which tends to `sqrt(N)` as `L->infinity`.
The case `N=2` instead approaches the vanishing-distance equality shape
and has the different dust-direction ratio from the inherited theorem;
one must not substitute `N=2` into (10) as an actual finite-row sharpness
claim.

Cheap hostile probes of unequal positive three-atom shapes and one
distinguished positive or negative atom did not improve the known
three-atom upper obstruction. They are bounded controls only. The next
precise obstacle is to retain integer tail packing while minimizing at
fixed top-positive coordinates; the present proof neither assumes that
the third root vanishes nor promotes a local Hessian constant.

## 6. Reproduction and audit manifest

[Source](../../04-computation/open_frontier_sep06_stability.py) and
[output](open_frontier_sep06_stability.out).

```bash
python3 -B 04-computation/open_frontier_sep06_stability.py
python3 -B -O 04-computation/open_frontier_sep06_stability.py
```

All **341 explicit gates** pass. Formal bivariate arithmetic in
`Q(sqrt(2))` checks the complete target, secant, derivative, domain
remainder and diagonal identities. The declared signed finite universe
is the inherited 49 rational prefixes, with all six exceptions printed
and all 43 eligible rows tested against the stronger bound through
literal ordinary-square coefficients and full root energy. Zero entries
and negative normalization sums are retained.

The source also checks 21 actual equal-macro-plus-dust families
(`N=2,...,8`, `L=2,8,100`), six signed three-macro moment controls, three
exact one-atom-boundary rows, and three relaxed diagonal controls. The
six macro controls test the square-normalized moment inequality (3);
they are not asserted to be literal square-coefficient controls without
the first-moment normalization. Compressed actual-family coefficients
are compared against individual-factor multiplication in the small
rows. Exact signs of the needed biquadratic expressions avoid floating
root or sign tests. No unbounded numerical search or finite-height
extrapolation is used in the proof.

The root agent independently read and audited Sections 1--5, including
the energy factor in (3), signed secant, both derivative factorizations,
the two boundary cases, sharpness for the relaxation, and the integer
`K_N` formula: **PASS**. Its scope clarification distinguishing a
necessary relaxation equality condition from finite `p_1=1` feasibility
is included above. Normal and optimized outputs agree byte for byte.

```text
source SHA256 c431c98cba8aef4aa1462b69bc5c5eda38712358ba3bfbfef376be2dcb881a63
output SHA256 d72fdcf02afd9cdd53a55972ab775a7f2dc6e751f62afcfc75fc810a52bd63f4
trace  SHA256 aa4036db4aaec90e3c438eb7e13fa1cf4367b67636fa30a073c0a3ba1cdcd391
```
