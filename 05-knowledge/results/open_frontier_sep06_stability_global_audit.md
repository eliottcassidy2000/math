# Third independent audit of the sharp global signed-root stability constant

**Status: INDEPENDENT ANALYTIC AND SOURCE AUDIT PASS.** This referee accepts
the combination of the
[complementary proof](open_frontier_sep06_stability_complement.md) and
[regional packing proof](open_frontier_sep06_stability_packing.md), including
global coverage, strictness for every finite eligible list, and actual
sharpness. No mathematical repair was required.

The exact best constant is

    K3 = 4 sqrt(3) / [3(1+sqrt(2))(1+sqrt(3))].

For finite real lists with `p1=p2=1`, `E=(1-p4)/2>0`, and the declared
top-two-positive-root distance d, the theorem states

    J-c_* > K3 d^2.

The strict inequality and sharp infimum are compatible. The infimum is
approached by actual finite lists; it is not attained by one.

## 1. The two regions have identical mathematical types

Write `u=sqrt(2)`, `v=sqrt(3)`, `r=1/v`, `h=1/u`,
`A=2-u`, `B=u-1`, `alpha=3K3/8`, and `gamma=u alpha`.
Both proofs use the same normalized quantities and the same coefficient

    C=C0-gamma(a+b)=A+alpha d^2,
    C0=A+2alpha=(1+u+v-uv)/2.

Independent substitution of the moment formula gives

    B-alpha d^2-p3+C p4
      = (3E/4)(J-c_*-K3 d^2).                       (1)

Thus neither region proves a different surrogate target before the union
of domains is taken. Their only surrogates are explicitly one-sided
objective bounds, and neither claims to preserve the first moment.

The normalization implies `a^2+b^2<=1` and `a+b<=sqrt(2)`, including lists
with fewer positive entries when padded with zero. Equality in the latter
inequality would force exactly two roots `1/sqrt(2)` and every other root
zero, contradicting p1=1. Hence d^2>0 on every eligible finite list.
The two regions `b<=r` and `b>=r` overlap on their common boundary and
exhaust the domain without a missing sequence, sign, or multiplicity case.

## 2. Complementary-region coverage and equality

On `b<=r`, the bounds `C>=A>0`, `C<=C0-2gamma b` and monotonicity of
`b(C0-2gamma b)` give `1-2Cb>=2-sqrt(3)>0`. Thus `x-Cx^2` is increasing
on `[0,b]`. For each positive tail root, its contribution to `p3-Cp4`
is bounded above by its square mass times `b-Cb^2`; every negative root
has a strictly smaller contribution than its absolute-value replacement.
This gives the stated signed secant envelope, with no finite-tail
restriction beyond the original domain.

The envelope derivative is exactly `(1-a^2)T`. The upper bound U satisfies

    U-T=gamma[(1-a^2-b^2)+2b(a-b)]>=0.

On the full interval `[0,r]`, U is increasing and U(r)<0. The latter sign
is equivalent to `sqrt(6)>1+sqrt(2)`, whose squared comparison reduces to
9>8. Hence the envelope strictly decreases with b when a<1.

For each fixed a, the largest permitted b is
`min(a,r,sqrt(1-a^2))`. The three ranges in the complementary proof are
therefore exhaustive: `a<=r`, `r<=a<=sqrt(2/3)`, and
`sqrt(2/3)<=a<=1`. Their factorizations have the stated signs. In the
fixed-r boundary, the cubic decreases to a positive right endpoint.
In the two-root boundary, the exact moment identity gives a quotient
strictly above 1/2, while K3<1/2. The positive denominator factors used in
that division hold throughout the open boundary; its a=1 endpoint is
treated separately.

The only zeros of this envelope in the complementary domain are `(r,r)`
and `(1,0)`. The second forces E=0 and is excluded. At `(r,r)`, the first
two positive roots already sum to more than one. The actual p1=1 list
must therefore contain a negative root, making the signed comparison
strict. Thus actual complementary lists have a strictly positive gap in (1).

## 3. Regional packing, domain geometry, and strictness

For `b>=r`, the entire tail square mass is `c^2=1-a^2-b^2<=1/3`, so
`0<=c<=r<=b`. Also `C<=C3=(3-sqrt(3))/2`. The function

    H_C(x)=x^(3/2)-C x^2

has second derivative at least `7sqrt(3)/4-3>0` on the relevant positive
interval. Its convexity and H_C(0)=0 imply
`sum H_C(x_i)<=H_C(sum x_i)` by the elementary chord inequality.
Replacing a negative root by its positive absolute value strictly raises
its cubic contribution. Since `a+b>=2/sqrt(3)>1`, every actual list in
this region contains a negative tail. Consequently the comparison with
the three-nonnegative-root surrogate is always strict.

At fixed c, the feasible interval of t=a+b is exactly

    r+sqrt(2/3-c^2) <= t <= sqrt(2(1-c^2)).

This follows by varying b from r to a at fixed `a^2+b^2=1-c^2`.
The second derivative of the surrogate objective was independently
rederived from its two moment expressions. Its c^2 coefficient is
negative. The resulting one-variable upper bound Q has Q'<0 on the
enclosing rational interval, and Q(2/sqrt(3))<0. The displayed rational
square comparisons prove these signs on the entire interval, not only
at sampled t values. Thus the objective is strictly concave in t.

The two boundary reductions are valid with their exact constraints. At
`b=r`, monotonicity of `x-Cx^2` lowers the objective to the already positive
fixed-r envelope. At `a=b=x`, the negative diagonal envelope is compensated
by the exact term `c^2(x-c)[1-C(x+c)]`. Using

    x-c=3(x-r)(x+r)/(x+c),
    h-x=c^2/[2(h+x)]

gives the stated positive-brace factorization. Its denominators are
strictly positive throughout the closed boundary. The inequalities
`x+c<=2r`, `1-C(x+c)>=2-sqrt(3)>1/4`, and `gamma<1/4` bound the braces
strictly above 1/2. Endpoints are retained through the undivided formula.

Strict concavity and the boundary analysis leave only the two relaxed
zero shapes: three equal positive roots, or two equal positive roots with
zero tail. Neither is an actual p1=1 list. More directly, the signed-tail
comparison was already strict, so every actual regional list satisfies (1)
strictly. Combining this with Section 2 proves the global inequality.

## 4. Actual finite sharpness, rather than a relaxed extremizer

For n>=3, the raw positive roots are `1+1/n,1+1/n,1-2/n`, with sum 3 and
square sum `Q_n=3+6/n^2`. Put `L=n^2` and append L negative roots -q_n,
where the positive expression in the packing theorem is the rationalized
small root of

    L(L-1)q_n^2-6Lq_n+9-Q_n=0.

It follows exactly that `S_n=3-Lq_n>0` and
`S_n^2=Q_n+Lq_n^2`. Dividing every root by S_n therefore gives both p1=1
and p2=1. There are several nonzero entries, so E>0.

The estimate `0<q_n<2/L` gives

    S_n^2<3+10/n^2<=3(1+1/n)^2.

Hence its two largest positive roots lie strictly above 1/sqrt(3) for
every n>=3. The actual sharpness sequence is inside the regional domain,
not merely approaching its boundary from the wrong side.

As n grows, S_n tends to sqrt(3), the three positive roots tend to
1/sqrt(3), and the dust square, cube and fourth-power contributions vanish.
The negative first moment remains the one needed for p1=1. Consequently

    E -> 1/3,
    d^2 -> 2-2sqrt(2/3)>0,
    J -> 3-4/sqrt(3),
    (J-c_*)/d^2 -> K3.

Both limiting denominator factors are positive, so there is no hidden
zero-over-zero sharpness inference. No larger constant can satisfy the
global inequality on all finite eligible lists. These facts establish the
claimed best constant with no finite equality case.

## 5. Source review and reproduction

The exact field arithmetic and formal-polynomial checks in both sources
were read. The two implementations use different quadratic-field towers
for their signs; each reduces an opposite-sign comparison to a squared
norm while retaining the original sign. This is valid over the stated real
embedding. The formal derivative and boundary identities retain all domain
relations, and the finite-row checks do not supply an extrapolation premise.

The following four runs were independently replayed. Each normal and
optimized output is byte-identical to its frozen output:

```
python3 -B 04-computation/open_frontier_sep06_stability_complement.py
python3 -B -O 04-computation/open_frontier_sep06_stability_complement.py
python3 -B 04-computation/open_frontier_sep06_stability_packing.py
python3 -B -O 04-computation/open_frontier_sep06_stability_packing.py
```

Complement: 90 exact gates; source SHA256
`bc1c612f840d5e9f1e8135fa17086f3244442d13962e8e37b633080c77c2d5a1`;
output SHA256
`c43d2e8da2947154415299d155766777a3384754221d4526dff92b00ed3dda9d`.

Packing: 72 exact gates; source SHA256
`a31fc847573b9fbc2799863057870c0be90c7645bf2c98b76da0c516c202c09b`;
output SHA256
`98bd68e9ffb7023bc838764768d6f56a3b840b21651381540f33ca9148ab654c`.

This audit concerns the signed-root optimization problem. It does not
claim an actual binomial-response identity, a trinomial separation theorem,
or an LRC consequence beyond already justified maps.
