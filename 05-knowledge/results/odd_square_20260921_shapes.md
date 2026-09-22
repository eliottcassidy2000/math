# Odd-square fibres and the measure on right-triangle shapes

**Status: PROVED elementary counting and change-of-variable theorems;
FINITE-EXACT controls; no priority claim.** Date: 2026-09-21. All limiting
probabilities below specify a finite sampling space first. None is a
statement about a Collatz trajectory or a single Berggren/Gaussian orbit.

## 1. Inheritance and the live board

The shape identities are inherited from
[arithmetic braids, geometry §§4 and 6](arithmetic_braids_20260917_geometry.md)
and [the semicircle continuation](collatz_mod6_20260917_pythagorean_semicircle.md).
The complete odd-square chart and its exact fibre size are already
[THM-3756, odd-square ordinal chart and affine descent](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md).
The closest counting mechanism is the Möbius/visible-point argument in
[the signed density note](glued_xor_20260921_density.md). Its canonical
hostile is that a gcd-preserving finite orbit has primitive frequency zero
or one even when the surrounding lattice has density `6/pi^2`.

The corrected geometric near miss is the literal `3-4-5` altitude: it is
`12/5`, or `12/25` after hypotenuse normalization, not `sqrt(2)`. The
least-used sidecar in this continuation is the **height used to sample**:
fixing or bounding the outer odd square and bounding the hypotenuse give
different shape measures on exactly the same triples.

| Lane | Object and question | Preserved coordinate / hostile probe |
|---|---|---|
| Anchor | `(e/d,l,theta)` on the unit-hypotenuse semicircle | One shape, three dependent coordinates; test `3-4-5` |
| Niche | A fixed odd-square fibre | Retain inner root `v`; compare primes and highly composite `u` |
| Wildcard | Change hypotenuse height to outer-root height | Same support, different limiting measure |
| Bridge | Coprime lattice points with parity | Separate unconditional and conditional densities at prime 2 |
| Boundary | Integer square-plus-one family | A thin subfamily can approach degeneration despite a uniform bulk law |

The new work is the explicit measure transport, uniform fibre discrepancy,
and the contrast of these two arithmetic heights. The chart itself is not
being rediscovered or renamed.

## 2. One shape interval and the epsilon constraints

Normalize the hypotenuse to one. Let `a<b` be the shorter and longer legs,
and drop the altitude from the right angle. Its foot splits the hypotenuse
into `e=a^2` and `d=b^2`. Write `l=ab` and let `theta` be the angle between
the altitude and the shorter leg. Similar right triangles give

```text
e+d=1,              l^2=ed,
r=e/d=tan(theta)^2, 0<theta<pi/4,
e=r/(1+r),         d=1/(1+r),
l=sqrt(r)/(1+r)=sin(2theta)/2.                         (OS1)
```

Thus `0<r<1` and `0<l<1/2`. These coordinates are strictly increasing
functions of each other, not three independently adjustable variables.
The isosceles endpoint is allowed for real triangles, but is excluded for
nondegenerate integer Pythagorean triangles. The third vertex `(e,l)` lies
on `(e-1/2)^2+l^2=1/4`: hypotenuse one is the **diameter**, with radius one
half. The speed of this vertex as a function of `theta` is one, so uniform
`theta` is also uniform arc length on the folded quarter circle.

For `3-4-5`, `(r,l)=(9/16,12/25)` and `theta=arctan(3/4)`.
At every interior point `l<theta`, since `sin(2theta)<2theta`.
Consequently the simultaneous strict bounds

```text
epsilon<r<1, epsilon<l<1/2, epsilon<theta<pi/4
```

are empty if `epsilon>=1/2`. For `0<epsilon<1/2`, they are equivalent to

```text
theta > Theta(epsilon),
Theta(epsilon)=max(arctan(sqrt(epsilon)),
                   arcsin(2epsilon)/2).                (OS2)
```

The angle cutoff is redundant once the altitude cutoff is imposed. The
binding coordinate switches at the unique root `epsilon_*` in `(0,1/2)`
of `epsilon(1+epsilon)^2=1`, approximately `0.4655712319`: below it the
ratio cutoff binds, above it the altitude cutoff binds. This follows by
evaluating `l` at `r=epsilon`; the cubic's derivative is positive on this
interval. A common epsilon here uses radians for the angle.

## 3. Why the constants are 6, 4, 2, and 8 over pi squared

All densities in the next table refer to large two-dimensional regions
whose dilates have area of order `R^2` and boundary length of order `R`.
The sectors and triangles needed below satisfy this condition.

| Predicate on integer pairs | Density among all pairs | Conditional gcd-one density |
|---|---:|---:|
| Coprime, no parity restriction | `6/pi^2` | — |
| Coprime and opposite parity | `4/pi^2` | `8/pi^2` among opposite-parity pairs |
| Coprime and both odd | `2/pi^2` | `8/pi^2` among odd-odd pairs |

In particular, two thirds of coprime pairs have opposite parity; one third
are both odd. Treating the parity of a coprime pair as two independent fair
bits would give the wrong constants.

**Proof with an error term.** A fixed parity class has one quarter of the
area, up to `O(R+1)`, by tiling with two-by-two squares and charging only
boundary squares. Opposite parity consists of two such classes. If the
parity is opposite, every common divisor is odd, so Möbius inversion gives

```text
sum_(d odd, d<=O(R)) mu(d)
  [area/(2d^2)+O(R/d+1)]
 = (4/pi^2) area + O(R log(2R)),                        (OS3)
```

because

```text
sum_(d odd) mu(d)/d^2
 = product_(p odd)(1-p^(-2)) = 8/pi^2.                 (OS4)
```

The omitted tail of the main term is `O(R)`. The same calculation with
`area/(4d^2)` proves the odd-odd row. Without a parity restriction, sum
over all `d` to obtain `6/pi^2`. These arguments also hold for a fixed
finite union of sectors, including choices of open or closed radial edges.
Boundary changes contribute to the stated error.

The equality with the odd-squarefree density `8/pi^2` has a precise local
meaning: at each odd prime the forbidden event has relative size `1/p^2`.
It does not identify individual squarefree integers with primitive triples,
and it does not manufacture an orbit sampling law. For primary-source
context, Carl Pomerance's author-hosted
[Euler phi lecture, slides 8–9](https://math.dartmouth.edu/~carlp/PDF/phitalk2.pdf)
records the ordinary visible-point constant `6/pi^2`. The parity refinements
and all error estimates used here are proved above, rather than imported
from a stronger unexamined theorem.

## 4. Hypotenuse height makes the smaller angle uniform

Let `P_X` be the finite set of positive primitive Pythagorean triples with
hypotenuse `c<=X`, counting each unordered pair of legs once. Equivalently,
retain the labels `(A,B,c)=(odd leg,even leg,hypotenuse)` and use the unique
Euclid pair

```text
(A,B,c)=(m^2-n^2,2mn,m^2+n^2),
m>n>0, gcd(m,n)=1, m-n odd.                            (OS5)
```

For fixed `0<=alpha_0<alpha_1<=pi/4`, the points with
`alpha=arctan(n/m)` in that interval lie in a circular sector of area
`X(alpha_1-alpha_0)/2`. Applying (OS3) gives

```text
# {P_X : alpha_0<alpha<alpha_1}
 = 2(alpha_1-alpha_0)X/pi^2 + O(sqrt(X) log(2X)),
# P_X = X/(2pi) + O(sqrt(X) log(2X)).                  (OS6)
```

The constants can be taken uniformly over such sectors because their
boundary lengths are uniformly bounded multiples of `sqrt(X)`.
Euclid's leg angle is `2alpha`, while the smaller angle is

```text
theta=min(2alpha,pi/2-2alpha).                         (OS7)
```

The inverse image of a `theta` interval of length `h` consists of two
`alpha` intervals, each of length `h/2`. Thus, for every fixed
`0<=t<=pi/4`,

```text
lim_(X->infinity) # {P_X : theta<=t}/#P_X = 4t/pi.      (OS8)
```

This is a bulk height theorem, not equidistribution of any one tree ray.
The limit probability laws of the other coordinates follow by monotone
change of variable:

| Coordinate | Limit cumulative distribution | Limit density |
|---|---|---|
| `r=e/d`, `0<r<1` | `4 arctan(sqrt(r))/pi` | `2/[pi sqrt(r)(1+r)]` |
| `l`, `0<l<1/2` | `2 arcsin(2l)/pi` | `4/[pi sqrt(1-4l^2)]` |
| `e`, `0<e<1/2` | `4 arcsin(sqrt(e))/pi` | `2/[pi sqrt(e(1-e))]` |

Each density integrates to one. Their endpoint singularities are
integrable; they are not atoms. In particular, (OS2) has limiting
probability `1-4Theta(epsilon)/pi`, for `0<epsilon<1/2`.

For comparison, choosing a real coordinate uniformly is a different
operation:

| Uniformly sampled coordinate | Resulting density in `theta` |
|---|---|
| `theta in (0,pi/4)` | `4/pi` |
| `r in (0,1)` | `2 tan(theta) sec(theta)^2` |
| `l in (0,1/2)` | `2 cos(2theta)` |

There is no measure-free meaning of a uniformly chosen shape.

## 5. Fixed odd-square fibres have a different, uniform slope law

Use the inherited odd-square roots

```text
u=sqrt(c+B)=m+n, v=sqrt(c-B)=m-n,
u>v>0 odd, gcd(u,v)=1,
(A,B,c)=(uv,(u^2-v^2)/2,(u^2+v^2)/2).                 (OS9)
```

For every fixed odd integer `u>=3`, the number of allowed `v` is
`N_u=phi(u)/2`. Pairing reduced residues `v` and `u-v` proves the count:
exactly one member of each pair is odd. This is THM-3756's fibre theorem.

The following quantitative refinement controls *every* large fibre,
including highly composite `u`. Let

```text
F_u(t)=#{v odd: 0<v<u, gcd(u,v)=1, v/u<=t}/N_u,
0<=t<=1.
```

Then

```text
sup_t |F_u(t)-t| <= 2^omega(u)/phi(u)
                  <= sqrt(15)/(2sqrt(u)).             (OS10)
```

**Proof.** If `O(y)` counts positive odd integers at most the real number
`y`, then `|O(y)-y/2|<=1/2`. For `0<=t<=1`, Möbius inversion gives

```text
N_u F_u(t) = sum_(d|u) mu(d) O(tu/d)
           = t phi(u)/2 + E,
|E| <= 2^(omega(u)-1).                                (OS11)
```

At `t=1`, including `v=u` in the inversion causes no problem: its total
weight is zero because `u>1`. Divide by `N_u`. To obtain the last bound,
factor `phi(u)/2^omega(u)`. For a prime power `p^a`, its factor divided by
`sqrt(p^a)` is `p^(a/2-1)(p-1)/2`. It is at least one for `p>=7`, at least
`1/sqrt(3)` for `p=3`, and at least `2/sqrt(5)` for `p=5`. Multiplication
gives `phi(u)/2^omega(u)>=2sqrt(u)/sqrt(15)`. Equality in this latter
inequality occurs at `u=15`. QED.

Thus `t=v/u` tends to the uniform distribution along **all** odd
`u->infinity`, with no primality or squarefreeness assumption on `u`.
Writing `beta=arctan(t)`, the smaller angle is again
`theta=min(2beta,pi/2-2beta)`. Consequently its limit CDF in a fixed
outer-root fibre is

```text
G(theta)=tan(theta/2)+1-tan(pi/4-theta/2),
G'(theta)=[sec(theta/2)^2+sec(pi/4-theta/2)^2]/2.        (OS12)
```

For the finite fibre, the angle CDF differs from `G` by at most
`2^(omega(u)+1)/phi(u)`: apply (OS10) at the two preimage endpoints,
using the left limit for the second endpoint. The same bound holds for
left limits because the comparison CDF is continuous.

This law is not uniform in angle: its density starts at `3/2` and ends
at `sec(pi/8)^2`, whereas the hypotenuse law has constant density `4/pi`.
For example, at the `3-4-5` angle `theta=arctan(3/4)`, (OS12) is `5/6`,
while (OS8) is `4 arctan(3/4)/pi`.

## 6. Aggregate outer-root height and the missing radial weight

Let `Q_U` contain each primitive triple whose `u` in (OS9) satisfies
`u<=U`. Odd-odd primitive pairs occupy the triangle `0<v<u<=U` with
density `2/pi^2`. Hence

```text
# Q_U = (1/2) sum_(3<=u<=U, u odd) phi(u)
      = U^2/pi^2 + O(U log(2U)).                       (OS13)
```

Cutting this triangle by `v<=tu` cuts its area by the proportion `t`.
Thus `v/u` is asymptotically uniform under this aggregate height also,
and its shape CDF is (OS12). This conclusion can equally be deduced from
(OS10), after discarding finitely many small fibres.

The different hypotenuse law is explained by an exact area calculation,
not by a change in which triangles exist. In coordinates `(u,t=v/u)`,
area is `u du dt`. The restriction `u<=U` gives a constant `t` weight
`U^2/2`. But `c<=X` gives `u<=sqrt(2X/(1+t^2))`, so its weight is
`X/(1+t^2)`. After normalization, the latter slope density is

```text
4/[pi(1+t^2)], 0<t<1.                                (OS14)
```

Changing to `beta=arctan(t)` gives uniform `beta`, and folding gives
(OS8). Formula (OS14) retains precisely the radial information lost by
calling all heights the same sampling scheme.

The integer square-plus-one family is a separate thin-set test: the
companion [geometry note](odd_square_20260921_triangles.md) identifies its
two odd-root boundary families. Its count is of order `sqrt(X)` under
`c<=X`, hence its fraction of all triples tends to zero by (OS6). Uniform
bulk angle distribution does not imply a thin family's angle distribution.

## 7. Exact controls and the connection contract

The companion
[standard-library Python experiment](../../04-computation/experiments/odd_square_20260921_shapes.py)
uses explicit `require` calls that remain active under `python -O`.
Its [JSON certificate](../../04-computation/experiments/odd_square_20260921_shapes.json)
declares each finite universe. Ordinary and optimized Python produce
byte-identical JSON, with **480,912 active checks** and all **159,139**
primitive triples through hypotenuse `10^6`. Independent odd-root enumeration
agrees through `10^4`, and direct integer-leg enumeration agrees through
`300`. Exact integer/Fraction checks cover the
shape identities, the two parametrizations, fibre sizes and discrepancy,
epsilon cutoff predicates, and parity-conditioned Möbius counts. Decimal
asymptotic comparisons are separately labelled diagnostics; no finite
approximation is substituted for a proof of a limit.

Reproduce from the repository root:

```text
python 04-computation/experiments/odd_square_20260921_shapes.py
python -O 04-computation/experiments/odd_square_20260921_shapes.py --output .scratch/odd_square_shapes_optimized.json
```

| Source -> target | Map | Preserves | Forgets / required sidecar |
|---|---|---|---|
| Primitive odd/even Euclid pair -> triangle | (OS5) | Integrality, primitivity, hypotenuse | Pair is recoverable with parity-labelled legs |
| Triangle -> `(r,l,theta)` | Normalize and sort legs | Unordered similarity class | Absolute scale; parity label of shorter leg |
| Odd-root fibre -> unit slope | `v -> v/u` | Order, exact shape with `u` | Denominator and coprimality need `u` |
| Outer-root height -> hypotenuse height | `c=u^2(1+t^2)/2` | Same arithmetic support | Radial cutoff changes the probability measure |
| Coprimality -> Euler product | Möbius inversion | Aggregate densities with parity | Individual prime pattern, chronology, orbit access |

**Boundary.** Primitive-triple shapes are countable and dense, and have
Lebesgue measure zero as a subset of the continuum. Their height-weighted
empirical measures nevertheless converge to the continuous laws proved
above. None of these statements says a single arithmetic orbit is typical,
terminates, or samples squarefree integers independently.

**Independent audit:** the root agent read and accepted (OS1)–(OS14),
including the parity factors, uniform all-odd-fibre discrepancy constant,
epsilon threshold, and both height laws. The proof is not inferred from
the finite certificate.
