# An explicit dimension-independent signed-duplication stability bound

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
New continuation of the independently audited
[uniform stability theorem](signed_uniform_stability_empty_core_next_sep06.md).
The constant below is explicit and nonsharp; the best global constant is
**OPEN**. No theorem ID or external-priority claim.

## 1. Statement and normalization

Let `r=(r_1,...,r_n)` be a finite real list with

```text
p_1=sum r_i=1,       p_2=sum r_i^2=1,
G(s)=prod_i(1+r_i s),
E=e_2(r_1^2,...,r_n^2)=(1-p_4)/2>0,
D=-[s^4]G(s)^2=(5-8p_3+3p_4)/6,
J=D/E,              c_*=(13-8sqrt(2))/3.
```

Zeros and both root signs are allowed. Pad with zeros and let `a>=b>=0`
be the two largest positive entries, replacing missing positive entries
by zero. The square distance to permutations of the two-atom vector is

```text
d(r)^2 := inf_pi ||r-pi(1/sqrt(2),1/sqrt(2),0,...)||_2^2
        = 2-sqrt(2)(a+b).                                  (1)
```

Then, independently of the list length,

```text
J-c_* >= K d(r)^2,
K=(6sqrt(2)-8)/3 = 4(3-2sqrt(2))/(3sqrt(2)) >0.             (2)
```

Thus `d(r)<=sqrt((J-c_*)/K)` is an explicit form of the inherited
qualitative stability statement. The normalized signed dust first moment
is still retained exactly: after choosing `a,b`, the remaining entries
have sum `1-a-b`. Nothing here asserts a one-sign or equal-multiplicity
dust pattern.

More explicitly, write `Delta=J-c_*` and let `eta^2` be the remaining
square mass after choosing `a,b`. Cauchy--Schwarz gives
`eta^2<=d^2-d^4/4<=d^2`, and the distance identity gives

```text
eta^2 <= Delta/K,
sum dust-(1-sqrt(2)) = d^2/sqrt(2) <= Delta/(sqrt(2)K).   (2a)
```

Thus both dust square mass and its signed first-moment error have an
explicit bound linear in the quotient deficit.

For an unnormalized list with `e_2=0` and positive energy, its nonzero
sum satisfies `p_1^2=p_2`; divide by that sum before interpreting the
positive atoms and the distance. A real polynomial scalar cancels from
`J`; monomials must first be removed with the coefficient index shifted.
These are the same real-gauge conventions as the inherited notes.

## 2. Recovery, map, and hostile controls

The closest proved mechanism is the sharp sphere inequality and its exact
two-support remainder in the uniform stability note. The closest hostile
boundary is the one-atom zero-energy configuration, at which the quotient
tends to `1` even though the unnormalized moment gap tends to zero.
The least-used sidecar is the top two *positive* coordinates: sorting by
absolute value would discard the sign information needed by the target
distance. The live board is: actual square coefficient; energy divisor;
largest positive roots; signed tails; unequal two-atom limits; three-atom
global obstruction.

A bounded repository search recovered the qualitative result, the sharp
finite-degree minima, and the earlier quartic quantitative remainder, but
no checked global estimate of the present form. Nearby primary literature
includes Acevedo--Blekherman--Debus--Riener,
[*The Wonderful Geometry of the Vandermonde map*](https://link.springer.com/article/10.1007/s10208-025-09718-6),
which studies finite- and growing-dimension images of the nonnegative
simplex under power-sum maps. The possible map here is `x_i=r_i^2`,
retaining square mass and unsigned moments with exponents `3/2,2`.
It destroys the signs and does not retain the top-positive-coordinate
distance. Consequently its moment-image results are context, not a
transported quantitative-stability theorem. The proof below is elementary
and self-contained once the displayed Newton identities are expanded.
No claim of a new general moment-image principle is made.

The decisive map used here keeps `a,b`, compresses the rest of the list
only through its square mass, and gives an upper envelope for the signed
third/fourth-moment objective. It loses the tail multiplicities and signed
first moment, but retains a one-sided inequality sufficient for (2).
The energy divisor is kept throughout; dropping it fails at one atom.

## 3. A signed tail envelope

Write

```text
A=2-sqrt(2), B=sqrt(2)-1, h=1/sqrt(2),
m=3-2sqrt(2)=1/(1+sqrt(2))^2,
f(x)=x-Ax^2,
H=p_3-Ap_4,          g=B-H.
```

Since `b<=h`, the function `f` is increasing and nonnegative on `[0,b]`.
Each positive tail entry `r_i<=b` contributes
`r_i^3-Ar_i^4=r_i^2 f(r_i)<=r_i^2 f(b)` to `H`. Every nonpositive
entry contributes at most zero, also at most `r_i^2 f(b)`. Keeping the
largest positive entry `a` separate proves the exact one-sided envelope

```text
H <= a^2 f(a)+(1-a^2)f(b).                                (3)
```

For any `0<=x<=y<=h`,

```text
f(y)-f(x)=(y-x)[1-A(x+y)] >= m(y-x),                     (4)
```

because `x+y<=sqrt(2)` and `1-A sqrt(2)=m`. Also `B=f(h)`.
Newton's identities give the exact connection to the requested ratio:

```text
J-c_* =4g/(3E).                                          (5)
```

## 4. Two cases give the same constant

First suppose `a<=h`. Rewrite (3) as

```text
g >= B-f(a)+(1-a^2)[f(a)-f(b)]
  >= m(h-a)+(m/2)(a-b)
   = (m/(2sqrt(2)))d(r)^2.                               (6)
```

Here `1-a^2>=1/2`, and (4) applies to both differences. Since `E<=1/2`,
equations (5)--(6) imply (2).

Now suppose `a>=h` and put `s=sqrt(1-a^2)`. Positive energy excludes
`a=1`, so `s>0`. We have `0<=b<=s<=h`. Set `t=a+s`, so
`1<t<=sqrt(2)`, and let `d_2^2=2-sqrt(2)t`. The exact two-atom gap is

```text
g_2=B-a^2f(a)-s^2f(s)
   = (t-1)^2(sqrt(2)-t)(A t+1)/2.                        (7)
```

At `t=sqrt(2)`, both `g_2,d_2^2` vanish. Otherwise divide (7) by
`s^2 d_2^2`, and use `(t-1)(t+1)=2as`, to obtain

```text
g_2/(s^2 d_2^2)
   = sqrt(2)a^2(A t+1)/(t+1)^2
  >= (1/sqrt(2))/(1+sqrt(2))^2
   = m/sqrt(2).                                         (8)
```

The inequality uses only `a^2>=1/2`, `A t+1>=1`, and
`t+1<=1+sqrt(2)`. Thus it holds without division at the endpoint as well.
Combining (3), (4), and (8) gives

```text
g >= g_2+s^2[f(s)-f(b)]
  >= (m/sqrt(2))s^2[d_2^2+sqrt(2)(s-b)]
   = (m/sqrt(2))s^2 d(r)^2.                              (9)
```

Finally `p_4>=a^4`, hence

```text
E=(1-p_4)/2 <=(1-a^4)/2=s^2(1+a^2)/2 <=s^2.             (10)
```

Substitute (9)--(10) into (5). This again yields (2), completing the
dimension-independent proof. Negative roots cause no missing case:
they were explicitly included in (3), while the target (1) uses only the
largest positive entries and permits zero padding.

## 5. The exponent is optimal; the global constant remains open

The quadratic distance power in (2) cannot be replaced by any smaller
positive power with a universal positive constant. To see this, use the
actual finite-degree extremizers from the inherited all-degree theorem.
For `L>=2`, put

```text
q_L=1/[L+sqrt(L(L+1)/2)],
u_L=(1,1,[-q_L repeated L]),
S_L=2-Lq_L,       r_L=u_L/S_L,
z_L=Lq_L^2/S_L^2.
```

These satisfy `p_1=p_2=1` and `e_2=0` exactly. Their two positive entries
are `sqrt((1-z_L)/2)`, their dust square mass is `z_L->0`, and their dust
third and fourth moments are respectively `o(z_L)` and `o(z_L)`.
Taylor expansion at `z_L=0` in the exact moment formulas yields

```text
d(r_L)^2 = z_L+O(z_L^2),
J(r_L)-c_* = [(28sqrt(2)-32)/3]z_L+o(z_L).               (11)
```

Thus a bound `J-c_*>=K_p d(r)^p` with `K_p>0` is impossible for `p<2`.
This determines the power, not the best constant at that power.

Unequal two-positive-atom limits give a stricter candidate constant than
(11). With `a^2+b^2=1`, `a,b>0`, and vanishing-square-mass signed dust
carrying sum `1-a-b`, the limiting ratio is determined by the two atoms.
Such limits are realized exactly: append `L` copies of `-q_L`, where
`q_L=2ab/[L(a+b+sqrt(1+2ab/L))]`, and divide by the total sum.
The cancellation quadratic is exact, and the total sum tends to one.
For `t=a+b in (1,sqrt(2))`, equations (5) and (7) give

```text
(J-c_*)/[2-sqrt(2)t]
   = [8/(3sqrt(2))](A t+1)/(t+1)^2.                     (12)
```

As `t->sqrt(2)`, this approaches
`K_two=(64-44sqrt(2))/3`, approximately `0.591534`.
This is a local two-atom obstruction, not the sharp global constant.

Indeed three equal macroscopic positive atoms give a smaller obstruction.
For `L>=2`, define

```text
q_L=2/[L+sqrt((L^2+2L)/3)],
v_L=(1,1,1,[-q_L repeated L]),
T_L=3-Lq_L.
```

The quadratic `6-6Lq_L+L(L-1)q_L^2=0` proves `e_2(v_L)=0` exactly.
After normalization by `T_L`, the three positive entries tend to
`1/sqrt(3)`, and the negative dust retains their missing first moment
while its square mass vanishes. Consequently

```text
J -> 3-4sqrt(3)/3,
d^2 -> 2-2sqrt(2/3),
(J-c_*)/d^2 ->
K_three=[3-4sqrt(3)/3-c_*]/[2-2sqrt(2/3)]
       approximately 0.3501345013 < K_two.               (13)
```

Therefore any valid best global constant is at most `K_three`, while (2)
proves it is at least `(6sqrt(2)-8)/3`, approximately `0.1617604581`.
Neither equality of the upper bound nor a complete global extremizer
classification is claimed. The failed inference is extending a local
two-atom Hessian constant to all root configurations; the missing
coordinate is the possible third macroscopic positive atom.

The separate one-atom hostile `(1,u,-u/(1+u))`, normalized by its sum,
has exact `J=1` and tends as `u->0+` to `(1,0,...)`. Its distance square
tends to `2-sqrt(2)`, so the energy boundary is compatible with a global
quadratic estimate and does not force a weaker power.

## 6. Exact controls and audit manifest

[Source](../../04-computation/quantitative_stability_empty_core_morning_sep06.py)
and [output](quantitative_stability_empty_core_morning_sep06.out).

```bash
python3 -B 04-computation/quantitative_stability_empty_core_morning_sep06.py
python3 -B -O 04-computation/quantitative_stability_empty_core_morning_sep06.py
```

All **434 explicit gates** pass. The exact rational universe starts from
all 49 prefixes `(1,x,y)` with
`x,y in {-2,-1,-1/2,0,1/2,1,2}`, tuning one final root to cancel `e_2`.
Five prefixes have zero sum and do not admit this tuning formula; one
has zero energy. All six exceptions are printed, leaving 43 eligible
signed controls. Normalization retains negative total sums correctly,
and zero padding is checked separately.

The algebraic controls consist of five actual finite-degree extremizers,
five three-positive-atom rows, six unequal-two-positive-atom rows, and
four one-atom-boundary rows. They retain all multiplicities, use exact
quadratic arithmetic and exact signs of the necessary biquadratic
expressions, compare literal ordinary-square coefficients with Newton's
identities and complete root-product energy, and cross-check compressed
binomial products by multiplication of individual factors in the bounded
small rows. The entire two-support factorization, the dust asymptotic
coefficient, and the strict three-atom/local-two-atom separation are also
checked exactly. These controls are not a census or a proof substitute.

The root agent independently audited the full analytic note: the signed
envelope, both largest-root regimes, cancellation in (8), the energy
divisor, the quadratic-power obstruction, and both two-/three-atom
constants: **PASS**. Independent referee `certificate_audit` separately
audited Sections 1--4, including the zero and signed tails, top-positive
distance, endpoint treatment and exact global constant: **PASS**. That
second referee did not audit Section 5 or replay the source; the audit
scope is recorded explicitly. Normal/optimized execution and frozen
hashes are recorded below.

```text
source SHA256 e33e95cee17e904834634633640df78146d6ca24889c07a994b4407f5e9afb63
output SHA256 2e85bb24118b3da6f5070f75d342beef10c6929aabbd91e251723c7496917d19
trace  SHA256 9a09c7b68404dccf609758ce028580b40a9913708becd671b9227b4692a0f98d
```
