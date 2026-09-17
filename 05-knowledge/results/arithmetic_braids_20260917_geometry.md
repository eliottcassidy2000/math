# Arithmetic braids: quadratic time reversal and right-triangle coordinates

**Status: PROVED elementary identities and scoped theorems + FINITE-EXACT
controls + CITED literature.** This is a research result note, not a new
canon theorem ID. No priority claim is made. Collatz, twin primes, and LRC(14)
are not consequences of the results below.

## Inheritance and concept board

The closest proved mechanism is Gaussian multiplication and primitive content
in [THM-3336, primitive Gaussian multiplication and content curvature](../../01-canon/theorems/THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md).
The directly relevant prior dynamics is [THM-4139, rational three-cycle and
order-six lift](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md)
and [THM-4146, universal three-cycle lift and fibre firewall](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md).
In particular, THM-4146 already proves the `3:4:5` signed-cycle rigidity; it
must not be presented as a discovery of this session.

The canonical hostile is `(3,4,5) star (7,24,25)=25(-3,4,5)`: general
primitive Gaussian multiplication needs a content sidecar. The corrected
near miss is THM-4139/4146's failed identification of a source three-point
divisor with one target elliptic fibre. The least-used relevant sidecar here
is the **chosen cycle trace**: a quadratic has two period-three cycles
generically, and its parameter alone does not identify the cycle.

| Lane | Object / representation | Predicate | Lost coordinate / next probe |
|---|---|---|---|
| Anchor | A chosen quadratic three-cycle / trace | Transport dynamics under reversal | Keep which cycle; test both Chebyshev cubics |
| Niche | Primitive right triangle / signed Gaussian point | New prime divisors under squaring | Keep the leg, hypotenuse, sign, and denominator separately |
| Wildcard | Semicircle foot / `(e/d,l,theta)` | Relate continuous shape to arithmetic | Keep scale and rationality; test `3-4-5` literally |
| Bridge | Seventh-root traces / doubling | Explain `-7/4` through a genuine map | Equality on the cycle is not global conjugacy |
| Hostile | Critical orbit / third numerator | Distinguish repetition from absent new primes | Specify start, denominator, and primitive-divisor convention |

## 1. Affine reversal of every quadratic three-cycle

**PROVED.** Work over any characteristic-zero field containing a chosen
cycle of three distinct points of `f_c(x)=x^2+c`. Write its trace as `sigma`.
Define

```text
L(x) = -x-1/2,
sigma' = -sigma-3/2,
c' = c-sigma-3/4.                                      (1)
```

Then `L` carries the **reversed** cycle of `f_c` to a forward cycle of
`f_(c')`. Applying this operation twice returns the original marked cycle
and parameter. Among nonconstant affine maps, `(1)` is the unique such
reversal into a centered monic quadratic.

**Proof.** On the three-point cycle, the sum of its three iterates gives

```text
f_c^2(x) = sigma-x-f_c(x) = sigma-x-x^2-c.               (2)
```

Substituting into `f_(c')(L(x))=L(f_c^2(x))` proves (1). For uniqueness,
replace `L` by `a x+b`, `a!=0`. Modulo the cycle polynomial, the difference
between the two sides is a polynomial of degree at most two. It vanishes
at three distinct points and therefore is zero. Its coefficients give

```text
a^2=-a,       2ab=-a,       b^2+c'=a(sigma-c)+b.
```

Consequently `a=-1`, `b=-1/2`, and `c'=c-sigma-3/4`. The transformed trace
is `sigma'=-sigma-3/2`, so `c'-sigma'-3/4=c`. Also `L^2=id`. QED.

This transports exact period, field of definition, and splitting of the
cycle polynomial. It does **not** assert a conjugacy of the two entire
quadratic maps: the source inverse is represented by `f_c^2` only on the
cycle. Its global degree is four.

THM-4146 Section 2 supplies the trace formulas

```text
c = -(sigma^2+sigma+2),
P_sigma(x) = x^3-sigma*x^2-(sigma^2+2sigma+3)*x
             +(sigma^3+2sigma^2+3sigma+1).              (3)
```

In the especially economical coordinate `eta=2sigma+1`, these become

```text
c = -(eta^2+7)/4,
time reversal: eta -> -eta-1,
the other cycle at the same c: eta -> -eta,
disc(P_sigma) = (eta^2+eta+7)^2.                         (4)
```

The two operations in (4) are different reflections. The first transports
an actual chosen cycle by `L`; the second selects the other cubic factor
of the third dynatomic polynomial. The second need not preserve rational
splitting. Composing them translates `eta` by one, but this does not generate
rational cycles automatically. The discriminant is invariant under time
reversal, as it must be under an affine reflection of the three roots.

### Why `-29/16` is the fixed parameter of this operation

Keeping the chosen parameter fixed requires `c'=c`, hence

```text
sigma=-3/4,       eta=-1/2,       c=-29/16.               (5)
```

Here the cycle is

```text
-7/4 -> 5/4 -> -1/4 -> -7/4,
```

and `L` interchanges `-7/4` and `5/4` while fixing `-1/4`. Thus this is
the unique fixed chosen-cycle trace of affine time reversal. This recovers
THM-4139/4146's unique arithmetic-progression cycle through a different
operation: a three-element set invariant under a reflection consists of its
center and one reflected pair.

This does not mean every cycle at `c=-29/16` is fixed: its other trace is
`-1/4`, and that chosen cycle reverses to the parameter `-37/16`.

### The rational chart has a genuine three-by-two symmetry

For a marked cycle `alpha -> beta -> gamma -> alpha`, put

```text
t=alpha+beta=(beta-gamma)/(alpha-beta),    t!=0,-1.
```

Subtracting the cycle equations gives the other two pair sums
`beta+gamma=-(t+1)/t` and `gamma+alpha=-1/(t+1)`. Solving these three
linear equations recovers the complete chart from THM-4139:

```text
p0(t)=(t^3+2t^2+t+1)/(2t(t+1)),
p1(t)=(t^3-t-1)/(2t(t+1)),
p2(t)=-(t^3+2t^2+3t+1)/(2t(t+1)),
sigma(t)=(t^3-3t-1)/(2t(t+1)),
eta(t)=(t^3+t^2-2t-1)/(t(t+1)).                          (5a)
```

Conversely these form an exact three-cycle for every rational `t!=0,-1`;
their pairwise differences have nonzero factor `t^2+t+1`. Thus the chart
is complete for rational marked cycles, not merely a construction.

The operations (with `rho` a backward cyclic change of marking)

```text
rho(t)=-1/(t+1),              tau(t)=1/t
```

satisfy `rho^3=tau^2=id` and `tau rho tau=rho^(-1)`. Their exact action is

```text
(p0,p1,p2)(rho(t))=(p2,p0,p1)(t),
(p0,p1,p2)(tau(t))=(L(p2),L(p1),L(p0))(t).                (5b)
```

Hence this is an actual `S3` action on chart parameters: three cyclic
markings of a cycle, paired with three markings of its affine-reversed
partner. It is not a tournament, and the partner generally lives at a
different quadratic parameter. Passing to cyclic markings leaves the
two-way reversal. The sole rational fixed cycle follows from

```text
eta(t)+1/2=(t-1)(t+2)(2t+1)/(2t(t+1)),                  (5c)
```

whose roots `{1,-2,-1/2}` are the three markings of the AP cycle at
`-29/16`. The distinguished starting point of a marked cycle need not
be fixed; the fixed object here is the cycle modulo cyclic remarking.

The numerator of `eta(t)` is exactly the seventh-root trace polynomial
used next. This gives another precise junction: the parabolic branch
`eta=0` occurs at cubic, rather than rational, chart parameters.

## 2. The exact seventh-root bridge to `-7/4`

Put `D(t)=t^2-2`. Its third dynatomic polynomial factors as

```text
(D^3(t)-t)/(D(t)-t)
 = (t^3+t^2-2t-1)(t^3-3t+1) = q7(t) q9(t).             (6)
```

The three roots of `q7` are `2cos(2pi/7)`, `2cos(4pi/7)`, and
`2cos(6pi/7)`. Indeed, if `z` is a nontrivial seventh root of unity, dividing
`1+z+...+z^6=0` by `z^3` and writing `t=z+z^(-1)` gives `q7(t)=0`.
The identity

```text
D(z+z^(-1)) = z^2+z^(-2)                               (7)
```

identifies `D` with exponent doubling after quotienting `z` by inversion.
These three distinct real roots make one exact three-cycle.

**PROVED transport.** The trace of this cycle is `sigma=-1`, so (1) gives
`c'=-7/4` and `sigma'=-1/2`. More explicitly,

```text
D^2(t) = -t^2-t+1                 modulo q7(t),
f_(-7/4)(L(t))-L(D^2(t)) = (t-1) q7(t).                 (8)
```

Therefore `L(t)=-t-1/2` turns this doubling cycle into the reversed
three-cycle of `x^2-7/4`. Its cycle polynomial is

```text
P(x) = x^3 + x^2/2 - 9x/4 - 1/8,
disc(P)=49,
8 * (product of its roots) = 1.                         (9)
```

The last quantity is the derivative multiplier of the cycle. Thus the
cycle is parabolic, and the full third dynatomic polynomial is exactly
`P(x)^2`. The two generic three-cycle traces collide at `eta=0`, explaining
why `-7/4` is also the branch parameter in (4).

The rational roots test makes `q7` irreducible. Its cycle and its image
are consequently cubic algebraic cycles, unlike the rational cycle at
`-29/16`.

**Hostile companion.** For the other factor in (6),

```text
D^2(t) = -t^2-t+2                 modulo q9(t).
```

Its trace is zero. The same reversal therefore gives parameter `-11/4`,
with cycle multiplier **19**, not 1. This identifies exactly which
three-cycle makes the `-7/4` bridge work; being a period-three doubling
orbit alone is insufficient.

The arithmetic `63=3^2*7=2^6-1` meets this bridge at an actual shared
coordinate: `ord_7(2)=3`. But no equality of all three dynamical systems
follows. As already shown in THM-4139 Section 6, the factors at exponent
six combine orders two and three and increased 3-adic depth; absence of a
new prime is not absence of new order. Repeating the numeral 7 or 63
without the map (7) would not establish a connection.

## 3. Critical orbits and the sign audit

The user's small integer graphs correspond to `x^2`, `x^2-1`, and
`x^2-2`, with minus signs in the latter two:

```text
x^2:     0 -> 0;            -1 -> 1 -> 1.
x^2-1:   1 -> 0 -> -1 -> 0.
x^2-2:   0 -> -2 -> 2 -> 2;  1 -> -1 -> -1.
```

**PROVED classification.** These are exactly the rational parameters
`c` for which the critical orbit of `x^2+c`, starting at zero, is finite.
For nonintegral `c=a/b` in lowest terms, the denominator of the nth
iterate is `b^(2^(n-1))`, with no cancellation, so it grows strictly.
For integer `c>=1` the orbit grows, while for `c<=-3` its second term is
at least 6 and exceeds `|c|`, after which it grows. The three remaining
integer parameters are checked above.

For `c=-7/4`, the critical orbit begins

```text
0 -> -7/4 -> 21/16 -> -7/256 -> ... .                    (10)
```

Its third **numerator**, in lowest terms, has no primitive prime divisor.
The rational value does not repeat and the critical orbit is infinite.
This is a different cycle from the algebraic parabolic orbit in (9).

There is a useful exact small-index arithmetic test. For nonintegral
`c=a/b`, `b>=2`, in lowest terms,

```text
f_c^3(0) = a Q(a,b)/b^4,
Q(a,b) = a^3+2a^2*b+a*b^2+b^3,
gcd(Q(a,b), a(a+b)b)=1.                                 (11)
```

The gcd follows by reducing `Q` modulo `a`, `a+b`, and `b`. Thus every
prime dividing `Q` is new at the third numerator, and the third numerator
has no new prime **if and only if** `Q(a,b)=1` or `-1`. At `(-7,4)`,
`Q=1`. An exact census over reduced `-2<a/b<-1`, `2<=b<=2000` finds
`(-7,4)` as its sole solution; this finite census is not an all-height
solution of the cubic Thue equations.

**CITED boundary.** Holly Krieger's
[*Primitive prime divisors in the critical orbit of z^d+c*](https://www.dpmms.cam.ac.uk/~hk439/PPD.pdf),
Theorem 1.1, bounds the number of exceptional indices by 23 for infinite
critical orbits. Theorem 6.1 gives the bound 3, with the third-index exception
exactly at `-7/4`, on its explicitly defined set of parameters outside
specified attracting-basin regions. The additional hypotheses matter;
they have not been suppressed or replaced here by a universal bound of 3.
The paper was read directly on 2026-09-17. Equations (10)-(11) need none
of its deeper machinery.

For the affine iteration `u -> 2u+1` starting at zero, `u_n=2^n-1`.
The first term `u_1=1` also has no prime; among indices `n>1`, 6 is the
classical unique Bang-Zsigmondy exception, as stated in the introduction of
[Zhi-Wei Sun's primary research paper](https://maths.nju.edu.cn/~zwsun/195g.pdf).
The immediate mechanism at
6 is `63=3^2*7`, with the primes already present at indices 2 and 3.
This is separate from the numerator unit in (11).

## 4. Primitive triples: the odd square is the even-leg sum

For a positive primitive Pythagorean triple, take

```text
A=m^2-n^2 (odd),  B=2mn (even),  C=m^2+n^2,
m>n>0, gcd(m,n)=1, m-n odd.
```

Then

```text
C+B=(m+n)^2,    C-B=(m-n)^2,                             (12)
C+A=2m^2,       C-A=2n^2.
```

Thus the user's square observation is correct with the **even leg**, and
both signs give odd squares. One half alone does not identify a triple:
`C+B=u^2` permits many choices. The faithful pair is

```text
u=m+n, v=m-n,   C=(u^2+v^2)/2,
B=(u^2-v^2)/2, A=uv,
```

where `u>v>0` are odd and coprime. This retains the missing coordinate.

THM-4146 Section 4.2 already gives the strongest direct bridge to `29`:
if positive `a,b,h` satisfy `a^2+b^2=h^2` and the signed cycle

```text
-(a+b) -> h -> -(b-a) -> -(a+b)
```

is generated by `(y^2-D)/b`, then necessarily

```text
(a,b,h)=(3k,4k,5k),  D=29k^2.
```

The present reversal theorem explains the same special parameter as an
involution's fixed point. The two mechanisms agree, rather than being
independent numerical coincidences.

## 5. Gaussian squaring: fixed hypotenuse support, fresh leg primes

**PROVED.** Start with a nondegenerate signed primitive triple
`A_0^2+B_0^2=C_0^2`, with `A_0` odd, `B_0` even, `C_0>0`. Define

```text
A_(n+1)=A_n^2-B_n^2,
B_(n+1)=2A_n B_n,
C_(n+1)=C_n^2.                                         (13)
```

Every iterate is again primitive of these parities, and

```text
C_n=C_0^(2^n),
B_n=2^n B_0 product_(j<n) A_j,
gcd(A_i,A_j)=1 for i!=j,
|A_j|>1 for every j.                                    (14)
```

**Proof.** The Pythagorean identity in (13) follows from Gaussian squaring.
Any odd prime common to `A^2-B^2` and `2AB` would divide both `A` and `B`;
2 cannot divide the odd first leg. Thus primitivity is preserved. The first
two formulas in (14) follow by induction. Every earlier `A_i` divides
`B_j`, so its gcd with `A_j` is one. A nondegenerate integer right triangle
cannot have leg of absolute value one, since
`1=C^2-B^2=(C-|B|)(C+|B|)` would force `B=0`. No iterate becomes degenerate:
`A^2=B^2` would force `C^2=2A^2`, impossible for nonzero integers. QED.

Consequently each odd leg contributes a prime absent from **all earlier
legs and the fixed hypotenuse prime support**. The hypotenuse itself
introduces no new prime after the initial triple. This is an exact, all-depth
separation between support growth and coordinate growth.

For example,

```text
(3,4,5) -> (-7,24,25) -> (-527,-336,625)
        -> (164833,354144,390625) -> ... .               (15)
```

The elementary gcd proof is stronger for this purpose than factoring a
long list of outputs. General multiplication does not inherit (14);
THM-3336's content hostile above explains the boundary of this theorem.

Under

```text
z_n=(A_n+i B_n)/C_n,       x_n=2A_n/C_n,
```

we have `z_(n+1)=z_n^2` and exactly

```text
x_(n+1)=x_n^2-2.                                       (16)
```

This is the promised actual right-triangle/quadratic map. Signs must be
retained for the displayed conjugacy; sorting the two positive leg lengths
at every step folds the angle and changes the scalar rule.

The reduced denominator of `x_n` is `C_0^(2^n)`. Hence no nondegenerate
primitive-triangle orbit is preperiodic under (16), even though all its real
coordinates lie in `[-2,2]`. Bounded real motion is not bounded arithmetic
height. More generally the only rational preperiodic points of `x^2-2`
are `{-2,-1,0,1,2}`: noninteger denominators square, and integers outside
this set escape.

## 6. Semicircle coordinates and the standard square-plus-one family

Normalize a right triangle's hypotenuse to one and let its shorter and
longer legs have lengths `a<=b`, so `a^2+b^2=1`. The altitude to the
hypotenuse has length `l=ab`; the foot splits the hypotenuse into
`e=a^2` and `d=b^2`. The angle between the altitude and the short leg is
`theta=arctan(a/b)`, in `(0,pi/4]`. Therefore

```text
r=e/d,
e=r/(1+r),   d=1/(1+r),
l=sqrt(r)/(1+r),
theta=arctan(sqrt(r)),
e=sin(theta)^2, d=cos(theta)^2, l=sin(2theta)/2.          (17)
```

These are three coordinates on **one** shape interval, not independent
parameters. They tend to `(1,1/2,pi/4)` at the isosceles boundary and to
`(0,0,0)` at degeneration. No nondegenerate integer Pythagorean triangle
attains the isosceles endpoint.

The correct folded dynamics can also be retained explicitly. Gaussian
squaring doubles the leg angle before sorting, so on this small-angle
coordinate its operation is the tent map

```text
theta' = min(2theta, pi/2-2theta),
s=4theta/pi:     s'=2 min(s,1-s).                        (17a)
```

Equivalently, if `t=tan(theta)=a/b`, its positive update is
`t'=min(2t/(1-t^2),(1-t^2)/(2t))`. This rational identity is checked on
every triple in the exact census. The folded shape operation is therefore
fully determined, but a sheet/quadrant sidecar is needed to restore the
signed Gaussian point. The primitive-triple shapes in this tent map never
repeat, because their reduced integral hypotenuses grow by squaring.

With endpoints `(0,0),(1,0)`, the third vertex is `(e,l)` and

```text
(e-1/2)^2+l^2=1/4.                                    (18)
```

Thus the hypotenuse is the semicircle's **diameter**; the radius is `1/2`.
Primitive triples give a countable dense set of shapes on this interval,
because rational Euclid slopes are dense and the shape map is continuous.
This is a statement about all triples, not density of one squaring orbit.

The square-minus/plus-one family is

```text
(k^2-1,2k,k^2+1),  k>=2.
```

For even integer `k` it is primitive; for odd integer `k` divide all sides
by 2. Its true unnormalized altitude and normalized altitude are

```text
h=2k(k^2-1)/(k^2+1),
l=2k(k^2-1)/(k^2+1)^2.                                (19)
```

At `k=2`, `(3,4,5)` has `h=12/5`, `l=12/25`, `r=9/16`, and
`theta=arctan(3/4)`. Its altitude is not `sqrt(2)`.

The user subsequently clarified that this standard family is intended.
It has a sharp shape restriction when `k` ranges over **integers**:

```text
max_(integer k>=2) l(k)=12/25,
with equality precisely at k=2 and k=3.                 (20)
```

Indeed differentiation gives

```text
l'(k)=-2(k^4-6k^2+1)/(k^2+1)^3,
```

so `l` increases up to `k=1+sqrt(2)` and decreases thereafter; direct
substitution gives `l(2)=l(3)=12/25`. These two parameters both reduce to
the `3-4-5` similarity class. More generally the only two parameter values
with the same unordered normalized shape are

```text
k  and  (k+1)/(k-1).                                   (21)
```

To see this, use the strictly increasing unsorted leg ratio
`(k^2-1)/(2k)`: equal unordered shapes give either equal ratios or reciprocal
ratios. Solving the latter gives (21). For integral `k>=2`, (21) is integral
only when `k-1` divides 2. Hence `k=2,3` is the unique repeated integer
shape.

The integer family consequently is not dense in the semicircle: its
normalized altitudes never enter `(12/25,1/2]`, and its only accumulation
shape is the degenerate one as `k` tends to infinity. In contrast the set
of all primitive triples is dense in the shape interval. Allowing rational
`k>1` in the same formula recovers all rational right-triangle shapes,
because `k=m/n` is the Euclid parameter after scaling. Thus it is the
integer restriction, not the square-plus-one formula itself, that creates
the thin subfamily.

There is a sharp obstruction to the literal proposed family: every positive
primitive integer right triangle has altitude `AB/C` with reduced denominator
exactly `C`, because `gcd(AB,C)=1`. It is rational but not an integer. If
`k` is an integer, `sqrt(k)` is either irrational or an integer. Therefore
**no positive primitive integer right triangle has altitude `sqrt(k)` for
integer `k`**, regardless of the hypotenuse condition. A different altitude
or scaling would be needed for that literal condition; the clarified
standard family instead has (19).

## 7. Connection ledger and verification

| Source -> target | Map | Preserved predicate | Lost information / sidecar |
|---|---|---|---|
| Chosen quadratic three-cycle -> reversed quadratic cycle | `L=-x-1/2`, `(1)` | Exact period, splitting, field, cubic discriminant | Parameter alone loses chosen trace; global quadratic dynamics is not preserved |
| Seventh roots -> Chebyshev three-cycle | `z -> z+z^-1` | Doubling after inversion quotient | Distinguish `z` from `z^-1` with a sheet sign |
| Primitive triple squaring -> `x^2-2` | `x=2A/C` | Exact forward scalar dynamics | Imaginary sign; keep `B/C` to reconstruct circle point |
| Triangle -> unit semicircle | `(A,B,C) -> (A^2/C^2,AB/C^2)` | Shape and altitude | Scale, integer content, leg label |
| Critical third numerator -> unit equation | `(11)` | Iff absence of a new numerator prime | Parameter height; finite census is not a Thue classification |

The concept board changed in two useful ways. First, `-2`, `-7/4`, and
`-29/16` now sit on one explicit chosen-cycle moduli operation, rather than
being linked by matching numerals. Second, primitive primes behave in
opposite ways in the two triangle coordinates: the same squaring operation
freezes hypotenuse support and forces fresh leg support. This is exactly why
an unlabelled prime-support graph cannot transport a conclusion by itself.

Reproduce with

```text
python 04-computation/experiments/arithmetic_braids_20260917_geometry.py
python -O 04-computation/experiments/arithmetic_braids_20260917_geometry.py
```

The script uses only standard-library exact integers and fractions. It checks
polynomial quotient identities by coefficient arithmetic independently of
the trace proof, 181 rational marked parameters yielding 91 distinct cycles,
all 1,314 primitive Euclid pairs with `m<=80` through five squaring stages,
the standard family at all `2<=k<=1000`, and all 1,216,587 reduced parameters
in `-2<a/b<-1`, `2<=b<=2000` for the third-numerator test. Positive controls
are the rational AP cycle and
`3-4-5`; hostile controls include the ninth-root cubic and the altitude
denominator. Normal and optimized execution both retain all checks.

An independent agent rederived the affine uniqueness coefficients, the
all-depth fresh-prime proof, and the rational-chart `S3` action and reported
no substantive mathematical issue. Its suggested clarification that `rho`
changes the marking backward has been included.

Final source SHA-256 (working-tree bytes):
`e7fc46feab6f343744dfb9c1a6c9565840cd769c315a4b6b68368e3c70f43f4c`.
The normal and optimized transcripts agree:

```text
polynomial quotient identities: PASS (q7 multiplier 1; q9 multiplier 19)
rational chart: 181 markings, 91 cycles; reversal PASS
chart reciprocity/rotation S3: PASS; fixed orbit [-2, -1/2, 1]
primitive Euclid universe m<=80: 1314 triples; 6570 depth checks PASS
signed triple square orbit: [(3, 4, 5), (-7, 24, 25), (-527, -336, 625), (164833, 354144, 390625)]
integer family 2<=k<=1000: maximum 12/25 at k=2,3; sole duplicate (2,3); PASS
critical orbit at -7/4: ['-7/4', '21/16', '-7/256']
FINITE-EXACT -2<a/b<-1, 2<=b<=2000: 1216587 parameters; units [(-7, 4, 1)]
All checks passed.
```

The next honest frontier is arithmetic, not an analogy: study integral
height and primitive-prime transport under the exact reciprocity action
(5b), or solve the explicit unit equations (11) at all heights with a
certified Thue solver/proof. Neither has been claimed solved by the bounded
computation. The rational-cycle chart and its reversal action themselves
are completely described above.
