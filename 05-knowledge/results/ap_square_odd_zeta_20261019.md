# Odd zeta values from AP square-sum triangles: the exact height series and its angular correction

**2026-10-19. Status: PROVED synthesis, with FINITE-EXACT checks and an
independent mathematical audit.** The Gaussian factor-count and special-function
identities are classical. No priority claim or new irrationality result is made.
The contribution is the precise bridge, its scale/axis/orientation normalization,
and an explicit angular correction for the positive-label subfamily.

## 1. The bridge, and what it does not identify

The starting classification is
[AP square-sum triangles](ap_square_pythagorean_20261019.md), Theorem 1.
Its closest existing arithmetic mechanism is
[THM-3334, Gaussian hypotenuse-fibre torsors](../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md),
§5: a fixed primitive hypotenuse has one conjugate Gaussian factor choice per
split prime. The current note reconnects this arithmetic with the original
questions about the odd values `ζ(5), ζ(7), ζ(9)`.

The live concepts are: AP square-sum labels; primitive Gaussian norms;
root versus label height; odd zeta/beta values; and the angular positivity filter.
The decisive hostile is the primitive triple `(3,4,5)`: it contributes to the
unrestricted Dirichlet series but gives a negative smallest graph label.

There is an exact, elementary statement:

\[
\boxed{\displaystyle
P(s):=\sum_{\substack{A>B>0,\ A^2+B^2=C^2\\\gcd(A,B,C)=1}} C^{-s}
=\frac12\left[\frac{\zeta(s)\beta(s)}{(1+2^{-s})\zeta(2s)}-1\right],
\qquad \Re s>1,}                                      \tag{1}
\]

where `β(s)=Σ_{j≥0}(-1)^j/(2j+1)^s`. This is also the root-height Dirichlet
series of the **base signed AP square-sum triangles** defined below. In particular,

\[
\begin{aligned}
P(5)&=\frac{4725}{32}\frac{\zeta(5)}{\pi^5}-\frac12,\\
P(7)&=\frac{8243235}{5504}\frac{\zeta(7)}{\pi^7}-\frac12,\\
P(9)&=\frac{334084125375}{22459904}\frac{\zeta(9)}{\pi^9}-\frac12.
\end{aligned}                                         \tag{2}
\]

Thus the odd-zeta values are recovered exactly from weighted counts of these
integer-labelled graph configurations. The map preserves the full triple and
the scale. Passing to `C` alone loses the Gaussian factor allocation and leg
angle; the angle is precisely the extra information needed for positivity.

This weighted-count realization is a different construction from the moduli-space
periods in [the odd-zeta/Petersen note](odd_zeta_boundary_graphs_20261019.md).
No identification of the two geometries, period chains, or irrationality criteria
is asserted. Formula (2) does not evaluate `P(5),P(7),P(9)` by simpler constants;
it supplies a faithful arithmetic realization of the zeta values.

## 2. Signed labels make the Pythagorean correspondence unrestricted

A **signed AP square-sum triangle** is an increasing integer triple `x<y<z`
such that `x+z=2y` and all three pair sums are strictly positive perfect squares.
“Triangle” refers to the three edges of a square-sum graph; the labels are not
Euclidean side lengths. In particular, a negative label is permitted here.

### Proposition 1: complete signed classification

Every such triple is uniquely

\[
(x,y,z)=2g^2(C^2-4AB,\ C^2,\ C^2+4AB),                \tag{3}
\]

with `g≥1` and one primitive Pythagorean triple with size-ordered positive legs
`A>B`. Conversely every choice in (3) gives a signed AP square-sum triangle.
Its pair-sum roots and root height are

\[
(p,q,r)=2g(A-B,C,A+B),\qquad
h=\sqrt{y/2}=gC.                                      \tag{4}
\]

Every pair of labels has gcd `2g²`; every pair of pair-sum roots has gcd `2g`.
The **base** case is `g=1`, so its label gcd is two, not one. The positive-label
subfamily is exactly the strict cone `C²>4AB`.

**Proof.** Let `p,q,r` be the positive roots of `x+y,x+z,y+z` in increasing
order. Then `p²+r²=2q²` and `q²=2y`. Since `y` is an integer, `q` is even;
the first equation modulo four forces `p,r` even. Set
`(u,v,w)=(p,q,r)/2`. The integers `u,w` have the same parity, so

\[
P_0=(w+u)/2,\quad Q_0=(w-u)/2,\quad P_0^2+Q_0^2=v^2
\]

are positive integers with `P_0>Q_0`. Divide by
`g=gcd(P_0,Q_0,v)` to obtain `(A,B,C)`. Solving the pair-sum equations gives
(3) and (4), and every step is reversible. The converse follows by expansion.
Primitive Pythagorean legs are coprime with opposite parity, and `C` is odd.
The three cores `C²−4AB,C²,C²+4AB` are pairwise coprime, as are
`A−B,C,A+B`; their only possible common odd prime would divide both `C` and
`AB`. This proves the gcd assertions, including for negative labels. Finally,
`y,z>0` automatically and `x>0` precisely in the stated cone. QED.

There is no zero-label boundary: `C²=4AB` is impossible since `C` is odd.
The first unrestricted base configuration is

\[
(A,B,C)=(4,3,5)\longmapsto(-46,50,146),
\quad(p,q,r)=(2,10,14).
\]

The first positive one by root height is
`(40,9,41) ↦ (482,3362,6242)`. This also uniquely minimizes the largest
positive label, by the predecessor's unbounded minimum proof.

## 3. Gaussian multiplicity, the axis, and the Euler product

Let `r(C)` count size-ordered primitive Pythagorean triples of hypotenuse `C`.
Then

\[
r(C)=\begin{cases}
2^{\omega(C)-1},&C>1\text{ and every prime divisor of }C\text{ is }1\pmod4,\\
0,&\text{otherwise}.
\end{cases}                                           \tag{5}
\]

Here `ω` counts distinct primes, without their exponents. For example,
`r(5)=r(25)=1`, `r(65)=2`, `r(1105)=4`, while `r(1)=r(10)=r(15)=0`.

**Derivation.** Euclid's unique parameters satisfy
`C=n²+m²`, `0<m<n`, `gcd(m,n)=1`, and `m+n` odd. Factor `n+im` in the
Euclidean domain `Z[i]`. Primitivity excludes inert primes `q≡3 mod4`;
opposite parity excludes the ramified prime over two. For each
`p^e || C`, with `p≡1 mod4`, all `e` copies must be assigned to one of the
two conjugate Gaussian prime factors, or a rational `p` would divide both
coordinates. There are `2^{ω(C)}` choices and four Gaussian units. Quotienting
the eight sign/interchange symmetries gives `2^{ω(C)−1}` for `C>1`.

At `C=1`, the only representations are the four axes. They have a different
stabilizer and do not describe a nondegenerate triple. One must subtract this
case before halving. Equivalently,

\[
1+2P(s)
=\prod_{p\equiv1\ (4)}\left(1+2\sum_{e\ge1}p^{-es}\right)
=\prod_{p\equiv1\ (4)}\frac{1+p^{-s}}{1-p^{-s}}.       \tag{6}
\]

The factor of two does not count a second graph triangle obtained by swapping
the legs: (3) is symmetric in `A,B`, and size ordering has already identified
that symmetry.

For `Re(s)>1`, compare prime factors in `ζ(s)β(s)/ζ(2s)`. A split prime has
factor `(1+p^−s)/(1−p^−s)`; an inert prime has factor one. The factor at two
is `1+2^−s`. Removing it gives (1). All these products and rearrangements
converge absolutely in the stated half-plane.

## 4. Why odd arguments 5, 7, 9 give exact rational coefficients

With Euler numbers defined by `sech(t)=Σ E_j t^j/j!`, the classical identities are

\[
\beta(2k+1)=\frac{|E_{2k}|\pi^{2k+1}}{4^{k+1}(2k)!},
\qquad
\zeta(2s)=\frac{|B_{2s}|(2\pi)^{2s}}{2(2s)!}.           \tag{7}
\]

For the first formula, set `x=1/2` in the Euler-polynomial Fourier series
[DLMF 24.8.4](https://dlmf.nist.gov/24.8.E4), and use
[DLMF 24.2.9](https://dlmf.nist.gov/24.2.E9),
`E_{2k}(1/2)=2^−2k E_{2k}`. The second formula is
[DLMF 25.6.2](https://dlmf.nist.gov/25.6.E2).

Consequently, for every odd integer `s≥3`,

\[
P(s)=K_s\frac{\zeta(s)}{\pi^s}-\frac12,
\qquad
K_s=\frac{|E_{s-1}|(2s)!}
{2^{2s+1}(2^s+1)(s-1)!|B_{2s}|}.                      \tag{8}
\]

The inputs for (2) are `E_4=5,E_6=−61,E_8=1385` and
`B_10=5/66,B_14=7/6,B_18=43867/798`. The reproducibility script derives
these by rational recurrences and simplifies (8), rather than storing the
answers as floating-point fits.

### The same Eisenstein family, with a different observer

The [E3/E4/E5 atlas](e_atlas_20261019.md), §2.1, uses
`E(z,s)=(1/2)Σ_{gcd(c,d)=1} y^s/|cz+d|^{2s}`. At the Gaussian CM point
`z=i`, the full nonzero Gaussian lattice sum is `4ζ(s)β(s)`: each nonzero
Gaussian ideal has four generators, all of the same norm, and its ideal series
is `ζ(s)β(s)`. Separating the positive common divisor of `(c,d)` contributes
`ζ(2s)`. Therefore, for `Re(s)>1`,

\[
E(i,s)=\frac{2\zeta(s)\beta(s)}{\zeta(2s)},\qquad
\boxed{P(s)=\frac{E(i,s)}{4(1+2^{-s})}-\frac12.}        \tag{8a}
\]

This connects the PPT height series to the **same real-analytic Eisenstein
family** as the atlas. The observer differs: (8a) evaluates at `i`, while the
atlas integrates along a horocycle to obtain the constant Fourier coefficient
`y^s+φ(s)y^{1−s}`, with `φ(s)` containing `ζ(2s−1)`. For example,
`E(i,3)=945ζ(3)/(16π³)`, whereas `φ(3)=2835ζ(5)/(8π⁵)`.
Thus the CM evaluations recovering `ζ(5),ζ(7),ζ(9)` have parameters `5,7,9`;
the atlas scattering parameters are `3,4,5`. Point evaluation and the constant
Fourier coefficient must not be interchanged. The script independently checks
the full and primitive Gaussian lattice coefficients, including both axis and
norm-two terms, to test all factors in (8a).

## 5. Scaling: the exponent depends on the actual height

The unique square-scale decomposition yields the following distinct series:

\[
\begin{aligned}
\text{base root height:}\quad &\sum C^{-s}=P(s),&&\Re s>1,\\
\text{all root heights:}\quad &\sum h^{-s}=\zeta(s)P(s),&&\Re s>1,\\
\text{all normalized middle labels:}\quad
&\sum(y/2)^{-s}=\zeta(2s)P(2s),&&\Re s>1/2,\\
\text{all literal middle labels:}\quad
&\sum y^{-s}=2^{-s}\zeta(2s)P(2s),&&\Re s>1/2.
\end{aligned}                                         \tag{9}
\]

Each sum is over triangles, retaining multiplicity when heights coincide.
The corresponding base middle-label series is `P(2s)` (or `2^−s P(2s)`).
For the separate scale variable `t=g²`, the scale-only series is `ζ(2s)`.
A two-variable bookkeeping form is

\[
\sum_{\text{base triples}}\sum_{g\ge1}C^{-u}g^{-v}
=P(u)\zeta(v),\qquad\Re u>1,\ \Re v>1.                \tag{10}
\]

Thus `ζ(2s)P(s)` is valid if one deliberately weights the hybrid quantity
`g²C`; it is neither the actual root-height series nor the middle-label series.

There is also an exact count at each root height. If `a(h)` is the number of
signed triangles of height `h`, then

\[
a(h)=\sum_{C\mid h}r(C)
=\frac12\left[\prod_{p\equiv1\ (4)}(2v_p(h)+1)-1\right].               \tag{11}
\]

The empty product is one. For instance, height `65` has two primitive shapes
and the rescaled shapes from `C=5` and `C=13`, giving four triangles.
Inert primes and powers of two can occur in `h` through `g`, despite being
forbidden in primitive `C`.

## 6. Counts and convergence boundaries

Let `N(X)=Σ_{C≤X}r(C)`. Counting coprime opposite-parity lattice points in the
Euclid sector `0<arg(n+im)<π/4`, `n²+m²≤X`, gives

\[
N(X)=\frac{X}{2\pi}+O(\sqrt X\log(2X)).                \tag{12}
\]

For completeness, opposite-parity points have density `1/2`; Möbius inversion
over odd common divisors changes this to
`(1/2)Σ_{d odd} μ(d)/d²=4/π²`. The sector area is `πX/8`.
Summed boundary errors are `O(√X log(2X))`, and the tail of the area sum is
`O(√X)`. This gives (12) without a Tauberian assumption.

All square scales counted by root height therefore satisfy

\[
\sum_{h\le X}a(h)
=\sum_{g\le X/5}N(X/g)
=\frac{X\log X}{2\pi}+O(X).                           \tag{13}
\]

Indeed, the main term is a harmonic sum, and the summed errors are bounded by
`√X Σ_{g≤X/5}g^−1/2 log(2X/g)=O(X)`. Partial summation establishes the
absolute convergence domains in (9). The primitive and all-root series diverge
at real `s=1`, and the middle-label series at real `s=1/2`; hence those are the
exact abscissae of absolute convergence. In particular, the condition is not
`Re(s)>2` merely because Euclid parameters have two coordinates: the height
`C=n²+m²` is quadratic in their radius.

As a consistency check, (1) has residue

\[
\mathop{\rm Res}_{s=1}P(s)
=\frac{\beta(1)}{2(1+1/2)\zeta(2)}=\frac1{2\pi},       \tag{14}
\]

and `ζ(s)P(s)` has leading double-pole coefficient `1/(2π)`.

## 7. Positivity is an angular condition, not a factor of one third

Define

\[
P_+(s)=\sum_{\substack{\text{primitive PPTs}\\C^2>4AB}}C^{-s}.
\]

In Euclid angle `φ=arg(n+im)∈(0,π/4)`, positivity means
`φ∈(0,π/24)∪(5π/24,π/4)`. These sectors have one third of the full angular
width. Consequently

\[
N_+(X)=\frac{X}{6\pi}+O(\sqrt X\log(2X)),\qquad
\sum_{\substack{\text{positive triangles}\\h\le X}}1
=\frac{X\log X}{6\pi}+O(X).                           \tag{15}
\]

This is a ratio among triples counted with multiplicity by a hypotenuse cutoff,
not a density assertion about distinct represented integers.
It does **not** imply `P_+(s)=P(s)/3`. The `C=5` coefficients are already zero
and one respectively. Even the augmented coefficients `δ_{C,1}+2r_+(C)` fail
multiplicativity: those at `5` and `13` vanish, but at `65` the coefficient is two.
The admitted `C=65` shape is `(63,16,65)`; the other primitive shape is excluded.

The exact scaling factors in (9) continue to hold after replacing `P` by `P_+`.
The loss occurs in the primitive angular filter, not in the scale decomposition.

### Theorem 2: an exact Gaussian angular correction

Use parity-oriented legs `(a,b)=(n²−m²,2nm)` and let
`θ=arg(a+ib)=2φ∈(0,π/2)`. Define, for `k≥0`,

\[
S_k(s)=\sum_{\text{primitive PPTs}}\cos(4k\theta)C^{-s},
\qquad S_0(s)=P(s).
\]

Size-ordering the legs replaces `θ` by `π/2−θ` when necessary and does not
change `cos(4kθ)`. Choose a Gaussian prime `π_p` of norm `p` for each
`p≡1 mod4`. Then, absolutely for `Re(s)>1`,

\[
\boxed{\displaystyle
S_k(s)=\frac12\left[
\prod_{p\equiv1\ (4)}
\frac{1-p^{-2s}}
{1-2\cos(8k\arg\pi_p)p^{-s}+p^{-2s}}-1\right].}        \tag{16}
\]

The choice of conjugate or associate `π_p` does not change this product.
The exact cone series is

\[
\boxed{\displaystyle
P_+(s)=\frac13P(s)+\frac2\pi\lim_{M\to\infty}
\sum_{k=1}^{M}\left(1-\frac{k}{M+1}\right)
\frac{\sin(k\pi/3)}{k}\,S_k(s),\qquad\Re s>1.}          \tag{17}
\]

The limit is a Fejér limit, locally uniform in this half-plane. No absolute
convergence assertion for the outer `k`-series is needed or intended.

**Proof of (16).** The Gaussian parameter has `θ=2 arg(n+im)`, so the weight
is the real part of `(α/|α|)^{8k}`. At `p^e || C`, the two allowed primitive
factor choices contribute `2cos(8ke arg π_p)`. They multiply independently over
the split primes. Summing over `e` gives the elementary identity

\[
1+2\sum_{e\ge1}\cos(eu)t^e
=\frac{1-t^2}{1-2\cos(u)t+t^2},\qquad |t|<1.
\]

The same unit/conjugation quotient and axis subtraction as in §3 produce
`1+2S_k`, proving (16). Absolute convergence follows from
`|cos(4kθ)|≤1` and convergence of `P(Re(s))`.

**Proof of (17).** Extend the indicator
`f(θ)=1_{|sin(2θ)|<1/2}` periodically with period `π/2`. Its mean is `1/3`,
and its cosine coefficient at frequency `4k` is

\[
\frac4\pi\left(\int_0^{\pi/12}+\int_{5\pi/12}^{\pi/2}\right)
\cos(4k\theta)\,d\theta
=\frac{2\sin(k\pi/3)}{\pi k}.
\]

The Fejér approximants have the coefficients in (17), lie between zero and
one, and converge to `f` at every non-boundary angle. No primitive triple lies
on the boundary: it would require `C²=4ab` with odd `C`.
Dominated convergence using `Σ C^−Re(s)<∞` now proves (17).
For compact subsets of `Re(s)>1`, one summable majorant works uniformly, giving
local uniformity. QED.

### Equivalent ideal-series normalization

Because `Z[i]` is principal, the multiplicative angular character on nonzero
Gaussian ideals is well defined by

\[
\chi_{8k}((\alpha))=(\alpha/|\alpha|)^{8k},\qquad
L_k(s)=\sum_{\mathfrak a\ne0}\chi_{8k}(\mathfrak a)N(\mathfrak a)^{-s}.
\]

Changing a generator by a Gaussian unit has no effect. This is the Gaussian
angular Hecke L-series, but only its explicit convergent ideal sum and Euler
factors are used here. A split prime gives the two denominator factors in
(16); an inert prime contributes `(1−p^−2s)^−1`. At the ramified prime,
`χ_{8k}((1+i))=exp(2πik)=1`, giving `(1−2^−s)^−1`. Thus

\[
1+2S_k(s)=\frac{L_k(s)}{(1+2^{-s})\zeta(2s)}.           \tag{18}
\]

At `k=0`, `L_0(s)=ζ(s)β(s)`, recovering (1). Equations (16)–(18) identify
the precise angular information missing from the one-third density heuristic.
No claim of a finite Hecke combination or elementary closed form for `P_+(5),
P_+(7),P_+(9)` follows.

## 8. Finite-exact controls and reproducibility

Run the standard-library script:

```bash
python3 04-computation/experiments/ap_square_odd_zeta_20261019.py
python3 -O 04-computation/experiments/ap_square_odd_zeta_20261019.py
```

Before repository integration, run the same file in its staging directory.
Both modes have identical output because the checks do not use removable
Python `assert` statements.

The frozen universes and independent paths are:

- `C≤100000`: **15,919** primitive triples, of which **5,305** are positive.
  Euclid enumeration agrees with (5) at every one of the 100,000 integers.
- Root height `h≤300`: **209** signed triangles, **33** positive.
  Direct enumeration of `p²+r²=8h²` agrees exactly with the square-scaled
  Pythagorean census; all inverse maps and pairwise gcds are checked.
- `n≤10000`: independent coefficient convolution of `ζ`, `β`,
  `1/ζ(2s)`, and `1/(1+2^−s)` agrees with `δ_{n,1}+2r(n)`.
- Gaussian norm `≤1000`: direct enumeration of all nonzero integer coordinate
  pairs checks the full Epstein coefficients against `4ζ(s)β(s)` and the
  primitive, half-counted coefficients against `2(1+2^−s)(1+2P(s))`.
- `C≤1000`, modes `k=0,1,2,3`: exact rational values
  `cos(4kθ)=Re((a+ib)^{4k})/C^{4k}` agree with the Gaussian-prime coefficient
  products from (16). This tests the doubled angular frequency and orientation
  quotient without trigonometric floating-point comparisons.
- Euler/Bernoulli recurrences give the rational constants in (2) exactly.
  Independent elementary zeta/beta series enclosures overlap both the Euler
  product and odd-value formulas. Machin's formula supplies a rational pi
  enclosure; floating-point special-function libraries are unnecessary.

There are **117,669 passing checks**, including the axis, split-prime powers,
inert/ramified exclusions, multiplicities, scale convolution, and the failed
positive-filter multiplicativity control.

For a rigorous tail bound, the number of all possible parameter pairs below
`C≤t` is at most `t/2`, hence `N(t)≤t/2`. For real `s>1`, partial summation
therefore gives

\[
0\le\sum_{C>X}r(C)C^{-s}
\le\frac{s}{2(s-1)}X^{1-s}.                            \tag{19}
\]

The same bound holds for the positive subset. At `X=100000`, these tails are
at most `6.25·10^−21`, `(7/12)·10^−30`, and `(9/16)·10^−40` for `s=5,7,9`.
The resulting positive/full ratios, with the bounds rounded outward, are

| s | Certified enclosing interval for P_+(s)/P(s) |
|---:|---|
| 5 | [3.429860722351164301·10^−5, 3.429860722353095892·10^−5] |
| 7 | [4.449350234251218669·10^−7, 4.449350234251218670·10^−7] |
| 9 | [6.236975535236050824·10^−9, 6.236975535236050825·10^−9] |

These are weighted-series ratios, so they need not resemble the large-cutoff
count ratio `1/3`. Higher `s` puts almost all mass near the first unrestricted
hypotenuse `5`, while positivity does not begin until `41`.

The companion audit independently rederives the normalizations, convergence
thresholds, counts, exact constants, and angular extension. Computation verifies
the finite universes; the unbounded statements rest on the proofs above.
