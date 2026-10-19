# A bounded rational-section atlas for minus-port contacts

**Status: PROVED**, with an independent proof-critical audit passed. The
companion computation is **FINITE-EXACT**. This addendum closes the missing
port in the bounded-section construction. Its cusp chart uses a different
sheet of the anchor quartic. Together with the plus-port note, it gives a
polynomial numerical-height and linear output-bit bound for either
prescribed mixed-sign anchor port. No novelty or optimal-exponent claim is
made.

## 1. Statement and inheritance

Let `T=(a,b,C)` be a primitive positive Pythagorean triple with `C²>4ab`, and
put `f_T=2ab/C² in (0,1/2)`.

**Theorem 1 (minus port).** There are effective absolute constants `A,B>0`
and a deterministic exact algorithm producing a distinct admissible
primitive triple `U=(a',b',C')` and a positive rational `r` such that

\[
 r^2=1-f_T+f_U,\qquad C'\le B C^A.                      \tag{1}
\]

The algorithm uses a bounded number of rational operations and comparisons
on `O(log C)`-bit operands. Its triple and square certificate have
`O(log C)` bits, and its bit complexity is polynomial in `log C`.

**Corollary 2 (both prescribed ports).** For either prescribed
`epsilon in {-1,+1}`, the same conclusions hold with

\[
 r^2=1+\varepsilon f_T-\varepsilon f_U.                 \tag{2}
\]

Take the maxima of the absolute constants here and in
[the plus-port bounded atlas](square_sum_bounded_section_contacts_20261019.md).
The negative input port is joined to the positive output port in (1);
the opposite signs occur for `epsilon=+1`.

The inherited results are:

- [Universal signed Pythagorean contacts](square_sum_ppt_bridges_20261019.md),
  Lemmas 2–3 and Theorem 5: the quartic/cubic maps, non-torsion translating
  point, primitive normalization, and signed-port criterion.
- [Exact rational recurrence budget](square_sum_contact_recurrence_bound_20261019.md),
  equation (18): a successful even signed translate satisfies
  `0<|q|<=2 floor(100(4+J)/7)` if `2^J f_0>=D`; for the minus port `D=2`.
- [The plus-port bounded atlas](square_sum_bounded_section_contacts_20261019.md):
  low-cone leg-swap normalization and the first prescribed-port height
  theorem.

Only one distinct neighbor is promised. Neither a prescribed spread
direction away from the cusp, an arbitrary spread window, nor avoidance
of an arbitrary finite forbidden set is supplied by this atlas.

## 2. Normalize to the low cone and select the negative sheet

Use

\[
 f(t)=\frac{4t(t^2-1)}{(t^2+1)^2},\qquad
 \iota(t)=\frac{t+1}{t-1}.
                                                               \tag{3}
\]

The involution exchanges normalized legs and preserves `f`. Starting with
`t_0=(C+a)/b`, select `t=min(t_0,iota(t_0))`. Then

\[
 1<t<\cot(5\pi/24)<4/3,\qquad 0<f(t)<1/2.              \tag{4}
\]

On this low cone, the positive rational square root is

\[
 k=\frac{-t^2+2t+1}{t^2+1}>0,\qquad K=k^2=1-f(t).      \tag{5}
\]

For reduced `t=n/m`, primitive normalization gives

\[
 n^2+m^2=\delta C,\quad \delta\in\{1,2\},\qquad
 H(t)=\max(n,m)\le\sqrt{2C}.                           \tag{6}
\]

The quartic and cubic are

\[
 z^2=K(u^2+1)^2-4u(u^2-1),\qquad
 E_t:Y^2=(X-K)(X^2-1).                                 \tag{7}
\]

Their ordinary coordinate maps are

\[
 X=\frac{kz+K+2u+Ku^2}{2u^2},\qquad
 Y=\frac{(X-K)(Xu-1)}k,                                \tag{8}
\]

\[
 u=\frac{kY+X-K}{X(X-K)},\qquad
 z=\frac{2u^2X-K-2u-Ku^2}{k}.                           \tag{9}
\]

For the cusp chart choose the **negative-sheet** anchor
`(u,z)=(-t,-(t²+1))`, which maps to

\[
 Q_- =\left(\frac{t^2-4t+1}{t^2+1},
              -\frac{8t(t-1)}{(t^2+1)^2}\right).        \tag{10}
\]

This sheet choice is essential to the construction. At `t=1`, the cubic
is nodal, but `Q_-` tends to the smooth 2-torsion point `(-1,0)`. Moreover,

\[
 X(2Q_-)=\frac{(t^2+1)(t^2-4t+5)}{4(t-1)^2}.            \tag{11}
\]

Thus `2Q_-` tends to the identity and `3Q_-=Q_-+2Q_-` gives a small
displacement of the anchor. Formula (11) describes the mechanism; the
uniform admissibility proof below uses exact polynomial inequalities.

## 3. An explicit minus-port cusp section

Define integer polynomials

\[
\begin{aligned}
 A_4&=t^4+2t^2-8t+1,\\
 B_6&=3t^6-12t^5+17t^4-23t^2+12t-5,\\
 C_5&=t^5-2t^4+2t^3+4t^2+t-2,\\
 D_6&=t^6-4t^5+3t^4+16t^3-29t^2+20t+1,
\end{aligned}
\]

and set `N=C_5 D_6`, `E=A_4 B_6`. The inverse quartic coordinates of
`R=3Q_-` are

\[
 \boxed{\quad u(t)=-\frac{N(t)}{E(t)},\qquad
                 z(t)=-\frac{W(t)}{E(t)^2}.\quad}       \tag{12}
\]

Here `W` has degree 22; its coefficients from degree 22 down to degree 0
are

```
1, -8, 27, -112, 743, -3368, 8293, -8000, -11894, 44912,
-39122, -48288, 150510, -132112, -7414, 100544, -66331,
24216, -26297, 21392, -6469, 824, 1.
```

The group law gives (12) directly. A separate integer-polynomial identity
certifies the required square independently of any group-law chart:

\[
 (t^2+1)^2W^2
  =(t^2-2t-1)^2(N^2+E^2)^2
     +4NE(N^2-E^2)(t^2+1)^2.                           \tag{13}
\]

Writing `v=-u=N/E`, division gives

\[
 \boxed{\quad
 \left(\frac{W}{N^2+E^2}\right)^2=1-f(t)+f(v(t)).
 \quad}                                                \tag{14}
\]

The checker derives (12) by rational-function doubling and addition and
also checks (13) as a polynomial identity.

### Lemma 3. The entire closed rational cusp chart is valid

For every real `1<t<=17/16`,

\[
 1<t<v(t)<5/4,\qquad 0<f(t)<f(v(t))<1/2.                \tag{15}
\]

All denominators in (12) and (14) are nonzero.

**Proof.** Put `h=t-1`, so `0<h<=1/16`. The denominator expands as

\[
\begin{split}
 E(1+h)={}&32+32h-80h^2-128h^3-16h^4+48h^5+40h^6\\
          &+64h^7+50h^8+18h^9+3h^{10}.
\end{split}                                             \tag{16}
\]

The sum of its absolute negative terms is at most

\[
 80/16^2+128/16^3+16/16^4=1409/4096<1.
\]

Hence `E>31`, excluding poles. Next, the coefficients of
`L(h)=(N(1+h)-(1+h)E(1+h))/h`, in ascending degree, are

```
64, 256, 480, 384, 96, 0, -48, -96, -60, -16, -2.
```

The absolute negative tail at `h=1/16` is
`1777281/549755813888<1`. Thus `L>63>0`, so `N-tE>0` and `v>t`.
Finally,

\[
\begin{split}
5E-4N={}&32-352h-1232h^2-1728h^3-1040h^4-272h^5\\
        &-152h^6+96h^7+178h^8+58h^9-5h^{10}-4h^{11}.
\end{split}                                             \tag{17}
\]

The absolute negative tail at `h=1/16` is
`119849021866005/4398046511104<28`. Consequently `5E-4N>4>0`, and,
since `E>0`, `v<5/4`.

The derivative
`f'(s)=-4(s⁴-6s²+1)/(s²+1)³` is positive on `(1,5/4)`.
Also `f(5/4)=720/1681<1/2`. These facts prove the spread inequalities.
The denominator in (14) is positive because `E>0`, and its radicand is
strictly greater than 1. Therefore the rational square root is nonzero.
QED.

The exact local expansions are

\[
 v(1+h)=1+3h+6h^2+14h^3+O(h^4),
\]

\[
 f(1+h)=2h-h^2-h^3+O(h^4),\qquad
 f(v(1+h))=6h+3h^2-35h^3+O(h^4).                        \tag{18}
\]

Thus the output/input spread ratio tends to 3. The cutoff cannot simply
be discarded: at `t=5/4`, this section has
`v=345349397/138407340` and `f(v)>1/2`, so it is inadmissible.

### Reciprocal involution links the two port sections

Let `N_+, E_+, Z_+` be the numerator, denominator, and square-certificate
polynomials of the plus-port cusp section, with `u_+=N_+/E_+`.
Polynomial reversal gives the exact identities

\[
 t^{11}N_+(1/t)=-E(t),\qquad
 t^{11}E_+(1/t)=-N(t),\qquad
 t^{22}Z_+(1/t)=W(t).
\]

Therefore

\[
 \boxed{\quad v(t)=\frac1{u_+(1/t)}.\quad}              \tag{19}
\]

Since `f(1/s)=-f(s)`, substituting `1/t` into the plus section's rational
square identity transports it exactly into (14), including the root
`W/(N²+E²)`. This is an algebraic port duality of rational sections.
It does not transport the earlier admissibility theorem: when `t>1`,
its reciprocal is below 1, outside that theorem's cusp domain. The new
uniform sign inequalities in Lemma 3 remain necessary. The checker
independently verifies all three reversals, (19), and the spread identity.

The group-law mechanism also agrees. Because
`k_+(1/t)=k_-(t)`, the fixed-curve quartic involution
`(u,z) ↦ (-1/u,-z/u²)` maps the plus anchor at `1/t` to the minus
negative-sheet anchor at `t`. Under (8)–(9), its action on the cubic is
`R ↦ -R-P`, where `P=(0,k)`. For example its new X-coordinate is
`K+u²(K-X)`; the cubic equation makes this equal to the chord formula for
`X(-R-P)`, and the Y-coordinate agrees as well. The checker verifies both
coordinate identities modulo the cubic equation, as well as preservation
of the quartic. Thus it sends `3Q+2P` to `-3Q-3P=3Q_-`, explaining the
pure tripling in (12) and the translation term in the plus section.

## 4. A constant fallback on the rest of the low cone

For `t>17/16`, monotonicity in the low admissible cone gives

\[
 f_0=f(t)>\frac{35904}{297025}>\frac1{10}.                \tag{20}
\]

Use `P=(0,k)` and now change back to the **positive-sheet** anchor
`Q_+`, the image under (8) of `(-t,t²+1)`. This is exactly the anchor
required by the inherited recurrence theorem. The subscripts here denote
quartic sheets; the anchor port remains minus throughout.

For the minus-port recurrence `D=2`; hence `D/f_0<20<32`. The least
nonnegative integer `J` with `2^J f_0>=D` is at most 5. Consequently
one of the 256 candidates

\[
 Q_++qP,\qquad
 q\in\{-256,-254,\ldots,-2,2,\ldots,254,256\}           \tag{21}
\]

has an ordinary inverse-chart image satisfying

\[
 u<-1,\quad z>0,\quad 0<f(-u)<1/2,\quad f(-u)\ne f_0.
                                                               \tag{22}
\]

Indeed the recurrence bound is
`2 floor(100(4+5)/7)=256`. Formula (7) then gives
`r=|z|/(u²+1)` and `r²=1-f_0+f(-u)`.

The algorithm tests these predicates exactly. It may use the smaller
input-dependent budget, but never needs more than (21). Ordinary
identity, inverse and doubling cases handle elliptic additions; skip
candidates outside the inverse chart. The inherited recurrence guarantees
at least one successful nonexceptional candidate. This includes inputs
arbitrarily close to the upper admissibility boundary `f_0=1/2`.

## 5. Primitive normalization, height, and exact complexity

For either chart take `v=-u=n'/m'>1` in lowest terms. Let `delta'=2` when
both numerator and denominator are odd, and `delta'=1` otherwise. Return

\[
 U=\frac1{\delta'}(n'^2-m'^2,\ 2n'm',\ n'^2+m'^2),       \tag{23}
\]

with legs reordered to odd/even convention, and the positive absolute
value of the corresponding rational root. Formula (23) is primitive,
positive, and preserves spread. Changed spread implies a different
unordered primitive petal, so the output remains distinct after leg
exchange.

The cusp polynomials have coefficient absolute sums
`||N||_1=446`, `||E||_1=760`; their degrees are 11 and 10. Homogenize both
to degree 11. Reduction and division by `delta'` only decrease the
hypotenuse, and (6) gives the explicit bound

\[
 \boxed{\quad C'\le(446^2+760^2)H(t)^{22}
                  \le776516(2C)^{11}.\quad}            \tag{24}
\]

Away from the cusp, the fixed set (21) is a finite rational-section atlas,
with finitely many ordinary exceptional group-law cases. Each case
requires only a bounded number of rational operations on the rational
input `t`. The starting coordinates have `O(log C)` bits by (6). The
height inequalities

\[
 H(xy),H(x/y)\le H(x)H(y),\qquad
 H(x+y)\le2H(x)H(y)
\]

show that any bounded number of field operations preserves `O(log C)`
bit size. This gives effective absolute polynomial numerical-height
constants for both the output and the certificate, without expanding all
256 sections. Equality tests, reductions, and selection of ordinary
exceptional cases obey the same uniform bound. This proves Theorem 1.

Taking constants large enough for both the present atlas and the plus
atlas proves Corollary 2. The fixed multiplier can be absorbed into a
larger exponent because admissible input hypotenuses are bounded away
from 1. With schoolbook integer arithmetic and a conservative Euclidean
analysis, the bounded procedure uses `O((1+log C)^3)` bit operations and
`O(1+log C)` working bits. Constants are large and no optimized complexity
claim is intended.

For either port, the graph transfer at `H=lcm(C,C')<=CC'` and common
center `M=2H²` therefore gives a glued two-petal square-sum cycle with
all five positive integer labels below
`4H²<=4B² C^{2(A+1)}`. This is a polynomial bound for one pair, not for
arbitrarily long distinct chains or interval-graph coverage.

## 6. Exact controls and scope boundary

The companion checker validates the rational-function group-law derivation,
the independent integer-polynomial square identity, all three uniform
sign-polynomial bounds, Taylor coefficients, the rational cutoff, the
fallback budget, and primitive normalization. It enumerates all 528
admissible primitive odd/even oriented triples with `C<=10000`, testing
swapped-leg invariance for every input. Separate probes include the
cutoff, both swapped-leg cusps, and parameters with 30-digit components.
The nonlocal failure at `t=5/4` is a retained hostile control.

Run, with SymPy installed:

```bash
python3 04-computation/experiments/square_sum_minus_bounded_section_contacts_20261019.py
python3 -O 04-computation/experiments/square_sum_minus_bounded_section_contacts_20261019.py
```

The normal and optimized executions have byte-identical output. They pass
39,092 exact rational checks and 28 exact symbolic checks. In the complete
528-input census, 119 anchors use the cusp chart and 409 use the recurrence
chart. The largest successful absolute recurrence index is 40; both
primitive divisors occur, with 269 divisor-1 and 259 divisor-2 outputs.
The largest output hypotenuse has 43,466 bits, illustrating that the
uniform asymptotic theorem is not a small-constant practical-size claim.

Finite controls support implementation correctness; the uniform assertion comes
from Lemma 3, the inherited recurrence theorem, and the bounded-height
argument. The connection contract is: source, an admissible PPT with
prescribed minus port; target, one distinct admissible primitive
positive-port neighbor; map, low-cone normalization, negative-sheet
tripling near the cusp, and positive-sheet recurrence elsewhere;
preserved predicate, the rational mixed-sign square; necessary sidecars,
the sheet distinction, cutoff, strict changed spread, ordinary inverse
chart, and primitive divisor. A prescribed target window and arbitrary
finite forbidden-set avoidance remain outside this bounded-atlas theorem.
