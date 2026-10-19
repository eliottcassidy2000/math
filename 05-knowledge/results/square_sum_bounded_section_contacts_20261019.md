# A bounded rational-section atlas gives polynomial-height plus-port contacts

**Status: PROVED**, with an independent proof-critical audit passed. The
companion computation is **FINITE-EXACT**. A new section near the nodal cusp,
combined with the already proved recurrence budget away from it, gives a
constant-size rational atlas for a distinct admissible plus-port neighbor.
Consequently the output hypotenuse is polynomial in the numerical input
hypotenuse, and its bit length is linear in the input bit length. No novelty
claim or optimized exponent is made.

## 1. Statement, inheritance, and scope

Let `T=(a,b,C)` be a primitive positive Pythagorean triple with `C²>4ab`.
Its normalized spread is `f_T=2ab/C² in (0,1/2)`.

**Theorem 1.** There are effective absolute constants `A,B>0` and a
deterministic exact algorithm returning a distinct admissible primitive
triple `U=(a',b',C')` and a positive rational number `r` such that

\[
 r^2=1+f_T-f_U,\qquad C'\le B C^A.                         \tag{1}
\]

The algorithm uses a bounded number of rational field operations and
comparisons, independently of the input, on `O(log C)`-bit operands.
Its output and contact certificate have `O(log C)` bits, and it has
polynomial bit complexity. Absorbing the fixed factor `B` into a larger
absolute exponent yields the equivalent numerical bound `C'<=C^{A'}`.

This theorem prescribes the anchor's **plus port** and the new petal's minus
port. The proof does not give a minus-port analogue, a prescribed spread
direction away from the cusp, arbitrary forbidden-set avoidance, or an
arbitrary requested spread window.

The inheritance chain is:

- [Universal signed Pythagorean contacts](square_sum_ppt_bridges_20261019.md),
  Lemmas 2–3 and Theorem 5: the quartic/cubic map, infinite-order point,
  primitive normalization, and exact port predicate.
- [Exact rational recurrence budget](square_sum_contact_recurrence_bound_20261019.md),
  equation (18): an even signed success with
  `0<|q|<=2 floor(100(4+J)/7)` when `2^J f_0>=2+f_0`.
- [Sharp cubic output length](square_sum_contact_sharp_output_height_20261019.md):
  the prescribed coset `Q+ZP` has cubic output-bit growth on the outer
  Berggren spines. That is a hostile example for this *coset*, not for
  arbitrary contact sections.

The new operation is `R=3Q+2P=Q+2(Q+P)`. It changes the coefficient of the
anchor section and supplies a nontrivial section tending to the anchor at
the cusp. The live concepts are signed ports, swapped-leg symmetry, nodal
specialization, rational sections, finite recurrence, and rational height.
The new small displacement does not come from a bounded nonzero multiple
of the old translating section `P`.

## 2. Put every plus-port anchor in one bounded parameter cone

Use the rational Pythagorean parameter and square root

\[
 f(t)=\frac{4t(t^2-1)}{(t^2+1)^2},\qquad
 k(t)=\frac{t^2+2t-1}{t^2+1},\qquad k(t)^2=1+f(t).       \tag{2}
\]

The involution

\[
 \iota(t)=\frac{t+1}{t-1}                                \tag{3}
\]

exchanges the normalized legs, satisfies `f(iota(t))=f(t)` and
`k(iota(t))=k(t)`, and maps each of the two positive admissible parameter
intervals to the other. Thus choose `t=min(t_0,iota(t_0))`, where for example
`t_0=(C+a)/b`. Then

\[
 1<t<\cot(5\pi/24)<4/3,\qquad 0<f(t)<1/2.               \tag{4}
\]

On this low interval `f` is strictly increasing: its derivative is
`-4(t⁴-6t²+1)/(t²+1)³`, whose first zero above 1 is `1+sqrt(2)`.
The branch selection in (3) is an exact rational comparison; no
trigonometric computation is required.

If the reduced selected parameter is `t=n/m`, primitive normalization gives

\[
 n^2+m^2=\delta C,\qquad \delta\in\{1,2\},\qquad
 H(t)=\max(n,m)\le\sqrt{2C}.                             \tag{5}
\]

This includes the odd/odd parameter arising when the legs are exchanged.

## 3. An explicit cusp section

Write `K=k(t)²` and use

\[
 E_t:Y^2=(X-K)(X^2-1),\quad P=(0,k),\quad
 Q=\left(\frac{t^2+4t+1}{t^2+1},
          \frac{8t(t+1)}{(t^2+1)^2}\right).              \tag{6}
\]

The anchor quartic point is `(u,z)=(t,t²+1)`. The group law gives

\[
 Q+P=\left(\frac{t^2-4t+1}{t^2+1},
           -\frac{8t^2(t-1)}{(t^2+1)^2}\right),          \tag{7}
\]

and in particular

\[
 X\bigl(2(Q+P)\bigr)
 =\frac{(t^2+1)(5t^2-4t+1)}{4t^2(t-1)^2}.               \tag{8}
\]

Define four integer polynomials

\[
\begin{aligned}
 A_4&=t^4-8t^3+2t^2+1,\\
 B_6&=5t^6-12t^5+23t^4-17t^2+12t-3,\\
 C_5&=2t^5-t^4-4t^3-2t^2+2t-1,\\
 D_6&=t^6+20t^5-29t^4+16t^3+3t^2-4t+1,
\end{aligned}
\]

and put `N=t A_4 B_6`, `E=C_5 D_6`. The inverse quartic parameter of
`R=3Q+2P` is the degree-11 rational function

\[
 \boxed{\quad u(t)=\frac{N(t)}{E(t)}.\quad}              \tag{9}
\]

For a completely explicit square certificate, let `Z(t)` be the degree-22
polynomial whose coefficients, from degree 22 down to degree 0, are

```
1, 824, -6469, 21392, -26297, 24216, -66331, 100544,
-7414, -132112, 150510, -48288, -39122, 44912, -11894,
-8000, 8293, -3368, 743, -112, 27, -8, 1.
```

Direct group-law simplification gives `z(t)=Z(t)/E(t)²`. Independently,
the integer-polynomial identity

\[
 (t^2+1)^2Z^2
 =(t^2+2t-1)^2(N^2+E^2)^2
       -4NE(N^2-E^2)(t^2+1)^2                           \tag{10}
\]

proves

\[
 \left(\frac{Z}{N^2+E^2}\right)^2=1+f(t)-f(u(t)).        \tag{11}
\]

Thus (9) and (11) can be evaluated without any intermediate elliptic
addition denominators. The checker verifies both the group-law derivation
and the independent polynomial identity.

### Lemma 2. The entire rational cusp chart is valid

For every real `1<t<=17/16`,

\[
 1<t<u(t)<5/4,\qquad 0<f(t)<f(u(t))<1/2.                \tag{12}
\]

All denominators in (9) and (11) are nonzero.

**Proof.** Put `h=t-1`, so `0<h<=1/16`. The shifted factors are

\[
\begin{aligned}
 A_4(1+h)&=h^4-4h^3-16h^2-16h-4<0,\\
 B_6(1+h)&=5h^6+18h^5+38h^4+72h^3+76h^2+40h+8>0,\\
 C_5(1+h)&=2h^5+9h^4+12h^3-8h-4<0,\\
 D_6(1+h)&=h^6+26h^5+86h^4+120h^3+92h^2+40h+8>0.
\end{aligned}                                           \tag{13}
\]

For `C_5`, the positive tail at `h=1/16` is
`1609/524288<1`. The first inequality is immediate from `h⁴<1`.
Consequently `E<0`, so there is no pole in the chart.

The coefficients of `-(N-tE)/h`, in ascending powers of `h`, are

```
64, 448, 1440, 2976, 4704, 6048, 6000, 4176, 1860, 476, 58, 2.
```

They are all positive. Thus `N-tE<0`, and division by `E<0` proves `u>t`.
For the other endpoint, expansion gives

\[
\begin{split}
4N-5E={}&32-160h-2000h^2-7392h^3-15632h^4-22512h^5\\
 &-23576h^6-17744h^7-8878h^8-2554h^9-293h^{10}+10h^{11}.
\end{split}                                             \tag{14}
\]

The sum of the absolute negative terms is bounded by its value at `1/16`,
namely `21856812814021/1099511627776<20`. Hence `4N-5E>12>0`.
Dividing by the negative `E` gives `u<5/4`. The function `f` is increasing
throughout `(1,5/4)`, and `f(5/4)=720/1681<1/2`, proving (12).
Finally `N²+E²>0`, while `t²+1>0`, proving the claimed denominator facts.
The radicand in (11) is positive since `f(t)>0` and `f(u)<1/2`, so its
rational root is nonzero. QED.

For orientation, the exact local expansions are

\[
 u(1+h)=1+3h+2h^3+O(h^4),
\]
\[
 f(1+h)=2h-h^2-h^3+O(h^4),\quad
 f(u(1+h))=6h-9h^2-23h^3+O(h^4).                         \tag{15}
\]

In particular the output/input spread ratio tends to 3. The proof of the
whole chart uses the polynomial inequalities, not an asymptotic remainder.
The cutoff matters: this same section fails admissibility at `t=5/4`.

## 4. A constant recurrence budget on the rest of the cone

Suppose `17/16<t<cot(5pi/24)`. Monotonicity and (2) give

\[
 f_0=f(t)>f(17/16)=\frac{35904}{297025}>\frac1{10}.
\]

With `D=2+f_0<5/2`, one has `D/f_0<25<32`. Therefore the least nonnegative
integer `J` with `2^J f_0>=D` is at most 5. The inherited exact rational
recurrence theorem supplies a new admissible quartic point among

\[
 Q+qP,\qquad q\in\{-256,-254,\ldots,-2,2,\ldots,254,256\}.
                                                               \tag{16}
\]

Indeed `2 floor(100(4+5)/7)=256`. Its guaranteed successful point has
`u>1`, `z>0`, `0<f(u)<1/2`, and `f(u)!=f_0`, in the ordinary inverse chart.
This covers anchors arbitrarily close to the upper admissibility wall:
there is no additional degeneration at `f_0=1/2`.

The finite atlas consists of the new section (9) on the cusp chart,
and the 256 signed even translates (16) on the rest. At any exceptional
addition or inverse chart, use the ordinary exact identity/inverse/doubling
branches or skip the unusable candidate. The inherited recurrence theorem
guarantees that one candidate lies in a nonexceptional successful chart.
No numerical real integration or undecidable endpoint test is used.

## 5. Primitive normalization and numerical height

Given an accepted `u=n'/m'>1` in lowest terms, let `delta'=2` if both
`n',m'` are odd, and let `delta'=1` otherwise. Return

\[
 U=\frac1{\delta'}(n'^2-m'^2,2n'm',n'^2+m'^2),             \tag{17}
\]

with its legs reordered to odd/even convention, and return the absolute
value of the rational root in (11) or its inverse-map analogue.
Formula (17) is a positive primitive triple and preserves the spread.
Different spreads imply different unordered primitive petals, so the
output is distinct from the anchor even after leg exchange.

For the cusp section, the coefficient absolute sums of `N,E` are
respectively 760 and 446. Homogenize both polynomials to degree 11 and
evaluate at coprime `(n,m)`. Reduction and (17) can only decrease the
hypotenuse, so (5) gives the explicit estimate

\[
 \boxed{\quad C'\le(760^2+446^2)H(t)^{22}
               \le776516(2C)^{11}.\quad}                \tag{18}
\]

For the rest of the atlas, each of the finitely many sections `Q+qP` has
rational-function coordinates in `t`, and so does its quartic inverse.
The finite degrees and integer coefficient norms are absolute constants.
For `g(t)=p(t)/q(t)` of degree at most `d`, homogeneous evaluation proves

\[
 H(g(t))\le \max(\|p\|_1,\|q\|_1)H(t)^d                 \tag{19}
\]

whenever the denominator is nonzero. Taking maxima over the finite atlas
and its finite exceptional group-law branches proves (1). Alternatively,
run the bounded exact arithmetic procedure directly: for rational height,
`H(xy),H(x/y)<=H(x)H(y)` and `H(x+y)<=2H(x)H(y)`. A bounded number of
operations on the `O(log C)`-bit initial coordinates therefore preserves
`O(log C)`-bit size. Equality and sign tests, reductions, and exceptional
branches have the same uniform bit bound. This argument also shows the
constants in the theorem are effective without symbolically expanding
all 256 translates.

Using schoolbook arithmetic, each bounded-size rational operation has
polynomial bit cost; even the conservative unamortized Euclidean analysis
gives `O((1+log C)^3)` total bit time and `O(1+log C)` working bits. These
bounds have large absolute constants and are not advertised as practical
or optimal. Output length `O(log C)` is a worst-case upper bound; no
matching per-input lower bound is claimed.

As an immediate graph consequence, choose `H=lcm(C,C')<=CC'` and the common
center `M=2H²`. The two AP petals and their mixed-sign bridge form the
five-vertex square-sum cycle supplied by the contact theorem, with every
vertex below `4H²<=4B² C^{2(A+1)}`. This is a polynomial numerical bound
for one glued AP pair. It says nothing about coverage of interval labels
or a polynomial bound for arbitrarily long distinct chains.

## 6. Why this escapes the old coset: two independent generic sections

The construction already proves Theorem 1 without this interpretation.
It is nevertheless useful to verify that its new operation is genuinely
outside the old generic translate coset.

**Proposition 3.** The sections `P,Q` are linearly independent over `Z` in
`E(Q(t))`.

**Proof.** At `t=1` the cubic is `Y²=(X-1)²(X+1)`. Its only singular point
is `(1,0)`, whereas the reductions

\[
 P_0=(0,1),\quad Q_0=(3,4),\quad (Q+P)_0=T=(-1,0)
\]

are smooth, and `T` has order 2 in the smooth nodal group. Over
`Q(sqrt(2))`, put `w=Y/(X-1)` and

\[
 \zeta=\frac{w-\sqrt2}{w+\sqrt2}.
\]

A direct chord calculation gives
`w(A+B)=(w(A)w(B)+2)/(w(A)+w(B))`; substitution therefore proves that
this identifies the smooth nodal group, with the point at infinity as
identity, with the multiplicative group. In particular
`zeta(P_0)=-(1+sqrt(2))²`, of infinite order, and `zeta(T)=-1`.
Reduction is a group homomorphism on sections meeting this smooth locus.
Thus a generic relation `aP+bQ=O` reduces to
`(a-b)P_0+bT=O`, forcing `a=b` by the absolute value of `zeta(P_0)`.

Let `S=2(P+Q)`. Equation (8) has a genuine double pole at `t=1`, and the
y-coordinate is

\[
 Y(S)=-\frac{(t^2+2t-1)(3t^2-2t+1)(t^4-6t^3+2t-1)}
                 {8t^3(t-1)^3(t^2+1)}.
\]

Consequently `S` is nonzero and specializes to `O`; its formal parameter
`v=-X(S)/Y(S)` satisfies `v=-(t-1)+O((t-1)²)`.
The formal group at the smooth identity over `Q[[t-1]]` has
`[n]v=nv+O(v²)` for every nonzero integer `n`. Since the residue field has
characteristic zero, the leading nonzero coefficient cannot vanish; no
nonzero point in this reduction kernel is torsion. Hence `S` is not
torsion. A nonzero relation with `a=b` would make `P+Q`, and therefore
`S`, torsion, a contradiction. QED.

The nodal multiplicative coordinate and the linear term of the formal
group explain the mechanism: `P` retains a nontrivial limiting
translation, whereas `2(P+Q)` tends to the identity and is still a nonzero
rational section. Since `P,Q` are independent, `3Q+2P` is not generically
`Q+qP` for any fixed integer `q`. No claim of independence at every
rational specialization follows, or is needed. For the standard smooth nodal
group, reduction homomorphism, and formal-law facts used here, see
[J. S. Milne, Elliptic Curves](https://www.jmilne.org/math/Books/ectext6.pdf),
Chapter II §§3–4 and Remark 4.3; the characteristic-zero kernel argument is
given explicitly above.

## 7. Exact controls and remaining boundary

The companion script checks the rational map and its square certificate
by two routes, the sign-polynomial proofs on the closed rational cusp
chart, the local Taylor coefficients, leg-swap invariance, primitive
normalization, and the fallback budget. It enumerates all 528 admissible
primitive odd/even oriented triples with `C<=10000`:

- 119 use the new cusp section;
- 409 use a fallback even translate, with successful signed indices
  `-6,-2,4,6,10` in this finite census;
- both normalization divisors occur: 269 outputs use divisor 1, and 259
  use divisor 2.

Selected exact probes cover the cutoff, both swapped-leg cusps, and
parameters with 30-digit components. The global algorithm is checked to
return the same petal and certificate after swapping the input legs.
The nonlocal failure at `t=5/4` is retained as a hostile control.

Run, with SymPy installed:

```bash
python3 04-computation/experiments/square_sum_bounded_section_contacts_20261019.py
python3 -O 04-computation/experiments/square_sum_bounded_section_contacts_20261019.py
```

The normal and optimized executions have byte-identical output. The controls
make no census-based extrapolation. Their 16,628 exact rational and 23 exact
symbolic checks are separate from the uniform proof.
The connection contract is: source, all admissible primitive PPTs with
selected plus port; target, one distinct admissible primitive minus-port
neighbor; map, leg-swap normalization followed by the finite section atlas;
preserved predicate, the rational mixed-sign square contact; indispensable
sidecars, the cusp cutoff, strict changed spread, exceptional charts, and
primitive divisor. The atlas gives one neighbor, not enough choices to
avoid an arbitrary finite forbidden set. That remains a separate problem,
as does an analogous bounded atlas for a prescribed minus port.
