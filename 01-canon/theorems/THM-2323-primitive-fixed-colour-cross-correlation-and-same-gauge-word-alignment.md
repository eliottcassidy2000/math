---
id: THM-2323
title: "Primitive fixed-colour cross-correlation and same-gauge word alignment"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. Let N be
  any positive integer divisible by seven. If 0<=f<=g are
  nonzero rational-valued step functions with rational breakpoints,
  supported in one open circle arc of length 1/7, then for every primitive
  N-character k every block of J_f J_g consecutive gauge indices contains
  one at which both Fourier coefficients f-hat(k+Nh) and
  g-hat(k+Nh) are nonzero. In particular there is a common nonnegative
  index 0<=h<=J_f J_g-1. The mechanism is a nonnegative rational
  fixed-colour cross-correlation Laurent polynomial: Galois straightens
  any primitive colour to the standard root, where its complete support
  lies in the acute sector |arg|<2*pi/7 and cannot cancel. Parseval and
  endpoint-product Vandermonde then land the common index. At N=91,
  applying the theorem to THM-2319's word source H_Q and bare source H_E
  gives, for every one of the 72 unit colours, a common index
  h<=12S^2-1 and common multiplier n<=1092S^2-1. This closes the
  bare/word same-gauge-index loss. Conditionally, for the middle owner,
  write c_2/gcd(c_2,c_3)=a and
  c_3/gcd(c_2,c_3)=d'. If gcd(a,91)=1, exposure at N=91d' and its
  K_7 x K_13 / K_6 x K_13 primitive fibre force two common bare/word
  atoms to form a unit-coloured c_3 edge, with no restriction on the
  unit cofactor of d'. If gcd(a,91)>1 this normalized sublattice cannot
  yield a unit-coloured edge. The theorem does not select target-plane
  gain, control terminal-component phase, or exclude a scalar row.
  LRC(14) remains open.
source: codex-2026-07-25-primitive-fixed-colour-cross-correlation
depends_on:
  - THM-2319-crt-unit-bispectrum-needle-and-mixed-polarization-no-go
related:
  - THM-2293-quadratic-root-energy-raises-the-ancestry-shell
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2321-prescribed-root-character-bispectrum-slice-positivity
script: 04-computation/lrc14_primitive_cross_correlation_same_gauge_thm2323.py
output: 05-knowledge/results/lrc14_primitive_cross_correlation_same_gauge_thm2323.out
script_sha256: 520e35ddd167e1087ee35cb9e34f3805fa6211e569c7a933a8123d594eb6dea2
output_sha256: c11dc72e33ade2f598b4f0551fc71e71ddde4709bcfb18cb4131a6bf32b16b0f
hash_basis: working-tree bytes (LF)
---

# THM-2323 -- primitive colour forces a bare/word diagonal

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2319 lands every primitive `91`-colour separately in the Fourier
spectra of the bare owner source and of its literal positive-word subset.
Separate landing does not formally give one gauge index at which the two
coefficients coexist. The complete primitive quadratic cross-energy cannot
repair this by positivity: THM-2319 gives an exact consecutive-needle
example for which its sum over all `72` colours is zero.

The missing operation is to freeze one primitive colour before taking the
cross-correlation. The physical support then makes the correlation a short
nonnegative rational Laurent polynomial. Galois conjugation straightens
the chosen primitive colour to the standard root; all physical shifts then
lie in one acute sector and cannot cancel. This is a phase-sensitive
algebraic nonvanishing statement, not a lower bound on the real part of the
original colour.

## 1. Root fibres and the fixed-colour cross-correlation

Let `N` be a positive multiple of seven and put

```text
R=N/7-1,
2*pi*R/N<2*pi/7<pi/2.                               (1)
```

Let `I` be one open circle arc of length `1/7`. Let `f,g` be
rational-valued step functions on the circle, with rational breakpoints,
such that

```text
0<=f<=g,
support(f),support(g) subset I,
integral f>0.                                        (2)
```

Boundary values are immaterial. In particular, (2) means nonzero on a set
of positive measure.

Write

```text
zeta=exp(2*pi*i/N),
x_r(y)=(y+r)/N mod 1,       r in Z/NZ,

M_(k,f)(y)=sum_r f(x_r(y))zeta^(-kr),                (3)
M_(k,g)(y)=sum_r g(x_r(y))zeta^(-kr).
```

Fix `k` with `gcd(k,N)=1`, and form the cross-correlation

```text
C_k=integral_0^1 M_(k,f)(y)
                    conjugate(M_(k,g)(y))dy.         (4)
```

For a residue `d modulo N`, put

```text
c_d=sum_r integral_0^1
       f(x_r(y))g(x_(r+d)(y))dy.                     (5)
```

Every `c_d` is a nonnegative rational number. Indeed, the integrand in
(5) is a rational-valued step function and all of its breakpoints are
rational.

If the integrand in (5) is nonzero, the two root samples lie in the same
open arc `I` and differ by `d/N` on the circle. Because the arc has length
exactly `1/7`, the least signed representative of `d` must satisfy

```text
|d|<N/7,
```

apart from boundary coincidences of measure zero. Since `N/7` is an
integer,

```text
c_d=0 unless -R<=d<=R.                              (6)
```

Expansion of (4) gives the short Laurent polynomial

```text
C_k=sum_(d=-R)^R c_d zeta^(kd).                     (7)
```

Its diagonal coefficient is strictly positive:

```text
c_0
 =sum_r integral f(x_r(y))g(x_r(y))dy
 =N integral f(x)g(x)dx
 >=N integral f(x)^2dx
 >0.                                                (8)
```

Thus

```text
Q(X)=sum_(d=-R)^R c_d X^d in Q[X,X^(-1)]           (9)
```

is a nonzero rational Laurent polynomial with nonnegative coefficients.

## 2. Galois straightens every primitive colour into an acute sector

Suppose `C_k=0`. Since `k` is a unit modulo `N`, the Galois automorphism

```text
sigma_(k^(-1)): zeta -> zeta^(k^(-1))
```

fixes every rational coefficient and sends `zeta^k` to `zeta`. Applying
it to (7) gives

```text
0=sigma_(k^(-1))(C_k)
 =sum_(d=-R)^R c_d zeta^d.                          (10)
```

But (1), nonnegativity of the `c_d`, and (8) give

```text
Re sum_(d=-R)^R c_d zeta^d
 =sum_(d=-R)^R c_d cos(2*pi*d/N)
 >=cos(2*pi*R/N) sum_(d=-R)^R c_d
 >cos(2*pi/7) sum_(d=-R)^R c_d
 >0.                                                (10a)
```

This contradicts (10). We have proved:

> **Fixed-colour cross-correlation lemma.** Under (1)--(2),
>
> ```text
> C_k!=0 for every k in (Z/NZ)^*.                    (11)
> ```

The rationality in (2) is load-bearing because it fixes every coefficient
under the Galois straightening. Nonnegativity is load-bearing in (10a),
and (8) makes the final inequality strict. The common short support is
load-bearing at (6). No totient lower bound is needed: arbitrary extra
prime factors of `N` do not weaken (11).

The lemma explains why THM-2319's aggregate hostile control is not a
counterexample. On sites `0` and `7`, take

```text
f=(1,0),
g=(1,12).
```

Then the complete sum over primitive colours is

```text
sum_(k unit mod 91) (1+12 zeta^(7k))
 =72+12 c_91(7)
 =72-72
 =0,                                                (12)
```

but every fixed summand is nonzero. Galois straightens it to
`1+12zeta^7`, whose two terms lie in the acute sector of (10a). The colour
sum can cancel even though no colour does.

## 3. Parseval gives a common gauge

Define the periodic gauges

```text
N_(k,f)(y)=exp(-2*pi*i*k*y/N)M_(k,f)(y),
N_(k,g)(y)=exp(-2*pi*i*k*y/N)M_(k,g)(y).            (13)
```

The root shift `y -> y+1` multiplies each `M_k` by `zeta^k`, so the
factor in (13) makes the gauges genuinely one-periodic. Direct branchwise
change of variables gives, for every integer `h`,

```text
(N_(k,f))_hat(h)=N f_hat(k+Nh),
(N_(k,g))_hat(h)=N g_hat(k+Nh).                     (14)
```

The gauge factors cancel in the physical inner product. Parseval and
(11) therefore give

```text
0!=C_k
 =sum_(h in Z)
    (N_(k,f))_hat(h)
    conjugate((N_(k,g))_hat(h)).                    (15)
```

The series in (15) is absolutely convergent by Cauchy--Schwarz. Hence at
least one integer `h` obeys

```text
f_hat(k+Nh)g_hat(k+Nh)!=0.                          (16)
```

At this point the index can have either sign. The endpoint recurrence
makes it nonnegative with a uniform finite bound.

## 4. Endpoint-product Vandermonde lands a nonnegative index

Let `J_f,J_g` be the numbers of nonzero jumps of `f,g`. The root amplitudes
in (3) have at most `J_f,J_g` distinct image-jump nodes respectively.
With

```text
alpha=k/N,
```

distributional differentiation of the covariant step gauges gives

```text
2*pi*i(h+alpha)(N_(k,f))_hat(h)
 =sum_x Delta_(k,f)(x)
      exp(-2*pi*i(h+alpha)x),                       (17)
```

and the analogous identity for `g`. The endpoint at the chosen circle
section is included with its covariant jump, so (17) has no omitted
boundary term.

Multiply (17) by the conjugate of its `g` analogue and combine equal
endpoint-difference nodes. The resulting sequence in `h` is an
exponential sum on at most

```text
L<=J_f J_g                                           (18)
```

distinct nonzero nodes. It is not the zero sequence: (16) supplies one
nonzero term, and `h+alpha` never vanishes for an integer `h` because
`0<k<N`.

A nonzero exponential sum on `L` distinct nodes cannot vanish at `L`
consecutive integers, by invertibility of the corresponding Vandermonde
matrix. The same argument applies to a block beginning at any integer.
Equations (14), (17), and (18) yield

> **Primitive same-gauge theorem.** For every `k` coprime to `N` and every
> integer `H`, some
>
> ```text
> H<=h<=H+J_f J_g-1
> ```
>
> satisfies
>
> ```text
> f_hat(k+Nh)g_hat(k+Nh)!=0.                         (19)
> ```

Taking `H=0`, the first common positive frequency obeys

```text
1<=k+Nh<=N J_f J_g-1.                               (20)
```

This conclusion is diagonal: it retains the same primitive colour and the
same gauge index in both functions. A separate landing theorem for the two
functions would not imply it. It also makes the common-index set
two-sided syndetic, with gauge gaps at most `J_fJ_g`; the first bounded
positive atom is only its smallest immediately useful consequence.

## 5. The arbitrary-prime stalk

There is no restriction on the primes or exponents in `N/7`. In
particular, take

```text
N_a=7*13^a,                 a>=1.
```

Then

```text
R_a=13^a-1,
2*pi*R_a/N_a<2*pi/7<pi/2.                          (21)
```

For every primitive `N_a`-colour, the same
pair `f<=g` has a common gauge satisfying (19), with frequency at most
`N_a J_fJ_g-1`.

This is a genuine stalk law: increasing the thirteen depth enlarges the
physical root needle from `13` to `13^a` teeth, while Galois always
straightens the selected primitive colour back into the same acute sector.
At each fixed depth and colour the common spectrum is a bounded-gap
toothpick. The theorem is levelwise and does not assert compatible gauge
indices between depths.

More generally, one may multiply `N_a` by any integer. The proof is
unchanged as long as the resulting modulus remains divisible by seven.
This arbitrary-prime freedom is load-bearing in Section 7.

## 6. Exact LRC word/bare specialization at `N=91`

Use THM-2319's notation for a positive canonical word `Q` at owner `j`:

```text
E_Q=E_j intersection T^(-(lambda_j+1))Q,

H_Q=P_(c_j)1_(E_Q),
H_E=P_(c_j)1_(E_j).                                 (22)
```

The literal inclusion `E_Q subset E_j` and positivity of the Perron
operator give

```text
0<=H_Q<=H_E.                                        (23)
```

Both functions are rational-valued step functions with rational
breakpoints, both are supported in `D_1`, and

```text
integral H_Q=measure(E_Q)>0.                         (24)
```

The jump ledger already proved in THM-2319 is

```text
J_(H_E)<=2S,
J_(H_Q)<=6S.                                        (25)
```

For `N=91`,

```text
R=12,
2*pi*R/91<2*pi/7<pi/2.
```

Apply (19) to (22). For every one of the `72` unit colours `k`, every
block of `12S^2` consecutive gauge indices contains a common index. In
particular there is one

```text
0<=h<=12S^2-1                                       (26)
```

such that, with

```text
n=k+91h<=1092S^2-1,
```

one has

```text
(H_E)_hat(n)(H_Q)_hat(n)!=0,

(1_(E_j))_hat(c_j n)
(1_(E_Q))_hat(c_j n)!=0.                            (27)
```

The second line uses the exact Perron identities from THM-2319. Hence the
same ordinary Fourier atom is simultaneously a bare-source atom and an
atom of the literal subset whose prescribed orbit realizes the complete
terminal word `Q`.

Given a nonzero shell character `kappa modulo 13`, choose any unit
`k modulo 91` satisfying

```text
u_j k congruent kappa mod 13.                        (28)
```

There are six such choices, and (26)--(27) hold for every one. Since `n`
is coprime to `91`, the common atom has exact owner grade
`nu_13(c_j)=lambda_j`, prescribed normalized residue `kappa`, and a whole
outside multiplier prime to seven and thirteen. For the middle owner
`j=2`, it is a grade-`b` marked bare vertex with the same positive-word
coefficient at the same gauge index.

This removes the same-gauge-index loss explicitly left in THM-2319.

## 7. A conditional unit-coloured incidence corollary

There is one exact arithmetic boundary on which the diagonal theorem
forces the still-missing edge incidence. This is a conditional corollary,
not a property of every live scalar row.

Work with the middle owner, so

```text
c_2=13^b u_2,
c_3=13^c u_3,                 c>b.
```

Put

```text
g=gcd(c_2,c_3),
a=c_2/g,
d'=c_3/g.                                           (29)
```

By the definition of `g`,

```text
gcd(a,d')=1,
13 divides d'.                                      (30)
```

The second fact uses `c>b`. Hence `13` does not divide `a`. Assume
additionally

```text
gcd(a,91)=1.                                        (31)
```

Condition (31) is explicit and load-bearing. It is automatic when seven
divides `d'`, by (30); otherwise it is exactly the requirement that the
unshared quotient `a` not contain a factor seven. It is not implied merely by
the existence of THM-2293's unit-coloured edge, because the endpoints of
that edge need not be integer multiples of `c_2`.

Set

```text
N=91d'.                                             (32)
```

Choose the representative `1<=K_0<d'` of a unit class modulo `d'` whose
residue modulo thirteen
gives any prescribed physical root character after multiplication by
`u_2`. Such a choice exists by CRT: prescribe the required nonzero
thirteen residue and a unit residue modulo every other prime-power factor
of `d'`.
For the integer representatives `z=0,...,90`, put

```text
K_z=K_0+d'z.                                        (34)
```

Then `1<=K_z<N`.
If seven divides `d'`, choose `K_0` nonzero modulo seven; every `K_z` is
then a unit modulo `N`. If seven does not divide `d'`, exactly one row
`z modulo 7` makes `K_z` divisible by seven; delete that row. At every
other prime of `N`, the residue of `K_z` is the fixed unit residue of
`K_0`. Thus all retained `K_z` are primitive modulo `N`, even when `d'`
has arbitrarily many additional prime factors.

Apply the primitive same-gauge theorem at modulus `N` to every retained
residue. Choose one positive common bare/word atom in that residue and
write it as

```text
q_z=K_z+N h_z,             0<=h_z<=J_(H_Q)J_(H_E)-1. (35)
```

Join two parameters when

```text
z_2-z_1 in (Z/91Z)^*.                               (36)
```

By CRT, the full graph when seven divides `d'` is

```text
Cay(Z/91Z,(Z/91Z)^*) = K_7 x K_13,                 (37)
```

the categorical product: two vertices are adjacent exactly when both CRT
coordinates differ. If seven does not divide `d'`, deleting the one
forbidden row leaves `K_6 x K_13`. Either graph has an edge. Choose any
adjacent retained pair. Equations (32), (34), and (35) give

```text
q_(z_2)-q_(z_1)
 =d'[(z_2-z_1)+91(h_(z_2)-h_(z_1))].                (38)
```

The bracket in (38) is coprime to `91` by adjacency. Multiplying by
`c_2=ga` and using `c_3=gd'` gives

```text
c_2 q_(z_2)-c_2 q_(z_1)
 =m c_3,

m=a[(z_2-z_1)+91(h_(z_2)-h_(z_1))],
gcd(m,91)=1.                                        (39)
```

Both endpoints in (39) are bare-source coefficients and literal
positive-word coefficients, at the same owner grade and prescribed root
character. Thus (39) is a THM-2293-type unit-coloured `c_3` edge whose
two endpoints are both marked by the same positive word.

Put `L=J_(H_Q)J_(H_E)<=12S^2`. Each selected atom lies between `1` and
`NL-1`. Hence the bracket `B` in (39) satisfies

```text
0<|B|=|q_(z_2)-q_(z_1)|/d'
     <91L,

|B|<=1092S^2-1,
|m|<=a(1092S^2-1).                                  (39a)
```

The normalized edge bound is independent of the size and prime
factorization of `d'`; only the unavoidable quotient `a` remains.

This argument uses a genuine graph but not a tournament. Its intrinsic
binary relation is unit difference in the CRT fibre, which is symmetric;
the coherent lift retains precisely the integer coordinate lost by reducing
a Fourier address modulo `N`.

The hypothesis (31) is also the exact failure boundary of this normalized
sublattice. If two `c_2`-multiple vertices satisfy

```text
c_2(q_2-q_1)=m c_3,
```

then coprimality of `a,d'` forces

```text
d' divides q_2-q_1,
m=a(q_2-q_1)/d'.                                   (40)
```

Thus every such multiplier is divisible by `a`. If `7` divides `a`, no
edge between two `c_2`-multiple vertices can have multiplier coprime to
`91`. A different affine coset or a non-`c_2`-multiple shell atom is then
genuinely necessary.

## 8. Exact remaining boundary

Equation (27) alone marks one vertex and does not prove that this particular
vertex is an endpoint of THM-2293's edge

```text
A-A'=m c_3,                 gcd(m,91)=1.             (41)
```

THM-2302's endpoint recurrence gives every marked vertex an incident
`c_3`-multiple edge only after forgetting the condition `gcd(m,91)=1`.
THM-2293 supplies a unit-coloured edge somewhere in the same character
graph, not necessarily at the vertex in (27). These are distinct
quantifiers at modulus `91`. Section 7 resolves the incidence by exposing
the larger modulus `91d'` whenever (31) holds. When seven divides `a`,
equation (40) proves that the complete `c_2`-multiple carrier is the wrong
affine lattice for a unit-coloured edge.

Nor does the fixed-colour correlation retain the target-plane gain or the
relative terminal-component phase of THM-2303. The complete connection
contract is

```text
source:
  THM-2319's separate primitive-colour spectra of a bare owner source and
  its literal positive-word subset;

target:
  a diagonal edge in the bipartite bare/word spectral-incidence graph;

map:
  freeze one primitive colour, form its physical cross-correlation, use
  rational Galois straightening plus the acute support sector to forbid
  cancellation, and land a syndetic gauge set by endpoint-product
  Vandermonde;

preserved:
  source owner, exact word, prescribed clock, primitive N-colour,
  ordinary Fourier frequency, same gauge index, owner grade, and at N=91
  the complete seven/thirteen unit condition; at N=91d' under (31), two
  word-marked endpoints and their unit-coloured c_3 edge;

destroyed or unselected:
  which terminal component carries the phase and exact target-plane gain;
  when seven divides a, the needed non-c_2-multiple affine shell coset;

needed sidecar:
  target-gain alignment and phase-tree transport on the edge in Section 7;
  on the seven-divisible-a branch, an affine-coset shift before applying
  the same fixed-colour mechanism;

cheapest decisive next test:
  couple the two word-marked endpoints in (39) to THM-2315's target-gain
  corolla, and separately test whether one legal owner/pivot translation
  escapes the exact seven-divisible obstruction (40).
                                                               (42)
```

The faithful object in this step is a bipartite spectral-incidence graph:
one side is the bare spectrum, the other is the word-restricted spectrum,
and edges mean equality of the full coloured gauge address. There is no
intrinsic pairwise comparison and therefore no tournament orientation.

No scalar profile is excluded. The exact scalar ledger remains `165`, and
LRC(14) remains open.

## 9. Exact companion

The companion verifies the `N=91` acute support sector, the arbitrary-prime
control `N=210` where the old totient-width test fails, exact construction
of the relevant cyclotomic polynomials, all primitive colours on
deterministic rational needles, the aggregate-zero hostile control (12)
together with fixed-colour nonvanishing, the
`12S^2-1` gauge bound, and the `1092S^2-1` frequency bound. Every
load-bearing check raises explicitly in ordinary and optimized Python.
It also checks the two CRT fibre graphs, primitive coherent lifts at
modulus `91d'` with a large composite unit cofactor, the factorization
(38)--(40), the normalized edge bound (39a), and both sides of the exact
seven-divisibility boundary.

Reproduce with

```bash
python3 04-computation/lrc14_primitive_cross_correlation_same_gauge_thm2323.py
python3 -O 04-computation/lrc14_primitive_cross_correlation_same_gauge_thm2323.py
```

The two transcripts must match

```text
05-knowledge/results/lrc14_primitive_cross_correlation_same_gauge_thm2323.out
```

byte-for-byte after LF normalization.
