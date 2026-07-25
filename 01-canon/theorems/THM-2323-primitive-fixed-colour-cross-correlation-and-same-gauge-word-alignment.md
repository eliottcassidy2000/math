---
id: THM-2323
title: "Primitive fixed-colour cross-correlation and same-gauge word alignment"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. Let N be
  divisible by seven and satisfy phi(N)>2(N/7-1). If 0<=f<=g are
  nonzero rational-valued step functions with rational breakpoints,
  supported in one open circle arc of length 1/7, then for every primitive
  N-character k every block of J_f J_g consecutive gauge indices contains
  one at which both Fourier coefficients f-hat(k+Nh) and
  g-hat(k+Nh) are nonzero. In particular there is a common nonnegative
  index 0<=h<=J_f J_g-1. The mechanism is a fixed-colour
  cross-correlation Laurent polynomial of width 2(N/7-1), whose positive
  diagonal coefficient and rationality prevent primitive-root vanishing;
  Parseval and endpoint-product Vandermonde then land the common index.
  The hypothesis holds for the complete stalk N=7*13^a, a>=1. At N=91,
  applying the theorem to THM-2319's word source H_Q and bare source H_E
  gives, for every one of the 72 unit colours, a common index
  h<=12S^2-1 and common multiplier n<=1092S^2-1. This closes the
  bare/word same-gauge-index loss. Conditionally, for the middle owner,
  write c_2/gcd(c_2,c_3)=a and
  c_3/gcd(c_2,c_3)=7^alpha*13^delta*r with gcd(r,91)=1. If
  gcd(a,91)=1 and either alpha>=1,r<=6 or alpha=0,r<=5, a
  K_7 x K_13 / K_6 x K_13 colouring argument forces two common
  bare/word atoms to form a unit-coloured c_3 edge. Outside this explicit
  arithmetic branch, the theorem does not force edge incidence, select
  target-plane gain, or exclude a scalar row. LRC(14) remains open.
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
script_sha256: a62af2e409d87570634e8c069652be54f9ec2d3f1f345b7c1c719fdc1e700649
output_sha256: ad146cf6eb5bafd1eb8301a1a2cbb67e05046798ed3219c8868ae28fc7e36d24
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
rational Laurent polynomial. Its width is smaller than the cyclotomic
degree, so it cannot vanish. This is a phase-sensitive algebraic
nonvanishing statement, not a lower bound on the real part.

## 1. Root fibres and the fixed-colour cross-correlation

Let `N` be a positive multiple of seven, put

```text
R=N/7-1,
```

and assume

```text
phi(N)>2R.                                           (1)
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
P(X)=sum_(d=-R)^R c_d X^(d+R) in Q[X]              (9)
```

is nonzero and has degree at most `2R`.

## 2. Cyclotomic width forbids fixed-colour cancellation

Suppose `C_k=0`. Since `zeta^k` is a primitive `N`th root, (7)--(9) give

```text
P(zeta^k)=zeta^(kR)C_k=0.
```

The minimal polynomial of `zeta^k` over `Q` is `Phi_N`, of degree
`phi(N)`. Therefore `Phi_N` divides `P`. This contradicts

```text
deg(P)<=2R<phi(N)                                   (10)
```

and the fact that `P` is nonzero. We have proved:

> **Fixed-colour cross-correlation lemma.** Under (1)--(2),
>
> ```text
> C_k!=0 for every k in (Z/NZ)^*.                    (11)
> ```

The rationality in (2) is load-bearing: it places `P` in `Q[X]`, so the
full cyclotomic minimal polynomial is forced by one primitive zero.
Nonnegativity is used at the diagonal (8). The common short support is
load-bearing at (6).

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

but every fixed summand is nonzero: the rational polynomial
`1+12X^7` has degree seven, below `deg Phi_91=72`. The colour sum can
cancel even though no colour does.

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

## 5. The all-depth `7*13^a` stalk

Take

```text
N_a=7*13^a,                 a>=1.
```

Then

```text
R_a=13^a-1,
phi(N_a)=72*13^(a-1),

phi(N_a)-2R_a
 =46*13^(a-1)+2
 >0.                                                (21)
```

Thus (1) holds at every depth. For every primitive `N_a`-colour, the same
pair `f<=g` has a common gauge satisfying (19), with frequency at most
`N_a J_fJ_g-1`.

This is a genuine stalk law: increasing the thirteen depth enlarges the
physical root needle from `13` to `13^a` teeth, while the cyclotomic degree
grows fast enough to dominate twice its correlation width at every level.
At each fixed depth and colour the common spectrum is a bounded-gap
toothpick. It does not assert compatibility of the selected gauge indices
between different depths.

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
2R=24<72=phi(91).
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

There is one arithmetic range in which the diagonal theorem does force
the still-missing edge incidence. This is a conditional corollary, not a
property of every live scalar row.

Work with the middle owner, so

```text
c_2=13^b u_2,
c_3=13^c u_3,                 c>b.
```

Put

```text
g=gcd(c_2,c_3),
a=c_2/g,
d'=c_3/g,

d'=7^alpha 13^delta r,
gcd(r,91)=1.                                        (29)
```

By the definition of `g`, one has `gcd(a,d')=1`, and here `delta>=1`.
Assume additionally

```text
gcd(a,91)=1,                                        (30)
```

and one of the two small-cofactor conditions

```text
alpha>=1 and r<=6,

or

alpha=0 and r<=5.                                   (31)
```

Condition (30) is explicit and load-bearing. It is not implied merely by
the existence of THM-2293's unit-coloured edge, because the endpoints of
that edge need not be integer multiples of `c_2`.

Set

```text
g_0=7^alpha 13^delta,
N=91g_0.                                            (32)
```

This modulus satisfies the abstract width hypothesis:

```text
phi(N)=72g_0,
2(N/7-1)=26g_0-2,
phi(N)-2(N/7-1)=46g_0+2>0.                          (33)
```

Choose an integer `K_0`, prime to `g_0`, whose residue modulo thirteen
gives any prescribed physical root character after multiplication by
`u_2`. For `z in Z/91Z`, use the coherent lift

```text
K_z=K_0+d'z.                                        (34)
```

If `alpha>=1`, every `K_z` is a unit modulo `N`: its residues modulo seven
and thirteen are the fixed nonzero residues of `K_0`. If `alpha=0`, its
residue modulo thirteen is still fixed and nonzero, while exactly one
row `z modulo 7` makes it divisible by seven. Delete that row.

Apply the primitive same-gauge theorem at modulus `N` to every retained
residue. Choose one positive common bare/word atom in that residue and
write it relative to the coherent lift as

```text
q_z=K_z+N H_z,                  H_z in Z.            (35)
```

The integer `H_z` may be negative because `K_z` is a coherent lift rather
than the canonical representative; the atom `q_z` itself is positive and
obeys the bound (20).

Join two parameters when

```text
z_2-z_1 in (Z/91Z)^*.                               (36)
```

By CRT, the full graph in the `alpha>=1` case is

```text
Cay(Z/91Z,(Z/91Z)^*) = K_7 x K_13,                 (37)
```

the categorical product: two vertices are adjacent exactly when both CRT
coordinates differ. Its chromatic number is seven. In the `alpha=0` case,
deleting one complete mod-seven row leaves `K_6 x K_13`, whose chromatic
number is six. The upper bounds colour by the first coordinate, and the
matching diagonal gives cliques of sizes seven and six, respectively.

Colour each retained vertex by

```text
H_z mod r.                                          (38)
```

Under (31), fewer colours are available than the chromatic number.
Therefore two adjacent vertices have the same colour. For their coherent
lifts,

```text
H_(z_2)-H_(z_1)=rw
```

for some integer `w`, and (29), (32), (34), and (35) give

```text
q_(z_2)-q_(z_1)
 =d'[(z_2-z_1)+91w].                                (39)
```

The bracket in (39) is coprime to `91` by adjacency. Multiplying by
`c_2=ga` and using `c_3=gd'` gives

```text
c_2 q_(z_2)-c_2 q_(z_1)
 =m c_3,

m=a[(z_2-z_1)+91w],
gcd(m,91)=1.                                        (40)
```

Both endpoints in (40) are bare-source coefficients and literal
positive-word coefficients, at the same owner grade and prescribed root
character. Thus (40) is a THM-2293-type unit-coloured `c_3` edge whose
two endpoints are both marked by the same positive word.

This argument uses a genuine graph but not a tournament. Its intrinsic
binary relation is unit difference in the CRT fibre, which is symmetric;
the chromatic obstruction retains precisely the lift coordinate lost by
reducing a Fourier address modulo `N`.

## 8. Remaining incidence boundary

Equation (27) marks one vertex. It does not prove that the vertex is an
endpoint of THM-2293's edge

```text
A-A'=m c_3,                 gcd(m,91)=1.             (41)
```

THM-2302's endpoint recurrence gives every marked vertex an incident
`c_3`-multiple edge only after forgetting the condition `gcd(m,91)=1`.
THM-2293 supplies a unit-coloured edge somewhere in the same character
graph, not necessarily at the vertex in (27). These are distinct
quantifiers in general. Section 7 combines them only under the explicit
arithmetic hypotheses (30)--(31).

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
  support width plus cyclotomic degree to forbid cancellation, and land a
  nonnegative gauge by endpoint-product Vandermonde;

preserved:
  source owner, exact word, prescribed clock, primitive N-colour,
  ordinary Fourier frequency, same gauge index, owner grade, and at N=91
  the complete seven/thirteen unit condition;

destroyed or unselected:
  which terminal component carries the phase, target-plane gain, a second
  shell vertex, and unit-coloured c_3-edge incidence;

needed sidecar:
  outside Section 7's low-cofactor branch, a unit-coloured marked-degree
  theorem at the common vertex, or an exact third-order coefficient
  coupling that vertex to a THM-2293 shell edge;

cheapest decisive next test:
  rerun the THM-2302 seven-periodic hostile mechanism with the newly forced
  marked frequency coprime to 91, and classify the endpoint-orbit patterns
  which could still concentrate every incident multiplier on 7Z or 13Z.
                                                               (42)
```

The faithful object in this step is a bipartite spectral-incidence graph:
one side is the bare spectrum, the other is the word-restricted spectrum,
and edges mean equality of the full coloured gauge address. There is no
intrinsic pairwise comparison and therefore no tournament orientation.

No scalar profile is excluded. The exact scalar ledger remains `165`, and
LRC(14) remains open.

## 9. Exact companion

The companion verifies the `N=91` support width and cyclotomic degree, the
all-depth identity (21), exact construction of `Phi_91`, all `72`
primitive colours on a deterministic rational needle, the aggregate-zero
hostile control (12) together with fixed-colour nonvanishing, the
`12S^2-1` gauge bound, and the `1092S^2-1` frequency bound. Every
load-bearing check raises explicitly in ordinary and optimized Python.
It also checks the two CRT fibre graphs and their chromatic thresholds,
the coherent-lift factorization (39)--(40), and the exact modulus gap (33).

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
