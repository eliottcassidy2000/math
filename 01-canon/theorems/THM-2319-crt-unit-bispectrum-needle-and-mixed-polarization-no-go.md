---
id: THM-2319
title: "CRT unit-bispectrum needle and mixed-polarization no-go"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. On Z/91,
  every nonzero rational root vector supported on thirteen consecutive
  sites has all 72 primitive-character transforms nonzero. For nonnegative
  vectors, the whole bispectrum restricted to unit characters k,l,k+l is
  at least (1305/98) times the cube of the mass. Its carry-zero and
  carry-one halves are conjugate and equally positive, so one carry-one
  pair has real cubic current at least (29/8624) times mass cubed. For a
  rational nonnegative step function supported in D_1, every unit colour
  lands at bounded Fourier height, and the cubic current gives a bounded
  Schur triangle A+B=C with all three terms coprime to 91. Applied
  separately to the bare LRC owner source and its literal positive-word
  subset, this supplies all unit colours at the base owner valuation; for
  the middle owner it gives an actual grade-b bare vertex in any prescribed
  shell character.
  The consecutive-needle geometry is load-bearing: an exact nonconsecutive
  seven-site vector has negative unit bispectrum. More decisively for
  LRC(14), even on eight consecutive sites the polarization with two bare
  source legs and one word-restricted leg can be negative. Thus the theorem
  supplies an all-unit word-current triangle but not its incidence with a
  bare-source shell edge. No scalar row is excluded and LRC(14) remains
  open.
source: codex-2026-07-25-crt-unit-bispectrum
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
related:
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2306-owner-normalized-disjoint-support-and-first-collision-shell
  - THM-2312-sparse-root-bispectrum-positive-word-current
  - THM-2313-biprime-pareto-collision-frontier-and-91-unit-current-shell
  - THM-2318-one-shot-three-prime-mobius-amplifier
script: 04-computation/lrc14_crt_unit_bispectrum_needle_thm2319.py
output: 05-knowledge/results/lrc14_crt_unit_bispectrum_needle_thm2319.out
script_sha256: 65ad1e46ab39b7a99ce97c7e4b611a37a5bb2d0f61c3b9b88b2e886bbcbe589d
output_sha256: cd8e23fb29ed99a518ca10226ad26ad81b101b91e17c0283d0ae23176893611e
hash_basis: working-tree bytes (LF)
---

# THM-2319 -- a thirteen-tooth needle has a positive 91-unit face

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2312 makes every nonzero thirteenth-root character available on one
positive blocker word, but its charge-zero contraction does not force the
septimal unit color required by THM-2302. THM-2313 instead obtains a
`91`-unit residual by a two-prime collision derivative, at the cost of two
uncontrolled valuation depths.

There is a third operation: expose all `91` roots in one shot and sum only
the character face which is nonzero modulo both seven and thirteen. The
owner-normalized interval has length

```text
1/7=13/91.
```

Every active `91`-root fibre is therefore a thirteen-consecutive-site
needle. Rationality already makes every primitive character nonzero, and
its CRT geometry makes the entire unit face quantitatively positive.

The gain is real but sharply scoped. Applying the theorem after restricting
to one blocker word puts all three legs on the word-restricted function.
The mixed polarization needed to replace two of those legs by bare-source
coefficients is false.

## 1. The whole unit-character face

Put

```text
N=91=7*13,
U=(Z/NZ)^*,
zeta=exp(2*pi*i/N).
```

For a real vector `v=(v_r)_(r in Z/NZ)`, define

```text
M_k(v)=sum_r v_r zeta^(-kr),                         (1)
```

and its whole unit bispectrum

```text
B_U(v)
 =sum_(k,l in U; k+l in U)
    M_k(v)M_l(v)conjugate(M_(k+l)(v)).              (2)
```

All character indices in (2) are modulo `91`.

For a prime `p`, define the local kernel

```text
K_p(A,B)
 =sum_(x,y nonzero mod p; x+y nonzero mod p)
    exp(2*pi*i(xA+yB)/p).                            (3)
```

Inclusion-exclusion on the three forbidden lines gives

```text
K_p(A,B)
 =p^2 [A=0,B=0]
  -p([A=0]+[B=0]+[A=B])
  +2.                                                (4)
```

Consequently its value depends only on the equality pattern of three
physical sites `r,s,t`, with `A=t-r`, `B=t-s`:

```text
r=s=t:                 (p-1)(p-2);
exactly two equal:      2-p;
all three distinct:     2.                           (5)
```

The Chinese remainder theorem factors both the unit conditions in (2) and
the character kernel. Expanding the three Fourier factors therefore gives

```text
B_U(v)
 =sum_(r,s,t)
   v_r v_s v_t
   K_7(t-r,t-s)K_13(t-r,t-s).                       (6)
```

In particular, `B_U(v)` is real. Formula (6), rather than an unsigned
absolute-value estimate, is the exact CRT current.

## 2. Quantitative positivity on a thirteen-consecutive needle

Assume that `v_r>=0`, that `v` is nonzero, and that its support lies in a
cyclic interval of thirteen consecutive residues. Translation invariance
of (2) lets us use the sites

```text
0,1,...,12.                                         (7)
```

They are pairwise distinct modulo thirteen. Modulo seven they occupy six
two-point fibres and one one-point fibre. For each septimal fibre write its
two possible masses as

```text
a_i,b_i>=0,
X_i=a_i+b_i,
V=sum_i X_i=sum_r v_r,             i in Z/7Z,       (8)
```

putting `b_i=0` on a one-point fibre.

Equations (5)--(6) give the complete coefficient table relevant to (7):

| physical-site pattern | mod-seven pattern | CRT coefficient |
|---|---|---:|
| all equal | all equal | `30*132=3960` |
| exactly two sites, same seven-fibre | all equal mod 7 | `30*(-11)=-330` |
| exactly two sites, different seven-fibres | exactly two equal mod 7 | `(-5)*(-11)=55` |
| three sites, one same-fibre pair | exactly two equal mod 7 | `(-5)*2=-10` |
| three sites, three seven-fibres | all distinct mod 7 | `2*2=4` |

The positive cross-fibre classes with coefficients `55` and `4` may be
discarded in a lower bound. For a fixed septimal fibre, the negative terms
are the six permutations of a same-fibre pair and one exterior site. Hence

```text
B_U(v)
 >=sum_i [
    3960(a_i^3+b_i^3)
   -990a_i b_i X_i
   -60a_i b_i(V-X_i)
   ].                                                (9)
```

Now

```text
a_i^3+b_i^3=X_i^3-3a_i b_i X_i,
a_i b_i<=X_i^2/4.
```

The first two terms in the bracket in (9) are at least

```text
3960X_i^3-12870(X_i^3/4)
 =(1485/2)X_i^3.                                    (10)
```

The last term is at least `-15X_i^2(V-X_i)`. Therefore

```text
B_U(v)
 >=15[(101/2)sum_i X_i^3-V sum_i X_i^2].            (11)
```

Put `Y=sum_i X_i^2`. Cauchy--Schwarz and the seven-fibre bound give

```text
Y^2<=V sum_i X_i^3,
Y>=V^2/7.                                           (12)
```

Writing `z=Y/V^2`, the right side of (11) is at least

```text
15V^3[(101/2)z^2-z].
```

This quadratic is increasing for `z>=1/7`. At `z=1/7` it equals `87/98`.
We have proved:

> **CRT unit-needle theorem.**
>
> ```text
> B_U(v)>=(1305/98)V^3>0.                           (13)
> ```

The displayed constant is a uniform certified floor, not a sharpness
claim. The proof deliberately spends every positive cross-fibre term.

There are

```text
(7-1)(7-2)(13-1)(13-2)=30*132=3960                (14)
```

ordered pairs in (2). Thus one satisfies

```text
Re[M_kM_l conjugate(M_(k+l))]
 >=(29/8624)V^3.                                   (15)
```

All three characters in (15) are coprime to `91`.

## 3. Carry one is equally positive

Represent `k,l` by integers in `{1,...,90}` and put

```text
c=(k+l) mod 91 in {1,...,90},
q=(k+l-c)/91 in {0,1}.                              (16)
```

The involution

```text
(k,l)->(91-k,91-l)                                  (17)
```

maps the `q=0` pairs bijectively to the `q=1` pairs. Since `v` is real, it
conjugates their cubic terms. The whole sum in (13) is real, so each half
has real sum `B_U(v)/2`. Both halves contain `1980` pairs. Consequently one
**carry-one** pair obeys the same bound (15).

This is stronger than merely finding some unit character. In the periodic
gauge

```text
N_k(y)=exp(-2*pi*i*k*y/91)M_k(y),                   (18)
```

the selected current has the exact spatial factor

```text
M_kM_l conjugate(M_c)
 =exp(2*pi*i*y)N_kN_l conjugate(N_c).               (19)
```

The carry is an integer frequency one. It does not by itself identify a
THM-2293 shell multiplier.

## 4. A bounded all-unit Schur triangle for interval-supported steps

Let `H` be a nonzero nonnegative real step function on the circle with

```text
support(H) subset D_1={x:||x||<1/14},
rho=integral H>0,                                   (20)
```

and let `J` be its number of nonzero jumps. Define its `91`-root fibres

```text
v_r(y)=H((y+r)/91),
M_k(y)=sum_(r=0)^90 v_r(y)zeta^(-kr).               (21)
```

Every fibre in (21) is supported on at most thirteen consecutive roots.

### 4.1 Rational primitive-character saturation

Suppose first that `H` is rational-valued. More generally, let `v` be any
nonzero rational vector supported on at most thirteen consecutive cyclic
sites. Translate those sites to `0,...,12` and put

```text
P_v(X)=sum_(r=0)^12 v_r X^r in Q[X].                (21a)
```

For every `k in U`, the number `zeta^(-k)` is a primitive `91`st root. Its
minimal polynomial over `Q` is `Phi_91`, of degree

```text
phi(91)=72.                                         (21b)
```

The nonzero polynomial in (21a) has degree at most twelve, so it cannot
vanish at `zeta^(-k)`. Therefore

```text
M_k(v)!=0                         for every k in U. (21c)
```

Nonnegativity is not needed for (21c); rationality and the
thirteen-consecutive support are the load-bearing hypotheses.

Apply this pointwise to (21). Wherever its root fibre is nonzero, every one
of the `72` primitive-character functions `M_k(y)` is nonzero. Consequently
each periodic gauge `N_k` from (18) is a nonzero function. If `H` has `J`
jumps, its step amplitude `M_k` has at most `J` image-jump nodes.
Distributional differentiation gives a nonzero exponential sequence

```text
2*pi*i(h+k/91)(N_k)_hat(h).
```

It cannot vanish at `J` consecutive nonnegative integers. Hence, for every
unit colour `k`, some

```text
0<=h_k<=J-1
```

satisfies

```text
H_hat(k+91h_k)!=0,
k+91h_k<=91J-1,
gcd(k+91h_k,91)=1.                                 (21d)
```

Thus primitive-character saturation is not only fibrewise: all 72 colours
land at bounded ordinary Fourier addresses.

### 4.2 Cubic carry and the Schur triangle

Apply Section 3 pointwise and integrate. If

```text
V(y)=sum_r v_r(y),
```

then

```text
integral V=91rho,
integral V^3>=(91rho)^3.                            (22)
```

One carry-one unit pair therefore has

```text
Re integral M_kM_l conjugate(M_c)
 >=(29/8624)(91rho)^3>0.                            (23)
```

The current in (23) forces an exact additive triangle in ordinary Fourier
support. Indeed, the gauge in (18) is periodic and

```text
(N_k)_hat(h)=91 H_hat(k+91h).                       (24)
```

Fejer approximation of the three bounded periodic factors in (23) shows
that some integers `h_1,h_2,h_3` satisfy

```text
h_3=h_1+h_2+1,
(N_k)_hat(h_1)(N_l)_hat(h_2)
 conjugate((N_c)_hat(h_3))!=0.                     (25)
```

There is a finite nonnegative landing. Distributional differentiation
writes

```text
2*pi*i(h+k/91)(N_k)_hat(h)
```

as an exponential sum over at most `J` image-jump nodes. The bivariate
product of the three differentiated sequences in (25) is a nonzero
exponential polynomial in `(h_1,h_2)`. It uses at most `J^2` distinct
first-coordinate nodes and at most `J^2` distinct second-coordinate nodes.
The tensor product of the two Vandermonde matrices shows that it cannot
vanish on the whole grid

```text
0<=h_1,h_2<=J^2-1.                                  (26)
```

Thus (25) holds with (26), and

```text
0<=h_3<=2J^2-1.                                     (27)
```

Put

```text
A=k+91h_1,
B=l+91h_2,
C=c+91h_3.                                         (28)
```

Equations (16), (25), and (27) give

```text
A+B=C,
H_hat(A)H_hat(B)conjugate(H_hat(C))!=0,
gcd(A B C,91)=1,                                   (29)

A,B<=91J^2-1,
C<=182J^2-1.                                       (30)
```

This is a bounded positive Schur triangle in the Fourier support of the
same interval-supported step function. It is an additive three-leg
relation, not a tournament orientation.

## 5. Exact LRC word-current specialization

Use THM-2305's positive word `Q` for owner `j`. Write

```text
c_j=13^lambda u_j,             13 does not divide u_j,
k_j=lambda+1,

E_Q=E_j intersection T^(-k_j)Q,
rho=measure(E_Q)>0.                                  (31)
```

Thus `E_Q` is the literal subset of the exclusive source whose prescribed
orbit lands in the complete terminal word. Put

```text
G=P^lambda 1_(E_j),
F_Q=G (1_Q after T).                                 (32)
```

The word factor in (32) is constant on every branch in the `P^lambda`
average: if `T^lambda z=x`, then

```text
1_Q(T^(lambda+1)z)=1_Q(Tx).
```

Therefore

```text
F_Q=P^lambda 1_(E_Q).                                (33)
```

This identity is load-bearing. Normalize by the remaining owner unit:

```text
H_Q=P_(u_j)F_Q=P_(c_j)1_(E_Q).                      (34)
```

It follows exactly that

```text
support(H_Q) subset D_1,
integral H_Q=rho,
(H_Q)_hat(n)=(1_(E_Q))_hat(c_j n).                  (35)
```

There is also a bare-source normalization

```text
H_E=P_(c_j)1_(E_j).                                 (35a)
```

Both `H_E` and `H_Q` are nonzero rational step functions supported in
`D_1`. Apply primitive-character saturation (21c)--(21d) separately. For
every `k in U`, there are indices

```text
0<=h_E<=2S-1,
0<=h_Q<=6S-1                                       (35b)
```

such that

```text
(1_(E_j))_hat(c_j(k+91h_E))!=0,
(1_(E_Q))_hat(c_j(k+91h_Q))!=0.                    (35c)
```

The first jump bound uses `#Jump(H_E)<=2S`; the second is proved below.
Both atoms have the base owner valuation and the same prescribed CRT
colour, although their gauge indices need not agree. In particular, given
any nonzero shell character `kappa modulo 13`, choose a unit `k modulo 91`
with

```text
u_j k congruent kappa mod 13.                       (35d)
```

There are six such CRT choices, one for each nonzero residue modulo seven.
Then the first atom in (35c) is a bare-source vertex of root character
`kappa` with a `91`-unit outside multiplier. If `j=2`, so
`nu_13(c_j)=b`, it is an actual grade-`b` bare vertex in the prescribed
THM-2293 shell-character graph.

Quantitatively,

```text
k+91h_E<=182S-1,
k+91h_Q<=546S-1.                                   (35e)
```

This removes “produce a prescribed-character bare vertex with a `91`-unit
outside multiplier” from the open target. What is not proved is a common
index `h_E=h_Q`, or incidence of either atom with the required unit-coloured
`c_3` edge.

The quantitative current in (23) now becomes

```text
sum_(all unit pairs) Re C_(k,l)
 >=(1305/98)(91rho)^3
  =(20069595/2)rho^3.                               (36)
```

The carry-one half contains `1980` pairs, so one carry-one pair satisfies

```text
Re C_(k,l)>=(445991/176)rho^3>0.                    (37)
```

Section 4 supplies three nonzero coefficients of the literal word-source
indicator:

```text
(1_(E_Q))_hat(c_j A)
(1_(E_Q))_hat(c_j B)
conjugate((1_(E_Q))_hat(c_j C))!=0,

c_j A+c_j B=c_j C,
gcd(A B C,91)=1.                                    (38)
```

Thus all three atoms retain the same exact source, prescribed clock, and
complete word ancestry. Their thirteen-adic grade is the base owner grade
`lambda`, and their outside multipliers are simultaneous seven- and
thirteen-units.

There is also a safe uniform height ledger. THM-2305 gives

```text
#Jump(G)<=2S,
#Jump(1_Q)<=2S.
```

Only the pullback jumps which meet `support(G)` can survive in the product
(32). For each of the at most `2S` jumps of `1_Q`, its thirteen preimages
under `T` contain at most two points in `closure(D_(u_j))`: multiplication
by the thirteen-unit `u_j` permutes the roots, and three roots would span
at least `2/13>1/7`. Thus the word factor contributes at most `4S`
effective product jumps. Together with the `2S` jumps of `G`, equation (32)
has at most `6S` jumps. Perron transport cannot increase this count. In
Section 4 one may therefore take

```text
J<=6S,

A,B<=3276S^2-1,
C<=6552S^2-1.                                      (39)
```

This improves THM-2312's root-only character pair to a literal `91`-unit
additive Fourier triangle at base owner valuation. It still does **not**
make any leg a coefficient of the full bare source `1_(E_j)`, much less a
second endpoint of a THM-2293 shell edge. Restriction to `E_Q` can create
Fourier atoms absent from `E_j`.

## 6. Mixed polarization is false

The missing replacement cannot follow from positivity by formal
polarization. For real root vectors `g,f`, define the mixed face

```text
B_U(g,f,g)
 =sum_(k,l,k+l in U)
    M_k(g)M_l(f)conjugate(M_(k+l)(g)).              (41)
```

On the eight consecutive sites `0,...,7`, take

```text
g=(62,842,894,609,819,261,207,761),
f=(62,0,0,0,0,0,0,0).                              (42)
```

Thus `0<=f<=g`; `f` is literally a one-site word restriction of `g`.
Exact evaluation of (6), with the middle mass replaced by `f_s`, gives

```text
B_U(g,f,g)=-3042309124<0.                           (43)
```

The unmixed current on the same source is strongly positive:

```text
B_U(g)=12219375808488
      >(1305/98)(4455)^3.                           (44)
```

So neither the unit restriction, the interval needle, nonnegativity, nor a
literal subfunction makes the required two-source/one-word current
positive.

The consecutive geometry in the positive theorem is also essential. On
the seven nonconsecutive sites

```text
R=(53,79,37,27,40,58,1) mod 91
```

with positive masses

```text
w=(283,57,16,270,174,4,196),
```

the unmixed unit face itself is

```text
B_U(w)=-4459059900<0.                               (45)
```

This is not a failure of the CRT formula. It is precisely the information
carried by the thirteen-consecutive Kakeya needle.

## 7. Connection and stopping boundary

The exact connection contract is:

```text
source:
  a nonnegative function on one owner-normalized interval, optionally
  already restricted to one complete THM-2305 blocker word;

map:
  expose 91 roots, sum the whole k,l,k+l unit face, split by radix carry,
  and use endpoint-Vandermonde landing;

preserved:
  every primitive colour for rational inputs, nonnegativity for the cubic
  floor, the exact word when present, all three unit residues, cubic
  complex multiplication, carry one, and a bounded Schur triangle;

destroyed or unselected:
  a uniform cubic character pair, a common gauge index between the bare
  and word-restricted spectra, c_3-edge incidence, terminal component
  phase, target-plane gain, and orbit holonomy;

sharp hostile mechanisms:
  the eight-consecutive mixed polarization (42) and the nonconsecutive
  seven-site current (45);

needed sidecar:
  a same-index theorem coupling (35c), followed by a unit-coloured
  c_3-edge incident to that common marked vertex; alternatively a
  three-prime collision cell which preserves the same colour and grade. (46)
```

The theorem answers the `Fano/chi_7` probe at the current level: seven is
useful through the full CRT equality kernel, not through a quadratic
character orientation. The resulting object is a marked additive triangle
inside a thirteen-tooth needle.

No scalar profile is excluded, and LRC(14) remains open.

## 8. Exact companion

The companion checks the prime kernels directly, the CRT coefficient
table, the `3960=1980+1980` pair counts, every rational constant in
(9)--(15), deterministic positive controls, the exact mixed hostile
(35)--(37), the nonconsecutive hostile (38), and the endpoint landing
bounds. Every load-bearing test raises explicitly under ordinary and
optimized Python.

Reproduce with

```bash
python3 04-computation/lrc14_crt_unit_bispectrum_needle_thm2319.py
python3 -O 04-computation/lrc14_crt_unit_bispectrum_needle_thm2319.py
```

The two transcripts must match

```text
05-knowledge/results/lrc14_crt_unit_bispectrum_needle_thm2319.out
```

byte-for-byte after LF normalization.
