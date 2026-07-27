---
source: codex-2026-07-24-affine-unit-bispectrum
status: >
  PROOF CANDIDATE + VERIFIED-EXACT CONSTANTS + INDEPENDENTLY AUDITED.
  On the strictly shallower selected-owner branch of THM-2302/2305, a
  nonnegative whole-character cubic face and the THM-2293 shell extractor
  force three actual Fourier atoms of one literal handoff subset whose
  signed sum is a gcd(m,91)=1 multiple of the deepest blocker. This is a
  unit-coloured affine-covariant third-order coefficient, not an ordinary
  pair edge, a THM-2303 component-phase edge, or a proof of LRC(14).
depends_on:
  - THM-2293-quadratic-root-energy-raises-the-ancestry-shell
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2312-sparse-root-bispectrum-positive-word-current
related:
  - THM-2286-endpoint-prony-lift-and-phase-no-go
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2304-deepest-boundary-cyclotomic-current-separation
  - THM-2306-owner-normalized-disjoint-support-and-first-collision-shell
---

# A literal handoff subset carries an affine unit bispectrum

## 1. Inheritance and the corrected target

THM-2302 reaches the exact boundary

```text
same owner + same root colour + prescribed service
  + a marked source atom
  + an uncoloured incident terminal edge,
```

but endpoint recurrence cannot force the load-bearing colour

```text
gcd(m,91)=1.                                         (1)
```

THM-2312 supplies a positive whole-character cubic face on one exact
THM-2305 word, but only at base frequency zero. It deliberately forgets
common affine phase and the multiplier `m`.

The two mechanisms compose in the opposite order:

```text
sum the complete root-character face
  -> obtain a nonnegative cubic density
  -> apply the deepest-comb covariance shell extractor
  -> select a nonzero affine base mode
  -> only then select one character triad.            (2)
```

The result is exactly the third-order alternative permitted by THM-2302.
It is not necessary to deconvolve a restricted quadratic coefficient into
an unrestricted pair edge.

## 2. Strict shallow-owner setup

Fix a positive THM-2305 word for a selected shallow owner `j`. Write

```text
c_j=13^lambda u_j,
k=lambda+1,

Q=Q_(j,sigma),             sigma in {{a},{b},{a,b}},

E_Q=E_j intersection T^(-k)Q,

rho_Q=measure(E_Q)>0.                               (3)
```

Put

```text
G=P^lambda 1_(E_j),

v_r(y)=G((y+r)/13),

M_k(y)=sum_(r=0)^12 v_r(y) zeta^(-kr).              (4)
```

Every fibre in (4) has at most two positive entries, each at most one.
Assume that the selected owner is strictly shallower than the deepest
blocker:

```text
s=c_3/13^(lambda+1) in Z.                           (5)
```

Because `E_j` is disjoint from `D_(c_3)`, every active final root fibre
lies outside `D_s`. In particular every function below is supported in
`D_s^c`.

Indeed, if `B_Q(y)>0`, then some `z in E_j` satisfies

```text
T^(lambda+1)z=y.
```

Since `c_3=13^(lambda+1)s`,

```text
||s y||=||c_3 z||>=1/14.                            (5a)
```

Condition (5) is automatic on THM-2302's strict same-label branch. It is
not automatic from THM-2296/2305 alone: on a repeated-first row the
selected owner can be deepest. Every repeated-first statement below is
therefore conditional on the same deeper-exclusion hypothesis, or applies
after replacing `c_3` by a named strictly deeper blocker.

## 3. The positive cubic density

Define the word-restricted whole-character face

```text
B_Q(y)
 =1_Q(y)
  sum_(k,l nonzero; k+l nonzero mod 13)
    M_k(y) M_l(y) conjugate(M_(k+l)(y)).             (6)
```

THM-2312 gives pointwise

```text
B_Q
 =1_Q[13^2 S_3-3*13 S_1 S_2+2S_1^3],

B_Q>=(99/4)1_Q S_1^3>=0.                            (7)
```

For two sheet masses `0<=a,b<=1`,

```text
B_13(a,b)
 =11[12(a^3+b^3)-3ab(a+b)]<=198.                   (8)
```

The cap is sharp at `a=b=1`. For fixed `b`, the sole nonnegative interior
critical point in `a` is a minimum, so the maximum reduces to the four
corners, with values `0,132,132,198`.

Every word in THM-2305 lies in a blocker danger set and hence has measure
at most `1/7`. Since

```text
integral_Q S_1=13 rho_Q,
```

Holder's inequality and (7) give

```text
mu_Q:=integral B_Q
 >=(10657647/4)rho_Q^3.                              (9)
```

## 4. The shell extractor supplies a unit base frequency

Apply THM-2293's conditional-covariance shell lemma to `B_Q`, with
pointwise cap `C=198`. Since

```text
B_Q 1_(D_s)=0,
```

its Fourier interaction with the deepest comb is a strictly negative sum
over exactly

```text
d=ms,             m!=0,             gcd(m,91)=1.    (10)
```

The tail estimate in THM-2293 says that a nonzero term occurs with

```text
0<|m|<=floor[33800*198/(961 mu_0)]+1,               (11)
```

where `mu_0` is the lower bound in (9). Substituting THM-2305's two word
floors gives

```text
strict:
  rho_Q>=39002430583/160481782761300,

  mu_0
   =59330091441448204464161088965287
     /1551228875536458713050164713004000000,

  0<|m|<=182078793;

repeated-first, conditional on (5):
  rho_Q>=13560199813/160481782761300,

  mu_0
   =2493436238631076936119894860797
     /1551228875536458713050164713004000000,

  0<|m|<=4332475508.                                (12)
```

Thus

```text
Fourier(B_Q,ms)!=0                                  (13)
```

for a unit-coloured multiplier.

## 5. Lift the character gauge before reducing modulo thirteen

For every integer `r` with `13` not dividing `r`, define

```text
N_r(y)
 =exp(-2*pi*i*r*y/13) M_(r mod 13)(y),

W_r=1_Q N_r.                                        (14)
```

The integer lift in (14) is load-bearing. For
`k,l in {1,...,12}` with `13` not dividing `k+l`, retain `k+l` as the
integer in `{2,...,24}`. Then

```text
W_k W_l conjugate(W_(k+l))
 =1_Q M_k M_l conjugate(M_(k+l mod 13)).            (15)
```

Reducing `k+l` before gauging would insert an erroneous carry factor
`exp(-2*pi*i*y)`.

The exact mask-transport identity is

```text
1_(T^(-1)Q) P^lambda 1_(E_j)
 =P^lambda 1_(E_Q).                                 (16)
```

Changing variables across the thirteen roots therefore gives, for every
integer `h`,

```text
Fourier(W_r,h)
 =13 Fourier(1_(E_Q),13^lambda(r+13h)).             (17)
```

This is the decisive deconvolution: a word-restricted rooted coefficient
is an actual Fourier atom of the **literal source handoff subset** `E_Q`,
not merely of a Perron density or product arrival.

## 6. Three actual atoms and their affine phase

Equations (6), (14), and (15) give the finite character sum

```text
B_Q
 =sum_(k,l nonzero; k+l nonzero mod 13)
    W_k W_l conjugate(W_(k+l)).                     (18)
```

From (13), some character pair has a nonzero coefficient at `d=ms`.
The Fourier-spectrum product inclusion on `Z`, justified by simultaneous
Fejer truncation, supplies integers `h_1,h_2,h_3` with

```text
h_1+h_2-h_3=d                                       (19)
```

and all three corresponding `W` coefficients nonzero. Equation (17)
therefore gives three actual atoms of the same indicator `1_(E_Q)`:

```text
A_1=13^lambda(k+13h_1),

A_2=13^lambda(l+13h_2),

A_3=13^lambda(k+l+13h_3).                           (20)
```

All three have exact thirteen-adic grade `lambda` and nonzero root
characters. Their signed sum is

```text
A_1+A_2-A_3
 =13^(lambda+1)d
 =m c_3,                    gcd(m,91)=1.             (21)
```

The nonzero coefficient

```text
Fourier(1_(E_Q),A_1)
Fourier(1_(E_Q),A_2)
conjugate(Fourier(1_(E_Q),A_3))                     (22)
```

transforms under translation by `tau` by the phase

```text
exp(-2*pi*i*m c_3 tau).                              (23)
```

It is therefore affine-covariant in exactly the deepest unit colour erased
by THM-2312's zero-mode bispectrum. The deepest-comb coefficient at
`-m c_3` supplies the opposite gauge if an invariant scalar certificate is
desired.

## 7. Sharp stopping boundary

The theorem proves

```text
selected shallow owner
  + prescribed expiration clock
  + exact pure/fork blocker word
  + literal source handoff subset
  + three nonzero rooted characters
  + an affine-covariant third-order current
  + gcd(m,91)=1 deepest colour.                     (24)
```

It does **not** prove a pair edge

```text
A-A'=m c_3.
```

There are three nonneutral characters in (21). A one-component handoff
subset can already support a diagonal self-triad while having no relative
component phase-tree edge. In a multi-component expansion, the surviving
term can likewise lie entirely in one component. Thus (22) is not a
THM-2303 same-frequency component-current ratio.

The failed quadratic deconvolution has a sharper hostile control. There
are nonnegative rooted gauges `N` and exact word restrictions `W=1_QN`
with

```text
Fourier(1_Q|N|^2,1)!=0,
```

but every contributing product

```text
Fourier(W,h) conjugate(Fourier(N,h-1))
```

has `Fourier(N,h)=0`. Spectral leakage from the word can therefore make
the restricted endpoint invisible in the unrestricted support. A sufficient
repair is

```text
support Fourier(W) subset support Fourier(N),
```

or the weaker quantitative condition that the leakage contribution is
smaller than the total covariance. Neither follows from the present LRC
data.

An explicit minimal mechanism is seven-periodicity. Put

```text
epsilon=1/1274,

G=1_(union_(r=0)^6 (r/7-epsilon,r/7+epsilon)).
```

Then `G` is inversion-symmetric, `1/7`-periodic, and `T(G)=D_7` with one
occupied root on every active fibre. For every nonzero root character `a`,

```text
Fourier(N_a,h)=13 Fourier(G,a+13h),
```

so the support of `Fourier(N_a,.)` lies in one residue class modulo seven.
Consequently no difference of two `N_a` frequencies has multiplier prime
to seven. Nevertheless, for

```text
Q=D_7 intersection D_1^c,
W=1_Q N_a,
F=W conjugate(N_a),
```

one has `F=1_Q` and

```text
Fourier(F,1)=-sin(pi/49)/pi!=0.
```

Thus the word-restricted cross-correlation has a unit mode although the
unrestricted rooted gauge has no unit-coloured edge anywhere. This is the
first failed implication; the cubic theorem does not use it.

The cubic theorem avoids that false implication by using the exact
transport (16). A subsequent quadratic audit, recorded in
`lrc14-literal-word-quadratic-unit-pair-opus-20260725.md`, does collapse
the restricted object to a same-character unit-coloured pair with much
smaller multiplier bounds. It also proves that this arity reduction is not
the missing ancestry step: the two restricted atoms need not survive in
the unrestricted owner spectrum.

The next decisive target is therefore one of:

```text
force one endpoint of the literal quadratic pair to be a common
  unrestricted-owner/current-service atom;

bound spectral leakage across the complementary words strongly enough to
  preserve one selected endpoint;

or couple the literal pair to THM-2304's middle-depth actual/virtual
  current balance without discarding the affine phase.               (25)
```
