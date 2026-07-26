---
id: THM-2410
title: "Full-coordinate projector, local Gram, and integrated-phase boundary"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING. The lawful
  independent C_13 translation bank of all six unit/guard safe factors,
  all three blocker safe factors, and one deep danger factor descends
  losslessly through the common-root gauge. Every all-nonzero character
  in the annihilator has a pointwise nonzero factorization and a uniform
  positive local Gram floor. There are exactly 396,907,776 such endpoint
  characters and 4,365,985,536 eligible endpoint/deep pairs. A finite
  exact translated-pattern hostile shows that this quadratic local
  positivity does not force the integrated linear coefficient to
  survive; the hostile is deliberately not a physical common-root
  trajectory. Conditional on one physical coefficient surviving, the
  existing m-then-X expansion gives an all-coordinate mod-13 relation
  residue and a 91-unit deep multiplier. The mod-seven word action
  remains coupled by convolution, so this is not an all-91-unit address
  theorem, a row exclusion, or a proof of LRC(14).
source: codex-2026-07-26-full-coordinate-projector
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2403-clean-toothpick-unequal-slope-target-axis-imbalance
related:
  - THM-2325-prescribed-target-gain-full-lattice-91-unit-needle-bank
  - THM-2397-clean-root-same-parent-charged-role-partition
  - THM-2400-clean-parent-root-gauge-quotient-and-target-slope-boundary
  - THM-2407-owner-or-source-deletion-target-current-dichotomy
  - THM-2409-unfiltered-septimal-source-completion-and-word-phase-boundary
script: 04-computation/lrc14_full_coordinate_projector_thm2410.py
output: 05-knowledge/results/lrc14_full_coordinate_projector_thm2410.out
script_sha256: e2e9ff1261c6aae2a497706d93f51d73a300db0bed3cf8d8ea760eae4e5936bd
output_sha256: cf38884ce00dc6c1b27180a3aaab7919feb13c3dd1553a4fe4bacb252ee13675
hash_basis: working-tree bytes (LF)
---

# THM-2410 -- the full-coordinate local Gram does not fix integrated phase

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2403 proves that one lawful unequal-slope target line fires every
nonzero first-target colour. The natural stronger experiment is to
translate all nine present factors independently and project directly
onto relation characters with every coordinate nonzero. That experiment
has an unexpectedly sharp answer:

```text
full independent C_13 projector
  -> lossless common-root quotient
  -> 396,907,776 all-nonzero endpoint characters
  -> pointwise nonzero factorization
  -> uniform positive local Gram

but

local Gram != integrated linear phase coherence.                 (1)
```

The first chain is proved. The failed implication in (1) is witnessed
by a finite exact translated-pattern hostile. The hostile is not passed
off as a physical LRC trajectory; it isolates precisely the sidecar
still missing from the physical packet.

## 1. The lawful ten-variable shift table

Put

```text
p=13,                    G=F_p.
```

Use the nine present labels in the order

```text
one guard, five ordinary q roles, three blocker roles.           (2)
```

Let their integer speeds be

```text
w=(w_1,...,w_9).
```

Modulo thirteen, exactly the first six coordinates are units and the
last three are zero. Let

```text
I_i:T -> {0,1}
```

be the corresponding present safe indicators. Thus the danger mask of
the guard has three or four consecutive roots on a `C_13` fibre, and
each other danger mask has one or two consecutive roots. Let

```text
d:T -> {0,1}
```

be the centred ordinary danger indicator and let `c=c_3`, so
`13|c`. Fix a positive-mass transported Boolean terminal word `Q(Rx)`
with

```text
R=13^h,                  h>=1.                                  (3)
```

Thus `mu(Q)>0`.

For

```text
r in G,                  s=(s_1,...,s_9) in G^9,
```

define the full coordinate-shift table

```text
H(r,s)
 =integral_T Q(Rx)d(cx-r/p)
    product_(i=1)^9 I_i(w_i x-s_i/p) dx.                        (4)
```

Every entry is a lawful coefficient of the globally shifted present
packet. No clean-set indicator, owner selector, or derived least-role
factor has been inserted.

## 2. The common-root gauge and the lossless quotient normalization

Write `wbar` for the reduction of `w` in `G^9`. Translation of the
integration variable by `t/p` gives

```text
H(r,s+t wbar)=H(r,s)                  for every t in G.          (5)
```

Indeed, the coordinate translations cancel exactly. The deep factor is
unchanged because `13|c`, and the word is unchanged because `13|R`.

Let

```text
zeta=exp(2*pi*i/p).
```

For `m in G` and `k in G^9`, define the normalized full transform

```text
B(m,k)
 =p^(-10) sum_(r in G,s in G^9)
     H(r,s) zeta^(m r+k.s).                                    (6)
```

Equation (5) gives

```text
B(m,k)=0 unless k.wbar=0.                                      (7)
```

When `k.wbar=0`, the summand is constant on the `p`-point gauge
orbits. Hence the full normalization is exactly the quotient
normalization:

```text
B(m,k)
 =p^(-9) sum_(r in G,[s] in G^9/<wbar>)
     H(r,[s]) zeta^(m r+k.s).                                  (8)
```

There is no hidden factor of thirteen in (8).

## 3. Pointwise factorization and the local Gram floor

Put

```text
A_(i,k)(x)
 =p^(-1)sum_(a in G)
    I_i(w_i x-a/p) zeta^(k a),

D_m(x)
 =p^(-1)sum_(a in G)
    d(cx-a/p) zeta^(m a).                                      (9)
```

Finite rearrangement in (6) gives the exact pointwise factorization

```text
B(m,k)
 =integral_T Q(Rx)D_m(x)
    product_(i=1)^9 A_(i,k_i)(x) dx.                           (10)
```

Fix `x` away from the strict interval endpoints. Every danger-root
set in (9) is a nonempty proper consecutive subset of `G`. For
`k!=0`, the safe coefficient is the negative of the danger
coefficient and has the form

```text
A_(i,k)(x)
 =-p^(-1) zeta^(k a)
    (1-zeta^(k n))/(1-zeta^k),                                 (11)
```

where `1<=n<=4`. The same formula without the minus sign holds for
`D_m`, with `1<=n<=2`. Since `p` is prime, these coefficients never
vanish when their characters are nonzero.

More quantitatively,

```text
|A_(i,k)(x)|, |D_m(x)|
 >=sin(pi/13)/13
 >2/169.                                                       (12)
```

For every

```text
m!=0,            all k_i!=0,            k.wbar=0,               (13)
```

put

```text
Phi_(m,k)(x)=D_m(x)product_i A_(i,k_i)(x).
```

The exact local Gram therefore obeys

```text
integral_T Q(Rx)|Phi_(m,k)(x)|^2 dx
 >=mu(Q)(sin(pi/13)/13)^20
 >mu(Q)(2/169)^20.                                             (14)
```

Here `Q` is Boolean, so `Q^2=Q`; Haar invariance gives
`mu(Q(R .))=mu(Q)`. Equation (14) says that every eligible local
channel has positive quadratic mass. It does **not** say that the
linear integral (10) is nonzero.

Without the deep factor, the corresponding endpoint-only exponent is
`18`.

## 4. The exact all-nonzero character census

Multiplying each of the six unit coordinates by its nonzero speed
reduces the annihilator equation to

```text
x_1+...+x_6=0,               x_i in G^*.                       (15)
```

The number of solutions is

```text
((p-1)^6+(p-1))/p
 =(12^6+12)/13
 =229,692.                                                       (16)
```

The three blocker characters are arbitrary nonzero elements of `G`.
Thus the number of all-coordinate-nonzero endpoint characters is

```text
229,692*12^3
 =396,907,776.                                                   (17)
```

For each such `k`, choose a nonzero deep colour `m` and form

```text
q=k+m e_c.                                                       (18)
```

Exactly one of the twelve nonzero choices, `m=-k_c`, makes the
`c` coordinate of `q` vanish. The number of pairs for which every
coordinate of `q` remains nonzero is therefore

```text
396,907,776*11
 =4,365,985,536.                                                 (19)
```

These are relation **characters** and eligible deep pairs. Equations
(17)--(19) do not assert that the corresponding integrated currents
are nonzero.

## 5. An exact finite integrated-phase hostile

The missing implication cannot be recovered from local nonvanishing or
from the Gram (14) alone.

Use the unit-speed residues

```text
wbar=(1,2,3,4,5,11,0,0,0)                                    (20)
```

and the all-one character

```text
k=(1,1,1,1,1,1,1,1,1).
```

Then `k.wbar=26=0` in `G`. On the first six coordinates choose danger
masks

```text
{9,10,11,12},
{0,1}, {2,3}, {3,4}, {5,6}, {7,8}.                             (21)
```

Every displayed mask is cyclic-consecutive. They cover `G`, have total
incidence `14`, and have their unique double root at `3`. On the last
three coordinates choose, for example,

```text
{1,2}, {3,4}, {5,6}.                                           (22)
```

Let `T(s)` be the rank-one Boolean safe tensor obtained by taking the
product of the nine safe-mask indicators. Put

```text
v=(1,1,1,1,1,1,0,0,0),

T_t(s)=T(s-t v),                    t in G.                      (23)
```

Retain the fixed cyclic-consecutive deep danger mask

```text
D(r)=1_(r in {0,1}),                m=1,

U_t(r,s)=D(r)T_t(s).                                             (23a)
```

Every one-coordinate nontrivial Fourier factor is nonzero by
prime-cyclotomic irreducibility, including `D_hat(1)`. Hence

```text
U_hat_t(1,k)
 =D_hat(1)zeta^(t k.v)T_hat_0(k)
 !=0
```

for every `t`. The unit masks in (21) cover `G`, so also

```text
T_t(0)=0                         for every t.                    (24)
```

Nevertheless

```text
k.v=6!=0,

(1/13)sum_t U_hat_t(1,k)
 =D_hat(1)T_hat_0(k)(1/13)sum_t zeta^(6t)
 =0.                                                            (25)
```

Thus an anchored ten-factor family can have a nonzero eligible
deep-inclusive local coefficient at every parameter value and still
have zero integrated coefficient.

The scope boundary is essential. The trajectory `v` in (23) is **not**
the physical common-root trajectory `wbar`. On the physical trajectory,
`k.wbar=0` makes the phase constant, exactly as THM-2400 predicts.
Consequently (25) is a universal inference hostile, not a counterexample
to the physical coefficient `B(m,k)` in (10). The physical integrated
nonvanishing problem remains open.

## 6. Conditional exact relation-residue consequence

Suppose an eligible physical coefficient in (13) actually satisfies

```text
B(m,k)!=0,             m!=-k_c.                                (26)
```

Embed (4) in the THM-2334/THM-2365 Poisson--Abel current with the same
shift convention. Its `m`-then-`X` expansion has exact target residue

```text
q=k+m e_c.                                                       (27)
```

Type the residue fibre explicitly. If `C_rho(a;X,n)` is THM-2334's
Poisson-smoothed exact-address coefficient, put

```text
C_13(q;X,n)
 :=lim_(rho->1-)
    sum_(a in Lambda; a=q mod 13) C_rho(a;X,n).                  (27a)
```

This is the **mod-thirteen Abel residue fibre** of THM-2334. It is not
being reinterpreted as an ordinarily absolutely convergent sum over the
infinite residue class.

With this typing, THM-2365's two convergence statements give the
absolutely iterated expansion

```text
B(m,k)
 =sum_(n=m mod 13)
    [sum_X C_13(q;X,n)].                                        (27b)
```

For each fixed `n`, the inner `X`-sum is absolutely convergent by
endpoint Parseval/Cauchy--Schwarz. After that inner sum is taken, the
outer `n`-sum is absolutely convergent by the deep-mode estimate:

```text
sum_(n=m mod 13)
 |sum_X C_13(q;X,n)|<infinity.                                  (27c)
```

No jointly absolute double sum and no ordinary absolute residue-class
sum is claimed. Consequently (26) supplies some exact deep multiplier
`n` and then some ordinary frequency `X` with

```text
n=m mod 13,

C_13(q;X,n)!=0.                                                 (28)
```

By (18), every coordinate of `q` is a unit modulo thirteen. Since
`m!=0`, also `13` does not divide `n`. The centred deep-danger Fourier
coefficient vanishes at every nonzero multiple of seven, so a live term
in (28) has

```text
7 does not divide n,

gcd(n,91)=1.                                                     (29)
```

This is the strongest unconditional implication of one surviving
physical full-coordinate coefficient:

```text
one nonzero B(m,k)
  -> one all-coordinate-unit mod-13 relation residue q
  -> one exact X and a 91-unit deep multiplier n.                (30)
```

It does not say that every coordinate of the exact relation address is
a unit modulo seven.

## 7. Why the projector does not split modulo seven

At modulus thirteen the word in (4) is target-neutral because

```text
R s_i/13 in Z.                                                   (31)
```

At modulus `91`, the same coordinate translation acts on the word by

```text
R ell_i/91=13^(h-1)ell_i/7.                                    (32)
```

Thus its septimal part is active. The present factor and transported
word are moved by the same coordinate parameter; their character
sequences combine by convolution. They cannot be projected
independently and multiplied after the fact.

The smallest exact stopping witness already occurs on `F_7`. Let

```text
f=1_{0},                  g=1-1_{0}.                             (33)
```

Both `f` and `g` have nonzero Fourier coefficients at every nontrivial
septimal character, but

```text
f g=0.                                                          (34)
```

Equivalently, the convolution of their two full nontrivial spectra
cancels identically. Therefore separate mod-seven support of the
present packet and the word does not produce an all-unit mod-seven
current.

The missing sidecar must control the **coupled** septimal convolution,
for example by a common phase cone, a polarized parent-labelled Gram,
or a genuinely joint word/present projector. Neither (14) nor the CRT
count in (19) supplies it.

Candidate THM-2409 reaches the same boundary from the complementary
source-deletion direction: with `Q=1` its seven translated-source
coefficients are all-or-flat, while a lawfully shifted delayed word can
destroy the partition identity completely. The two candidates therefore
agree on the next object--a joint present/word septimal phase gauge--from
opposite coordinate systems.

## 8. Exact scope

What is proved:

- the ten-variable table (4) is a lawful independent translation bank;
- the common-root gauge quotient and normalization (8) are lossless;
- every eligible local factor is nonzero and has the Gram floor (14);
- the exact endpoint and endpoint/deep censuses are (17) and (19);
- local nonvanishing and positive Gram alone do not imply integrated
  phase survival;
- any surviving physical coefficient has the conditional exact
  consequence (30); and
- mod-seven word coupling blocks a formal CRT upgrade.

What is not proved:

- that any of the `396,907,776` physical endpoint coefficients is
  nonzero after integration;
- that the finite hostile (25) is a physical scalar-row trajectory;
- that an exact relation address is all-unit modulo `91`;
- that a positive shallow owner, a preselected marked triangle, or a
  canonical terminal phase is retained; or
- that any scalar row is excluded.

The scalar ledger remains `165`, and LRC(14) remains open.

## 9. Exact companion

The dependency-free exact companion:

- verifies prime-cyclotomic nonvanishing for every one of the `8,190`
  nonempty proper masks of `F_13` and all twelve nonzero characters;
- checks the annihilator count (16), the blocker factor, and the
  `11`-colour deep refinement;
- verifies the full-to-quotient normalization;
- reconstructs the cover/incidence hostile (20)--(25) exactly in
  `Q[zeta_13]`, including its fixed nonzero deep mask/colour;
- records the rational Gram floors with exponents `18` and `20`; and
- checks the septimal convolution hostile (33)--(34).

Run

```bash
python3 04-computation/lrc14_full_coordinate_projector_thm2410.py
python3 -O 04-computation/lrc14_full_coordinate_projector_thm2410.py
```

Both transcripts must byte-match, after LF normalization,

```text
05-knowledge/results/lrc14_full_coordinate_projector_thm2410.out
```

Every truth-bearing finite check raises explicitly, so optimized mode
executes the same audit.
