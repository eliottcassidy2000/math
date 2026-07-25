---
id: THM-2300
title: "Small owner multipliers force same-character relation multiples"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. A nonzero
  nonnegative rational step density supported in D_u, with u a
  thirteen-unit, has the following sharp property. For every 1<=m<=6,
  some exact atom at frequency n*m*u is nonzero with n congruent 1 modulo
  13; if the density has at most 2S jumps, one may take
  0<|n|<=13S-1. Applied separately to THM-2276's strict c_1-private source
  and its image support, every pair multiplier m<=6 has bounded exact
  source/ancestry and image atoms at carries of relation multiples that
  remain in the primitive carry's mod-13 character. One multiplier lands
  simultaneously in the original source and its ancestry multiplicity; the
  separately selected image-support multiplier need not agree with it.
  Thus the primitive/base cancellation bank remains 696,
  but the no-forced-same-character-relation-multiple bank has exactly 693
  values. The cutoff is sharp from m=7 onward using an interval whose
  m-fold Perron image is constant. No scalar profile is excluded and
  LRC(14) remains open.
source: codex-2026-07-25-small-owner-same-character-multiple
depends_on:
  - THM-2276-shallow-owner-residue-aligned-crossing
  - THM-2286-endpoint-prony-lift-bank-and-sharp-owner-multiplier-landing
related:
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
  - THM-2278-two-shallow-proper-root-spectrum-and-gap-ancestry-activation
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
script: 04-computation/lrc14_small_owner_same_character_multiple_thm2300.py
output: 05-knowledge/results/lrc14_small_owner_same_character_multiple_thm2300.out
script_sha256: 99ab5281e75a600b3c7181ba521b374de15441363ec8eaf40fdfad9fd83d6d56
output_sha256: 1993e869d37153f147d40562123d97345d0e57e8d40273a4c95f197b9babbc6e
hash_basis: working-tree bytes (LF)
---

# THM-2300 -- small owner multipliers land in the same root character

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2286 proves that an owner multiplier `m<=6` has some nonzero
relation-multiple atom. Its support argument does not control that multiple
modulo thirteen. THM-2278 proves that proper rational root words fire in
every nonzero character. Composing the two mechanisms after normalizing both
the owner and the multiplier supplies the missing congruence:

```text
f supported in D_u
  -> P_m P_u f supported in an arc shorter than twelve root gaps
  -> every nonzero C_13 character fires
  -> the character 1 gives n congruent 1 mod 13
  -> f_hat(nmu)!=0.                                  (1)
```

The factor `n` is a thirteen-unit and preserves the primitive pair carry's
root character. It is not forced to equal one, so the primitive/base atom
can still vanish.

## 1. Weighted proper-root lemma

Let

```text
zeta=exp(2*pi*i/13)
```

and let

```text
(a_0,...,a_12) in Q_(>=0)^13                       (2)
```

be nonzero but have at least one zero entry. Then, for every
`k=1,...,12`,

```text
sum_(r=0)^12 a_r zeta^(kr)!=0.                      (3)
```

Indeed, if the rational polynomial

```text
A(X)=sum_(r=0)^12 a_r X^r                           (4)
```

vanished at `zeta^k`, the minimal polynomial

```text
Phi_13(X)=1+X+...+X^12                              (5)
```

would divide `A`. Since both degrees are at most twelve,

```text
A=c Phi_13
```

for some rational `c`, with the zero polynomial included when `c=0`.
All thirteen coefficients would be equal. That contradicts (2).

The rational-weight hypothesis is load-bearing. Arbitrary positive real
weights can solve two real linear cancellation equations.

## 2. Normalizing the owner and multiplier

Let `f` be a nonzero nonnegative rational-valued step function and assume

```text
support(f) subset D_u,
13 does not divide u.                               (6)
```

For positive integers `a`, use the normalized Perron operator

```text
P_a f(y)=(1/a)sum_(r=0)^(a-1)f((y+r)/a).            (7)
```

It preserves mass and obeys

```text
(P_a f)_hat(q)=f_hat(aq).                           (8)
```

Put

```text
h=P_m P_u f.                                        (9)
```

The first push gives support in `D_1`. When `1<=m<=6`, the second gives

```text
support(h) subset mD_1
 =(-m/14,m/14) mod 1,                               (10)

measure(mD_1)=m/7<=6/7<12/13.                       (11)
```

The shortest circular arc containing all thirteen equally spaced roots has
length `12/13`. Therefore every occupied thirteen-root fibre of `h` has at
least one zero entry. The values of `h` are rational because (7) averages
rational values. The weighted proper-root lemma applies pointwise.

For

```text
M_k(y)=sum_(r=0)^12 h((y+r)/13)zeta^(-kr),          (12)
```

every nonzero character is nonzero on `T(support(h))`, a set of positive
measure. The root/Fourier identity gives

```text
integral |M_k(y)|^2dy
 =169 sum_(q congruent k mod 13)|h_hat(q)|^2
 >0,                                    1<=k<=12.  (13)
```

In particular the residue class `k=1` contains a nonzero exact atom.

## 3. Endpoint-Prony in the character one progression

Suppose `f` has at most `2S` jumps. Perron transport maps jumps to their
images and cannot create jump locations, so `h` also has at most `2S`
jumps. For `q=1+13ell`, distributional differentiation writes

```text
2*pi*i(1+13ell)h_hat(1+13ell)
 =sum_(x in Jump(h))
    Delta_h(x)exp(-2*pi*i*x)
               [exp(-2*pi*i*13x)]^ell.              (14)
```

After collisions, this is a nonzero exponential sum with at most `2S`
nodes; it is nonzero by (13) with `k=1`. It cannot vanish at the `2S`
consecutive indices

```text
ell=-S,-S+1,...,S-1.                               (15)
```

For one of them,

```text
n=1+13ell congruent 1 mod 13,

0<|n|<=13S-1,

h_hat(n)=f_hat(nmu)!=0.                             (16)
```

This proves the general small-multiplier statement. The exact residue
condition, not merely `13` not dividing `n`, is what preserves a selected
root character after multiplying a relation.

## 4. Application to the THM-2276 shallow pair

On each of the `120` interior first-depth-one scalar rows, THM-2276 supplies
a primitive pair relation `p`, a thirteen-unit `m`, and

```text
K=A(p)=plus or minus m c_1,
kappa=K/13=plus or minus m u_1,
1<=m<=757.                                          (17)
```

Let `E_1` be its full strict `c_1`-private locus and put

```text
F=T(E_1),
f_F=1_F,
G=P1_(E_1),
g=13G
 =sum_(r=0)^12 1_(E_1)((y+r)/13).                   (18)
```

Both `G` and `1_F` are nonzero rational step functions supported in
`D_(u_1)`. If

```text
S=H+sum_i q_i+sum_j c_j,                            (19)
```

each has at most `2S` jumps: `E_1` has boundaries only in the nine scalar
boundary banks, Perron transport does not create jump locations, and the
image support has boundary contained in the image of the source boundary.

Assume `m<=6`. Apply (16) first to `(f,u)=(1_F,u_1)`.
There is an integer `n_F` such that

```text
n_F congruent 1 mod 13,
0<|n_F|<=13S-1,

(1_F)_hat(n_F kappa)!=0.                            (21)
```

The sign in (17) is handled by conjugation.

Apply (16) separately to `(f,u)=(G,u_1)`.
There is, generally different, `n_G` such that

```text
n_G congruent 1 mod 13,
0<|n_G|<=13S-1,

G_hat(n_G kappa)!=0.                                (22)
```

The Perron Fourier law and `K=13kappa` give, for this **same** `n_G`,

```text
g_hat(n_G kappa)=13G_hat(n_G kappa)!=0,

(1_(E_1))_hat(n_G K)
 =G_hat(n_G kappa)
 !=0.                                               (23)
```

Thus the relation multiples `n_F p` and `n_G p` land exactly in the image
and source/ancestry spectra, respectively. Their carries retain valuation
one and the primitive root character:

```text
n_F kappa congruent kappa mod 13,
n_G kappa congruent kappa mod 13.                   (24)
```

Their heights are at most

```text
(13S-1)9841.                                        (25)
```

Equations (21) and (22) do **not** assert `n_F=n_G`. The source and its
ancestry multiplicity use the same `n_G`, but no simultaneous
source/image-support atom has been proved.

## 5. Exact bank refinement

Among `1,...,757`, exactly

```text
699
```

integers are thirteen-units. THM-2276 forces the primitive/base atom for

```text
m=1,2,3
```

and has exact local primitive/base cancellations at the remaining

```text
696
```

values.

The present theorem does not change that primitive count. It proves a
different and useful refinement:

```text
m=1,...,6:
 a bounded same-character relation multiple is forced;

7<=m<=757, 13 does not divide m:
 693 values remain in the no-forced-same-character-multiple bank.        (26)
```

Thus the three multipliers `4,5,6` are semi-hard: their primitive atom can
cancel, but not every same-character relation multiple can cancel.

## 6. The cutoff at seven is sharp

For every `m>=7`, take

```text
f_m=1_[-1/(2m),1/(2m)).                            (27)
```

Up to its null endpoints,

```text
support(f_m) subset D_1.                            (28)
```

Multiplication by `m` maps the interval in (27) bijectively almost
everywhere onto the circle. Therefore

```text
P_m f_m=1/m                                         (29)
```

almost everywhere, and

```text
(f_m)_hat(nm)=0                    for every n!=0.  (30)
```

This kills every nonzero relation multiple, including all factors
`n congruent 1 mod 13`. At `m=7`, (27) is exactly the owner interval
`D_1` up to endpoints, so the `6/7` boundary in (11) cannot be improved.

The hostile interval is an owner-support control. It is not a strict
private locus or a global scalar cover.

## 7. Connection and loss ledger

```text
source:
  THM-2276's exact shallow pair and positive private source/image sets;

map:
  normalize first by the owner u and then by the absolute multiplier m,
  retain the rational weighted root word, select character 1, and use
  endpoint-Prony on that exact progression;

preserved:
  exact relation multiple, owner label, absolute multiplier, primitive
  mod-13 character, valuation one, source ancestry multiplicity, and an
  explicit finite bound;

new output:
  one bounded simultaneous source/ancestry atom and a separately bounded
  image-support atom for m<=6;

not preserved:
  primitivity, one common multiplier for source and image, current blocker
  service, coefficient sign in a cover convolution, or a uniform lower
  bound for the selected exact atom;

sharp boundary:
  m=7, where one owner interval makes P_m constant and kills all multiples.
                                                                    (31)
```

The theorem removes no scalar row. A later closure still has to pair one
landed relation multiple with the same current-service or signed-cover term.

## 8. Exact verification

The companion uses exact integer and `Fraction` arithmetic. It checks the
weighted-root degree ledger, the sharp `6/7<12/13` threshold, all endpoint
windows, the `699=3+3+693` census, the relation-multiple heights and residues,
and every hostile interval from `m=7` through `757`. Reproduce with

```bash
python3 04-computation/lrc14_small_owner_same_character_multiple_thm2300.py
python3 -O 04-computation/lrc14_small_owner_same_character_multiple_thm2300.py
```

Every load-bearing check raises explicitly, so optimized mode executes the
same audit. QED.
