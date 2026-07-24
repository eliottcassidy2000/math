---
id: THM-2145
title: "Two-block spectral crossing and the exact low-core relation carry"
status: >
  PROVED from THM-1234, THM-1221, a positive Jackson smoothing, and a finite
  exact 177-cell core sweep. Every zero-measure-safe 6+7 split has a genuinely
  crossing relation of coefficient height 584. If the seven retained speeds
  form a subset of {1,...,13}, height 298 suffices and its retained carry is
  at most 20860. A support-two output can still be a tautological reduced
  cross-pair relation, so this does not by itself close defect six or LRC(14).
source: codex-2026-07-24-relation-carry-spectrum
depends_on:
  - THM-1221
  - THM-1234
related:
  - THM-2051
  - THM-2052
  - THM-2054
  - THM-2144
script: 04-computation/lrc14_two_block_crossing_referee_codex_20260724.py
output: 05-knowledge/results/lrc14_two_block_crossing_referee_codex_20260724.out
script_sha256: baea8a3d8cd61a5bd987119bc71d81e993f65a2e374b8b65df0b51d584fd2e60
output_sha256: fecf29eb6ebd86fcf6534ea070ec50ad0841faef41bce84d80b6b7f57b458aa3
hash_basis: working-tree bytes (LF)
---

# THM-2145 -- two-block spectral crossing and low-core carry

At radius `1/14`, put

```text
G={x in R/Z: ||x||_(R/Z)>=1/14},                       (1)
G_A={t: at in G for every a in A}.                     (2)
```

Sets and rows below are labelled when coefficients are extracted. Repeated
speeds are allowed in the abstract lemma; the LRC specializations use
distinct positive integer speeds.

## 1. Abstract two-block crossing lemma

Let `A=(a_1,...,a_r)` and `B=(b_1,...,b_s)` be finite positive-integer rows.
Suppose `q` is a real trigonometric polynomial of degree at most `H` such
that

```text
0<=q<=1,              ||q-1_G||_1<=eta.                (3)
```

Assume

```text
mu(G_A)>=alpha,       mu(G_B)>=beta,
mu(G_A intersect G_B)=0,                               (4)
alpha-r eta>0,        beta-s eta>0,
(alpha-r eta)(beta-s eta)>(r+s)eta.                    (5)
```

Then there are coefficient vectors

```text
u in [-H,H]^r intersect Z^r,
w in [-H,H]^s intersect Z^s                          (6)
```

such that

```text
sum_i u_i a_i + sum_j w_j b_j=0,
sum_i u_i a_i !=0,
sum_j w_j b_j !=0.                                    (7)
```

Thus the relation genuinely crosses the labelled cut: both restrictions are
nonzero.

### Proof

Define

```text
Q_A(t)=product_i q(a_i t),       Q_B(t)=product_j q(b_j t). (8)
```

Integer dilation preserves Haar measure. The product telescope for
`[0,1]`-valued factors therefore gives

```text
integral Q_A>=alpha-r eta,
integral Q_B>=beta-s eta,                              (9)
integral Q_A Q_B<=mu(G_A intersect G_B)+(r+s)eta
                =(r+s)eta.                            (10)
```

By (5),

```text
(integral Q_A)(integral Q_B)>integral Q_A Q_B.         (11)
```

If the nonzero Fourier supports of `Q_A` and `Q_B` were disjoint after
reflection, finite Fourier orthogonality would instead give equality in
(11). Hence some nonzero integer frequency `nu` satisfies

```text
Fourier(Q_A,nu)!=0,       Fourier(Q_B,-nu)!=0.          (12)
```

Expand the two nonzero coefficients as finite convolutions. At least one
nonzero summand on each side supplies indices

```text
sum_i u_i a_i=nu,       sum_j w_j b_j=-nu,
|u_i|,|w_j|<=H.                                           (13)
```

Equations (12)--(13) prove (7). QED.

This is a common-frequency theorem, not an independence estimate. It loses
the number of convolution representations and any linear independence among
the relations it extracts.

## 2. The Jackson polynomial and its exact error certificate

For `N>=2`, use the normalized Fejer kernel

```text
F_N(t)=(1/N)(sin(pi N t)/sin(pi t))^2                  (14)
```

with its continuous values at integers, and set

```text
J_N=F_N^2/integral F_N^2,       q_N=J_N*1_G.           (15)
```

Then `J_N>=0`, `integral J_N=1`, and

```text
0<=q_N<=1,         degree(q_N)<=H=2N-2.                (16)
```

Translation of the interval `G` by `x` changes its indicator in `L^1` by
`2 min(||x||,1/7)`. Consequently

```text
||q_N-1_G||_1
 =integral 2 min(||x||,1/7) J_N(x) dx
 <=2 integral ||x|| J_N(x) dx=:eta_N.                  (17)
```

The Fourier coefficients can be kept integral. Put

```text
C_0=N(2N^2+1)/3,                                       (18)

C_k=(4N^3-6Nk^2+2N+3k^3-3k)/6,       0<=k<=N,
C_k=((2N-k)^3-(2N-k))/6,              N<k<=2N-2.       (19)
```

Then `Fourier(J_N,k)=C_|k|/C_0`, all these coefficients are positive, and the
Fourier series of circle distance gives

```text
eta_N
 =1/2-[4/(pi^2 C_0)]
       sum_(1<=k<=2N-3, k odd) C_k/k^2.                (20)
```

The classical rational inequality

```text
pi<355/113                                               (21)
```

turns (20) into a strict rational upper bound.

For `k!=0`,

```text
Fourier(1_G,k)=-sin(pi k/7)/(pi k).                    (22)
```

Because the Jackson coefficient is positive throughout its support,

```text
Fourier(q_N,k)!=0
iff
0<|k|<=H and 7 does not divide k.                       (23)
```

The crossing relation from Section 1 may therefore be chosen with every
coefficient either zero or not divisible by `7`.

## 3. Universal 6+7 crossing at height 584

Let `A` and `B` be disjoint labelled blocks of six and seven distinct
positive integer speeds, and suppose

```text
mu(G_A intersect G_B)=0.                               (24)
```

THM-1234, equation (27), gives

```text
mu(G_A)>=1-212/273=61/273.                             (25)
```

THM-1221 gives

```text
mu(G_B)>=15/154.                                       (26)
```

For `N=293`, hence `H=584`, the exact odd-mode sum in (20), together with
(21), proves

```text
eta_N<1439/1000000.                                    (27)
```

At the right side of (27), the crossing margin is

```text
(61/273-6eta)(15/154-7eta)-13eta
 =548679648961/10510500000000000
 >0.                                                   (28)
```

The two factors are positive. The abstract lemma proves:

> Every zero-measure-safe labelled `6+7` split of thirteen distinct positive
> integer speeds has a crossing integer relation of coefficient height at
> most `584`; every nonzero coefficient can be chosen prime to `7`.

The exact referee checks that `N=292` does not close (28) with this same
Jackson moment and rational-`pi` ledger. Thus `N=293` is certificate-minimal
for the stated ledger, not optimal among all positive kernels.

## 4. Exact seven-core floor inside `{1,...,13}`

For retained cores one can replace the universal floor (26) by a much larger
finite exact value.

> For every seven-set `E subset {1,...,13}`,
>
> ```text
> mu(G_E)>=45107/229320.                                (29)
> ```
>
> Equality occurs uniquely at
>
> ```text
> E=(1,5,7,8,9,11,13).                                 (30)
> ```

Here is the exhaustive certificate. Form the boundary set

```text
B={0,1}
  union {(14j +/- 1)/(14v) mod 1:
         1<=v<=13 and 0<=j<v}.                         (31)
```

Exact fraction deduplication gives `178` points and `177` open cells. On
each cell the thirteen safe/danger indicators are constant. Summing cell
lengths for all `binomial(13,7)=1716` seven-sets proves (29)--(30).

An independent interval-union merge for (30) gives

```text
mu(union_(e in E) D_e)=4052686/5045040
                      =184213/229320,
mu(G_E)=992354/5045040=45107/229320.                   (32)
```

Boundaries have measure zero, so the midpoint-cell and interval-union
computations use the same event.

## 5. Defect-six specialization: height 298 and carry 20860

Let `E` be any seven-subset of `{1,...,13}` and let `F` be a disjoint block
of six distinct positive speeds. Suppose the full thirteen-speed safe set has
measure zero. Use

```text
alpha=61/273,       beta=45107/229320.                  (33)
```

For `N=150`, `H=298`, equations (18)--(21) give

```text
C_0=2250050,
sum_(1<=k<=297, k odd) C_k/k^2>2760290,                (34)

eta_N
 <159340717/56712510250
 <281/100000.                                          (35)
```

At `eta=281/100000`,

```text
(61/273-6eta)(45107/229320-7eta)-13eta
 =322476924229/7825545000000000
 >0.                                                   (36)
```

Therefore there are coefficients `a_f,b_e`, each zero or a nonzero
`7`-unit, such that

```text
sum_(f in F) a_f f + sum_(e in E)b_e e=0,
0<|a_f|,|b_e|<=298 when nonzero,                        (37)
sum_(f in F)a_f f=-sum_(e in E)b_e e!=0.               (38)
```

Since the seven largest elements of `{1,...,13}` sum to `70`, the retained
partial sum, or **relation carry**, satisfies

```text
0<|sum_(e in E)b_e e|<=298*70=20860.                   (39)
```

The exact referee checks that `N=149` does not close this same
Jackson/rational-`pi` ledger. Hence `298` is certificate-minimal within this
specialized ledger.

## 6. The support-two branch and the honest frontier

Every relation in (7) uses at least one coefficient on each side. If its
total support is two, for some `f in F`, `e in E`, and `g=gcd(f,e)`, its
primitive form is

```text
(e/g)f-(f/g)e=0.                                      (40)
```

Thus the theorem gives the exact dichotomy

```text
some cross pair has e/g<=H and f/g<=H,
with both reduced coefficients prime to 7,
or
there is a support-at-least-three crossing relation.   (41)
```

This caveat is load-bearing. Even height `298`, though below the raw
first-replacement ceiling `3003/8`, does not automatically exclude (40):
a common gcd can make `f/g` small. The result therefore does **not** by itself

1. force a new independent relation;
2. finite-bound the six replacement speeds;
3. produce a complete defect-six dynamic program; or
4. prove LRC(14).

It does provide two rigid sidecars that the earlier crossing proposal lacked:
all nonzero coefficients are `7`-units, and the normalized retained carry is
in the explicit interval `[-20860,20860] minus {0}`. A successful next step
must either exclude the support-two commensurability branch, harvest relation
rank, or add a magnitude/termination sidecar to the carry state.

The companion script uses exact `Fraction` arithmetic, runs the 177-cell
sweep and an independent interval-union check, evaluates the Jackson sums,
and scans the adjacent `N` boundaries. Normal and optimized Python execute
the same explicit raising checks and produce identical output.
