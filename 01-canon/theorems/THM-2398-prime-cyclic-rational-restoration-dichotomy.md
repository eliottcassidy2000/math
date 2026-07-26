---
id: THM-2398
title: "Prime-cyclic rational restoration dichotomy"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING. A
  nonzero nonnegative rational translation-covariant operator on C_p,
  p prime, is either a scalar multiple of uniform rank-one averaging or is invertible
  on every character. A common nonflat kernel applied to two currents
  preserves their exact charged-product phase, with a cyclotomic norm
  floor. A rational common-fibre cross-correlation with zero base
  overlap and positive two-sided mass fires every nonzero character.
  Applied to a fixed singleton/adjacent mask and one disjoint two-root
  role pattern at p=13, all twelve colours survive on one cell with
  explicit uniform floors. The theorem does not construct the common
  terminal intertwiner, a terminal word, or a proof of LRC(14).
source: codex-2026-07-26-prime-cyclic-restoration
depends_on:
  - THM-2392-clean-toothpick-or-bounded-cross-ancestor-cage
  - THM-2396-common-core-forty-nine-orbit-word-incompatibility
related:
  - THM-2300-small-owner-multipliers-force-same-character-relation-multiples
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2368-owner-pivot-root-fibre-radon-invertibility
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2397-clean-root-same-parent-charged-role-partition
script: 04-computation/lrc14_prime_cyclic_rational_restoration_thm2398.py
output: 05-knowledge/results/lrc14_prime_cyclic_rational_restoration_thm2398.out
script_sha256: e2057ec13d6468685de21b0af9dec540a1dea0412523359b7e7fe73a592eab90
output_sha256: 919ef7aeb97650fd08cf360ddee8c1c4c559d453144980d13d42a03ef433f3c4
hash_basis: working-tree bytes (LF)
---

# THM-2398 -- prime-cyclic rational restoration dichotomy

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

The terminal-transport obstruction after a charged root carrier can be
phrased as a question about a translation-covariant restoration operator.
Positivity, mass preservation, reversibility, and even strict
positivity-improvement do not suffice: uniform averaging has all those
properties and kills every charged character.

At a prime cyclic fibre, rationality makes that hostile unique. A rational
kernel which kills one nonzero character is forced by cyclotomic
irreducibility to be exactly uniform. Consequently:

```text
common rational target-covariant restoration
  + one unequal transition weight
  -> every charged target survives;

the same oriented restoration on two currents
  -> their charged-product phase is unchanged.                 (1)
```

The theorem also gives a different, often cheaper test. Two disjoint
rational packets in one common root fibre have a nonnegative
cross-correlation whose zero shift vanishes. If the packets coexist on
positive base mass, that correlation is nonzero and hence nonflat, so all
nonzero characters survive.

## 1. Prime-cyclic all-or-flat theorem

Let `p` be an odd prime, let

```text
G=F_p,

zeta=exp(2*pi*i/p),

f_hat(k)=(1/p)sum_(x in G)f(x)zeta^(-kx).            (2)
```

For a nonzero kernel

```text
kappa=(kappa_s)_(s in G),        kappa_s in Q_(>=0), (3)
```

define the translation-covariant operator

```text
(Kf)(x)=sum_(s in G)kappa_s f(x-s).                  (4)
```

Its character multiplier is

```text
m_k=sum_s kappa_s zeta^(-ks),

(Kf)_hat(k)=m_k f_hat(k).                            (5)
```

Write `kappa_s=a_s/D`, with `D>=1` and `a_s` nonnegative integers, and
put

```text
P(X)=sum_(s=0)^(p-1)a_s X^s.                         (6)
```

For `k!=0`, `zeta^(-k)` is a primitive `p`-th root. If `m_k=0`, then
the minimal polynomial

```text
Phi_p(X)=1+X+...+X^(p-1)                             (7)
```

divides `P` over `Q`. Both polynomials have degree at most `p-1`, hence

```text
P=c Phi_p                                             (8)
```

for an integer `c`. Thus all `a_s`, and therefore all `kappa_s`, are
equal. The converse is immediate from the sum of all nontrivial
characters. We have proved:

> **Prime-cyclic rational dichotomy.** For every `k!=0`,
>
> ```text
> m_k=0
> iff
> kappa_0=kappa_1=...=kappa_(p-1).                   (9)
> ```

Since the kernel in (3) is nonzero and nonnegative,

```text
m_0=sum_s kappa_s>0.                                 (10)
```

Therefore exactly one of the following occurs:

```text
kappa is uniform:
  K is a positive rank-one averaging operator;

kappa is nonuniform:
  every m_k is nonzero and K is invertible on C[G]. (11)
```

The proof of (9) works for rational signed kernels as well. Nonnegativity
is used only in (10) and in the physical corollaries below.

## 2. An explicit cyclotomic norm floor

Let

```text
S=sum_s a_s>0,

g=gcd(a_0,...,a_(p-1)).                             (12)
```

If the kernel is nonuniform, the algebraic-integer norm

```text
N=Norm_(Q(zeta)/Q)(P(zeta))
 =product_(j=1)^((p-1)/2)|P(zeta^j)|^2              (13)
```

is a nonzero integer divisible by `g^(p-1)`. Fix `k!=0`. Every other
conjugate has modulus at most `S`, so

```text
g^(p-1)
 <=|P(zeta^k)|^2 S^(p-3).                            (14)
```

Equations (5)--(6) give

```text
|m_k|
 >=g^((p-1)/2)/(D S^((p-3)/2)).                     (15)
```

At `p=13`,

```text
|m_k|>=g^6/(D S^5).                                 (16)
```

If `K` is Markov and the displayed denominator is chosen so that
`S=D` and the fractions are primitive, this becomes

```text
|m_k|>=D^(-6).                                      (17)
```

There is a complementary variance form. Put

```text
V=13 sum_s a_s^2-S^2
 =sum_(k=1)^12|P(zeta^k)|^2.                        (17a)
```

The six conjugate-pair squared moduli have sum `V/2`. Applying AM--GM
to the five pairs other than a fixed `k` gives

```text
|m_k|^2
 >=(g^12/D^2)(10/V)^5.                              (17b)
```

For the normalized joint correlation in Section 4, the corresponding
quantity is divided by `13^2`. Also

```text
V=sum_(r<s)(a_r-a_s)^2>=12g^2                      (17c)
```

for a nonflat integer kernel, and Parseval shows that some normalized
joint mode has magnitude at least `g/(169D)`. The last bound is sharp for
one coefficient differing from the other twelve by `g`.

This is a denominator-sensitive floor. There is no positive
denominator-free floor: the nonuniform Markov kernels

```text
(1-1/N)Pi+(1/N)delta_0                              (18)
```

have every nontrivial multiplier equal to `1/N`.

Two small useful variants require no rationality.

1. If a nonnegative kernel has total mass `tau` and support on at most two
   shifts, then, for `k!=0`,

   ```text
   |m_k|>=tau sin(pi/(2p)).                          (19)
   ```

   Indeed, two distinct `p`-th roots subtend an angle at most
   `pi-pi/p`; at fixed total mass the smallest chord combination occurs
   with equal weights at that largest angle.

2. Translating a nonzero kernel through all `p` shifts rotates `m_k`
   through every `p`-th root. Some translate therefore obeys

   ```text
   Re m_k>=cos(pi/p)|m_k|.                           (20)
   ```

Equation (20) is only a selectable phase sidecar when the whole translation
bank is lawful.

## 3. A common kernel preserves a charged edge

Let `d,w:G->C` be two currents, and apply the same oriented operator `K`
to both. Equation (5) gives, character by character,

```text
(Kd)_hat(k) conjugate((Kw)_hat(k))
 =|m_k|^2 d_hat(k) conjugate(w_hat(k)).              (21)
```

Thus a common nonuniform rational kernel:

- preserves exact same-target support;
- preserves the complex argument and the sign of the real charged
  product; and
- at `p=13` sends a magnitude floor `eta` to

```text
eta/(D^2 S^10).                                     (22)
```

For a denominator-`D` Markov kernel, (22) is `eta D^(-12)`.

The word **common** is load-bearing. Two different nonuniform kernels
preserve support but can rotate the cross phase. For example, with `Pi`
the uniform projection,

```text
K=(1-epsilon)Pi+epsilon delta_0,

L=(1-epsilon)Pi+epsilon delta_s,                    (23)
```

the nonzero-character cross multiplier is

```text
epsilon^2 zeta^(ks),                                (24)
```

whose phase is arbitrary as `s` varies. A common but oppositely oriented
kernel has the same problem. The operator must act in one common endpoint
gauge and with one orientation.

The adjoint square `K^*K` has multiplier `|m_k|^2` and automatically
repairs phase, but current LRC canon supplies no lawful reverse terminal
pass.

## 4. Two-sided survival makes every colour nonzero

There is a direct overlap version which does not posit an operator.
Let `Y` be a finite union of intervals with rational endpoints. Let

```text
E,F:Y x F_p->{0,1}                                  (25)
```

be step functions whose breakpoints in `Y` are rational. Define the
unnormalized cyclic cross-correlation

```text
C(s)
 =integral_Y sum_(r in F_p)
    E(y,r+s)F(y,r) dy.                              (26)
```

Every `C(s)` is a nonnegative rational number. If the packets are
disjoint in the common base gauge,

```text
E(y,r)F(y,r)=0             a.e.,                    (27)
```

then

```text
C(0)=0.                                             (28)
```

Summing (26) over all shifts gives the exact mass identity

```text
sum_s C(s)
 =integral_Y |E_y| |F_y| dy.                        (29)
```

Suppose the right side is positive. Then `C` is nonzero, but (28) says
it is not constant. Applying (9) to its rational weights proves

```text
C_hat(k)!=(0)              for every k!=0,          (30)
```

where `C_hat` uses the normalized convention (2).

For

```text
e_y(k)=(1/p)sum_r E(y,r)zeta^(-kr),

f_y(k)=(1/p)sum_r F(y,r)zeta^(-kr),                 (31)
```

direct finite inversion gives the exact typing

```text
C_hat(k)
 =p integral_Y e_y(k)conjugate(f_y(k))dy.           (32)
```

Thus (30) is a same-parent joint charged overlap at every nonzero colour,
not merely separate energy in the two packets.

The simultaneous-fibre condition in (29) is essential. Positive `E` and
positive `F` supported on disjoint base subsets give `C=0`. Likewise,
currents in different endpoint gauges cannot be inserted into (26).

## 5. The fixed singleton/two-root quantitative corollary

At `p=13`, suppose on a base cell of mass `rho`:

```text
A is a fixed singleton or fixed two-root mask;

F_y is a two-root mask disjoint from A for every y.  (33)
```

For arbitrary two-root masks there are at most

```text
binom(13,2)=78                                      (34)
```

possibilities. In the fixed-unit-comb application the two-root word is a
cyclic translate of one fixed shape in the unit-scaled root order. It
therefore has at most thirteen possibilities. Refine by that translate
and retain a cell `Z` with

```text
mu(Z)>=rho/13.                                      (35)
```

On `Z`, both masks are fixed. A nonempty proper Boolean mask on `C_13`
has no vanishing nonzero Fourier coefficient by (9). More quantitatively,
the smallest normalized coefficient of a two-root mask is

```text
2 sin(pi/26)/13>2/169,                              (36)
```

where the strict inequality uses `sin x>2x/pi` on `(0,pi/2)`. A singleton
coefficient has modulus `1/13>2/169`. Therefore, simultaneously for all
`k!=0`,

```text
|integral_Z a_k conjugate(f_k)dy|
 =mu(Z)|a_k f_k|
>4rho/(13*13^4)
 =4rho/371293.                                     (37)
```

This finite-pattern form does not need rational integration: the two masks
are constant on `Z`. It retains one actual two-root label and its complete
root pattern. It does not partition by a full six-role word. Without the
fixed-translate hypothesis, the same proof remains valid with the coarser
factor `78` in (34).

THM-2392 proves that each fixed ordinary lower role has exactly the
two-root translate structure in (33)--(35), disjoint from the exclusive
top mask. THM-2396 supplies the uniform mass inputs. Thus (37) gives the
proved exact floors

```text
rho>=1/1391208:
  every colour >1/129136447986;

rho>=33/115934:
  every colour >66/21522741331;

rho>=33/753571:
  every colour >132/279795637303.                   (38)
```

If the last cell also carries a fixed septimal coefficient of modulus
`1/7`, its derived joint `C_7 x C_13` tensor floor is

```text
132/1958569461121.                                 (39)
```

This application does not use THM-2397. It retains one actual named
lower role on all twelve root colours. The current THM-2397 proof
candidate gives complementary signed and aggregate information through
derived least-role selectors; it is related rather than a dependency.
The coefficients in (37)--(39) are fibrewise joint mixed coefficients.
A product of separately integrated aggregate currents would contain one
additional factor of `mu(Z)`.

## 6. Higher-dimensional and composite boundaries

For a rational kernel on `F_p^d`, fix `q!=0` and collect the kernel mass
on the `p` character fibres:

```text
w_r=sum_(
  s:<q,s>=r
) kappa_s.                                           (40)
```

Then

```text
m_q=sum_r w_r zeta^(-r).                            (41)
```

The same minimal-polynomial proof gives

```text
m_q=0
 iff
w_0=w_1=...=w_(p-1).                               (42)
```

Thus global nonuniformity in `F_p^d` is insufficient; the projection along
the selected character must be nonuniform. This is the correct form for a
two-dimensional target bank.

Primality is also essential. On `C_4`, the nonnegative anchored kernel

```text
(0,1/2,0,1/2)                                      (43)
```

is nonzero and nonuniform but has zero character `k=1`. A strictly
positive rational example exists on `C_91`: mix full uniform measure with
uniform measure on the order-seven subgroup

```text
{0,13,26,39,52,65,78}.                              (44)
```

The result is nonuniform and strictly positive, while its primitive
character `k=1` vanishes.

Rationality cannot be weakened to arbitrary real positivity. For
`0<epsilon<1`, the reversible, strictly positive, nonuniform Markov kernel

```text
kappa_s
 =(1+epsilon cos(4*pi*s/13))/13                    (45)
```

has `m_1=0` by ordinary Fourier orthogonality. This is the exact first
failed implication in any attempted Perron--Frobenius replacement of
(9).

The rationality boundary persists even with the pinned anchor and positive
two-sided survival from Section 4. On three or more irrational base
strata one can realize the following exact correlation:

```text
C(0)=0,

C(s)=(11+2cos(2*pi*s/13))/130,       s!=0.          (45a)
```

All twelve displayed weights are positive and sum to one. Yet

```text
sum_s C(s)zeta^(-s)
 =11(-1)+(12-1)=0.                                  (45b)
```

Indeed, on the stratum of mass `C(s)` take `E=delta_s` and
`F=delta_0`. The packets are disjoint and coexist on every base fibre;
only rational chamber lengths have been removed.

## 7. The LRC interface and the surviving debt

The mechanism connects to proved canon as follows.

- THM-2300 §1 is the weighted proper-root precursor: a rational kernel
  with a zero entry cannot be flat.
- THM-2365 §6 uses the same rational cyclotomic mechanism on its
  anchored drift profile.
- THM-2368 §3 proves the Boolean proper-word unit special case.
  Its later pointwise-versus-integrated warning remains valid.
- THM-2370 §6 permits a linear decomposition only while a common right
  endpoint is retained. Replacing it by the canonical fully masked
  endpoint introduces cross-layer terms.
- THM-2380 identifies the charged product in (21), but does not supply
  the common terminal operator.
- THM-2392 and THM-2396 supply the clean-root application in
  (33)--(39), including its uniform mass.

Accordingly, the exact positive bridge would be:

```text
clean-parent charged edge
  + one common rational C_13-equivariant terminal operator
  + any proof that its 13 transition weights are not all equal
  -> terminal charged edge.                         (46)
```

The cheapest decisive test is one `13 x 13` transition table on a fixed
owner/status cell:

1. verify that both currents use the same oriented circulant table;
2. verify rationality;
3. find one forbidden shift or any unequal pair of transition weights.

The sole remaining algebraic hostile is the flat mean projection. The
current canon has not constructed the common operator in (46).

Two exact hostiles isolate the missing geometric graft.

1. Separate terminal masses and separate full spectra do not give
   co-support. Split a rational base in halves; take

   ```text
   (E,F)=(delta_0,0)
   ```

   on one half and `(0,delta_1)` on the other. Both aggregate packets
   have every nonzero root character, but the joint correlation is zero.

2. Positive co-support without the pinned common anchor can still be flat.
   On thirteen equal rational cells indexed by `t`, take

   ```text
   E=delta_t,          F=delta_0.
   ```

   Then `C(s)=1/13` for every shift, so every charged mode vanishes.

Thus the precise missing service is:

```text
one common oriented endpoint/root gauge
  + one common Boolean terminal filter
  + positive simultaneous descendant mass.          (47)
```

Common filtering preserves the clean-parent disjointness, hence `C(0)=0`;
(29)--(32) then produce all twelve terminal colours automatically. No
proved result among THM-2305, THM-2370, THM-2380, or THM-2392--2397
supplies the positive aligned descendant mass in (47). The
finite-pattern coefficient in (37) remains a clean-parent joint current,
not a canonical terminal word or scalar-frequency atom. No scalar profile
is excluded. The ledger remains `165`, and LRC(14) remains open.

## 8. Exact companion

The dependency-free exact companion:

- exhausts coefficient alphabets `{0,1,2}` at primes `3,5,7` and all
  `8,192` Boolean kernels at prime `13`, finding charged zeros exactly
  on the flat kernels;
- exhausts all `1,932` nonempty disjoint two-packet assignments on
  `C_7`;
- checks all `121` singleton/adjacent versus disjoint-two-root patterns
  used in (33)--(37);
- evaluates five exact `12 x 12` cyclotomic multiplication determinants
  and checks the gcd, variance, Parseval, and paired norm floors in
  (15)--(17c);
- verifies the separate-spectrum and flat-restoration terminal hostiles;
- verifies the composite anchored hostile and numerically replays the
  explicitly analytic strict-positive hostiles, including (45a); and
- checks every reduced fraction in (37)--(39).

Run

```bash
python3 04-computation/lrc14_prime_cyclic_rational_restoration_thm2398.py
python3 -O 04-computation/lrc14_prime_cyclic_rational_restoration_thm2398.py
```

Both transcripts must byte-match, after LF normalization,

```text
05-knowledge/results/lrc14_prime_cyclic_rational_restoration_thm2398.out
```

All finite truth-bearing assertions use explicit raising checks, so optimized
mode executes the same audit.
