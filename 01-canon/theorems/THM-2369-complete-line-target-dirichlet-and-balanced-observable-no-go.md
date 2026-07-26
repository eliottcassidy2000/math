---
id: THM-2369
title: "Complete-line target Dirichlet decomposition and balanced-observable no-go"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. Averaging
  THM-2350's magnetic defects over the two complete coordinate lines
  gives D_1=2(E_10+E_11), D_2=2(E_01+E_11), hence recovers the two pure
  target-axis energies and the fork energy exactly. The all-direction
  Cayley average equals total nonzero-target energy. The same formulas
  are physical finite-difference observables on THM-2365's lawful
  tensor, and line circulation transfers at least 1/13 of every positive
  weighted energy to nonzero deep colours. Conversely every balanced
  target-family monomial, including all intensities, autocorrelations,
  and ordinary bispectra, is blind to character slope: constant and
  inverse-character currents have identical balanced data but land at a
  successful target and zero respectively. A charged edge defect is
  sufficient; two independent charged directions fix the modulation
  gauge conditionally on nonzero coefficients. Canon supplies neither
  charged pair service. No scalar-row exclusion, profile decrement, or
  LRC(14) consequence follows.
source: codex-2026-07-25-complete-line-target-dirichlet
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2343-deep-comb-affine-target-catalyst
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
related:
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2361-familywise-fixed-colour-cone-and-offdiagonal-phase-boundary
  - THM-2366-retained-probe-target-covariance-and-subthirteen-budget-bridge
script: 04-computation/complete_line_target_dirichlet_thm2369.py
output: 05-knowledge/results/complete_line_target_dirichlet_thm2369.out
script_sha256: 5b3f52ddef48fe93cb08af3f044e67fb4df4aa70d985ca03ac61fda1c2c44cae
output_sha256: 2e9413e3377eba62f365d413642b84dd775217e7f343dee9fd101755a975bc54
hash_basis: working-tree bytes (LF)
---

# THM-2369 -- complete target lines separate pure and fork energy

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2350 turns target landing into magnetic Dirichlet energy on the
owner-pivot twist torus. Nearest-neighbour edges decide whether that
energy is positive, but their spectral constants are irrational and do
not distinguish a pure target axis from a fork.

Complete coordinate lines remove both limitations. They give integer
weights `0` and `2`, recover all three target-mask energies exactly, and
remain physical overlap-difference observables. They also reveal a sharp
invariant-theoretic boundary: another balanced moment cannot determine
the missing character slope.

## 1. Gauge away the deepest character

Put

```text
G=F_13^2,

Gamma=G^,

|G|=|Gamma|=169.
```

Choose a basis `epsilon_1,epsilon_2` of `Gamma`, and write

```text
q=(q_1,q_2) in G
```

in the dual target coordinates. Let `p_0 in G` be the known
deepest-comb charge from THM-2350. For a phase-free target response
`K:Gamma->C`, define

```text
H(ell)=chi_ell(p_0)K(ell),                         (1)

A(q)=1/169 sum_(ell in Gamma)
     conjugate(chi_ell(q))H(ell).                  (2)
```

Thus

```text
H(ell)=sum_q chi_ell(q)A(q),
```

and `A` is the actual target array after the lawful deepest modulation.
Put

```text
E_tar=sum_(q!=0)|A(q)|^2.                          (3)
```

The bad THM-2343 inverse character is precisely

```text
K(ell)=c chi_ell(-p_0),
```

for which `H` is constant and `E_tar=0`.

## 2. Complete coordinate-line energies

For `j=1,2`, define the covariant complete-line energy

```text
D_j
 =1/(13*169)
  sum_(ell in Gamma)sum_(a in F_13)
  |chi_(a epsilon_j)(p_0)K(ell+a epsilon_j)
    -K(ell)|^2.                                    (4)
```

Multiplication by `chi_ell(p_0)` turns (4) into the ordinary line
difference of `H`. For a target character `q`, the exact multiplier is

```text
1/13 sum_a |chi_(a epsilon_j)(q)-1|^2
 =0                    if q_j=0,
 =2                    if q_j!=0.                 (5)
```

Define the three disjoint target-mask energies

```text
E_10=sum_(q_1!=0,q_2=0)|A(q)|^2,

E_01=sum_(q_1=0,q_2!=0)|A(q)|^2,

E_11=sum_(q_1!=0,q_2!=0)|A(q)|^2.                 (6)
```

Normalized Parseval and (5) give the exact identities

```text
D_1=2(E_10+E_11),

D_2=2(E_01+E_11),                                  (7)

D_line:=D_1+D_2
 =2E_tar+2E_11.
```

Consequently

```text
2E_tar<=D_line<=4E_tar,                            (8)

E_10=E_tar-D_2/2,

E_01=E_tar-D_1/2,

E_11=(D_1+D_2)/2-E_tar.                            (9)
```

The lower equality in (8) means there is no fork target. The upper
equality means there is no pure-axis target. Unlike a tournament
encoding, `(D_1,D_2,E_tar)` preserves amplitudes, both labelled axes,
and the pure/fork split.

## 3. The complete Cayley bank is exactly flat

Averaging every nonzero direction removes even the factor-two
conditioning:

```text
E_tar
 =1/(169*338)
  sum_(ell in Gamma)sum_(d in Gamma;d!=0)
   |H(ell+d)-H(ell)|^2.                            (10)
```

For `q!=0`,

```text
sum_(d!=0)|chi_d(q)-1|^2=338,
```

while the constant target has multiplier zero. Thus (10) is exact, not
an inequality.

The coordinate-line pair is normally cheaper than the full bank: it
requires `2*13*169` differences and additionally reports target-mask
type. The full bank is useful as a symmetry-complete hostile control.

## 4. Physical form on THM-2365's tensor

Let `mathcal H(r,s,t)` be THM-2365's nonnegative lawful overlap tensor,
with normalized squared norm on `F_13^3`. Its target action is

```text
(T_(u,v)mathcal H)(r,s,t)
 =mathcal H(r+v,s+u,t+v).                          (11)
```

Define

```text
L_1=1/13 sum_u
    ||mathcal H-T_(u,0)mathcal H||_2^2,

L_2=1/13 sum_v
    ||mathcal H-T_(0,v)mathcal H||_2^2.            (12)
```

THM-2365's actual target is `q=(b,a+h)`. Therefore

```text
L_1=2(E_10+E_11),

L_2=2(E_01+E_11),                                  (13)

2D_H<=L_1+L_2<=4D_H.
```

These are finite rational mass observables: no Fourier coefficient or
terminal phase needs to be selected to evaluate them.

For every fixed nonzero target `q`, the THM-2365 deep-colour entries
sum to zero. Cauchy--Schwarz gives

```text
energy at a=0
 <=12*(energy at a!=0).
```

Since every weight in (5) depends only on `q`, at least `1/13` of any
positive `D_1`, `D_2`, or `D_line` energy is carried by deep colours
`a!=0`. THM-2365 then extracts a `91`-unit deep multiplier.

The identities remain true after retaining target-neutral THM-2364
probe colours: apply them at each fixed probe colour and sum. This does
not prevent cancellation when those colours are collapsed back to the
canonical current, exactly as recorded in THM-2366.

## 5. Every balanced observable loses character slope

Let the target-character modulation group act by

```text
K(ell)->chi_ell(t)K(ell),             t in G.      (14)
```

A monomial

```text
product_i K(ell_i)
product_j conjugate(K(m_j))
```

acquires the charge

```text
chi_(sum_i ell_i-sum_j m_j)(t).                    (15)
```

It is invariant whenever the index charge in (15) is zero. Thus all
intensities, autocorrelations, closed-index moments, and ordinary
bispectra are blind to character slope.

The sharp pair is

```text
K_good(ell)=c,

K_bad(ell)=c chi_ell(-p_0),             c!=0.      (16)
```

After the deepest modulation (1),

```text
K_good -> A=c delta_(p_0),

K_bad  -> A=c delta_0.                             (17)
```

The first lands on the prescribed deepest target; the second is the
zero-only hostile. Yet every balanced monomial agrees between them. In
particular, for every `ell,m`,

```text
K(ell)K(m)conjugate(K(ell+m))
 =|c|^2 c                                           (18)
```

in both cases. Even the complete target-family bispectrum bank cannot
distinguish success from failure.

This is the same structural lesson as the Gaussian
moment/nullcone interface: weight-zero invariants classify orbit
closure but do not choose a nonzero covariant weight. Here target
landing requires a charged covariant, not another balanced moment.
This is a mechanism analogy, not an import of a Gaussian theorem into
LRC.

## 6. The minimal charged repair

The simplest charged quantity is one directed edge

```text
P_j(ell)=K(ell+epsilon_j)conjugate(K(ell)).         (19)
```

Its exact covariant defect is

```text
|chi_(epsilon_j)(p_0)K(ell+epsilon_j)-K(ell)|^2

 =|K(ell+epsilon_j)|^2+|K(ell)|^2
  -2 Re(
    chi_(epsilon_j)(p_0)P_j(ell)
   ).                                              (20)
```

One edge with positive defect proves `E_tar>0`. A merely nonzero edge
does not: every edge of `K_bad` is nonzero. One vanishing edge also does
not prove the hostile; vanishing defects on a connected spanning tree
do.

Balanced data leave a two-dimensional modulation gauge `t in G`.
Conditionally on nonzero reconstructed coefficients, one complex
charged edge in each of two independent directions determines the two
coordinates of `t`, and is necessary in rank. This is a conditional
gauge-fixing statement, not a claim that canon supplies either edge.

## 7. Why a diagonal reference is not yet the edge

Suppose a lawful nonzero reference `R` and cyclic pair-twist energies

```text
Q_ell(z)=|R+zeta^z K(ell)|^2
```

were available. With the convention

```text
Qhat_ell(1)
 =1/13 sum_z Q_ell(z)zeta^(-z),
```

finite Fourier inversion gives

```text
Qhat_ell(1)=conjugate(R)K(ell),                    (21)
```

so two such values manufacture (19), up to the known factor `|R|^2`.

THM-2361 proves nonvanishing for a diagonal fixed-colour reference, but
does not supply the same-word, same-target pair-twist energies in (21).
Its exact `N=169` hostile has diagonal current `1/2` and raw
off-diagonal current zero. Therefore reference nonvanishing alone is
not a charged-edge service.

## 8. Scope and next test

The theorem gives three exact upgrades:

```text
complete line differences
  <-> pure/fork target-energy decomposition;

all-direction differences
  <-> total target energy;

balanced-observable bank
  -> exact character-slope blindness.              (22)
```

It does not prove a positive line difference on a canonical row, realize
the pair-twist service, prevent probe-colour cancellation, mark a prior
triangle, exclude a scalar profile, or prove LRC(14). The ledger remains
`165`.

The cheapest decisive physical tests are now the `338` nearest-neighbour
dipole defects of THM-2350 or, when target-mask typing matters, the two
complete line energies `(L_1,L_2)`. The minimal phase sidecar is a
same-word charged edge in each independent dipole direction.

## 9. Exact companion

The dependency-free companion uses `Fraction` arithmetic and an exact
`Q[zeta_13]` representation to:

- verify all `338` coordinate-line character multipliers;
- reconstruct `E_10,E_01,E_11` from exact Gaussian-rational target
  arrays;
- verify all `169` full-Cayley multipliers and the factor `338`;
- check the sharp `1/13` deep-colour transfer and equality profile;
- exhaust the constant/inverse-character balanced-bispectrum hostile;
- verify the rank-two charged-edge modulation map and spanning-tree
  boundary; and
- check the cyclic reference polarization (21).

Run

```bash
python3 04-computation/complete_line_target_dirichlet_thm2369.py
python3 -O 04-computation/complete_line_target_dirichlet_thm2369.py
```

Both transcripts must match

```text
05-knowledge/results/complete_line_target_dirichlet_thm2369.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independent audit is pending. QED.
