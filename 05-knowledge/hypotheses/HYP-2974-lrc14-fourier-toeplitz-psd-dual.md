---
id: HYP-2974
title: LRC14 Fourier-Toeplitz PSD dual
status: PROOF-INTERFACE / finite Fourier dual evidence, not a proof
source: codex-2026-06-24-S157
script: 04-computation/lrc14_fourier_toeplitz_psd_dual_codex_s157.py
result: 05-knowledge/results/lrc14_fourier_toeplitz_psd_dual_codex_s157.out
related:
  - HYP-2973
  - HYP-2972
  - HYP-2971
  - HYP-2970
  - HYP-2969
  - HYP-2968
  - HYP-2966
  - HYP-2965
  - HYP-2963
  - HYP-2956
  - HYP-2901
  - THM-406
  - THM-534
  - THM-572
  - OPEN-Q-108
---

# HYP-2974: LRC14 Fourier-Toeplitz PSD Dual

This completes the reserved Fourier version of the multiplicity-dual route.
For a 13-speed row `S`, let

```text
C_S(t) = #{v in S : ||v t|| < 1/14}.
```

A strict LRC14 counterexample would force the open danger arcs to cover the
circle, hence

```text
F_S(t) = C_S(t) - 1 >= 0
```

almost everywhere.  Therefore every finite Toeplitz moment matrix built from
the Fourier coefficients of `F_S` must be positive semidefinite:

```text
T_d(S) = [hat F_S(i-j)]_{0<=i,j<=d} >= 0.
```

If any low-degree `T_d(S)` has a negative eigenvalue, the associated
trigonometric square `|P(e^{2 pi i t})|^2` gives a dual certificate:

```text
integral F_S(t) |P(e^{2 pi i t})|^2 dt < 0,
```

so `F_S` is not nonnegative and `S` has a strict safe interval.

## Explicit Coefficients

Let `B(t)=1_{||t||<1/14}`.  Its Fourier coefficients are

```text
Bhat(0) = 1/7,
Bhat(l) = sin(pi*l/7)/(pi*l),  l != 0.
```

Since `1_{||s t||<1/14}` is the pullback of `B` by multiplication by `s`,

```text
hat F_S(0) = |S|/7 - 1 = 6/7,
hat F_S(k) = sum_{s in S, s|k} sin(pi*(k/s)/7)/(pi*(k/s)),  k != 0.
```

This is a closed-form harmonic test, not a grid sample.

## Dual Certificate

If `T_d(S)` has vector `c` with

```text
c^* T_d(S) c < 0,
```

set

```text
P(t) = sum_{j=0}^d c_j exp(2*pi*i*j*t).
```

Then

```text
int_0^1 F_S(t) |P(t)|^2 dt < 0.
```

If `F_S>=0` a.e. this is impossible.  Therefore `F_S<0` on positive
measure.  Since `C_S` is integer-valued off finitely many endpoints, this
means `C_S=0` on a positive-measure open set, so `M(S)>1/14`.

This is the clean Fourier/Farkas certificate suggested by the prompt.

## Computation

Script:

```text
04-computation/lrc14_fourier_toeplitz_psd_dual_codex_s157.py
```

Stored output:

```text
05-knowledge/results/lrc14_fourier_toeplitz_psd_dual_codex_s157.out
```

Default run scans named hard rows through degree `160` with eigenvalue
tolerance `1e-9`.

Readout:

```text
rows tested                         = 15
no PSD failure through degree 160   = AP, GW 12->24
positive hard rows with PSD failure = 13
```

First negative Toeplitz degree examples:

```text
repair drop(2,6)->add(17,42)   degree 32
petal 13->26                   degree 37
covering 12->84                degree 51
covering 12->168               degree 51
hard one-swap 12->48           degree 53
hard one-swap 6->69            degree 56
repair drop(4,6)->add(19,42)   degree 57
loose q-witness 12->26         degree 59
petal 10->20                   degree 70
near K33 12->36                degree 101
P10+GW two-swap                degree 160
```

AP and GW remain PSD to numerical precision; their least eigenvalues approach
zero from above on the tested range.  This matches their role as closed
boundary cover atoms rather than strict open covers.

## Relationship To HYP-2973 And HYP-2971

HYP-2973 keeps only the distribution of the danger-count random variable
`C_S(t)`.  HYP-2971 uses the scalar multiplicity histogram of
`X_S(t)=C_S(t)` and searches admissible polynomial barriers.  HYP-2974 keeps
harmonic location through the Toeplitz matrix of `C_S-1`.

They are complementary:

```text
HYP-2973: S -> law(C_S) -> polynomial lower bound for safe_mu
HYP-2971: S -> law(X_S) -> admissible sign barrier
HYP-2974: S -> hat F_S -> Toeplitz T_d -> |P(t)|^2
```

This is the curried-function form of the Fourier proof:

```text
S |-> (P |-> int (C_S-1)|P|^2).
```

The scalar moment barriers can be lower degree for some rows (`12->36` is
certified by HYP-2971 before Toeplitz degree `101`).  The Toeplitz route buys
localization.  The negative eigenvector should point toward the safe interval
in harmonic space, so it can be reattached to endpoint-credit, twist-ladder,
or state-lift packets.

## Theorem Target

The Fourier theorem target is:

```text
Toeplitz PSD rigidity for LRC14.

If S is a primitive 13-speed row and every finite Toeplitz section of
F_S=C_S-1 is PSD, then S is an AP/Goddyn-Wong closed boundary atom or the row
emits a retained K33/H=7 state-lift packet.
```

Equivalently, every non-AP/GW row should have a finite negative Toeplitz
section unless it routes to the HYP-2908/THM-572 endpoint.

This route is phase-sensitive in a way HYP-2973 deliberately is not.  HYP-2973
keeps only the distribution of the count random variable `C_S`; HYP-2974 keeps
the Fourier locations of that count function through the Toeplitz matrices.

## Tournament Analysis

Assumption challenge: the vertices should not be runners.  The useful vertices
are harmonic proof carriers.

Candidate vertex sets considered:

```text
runners,
danger arcs,
endpoints,
safe components,
Fourier modes,
Toeplitz rows,
negative eigenvectors,
trigonometric-square weights,
moment barriers,
twist ladders,
proof obligations.
```

Chosen vertices:

```text
Toeplitz/Fourier proof carriers:
  PSD finite section,
  negative eigenvector,
  Fourier coefficient ledger,
  multiplicity moment barrier,
  endpoint-credit graph,
  boundary-moment packet ledger,
  raw runner set.
```

Pair observable:

```text
which carrier preserves counterexample necessity, harmonic degree,
localization, compatibility with HYP-2971 barriers, and state-lift visibility.
```

Switch/gauge:

```text
Toeplitz certificate
> eigenvector localization
> Fourier coefficient ledger
> multiplicity histogram
> endpoint packet
> raw row.
```

The proof-carrier tournament is transitive in this run:

```text
directed_3_cycles = 0,
Hamiltonian paths = 1.
```

Destroyed information:

```text
endpoint ownership and exact component labels until the eigenvector is decoded.
```

Preserved LRC predicate:

```text
failure of PSD is a rigorous positive-Haar lonely-time certificate.
```

## Next Work

The next step should decode the negative eigenvectors:

1. recover the trigonometric square `|P(t)|^2`;
2. locate where it concentrates relative to the exact safe intervals;
3. attach endpoint-owner labels at that location;
4. compare the Toeplitz degree to HYP-2971 barrier degree, HYP-2973 count-dual
   degree, and HYP-2972 first twist denominator;
5. test whether the only rows with PSD through a growing degree bound are
   AP/GW boundary atoms.

If this stabilizes, the Fourier route gives a compact dual certificate family
for the NORK bucket: no positive-winding cover can survive all low-degree
Toeplitz tests unless it is already AP/GW or state-lift-owned.
