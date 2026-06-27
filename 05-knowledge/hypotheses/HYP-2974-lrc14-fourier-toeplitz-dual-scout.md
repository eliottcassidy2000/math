---
id: HYP-2974
title: LRC14 Fourier-Toeplitz dual scout, PSD refinement, and Fejer full-bank extension
status: PROOF-INTERFACE / harmonic dual necessary condition with degree-280 Fejer full-bank extension; not a proof
source: codex-2026-06-24-S157 + codex-2026-06-24-S156
related:
  - HYP-2977
  - HYP-2976
  - HYP-2975
  - HYP-2973
  - HYP-2972
  - HYP-2971
  - HYP-2970
  - HYP-2969
  - HYP-2968
  - HYP-2967
  - HYP-2966
  - HYP-2965
  - HYP-2964
  - HYP-2963
  - HYP-2961
  - HYP-2956
  - HYP-2954
  - HYP-2953
  - HYP-2901
  - HYP-2908
  - THM-406
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_fourier_toeplitz_dual_scout_codex_s156.py
  - 05-knowledge/results/lrc14_fourier_toeplitz_dual_scout_codex_s156.out
  - 04-computation/lrc14_fourier_toeplitz_psd_dual_codex_s157.py
  - 05-knowledge/results/lrc14_fourier_toeplitz_psd_dual_codex_s157.out
  - 04-computation/lrc14_fourier_toeplitz_fejer_fullbank_codex_s157.py
  - 05-knowledge/results/lrc14_fourier_toeplitz_fejer_fullbank_codex_s157.out
  - 04-computation/lrc14_packet_fejer_interval_scaffold_codex_s162.py
  - 05-knowledge/results/lrc14_packet_fejer_interval_scaffold_codex_s162.out
---

# HYP-2974: LRC14 Fourier-Toeplitz Dual Scout, PSD Refinement, And Fejer Full-Bank Extension

This hypothesis fills the S157 Fourier-Toeplitz PSD stub with the S156
computation and the later S157 full-bank Fejer extension.  It is the Fourier/Toeplitz member of the current
dual-certificate cluster, alongside HYP-2970 endpoint-credit winding cycles, HYP-2971
multiplicity-moment barriers, HYP-2972 twist ladders, and HYP-2973
danger-count moment duals.  After rebasing over HYP-2977, distinguish the two
Fourier lanes this way: HYP-2974 tests the nonnegativity forced by danger-cover
via Toeplitz moments of `C_S-1`, while HYP-2977 tests positive strict-safe mass
through Fourier/Fejer shadows of `1_U`.  Instead of classifying the residual
first, HYP-2974 looks directly at the covering formulation.

For a speed set `S`, define the danger multiplicity

```text
D_S(t) = sum_{v in S} 1_{||v t|| < 1/14}.
```

If `S` were a strict LRC14 counterexample, then the danger arcs cover the circle
and `F_S(t)=D_S(t)-1 >= 0` almost everywhere.  Hence every Toeplitz moment matrix

```text
T_K(S) = (hat F_S(i-j))_{0 <= i,j <= K}
```

must be positive semidefinite.  Any negative eigenvalue is a Fourier/Farkas
dual certificate that `S` has a safe open interval.  This is a necessary
condition for counterexamples, independent of endpoint-owner packet labels.

The coefficients used by the S157 refinement are explicit.  For
`B(t)=1_{||t||<1/14}`,

```text
Bhat(0) = 1/7,
Bhat(l) = sin(pi*l/7)/(pi*l),  l != 0.
```

Since `1_{||s t||<1/14}` is the pullback of `B` by multiplication by `s`,

```text
hat F_S(0) = |S|/7 - 1 = 6/7,
hat F_S(k) = sum_{s in S, s|k} sin(pi*(k/s)/7)/(pi*(k/s)),  k != 0.
```

If `T_d(S)` has vector `c` with `c^*T_d(S)c<0`, then

```text
P(t) = sum_{j=0}^d c_j exp(2*pi*i*j*t)
```

satisfies `int F_S(t)|P(t)|^2 dt < 0`.  If `F_S>=0` a.e. this is impossible;
therefore `C_S=0` on positive measure, and the row has a strict lonely
interval.

## Evidence

The default scout audits `52` rows: curated AP/GW, K33, petal, covering,
few-apex, and AP-mutation rows plus deterministic qdiv>=14 samples.

```text
rows audited                 52
boundary-only rows           2
positive-open rows           50
dual-certified rows          48
PSD-through-degree-90 rows    4
```

The two boundary-only rows, AP and Goddyn-Wong `12->24`, remain PSD through
degree `90`.  That is exactly what the dual condition should do: it does not
falsify known equality atoms.

The two positive-open rows not certified by degree `90` are also meaningful:
`K33 near 12->36` and `two-hole P10+GW`.  Those are already named exits in the
existing packet stack.  Thus the new harmonic lens appears to discharge many
covering/migration rows while naturally handing K33/petal exceptions back to
the state-lift and petal-rigidity machinery.

The S157 named-row refinement pushes the Toeplitz sections to degree `160` on
the hard rows singled out by the user prompt and the HYP-2971/HYP-2972 dual
cluster:

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

Thus the degree-90 scout's two harmonic-invisible positive rows are not
invisible to the full named-row Toeplitz test: `12->36` turns negative at
degree `101`, and `P10+GW` at degree `160`.  AP and GW still remain PSD to
numerical precision on this range, as expected for closed boundary atoms.

The S157 full-bank extension avoids NumPy and tests an explicit Fejer vector
instead of computing eigenvalues.  For each positive-open row, choose a point
`x` in the largest exact safe component and use

```text
p_j = exp(-2*pi*i*j*x) / sqrt(d+1),  0 <= j <= d.
```

Then

```text
Q_d(x) = c_0 + 2*sum_{k=1}^d (1-k/(d+1))*c_k*cos(2*pi*k*x)
```

is the Toeplitz quadratic form of that vector.  A negative value is still a
valid PSD violation, though it is only one possible test vector and not the
least eigenvalue.

Default full-bank audit over the HYP-2963 labelled-packet bank:

```text
rows audited                 21913
zero-safe rows               2
positive-open rows           21911
Fejer PSD-vector hits        21911
misses at cap degree 512     0
max first degree             280 at P10+GW
```

The named hard rows are now all hit by explicit Fejer vectors:

```text
near/K33 12->36       degree 159
petal 10->20          degree 115
petal 13->26          degree 65
P10+GW                degree 280
P10+K33               degree 124
covering 12->84       degree 64
covering 12->168      degree 63
few-apex 6->14        degree 38
few-apex 6->28        degree 106
```

This does not replace the eigen-scan.  It gives a stronger computational
message for the labelled-packet bank: every positive row tested has an explicit
trigonometric-square dual witness at modest degree, while AP/GW are exactly the
two equality atoms.

The S162 scaffold begins the formalization step requested after the floating
audit.  It rewrites selected Fejer quadratic forms as rational interval
certificates and attaches each certificate to the labelled packet fiber

```text
P(S) = (route, family, q_class, packet_route, state_lift, q_threshold).
```

The first selected packet-anchored interval certificates all have upper endpoint
strictly below zero:

```text
near/K33 12->36                    degree 159  P(S)=K33-STATE-LIFT
P10+GW                             degree 280  P(S)=BOUNDARY-PETAL-SPORADIC
covering 12->168                   degree  63  P(S)=COVERING-MOMENT
two drop(12,13)->add(14,29)        degree  41  P(S)=Q-WITNESS
single swap 6->63                  degree 266  P(S)=COVERING-MOMENT
```

The interval backend is exact-rational around a rational enclosure of `pi` and
Taylor enclosures for `sin(pi*m/7)` and `cos(pi*rational)`.  This is still a
scaffold, not a final formal certificate: the hard-coded `pi` enclosure should
be replaced by a formally sourced backend, and the full bank should be emitted
by packet fiber rather than by hand-selected rows.

## Proof Target

A possible theorem shape is:

```text
After AP/GW boundary atoms and named K33/petal exits are removed,
every primitive qdiv>=14 LRC14 residual has a bounded-degree negative
Toeplitz moment certificate for F_S=D_S-1.
```

The bounded-degree phrase is currently empirical.  The proof would need a
structured harmonic argument, likely tracking the unit-apex residue band that
dominates the observed negative eigenvectors.

The sharper S157 target is Toeplitz PSD rigidity:

```text
If every finite Toeplitz section of F_S=C_S-1 is PSD, then S is an AP/GW
closed boundary atom or emits a retained K33/H=7 state-lift packet.
```

Equivalently, every non-AP/GW row should have a finite negative Toeplitz
section unless it routes to the HYP-2908/THM-572 endpoint.  The full-bank
Fejer audit and S162 interval scaffold sharpen the next formal task: replace
all floating evaluations of `Q_d(x)` by interval-enclosed trigonometric
certificates, grouped by labelled packet family and endpoint-owner data.

## Tournament Analysis

Vertices are harmonic proof lenses rather than runners, arcs, or packets.
The pairwise observable is which lens preserves the implication

```text
danger cover => F_S >= 0 => Toeplitz PSD
```

while retaining enough labelled information for known exits.  The default
fingerprint is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
c3=0
```

The challenged assumption is that the endgame must be packet-first.  HYP-2974
says a counterexample might first be ruled out by a trigonometric moment
obstruction, with packets reattached only for the few harmonic-invisible
exceptions.

S157 uses the same assumption challenge but makes the carrier tournament more
concrete: vertices are Toeplitz/Fourier proof carriers, not runners.  The
transitive switch/gauge is

```text
Toeplitz certificate
> eigenvector localization
> Fourier coefficient ledger
> multiplicity histogram
> endpoint packet
> raw row.
```

The next useful computation is to decode the negative eigenvectors into
localized trigonometric squares and attach endpoint-owner labels at their
mass concentration sites.  The Fejer full-bank extension offers a simpler
parallel target: certify the explicit Fejer sums by interval arithmetic first,
then use eigenvector localization only where the labelled packet theorem needs
shorter or more structural witnesses.
