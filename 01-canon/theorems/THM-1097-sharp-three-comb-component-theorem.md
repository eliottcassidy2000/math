---
id: THM-1097
title: Sharp three-comb component theorem; the four-killer clustered stratum is uniformly closed
status: PROVED FINITE-EXACT + ANALYTIC — the sharp periodic-comb discrepancy and component count reduce every scale to an explicit guarded bank; an exact rational endpoint referee checks all 39,778,595 triples over all 220 nine-speed cores with zero failures, and an independent Fraction replay verifies the atlas, guard cardinality/boundaries, and hardest finite-bank row
source: codex-2026-07-18-S73 r4 continuation
depends_on: [THM-1094]
related: [THM-1042, THM-1061, THM-1071, THM-1081, THM-1101, MISTAKE-163, MISTAKE-164]
script:
  - 04-computation/lrc14_r4_three_comb_component_exact_codex_S73.cpp
  - 04-computation/lrc14_r4_three_comb_atlas_replay_codex_S73.py
  - 04-computation/lrc14_r4_three_comb_independent_guard_replay_codex_S73audit.py
  - 04-computation/lean/TournamentH7/TournamentH7/LRCSharpCombArithmetic.lean
result:
  - 05-knowledge/results/lrc14_r4_three_comb_component_exact_codex_S73.out
  - 05-knowledge/results/lrc14_r4_three_comb_atlas_replay_codex_S73.out
  - 05-knowledge/results/lrc14_r4_three_comb_independent_guard_replay_codex_S73audit.out
---

# THM-1097 — the sharp three-comb component theorem

For a positive integer `k`, put

```text
D_k = {t in [0,1] : ||kt|| < 1/14}.
```

For a nine-element core `P subset {1,...,12}`, let

```text
S(P) = [0,1] minus union_{p in P} D_p.
```

Open versus closed danger endpoints changes no interval length below.

## Theorem

For every nine-element `P subset {1,...,12}` and every ordered triple

```text
13 max(P) < k1 < k2 < k3,                              (1)
```

the remainder

```text
S(P) minus (D_k1 union D_k2 union D_k3)                (2)
```

has a positive-length component `I` satisfying

```text
|I| > 1/(7k3).                                         (3)
```

Consequently every fourth killer `k4>k3` leaves a point of `I` safe.  Thus
the four-killer clustered stratum is uniformly closed, with no covering
hypothesis and no sampled scaling inference.

## 1. The sharp one-comb discrepancy

Let `J=[u,u+ell]`.  Scale by `k` and write `k ell=m+s`, where `m` is a
nonnegative integer and `0<=s<1`.  The lifted danger comb is 1-periodic and
has one circular danger arc of length `1/7` per period.  Every interval of
integer length `m` therefore has occupied length exactly `m/7`; the remaining
length-`s` interval meets the danger arc in at most `min(s,1/7)`.  Hence

```text
|J intersect D_k|
  <= (m/7+min(s,1/7))/k
  <= ell/7+6/(49k).                                    (4)
```

Indeed the excess `min(s,1/7)-s/7` is `6s/7` for `s<=1/7` and
`(1-s)/7` for `s>=1/7`, so its exact maximum is `6/49`.  The periodic-lift
argument includes the half-teeth at `0` and `1`; there is no boundary caveat.

A tooth meeting `J` has centre `j/k` in the radius-`1/(14k)` enlargement of
`J`.  After scaling, the centre interval has length `k ell+1/7`, so

```text
N_k(J) <= floor(k ell+8/7) <= k ell+8/7.               (5)
```

This is the exact incidence side of the same one-dimensional needle picture.

## 2. The 220-core exact atlas

Every endpoint of `S(P)` belongs to the rational arrangement

```text
(14j-1)/(14p), (14j+1)/(14p),    p in P.               (6)
```

Classifying the cells between consecutive endpoints gives `S(P)` exactly.
If `ell(P)` is its largest component length, the complete
`C(12,9)=220`-core atlas gives

```text
ell(P) >= 11/1008.                                     (7)
```

Equality occurs only for

```text
Pmin={1,2,3,7,8,9,10,11,12}.                           (8)
```

Its two longest components are the reflected pair

```text
[29/126,27/112], [85/112,97/126],                      (9)
```

both of length `11/1008`.  The C++ referee obtains this atlas by successive
exact tooth subtraction.  The independent Python replay instead builds the
full arrangement (6) and classifies each cell by a rational midpoint.

## 3. Three removals: mass, components, and the analytic tail

Fix a largest core-safe interval `J`, of length `ell`.  Put
`a=k1`, `b=k2`, and `c=k3`.  Applying (4) to the three combs gives surviving
mass at least

```text
A = 4ell/7-(6/49)(1/a+1/b+1/c).                        (10)
```

By (5), the number of positive-length survivor components is at most

```text
C <= N_a(J)+N_b(J)+N_c(J)+1
  <= ell(a+b+c)+31/7.                                  (11)
```

Thus some component has length at least `A/C`.  The sufficient inequality
`A/C>1/(7c)` is exactly

```text
Q(ell;a,b,c)
 := 21ell c-7ell(a+b)-6(c/a+c/b)-37 > 0.              (12)
```

The all-scale tail is elementary.  Write `r=b/a>1`, `s=c/b>1`, and
`x=ell a`.  At `x=7`, the left side of (12) is

```text
s(141r-6)-49r-86 > 92(r-1) > 0,                       (13)
```

where the first strict inequality uses `s>1`.  Since the coefficient of
`x` is positive, every triple with

```text
ell k1 >= 7                                             (14)
```

is automatic.  The kernel-pure Lean theorem `threeComb_ratio_tail` checks the
equivalent ratio algebra.

## 4. Exact guard surfaces and exhaustive coverage

Below (14), `Q` is strictly increasing in `c`: its coefficient is

```text
21ell-6/a-6/b > 0.                                     (15)
```

Positivity of every denominator follows already from (1) and (7).  At the
formal boundary `c=b`, condition `Q>0` becomes

```text
b(14ell-6/a) > 7ell a+43.                              (16)
```

For a remaining pair `(a,b)`, solving (12) for `c` gives

```text
c > (7ell(a+b)+37)/(21ell-6/a-6/b).                    (17)
```

The exact referee therefore scans the deliberately guarded superset

```text
13 max(P) < a <= floor(7/ell)+1,                       (18)
a < b <= floor(a(7ell a+43)/(14ell a-6))+1,           (19)
b < c <= floor((7ell(a+b)+37)
               /(21ell-6/a-6/b))+1.                   (20)
```

The `+1` layers are redundant by design.  Every omitted `a` satisfies (14);
every omitted `b` satisfies (16), hence all `c>b` satisfy (12); and every
omitted `c` satisfies (17).  The independent Fraction replay checks all
`23,589` first-level and `1,238,741` second-level guard boundaries.  Its
smallest exact first-omitted margins are

```text
ell*a_omitted-7              = 17/1008,
Q(ell;a,b_omitted,b_omitted) = 31/252,
Q(ell;a,b,c_omitted)         = 28429/182160.            (21)
```

Thus (18)--(20) and the analytic tail cover every integer triple in (1), not
a sample or extrapolated residue window.

The three guard ceilings are all `642`; the largest coordinates that occur
in an actually counted triple are only `(639,640,641)`.  Exact summation of
the guard bank gives

```text
cores                     220
guarded finite triples    39,778,595
```

## 5. The exact endpoint certificate

For every guarded row, the standard-library C++ referee:

1. constructs `S(P)` by exact rational tooth subtraction;
2. removes the `a`, `b`, and `c` teeth in endpoint order;
3. computes the longest surviving component by exact endpoint differences;
4. rejects if the exact integer comparison `7cL<=1` holds.

All endpoint comparisons, cutoff floors, and acceptance decisions use signed
`__int128` cross-products.  No floating-point arithmetic occurs.  The result
is

```text
39,778,595 / 39,778,595 accepted,
0 failures.                                             (22)
```

The hardest row *inside the finite bank* is

```text
P={1,2,3,6,7,8,10,11,12},
(a,b,c)=(162,168,174),
L=1/522,
7cL=7/3,
R:=1/(7cL)=3/7.                                        (23)
```

Its longest components are

```text
[53/252,517/2436], [1919/2436,199/252].                (24)
```

The Python replay reconstructs (23)--(24) from the complete breakpoint
arrangement of all twelve speeds, rather than calling the C++ subtraction
logic.  A separate independent audit replay checks all omitted guard
boundaries, all 220 primary atlas/count rows, the 200 nonempty cores' reported
hardest rows, and the 20 analytic-only rows.  We do not call (23) an
infinite-tail extremum; the tail proves the required strict inequality
without optimizing its margin.

Combining Sections 3--5 proves (3).

## 6. The fourth killer

Let `I` be the component from (3), with length `L`, and let `k4>k3`.  Then

```text
L > 1/(7k3) > 1/(7k4).
```

Applying the sharp discrepancy (4) once more gives

```text
|I intersect D_k4| <= L/7+6/(49k4) < L.               (25)
```

Therefore `D_k4` cannot cover `I`; a positive-measure set of times is safe
for all thirteen speeds.  This last implication is also kernel-checked in
Lean as `sharpFinalComb_mass_lt`.

## 7. Toothpick self-similarity and the r=5 wall

The useful self-similarity is not a literal repetition of the full endpoint
word.  It is the scale-free pair `(occupied mass, tooth incidence)` in
(4)--(5): dilation changes tooth locations and endpoint order, but preserves
the `1/7` density and leaves only a bounded endpoint discrepancy.  The guard
surfaces (18)--(20) isolate the finite region where phase order matters; past
them, mass and component slope alone decide the proof obligation.

This mechanism has an exact phase transition.  After `m` removals at nearly
equal large scales, the survivor-mass slope is `(7-m)/7`, while the component
slope is `m`.  Beating the final-comb target `1/(7k)` requires

```text
(7-m)/m > 1  iff  m < 7/2.                             (26)
```

Thus it reaches `m=3` removals, exactly the present `r=4` theorem, but fails
at `m=4`, the `r=5` problem.  The Lean theorem
`toothpick_phase_transition` verifies (26).  This explains the honest wall
in THM-1101: it is a failure of this sufficient mass/component quotient, not
a counterexample to the lonely runner conjecture.  Uniform `r=5` needs new
endpoint self-similarity, overlap, or arithmetic information.

The Kakeya-needle analogy is one-dimensional and metric: periodic teeth of
several scales pierce a fixed core interval, (4) controls occupied length,
and (5) controls the number of cuts.  Direction-only or order-only data
cannot replace either metric sidecar.

## 8. Carrier and Tournament Analysis audit

The predicate-preserving carrier is the ordered list of core-safe components,
rational boundary coordinates, endpoint owners, and surviving lengths.  The
three combs refine this carrier by wall-crossing events.  Candidate vertices
considered here include runners, teeth, core gaps, section boundaries,
residues, exposed endpoints, and proof obligations.  A runner tournament
forgets which teeth meet which component; a residue tournament forgets
metric length; a plain endpoint tournament remembers only order.

For exposed endpoints, the natural pair observable is left/right order, the
coordinate is the gauge, and exact equality is the tie.  The induced
tournament is transitive: score multiset `{0,1,...,q-1}`, no directed cycles,
singleton SCCs, and one Hamiltonian path (the sorted endpoint word).  An edge
flips only at an endpoint coincidence.  Consequently the tournament has no
information beyond the order and is not faithful without coordinate and
owner sidecars.  The exact interval carrier preserves (3); its tournament
quotient alone does not.

## 9. Reproducibility and formalization boundary

The final source and output hashes are

```text
C++ source SHA-256     afc8875054aa3744dbbb5dadcb886826674a5be7ba7834b2a5c46a1de3b4359e
C++ output SHA-256     c030f563abf659f66acf44f4ce6cae84e493e0dcf4f40aed67f02a10a6591c09
Python source SHA-256  cb588f0ea994d23ed99cb6750437ecc483e361ab598a220ca4ea85f6ea65e4cf
Python output SHA-256  cc06e14d9eb5e4fb006b174ba971f5d01a91c5cbdd4ca391c8f3986c62fb6188
Audit source SHA-256   24185c6ec9d71eefa7c124812f3a1b489519135f002b0c1422d70987a12861cc
Audit output SHA-256   b42da6f169333ae261b4868b25898f3fcb3d29e724e00b6869645cdc7b805832
```

Unoptimized, optimized, and AddressSanitizer+UndefinedBehaviorSanitizer C++
builds produce byte-identical frozen output.  A separate `-Werror` optimized
build run by an independent agent produces the same zero-failure result and
the same hardest row.

The Lean file kernel-checks the ratio-tail algebra, final-killer inequality,
and phase transition with no `sorry`, `native_decide`, or new axioms.  It
does not yet internalize the sharp interval-measure lemma, the 220-core atlas,
or the 39,778,595-row reflected certificate.  Those are the explicit
formalization frontier; the theorem is computer-assisted finite-exact, not a
claim that the finite bank has already been reduced in the Lean kernel.

## Scope

THM-1097 uniformly closes exactly the four-killer clustered stratum.  It
supersedes sampled all-scale claims for that stratum without rehabilitating
them, and it leaves MISTAKE-163/MISTAKE-164 intact.  It does not close the
four-removal `r=5` bridge or prove global LRC(14).
