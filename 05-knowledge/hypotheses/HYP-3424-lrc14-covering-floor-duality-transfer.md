---
id: HYP-3424
title: LRC14 covering-floor duality transfer
status: SYNTHESIS / proof-routing scout; not an LRC14 proof
source: codex-2026-06-28 pass over HYP-3238/HYP-3234/HYP-2272/HYP-2128/HYP-2129 after S259/HYP-3418, integrating HYP-3423 topology-to-magnitude guardrail, HYP-3422 interval relocation, and HYP-3421 off-grid/Rprime branch
tangent: T1385
technique: LTI-385
tournament_technique: LTT-285
script: 04-computation/lrc14_covering_floor_duality_transfer_codex_20260628.py
result: 05-knowledge/results/lrc14_covering_floor_duality_transfer_codex_20260628.out
reflection: 07-reflections/lrc14-covering-floor-duality-transfer-codex-20260628.md
related:
  - HYP-3423
  - HYP-3422
  - HYP-3421
  - HYP-3420
  - HYP-3419
  - HYP-3418
  - HYP-3417
  - HYP-3416
  - HYP-3415
  - HYP-3238
  - HYP-3234
  - HYP-3137
  - HYP-2272
  - HYP-2129
  - HYP-2128
  - THM-414
  - THM-523
  - OPEN-Q-108
---

# HYP-3424: LRC14 Covering-Floor Duality Transfer

## Claim

The older even/odd, odd/even, addition/multiplication, and signed-chart
dualities should now be read through the S259/HYP-3418 correction:

```text
the covering floor is 2-adic;
even speeds are the binding obstruction;
odd/coprime data is a phase reservoir, not a standalone witness.
```

So the productive duality is not a new terminal scalar.  It is a transfer
protocol:

```text
even fold / 2-adic descent
-> signed SPEC decorrelation floor
-> odd phase-cover debt or finite owner-current sidecar
-> HYP-3423 topology-to-magnitude legality guardrail
-> HYP-3422 interval relocation branch
-> HYP-3421 off-grid/Rprime transparency branch
-> off-path filter for 7-adic census data.
```

This is the clean merger of HYP-3238 and S259.  HYP-3238 said:

```text
even/odd          = 2-adic fold versus coordinate resurrection
positive/negative = observer-blind gauge versus pair-visible cut/orientation
```

HYP-3418 says the covering floor binds on the even side at the actual optimum.
Together they say that every floor proof must keep the even fold as a
floor-bearing object and resurrect the odd side only as positional phase-cover
debt.  The old odd/coprime witness at `t=1/2` is exactly the illegal
scalarization: it forgets that all even speeds die there.

## Executable Readout

Script:

```text
04-computation/lrc14_covering_floor_duality_transfer_codex_20260628.py
```

Stored output:

```text
05-knowledge/results/lrc14_covering_floor_duality_transfer_codex_20260628.out
```

Exact anchors retained by the scout:

```text
C(14,2)=91
8*C(14,2)+1=729=27^2
HYP-3418 binding counts at optima: even=376 odd=279
HYP-3418 coprime-transparency test: 0/400
S558 even-fold wall: AP/V* safe slack=0; random min safe slack about 0.10448
HYP-3423 topology guardrail: q-uniform topology cannot close q-specific magnitude floors
HYP-3422 interval relocation: E_safe half-lift branch overlap holds 24/24 on audited rows
HYP-3421 off-grid/Rprime branch: resonance transparency feeds the signed SPEC floor
HYP-3420 owner/chiral sidecar: residue+owner_chiral_class leaves 0 mixed fibers on scanned banks
HYP-3417 frontier current: {2:g2,11:g1,13:g1}
```

The carrier ranking is proof-facing, not aesthetic.  The top path is:

```text
D00 two_adic_even_floor_descent
-> D01 signed_SPEC_decorrelation_floor
-> D02 even_good_odd_phase_cover
-> D03 recursive_quotient_sidecar_router
-> D04 owner_current_even_cover_sidecar
-> D05 additive_multiplicative_energy_transfer
-> D06 parity_redei_odd_sector_guard
-> D07 c3_galois_residue_skeleton_filter
-> D08 raw_7adic_or_scalar_shadow
```

The minimum cover of the tracked proof obligations has size `3`; the first two
minimum covers are:

```text
D00 two_adic_even_floor_descent
D01 signed_SPEC_decorrelation_floor
D03 recursive_quotient_sidecar_router

D00 two_adic_even_floor_descent
D03 recursive_quotient_sidecar_router
D05 additive_multiplicative_energy_transfer
```

This is the useful negative result: the parity/Galois/additive imagery alone
does not cover the floor obligations.  It becomes useful only when paired with
the two-adic descent, SPEC inequality, quotient legality, and finite owner
sidecars.

## Transfer Rules

### 1. Odd/Coprime Becomes Phase Debt

The S259 test `0/400` kills the naive reduction:

```text
R_nr lonely at t=1/2 + resonant transparency
```

The reason is duality-exact: `R_nr` is odd, but every even speed vanishes at
`t=1/2`.  Thus the odd sector is not a witness; it is the phase-cover debt that
can blanket the even-good window.

### 2. Even Fold Is Not "Free"; It Is The Floor Carrier

S558 correctly reduces LRC14 to whether odd danger arcs cover the even-good set
`G`.  S259 corrects the proof emphasis: at the full optimum the even speeds are
the binders.  Therefore `G` is not a solved half to forget.  It is the object on
which the two-adic floor must be proved.

### 3. Relocation And SPEC Are The Two Floor Companions

Incoming HYP-3422 supplies the interval form of the two-adic branch:

```text
E_safe(1/14) cap (odd_branch_0_good union odd_branch_1_good) != empty.
```

Read through HYP-3424, this is the concrete half-lift obligation behind
"even fold / two-adic descent"; it prevents the duality language from staying
at slogan level.

HYP-3421 supplies the companion route on the SPEC side: resonant danger is
grid-local, checked full optima live off-grid, and the canonical `84m` branch
feeds the Rprime chase.  Therefore HYP-3424 treats off-grid transparency as a
signed-SPEC exit, not as permission to revive the failed odd/coprime
reduction.

### 4. Additive Energy Must Become SPEC Penalty

HYP-2272 proves the additive pair-sum face is multiplicative energy of roots of
unity.  HYP-2128/HYP-2129 add the exact triangular/parity bridge

```text
8*C(14,2)+1=27^2
addition certifies multiplication mod 2
```

For the covering floor, this is useful only after translation:

```text
additive-pair abundance
-> multiplicative-energy concentration
-> low-mode SPEC penalty or named odd phase debt.
```

Without that final arrow it is another equality/census shadow.

### 5. Owner Cuts Are Finite Shadows Of The Same Transfer

HYP-3417/HYP-3419 give the local finite shadow:

```text
10->20 frontier current = {2:g2,11:g1,13:g1}
```

Incoming HYP-3420 adds the adjacent chiral-owner readout: on the scanned
expanded banks, `residue_plus_owner_chiral_class` and
`residue_plus_owner_support` leave no mixed theorem-exit fibers, and the
residue-mixed families have size-one owner cuts.  Read through HYP-3424, these
finite owner/chiral facts feed the D04 sidecar branch.  They are important
because the `10->20` frontier current has one even-cover label plus two binding
labels; they are not a substitute for the SPEC/two-adic floor.

### 6. Topology/Galois Is A Legality Filter, Not The Floor

Incoming HYP-3423 supplies the legality rule needed for the user's
`C3/Q(sqrt(-7))/C6 = C2 x C3` synthesis.  C2/Borsuk-Ulam and C6 orbit data are
q-uniform: they can certify residue, symmetry, or equioscillation obligations.
The covering floor is q-specific magnitude data.  Therefore any route that
uses topology, Galois orbits, C3 traces, or apex-7 census data must restore at
least one floor sidecar before concluding:

```text
HYP-3413 q-mod-3 arithmetic,
HYP-3417/HYP-3420 labelled owner-current packets,
S259/HYP-3418 two-adic descent,
HYP-3415/HYP-3421 signed-SPEC/Rprime input,
or HYP-3422 interval relocation.
```

This is why HYP-3424 treats the `C3` skeleton as a routing filter rather than
a terminal covering-floor proof.

### 7. C3/Galois Is A Filter, Not The Floor

The C3/`Q(sqrt(-7))`/apex-7 material organizes the tight/equality/census
regime.  HYP-3415 and HYP-3418 demote it on the covering floor unless it feeds
SPEC, two-adic descent, odd phase debt, finite sidecars, or the HYP-3423
topology-to-magnitude guardrail.

## Candidate Lemma

For every primitive LRC14 covering packet after q-witness and LRC<=13
induction split, one of the following must hold:

1. the even-speed two-adic descent supplies a signed-SPEC floor packet with
   `|SPEC| < product`;
2. HYP-3422 interval relocation supplies an `E_safe` half-lift overlap;
3. the even-good window survives but the odd danger arcs emit a positional
   phase-cover debt;
4. the first mixed quotient is discharged by owner-current/Menger sidecar,
   with even-cover labels explicitly named;
5. HYP-3421 off-grid transparency routes the signed-SPEC/Rprime branch;
6. additive/multiplicative energy concentration converts to a low-SPEC penalty
   or named odd debt;
7. HYP-3423 blocks q-uniform topology/Galois data from closing a q-specific
   floor;
8. the packet is off-path 7-adic/census data and is demoted unless it feeds one
   of the previous exits.

This lemma would not prove LRC14 by itself.  It is a proof-obligation router
whose purpose is to keep future work from proving the wrong dual shadow.

## Tournament Analysis

Vertices are proof obligations and transfer gates, not runners, arcs, residues,
or named constants.

Pairwise observable:

```text
floor-feed + 2-adic payload + odd resurrection + quotient/topology legality + finite sidecar
```

Fingerprint from the scout:

```text
vertices = 9
directed_3cycles = 0
hamiltonian_path =
  D00 -> D01 -> D02 -> D03 -> D04 -> D05 -> D06 -> D07 -> D08
```

## Assumption Challenge

Alternate vertices considered:

```text
runners, gaps, fixed circle sections, section boundaries, wall-crossing events,
residues, cover arcs, Fourier modes/SPEC bins, parity sectors, additive-pair
shells, multiplicative-energy atoms, owner currents, quotient states, proof
obligations.
```

The chosen vertices are proof obligations and duality-transfer gates.  This
preserves the LRC covering-floor predicate only when the quotient feeds the
SPEC floor, the two-adic descent, odd phase debt, finite owner sidecars, or an
explicit HYP-3423 legality/off-path filter.  It destroys raw row identity,
terminal 7-adic numerology, and scalar shadows.
