---
id: THM-3518
title: "U_full all-role phase-normalized dual contraction and cut-cycle certificate"
status: >
  PROVED FINITE-EXACT CERTIFICATE + SECOND INDEPENDENT AUDIT OF THM-3515.
  The five U_full refined role-class bucket tensors reconstructed by direct
  three-character tau slices agree entrywise with a separately contracted
  canonical guard-kernel route.  All 3042 atom-pair/frequency translation
  identities carry the required zeta^(q_t a_L) phase; explicit opposite-sign
  and 13^-2-normalization hostiles disagree.  A rank-five minor certifies the
  5x52 tensor, rank-four minors certify every drift and Fourier slice, and
  explicit sixteen-tree, reduced-Laplacian, and full-graph determinants agree
  in every bank and chart.  The graph incidence rank is seven, its cycle rank
  is six, the forced bridge lies in no cycle, and all 56,592 tested cycle
  pairings vanish.  This strengthens the audit certificate but not the scope:
  endpoint B1 only, with no ancestry, current, absolute H1 flux, bispectrum,
  scalar-row exclusion, or LRC(14).
source: codex/ufull-allrole-phase-audit/2026-08-16
audit: >
  Candidate code and its bucket/role helpers are not imported.  The audit
  imports only the pinned THM-3514 atom-table engine, then independently
  materializes tau/chamber/drift slices, contracts translated guard kernels,
  performs forward/inverse drift transforms, exhibits rank minors, enumerates
  each K4's sixteen spanning trees, evaluates independent Laplacian
  determinants, constructs six fundamental cycles, checks every cycle
  pairing, reconstructs the complete THM-3515 candidate semantic digest, and
  passes normal/optimized/stored replay and documentation gates.
depends_on:
  - THM-3514-ufull-owner-boundary-k4xf13-endpoint-factorization-and-walsh-spectrum
  - THM-3515-ufull-all-role-endpoint-weighted-tree-spectral-closure-in-b1
related:
  - THM-3479-literal-half-twist-relation-current-two-transplant-certificate
  - THM-3482-private-count-gradient-weighted-spectral-closure-without-absolute-h1-flux
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
script: 04-computation/lrc_ufull_all_role_modewise_spectral_independent_audit_20260816.py
output: 05-knowledge/results/lrc_ufull_all_role_modewise_spectral_independent_audit_20260816.out
script_sha256: 975b48c148b88a04b9fea07dcb6352cfaa98591f14eace1a38cd581a5e5876dc
output_sha256: 6eb3f28c07e99c352f1d88615d6beddeb82c264892a9f4735a420b75bf932c19
semantic_sha256: b745ab5d95936b9a134cfe8b6ffde032e72020d950670ab8b1a4643dd58c8be6
hash_basis: LF-normalized UTF-8 bytes; exact finite-field semantic ledger
---

# THM-3518 -- the all-role endpoint spectrum has a phase-normalized dual reconstruction

**PROVED FINITE-EXACT CERTIFICATE + SECOND INDEPENDENT AUDIT OF THM-3515.**

This theorem strengthens the verification record of THM-3515.  It does not
enlarge THM-3515's endpoint-only mathematical scope.

## 1. Fixed atom bank and five characters

Use THM-3514's unguarded endpoint atom tables

```text
A_(alpha,beta)(a,C), B_(alpha,beta)(a,C),
a in F_13, C in {L,M,R},                           (1)
```

over the certified field

```text
P=572252886246508880869,
zeta=505438565698892403012,       ord(zeta)=13.      (2)
```

The speed-thirteen owner gate makes all middle-chamber tables zero.  The
five character classes inherited from THM-3479 and THM-3515 are

```text
Q={(0,0,0),(0,1,0),(1,0,0),(1,0,1),(1,12,0)}.      (3)
```

Their eight role labels and normalized global values are exactly

```text
c1          -> (0, 0,0) -> 405336876493642499425,
c3          -> (0, 1,0) -> 518539850465495448196,
c2,q3,q4,q5 -> (1, 0,0) -> 503604956476841920373,
H           -> (1, 0,1) -> 320618948602619577408,
q2          -> (1,12,0) ->  15703541686881447885.  (4)
```

No value in `(4)` is injected as a bucket coefficient.  Both routes below
reconstruct every bucket from the atom tables and recover `(4)` only after
summing the buckets.

## 2. Direct tau-slice route

Write the guard-safe indicator on chamber `C` and sheet `a` as

```text
g_C(a+tau) in {0,1}.                                (5)
```

For chambers `C,D`, drift `d`, and one shift `tau`, materialize first

```text
S_(alpha,beta,tau)(C,D,d)
 = sum_(b-a=d)
     A_(alpha,beta)(a,C) B_(alpha,beta)(b,D)
     g_C(a+tau) g_D(b+tau).                         (6)
```

This is formed before any character contraction.  For
`q=(q_a,q_b,q_t)` define

```text
B^dir_q(C,D,d)
 =13^(-3) sum_(alpha,beta,tau in F_13)
   zeta^(beta-alpha*q_a-beta*q_b-tau*q_t)
   S_(alpha,beta,tau)(C,D,d).                       (7)
```

The leading `zeta^beta` is the fixed endpoint phase in THM-3479's refined
bank; the remaining three signs are the normalized inverse-transform signs.
Summing `(6)` over all chamber/drift buckets recovers the full guarded
endpoint product at each `(alpha,beta,tau)`.  Consequently the complete
`13^3` array has the frozen digest

```text
1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682. (8)
```

The normalization in `(7)` is forced by inversion on `F_13^3`.  Replacing
it with `13^(-2)` changes the first active coefficient from

```text
60583331018602792757  to  215330416995327424972.    (9)
```

Thus the normalization gate is not a support-only check.

## 3. Canonical guard-kernel route and the left-sheet phase

Independently define the translated kernel

```text
K_(C,D,d)(k)
 =sum_(s in F_13) g_C(s)g_D(s+d)zeta^(-ks).         (10)
```

For an actual atom pair `(a,C),(b,D)` with `d=b-a`, the substitution
`s=a+tau` gives the exact identity

```text
sum_tau g_C(a+tau)g_D(b+tau)zeta^(-k tau)
 =zeta^(ka) K_(C,D,d)(k).                           (11)
```

The audit checks `(11)` independently for all

```text
39^2 * 2 = 3042                                    (12)
```

atom-pair/frequency instances with `k in {0,1}`.  It then contracts

```text
B^ker_q(C,D,d)
 =13^(-3) sum_(alpha,beta)
   zeta^(beta-alpha*q_a-beta*q_b)
   sum_(b-a=d) A_(alpha,beta)(a,C)B_(alpha,beta)(b,D)
                  zeta^(q_t a)K_(C,D,d)(q_t).       (13)
```

The result is the entrywise identity

```text
B^dir_q(C,D,d)=B^ker_q(C,D,d)                       (14)
```

for all five classes and all `117=3*3*13` chamber/drift types.  In
particular, the phase `zeta^(q_t a)` in `(13)` has positive sign.  It is not
an optional origin convention once the direct route `(7)` is fixed.

Three explicit hostile witnesses distinguish the conventions:

```text
inverse target sign, q=(0,1,0), (L,L,0):
  correct 26461852298079162167,
  reversed 119509360393971210577;

tau sign, q=(1,0,1), (L,L,0):
  correct 309477751434222963513,
  reversed 13758465298279769122;

drift DFT sign, q=(0,0,0), frequency 1:
  correct 83824780274804594487,
  reversed 468665068434428720737.                  (15)
```

Hence neither complex-conjugating a character nor reversing the drift
transform silently reproduces the certified tensor.

## 4. Support, inversion, and rank minors

Both routes in `(14)` give exactly `65` zero middle-involving types and all
`52=2*2*13` boundary types nonzero for each of the five classes.  The
resulting `5 x 52` matrix has rank five.  In the canonical order

```text
(L,L,0),(L,L,1),...,(L,L,12),(L,R,0),...,           (16)
```

its first five columns already give the nonzero determinant

```text
215017650918161056500 mod P.                         (17)
```

Thus rank five is witnessed by one displayed minor rather than only by a
Gaussian-elimination counter.

At every fixed drift, the `5 x 4` class/chamber matrix has rank four.  The
same holds at all thirteen Fourier frequencies.  The audit exhibits a
nonzero `4 x 4` row minor for each slice; the two certificate-table digests
are

```text
drift:   7d9f7bc2091007a57467d5e1a9bf0f132aef9018a769fbaf5fe966a4f1f318c6,
Fourier: 2ef04676425e8dbebd1957c84e851a30426c17618c5fbc8b3d6008ea7abffe59. (18)
```

For every class, and separately for every chamber pair, the audit checks
the full inversion identity

```text
f(d)=13^(-1) sum_(k in F_13) fhat(k)zeta^(kd),
fhat(k)=sum_(d in F_13) f(d)zeta^(-kd).              (19)
```

Thus the drift and Fourier rank statements refer to inverse-equivalent
arrays with pinned signs and normalizers.

## 5. Independent tree and cut-cycle certificates

Use the THM-3515 graph with eight vertices and thirteen increasingly
oriented edges

```text
03 04 05 12 14 17 24 27 34 35 45 46 47.            (20)
```

For every address or transformed mode and every one of the `72` lawful
charts, assign the eight role responses as vertex potentials `p_v` and set

```text
w_(uv)=p_u-p_v.                                     (21)
```

The audit evaluates each `K4` tree polynomial by two independent routes:

1. explicit enumeration and summation of its sixteen spanning-tree
   monomials; and
2. a reduced `3 x 3` weighted-Laplacian determinant.

It also evaluates the full reduced `7 x 7` determinant directly.  At every
bank entry and chart, all three routes satisfy

```text
det L_red=(H-q5) Tau_L(w) Tau_R(w).                  (22)
```

The four pointwise, drift, Fourier, and chamber-Fourier chart digests are
exactly those in THM-3515.  More strongly, these independently reconstructed
coefficient arrays regenerate the complete candidate semantic digest

```text
2c9495fb8bcb731361ba331d9ca4b84a60f21551dc49b16e6519c1fc4f2e9f97. (23)
```

There is one transcript-only defect in the submitted candidate output.  Its
`pointwise_52x72_chart_bank` line serializes `71`, rather than `72`, zero
entries in `product_zeros_by_chart`.  The candidate code iterates over all
`72` charts, its label and digest assert the 72-chart object, and this audit's
direct reconstruction produces the missing seventy-second zero.  No tensor
coefficient, digest, rank, or nonvanishing claim changes; the submitted
candidate files are left untouched.

The signed incidence matrix of `(20)` has rank seven, while its cycle space
has dimension six.  The audit constructs three fundamental triangles in
each `K4`, obtaining a rank-six cycle basis.  Edge `46`, the forced bridge,
has coefficient zero in every basis cycle.  Since `(21)` is a gradient, the
pairing with every cycle vanishes.  This is checked at all chart instances:

```text
pointwise:          22,464 pairings,
drift marginal:      5,616 pairings,
drift Fourier:       5,616 pairings,
chamber Fourier:    22,464 pairings,
global:                432 pairings,
total:              56,592 pairings, all zero.       (24)
```

Equations `(22)` and `(24)` hold simultaneously.  The determinant is a
polynomial in a chosen coboundary representative; it does not descend to an
absolute `H^1` class.  The unique bridge contributes essentially to every
spanning tree while contributing to no cycle.

## 6. Hostiles and exact scope

Flat role potentials kill every edge and all four factors in `(22)`.
Replacing the `q5` potential by the `H` potential kills the forced bridge
and every full determinant, independently of the two tetrahedral factors.
Together with `(9)` and `(15)`, these controls distinguish constant,
bridge-killed, wrong-character, wrong-drift, and wrong-normalization
constructions.

The result remains an endpoint theorem.  The source-to-target contract is

```text
source:       five refined endpoint character responses;
target:       52 address tensors and their drift transforms;
map:          atomwise Cartesian pair contraction, then role-potential
              differences on the private-support graph;
preserved:    guard covariance, chamber/drift address, all five classes,
              all 72 charts, weighted-tree nonvanishing;
destroyed:    common Boolean base, collision root, owner chronology, word,
              source/arrival sheets, horizon, and physical pairing;
needed sidecar: one character-independent common-stalk ancestry relation;
first test:    recover the frozen H-q5 bridge before any tree factor.       (25)
```

Therefore this theorem proves no common ancestry, lawful Boolean or physical
current, grouped exact-address coefficient, all-unit projector, nonzero
absolute graph `H^1`, LRC `7 x 13` bispectrum theorem, scalar-row exclusion,
or LRC(14).

## 7. Exact audit package

Run from the repository root:

```text
python -B 04-computation/lrc_ufull_all_role_modewise_spectral_independent_audit_20260816.py
python -B -O 04-computation/lrc_ufull_all_role_modewise_spectral_independent_audit_20260816.py
```

Normal and optimized executions reproduce the stored output exactly.  The
LF-normalized script and output hashes are

```text
975b48c148b88a04b9fea07dcb6352cfaa98591f14eace1a38cd581a5e5876dc,
6eb3f28c07e99c352f1d88615d6beddeb82c264892a9f4735a420b75bf932c19, (26)
```

and the independent semantic digest is

```text
b745ab5d95936b9a134cfe8b6ffde032e72020d950670ab8b1a4643dd58c8be6. (27)
```
