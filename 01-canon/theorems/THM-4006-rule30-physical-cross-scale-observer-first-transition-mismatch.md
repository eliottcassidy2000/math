---
id: THM-4006
title: "Rule 30 physical cross-scale observer and first finite transition mismatch"
status: >
  FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the
  canonical highest-shell universe n=2^m+t, 0<=m<=9, 0<=t<2^m, the exact
  physical phase owner, four odd units, ordinary carry, three projective
  levels, routed chains, center, and next selected bank are reconstructed for
  every 1<=n<=1023. The finite depth-four/mod-16/projective observer has 1019
  fibres and four repeated states; the first repeat is the harmless pair
  128/132, while 943/951 is the unique and first target/transition mismatch.
  A phase-owned base chain resolves it. The first cross-scale mismatch before
  projective history is 20/574, and that history resolves it. General-phase
  routing must start at the actual d-bit phase-tail prefix a: replacing
  a+r*2^d by r*2^d first fails at n=6 and gives the wrong selected bank on
  676 of the 684 nonzero extensions in this universe. These are bounded
  finite statements, not an all-scale finite observer, minimality result,
  gap bound, center-trace theorem, or Rule 30 prize consequence.
source: root + rule30_cross_scale_first_mismatch + independent no-import audit, 2026-08-24
audit: >
  PASS (independent no-import hostile audit, 2026-08-24). The verifier rebuilds
  the Mealy right action, signalizer gaps, Rule 30 rows, physical phase
  sections, phase-tail routing, normalized units, ordinary carry, projective
  ratios, every quotient fibre, and chronological first-mismatch tests. It
  checks 54,313 gates, including 2,036 THM-3824 controls, and distinguishes
  the first state repeat from the first target mismatch. Primary and audit
  normal/optimized streams match their frozen outputs after LF normalization.
depends_on:
  - THM-3511-rule30-orbit-signalizer-gap-renormalization-and-shallow-portrait-hostile
  - THM-3516-rule30-marked-van-der-put-carry-and-power-section-bridge
  - THM-3824-rule30-fixed-division-tariff-and-physical-phase-separation
related:
  - THM-3512-rule30-van-der-put-haar-cocycle-and-profinite-automaton-boundary
  - THM-3804-rule30-all-period-amplitude-lattice-smith-law
script: 04-computation/rule30_physical_cross_scale_observer_thm4006.py
output: 05-knowledge/results/rule30_physical_cross_scale_observer_thm4006.out
script_sha256: dd2557860f92421121e808d879353429aee0dcb1029744e1dd8f215bff834cdd
output_sha256: 7e8d71b6b6dab9c4ec1a1aa9578da5478f1212ea1ad912db83e8aa61434238f2
semantic_sha256: 6cd1995718b50ab1225e0b6d4725af00357e08eb6ca4304fcaeacc6640dff09c
independent_audit_script: 04-computation/rule30_physical_cross_scale_observer_independent_audit_thm4006.py
independent_audit_output: 05-knowledge/results/rule30_physical_cross_scale_observer_independent_audit_thm4006.out
independent_audit_script_sha256: 5e71d754306ea51ec64f14b88f2808eaf5ffcc43d37e30aebdbdb6bced406288
independent_audit_output_sha256: 51d1936b605c457e5be296535a926c6eb9758e07139323c94527cb088b8b8eec
independent_audit_semantic_sha256: 47c79b0e6ab72cc6519b32602608b77ddb9da7b9043c4579b97d2e159cc2231b
independent_record_stream_sha256: d5046b3c1a7bfc8b921b7bcd6ece5368a24a3feb13de43d4811883769203f084
hash_basis: raw LF bytes
---

# THM-4006 -- physical cross-scale observer through `n=1023`

**FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Retain the
packed Rule 30 orbit and notation of THM-3511, THM-3516, and THM-3824. The
finite universe is

```text
n=2^m+t,             0<=m<=9,             0<=t<2^m.   (1)
```

Thus every integer `1<=n<=1023` occurs exactly once. The result below is an
exact census in this universe, not an extrapolation.

## 1. Physical records and their observers

At scale `m`, put `q=2^m`, let `v_m` be the innovation valuation, and put

```text
d_m=v_(m+1)-v_m,
U_m(t)=(R_(t+q)-R_t)/2^v_m.                            (2)
```

Every unit used below is checked integral and odd. For each physical phase
`t`, the record retains:

1. the exact gap `d_m`;
2. the phase owner `s_(m,t)`, its selected rays `(0,1,2)`, and its complete
   depth-four portrait;
3. the four consecutive normalized units
   `U_m(t+jq)`, `0<=j<4`, both exactly and modulo `16`;
4. THM-3516's ordinary center-carry decomposition, retaining the direct
   target parity and the lower-residue carry parity separately;
5. the last three projective scale entries

   ```text
   (d, G mod 16, Z mod 16),
   G=-U_m(t+q)/U_m(t),          Z=(1-G)/2^d,            (3)
   ```

   together with the exact reduced rational pairs; and
6. the physical center `c_n` and the next phase-owner bank on rays
   `(0,1,2)`.

For routed chains, let `a` be the first `d=d_m` bits of the actual phase
tail. The three physical inputs are

```text
x_r=a+r*2^d,                 r in {0,1,2}.              (4)
```

Apply the owner once and twice to each `x_r`. Stripping the fixed `d` bits
from the second images exactly recovers the next selected owner bank. This is
the finite query form of THM-3511's lawful phase transition

```text
(s,eta) -> ((s^2)|_a, shift^d eta).                    (5)
```

The **finite projective observer** is the tuple consisting of `d_m`, the
selected current bank, full depth-four portrait, four units modulo `16`,
ordinary carry state, and three projective entries in `(3)`. The **exact
faithful state** instead retains the exact owner word, exact units, all carry
terms, and exact projective rational pairs. The exact faithful state has
`1023` distinct values. This is not a finite-state theorem: its owner words
and integers are unbounded coordinate types.

## 2. Complete refinement census

For each observer level, the exact row below records

```text
(number of fibres, nontrivial fibres, target-mismatch fibres,
 same-scale mismatch fibres, cross-scale mismatch fibres).          (6)
```

The target is the pair `(center,next selected bank)`. The monotone refinement
is

```text
branch       (4,    4,   4,   4,  2)
bank         (246,185, 185, 176, 58)
portrait     (395,229, 226, 224, 15)
odd          (753,190, 183, 181,  2)
carry        (870,128, 109, 107,  2)
projective  (1019,  4,   1,   1,  0)
base        (1022,  1,   0,   0,  0)
one         (1023,  0,   0,   0,  0)
two         (1023,  0,   0,   0,  0)
exact       (1023,  0,   0,   0,  0).                 (7)
```

Here `base` adds the phase-owned chain at ray zero, while `one` and `two`
add the other routed chains. Thus the base chain removes every target
mismatch in this universe; either one additional off-ray chain makes the
whole bounded record injective. No minimality or all-scale sufficiency is
claimed.

## 3. First repeated state and first target mismatch

The four repeated finite projective states are exactly

```text
(128,132),        (256,260),        (515,519),
(943,951).                                                (8)
```

The first three are harmless: both records in each pair have the same center
and next bank. Therefore `128/132`, not `943/951`, is the first state
collision. The last pair is the unique target mismatch and hence the first
transition failure:

```text
(n,m,t)=(943,9,431), (951,9,439),
gap=2,
selected current bank=(15,6,1),
four-unit shadow=(5,7,5,15),
carry state=(visible,target 1,carry 0),
projective shadow=((8,1,5),(1,3,15),(2,5,11)),
centers=1,1,
next banks=(7,10,13),(11,14,1).                        (9)
```

Their phase-owned base chains begin

```text
(input,first image,second image)=(1,54,29), (1,22,45), (10)
```

and therefore separate the transition. The exact THM-3824 free defects at
phases `431` and `439` are also distinct. The mismatch in `(9)` is a
finite truncation collision, not an exact same-scale Smith collision.

## 4. First cross-scale mismatch and its repair

After gap, bank, portrait, odd shadow, and ordinary carry are retained, the
first cross-scale target mismatch is

```text
(n,m,t)=(20,4,4), (574,9,62),
gap=2,
selected current bank=(1,12,11),
four-unit shadow=(3,1,3,9),
carry state=(visible,target 0,carry 0),
centers=0,0,
next banks=(5,12,7),(5,8,11).                          (11)
```

Their three-level projective histories differ, so the projective sidecar
repairs this first cross-scale hostile. No cross-scale target mismatch
survives the complete finite projective observer through `n=1023`.

As a same-scale control, the companion independently checks all `2036`
strict inequalities

```text
Delta_m(t)=U_m(t+2^m)-U_m(t)>0,
Delta_m(t+1)>Delta_m(t),                               (12)
```

in this universe. THM-3824 already proves `(12)` at every scale; the finite
replay guards the translation into the present records.

## 5. The phase-owned base is load-bearing

Replacing the physical inputs in `(4)` by the marked-origin zero-base inputs
`r*2^d` is false. The first failure is

```text
(n,m,t,a)=(6,2,2,1),
naive selected bank=(11,12,9),
physical selected bank=(15,8,5),
physical chains=(1,16,61),(5,4,33),(9,24,21).          (13)
```

This correction affects routed chain coordinates and the next bank recovered
from them. It does not change the independently defined odd, projective,
carry, or center coordinates. The fact that `n=6` is also the first record
with ordinary carry one is an unrelated coincidence.

There are `684` records with nonzero phase-tail prefix in this universe. The
zero-base selected bank is wrong on `676` of them. The eight accidental
selected-bank equalities occur at

```text
307,319,371,387,433,593,756,955.                       (14)
```

Thus neither “the zero base always works” nor its converse “every nonzero
extension changes this selected bank” is valid. The lawful basepoint in
`(4)--(5)` is the needed coordinate.

## 6. Scope and reproduction

The result stops at scale nine, depth-four owner portraits, four odd units,
four retained arithmetic bits, three projective levels, and the three
selected rays. It does not prove:

1. any statement for `n>=1024`;
2. an all-scale finite observer or minimal symbolic/Mealy representation;
3. bounded gaps or eventual repetition of exact states;
4. periodicity, nonperiodicity, balance, or complexity of the center trace;
5. any of the Rule 30 prizes.

Reproduce from the repository root:

```text
python -B 04-computation/rule30_physical_cross_scale_observer_thm4006.py
python -B -O 04-computation/rule30_physical_cross_scale_observer_thm4006.py
python -B 04-computation/rule30_physical_cross_scale_observer_independent_audit_thm4006.py
python -B -O 04-computation/rule30_physical_cross_scale_observer_independent_audit_thm4006.py
```

The independent audit passes `54,313` exact gates. Its semantic SHA-256 is
`47c79b0e6ab72cc6519b32602608b77ddb9da7b9043c4579b97d2e159cc2231b`,
and its exact record-stream SHA-256 is
`d5046b3c1a7bfc8b921b7bcd6ece5368a24a3feb13de43d4811883769203f084`.
