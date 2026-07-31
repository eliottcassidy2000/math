---
id: THM-2751
title: "Frozen source-clock-one root-zero wing spectrum and positive quotient decoder"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the canonical
  rail-8 source_clock=1 fibre, the physical left
  wing is nonempty but its marked delayed coefficient is identically zero,
  while the normalized right-wing target-window shape Q has a ten-label
  91-unit spectrum.  The raw right-wing coefficient has nonunit scalar G.
  No linear wing decoder exists.  Nevertheless a positive coefficient-derived
  convolution transports the full fixed-clock source profile to the target
  profile modulo the uniform target-null line.  This is not the full
  unclocked Mayer--Vietoris decomposition, a physical packet action, an
  endpoint current, a row exclusion, or LRC(14).
source: root/fixed-clock-root-zero-wing-spectrum-2026-07-28
audit: >
  thm2751-final-referee-2026-07-28 (independent immutable carrier, prefix,
  direct-target, terminal-annihilation, resultant, convolution,
  raw-versus-normalized coefficient-ring scope, normal/-O/stored and hash
  replay: ACCEPT)
depends_on:
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
related:
  - THM-830-b3-deletion-deck-mirror-current-calculus
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2750-arm-blind-clutch-no-go-and-minimal-marked-leakage
script: 04-computation/lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.py
output: 05-knowledge/results/lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.out
script_sha256: 25cbed38026d61891173c687006250a69fe38aea56d67439406bd8bb60fa2552
output_sha256: 1f0facaa6184370ee38972e04c57557b79f5ab9e3cff93b52ec3d1ad5c20dfdb
hash_basis: LF-normalized bytes
---

# THM-2751 -- the frozen source-clock-one wing is coefficientally one-sided

> **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
> repairs the fixed-clock boundary paragraph of promoted THM-2749 and the
> numerical premise of the former THM-2751 provisional body.  It does not change the
> proved two-sided common-section theorem.

## 1. First failed implication: the legacy carrier omitted the source-one factor

Promoted THM-2749 defines

```text
A_t = E3 intersect F_(1,0,t),
M_t = A_t intersect T_tau^(-1) A_t.
```

The canonical constructor for `F_(1,0,t)` is
`two.source_present_section(..., source_clock=1, ...)`.  Its first operation is

```python
intersect_sorted(source_intervals, clock_comb[source_clock])
```

and therefore inserts the load-bearing source-one danger factor
`d_(c1,1)`.  By contrast, the historical comparator
`lrc14_semantic_root_zero_clutch_refinement_probe_20260728.py` calls
`restrict_e3_and_sheet`, whose signature has no source-clock argument and
whose body contains only the `E3` and four `U_(s,t)` cuts.  It never intersects
`clock_comb[1]`.

This is the entire constructor mismatch.  It is **not** a delayed-word
mismatch: the companion imports the hash-pinned historical comparator and
requires the historical `build_q3_pair_prefixes[ell][6][kappa]` to equal
THM-2749's `marked_prefixes(..., deepest_fork)[ell][kappa]` on all fourteen
`(ell,kappa)` cells before printing the audit line.  Their masses are zero for
`ell=0` and `206986279500` for both half-digits at every `ell=1,...,6`.

Consequently the historical amplitudes

```text
1812281403506324508838080,
1826551436254490256030720
```

belong to the clock-unrestricted carrier.  They cannot be subtracted from
THM-2749's source-one common amplitude.  The old wing gains `2` and `7` in the
THM-2751 reservation are therefore retracted premises, not candidate facts.

## 2. Correct fixed-clock common-coordinate decomposition

Fix the promoted rail `8`, `(source_clock,s)=(1,0)`, the actual terminal fork
`D^6 Q_(3,{1,2})`, both relative-present factors, and the source seam
`Sigma_-=(169,181)/182`.  In one source coordinate put

```text
A_t = rail_common intersect E3 intersect F_(1,0,t),
B_t = rail_common intersect T_tau^(-1)(E3 intersect F_(1,0,t)),
M_t = A_t intersect B_t,
L_t = A_t minus B_t,
R_t = B_t minus A_t.
```

The script constructs these five interval unions first, checks
`A_t=M_t disjoint-union L_t` and `B_t=M_t disjoint-union R_t` exactly, and
only then applies the common terminal/relative/seam functional.  At the
nonzero-window representative `t=3` it also pushes `B_t` and `M_t` to the
forward target coordinate and recomputes their target coefficients directly;
the exceptional `B_12` is checked the same way.  The interval
lengths below are exact numerator lengths on the canonical `T`-grid before
the clock-dependent relative cut; divide by `T=297836897838480` for Haar
measure.

| target labels | `|A|` | `|B|` | `|M|` | `|L|` | `|R|` |
|---|---:|---:|---:|---:|---:|
| `t=0,1,2` | 0 | 0 | 0 | 0 | 0 |
| each `t=3,...,11` | 13751337600 | 13808634840 | 6320326320 | 7431011280 | 7488308520 |
| `t=12` | 7404566400 | 7435418760 | 0 | 7404566400 | 7435418760 |

Thus both wings are genuine nonempty physical interval unions.  Their marked
coefficient functionals are radically asymmetric.

Let

```text
C = 339633525654239542165440,
G = C/119 = 2854063240791928925760,
W = z^3+...+z^11,
U = z^3+...+z^12,
Q = 2(z^3+...+z^11)+121 z^12.
```

Every nonzero seven-clock vector has shape `(0,c,c,c,c,c,c)`.  The exact full
target profiles are

```text
S(A)  = C W,
T(B)  = 121 G U,
S(M)  = T(M) = C W,
S(L)  = 0,
T(R)  = G Q.                                             (1)
```

In particular, for `t=3,...,11`, the right-wing amplitude is
`2G=5708126481583857851520`; at `t=12` it is
`121G=345341652135823400016960`.  The left wing is physically nonempty but
lies in the kernel of this delayed coefficient functional for every target
label.  The raw fixed-clock shear is one-sided:

```text
T(B)-S(A) = T(R), not T(R)-S(L) with two charged endpoints. (2)
```

The null left coefficient has a sharper support explanation.  For every
`t=3,...,11` and all seven delayed clocks, `L_t` survives the relative-present
and seam cuts with exact interval mass `26444880`, but its intersection with
both actual delayed terminal half-prefixes is empty.  At `t=12`, `L_t` is
already empty after the relative-present/seam cuts.  Thus “physical left
wing” refers to the preterminal interval carrier; the fully marked terminal
left packet is empty, rather than a nonempty signed current cancelling to
zero.

## 3. Root units and target characters

After content `26`, root normalization, and reduction modulo `Phi_7`, the
per-label profiles are:

| channel | labels | root | reduced constant | determinant |
|---|---|---:|---:|---:|
| `A` | `3..11` | source `12` | 9 | 1 |
| `B` | `3..12` | target `1` | 8 | 12 |
| `M_source` | `3..11` | source `12` | 9 | 1 |
| `M_target` | `3..11` | target `1` | 4 | 1 |
| `L` | none | source `12` | 0 | 0 |
| `R` | `3..11` | target `1` | 4 | 1 |
| `R` | `12` | target `1` | 8 | 12 |

Hence the common gain remains `4/9=12=-1`.  On common fixed-clock labels the
source-to-target ratio is `8/9=11`, but it cannot extend globally because
`B` has the extra nonzero label `12`.  Target rotations and reflections
preserve support cardinality, so no scalar/dihedral map sends the nine-label
source profile to the ten-label target profile.

The primitive cyclotomic norms are

```text
Norm(W)=1,
Norm(U)=1,
Norm(Q)=8492431042211308167354471.                       (3)
```

In particular `Norm(Q)=1 mod91`, so the **normalized target-window shape**
`Q` is a `91`-unit in every primitive target-character fibre.  This does not
make the raw polynomial `T(R)=GQ` a `91`-unit:

```text
gcd(G,91)=91,
gcd(G/26,13)=1,                  gcd(G/26,91)=7.          (3a)
```

Canonical content/root normalization therefore leaves an invertible scalar
over `F_13`, but it still leaves septimal content.  No raw or `C_91`-unit
right-wing current is asserted.

Therefore all twelve primitive target characters survive for `A`, `B`, both
common rows, and `R`; every target character of `L` vanishes.  Factoring the
raw amplitudes gives

```text
Res(Phi13,S(A)) = C^12,
Res(Phi13,T(B)) = (121G)^12,
Res(Phi13,S(M)) = Res(Phi13,T(M)) = C^12,
Res(Phi13,S(L)) = 0,
Res(Phi13,T(R)) = G^12 Norm(Q).                          (4)
```

The proposed cross-wing character bank is identically zero:

```text
sum_t chi(t) S(L_t) conjugate(T(R_t)) = 0               (5)
```

for the trivial and all twelve primitive characters.  This is not a
cyclotomic cancellation: the left coefficient is already zero labelwise.

## 4. Decoder boundary

There is no scalar, dihedral, convolutional, or positive **wing decoder**
from `L` to `R`: every linear operator sends the zero left-wing coefficient
profile to zero, while `R` is nonzero.

There are useful coefficient-algebra decoders after forgetting the wing
decomposition:

```text
W^(-1) = z^2+z^6+z^10,      W W^(-1)=delta_0+2N,
U^(-1) = z+z^4+z^7+z^10,    U U^(-1)=delta_0+3N,
```

where `N=1+z+...+z^12`.  With

```text
K=U W^(-1)=(3,3,2,2,2,3,2,2,2,3,2,2,2),
```

one has `W K=U+20N`, hence

```text
121 [S(A)*K] = 119 [T(B)]  modulo the uniform target-null line. (6)
```

This is a positive coefficient-section recombination, not a physical packet
action: the marked `t`-sections were not proved target-covariant as whole
Boolean packets, and it erases precisely the one-sided wing data.

The right wing alone also has a positive integer adjugate decoder because
`Norm(Q)` is nonzero.  The exact vector and convolution are printed by the
companion.  That fact preserves one-sided target activity but does not repair
the zero cross-wing pairing.

More precisely, the printed nonnegative decoder `K_Q` satisfies

```text
Q*K_Q = delta_Q delta_0 + cN,
delta_Q=943603449134589796372719=81 mod91.               (7)
```

The coefficient ring matters here.  In the integral augmentation quotient

```text
A_Z = Z[C_13]/<N>,
```

multiplication by `Q` has determinant and cokernel order `Norm(Q)>1`; hence
`Q` does **not** generate all of `A_Z`.  Identity `(7)` proves only
`delta_Q A_Z subset (Q)`.  After rationalization (or localization at
`delta_Q`) `Q` is invertible.  Modulo `91`, `delta_Q=81` and `81^(-1)=9`, and
the companion checks the stronger typed identity

```text
Q*(9 K_Q) = delta_0 + 10N  modulo 91.                  (8)
```

Thus the **shape `Q`**, by itself, generates the rational/localized and
mod-`91` `C_13` augmentation quotients coefficientwise.  This shape statement
does not pass through the physical scalar.  The raw right cofiber is `GQ`
with `gcd(G,91)=91`; after canonical content division its scalar is `G/26`,
which is a unit modulo `13` but remains divisible by `7`.  Therefore the
canonically normalized **physical coefficient class** generates the
`F_13[C_13]` augmentation quotient, but not the `F_7` or `C_91` quotient.  It
can synthesize a correction of the form `-7a` only in that explicitly
`F_13`-typed sense.  This is the strongest holotopy survivor.  Its connection
to the common-ancestry vertical-edge debts of THM-2542 and THM-2591 is
conditional on a new physical attachment theorem selecting one
common-ancestry semantic vertical edge.  Nothing here realizes `K_Q` as a
whole-packet action, retains an external arm or endpoint phase, or supplies
that attachment.

## 5. Fixed fibre versus the full unclocked cover

Every result above is on the frozen physical-present fibre
`source_clock=1`.  It must not be augmented into a claim about the unclocked
one-sided sheet.  The independently exact reflection
`07-reflections/lrc14-full-present-intersection-wing-rank-drop-hostile-audit-20260728.md`
(commit `bd53dc4c5`) constructs the full cover as

```text
M_full = disjoint-union_e M_e
```

with nonzero same-clock pieces at `e=1,2,3` and empty cross-clock pieces.
On that larger object the source wing is nonzero and the target wing loses a
unit only **after augmentation** in the present-clock coordinate.  That is a
different rank-drop theorem.  The reflection's cited full-unclocked companion
was not present when this candidate was written, so its larger exact table is
related hostile evidence rather than a dependency of THM-2751.

## 6. Scope and reproduction

The result is theorem-grade as a **repair/no-go**:

- it preserves the proved THM-2749 two-sided common clutch and its norm-one
  target window;
- it retracts only the legacy single-sheet comparison and the reserved
  THM-2751 numerical premise;
- it proves that the fixed-`e=1` defect is a coefficient-null physical left
  wing paired with a charged right wing, and identifies terminal-support
  annihilation as the exact reason the left coefficient vanishes;
- it supplies a positive coefficient-derived transporter on the target
  augmentation quotient, and a normalized `F_13` right-cofiber correction
  generator whose target-window shape `Q` is separately a rational/localized
  and mod-`91` unit, but not an integral unit;
- it does **not** supply a raw, canonically normalized `F_7`, or `C_91`-unit
  physical correction;
- it supplies no endpoint current, physical target action, owner/root
  provenance, row exclusion, or LRC(14) decrement.

Run from the repository root:

```powershell
python 04-computation/lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.py
python -O 04-computation/lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.py
```

The companion imports the promoted THM-2749 companion and its pinned proved
dependencies.  It also hash-verifies the historical comparator before loading
it solely to reconstruct and compare the fourteen terminal prefixes; that
comparator is evidence for the correction boundary, not a proved dependency
or truth source.

Both modes must byte-match
`05-knowledge/results/lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.out`.
No truth-bearing Python `assert` is used.

QED.
