---
id: THM-4423
title: "LRC14 height-499 atlas and low-carrier shell ladder"
status: >
  FINITE-EXACT THROUGH HEIGHT 499 + PROVED ALL-HEIGHT N<=112 COROLLARY
  RELATIVE TO THM-4414/4418/4422/4425 + VERIFIED-EXACT. All 753,853 eligible
  triples through height 499 satisfy the degree-zero raw projection target,
  uniquely sharply at (1,5,11). Of the 4,599 nonautomatic rows, 4,598 are
  signed norm-four rank-one rows and the only rank-two row is the A2 hexagon
  at (19,23,29). A 38-shell continuation classifies all 482
  nonautomatic rows with even carrier count 92 through 112 as signed
  norm-four rows, already closed analytically by THM-4422. Any remaining
  hostile has at least 114 carriers, largest speed at least 503, is genuinely
  rank two, and lies below the THM-4418 ratio tail. The universal inequality,
  entry, synchronization, and LRC(14) remain open.
source: root + network_universal + projection_inequality + dense_geometry_dichotomy / LRC14 continuation session, 2026-09-05
depends_on:
  - THM-4414-lrc14-six-separated-contact-capacity-collapse
  - THM-4418-lrc14-sharp-pair-arithmetic-and-forty-four-thirteenths-tail
  - THM-4422-lrc14-projection-deficit-and-beatty-row-reduction
  - THM-4425-lrc14-all-height-rank-one-carrier-closure
related:
  - THM-4420-lrc14-near-doubling-ray-network-closure
script: 04-computation/lrc14_height499_projection_atlas_thm4423.py
output: 05-knowledge/results/lrc14_height499_projection_atlas_thm4423.out
script_sha256: 981911d28c223d588eace6e262e24e9874432980b9bec3fd74773bef7c0e65e0
output_sha256: f47c6c89a4fbadf9df66b2d169855175bfda53b26ba2e0194904ba99a9f63a18
rows_sha256: a6d482622d2c96a9981239c66aa40aeb21cd570e46fd1cd6d5239d831eb72144
ladder_script: 04-computation/lrc14_low_carrier_shell_ladder_thm4423.py
ladder_output: 05-knowledge/results/lrc14_low_carrier_shell_ladder_thm4423.out
ladder_script_sha256: ed9ea158a7aa8ff6014955ec7dbdec69122013fe5b02c017abc31f1960ad3d00
ladder_output_sha256: 63932f2a166c2c2731de0df645d2c41a53db505c398f8d8e0f1e08161dd4db4f
ladder_rows_sha256: c66c18cdf85cc7128703d7f758f383fa6013ff4feb30feda624526779b73628e
dense_classifier_source: 04-computation/lrc14_height499_dense_carrier_geometry_classifier.cpp
dense_classifier_output: 05-knowledge/results/lrc14_height499_dense_carrier_geometry_classifier.out
dense_classifier_source_sha256: 3e7882e9df72b6bc3429af7ebbde7ef10d2aae1fe11d9516615436bab136ccde
dense_classifier_output_sha256: 107f45ad118c1c8d16b2c30161b1d096bc98b26e79544c1d5f131f0e5e990891
hash_basis: raw LF bytes
audit: >
  PASS. The integer-only verifier performs 47,467,792 explicit gates over
  753,853 triples. A definition-level relation-lattice enumeration independently
  reconstructs all 2,910 rows through height 79 plus eight outer controls.
  The canonical deterministic full replay and the frozen output byte-match.
  A separate exact C++ classifier reconstructs every one of the 4,599 dense
  rows from the relation-lattice box before classifying its directions; full
  height-499 GCC -O0 and -O3 outputs byte-match the frozen LF artifact. The
  38-shell compiler likewise reconstructs all 482 retained rows by both
  kernels, freezes their full semantic digest, and performs 30,402,720
  explicit gates.
---

# THM-4423 -- LRC14 height-499 atlas and low-carrier shell ladder

**FINITE-EXACT THROUGH HEIGHT 499 + PROVED ALL-HEIGHT `N<=112` COROLLARY
RELATIVE TO THM-4414/4418/4422/4425 + VERIFIED-EXACT.** This is a local three-speed
degree-zero theorem. It does not prove arbitrary chart entry, synchronization
with eleven other speeds, or `LRC(14)`, which remains **OPEN**.

## 1. Exact bounded atlas

For every sorted primitive distinct positive triple

```text
w=(a,b,c),       1<=a<b<c<=499,                         (1)
```

with all speeds odd and nonzero modulo three, let `Lambda(w)` and the three
raw projections `S_i(w)` be those of THM-4414/4422. Then

```text
min_i S_i(w)<=6/77.                                     (2)
```

The atlas contains exactly

```text
eligible speed values                                  167
primitive triples                                  753,853
hostiles                                                  0
equality rows                                               1.             (3)
```

The unique equality row is `(1,5,11)`, where all three projections are
`6/77`. The strongest strict row with `c>79` is

```text
w=(1,67,133),       min_i S_i=60/931.                   (4)
```

The automatic threshold from THM-4422 divides this finite atlas much more
rigidly than the projection values alone reveal. Exactly 4,599 rows have
`11N>2c`. An independent two-definition carrier classifier finds

```text
signed norm-four and rank one                         4,598
other rank-one rows                                       0
rank-two rows                                             1.
```

The unique rank-two dense row is

```text
w=(19,23,29),       N=6,
Lambda=+/-{(1,8,-7),(10,-7,-1),(11,1,-8)},
(11,1,-8)=(1,8,-7)+(10,-7,-1).
```

Thus the first dense non-ray object is exactly an `A_2` additive hexagon.
This classification is **FINITE-EXACT THROUGH HEIGHT 499**; it is not an
all-height rank dichotomy.

Across the full atlas the carrier count ranges from zero to 142. The numbers
of triples on which the three labelled projections individually succeed are

```text
750,105,       751,734,       751,687,                  (5)
```

while the first minimizing projection occurs respectively

```text
15,931,         76,291,       661,631                   (6)
```

times. These labels are diagnostic only; `(2)` uses their minimum.

## 2. Exact compiler and hostile controls

In raw `y=3x` coordinates, the script sweeps the exact intervals

```text
|w_i y-N_i|<3/14,
```

retains precisely the owner triples whose sheet labels are a permutation of
`0,1,2`, and forms `C=w cross N`. No floating-point comparison occurs. For one
carrier the three THM-4414 summands have common denominator

```text
D=14abc                                                (7)
```

and integer numerators

```text
min(c[3(a+b)-14|C_3|],6ab),
min(b[3(a+c)-14|C_2|],6ab),
min(a[3(b+c)-14|C_1|],6ab).                            (8)
```

The program sums `(8)`, selects the least projection without consulting the
physical three-facet minimum, and compares `77 numerator` with `6D`.

The interval sweep is checked by a separate kernel-lattice box enumeration on
every triple through height 79 and on eight fixed outer controls, for 2,918
definition-level reconstructions. It reproduces the canonical height-79 pair
counts and equality row. The full ordered semantic table has digest

```text
a6d482622d2c96a9981239c66aa40aeb21cd570e46fd1cd6d5239d831eb72144. (9)
```

The false unweighted-average shortcut is retained as a hostile:

```text
w=(1,5,7),
(S_1,S_2,S_3)=(4/35,6/49,8/245),
S_1+S_2+S_3=66/245>18/77.                              (10)
```

Thus the bounded result is a selected-projection certificate, not an average
or a physical-mass computation.

## 3. Finite-to-global carrier-count splice

Let `N=|Lambda(w)|`. THM-4422 proves

```text
N<=2c/11  implies  min_i S_i<=6/77.                    (11)
```

The live carrier set is centrally symmetric under `C -> -C`, and the zero
carrier is excluded by the residue gate. Thus `N` is even. Every row with

```text
N<=90                                                   (12)
```

satisfies the target at every height. If `c<=499`, this is `(2)`. If `c>499`,
then oddness and the ternary-unit condition make the first possible value
`c=503`; hence

```text
2c/11>=1006/11>91>90,                                  (13)
```

and `(11)` applies.

For a fixed even count `N_0>90`, criterion `(11)` is automatic at every
height `c>=11N_0/2`. Therefore closing a bounded count band requires only
the eligible height shells strictly between 499 and that linear threshold.
One consolidated continuation scans all 38 such shells needed for
`N_0<=112`. It retains precisely the rows with

```text
92<=N<=112,       N even,       11N>2c.                (14)
```

Every retained carrier dictionary is independently reconstructed by both
the interval sweep and the direct relation-lattice box. The exact count
ledger is

| `N` | `92` | `94` | `96` | `98` | `100` | `102` | `104` | `106` | `108` | `110` | `112` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| rows | 6 | 8 | 16 | 22 | 56 | 18 | 77 | 51 | 56 | 60 | 112 |

There are 482 rows in total. Every one lies in one of THM-4422's three
signed norm-four relation families; there is no `other` row. Hence their
projection target follows analytically from THM-4422. The direct exact sums
also satisfy the target and provide a hostile control on the classifier.
Together with `(12)--(13)`, this proves the all-height band

```text
N<=112.                                                 (15)
```

The next generated task is exact: `N=114` has 42 nonautomatic eligible
shells, beginning at `(N,c)=(114,503)`. The first eligible height where
`N=114` becomes automatic is `c=629`; the real threshold `627` is excluded
modulo three.

## 4. Consolidated residual frontier

THM-4418 independently closes every triple with

```text
c/b>=44/13.                                             (16)
```

THM-4425 independently closes every rank-one live carrier set. Combining
that theorem, `(15)`, and `(16)`, any counterexample to the remaining THM-4414
projection inequality must lie in

```text
c>=503,       N>=114,       c/b<44/13,                 (17)
```

and must contain at least two nonparallel live carrier directions. This is a
precise residual region for the local certificate, not a counterexample
region for the full lonely-runner conjecture.

Run the full deterministic replay with

```powershell
python -B 04-computation/lrc14_height499_projection_atlas_thm4423.py --height 499
python -B 04-computation/lrc14_low_carrier_shell_ladder_thm4423.py --workers 4
g++ -std=c++20 -O3 04-computation/lrc14_height499_dense_carrier_geometry_classifier.cpp -o .scratch/lrc14_h499_dense.exe
.scratch/lrc14_h499_dense.exe 499 --audit-dense
```

The exact atlas is intentionally expensive; its 47 million gates take on the
order of tens of minutes on a desktop machine. `--height 79` gives the fast
canonical regression.
