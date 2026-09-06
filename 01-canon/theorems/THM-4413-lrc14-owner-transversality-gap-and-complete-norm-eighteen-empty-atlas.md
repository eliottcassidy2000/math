---
id: THM-4413
title: "LRC14 owner transversality gap and complete norm-eighteen empty atlas"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4386/4393 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. For every odd ternary-unit three-speed comb, the
  distinct-owner gate excludes roof tangency and gives a sharp height-sensitive
  positive component floor. In the complete minimal norm-18 shell, exactly
  three physical combs are empty. This is local three-speed structure, not
  arbitrary entry, synchronization, or LRC(14).
source: root + gap_empty_referee / LRC14 continuation session, 2026-09-05
depends_on:
  - THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence
  - THM-4393-lrc14-minimal-ternary-unit-norm-eighteen-shell
related:
  - THM-4394-lrc14-minimal-ternary-unit-norm-twenty-shell
  - THM-4398-lrc14-one-zero-relation-residue-dichotomy-and-small-norm-atlas
  - THM-4402-lrc14-norm-eighteen-vanishing-live-carrier-gap
primary_script: 04-computation/lrc14_owner_transversality_norm18_empty_atlas_thm4413.py
primary_output: 05-knowledge/results/lrc14_owner_transversality_norm18_empty_atlas_thm4413.out
primary_script_sha256: 41323519efc4570e6a764e3f0dd7b3a736203a9f300d798e7f0c447c740355cd
primary_output_sha256: 39c090511245370530cc13b2f4ae2f74eaabc1e8ce65e1116882524921ac1af0
independent_audit_script: 04-computation/lrc14_owner_transversality_norm18_empty_atlas_thm4413_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_owner_transversality_norm18_empty_atlas_thm4413_independent_audit.out
independent_audit_script_sha256: 25430c2332858c70962117f0c3abc9eef59f9d95c7df5a9a092b53e6b5a5c285
independent_audit_output_sha256: 326549ec026bf509115ff4cb10161fcfd0a26e38c4ecb0c4d0452498a3e75956
hash_basis: raw LF bytes
audit: >
  PASS. The primary independently intersects the rank-two carrier lattice and
  sweeps every rational nearest-integer event on the physical circle. The
  second implementation instead builds exact affine defect fibres from
  Bezout sections and checks them against a separate literal interval sweep.
  They agree on all 180 minimal norm-18 triples in the complete height-73
  core. Normal, optimized, and fixed-seed outputs match; all gates remain live.
---

# THM-4413 -- LRC14 owner transversality gap and complete norm-eighteen empty atlas

**PROVED ELEMENTARY RELATIVE TO THM-4386/4393 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.** This theorem concerns one local scale-three comb on
three speeds. It proves no arbitrary nonresonance, chart entry, seam passage,
or synchronization with the other eleven speeds. `LRC(14)` remains **OPEN**.

## 1. Universal arithmetic boundary gap

Let

```text
w=(w1,w2,w3),       1<=w1<w2<w3,
```

be distinct positive odd integers, each nonzero modulo three. Let
`C in Z^3` be a raw carrier with

```text
C dot w=0,                 3 does not divide C1*C2*C3.       (1)
```

By THM-4386, the second condition is exactly the distinct-owner gate. At
radius `r=3/14`, define the pair margins and obstruction gap

```text
p_i(w,C)=3(w_j+w_k)-14|C_i|,       {i,j,k}={1,2,3},
D_w(C)=-min_i p_i(w,C).                                  (2)
```

The exact component roof is

```text
L_w(C)=max(0,min(
  3/(7w1), 3/(7w2), 3/(7w3),
  p_1/(14w2w3), p_2/(14w1w3), p_3/(14w1w2))).            (3)
```

Every `p_i` is even because `w_j+w_k` is even. Moreover

```text
p_i = -14|C_i| = |C_i|  (mod 3),                          (4)
```

so owner admissibility makes every margin a nonmultiple of three. Thus

```text
p_i in 2Z minus 6Z.                                       (5)
```

In particular the owner lattice never meets a roof tangent. Since the cap
terms in `(3)` are positive, the exact sharp dichotomy is

```text
L_w(C)>0  iff every p_i>0,

live owner carrier: D_w(C)<=-2,
dead owner carrier: D_w(C)>=+2.                           (6)
```

Both constants occur in the minimal norm-eighteen shell:

```text
live: w=(29,73,77), c=(7,-7,4), C=(32,1,-13),
      p=(2,304,124), D=-2;

dead: w=(7,25,29), c=(13,1,-4), C=(4,7,-7),
      p=(106,10,-2), D=+2.                                (7)
```

The hypotheses are load-bearing. If owner admissibility is dropped, then

```text
w=(1,5,13), C=(-2,3,-1), p=(26,0,4)                       (8)
```

is tangent at the owner-forbidden coordinate `C_2=3`. If oddness is dropped,

```text
w=(1,5,14), C=(-4,-2,1), p=(1,17,4)                       (9)
```

is owner-admissible with `L=1/980`; the parity quantum disappears.

## 2. Sharp height-sensitive component floor

For a live carrier put `p_*=min_i p_i`. Comparing every term of `(3)` over
the largest pair denominator `14w2w3` gives

```text
L_w(C) >= min(p_*,6w2)/(14w2w3)
       >= 1/(7w2w3).                                     (10)
```

The live witness in `(7)` has its first margin equal to two and satisfies

```text
L_(29,73,77)(32,1,-13)=1/(7*73*77)=1/39347.              (11)
```

Hence the second bound in `(10)` is sharp, even within the minimal
norm-eighteen shell. It depends on height and therefore does not contradict
THM-4402: no positive component quantum depending only on the discrete
relation pattern exists.

## 3. All-height reduction for minimal norm eighteen

Now suppose that a primitive full-support ternary-unit relation

```text
c dot w=0,             ||c||_1=18                         (12)
```

is shortest among such relations. Its magnitude pattern is one of

```text
(1,1,16), (1,4,13), (1,7,10),
(2,5,11), (4,7,7), (5,5,8).                              (13)
```

The nearest zero-defect owner carriers are `C=+c` and `C=-c`. If they are
dead at coordinate `i`, then `(6)` and the relation equation give

```text
3(w_j+w_k)<14|c_i|,
|c_i|w_i<=max(|c_j|,|c_k|)(w_j+w_k).                     (14)
```

Thus every speed is below `14M/3`, where `M=max_i|c_i|`. The largest odd
three-unit heights below the resulting strict bounds are

| pattern | real bound `14M/3` | analytic cutoff |
|---|---:|---:|
| `(1,1,16)` | `224/3` | 73 |
| `(1,4,13)` | `182/3` | 59 |
| `(1,7,10)` | `140/3` | 43 |
| `(2,5,11)` | `154/3` | 49 |
| `(4,7,7)` | `98/3` | 31 |
| `(5,5,8)` | `112/3` | 37 |

If the physical comb were empty, every available norm-eighteen chart would
have dead central carriers, so this table is an all-height completeness
reduction, not an empirical search limit.

## 4. Complete exact finite core

The exact universe

```text
w3<=73, sorted distinct positive odd ternary-unit speeds, gcd(w)=1           (15)
```

contains `2,289` triples. Direct coefficient-box generation gives `232`
norm-eighteen relation incidences on `209` triples. Removing every triple with
one of the twenty smaller even ternary-unit patterns leaves

```text
190 relation incidences on 180 physical triples.                           (16)
```

The primary rank-two-lattice dictionary and literal rational-circle sweep
agree carrier-by-carrier on all 180 triples, comprising 914 positive grouped
components. The patternwise census inside the analytic cutoffs is

| pattern | minimal triples | relation incidences | central-dead incidences |
|---|---:|---:|---:|
| `(1,1,16)` | 14 | 14 | 2 |
| `(1,4,13)` | 16 | 17 | 4 |
| `(1,7,10)` | 8 | 9 | 1 |
| `(2,5,11)` | 12 | 13 | 2 |
| `(4,7,7)` | 1 | 1 | 1 |
| `(5,5,8)` | 5 | 5 | 0 |

## 5. Central-dead and physical-empty atlases

For a chart `c`, every live carrier has intrinsic defect

```text
c cross C=delta*w,             delta in {-3,0,3}.          (17)
```

Exact minimization on each affine defect line gives the complete list of
central-dead presentations:

| `w` | `c` | pattern | `D_0` | best owner `D_3` | live counts `(-3,0,3)` | mass |
|---|---|---|---:|---:|---:|---:|
| `(1,17,37)` | `(2,-11,5)` | `(2,5,11)` | 40 | -8 | `(1,0,1)` | `8/4403` |
| `(1,17,55)` | `(1,-13,4)` | `(1,4,13)` | 14 | -34 | `(1,0,1)` | `2/385` |
| `(1,19,35)` | `(16,1,-1)` | `(1,1,16)` | 62 | -4 | `(1,0,1)` | `64/4655` |
| `(1,25,41)` | `(16,1,-1)` | `(1,1,16)` | 26 | -22 | `(2,0,2)` | `172/7175` |
| `(5,13,41)` | `(1,-13,4)` | `(1,4,13)` | 44 | -22 | `(1,0,1)` | `22/3731` |
| `(7,11,43)` | `(5,-11,2)` | `(2,5,11)` | 4 | 34 | `(0,0,0)` | 0 |
| `(7,11,47)` | `(13,-4,-1)` | `(1,4,13)` | 8 | 20 | `(0,0,0)` | 0 |
| `(7,25,29)` | `(4,7,-7)` | `(4,7,7)` | 2 | 20 | `(0,0,0)` | 0 |
| `(7,25,29)` | `(13,1,-4)` | `(1,4,13)` | 20 | 2 | `(0,0,0)` | 0 |
| `(13,23,31)` | `(1,-10,7)` | `(1,7,10)` | 8 | -22 | `(1,0,1)` | `22/4991` |

Defect `+/-3` rescues six of the ten presentations. The four remaining rows
represent exactly three physical triples because `(7,25,29)` has two charts.
Consequently the complete minimal norm-eighteen shell has exactly the empty
combs

```text
{7,11,43}, {7,11,47}, {7,25,29}.                         (18)
```

The first empty height is 29. The sharp global empty height ceiling is 47,
attained by `(7,11,47)`. The sharp central-dead presentation ceiling is 55,
attained by `(1,17,55)`; its comb is rescued off the central fibre.

## 6. Why those three combs are empty

If the owner gate is removed, each triple in `(18)` has eleven positive
geometric lattice carriers, including `C=0`. Each also has a one-zero
complement `q` satisfying `c cross q=w`:

| `w` | `q` | positive deleted carrier | roof length |
|---|---|---|---:|
| `(7,11,43)` | `(3,2,-1)` | `3q=(9,6,-3)` | `18/3311` |
| `(7,11,47)` | `(2,3,-1)` | `3q=(6,9,-3)` | `18/2303` |
| `(7,25,29)` | `(-3,2,-1)` | `3q=(-9,6,-3)` | `18/5075` |

The displayed carrier occupies the defect-three real roof, but every
coordinate is zero modulo three. It is deleted solely by the distinct-owner
gate. All neighboring owner-admissible lattice samples lie on the dead side
of the two-unit gap `(6)`. Thus these are residue-sampling failures of a
nonempty real roof, not geometric emptiness before discretization.

## 7. Reproduction and frontier

Run from the repository root:

```powershell
python -B 04-computation/lrc14_owner_transversality_norm18_empty_atlas_thm4413.py
python -B -O 04-computation/lrc14_owner_transversality_norm18_empty_atlas_thm4413.py
python -B 04-computation/lrc14_owner_transversality_norm18_empty_atlas_thm4413_independent_audit.py
python -B -O 04-computation/lrc14_owner_transversality_norm18_empty_atlas_thm4413_independent_audit.py
$env:PYTHONHASHSEED='4408'
python -B -O 04-computation/lrc14_owner_transversality_norm18_empty_atlas_thm4413_independent_audit.py
```

The primary has 6,993 explicit gates; the broader affine-fibre audit has
1,320,358. Neither uses `assert`. The next useful consumer is a
height-sensitive count-versus-floor estimate or the exact component-network
projection; `(10)` alone is a lower bound and does not close the universal
three-speed `6/77` target. Arbitrary entry and `LRC(14)` remain **OPEN**.
