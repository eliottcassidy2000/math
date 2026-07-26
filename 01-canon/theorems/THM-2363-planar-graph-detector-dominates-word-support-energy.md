---
id: THM-2363
title: "Planar graph detector dominates word-support energy"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. For any
  joint coefficient array A(q,z) on F_169^2, let A_0(q)=sum_z A(q,z)
  and let D_graph be the sum of THM-2356's 169 planar-graph detectors.
  Then D_graph is the nonzero-target joint l2 energy plus a nonnegative
  graph-sum term, and
  D_graph >= (1/169)sum_(q!=0)|A_0(q)|^2
          >= E_sigma/169
  for each of THM-2337's three word-support masks. The constant 1/169
  is sharp for every mask, with a complete equality mechanism.
  D_graph can be positive when every A_0(q) vanishes, so it strictly
  refines the old aggregate energies. This is a coefficient-level
  comparison; it does not realize the missing graph-channel phase
  ratio, decrement an LRC profile, or prove LRC(14).
source: codex-2026-07-25-graph-detector-word-energy
depends_on:
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing
related:
  - THM-2340-owner-word-anova-target-landing
script: 04-computation/planar_graph_detector_word_energy_thm2363.py
output: 05-knowledge/results/planar_graph_detector_word_energy_thm2363.out
script_sha256: 4475c977bf066ace6223944494a89724123b3e784d091f3a2c9e3d3df4602db2
output_sha256: 2e55b91766379d33c70cdd4b2c2ad3c41b48d1f1bb3ebd4fdec6a66293ec8b43
hash_basis: working-tree bytes (LF)
---

# THM-2363 -- every word-support energy pays the graph detector

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2337's word-support energies and THM-2356's planar-graph detector
look like different remaining debts. They are in fact ordered positive
quadratic forms. The graph detector dominates every word-support energy,
with a sharp constant.

## 1. Joint and graph coordinates

Fix the chosen identifications used in THM-2356 and write

```text
K=F_169,                    N=|K|=169.
```

Let

```text
A(q,z) in C,                (q,z) in K x K,          (1)

A_0(q)=sum_(z in K)A(q,z).                           (2)
```

The aggregate (2) is exactly the uncoloured target coefficient denoted
`A(q)` in THM-2337. For each `c in K`, put

```text
Z_c(q)=A(q,q^2/2+c).                                (3)
```

For a fixed `q`, the map

```text
c -> z=q^2/2+c
```

is a bijection of `K`. Hence the `169` graphs in (3) partition all
`169^2` joint target--jet cells.

THM-2356 defines the coefficient-level graph detector

```text
D_c
 =sum_(q!=0)|Z_c(q)|^2
  +|sum_(q!=0)Z_c(q)|^2.                            (4)
```

Put

```text
D_graph=sum_(c in K)D_c.                            (5)
```

## 2. Exact identity and domination

Summing (4) and using the graph partition gives

```text
D_graph
 =sum_(q!=0,z)|A(q,z)|^2
  +sum_c |sum_(q!=0)A(q,q^2/2+c)|^2.               (6)
```

The second term is nonnegative. For each fixed nonzero `q`,
Cauchy--Schwarz gives

```text
|A_0(q)|^2
 <=N sum_z |A(q,z)|^2.                              (7)
```

Therefore

```text
D_graph
 >=sum_(q!=0,z)|A(q,z)|^2
 >=1/N sum_(q!=0)|A_0(q)|^2.                        (8)
```

Let `P_sigma`, for

```text
sigma in {{a},{b},{a,b}},
```

be any of THM-2337's exact word-support masks. Each mask takes values in
`{0,1}` and vanishes at `q=0`. Its landing energy is

```text
E_sigma=sum_q P_sigma(q)|A_0(q)|^2.                 (9)
```

Combining (8)--(9) proves the sharp common comparison

```text
D_graph
 >=1/169 sum_(q!=0)|A_0(q)|^2
 >=E_sigma/169.                                    (10)
```

In particular,

```text
E_sigma>0  =>  D_graph>0.                           (11)
```

Thus a positive answer to any old word-support landing test automatically
pays the graph-channel coefficient debt. The three debts need not be
attacked separately.

## 3. Equality mechanism and sharpness

Equality in the first comparison in (10) has a complete description:

1. for every nonzero `q`, the row `z -> A(q,z)` is constant; and
2. for every graph label `c`,

   ```text
   sum_(q!=0)A(q,q^2/2+c)=0.                        (12)
   ```

The first condition is exactly equality in every instance of (7), and
the second kills the final nonnegative term in (6). Equality in the
masked bound also requires `A_0(q)=0` off the selected mask among
nonzero targets.

The constant is attained for every one of the three masks. Choose two
distinct points `q_1,q_2` in the support of `P_sigma` and
`alpha!=0`. Set, for every `c in K`,

```text
A(q_1,q_1^2/2+c)= alpha,

A(q_2,q_2^2/2+c)=-alpha,                            (13)
```

and set all other coefficients to zero. Every nonzero row in (13) is
constant in `z`, while the two entries on each graph cancel. Consequently

```text
D_graph=2N|alpha|^2,

E_sigma=2N^2|alpha|^2,

D_graph/E_sigma=1/N=1/169.                         (14)
```

The supports of the three masks have sizes `12,12,144`, so the required
pair always exists. No larger uniform constant is possible.

## 4. The graph detector is strictly finer

The reverse implication in (11) is false. Fix `q!=0`, distinct jet
values `z_1,z_2`, and `alpha!=0`, and put

```text
A(q,z_1)=alpha,

A(q,z_2)=-alpha,                                   (15)
```

with all other coefficients zero. Then

```text
A_0(q)=0,

E_sigma=0             for all sigma,               (16)
```

but the first term in (6) is `2|alpha|^2`, so

```text
D_graph>0.                                         (17)
```

Thus `D_graph` sees jet-resolved nonzero-target survival that is erased
by the old sum over `z`. The exact hierarchy is

```text
joint nonzero-target coefficient
  -> D_graph>0;

word-support aggregate survives
  -> D_graph>0;

D_graph>0
  -/-> word-support aggregate survives.             (18)
```

This identifies one common positive target rather than two unrelated
quadratic obligations.

## 5. Scope

The comparison is independent of the chosen coordinate presentation
once the graph foliation and the word masks have been fixed, but those
identifications remain noncanonical as in THM-2356.

Most importantly, `D_c` is still a derived coefficient functional. Its
physical reconstruction needs the cross-frequency term

```text
F_c(b)conjugate(F_c(0)),
```

or an equivalent graph-channel pair twist. Separate chirp intensities do
not determine it. Therefore (10) does not realize a new physical probe,
exclude the vertical tensor hostile, decrement any of the `165` LRC
profiles, force a scalar-row contradiction, or prove LRC(14).

It does sharpen the missing sidecar specification: a lawful realization
of the summed graph detector would simultaneously settle every positive
THM-2337 word-support energy, while retaining extra jet-resolved target
information.

## 6. Exact companion

The dependency-free companion works in

```text
F_169=F_13[t]/(t^2-2)
```

and uses Gaussian-integer coefficient arrays. It:

- exhausts the `169` graph-coordinate bijections;
- checks (6)--(10) on `80` deterministic hostile arrays;
- realizes equality (13)--(14) separately on all three masks; and
- realizes the strict jet-cancellation boundary (15)--(17).

Run

```bash
python3 04-computation/planar_graph_detector_word_energy_thm2363.py
python3 -O 04-computation/planar_graph_detector_word_energy_thm2363.py
```

Both transcripts must match

```text
05-knowledge/results/planar_graph_detector_word_energy_thm2363.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independent audit is pending. QED.
