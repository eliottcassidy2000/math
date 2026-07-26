---
id: THM-2405
title: "Two-level septimal sheet independence and middle-depth-two cage elimination"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In
  THM-2392's b=2 same-septimal-layer cage, the high quotient band
  D_C3 minus D_c3 outside the top word has exact mass 72/637. Every
  primitive low danger comb occupies exactly one seventh of this
  invariant set, so the two-low high-band contribution is at most
  144/4459. Consequently delta is at least 324/4459-mu(U_(a,b)).
  This forces explicit positive clean-mass floors for the oriented
  ratios (4,1), (3,4), and (4,3), eliminates the sole M=2 no-clean
  ratio, and leaves only seven oriented ratios at M=1. This is a
  direct-cage sharpening; it is not a profile-row decrement,
  canonical target landing, or proof of LRC(14).
source: codex-2026-07-26-two-level-septimal-cage
depends_on:
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
  - THM-2392-clean-toothpick-or-bounded-cross-ancestor-cage
related:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2388-thirteen-root-multiplicity-reflection-and-blocker-caged-toothpick-law
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
script: 04-computation/lrc14_two_level_septimal_cage_thm2405.py
output: 05-knowledge/results/lrc14_two_level_septimal_cage_thm2405.out
script_sha256: e0b9f5af824500eb57aa2a17706d786749671acff2cfc3080b3b88d40eeddd89
output_sha256: 812f400c611901c9c59b4be0d9012eec190b2311a21e6c7faa663b5d59d1a42a
hash_basis: working-tree bytes (LF)
---

# THM-2405 -- two septimal sheet counts eliminate the middle-depth-two cage

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2392 reduces the last fifteen strict profile rows to a clean-hole
alternative and an exact bank of ten oriented compatible low ratios.
THM-2391 further says that a no-clean packet has septimal depth

```text
M=1: one of the ten ratios;

M=2: only (a,b)=(4,3);

M>=3: impossible.                                      (1)
```

The cap used there treats the two pieces involving `C_3` separately.
That loses two exact coordinates: both high pieces live in the same
parent/child band, and both primitive low combs occupy one sheet of
the same seven-cover. Retaining those coordinates gives

```text
high band outside q_*:                    72/637;

one primitive low sheet in that band:     72/4459;

two primitive low sheets:                <=144/4459.   (2)
```

This small change crosses the `M=2` boundary. It also removes three
of the ten `M=1` orientations.

## 1. Exact one-sheet independence

Write

```text
D_v={x in R/Z: ||vx||<1/14}.                          (3)
```

All identities below are almost everywhere. The strict-open endpoints
form finite null sets and are never assigned a sheet.

> **Seven-sheet lemma.** Let `r>=0`, let
>
> ```text
> v=7^r u,                         7 does not divide u,
> ```
>
> and let `E` be a measurable subset of the circle invariant under
> translation by `1/7^(r+1)`. Then
>
> ```text
> mu(E intersection D_v)=mu(E)/7.                     (4)
> ```

Indeed, for almost every `x`,

```text
sum_(s=0)^6 1_(D_v)(x+s/7^(r+1))=1.                  (5)
```

After multiplication by `v`, the seven arguments differ by the
complete grid `us/7`; the strict interval of length `1/7` contains
exactly one of them. Invariance of `E`, translation invariance of Haar
measure, and integration of (5) prove (4).

This is literal one-sheet counting, not probabilistic independence.
It remains valid for an arbitrary invariant measurable `E`.

## 2. The common high band

Retain THM-2391/2392 notation in the last septimal lane:

```text
c_j=13C_j,

nu_7(q_*)=M<nu_7(C_3)=nu_7(c_3),

nu_7(c_1)=nu_7(c_2)=0.                                (6)
```

Define the same-line high band

```text
H_3=D_(C_3) intersection D_(c_3)^c.                  (7)
```

Multiplication by `C_3` preserves Haar measure. The intersection
`D_1 intersection D_13` is exactly the central interval of total
length `1/91`; the two neighboring teeth meet the strict endpoints
only. Therefore

```text
mu(H_3)=1/7-1/91=12/91.                               (8)
```

Because both `C_3` and `c_3` are divisible by `7^(M+1)`, the band
`H_3` is invariant under translation by `1/7^(M+1)`. Apply (4) with
`r=M` and `v=q_*`:

```text
mu(H_3 intersection D_(q_*))=(12/91)/7=12/637.       (9)
```

Thus the part of the high band outside the top word,

```text
A=H_3 intersection D_(q_*)^c,                        (10)
```

has the exact mass

```text
mu(A)=12/91-12/637=72/637.                            (11)
```

Every speed in (10) is divisible by seven, so `A` is invariant under
translation by `1/7`. The original low blockers `c_1,c_2` are
seven-units by (6). A second application of (4), now with `r=0`,
gives

```text
mu(A intersection D_(c_i))=mu(A)/7=72/4459,
                                                     i=1,2.          (12)
```

Consequently

```text
mu(A intersection (D_(c_1) union D_(c_2)))
 <=144/4459.                                         (13)
```

The constant in (13) is sharp at the level of its stated invariant-set
hypotheses. For example, the exact strict-open packet

```text
q_*=7,       C_3=49,       c_3=637,

(c_1,c_2)=(27,29)                                    (14)
```

has disjoint selected low sheets throughout `A` and attains
`144/4459`. This is an abstract positive control, not a canonical
LRC packet with the correlated four-speed roles below.

## 3. Direct sharpening of the cage

In the `b=2` same-septimal boundary, divide by the common low scale and
orient as in THM-2392:

```text
(C_1,C_2,c_1,c_2)=h(b,13a,13b,169a),

gcd(a,b)=1,                    gcd(ab,91)=1.           (15)
```

Put

```text
L_Q=D_(C_1) union D_(C_2),

L_O=D_(c_1) union D_(c_2),

U_(a,b)=L_Q intersection L_O.                         (16)
```

THM-2392's geometric cage is, modulo endpoints,

```text
Gamma
 =(L_Q union D_(C_3))
   intersection L_O
   intersection D_(q_*)^c
   intersection D_(c_3)^c.                           (17)
```

Split (17) according to whether the quotient-blocker contribution
comes from `L_Q` or `D_(C_3)`. The first piece is contained in
`U_(a,b)`. The second is contained in `A intersection L_O`, so (13)
gives

```text
mu(Gamma)<=mu(U_(a,b))+144/4459.                      (18)
```

If `delta` is THM-2392's clean-hole mass, its exact cage invoice is

```text
mu(Gamma)>=36/343-delta.                              (19)
```

Combining (18)--(19), and using `4459=13*343`, yields the central
inequality

```text
delta
 >=36/343-144/4459-mu(U_(a,b))
 =324/4459-mu(U_(a,b)).                               (20)
```

No pairwise independence, union heuristic, or unproved moving-root
alignment occurs in (20).

## 4. Exact ratio table and the three new positive branches

Multiplication by the common scale `h` preserves Haar measure, so the
compatible union in (16) is evaluated at

```text
(b,13a,13b,169a).                                    (21)
```

The complete THM-2392 bank and its exact union masses are:

| `(a,b)` | `mu(U_(a,b))` | consequence of (20) |
|---|---:|---:|
| `(1,1)` | `193/1183` | retained |
| `(1,2)` | `114/1183` | retained |
| `(2,1)` | `239/2366` | retained |
| `(1,3)` | `263/3549` | retained |
| `(3,1)` | `95/1183` | retained |
| `(4,1)` | `331/4732` | `delta>=629/231868` |
| `(2,3)` | `43/546` | retained |
| `(3,2)` | `95/1183` | retained |
| `(3,4)` | `491/7098` | `delta>=1213/347802` |
| `(4,3)` | `331/4732` | `delta>=629/231868` |

Thus a no-clean cage must satisfy

```text
mu(U_(a,b))>=324/4459,                               (22)
```

and its oriented ratio belongs to the seven-element bank

```text
(1,1), (1,2), (2,1), (1,3), (3,1), (2,3), (3,2).
                                                               (23)
```

For the three positive branches, THM-2392's literal top-labelled
charged cell has mass at least `delta/52`. Hence:

```text
(4,1),(4,3):
  rho_charged>=629/12057136;

(3,4):
  rho_charged>=1213/18085704.                         (24)
```

Its owner-resolved `C_7 x C_13` tensor cell has mass at least
`delta/338`, giving

```text
(4,1),(4,3):
  rho_7x13>=629/78371384;

(3,4):
  rho_7x13>=1213/117557076.                           (25)
```

All nonzero target colours and all septimal colours on those cells
retain the exact nonvanishing statements of THM-2392.

## 5. The depth-two boundary is empty

THM-2392 proves that a no-clean packet at septimal depth `M=2` must
have the unique oriented ratio

```text
(a,b)=(4,3).                                         (26)
```

But (20) gives

```text
delta>=629/231868>0                                  (27)
```

on that ratio. Therefore:

```text
M=2 no-clean boundary: empty;

M>=3 no-clean boundary: already empty by THM-2392;

remaining no-clean boundary:
  M=1 and one of the seven ratios in (23).             (28)
```

This is the promised middle-depth-two elimination. It is stronger
than a finite congruence scan: the obstruction is the exact common
high-band sheet budget.

## 6. Equality and failure boundaries

The hypotheses in the sheet argument are load-bearing.

1. If the selected low speed is divisible by seven, (4) can fail
   maximally: take `E=D_7` and `v=7`. Then

   ```text
   mu(E intersection D_v)=1/7,       not 1/49.         (29)
   ```

2. If `q_*` is not strictly below the high band in septimal depth,
   (9) can fail maximally. Taking `q_*=C_3=49` makes

   ```text
   mu(H_3 intersection D_(q_*))=12/91,
   ```

   rather than `12/637`.

3. Equality in (13) means that the unique `c_1` and `c_2` sheets are
   distinct almost everywhere on `A`. The witness (14) shows that
   the factor two cannot be reduced from invariant-set information
   alone.

4. The seven residual ratios in (23) genuinely survive this cap.
   Further progress must use their carry/translate data, the
   thirteen-drift identity, or a lawful endpoint/owner sidecar; the
   scalar mass invoice has been exhausted.

## 7. Scope

This theorem sharpens the **direct THM-2392 cage**, not the separate
global common-core clean-mass floor of THM-2396.

It proves positive clean-hole mass in the three displayed ratio
branches and removes the no-clean `M=2` branch. It does not:

- exclude any of the fifteen thirteen-adic profile rows;
- transport the charged root word to a canonical expiration target;
- align the transverse `C_7` and predecessor `C_13` fibres physically;
- settle any of the seven `M=1` oriented ratios; or
- prove LRC(14).

The scalar row ledger therefore remains `165`.

## 8. Exact companion

The dependency-free exact companion:

- checks `7,344` nonendpoint seven-sheet identities for seven-unit
  speeds;
- reconstructs the exact same-line mass `1/91` and high-band mass
  `12/91`;
- verifies `mu(A)=72/637` in four independent depth/scale controls;
- verifies `48` exact one-low-sheet intersections `72/4459`;
- attains the abstract two-low cap `144/4459`;
- supplies explicit nonunit-low and same-depth hostiles;
- independently reconstructs all ten compatible union masses in the
  table above;
- verifies the seven-ratio residual and all clean/charged/tensor
  fractions; and
- checks that the sole `M=2` ratio is eliminated.

Run

```bash
python3 04-computation/lrc14_two_level_septimal_cage_thm2405.py
python3 -O 04-computation/lrc14_two_level_septimal_cage_thm2405.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_two_level_septimal_cage_thm2405.out
```

after LF normalization. Every executable check raises explicitly under
optimized Python.

## 9. Independent hostile audit

An independent reconstruction audited the theorem at three levels.

1. Starting from THM-2392's literal definition of `Gamma`, it recovered
   (17) after using `T^(-1)D_(C_i)=D_(c_i)`, and hence recovered the
   containment in `U_(a,b) union (A intersection L_O)` without importing
   this proof's interval decomposition.
2. Equal-cell lattice counting, rather than the companion's endpoint
   sweep, reproduced all ten exact `mu(U_(a,b))` values. In particular,
   `(4,3)` gives `331/4732`, so

   ```text
   delta>=629/231868>0.
   ```

3. A separate congruence scan recovered the unique `M=2` bank
   `{(4,3)}` and the empty `M>=3` bank. Normal, optimized, and stored
   transcripts byte-match.

The audit also compared a proposed weighted-dual shift-matrix route.
That route would require additional blocker-address shift coverage; the
common-band proof above does not. It is therefore unnecessary, not an
implicit dependency of this theorem.
