---
id: THM-542
title: LRC14 one-tail mouth-retention lemma - the only AP drop-6 one-tail row below the AP one-hole second value is the (6,10)->20 tail
status: PROVED (periodic-comb rational cutoff plus finite exact interval certificate)
source: codex-2026-06-19-S35
depends_on:
  - THM-541
  - HYP-2654
related:
  - HYP-2651
  - HYP-2650
  - HYP-2652
  - HYP-2569
external: Lonely Runner Conjecture, n=14
---

# THM-542 - One-Tail Mouth-Retention Lemma

Let

```text
C_{h,r} = ({1,2,...,13} \ {6,h}) union {r},
```

where `h in {1,...,13}\{6}` and `r >= 14` is an integer.  Define

```text
G_{h,r} = {t in [0,1): ||c t|| > 1/14 for every c in C_{h,r}}.
```

Then

```text
meas(G_{h,r}) < 426/35035
```

if and only if

```text
(h,r) = (10,20).
```

In that exceptional case

```text
C_{10,20} = (1,2,3,4,5,7,8,9,11,12,13,20),
```

and

```text
meas(G_{10,20}) = 3859/420420 = 7/858 + 1/980.
```

Moreover, the four old drop-6 mouth intervals from THM-541 survive intact:

```text
[29/182, 9/56]
[29/168, 27/154]
[127/154, 139/168]
[47/56, 153/182].
```

## Certificate

Script:

```text
04-computation/lrc14_one_tail_mouth_retention_theorem_codex_s35.py
```

Stored output:

```text
05-knowledge/results/lrc14_one_tail_mouth_retention_theorem_codex_s35.out
```

The script uses exact rational interval arithmetic.  No floating-point
comparison is used in the proof checks.

## Periodic-Comb Cutoff

Let `B_h={1,...,13}\{6,h}` and let `G_h` be its level-`1/14` safe set.  Suppose
`G_h` has measure `M_h` and `c_h` interval components.  For the speed-`r` danger
comb `D_r`, each period has length `1/r` and danger length `1/(7r)`.

For a single interval `I`, all full `1/r` periods remove exactly a `1/7` share
of `I`.  The two endpoint partial periods can add at most one danger tooth each,
so

```text
meas(I cap D_r) <= meas(I)/7 + 2/(7r).
```

Summing over the `c_h` components gives the exact lower bound

```text
meas(G_h \ D_r) >= (6/7)M_h - 2c_h/(7r).
```

Therefore it is enough to check

```text
r < R_h = ceil( 2c_h / (6M_h - 7*(426/35035)) ).
```

The exact cutoff table is:

```text
h | M_h          | c_h | 6M_h - 7Q      | R_h
1 | 239/3003     | 6   | 1964/5005      | 31
2 | 461/12012    | 6   | 1453/10010     | 83
3 | 1159/18018   | 8   | 4517/15015     | 54
4 | 389/12012    | 8   | 1093/10010     | 147
5 | 379/8580     | 12  | 1801/10010     | 134
7 | 4739/51480   | 10  | 2551/5460      | 43
8 | 2243/42042   | 8   | 8233/35035     | 69
9 | 11191/168168 | 10  | 44027/140140   | 64
10| 313/9702     | 8   | 11399/105105   | 148
11| 1951/25480   | 12  | 807/2156       | 65
12| 1546/35035   | 6   | 6294/35035     | 67
13| 907/17640    | 14  | 93917/420420   | 126
```

For every `r >= R_h`, the comb bound proves

```text
meas(G_{h,r}) >= 426/35035.
```

## Finite Residue

The remaining exact finite check covers `863` rows with `14 <= r < R_h`.
The only row below `426/35035` is:

```text
h=10, r=20,
safe=3859/420420,
old_survivor=7/858,
new_mass=1/980.
```

All assertions are made in the script as `Fraction` equalities/inequalities.
This proves the theorem.

## Proof Use

THM-542 proves the first infinite mouth-retention subcase of HYP-2654.  In the
one-tail AP-drop-6 family, being below the AP one-hole second value does not
create new behavior: the row must be the `(6,10)->20` tail, and that tail is
exactly the old drop-6 collar plus a new `1/980` mouth contribution.

This is an arithmetic statement about rational wall combs.  The moving object is
the integer speed `r`, and the proof is controlled by exact denominators and
finite integer cutoffs, not by approximating real-valued runner motion.

## Tournament Analysis

Vertices: proof gates for the one-tail family.

Pairwise observable: which gate removes more candidate rows exactly.

Switch/gauge: use the periodic-comb denominator cutoff before the finite wall
scan.

Hamiltonian path:

```text
comb_cutoff > finite_exact_scan > old_mouth_survivor > new_mouth_mass > raw_tail_size
```

Fingerprint: transitive proof-obligation order; directed `3`-cycles `0`.
