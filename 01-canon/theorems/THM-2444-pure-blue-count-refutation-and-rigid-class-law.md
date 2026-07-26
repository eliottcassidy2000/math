---
id: THM-2444
title: "Pure-blue count formula refuted at n=9; the rigid-class law"
status: >
  REFUTED (HYP-4997/THM-791(d) closed form pure-blue(n) =
  floor((n+1)/2) - [n even] fails at its first untested point:
  pure-blue(9) = 6, predicted 5) + PROVED (rigidity lemma: every
  self-converse class with H = |Aut| is pure-blue) + FINITE-EXACT
  (exhaustive blue-sub-cube census n = 5..9 with n <= 8 canon
  controls reproduced: 3, 2, 4, 3). The six pure-blue classes at
  n = 9 are the five rigid classes H = |Aut| in {1, 3, 9, 9, 27}
  (all powers of three) plus the unique nonrigid class
  (H, |Aut|, tc) = (15, 5, 3). The theorem does not give a closed
  form for the rigid-class count and does not decide n = 10.
source: kind-pasteur-2026-07-26-S131
depends_on:
  - THM-643-gridsym-tiling-line-structure
  - THM-644-gridsym-fiber-law-and-anti-redei
related:
  - THM-790-blue-parity-law-proved
  - THM-791-H-companion-laws-to-the-transitivity-flow
  - THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay
script: 04-computation/metagraph_pureblue_n9_kps_S131.py
output: 05-knowledge/results/metagraph_pureblue_n9_kps_S131.out
script_sha256: f1a65a545b62ae0f47cf424205b264907fa43d5777de02bccfc2a99e008749b6
output_sha256: b3cd2a0ee09f7df035731152ae5f4d3ce96265bacdc9b1652c44467a120cf968
hash_basis: working-tree bytes (LF)
---

# THM-2444 -- pure-blue(9) = 6: the interleave formula was a rigid-count shadow

**REFUTED + PROVED + FINITE-EXACT** as itemized in the status.

## 1. The refutation

HYP-4997 (kps-S66 PB1, restated as THM-791(d)) conjectured

```text
pure-blue(n) = floor((n+1)/2) - [n even]     (2,1,3,2,4,3,5,...)  (1)
```

verified at n = 3..8 (six points, the n = 8 value 3 from opus-S305's
certified atlas). The first untested point refutes it:

```text
pure-blue(9) = 6,   (1) predicts 5.                              (2)
```

Census (blue-sub-cube method of kps-S66, upgraded WL-bucket
canonicalizer; universe = all grid-symmetric tilings, 2^16 = 65,536
at n = 9):

```text
n : touched SC classes : pure-blue : class inventory (H, |Aut|, tc)
5 :    8               : 3         : (1,1,1) (3,3,1) (15,5,3)
6 :   12               : 2         : (1,1,1) (9,9,1)
7 :   88               : 4         : (1,1,1) (3,3,1) (9,9,1) (15,5,3)
8 :  176               : 3         : (1,1,1) (9,9,1) (9,9,1)
9 : 2752               : 6         : (1,1,1) (3,3,1) (9,9,1) (9,9,1)
                                     (27,27,1) (15,5,3)
```

Controls: n = 5..8 reproduce canon (3, 2, 4, 3) exactly, including
THM-791(d)'s four n = 7 classes; the touched-class counts equal the
self-converse counts SC(n) = 8, 12, 88, 176, 2752 (A002785), as
THM-643 T2 requires (every SC class has an odd >= 1 number of blue
tilings, and blue lines only touch SC). The known-buggy
`purity_n8.out` PB row is not used anywhere.

## 2. The rigidity lemma (PROVED)

**Lemma.** Every self-converse class `C` with `H(C) = |Aut(C)|` is
pure-blue.

*Proof.* By LEM-003 the (unmerged) tiling fibre of `C` has exactly
`tc = H/|Aut| = 1` element, say `t`. The grid reflection `sigma`
sends the fibre of `C` to the fibre of `C^op` (THM-644(a)); for
self-converse `C` these coincide, so `sigma` permutes a one-element
set: `sigma(t) = t`. Hence the entire fibre is grid-symmetric, which
is the definition of pure-blue. QED.

So rigid SC classes (`tc = 1`) are pure-blue for free, and the
census splits

```text
pure-blue(n) = rigid-SC(n) + nonrigid-pure-blue(n),
rigid-SC(n)          = 2, 1, 2, 2, 3, 3, 5   (n = 3..9)
nonrigid-pure-blue(n)= 0, 0, 1, 0, 1, 0, 1   (n = 3..9).
```

The interleaved formula (1) was numerically correct through n = 8
only because `rigid-SC(n)` happened to match
`floor((n+1)/2) - [n even] - [n odd, n >= 5]` there; at n = 9 the
rigid count jumps by two (the second `(9,9,1)` class and the new
`(27,27,1)` class) while the formula steps by one. The failure is in
the rigid stratum, not in the nonrigid one.

## 3. Structure of the witnesses

1. All rigid pure-blue classes observed have `H = |Aut|` a **power of
   three** (`1, 3, 9, 27`), consistent with THM-643 C1's 3-power cap
   on `H_sym` and with the odd-`H`/odd-`|Aut|` laws (THM-643 T1). The
   first `3^3` appears exactly at n = 9, the first n admitting an
   iterated three-block structure `3x3` -- the natural candidates are
   iterated `C_3`-lexicographic products, whose converse-invariance
   and regular automorphism action force `H = |Aut|`.
2. The unique nonrigid pure-blue class is `(15, 5, 3)` at every odd
   `n in {5, 7, 9}` and absent at every even `n in {4, 6, 8}` -- at
   n = 5 it is the Paley/quadratic-residue circulant `T_5`; at
   n = 7, 9 its `H = 15, |Aut| = 5` signature persists (a `T_5`-core
   augmentation), while the regular classes are not pure-blue
   (THM-791(d) n = 7 refinement, reconfirmed here).

## 4. Repaired form and new questions

Repaired working conjecture (OPEN):

```text
pure-blue(n) = rigid-SC(n) + [n odd, n >= 5],                    (3)
```

with the entire growth carried by `rigid-SC(n)` = the number of
self-converse classes with a one-tiling fibre. New questions:

- Closed form or growth law for `rigid-SC(n)` (data
  `2, 1, 2, 2, 3, 3, 5`): is it the count of converse-invariant
  iterated odd-cyclic lexicographic towers on `n` vertices?
- Is `(15, 5, 3)` truly the **only** nonrigid pure-blue class at
  every odd n, and why does no nonrigid class survive at even n?
  (THM-648's even/odd blue self-loop asymmetry is the obvious
  suspect mechanism.)
- n = 10 decides between (3) (predicting `rigid-SC(10)`) and any
  surviving interleave variant; the same script at n = 10 is a
  2^20-tiling run.

## 5. What breaks downstream

THM-791(d)'s clause "(1) confirmed through n = 7" remains true as a
finite statement; its extrapolation is dead. kps-S66's PB2 pendant
law and blue-mass bookkeeping are untouched (they are per-class
statements). No proved theorem depended on (1) for n >= 9.

## 6. Reproduction

```bash
python 04-computation/metagraph_pureblue_n9_kps_S131.py 5 6 7 8 9
```

Exhaustive over each blue sub-cube; exact integer H, |Aut|, tc;
asserts the n <= 8 canon controls and the fibre-sum identity
`sum blue-mult = 2^e`; prints `ALL CONTROLS PASSED`.
