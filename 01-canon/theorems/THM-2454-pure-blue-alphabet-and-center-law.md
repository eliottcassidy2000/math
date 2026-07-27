---
id: THM-2454
title: "The pure-blue alphabet: palindromic stacks with rigid halves and atomic centers"
status: >
  FINITE-EXACT (pure-blue(11) = 9 by exhaustive 2^25 blue-cube
  census, C++ engine with n = 9, 10 controls reproducing canon; the
  nine census points n = 3..11 all match the closed form
  pure-blue(n) = c(n/2) for even n and
  c((n-1)/2) + c((n-3)/2) + c((n-5)/2) for odd n, c = A000930) +
  PROVED (center law: a palindromic stack is pure-blue iff its half
  parts are rigid and its center is pure-blue; the sigma-reversal
  fibre factorization is verified by the decisive mixed witness
  stack(T5,T5): blue-mult 3 vs tc 9) + CONFIRMED (THM-2453's locked
  rigid prediction at n = 11: exactly 7 rigid classes, H-multiset
  [1,3,9,9,9,27,27]) + REFUTED (the "+1 per odd n" clause of
  THM-2444's repaired form: n = 11 has TWO nonrigid pure-blues,
  (15,5,3) = (1,1,1,T5,1,1,1) and the new (135,45,3) = (C3,T5,C3)).
  The completeness of the atom alphabet {1, C3, T5} beyond n = 11
  and the strong-atom classification are OPEN.
source: kind-pasteur-2026-07-26-S132
depends_on:
  - THM-2453-palindromic-narayana-law-for-the-rigid-stratum
  - THM-2444-pure-blue-count-refutation-and-rigid-class-law
related:
  - THM-2450-rigid-self-converse-classes-are-cyclic-ternary-towers
  - THM-643-gridsym-tiling-line-structure
script: 04-computation/metagraph_pureblue_alphabet_kps_S132.py
output: 05-knowledge/results/metagraph_pureblue_alphabet_kps_S132.out
census_script: 04-computation/metagraph_pureblue_n11_kps_S132.cpp
census_output: 05-knowledge/results/metagraph_pureblue_n11_kps_S132.out
script_sha256: 07b6398f3bef5446b69b522d26be5b399d580bee5a0a3f5e4bc4f82da55f6436
output_sha256: 17e32bfba9e056534a13c5fd7324672879de06a6a5a5b4c58cc723356d176bb5
census_script_sha256: e4041a9a8aa76aa2364747e3499df5eb8bdccbb979ca162a3cd683e4dfce75fa
census_output_sha256: 79846933b6a6524c37172f3ce5d8b87619cae50f825b31946fcc488339079dc0
hash_basis: working-tree bytes (LF)
---

# THM-2454 -- pure-blue is spelled in a three-letter alphabet

**FINITE-EXACT + PROVED + CONFIRMED + REFUTED** as itemized.

## 1. The n = 11 census

The 2^25 blue-cube census (C++ engine; n = 9, 10 controls reproduce
canon 6 and 4 with exact inventories; touched classes
`2752, 8784, 279968` = the self-converse counts, extending A002785)
gives

```text
pure-blue(11) = 9:
rigid    (7): H = |Aut| in [1, 3, 9, 9, 9, 27, 27]
nonrigid (2): (H,|Aut|,tc) = (15, 5, 3) and (135, 45, 3).
```

The rigid stratum confirms THM-2453's pre-registered prediction
exactly. The two nonrigid classes identify as palindromic stacks
with a T5 center: `(1,1,1,T5,1,1,1)` (`15 = H(T5)`, `5 = |Aut T5|`)
and `(C3,T5,C3)` (`135 = 3*15*3`, `45 = 3*5*3`). This refutes the
"+1 nonrigid per odd n" clause of THM-2444 SS4 at its first
untested point -- and reveals the complete structure.

## 2. The center law (PROVED)

The grid involution acts on a stack's tiling fibre by reversing the
component order and reflecting each component (THM-644(a) plus the
canonical strong decomposition). Hence a sigma-fixed tiling of a
palindromic stack is exactly: an arbitrary tiling on each first-half
part, its mirror on the matching second-half part, and a sigma-fixed
tiling on the center. Therefore

```text
blue-mult = prod_{half} tc_i * blue-mult(center),
tc        = prod_{half} tc_i^2 * tc(center),
pure-blue <=> every half part is rigid (tc_i = 1)
              and the center is pure-blue.                       (1)
```

Decisive computational witness: `stack(T5, T5)` at n = 10 has
`blue-mult = 3 = tc(T5)` (direct 2^20 blue-cube enumeration) against
`tc = 9` -- mixed, exactly as (1) demands and as the n = 10 census
(no nonrigid pure-blue) required.

## 3. The complete law (FINITE-EXACT n = 3..11)

Through n = 11 the pure-blue strong atoms are `{1, C3, T5}` (T7 and
T11 are not pure-blue: their signatures are absent from the n = 7
and n = 11 censuses). With (1), pure-blue classes = palindromic
stacks with `{1, C3}` halves and center in `{none, 1, C3, T5}`:

```text
pure-blue(n) = c(n/2)                                   (n even)
pure-blue(n) = c((n-1)/2) + c((n-3)/2) + c((n-5)/2)     (n odd)    (2)
```

with `c = A000930`. Values `2,1,3,2,4,3,6,4,9` -- all nine censuses
match. Predictions: `pure-blue(12) = 6`, `13 -> 13`, `14 -> 9`,
`15 -> 19`, `16 -> 13`.

## 3b. n = 12 confirmation (same session)

The pruned survivor engine
(`04-computation/metagraph_pureblue_n12_pruned_kps_S133.cpp`, sha
`22a96e45...`; output `..._n12_pruned_kps_S133.out`, sha
`d0b2e6d1...`; sampler soundness: only false SURVIVORS possible,
eliminated by exhaustive path enumeration; controls n = 9, 10, 11
reproduce the exact known inventories) settles the even-line test:

```text
pure-blue(12) = 6, all rigid, H-multiset [1, 9, 9, 9, 9, 81]  (3)
```

exactly as predicted by (2), including the first `3^4 = 81` class
(the palindrome `C3 C3 C3 C3`) and zero nonrigid classes (the
center law's even-n exclusion). The 2^30 cube reduced to 163 raw
survivors in 89 classes. The complete law now stands at TEN census
points, with both pre-registered predictions (n = 11 rigid
inventory; n = 12 full inventory) confirmed.

## 3c. n = 13 confirmation (S134): the atom test passes

The threaded pruned engine
(`04-computation/metagraph_pureblue_n13_threaded_kps_S134.cpp`, sha
`149e60b8...`; output sha `94fcdfe4...`; n = 11, 12 controls exact;
2^36 cube -> 21 raw survivors in 15 classes) delivers

```text
pure-blue(13) = 13:
rigid   (10): [1, 3, 9, 9, 9, 9, 27, 27, 27, 81] = pal13(13)
nonrigid (3): (15,5,3) + 2 x (135,45,3)
            = T5-centered palindromes with {1,3} halves of sum 4  (4)
```

-- the pre-registered inventory (`..._n13_prediction_kps_S134.out`,
committed before the run) confirmed CLASS FOR CLASS. The law (2)
now stands at ELEVEN exhaustive census points (n = 3..13) with
three consecutive pre-registered predictions confirmed exactly
(n = 11 rigid; n = 12 full; n = 13 full), and **no pure-blue strong
atom beyond {1, C3, T5} exists through size 13** -- the first
genuinely open odd window for a new atom is n = 15 (as a size-15
strong atom) or any odd n via larger centers.

## 4. Open

- Completeness of the atom alphabet: does any new pure-blue strong
  atom appear at n >= 12? (A new atom of size s first affects (2)
  at n = s, odd s, or via centers at n = s + 2k.) The next census
  (n = 12, predicted all-rigid count 6, first 81-class) tests the
  even line cheaply.
- Classification of pure-blue strong tournaments: why exactly the
  rotational T5 and not T7/T11. Exact data
  (metagraph_rotational_bluemult_kps_S133, sha `a269a9ea...`/
  `bcb27843...`): the rotational family has

  ```text
  p   H       |Aut|  tc     blue-mult   status
  3   3       3      1      1           pure
  5   15      5      3      3           pure (all 3 regular blue
                                         tilings at n=5 are T5's)
  7   189     21     9      3           impure
  11  95095   55     1729   37          impure
  ```

  Note `tc(T11) = 1729`, the taxicab number. The `tc = 3^{(p-3)/2}`
  guess dies at `p = 11`, and blue-mult does not stall at 3 either;
  the honest law is that the purity ratio `bm/tc` collapses
  (`1, 1, 1/3, 37/1729`), so T5 is pure because `p = 5` is the last
  point where the regular blue stratum is small enough to coincide
  with the fibre. A structural proof of `bm(T_p) < tc(T_p)` for all
  `p >= 7` is the remaining open item.
- The arc THM-4997 -> 2444 -> 2453 -> 2454 is a candidate
  META-PATTERNS case study: two refutations, each at the first
  untested point, each converting a numerological formula into a
  structural law.

## 5. Reproduction

```bash
g++ -O2 -std=c++17 -o pb11 04-computation/metagraph_pureblue_n11_kps_S132.cpp && ./pb11 9 10 11
python 04-computation/metagraph_pureblue_alphabet_kps_S132.py
```

Both hard-fail on any mismatch; final lines `ALL CHECKS PASSED`.
