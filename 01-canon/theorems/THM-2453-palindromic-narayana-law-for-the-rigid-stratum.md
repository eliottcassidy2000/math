---
id: THM-2453
title: "The palindromic Narayana law for the rigid stratum"
status: >
  PROVED (every palindromic stack of singletons and C3's is a rigid
  self-converse class with H = |Aut| = 3^t, t = #C3 parts; their
  count is pal13(n) = c(n/2) for even n and c((n-1)/2) + c((n-3)/2)
  for odd n, with c = A000930, Narayana's cows) + FINITE-EXACT
  (rigid-SC(n) EQUALS pal13(n) with matching H-multisets for every
  n = 3..10 against the THM-2444 exhaustive census; strong towers
  C3[C3,1,1] and C3[C3,C3,C3] verified nonrigid, tc = 5 and > 1).
  Predictions locked before the 2^25 census: rigid-SC(11) = 7 with
  H-multiset [1,3,9,9,9,27,27] (hence pure-blue(11) = 8 under the
  THM-2444 repaired law), and rigid-SC(12) = 6 with
  [1,9,9,9,9,81]. The all-n converse (no rigid class outside the
  palindromic stacks) remains OPEN.
source: kind-pasteur-2026-07-26-S132
depends_on:
  - THM-2450-rigid-self-converse-classes-are-cyclic-ternary-towers
  - THM-2444-pure-blue-count-refutation-and-rigid-class-law
related:
  - THM-643-gridsym-tiling-line-structure
  - HYP-3003-summand-multiplicand-farey-basis-merge
script: 04-computation/metagraph_palindromic_narayana_kps_S132.py
output: 05-knowledge/results/metagraph_palindromic_narayana_kps_S132.out
script_sha256: d8a239cf062bbf850e07a7497dc4287a89756e0baee534d171a2d7b342be53f6
output_sha256: 1e034fac074e20ad7b70a8a194b4498706654fe6ab61c07a4b94c96c8b645c8b
hash_basis: working-tree bytes (LF)
---

# THM-2453 -- the rigid stratum is Narayana's palindromes

**PROVED + FINITE-EXACT** as itemized in the status.

THM-2450 identified the rigid stratum with tower-grammar classes
through n = 10. This theorem closes the structure completely on the
positive side: the rigid classes observed are not arbitrary towers
but exactly the **palindromic stacks over the alphabet
{single vertex, C3}**, and their count is a Narayana-cows palindrome
count.

## 1. The palindromic stack lemma (PROVED)

Let `w = (w_1, ..., w_m)` be a word in `{1, C3}` with sizes summing
to `n`, and `T(w)` the stack (transitive join) of its parts.

- **Self-converse iff palindrome.** The strong-component sequence of
  a tournament is canonical; the converse of `T(w)` is
  `T(reverse(w))` with each part conversed, and `C3^op = C3`,
  `1^op = 1`. Hence `T(w)` is self-converse iff `w` is a palindrome.
- **Always rigid.** By stack multiplicativity (THM-2450 SS2),
  `H = prod H(w_i) = 3^t` and `|Aut| = prod |Aut(w_i)| = 3^t`,
  `t = #{i : w_i = C3}`. So `tc = 1`: every palindromic stack is a
  one-tiling-fibre class, hence pure-blue (THM-2444 rigidity lemma).

Counting palindromes by their half-word: even `n` needs a half of
size `n/2`; odd `n` needs a centre `1` or `C3` with half
`(n-1)/2` or `(n-3)/2`. Compositions of `s` into parts `{1,3}`
number `c(s)` with `c(s) = c(s-1) + c(s-3)` -- **Narayana's cows**
(A000930). Hence

```text
pal13(n) = c(n/2)                     (n even)
pal13(n) = c((n-1)/2) + c((n-3)/2)    (n odd).                  (1)
```

## 2. Identification with the census (FINITE-EXACT n <= 10)

For every `n = 3..10`, `pal13(n)` and the palindromes' H-multiset
`{3^t}` match the THM-2444 exhaustive census of rigid classes
exactly: counts `2,1,2,2,3,3,5,4` and multisets `[1,3]`, `[1]`,
`[1,3]`, `[1,9]`, `[1,3,9]`, `[1,9,9]`, `[1,3,9,9,27]`,
`[1,9,9,9]`. So through `n = 10` **every** rigid self-converse
class is a palindromic `{1, C3}`-stack. Spot proofs that the
inclusion is strict inside the tower grammar: the strong towers
`C3[C3,1,1]` (`H = 15`, `|Aut| = 3`, `tc = 5`) and `C3[C3,C3,C3]`
are not rigid; the census data contain no rigid strong tower of
size `> 3`.

## 3. Predictions (locked before the 2^25 census run)

```text
rigid-SC(11) = 7,  H-multiset [1, 3, 9, 9, 9, 27, 27]
rigid-SC(12) = 6,  H-multiset [1, 9, 9, 9, 9, 81]
```

With the THM-2444 repaired law (`+1` nonrigid `(15,5,3)`-type class
at odd `n`), **pure-blue(11) = 8**. The `n = 11` blue-cube census
(C++ engine, in flight this session) decides. The first `3^4 = 81`
class is predicted at `n = 12` (the palindrome `C3 C3 C3 C3`).

## 4. The carry-family reading

`A000930` is the `d = 3` rung of the repo's Pascal-slope carry
family (HYP-3003): `d = 2` gives Fibonacci/Zeckendorf -- the spine
of the owner's summand-graph puzzle -- and `d = 3` gives Narayana's
cows, now realized as the metagraph's rigid stratum. The two
threads the owner dispatched separately (summand graphs; merged
metagraph) meet on consecutive rungs of one additive ladder: parts
of size `{1, 2}` for the integers' minimal summand chain, parts of
size `{1, 3}` for the tournaments' minimal symmetric atoms (the
singleton and the 3-cycle -- there is no size-2 tournament atom,
which is exactly why the metagraph rung is `d = 3`).

## 5. Open

- The all-n converse: no rigid-SC class outside the palindromic
  stacks (equivalently: no strong tower of size > 3 is rigid, and
  no non-tower class is rigid). The first half looks provable from
  the run-transfer structure of `H(C3[A,B,C])` (the ordered-path-
  system expansion makes `H` strictly exceed `3 |Aut A||Aut B|
  |Aut C|` once any block exceeds a point -- not yet written down).
- Whether the nonrigid pure-blue `(15,5,3)`-type class admits a
  matching closed law at all odd n (THM-2444 SS4).

## 6. Reproduction

```bash
python 04-computation/metagraph_palindromic_narayana_kps_S132.py
```

Exhaustive palindrome enumeration vs closed form vs census;
explicit `H`/`|Aut|` spot proofs; hard failures; final line
`ALL CHECKS PASSED`.
