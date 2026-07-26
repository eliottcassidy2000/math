---
id: THM-2450
title: "Rigid self-converse classes are cyclic-ternary towers through n=10"
status: >
  FINITE-EXACT (through n = 10: the tower grammar
  T ::= 1 | stack(T1,T2) | C3[T1,T2,T3] generates rigid self-converse
  classes matching THM-2444's exhaustive census exactly, in count
  (2,1,2,2,3,3,5,4) and in every H-multiset, e.g. [1,3,9,9,27] at
  n=9) + PROVED (stack multiplicativity: H(stack(A,B)) = H(A)H(B)
  and Aut(stack(A,B)) = Aut(A) x Aut(B); hence a stack is rigid iff
  its strong components are rigid) + VERIFIED (controls: 132
  exhaustive multiplicativity pairs; the Paley circulant T_5 --
  the unique nonrigid pure-blue class -- has (H,|Aut|) = (15,5) and
  lies OUTSIDE the grammar, so the grammar does not overgenerate).
  The all-n identification of rigid-SC classes with converse-stable
  towers remains OPEN; the theorem does not prove it, and does not
  give a closed form for the tower counts.
source: kind-pasteur-2026-07-26-S131b
depends_on:
  - THM-2444-pure-blue-count-refutation-and-rigid-class-law
related:
  - THM-643-gridsym-tiling-line-structure
  - THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay
script: 04-computation/metagraph_rigid_towers_kps_S131b.py
output: 05-knowledge/results/metagraph_rigid_towers_kps_S131b.out
script_sha256: 7d7295eb6e444793a7cc079d8ee92156a700272f19dfe46715d80f578589b477
output_sha256: d3b93ce80e8fba8e83c99ba348c616b361339d339524dfa9e2cc28246b2f47db
hash_basis: working-tree bytes (LF)
---

# THM-2450 -- the rigid stratum is the prime 3's tower

**FINITE-EXACT + PROVED + VERIFIED** as itemized in the status.

THM-2444 showed the pure-blue count is carried by the rigid stratum
(self-converse classes with `H = |Aut|`, i.e. one-tiling fibres) and
left its structure open, noting all observed `H` values are powers
of three. This theorem identifies the stratum exactly through
`n = 10`.

## 1. The grammar and the match

Let the **cyclic-ternary tower grammar** be

```text
T ::= 1                (single vertex)
   |  stack(T1, T2)    (transitive join: all arcs T1 -> T2)
   |  C3[T1, T2, T3]   (three blocks arranged in a directed 3-cycle)
```

The isomorphism-deduplicated closure has
`1, 1, 2, 4, 9, 24, 65, 182, 528, 1558` classes for `n = 1..10`.
Testing every class exactly for `H`, `|Aut|`, and self-converseness:

```text
n  : rigid-SC in grammar : H multiset        : THM-2444 census
3  : 2                   : [1, 3]            : 2  match
4  : 1                   : [1]               : 1  match
5  : 2                   : [1, 3]            : 2  match
6  : 2                   : [1, 9]            : 2  match
7  : 3                   : [1, 3, 9]         : 3  match
8  : 3                   : [1, 9, 9]         : 3  match
9  : 5                   : [1, 3, 9, 9, 27]  : 5  match
10 : 4                   : [1, 9, 9, 9]      : 4  match
```

Since THM-2444's census is exhaustive over all classes (via the blue
sub-cube), the match in count and `H`-multiset at every
`n <= 10` establishes: **every rigid self-converse class through
n = 10 is a tower, and every rigid-SC tower appears in the census.**

## 2. Stack multiplicativity (PROVED)

In `stack(A, B)` every arc crosses forward, so a Hamiltonian path
must exhaust `A` before entering `B`: `H = H(A) H(B)`. The strong
components of a tournament and their linear order are canonical, so
every automorphism preserves the stack cut:
`Aut = Aut(A) x Aut(B)`. Hence rigidity (`H = |Aut|`) is stack-
multiplicative, and a stack is rigid iff its strong components are.
(Both laws additionally re-verified exhaustively on all 132 grammar
pairs with `|A| + |B| <= 7`.) The content of the census is therefore
concentrated in the strongly connected case: which `C3`-blocks of
towers are rigid.

## 3. Why the prime 3

The `H` values of rigid classes are powers of three because the only
strongly connected atoms available to the grammar are directed
3-cycles of blocks; each `C3`-layer contributes its odd cyclic
symmetry to both `H` and `|Aut|` in lockstep. This is the metagraph
face of THM-643 C1's `3^floor((n-2)/2)` cap on `H_sym`, and the
prime-3 counterpart of the prime-5 phenomenon on the nonrigid side:
the unique nonrigid pure-blue class is the Paley circulant `T_5`
(`H = 15`, `|Aut| = 5`, quadratic-residue character mod 5), which
the control confirms lies outside the grammar. In the harmonics
reading (S131b reflection): the rigid stratum is owned by the
3-tower, the nonrigid exception by the 5-character.

## 4. Open

- Prove the identification for all `n` (both directions: every
  rigid-SC class is a tower; a converse-stability criterion for
  which towers are rigid-SC).
- Closed form for the rigid-SC tower counts `2,1,2,2,3,3,5,4` and
  for the grammar closure counts `1,2,4,9,24,65,182,1558`
  (OEIS-absence not yet checked for the latter).
- Predict pure-blue(11) = rigid-SC-towers(11) + 1 (the `T_5`-type
  class at odd n) and test against a future blue-cube run at
  `2^25`.

## 5. Reproduction

```bash
python 04-computation/metagraph_rigid_towers_kps_S131b.py
```

Deterministic; hard failure on any mismatch or control violation;
final line `ALL CHECKS PASSED`.
