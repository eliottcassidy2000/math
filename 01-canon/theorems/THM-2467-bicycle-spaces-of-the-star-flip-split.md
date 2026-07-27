---
id: THM-2467
title: "The bicycle spaces of the star-flip split: even-cut law on K_n, mod-12 law on the tile graph"
status: >
  PROVED (dim bicycle(K_n) = (n-2)[n even]: a cut delta(S) of K_n
  is a cycle vector iff n is even and |S| is even, and even-size
  cuts span a codimension-one subspace of the cut space; for odd n
  the CLAUDE.md-era Cut (+) Cycle split of GF(2)^{E(K_n)} is direct,
  quantifying the zoo table's 'direct iff n odd') + FINITE-EXACT
  (tile graph K_n minus the base path: dim bicycle in {0,1} for all
  4 <= n <= 30, and = 1 exactly when n = 2, 3, 6, 9, 10 (mod 12) --
  a palindromic residue set closed under n -> 12 - n; the competing
  mod-9 fit matches through n = 18 and BREAKS at n = 21, the
  MISTAKE-055 five-point-fit motif, kept as a hostile control).
  This computes klein-S399's 'named, never computed' top-pick
  object (TOURNAMENT-INVARIANT-ZOO-ATLAS SS II.e #1). UPGRADE
  (S135): the mod-12 law is now PROVED for all n >= 5 -- the parity
  system linearizes to the two-term recurrence
  x_{v+1} = x_{v-1} + sigma + n x_v over GF(2) with endpoint and
  sum-consistency constraints, so solvability depends only on
  n mod 12; the function of the residue is evaluated by the
  exhaustive table n = 5..60 (56 points, every residue class hit
  at least four times).
source: kind-pasteur-2026-07-26-S134
depends_on:
  - THM-1405-star-quotient-is-the-cycle-space-transverse-to-isomorphism
related:
  - THM-1382-map-graph-star-flips-nested-cut-cycle
  - THM-2444-pure-blue-count-refutation-and-rigid-class-law
script: 04-computation/metagraph_bicycle_spaces_kps_S134.py
output: 05-knowledge/results/metagraph_bicycle_spaces_kps_S134.out
script_sha256: 8bf895460a726b965d8a80f1554119546a01e289d7a859e479b7c4d3929eb863
output_sha256: 4a6a1915e757d192064a82df892d8882733d7e6affe0909425cb54eaf3740237
hash_basis: working-tree bytes (LF)
---

# THM-2467 -- the two bicycle spaces, finally computed

**PROVED + FINITE-EXACT** as itemized in the status.

klein-S399's invariant zoo names the bicycle space `Cut intersect
Cycle` of the star-flip GF(2) split as its top never-computed
object. There are two ambient graphs in play (THM-1405's
distinction), and they have different answers.

## 1. K_n with the base-path spanning tree (PROVED)

`delta(S)` has degree `n - |S|` at vertices of `S` and `|S|`
outside, so it is a cycle vector (all degrees even) iff `n` is even
and `|S|` is even. Even-size cuts form a codimension-one subspace
of the (n-1)-dimensional cut space, hence

```text
dim bicycle(K_n) = (n - 2) [n even]     (verified n = 4..18).   (1)
```

For odd `n` the split `GF(2)^{E(K_n)} = Cut (+) Cycle` is direct --
the zoo's "direct iff n odd", now with the even-n defect quantified
as a full `(n-2)`-dimensional degeneracy, NOT a small correction.
Any THM-1382-era argument that pairs cut and cycle coordinates at
even n must quotient these `n-2` dimensions first.

## 2. The tile graph K_n minus the base path (FINITE-EXACT)

```text
n            : 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 ... 30
dim bicycle  : 0 0 1 0 0 1 1  0  0  0  1  1  0  0  1  ... 1

dim in {0, 1} throughout, and dim = 1  iff
n = 2, 3, 6, 9, 10 (mod 12).                                  (2)
```

The residue set `{2,3,6,9,10}` is closed under `n -> 12 - n`.
Hostile control kept in the companion: the mod-9 fit
`{0,1,5,6} mod 9` reproduces all data through `n = 18` and breaks
at `n = 21` -- the recorded five-point-fit failure mode
(MISTAKE-055) caught by extension before adoption.

## 2b. Proof of the mod-12 law (S135 upgrade)

Write `x_v = [v in S]`, `sigma = |S| mod 2`. In the tile graph,
`N(v) = V minus {v-1, v, v+1}`, so `delta(S)` has all degrees even
iff

```text
interior v : x_{v-1} + x_{v+1} = sigma + n x_v        (mod 2)
v = 1      : x_2 = sigma + (n+1) x_1;   v = n symmetric.       (3)
```

Given `(x_1, sigma)` the recurrence determines the whole sequence.
For even `n` it reads `x_{v+1} = x_{v-1} + sigma`: period 2
(`sigma = 0`) or 4 (`sigma = 1`). For odd `n` it is the
Fibonacci-type `x_{v+1} = x_v + x_{v-1} + sigma`: period dividing 6
(Pisano pi(2) = 3, doubled by the affine constant). The terminal
condition at `v = n` and the consistency `sigma = sum x_v` are
therefore functions of `n mod lcm(4, 6) = 12` only, for `n >= 5`
(below 5 the endpoint and interior regimes overlap). Hence
`dim bicycle` is a function of `n mod 12` on `n >= 5`, and the
exhaustive table `n = 5..60` evaluates it: `1` on residues
`{2, 3, 6, 9, 10}`, `0` elsewhere, dimension never exceeding one
(at most the four seeds `(x_1, sigma)`, halved by
`delta(S) = delta(V minus S)`, minus the trivial solutions). QED

## 3. Reading and open

At the `n = 2, 3, 6, 9, 10 mod 12` positions there is exactly one
nonzero GF(2) vector that is simultaneously a star-flip combination
and a tile-cycle holonomy -- a single "gauge-invariant Wilson loop"
that THM-1405's transversality argument must except; everywhere
else the star quotient and the cycle holonomies are cleanly
complementary. Open: identify the exceptional vector's tournament meaning at the
first hit `n = 6` (the same n where blue self-loops first appear,
THM-648 -- possibly not a coincidence, unverified).

## 4. Reproduction

```bash
python 04-computation/metagraph_bicycle_spaces_kps_S134.py
```

Exact GF(2) ranks; hard failures; `ALL CHECKS PASSED`.
