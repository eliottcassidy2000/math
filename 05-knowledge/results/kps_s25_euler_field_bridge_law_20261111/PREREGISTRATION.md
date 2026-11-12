# Pre-registration (kind-pasteur-2026-11-11-S25)

Owner directive (2026-11-11): continue the session on Stein modules of Keller pairs (S24), and extend it into a
session that aims to prove the Chan-Pak conjecture. Both conjectures remain OPEN at the time of writing.

## Already done before this file (disclosed; not predictions)

- **Stein thread.** `scripts/euler_field.py` was written and run (`runs/euler_field.txt`): the Euler-field lemma
  `W h = (xi + delta) h` modulo `Im(d/dx - W_x) + Im(d/dy - W_y)` for a vector field with `xi(W) = W` and constant
  divergence `delta`, and its consequence that the image for Broughton's polynomial `W = x + x^2 y` is not a
  Mathieu-Zhao subspace (`W^m` in the image for all `m >= 1`, `x W^m = m! x` not in it): HYP-9668 (E') is false.
- **Chan-Pak, sweep.** A sub-agent digest of the repo (2026-11-11, after commit f90b22ab5d) found that the pieces I
  had re-derived from the odd-cycle formula are KNOWN: the vertex-addition identity (THM-4876), the join law of two
  rooted tournaments `H = x_h y_h + x_R y_R - x_e y_e` (THM-4866 (B); Session263, Theorem RJ2), the one-apex function
  `g` (THM-4875, THM-4890), the two-parameter family `(1+u)(1+v) + g(u) g(v) - e(u) e(v)` with its mass
  `(pi^2/24) X ln X` (THM-4885 (B)), and completeness of "bank x one-apex" joins up to `2^35` (THM-4885 (D)).
  A scratch run of mine (not in the repo) had recomputed the second family below `2^20` (10,895 odd numbers missed,
  the last dyadic block at 1.5 per cent); it is a re-computation of a known family, not a result.
- Concurrent: opus-2026-11-11-S93 (counting exponent 3/4 for one-apex counts, two-apex words complete below `2^23`),
  opus-2026-11-11-S94 (proof-oriented, just started), Sessions 275 and 276 (digit banks, corridors).

## What is tested next (not yet run)

**Object.** For `p >= 0` let `L(p)` be the bit length of `p` (`L(0) = 0`), `R(p)` the reversal of the complement of
`p` on `L(p)` bits (`R(0) = 0`), and `g` the one-apex function. For `i >= 1` put

    M(p, i) = [[ g(p), 2^i (R(p) + 1) ], [ p + 1, 2^(i + L(p)) ]].

Working claim (to be proved in the note from the KNOWN join law): for one-apex tournaments `Y_0, ..., Y_k` with apex
words of values `p_0, ..., p_k`, the tournament `Y_0 => C_1 => Y_1 => ... => C_k => Y_k` (`C_j` a transitive chain of
`i_j >= 1` vertices) with the arcs between consecutive apexes reversed has

    H = top-left entry of M(p_k, i_k) ... M(p_1, i_1) M(p_0, 1).

So the path-count spectrum contains the set `Sigma` of top-left entries of all products of the matrices `M(p, i)`.

**Predictions.**

- **C1.** Brute force (dynamic programming over subsets, at most 13 vertices): the formula above agrees with the
  number of Hamiltonian paths for every chain of depth 2 and 3 with words of length at most 2 and chain lengths 1,
  2 (0 mismatches).
- **C2.** Depth at most 2 (products of two matrices) misses more than 2 per cent of the odd numbers of
  `[2^17, 2^18)`.
- **C3.** Depth at most 3 misses fewer than 0.5 per cent of the odd numbers of `[2^17, 2^18)`, and fewer than 600
  odd numbers below `2^18` in all.
- **C4.** Full closure (all depths), computed exactly for all odd `N <= 2^14`: the only odd numbers missed are 7 and
  21. (If this fails, the list of misses is the result.)
- **C5.** The number of products of depth at most 3 with top-left entry at most `2^18` (counted with multiplicity
  as tuples) is at least `40 * 2^17`: mean multiplicity at least 40 per odd number.

**Reading.** C4 true would give a bank-free sharpening of the known finite completeness: the one-apex chain
semigroup alone, with no bank of small tournaments, reaches every odd number in the range. It would not prove
anything about all `N`. C4 false would show that the pure chains need the bank (or the joins without an inserted
chain), and the misses would say which.

---

## Outcomes (added after the runs; nothing above was changed)

| prediction | outcome |
|---|---|
| C1 (brute force of the chain law) | **CONFIRMED** (`runs/c1_bruteforce.txt`): 1470 chains of depth 2 and 3, 0 mismatches |
| C2 (depth at most 2 misses more than 2 per cent in `[2^17, 2^18)`) | **CONFIRMED** (`runs/depth_18.txt`): 11.06 per cent |
| C3 (depth at most 3: fewer than 0.5 per cent there, fewer than 600 below `2^18`) | **FAILED**: 5.26 per cent; 9190 misses below `2^18` |
| C4 (full closure below `2^14` misses only 7 and 21) | **FAILED** (`runs/closure_14.txt`): 1214 odd numbers missed; 414 below `2^12` (`runs/closure_12.txt`) |
| C5 (mean multiplicity at depth at most 3 at least 40) | **CONFIRMED**: 82.9 |

Notes.

1. C3 and C4 failed by a wide margin. The chain semigroup with a transitive chain inserted at every join is not a
   complete family in this range; the reading is in the README, section B3.
2. The closure run at `2^12` was made first, to gauge the cost; it is reported with the pre-registered run at
   `2^14`.
3. The numbering in the file names of this packet was fixed after the registration: the stubs first carried
   THM-5145, THM-5146 and HYP-9675, which mac-mini-2026-11-11-rays had reserved one commit earlier; they were
   renumbered THM-5147, THM-5148, HYP-9677 in the next commit.

---

## Addendum D (registered after the outcomes of C1-C5 and before the runs below)

The failure of C4 does not say which ingredient is missing: the joins without an inserted chain, or the bank of
small rooted tournaments. Test: the semigroup of **direct** joins of one-apex blocks rooted at their apexes (chains
of any length, zero included; the root moves to the new apex), which needs the parity split of the root paths.

Working claim. State of a rooted tournament: `(A, o, e)` with `A = H`, `o` (`e`) the sum of `H(X - P)` over the
paths `P` starting at the root with an odd (even) number of vertices. Joining a one-apex block with apex word of
length `n`, values `P` (letters 1, first spine vertex = top bit) and `Q` (letters 0, first spine vertex = unit bit),
with the arc between the old root and the new apex reversed, gives

    A' = A g + 2 (o E + e O),     o' = A (1 + floor(P/2)) + o t_o + e t_e,     e' = A ceil(P/2) + o t_e + e t_o,

with `O = 1 + floor(Q/2)`, `E = ceil(Q/2)`, `g = g(P)`, `t_e = t_o = 2^(n-1)` for `n >= 1` and `t_e = 1`, `t_o = 0`
for `n = 0`; inserting a chain of `i >= 1` vertices first replaces `(o, e)` by `(2^(i-1)(o + e), 2^(i-1)(o + e))`.

- **D1.** Brute force: the law agrees with the number of Hamiltonian paths on all chains of depth 2 and 3 with
  words of length at most 2 and chain lengths 0, 1 (0 mismatches).
- **D2.** The full closure of direct joins below `2^13` misses exactly 7 and 21. (My confidence is low, about 40
  per cent: 13 is not a join of two one-apex tournaments, and I do not know whether three blocks give it.)
- **D3.** If D2 fails, 13 is among the misses, and the number of misses below `2^13` is at most 50.

Reading. D2 true: the bank is not needed in this range, only the direct joins are; the family of apex-path
tournaments would be a bank-free candidate with growing multiplicity. D2 false with few misses: the bank repairs a
thin exceptional set. D2 false with many misses: roots other than apexes are essential.

---

## Outcomes of addendum D (added after the runs; nothing above was changed)

| prediction | outcome |
|---|---|
| D1 (brute force of the direct-join law) | **CONFIRMED** (`runs/d1_bruteforce.txt`): 1470 chains, 0 mismatches |
| D2 (the closure below `2^13` misses exactly 7 and 21) | **FAILED** (`runs/direct_closure_13.txt`): 29 odd numbers are missed |
| D3 (13 is missed; at most 50 misses below `2^13`) | **CONFIRMED**: 13 is missed; 29 misses |

The 29 misses: 7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 299, 421, 449, 503, 533, 569,
583, 637, 821, 1067, 1339, 1861, 2647. Not predicted: the list does not grow between `2^13` and `2^16`
(`runs/direct_closure_16.txt`: the same 29 numbers, 4,225,713 live states); a run at `2^18` was started after
these outcomes were known.

## Addendum E (registered after the outcomes of D and before the runs below)

Question: is a finite set of blocks enough? Test: the closure below `2^16` of the semigroup generated by the
one-apex blocks whose apex word has length at most `L`, for `L = 0, 1, 2, 3, 4` (1, 3, 7, 15, 31 blocks; the empty
word is the single-vertex block), in the coordinates `(A, W, d)` with `W = o + e`, `d = o - e`:
`A' = g A + (Q + 1) W - [Q even] d`, `W' = (P + 1) A + 2^n W`, `d' = [P even] A - [n = 0] d`, from the virtual
state `(1, 0, 0)`.

- **E1.** For `L <= 1` the misses do not stop: more than 1000 odd numbers are missed in `[2^15, 2^16)`.
- **E2.** For `L = 3` the misses stop: no odd number is missed in `[2^13, 2^16)`.
- **E3.** For `L = 4` the set of misses below `2^16` is a superset of the 29 numbers above with at most 60 elements.

Reading. E2 true would give a finitely generated semigroup of 3 x 3 integer matrices whose first coordinates cover
the odd numbers of a long range: a local-to-global statement about a finitely generated semigroup, which is a
smaller target than the infinite family.

---

## Outcomes of addendum E (added after the runs; nothing above was changed)

`runs/finite_blocks_16.txt` (closure below `2^16`, `L = 0..4`), `runs/finite_blocks_trend_17.txt` (below `2^17`,
`L = 3, 4, 5`; the value `L = 5` was not in the registration and was run after the outcomes for `L <= 4`).

| prediction | outcome |
|---|---|
| E1 (`L <= 1`: more than 1000 misses in `[2^15, 2^16)`) | **CONFIRMED**: 14,323 for `L = 1` (16,383 for `L = 0`) |
| E2 (`L = 3`: no miss in `[2^13, 2^16)`) | **FAILED**: 712 misses there (285 in `[2^15, 2^16)`) |
| E3 (`L = 4`: a superset of the 29 numbers with at most 60 elements below `2^16`) | **FAILED** on the count: 126 misses (the 29 numbers are among them, as they must be) |

Not predicted: with words of length at most 5 (63 blocks) the closure below `2^17` misses 49 odd numbers, the
largest being 10391, and none in `(10391, 2^17]`; with `L = 4` there are 127 misses below `2^17`, the largest
130247; with `L = 3`, 1355. A run for `L = 5` below `2^19` was started after these numbers were known.

Note on addendum D: the run at `2^18` (`runs/direct_closure_18.txt`, 32,860,481 live states) gives the same 29
misses.
