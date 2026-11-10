---
id: THM-5126
title: "Unit-offset rays and the 2-adic lines of the Hamiltonian-path spectrum. (A) Offset formula: for a rooted tournament (L, v), B = L - v and O = N+(v), the fixed-root tail ray is H(L o 0^t) = c0 + M 2^t with c0 = H(B) - W(O) and M = U_B(O), where W(O) sums the two-path-cover couplings W_ef over the pairs of O and U_B(O) sums the responses. (B) Lines: every odd N > 1 lies on one line L_c = {1 + c 2^s}, c the odd part of N - 1; a unit-offset root (W(O) = H(B) - 1) or a doubling vertex puts a whole ray of a line into the spectrum. Every odd N > 1 whose N - 1 has odd part at most 129 is a Hamiltonian-path count, except 7 and 21: all 65 lines with c <= 129 lie in the spectrum, L_3 misses exactly 7 and L_5 misses exactly 21. In particular every number of binary weight three 1 + 2^b + 2^(b+d) with d <= 7 is a path count except 7 and 21. (C) Extension-mean law: the mean of H over all extensions of T by k vertices is k! 2^(1-k) sum_j C(k+1, j) 2^(-j) P_j(T), minimised exactly by transitive T. (D) Trajectory formula for one-apex counts."
status: >
  (A) PROVED: the tail law itself is KNOWN (O/Z tails of Sessions 195-205, THM-4888 (B), THM-5085 at x = 2, the
  offsets of Session236); the proof is repeated in the present coordinates; new here are the identities
  c0 = d - delta' = H(B) - W(O) and M = H - E + delta' = U_B(O), which express the offset through the cut matrix W of
  THM-4114 / THM-5122. (B) PROVED from (A) and the KNOWN module law (THM-2195) + FINITE-EXACT certificates: 65
  explicit seeds (61 rooted tournaments of order <= 10 with unit offset, 4 doubling vertices of order <= 8), each
  recounted literally (verify_lines.py), and for the 119 points below the rays the complete list of path counts of
  order <= 9 (117 points) and two order joins (5313 = 3 * 1771, 3425 = 5 * 685); the exclusions of 7 and 21 are
  THM-1370 (KNOWN). The statement is about the 65 lines with c <= 129, a set of integers of density zero; whether
  every odd c has a seed is HYP-9658 (A) (2,016 of the 2,049 odd c <= 4097 have a seed of order <= 10; the least c
  without one is 245). (C) PROVED, elementary (new as an identity in this form; the one-vertex cut form is
  THM-4114 (32) and its strict positivity is THM-5122). (D) PROVED: an exact restatement of the KNOWN one-apex
  recursion g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2) (Session 193; THM-4875 (D)). Remark (E), the module face
  of THM-5124 (1), is credited to opus-2026-11-10-S90, with a second three-line proof. Independent audit: section
  "Audit". The Chan-Pak conjecture (every odd number other than 7 and 21 is a Hamiltonian-path count of a
  tournament) remains OPEN; nothing here bounds the spectrum from below in density.
session: kind-pasteur-2026-11-10-S21
source: 05-knowledge/results/kps_s21_two_adic_lines_20261110/README.md
scripts:
  - 05-knowledge/results/kps_s21_two_adic_lines_20261110/scripts/verify_lines.py (+ verify_lines.out)
  - 05-knowledge/results/kps_s21_two_adic_lines_20261110/scripts/verify_s21.py (+ verify_s21.out)
  - 05-knowledge/results/kps_s21_two_adic_lines_20261110/scripts/seeds.c, unitsets.c, unitpairs.py, dbl.c, modfam.c, finalize.py
  - 05-knowledge/results/kps_s21_two_adic_lines_20261110/line_seeds.txt, below_rays.txt, unit_offset_roots_orders_2_to_9.txt
related:
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md (the exclusions)
  - 01-canon/theorems/THM-4888-*, THM-5085-*, THM-5109-* (fixed-root tails, offsets, dyadic cones; this machine)
  - 01-canon/theorems/THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature.md (cut form)
  - 01-canon/theorems/THM-5122-shared-core-exchange-intervals-and-negative-inertia-capacity.md (the matrix W, W_ij >= 2 H(C - i - j) > 0)
  - 01-canon/theorems/THM-5124-additive-rigidity-of-the-path-count-on-presence-cubes.md (pair lemma)
  - 01-canon/theorems/THM-2195-transitive-quotients-exactly-control-universal-substitution-products.md (module law)
  - 01-canon/theorems/THM-4875-* (one-apex counts g, the join family F1)
  - 05-knowledge/hypotheses/HYP-9658-unit-offset-seeds-of-the-two-adic-lines.md
---

# THM-5126 - unit-offset rays, the 2-adic lines of the path spectrum, and the extension-mean law

## Setting

`H(T)` is the number of Hamiltonian paths of a tournament `T`; `H` of the empty tournament is `1`. `S` is the set
of values of `H` (the spectrum). For a rooted tournament `(L, v)` put `B = L - v` and `O = N+(v)`, and:

- `H = H(L)`, `d = H(B)`, `E` = number of Hamiltonian paths of `L` ending at `v`;
- `U` = number of ordered pairs `(P1, P2)` of vertex-disjoint paths covering `L` with `P2` starting at `v`
  (`P1` may be empty);
- a cover `(P1, P2)` counted by `U` is *crossed* if `P1` and `P2 - v` are both nonempty and `v` beats the last
  vertex of `P1`; `delta'` is the number of crossed covers;
- for a tournament `B`: `Q_B(e, f)` = number of ordered pairs of disjoint nonempty paths covering `B`, the first
  ending at `e` and the second starting at `f`; `W_ef = Q_B(e, f) + Q_B(f, e)` (the cut matrix of THM-4114 and
  THM-5122); `U_B(y)` = number of ordered pairs `(P1, P2)` covering `B` with `P2` starting at `y`, `P1` possibly empty.

For a binary word `y = y_1 ... y_t`, `L o y` is `L` with a new chain `z_1 -> ... -> z_t` (`z_i -> z_j` for
`i < j`), every vertex of `B` beating every `z_j`, and `v -> z_j` iff `y_j = 1`. `C(y) = sum_j (1 - y_j) 2^(j-1)`,
`eps(y) = 1 - y_1`, and `g(y)` is the count of the one-apex tournament `{v} o y`.

## (A) The tail law and the offset formula

**Tail law (KNOWN).** `H(L o y) = H + (U - d) C(y) + d (g(y) - 1) + (d - E) eps(y)`.

*Proof.* Delete `v` from a Hamiltonian path of `L o y`. It leaves `P1` (before `v`) and `P2` (after `v`), paths of
`B => Z`, so each is a path of `B` followed by an increasing run of `Z = {z_1, ..., z_t}`. Write `Z_i = Z ∩ P_i` and
`L_i` for the part of `P_i` in `B`.

1. `Z_1` empty. If `L_2` is nonempty, `L_1 v L_2` is a Hamiltonian path of `L` not ending at `v`, followed by all of
   `Z`: `H - E` paths. If `L_2` is empty, `L_1 v` is a Hamiltonian path of `L` ending at `v` and `v -> z_1` is
   needed: `E y_1` paths.
2. `Z_1` nonempty, with largest element `z_j`. Then `z_j -> v`, so `y_j = 0`. The pair `(L_1, v L_2)` is a cover
   counted by `U`.
   - `L_2` nonempty: `U - d` covers (the `d` covers with `L_2` empty are a Hamiltonian path of `B` followed by
     `(v)`). All junctions into `Z` are automatic. The splittings `(Z_1, Z_2)` are a choice of `j` with `y_j = 0`
     and of a subset of `{z_1, ..., z_(j-1)}`: `C(y)` of them.
   - `L_2` empty: `d` choices of `L_1`, and a number `C'(y)` of admissible `(Z_1, Z_2)` depending on `y` only.

So `H(L o y) = H - E (1 - y_1) + (U - d) C(y) + d C'(y)`. For `L = {v}` this reads `g(y) = y_1 + C'(y)`. Substituting
`C'(y) = g(y) - y_1` gives the law. QED

**Rays.** `g(0^t) = 1`, `C(0^t) = 2^t - 1`, `eps = 1`, so for `t >= 1`

    H(L o 0^t) = c0 + M 2^t,      c0 = H - E + 2d - U,      M = U - d.                              (1)

**Offset formula (new form).**

    c0 = d - delta' = H(B) - sum over pairs {e, f} of O of W_ef,        M = H - E + delta' = sum over y in O of U_B(y).   (1')

In particular `c0 <= H(B)`, with equality iff `|O| <= 1`; `c0` is odd; and `c0` is the 2-adic limit of the ray.

*Proof.* Split the covers counted by `U - d` (those with `P2 != (v)`) by `P1`: `P1` empty gives the Hamiltonian
paths of `L` starting at `v`; `P1` nonempty with last vertex beating `v` gives, by concatenation, the Hamiltonian
paths of `L` in which `v` is neither first nor last; the rest are the crossed covers. So `M = (H - E) + delta'`,
and `c0 = H - E + d - M = d - delta'`. A crossed cover is `(P1, v P2')` with `P1` ending at some `e in O` and `P2'`
starting at some `f in O`, `e != f`, and `(P1, P2')` an ordered two-path cover of `B`: `delta' = sum over ordered
pairs (e, f) of O of Q_B(e, f) = sum over pairs of W_ef`. A cover with `P2 != (v)` is `(P1, v P2')` with `P2'`
starting at an out-neighbour `y` of `v` and `P1` possibly empty: `M = sum over y in O of U_B(y)`. `W_ef` is even and
positive (THM-5122), which gives the parity and the equality case. QED

For a one-apex tournament with word `u` (binary value `A(u)`), `B` is transitive, `d = 1`, and `delta'` is the
weighted number of pairs of one-bits of `A(u)` (Session 16's weighted inversions; the `q` of the S90 note), so
`c0 = 1 - delta'` and the only unit-offset one-apex roots are the fans (`A(u)` a power of two).

## (B) Lines, seeds, and the isolation of 7 and 21

For odd `c >= 1` let `L_c = {1 + c 2^s : s >= 1}`. Every odd `N > 1` lies on exactly one line (`c` is the odd part
of `N - 1`), so Chan-Pak is the statement that every line lies in `S` except for `7 = 1 + 3 * 2` on `L_3` and
`21 = 1 + 5 * 4` on `L_5`.

**Definitions.** A *unit-offset root* is a rooted tournament with `c0 = 1`; by (1') this is the condition
`sum over pairs of O of W_ef = H(B) - 1`. A *doubling vertex* of `K` is a vertex `x` such that no Hamiltonian walk
of `K` visits `x` three or more times and the number `c_2` of Hamiltonian walks visiting `x` exactly twice is
`(H(K) - 1)/2`. (A Hamiltonian walk visits every other vertex once; `c_j = H(K[I_j at x]) / j!` with `I_j` an
independent set of `j` twins. A walk visits `x` at most `min(d+(x), d-(x)) + 1` times.)

**(B1) Ray lemma.**
- A unit-offset root with `M = c 2^j` (`c` odd) gives `1 + c 2^s in S` for every `s >= j + 1`, by (1).
- A doubling vertex gives `H(K[TT_a at x]) = 1 + (H(K) - 1) 2^(a-1)` for every `a >= 1`, and more generally
  `H(K[M at x]) = H(M) + (H(K) - 1) pi_M(2) / 2` for every tournament `M`, where `pi_M(2) = 2 H(M) + P_2(M)` and
  `P_2` counts ordered pairs of disjoint nonempty paths covering `M`. This is the KNOWN module law
  `H(K[M at x]) = sum_j c_j(K, x) P_j(M)` with `c_1 = H(K)`, `c_2 = (H(K) - 1)/2`, `c_j = 0` for `j >= 3`, and
  `P_2(TT_a) = 2^a - 2`.

**(B2) The two gaps are isolated on their lines.** `L_3 ∩ S = L_3 - {7}` and `L_5 ∩ S = L_5 - {21}`.

*Proof.* `7, 21 not in S` is THM-1370.
- `L_3`: let `L` be the strong tournament of order four (arcs `0->2, 1->0, 2->1, 2->3, 3->0, 3->1`) rooted at `v = 3`.
  Then `B` is a 3-cycle, `O = {0, 1}`, `W_01 = 2 = H(B) - 1`; `(H, d, E, U) = (5, 3, 1, 9)`, `c0 = 1`, `M = 6`, so
  `1 + 3 * 2^s in S` for all `s >= 2`.
- `L_5`: let `L` have out-lists `0: [3]; 1: [0]; 2: [0, 1]; 3: [1, 2, 4]; 4: [0, 1, 2, 5]; 5: [0, 1, 2, 3]`, rooted
  at `v = 2`. Then `(H, d, E, U) = (15, 9, 3, 29)`, `c0 = 1`, `M = 20`, so `1 + 5 * 2^s in S` for all `s >= 3`; and
  `11 = H` of a tournament of order five.

The four-tuples are recomputed by literal enumeration, and the first ray points are recounted on the literal
tournaments (verify_lines.py). QED

**(B3) All lines with `c <= 129`.** For every odd `c <= 129`, `L_c - {7, 21} ⊆ S`. Equivalently: every odd `N > 1`
such that the odd part of `N - 1` is at most `129` is a Hamiltonian-path count of a tournament, except `7` and `21`.

*Proof.* `line_seeds.txt` gives, for each odd `c <= 129`, a unit-offset root of order `<= 10` with
`M = c 2^(s0 - 1)`, or a doubling vertex of a tournament of order `<= 8` with `H - 1 = c 2^s0`; by (B1) the ray
`s >= s0` is in `S`. The least `s0` found:

    c    1  3  5  7  9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65
    s0   1  2  3  5  1  1  5  1  2  8  4  5  3  6  1  2  6  5  1  3  6  2  3  1  1  1  5  1  2  2  2  4  3
    c   67 69 71 73 75 77 79 81 83 85 87 89 91 93 95 97 99 101 103 105 107 109 111 113 115 117 119 121 123 125 127 129
    s0   3  3  6  2  3  1  5  1  7  3  2  1  2  1  1  4  1   4   6   3   6   1   1   5   2   3   3   4   1   1   1   1

The `121` points `1 + c 2^s` with `s < s0` are listed in `below_rays.txt`: two are the gaps `7` and `21`; `117` are
values of `H` at order `<= 9` (complete class list, 191,536 classes at order nine); and `5313 = 3 * 1771`,
`3425 = 5 * 685` are products of two such values, hence path counts of order joins. The seeds are recounted
literally by verify_lines.py (states, offset one, first ray points); the eleven seeds of order ten
(`c = 19, 33, 71, 79, 83, 97, 101, 103, 107, 113, 115`) were found through (1') by unitsets.c. The value list is the
output of classgen.c (validated against A000568). QED

**(B4) The weight-three stratum.** An odd number has binary weight three iff it is `1 + 2^b + 2^(b+d)`, i.e. a point
of the line `c = 2^d + 1`. By (B3): for `d <= 7` and every `b >= 1`, `1 + 2^b + 2^(b+d) in S` except `7`
(`d = b = 1`) and `21` (`d = b = 2`).

**Scope.** (B3) and (B4) concern the 65 lines with `c <= 129` only; their union has density zero. Whether every odd
`c` has a seed is HYP-9658 (A); a seed gives a ray, and the points below the ray remain a finite question for each
line.

## (C) The extension-mean law

Let `T` have `n >= 1` vertices and let `P_j(T)` be the number of ordered `j`-tuples of pairwise disjoint nonempty
paths covering `T` (`P_1 = H`). Add `k >= 1` new vertices and orient every pair containing a new vertex in all
`2^(kn + C(k,2))` ways.

**Theorem.** The mean of `H` over these extensions is

    E_k(T) = k! 2^(1-k) sum_(j=1)^(k+1) C(k+1, j) 2^(-j) P_j(T).                                      (2)

Moreover `P_j(T) >= j! S(n, j)` (Stirling numbers of the second kind) for every `j`, with equality for all `j` iff
`T` is transitive. Hence, among all tournaments of order `n`, the transitive one is the unique minimiser of
`E_k` for every `k`. For `k = 1`: `sum over the 2^n extensions of H = 2^(n-2) (4 H(T) + P_2(T)) >= 2^(2n-2) + 2^(n-1)`.

*Proof.* Deleting the new vertices from a Hamiltonian path of an extension leaves an ordered `j`-tuple of disjoint
paths covering `T`, `1 <= j <= k + 1`, no two consecutive in the path. Conversely an ordered cover, an order of
the new vertices (`k!`) and a choice of `j` of the `k + 1` gaps (`C(k+1, j)`) give a sequence with `j + k - 1`
junctions, each containing a new vertex, on distinct pairs; it is a path for a fraction `2^(-(j+k-1))` of the
orientations. This is (2). For the bound: every ordered partition of `V(T)` into `j` nonempty blocks carries at
least one tuple of Hamiltonian paths of its blocks (Redei), and distinct ordered partitions give distinct covers.
Equality for all `j` forces `H(T) = P_1 = 1`, i.e. `T` transitive, and holds for transitive `T`. The coefficient
of `P_1` in (2) is positive, so the minimiser is unique. QED

**Reading.** For `k = 1` the count `H(T + v_O)` is the capacity of a directed cut (THM-4114 (32)), and (2) says its
mean over the `2^n` cuts is `H + P_2/4 >= 2^(n-2) + 1/2`: a one-vertex extension family has `2^n` members and mean
at least `2^(n-2)`, and the one-apex family (transitive base) has the least mean. For the `k`-apex tape over
`TT_n`, (2) gives `E_k ~ k! (k+1)^n 4^(-k)` against `2^(kn + C(k,2))` labelled members: word-growth exponent
`k / log2(k + 1)` (`1, 1.26, 1.5, 1.72, ...`), in agreement with the flag eigenvalues `(k+1)^l` of THM-5123.

## (D) One-apex counts along the (x + 1)-trajectory

Let `G(x) = g(x - 1)` be the count of the one-apex tournament whose word has binary value `x - 1` (`g(0) = 1`,
`g(2v+1) = g(v)`, `g(2v) = g(v) + v + (v mod 2)`). Then `G(2x) = G(x)`, and for odd `m`

    G(m) = 1 + sum over o in the trajectory of m, o != 1, of 2 floor((o + 1)/4),

where the trajectory is `m -> (m + 1) / 2^(v_2(m+1)) -> ... -> 1`.

*Proof.* `G(2x) = g(2x - 1) = g(x - 1)`. For `x >= 2`, `G(2x - 1) = g(2(x-1)) = G(x) + 2 floor(x/2)`. With
`m = 2x - 1` and `G(x) = G(odd part of x)` this is `G(m) = G(m') + 2 floor((m+1)/4)`, `m'` the next point. QED

## (E) Remark: two modules always interact

If `x != y` both lie on cycles of `K`, some Hamiltonian walk of `K` visits `x` twice and `y` twice (and every other
vertex once); hence a family `H(K[TT_a at x, TT_b at y]) = al + be 2^a + ga 2^b + de 2^(a+b)` (exact when no walk
visits `x` or `y` three times, where `de` is the number of such walks) has `de >= 1` as soon as `x` and `y` lie on
cycles. This is the module face of THM-5124 (1) (opus-2026-11-10-S90, first). Direct proof: if `x, y` are in one
strong component `Q`, take a Hamiltonian cycle of `Q` through `x` (Camion), read as a closed walk from `x` to `x`,
and put the second copy of `y` before it if `y -> x`, after it if `x -> y`; if they are in different components,
concatenate two closed Camion walks. In the complete scan of orders 3-8 (406,412 triples `(K, x, y)`) no exact
bilinear family has `de = 0` with `be, ga != 0`.

## Controls

- verify_lines.py: the tail law on 15,570 triples `(L, v, y)` by literal counting; the 65 seeds (states, `c0 = 1`,
  first ray points on the literal tournaments); the certificates below the rays; gate: no ray through 7 or 21.
- unitpairs.py: (1') on the 336 unit-offset roots of order `<= 8`. unitsets.c reproduces the sets of
  (out-degree, `M`) of the direct runs of seeds.c at orders 8 and 9 (140 and 1,007 pairs).
- verify_s21.py: (2) for `k = 1` on every labelled tournament of order `<= 5` and for `k = 2` on order `<= 4`; the
  bound `P_j >= j! S(n, j)` with its equality case; (D) for all odd `m < 2^15`; the doubling laws of (B1) on the
  2,256 doubling vertices of orders 3-5; (E) on all labelled tournaments of orders 3-5; a negative control.
- Pre-registered runs (PREREGISTRATION.md): order nine P2, P3 CONFIRMED, P1, P4 FAILED; order ten Q1-Q4 CONFIRMED;
  both gates hold.

## Audit

See the source note, section "Audit", and `05-knowledge/results/kps_s21_two_adic_lines_20261110/audit/`.
