---
id: THM-4790
title: "Presentation-count fibres 7 and 21 are absent for every tournament on at most 11 vertices. (a) Reduction, all n: a tournament with tc = H/|Aut| in {7, 21} has a strong component S with tc(S) in {7, 21}, and then |Aut(S)| >= max(3, f(m)/21), where f(m) is the minimum number of Hamiltonian paths of a strong m-tournament (Moon 1972, Busch 2006). (b) Automorphism-type lemma, n <= 11: an odd-order subgroup of S_n of order at least L_n (L_7 = L_8 = 3, L_9 = 5, L_10 = 7, L_11 = 11) contains an element of one of the listed cycle types. (c) Every tournament on n <= 6 vertices, and every strong tournament on 7 <= n <= 11 vertices invariant under a representative of each listed type, was checked (non-strong ones follow from (a)): none has tc in {7, 21}. Hence HYP-9352 (2) holds through n = 11 (the THM-4741 census had n <= 9). In the graph chart, tc(T) is the number of labelled graphs G on [n-1] whose path tournament T_G is isomorphic to T, so the path-presentation map from graphs on [n-1] to tournament classes has no fibre of size 7 or 21 for n <= 11."
status: >
  PROVED: (a) from THM-2450 sec. 2 (stack multiplicativity), LEM-003 (Aut acts freely on Hamiltonian paths, so tc is an
  odd integer) and the strong minimum f(m) (KNOWN, Moon 1972; Busch 2006; written out in THM-4744); (b) Cauchy, the
  structure of groups of order 9 and 25, the Sylow 3-subgroup of S_11, Sylow's counting theorem and Jordan's theorem
  (proof below).
  FINITE-EXACT: (c), tc_fibres_upto11_20261026.py (numba; self-tests against C3, T5, P7 and brute force on 60 random
  tournaments; Moon/Busch minimum reproduced exhaustively for n <= 6; 256 s).
  INDEPENDENTLY AUDITED (05-knowledge/results/s57_untouchables_root27_20261026/AUDIT.md, MISTAKE-672; no error):
  - every case of (b) re-checked, including the full subgroup lattice of Z_3 wr Z_3;
  - an independent re-implementation (different conjugates, pull Held-Karp, refined-invariant automorphism counter,
    no candidate filter) reproduced every count in (c), and the two codes agree on 25,674 strong invariant
    tournaments with nontrivial Aut;
  - runs bypassing (b) found no hit: n = 9 type (3) (2^22), and n = 11 types (5) (2^23) and (3)(3) (2^25). So
    n <= 11 also holds using only Cauchy and the order-9 argument.
  CITED: THM-1370 (the same reduction for H, with the insertion lemma; H never 7 or 21 for any n), MathOverflow
  232751 and Grinberg's notes (the H version).
  OPEN: n >= 12 (the 3-group case needs the type (3)(3)(3) on 12 points, 2^24 invariant tournaments); HYP-9352 (2) in general.
session: opus-2026-10-26-S57
source: 05-knowledge/results/s57_untouchables_root27_20261026/README.md
depends_on:
  - 01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md
  - 01-canon/theorems/THM-2450-rigid-self-converse-classes-are-cyclic-ternary-towers.md
  - 01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md
related:
  - 01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
  - 05-knowledge/hypotheses/HYP-9352-presentation-count-laws-strong-spectrum-forbidden-fibre-sizes-and-the-three-move-law-at-n-11.md
scripts:
  - 04-computation/untouchables_s57/tc_fibres_upto11_20261026.py (+ tc_fibres_upto11.json)
  - 04-computation/untouchables_s57/tc_fibres_n10_20261026.py (+ tc_fibres_n10.json; first, pure-Python run of the n = 10 case)
---

# THM-4790: no presentation fibre of size 7 or 21 through eleven vertices

**Setting.** For a tournament `T` on `n` vertices, `H(T)` is its number of Hamiltonian paths and
`tc(T) = H(T)/|Aut(T)|`. By LEM-003 the automorphism group acts freely on Hamiltonian paths, so `tc(T)` is an integer;
it is odd because `H` (Rédei) and `|Aut|` are odd. THM-4741 identifies `tc(T)` with the number of marked-path
presentations. Concretely, a marked Hamiltonian path `v_1 -> ... -> v_n` turns `T` into the **path tournament** `T_G`
of a labelled graph `G` on `[n-1]`: the non-consecutive pair `{v_i, v_(j+1)}` (`i < j`) is backward iff `{i, j}` is an
edge of `G`. Every graph on `[n-1]` occurs, and two marked paths give the same `G` iff an automorphism carries one to
the other. So

    tc(T) = #{ G a labelled graph on [n-1] : T_G is isomorphic to T },     sum over classes T of tc(T) = 2^binom(n-1,2).

The pair `(G, complement of G)` is the pair of not necessarily connected graphs (backward arcs, forward non-path arcs)
that a marked path assigns to `T`. An isomorphism class receives exactly `tc(T)` such ordered pairs. HYP-9352 (2)
says no class receives 7 or 21 of them. THM-1370 says `H` itself is never 7 or 21.

**(a) Reduction.** Let the strong components of `T` be `S_1 => ... => S_r`. Every Hamiltonian path traverses them in
order, and every automorphism fixes each component (the condensation order is total), so `H(T) = prod H(S_i)` and
`Aut(T) = prod Aut(S_i)` (THM-2450 sec. 2). Hence `tc(T) = prod tc(S_i)`. Each factor is an odd positive integer, so
`tc(T) = 7` forces a factor 7, and `tc(T) = 21 = 3 * 7` forces a factor 7 or 21. A minimal counterexample is therefore
strong. For a strong `S` on `m` vertices, `H(S) >= f(m)`, with `f = 3, 5, 9, 15, 25, 45, 75, 125, 225` at
`m = 3, ..., 11` (Moon 1972; Busch 2006; THM-4744). If `|Aut(S)| = 1`, then `tc = H >= 25 > 21` for `m >= 7`.
Otherwise `|Aut(S)| >= 3` and `|Aut(S)| = H/tc >= f(m)/21`. The resulting odd lower bounds are
`L_7 = L_8 = 3`, `L_9 = 5`, `L_10 = 7`, `L_11 = 11`.

**(b) Automorphism-type lemma.** Let `G <= S_n` have odd order `|G| >= L_n`. Then `G` contains an element whose cycle
type (fixed points omitted) is in the list

| `n` | types |
|---|---|
| 7, 8 | `(3)`, `(3)(3)`, `(5)`, `(7)` |
| 9 | `(3)(3)`, `(3)(3)(3)`, `(5)`, `(7)` |
| 10 | `(3)(3)`, `(3)(3)(3)`, `(5)`, `(5)(5)`, `(7)` |
| 11 | `(3)(3)(3)`, `(5)(5)`, `(5)(3)`, `(5)(3)(3)`, `(7)`, `(11)` |

*Proof.*
- If `p` in `{5, 7, 11}` divides `|G|`, Cauchy's theorem gives an element of order `p`, of type `(p)^j`; on
  `n <= 11` points these are `(5)`, `(5)(5)`, `(7)`, `(11)`. For `n <= 8` every nontrivial odd-order group has an
  element of prime order 3, 5 or 7, and the order-3 types on at most 8 points are `(3)` and `(3)(3)`.
- For `n = 9, 10` the remaining case is a 3-group. Since `L_n >= 5` it has order at least 9 and contains a subgroup
  `K` of order 9. If `K` is cyclic, its generator is a 9-cycle, whose cube has type `(3)(3)(3)`. If `K` is `Z_3^2`
  and every nonidentity element were a single 3-cycle, take two generators: commuting 3-cycles that are not powers of
  one another are disjoint, so their product has type `(3)(3)`, a contradiction.
- For `n = 11` the 3-groups have order 27 or 81 (`L_11 = 11`, and the 3-part of `11!` is `3^4`). Each lies in a
  conjugate of the Sylow subgroup `P = Z_3 wr Z_3` acting on 9 points. A subgroup of order 27 of `P` either is the base
  group, which contains `(3)(3)(3)`, or contains an element `(b, t)` with `t != 1`. Such an element moves all nine
  points, so it has type `(3)(3)(3)` or `(9)`, and the cube of a 9-cycle has type `(3)(3)(3)`. The case of order 81
  is `P` itself.
- For `n = 11`, if `25 | |G|`, then a subgroup of order 25 is `Z_5^2` (no 25-cycle exists on 11 points). Two
  generators that are single 5-cycles commute without being powers of one another, so they are disjoint and their
  product has type `(5)(5)`.
- For `n = 11` with `5 || |G|` and `7, 11` not dividing `|G|`, `|G| = 5 * 3^j` with `j >= 1`. The number of Sylow
  5-subgroups divides `3^j` and is `1 (mod 5)`, so it is 1 or 81.
  - If it is 81, then `|G| = 405` and `G` contains a full Sylow 3-subgroup of `S_11`. That subgroup's 9-point orbit
    extends to a `G`-orbit whose size divides 405 and lies in `[9, 11]`, hence is 9; the other two points are fixed.
    So `G` acts faithfully and transitively on 9 points and contains a 5-cycle there.
  - A primitive group of degree 9 containing a 5-cycle contains `A_9`, by Jordan's theorem (`5 <= 9 - 3`), which has
    even order.
  - An imprimitive one has three blocks of size 3. A 5-cycle acts trivially on the set of three blocks and then
    trivially inside each block, which is impossible.
  - So the Sylow 5-subgroup is normal. A Sylow 3-subgroup acts on it through `Aut(Z_5) = Z_4`, hence trivially, and
    `G` contains an element of order 15. On 11 points its type is `(5)(3)` or `(5)(3)(3)`.
- Elements of order 7 or 11 on 11 points have types `(7)` and `(11)`. `QED`

**(c) Computation.** If a strong counterexample `S` exists, then `Aut(S)` contains an element `sigma` of a listed type.
Relabelling puts `sigma` in its standard form, and `S` is invariant under it. For each type, all `2^(#pair orbits)`
invariant tournaments were enumerated. The script keeps the strong ones, computes `H` exactly (Held–Karp), and
computes `|Aut|` by backtracking whenever `H = 7a` or `21a`, with `a >= L_n` divisible by the order of `sigma`.

| `n` | method | invariant tournaments (strong) per type | `|Aut|` computations | `tc` in `{7,21}` |
|---|---|---|---:|---|
| 3–6 | all labelled tournaments | `2^3, 2^6, 2^10, 2^15` | all | none |
| 7 | `(3), (3)(3), (5), (7)` | `2^11 (1088), 2^7 (56), 2^5 (8), 2^3 (8)` | 12 | none |
| 8 | same | `2^16 (44,640), 2^10 (528), 2^8 (96), 2^4 (0)` | 6,240 | none |
| 9 | `(3)(3), (3)(3)(3), (5), (7)` | `2^14 (10,336), 2^12 (3,760), 2^12 (2,176), 2^6 (16)` | 1,660 | none |
| 10 | `(3)(3), (3)(3)(3), (5), (5)(5), (7)` | `2^19 (391,104), 2^15 (23,520), 2^17 (89,280), 2^9 (480), 2^9 (192)` | 69,616 | none |
| 11 | `(3)(3)(3), (5)(5), (5)(3), (5)(3)(3), (7), (11)` | `2^19 (392,000), 2^11 (992), 2^13 (4,352), 2^9 (224), 2^13 (4,352), 2^5 (32)` | 57,394 | none |

The exhaustive part reproduces the strong `tc` minima `1, 5, 3, 9` at `n = 3..6` (HYP-9352 (1)) and the minimum
`f(n)` for `n <= 6`. Every strong invariant tournament has `H >= f(n)`, as it must.

**Scope.** The proof is complete for `n <= 11`; nothing is claimed for `n >= 12`. There `L_12 = 19`, and the 3-group
case needs the types `(3)(3)(3)` on 12 points (`2^24` invariant tournaments) or `(3)(3)(3)(3)` (`2^22`), which is the
natural next computation. The argument uses the automorphism group where THM-1370 uses the insertion lemma. Neither
proof explains *why* 7 and 21 are missing, and the hypothesis that every other odd number is attained (THM-1370;
HYP-9352 (2) companion) is untouched.
