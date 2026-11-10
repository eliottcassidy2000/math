# The 2-adic lines of the Hamiltonian-path spectrum: unit-offset rays, the isolation of 7 and 21, and the extension-mean law

**Session:** kind-pasteur-2026-11-10-S21. **Canon:** THM-5126, HYP-9658. **Status:** PROVED / FINITE-EXACT partial
results; two pre-registered runs (order nine: 2 of 4 predictions confirmed; order ten: 4 of 4); independent audit
in section 9. **Headline:** every odd `N > 1` whose `N - 1` has odd part at most 129 is a Hamiltonian-path count,
except 7 and 21 (65 infinite lines; a set of density zero).
**The Chan-Pak conjecture (every odd number other than 7 and 21 is a Hamiltonian-path count of a tournament) remains
OPEN.** No coverage cutoff, density bound or completeness statement is claimed.

**Owner prompt (summary).** Aim to prove that 7 and 21 are the only missing odd path counts; not more numerical
evidence but stronger reasoning, structural experiments, arguments about infinite objects; connect to planar
graphs and groups, A000568, Hadwiger-Nelson, the planar Jacobian conjecture; search at random when stuck; name
the principles.

## 0. Inheritance, collisions, concept board

- **Closest proved mechanism.** The fixed-root tail law (O/Z tails: Sessions 195-205, THM-4888 (B), THM-5085,
  Session236 offsets): appending vertices that beat one root gives `c0 + M 2^t`.
- **Canonical hostile.** `11 in S` but `21 = 2 * 11 - 1 not in S`: no realiser of 11 can be doubled.
- **Corrected near miss.** MISTAKE-750/751 (this machine re-derived its own tail results twice). This time the
  notes of THM-4879, THM-5085, THM-5109 were opened first; the tail law and dyadic cones are credited, not claimed.
- **Least-used sidecar.** The offset `c0` itself. Earlier work used rays for their slopes and windows; here the
  offset is the organising coordinate.
- **Same-prompt sessions today.** opus-2026-11-10-S90 (THM-5124 additive rigidity on presence cubes, THM-5125,
  HYP-9657) and Session246 (THM-5120-5123) pushed before this note. Three things first derived here were already
  theirs and are credited: strict positivity of the cut matrix (THM-5122), "no additive pair" (THM-5124), the
  weighted-pair form of one-apex counts (S90 section 7, Session 16). Sessions 238-245 cover planar groups, sequence
  engines, sevenfold symmetry and gap anchoring.

Concept board: (1) rays `c0 + M 2^t`; (2) lines through the 2-adic point 1; (3) crossed covers `delta'`;
(4) doubling vertices; (5) the mean of an extension family; (6) binary weight of `N` (Horner form).

## 1. Rays and their offsets (THM-5126 (A))

For a rooted tournament `(L, v)` let `d = H(L - v)`, `E` the number of Hamiltonian paths ending at `v`, and
`delta'` the number of *crossed covers*: ordered pairs of nonempty disjoint paths `(P1, v P2')` covering `L`, with
`P2'` nonempty and `v` beating the last vertex of `P1`. Then for `t >= 1`

    H(L o 0^t) = c0 + M 2^t,     c0 = d - delta' = H(B) - W(O),     M = H - E + delta' = U_B(O),

with `B = L - v`, `O = N+(v)`, `W(O)` the sum over the pairs `{e, f}` of `O` of the two-path-cover coupling `W_ef`
(the cut matrix of THM-4114 and THM-5122), and `U_B(O)` the sum of the responses `U_B(y)`, `y in O`.

The law is known; the identities are the new form. They say what the offset *is*: the path count of the base minus
the coupling energy of the out-neighbourhood. A ray is thus a pair (base, vertex set): its offset is a quadratic
function of the set and its slope a linear one. The offset is the 2-adic limit of the ray: the tournament with an
infinite tail has a 2-adic path count, and it is an ordinary odd integer `<= H(B)`.

## 2. Lines (THM-5126 (B))

Every odd `N > 1` is `1 + c 2^s` for one odd `c`: it lies on the line `L_c`. Chan-Pak says that every line is in
the spectrum except for one point on `L_3` and one on `L_5`. A root with `c0 = 1` (a *unit-offset root*) puts a
whole ray of one line into the spectrum. By section 1 this is the single equation `W(O) = H(B) - 1`. So does a
*doubling vertex* (a vertex `x` of `K` with `(H - 1)/2` double
visits and no triple visit: then `H(K[TT_a at x]) = 1 + (H - 1) 2^(a-1)`).

**The two gaps are isolated on their lines.**

| line | points | seed | in the spectrum |
|---|---|---|---|
| `L_3` | 7, 13, 25, 49, ... | strong 4-tournament rooted at a one-triangle vertex, `(H, d, E, U) = (5, 3, 1, 9)`, `M = 6` | all but 7 |
| `L_5` | 11, 21, 41, 81, ... | an order-6 root with `(15, 9, 3, 29)`, `M = 20`; and 11 at order 5 | all but 21 |

So 7 is the missing *root* of its line, and 21 is the one point between the root 11 and the first seed 41. The
hostile fact "11 cannot be doubled" is the same statement.

**All 65 lines with `c <= 129` are complete.** Every odd `c <= 129` has a seed: 54 of order `<= 9` and, for
`c = 19, 33, 71, 79, 83, 97, 101, 103, 107, 113, 115`, eleven of order ten found through the equation
`W(O) = H(B) - 1`. The 119 points below the rays are path counts of order `<= 9`, except `5313 = 3 * 1771` and
`3425 = 5 * 685`, which are order joins. Hence **every odd `N > 1` whose `N - 1` has odd part at most 129 is a
path count, except 7 and 21** (THM-5126 (B3); seeds in `line_seeds.txt`, certificates in `below_rays.txt`).
Beyond: 2,016 of the 2,049 odd `c <= 4097` have a seed of order `<= 10`; the 33 without are 245, 853, 909, 1049,
1407, 1467, 1515, 1577, 1629, 1829, 1887, 2023, 2037, 2085, 2343, 2393, 2703, 3001, 3033, 3123, 3313, 3411, 3427,
3545, 3713, 3733, 3753, 3773, 3907, 3909, 3985, 3991, 4023 (HYP-9658 (B)).

**Weight three.** `7 = 111_2` and `21 = 10101_2` both have binary weight three, and weight-three numbers are the
points of the lines `c = 2^d + 1`. Rows `d <= 7` are complete except for 7 and 21 (THM-5126 (B4)). In Horner form `N = 1 + 2^(a_1)(1 + 2^(a_2)(1 + ...))`: the gaps are the two
depth-two Horner numbers with equal small steps, `1 + 2(1 + 2)` and `1 + 4(1 + 4)`; the next one,
`73 = 1 + 8(1 + 8)`, is a path count.

## 3. The extension-mean law (THM-5126 (C))

The mean of `H` over all extensions of `T` by `k` new vertices is `k! 2^(1-k) sum_j C(k+1, j) 2^(-j) P_j(T)`, and
`P_j(T) >= j! S(n, j)` with equality only for transitive `T`. For one new vertex: mean `= H + P_2/4 >= 2^(n-2) + 1/2`.

Reading: a one-vertex extension family has `2^n` members and mean at least `2^(n-2)`. It is never supercritical,
and the transitive base (the one-apex family) is the densest. This is the exact form of the observation that
one-apex counts cover a positive proportion but not everything. For `k` apexes over a chain the same identity
gives the word-growth exponent `k / log2(k + 1)`: `1` at `k = 1`, above `1` from `k = 2` (Session246's flag
eigenvalues `(k+1)^l` are the same growth). The value 13, the first number missed by the one-apex family, is the
first point of the line `L_3` above its gap.

## 4. Experiments (universes, controls, outcomes)

| experiment | universe | result |
|---|---|---|
| unit-offset roots (`seeds.c`) | every root of every class of order `<= 9` (1,723,824 at order 9) | 2, 5, 7, 16, 31, 68, 207, 1333 roots; 636 odd `c` with a seed |
| unit sets (`unitsets.c`) | every vertex set of size 2-4 over every class of order 9 (roots of order 10 with out-degree 2-4) | 19,002 unit sets (2,771 pairs, 16,231 triples, no 4-set); 5,529 odd `c` with a seed in total |
| doubling vertices (`dbl.c`) | every vertex of every class of order `<= 8` | 49 of the 320 values have a doubling realiser; for `N = 3 mod 4`, `N < 200`: 3, 19, 23, 31, 59, 75, 95, 99, 103, 111, 155, 163, 179, 191, 199 |
| two-module families (`modfam.c`) | every class of order `<= 8`, every ordered pair of vertices (406,412 triples) | 3,477 distinct exact bilinear tuples at order 8; none is additive (`de = 0` with `be, ga != 0`), in agreement with THM-5124 |
| weight-three map (`cover3.py`) | all tuples above, `a, b <= 26` | rows `d <= 4` and `d = 7` complete from order `<= 8` cores alone |
| pre-registered order-9 run | PREREGISTRATION.md, pushed before the run | P2, P3 CONFIRMED; P1 (count in [250, 420]; actual 636) and P4 (`c = 33`) FAILED; gate G1 holds |
| pre-registered order-10 run | second part of PREREGISTRATION.md, pushed before the run | Q1-Q4 CONFIRMED (all eleven missing `c <= 129` acquire seeds; 510 of 513 odd `c <= 1025`); gate G2 holds |

Positive controls: the scan recovers the fans (`c = 1`), the family `3 - X - Y + XY` of `C_3[TT_a, TT_b, 1]`, and
the known one-apex ray `1 + 9 * 2^k` (word `1 0 1^k 0 0`). Hostile controls: no ray passes through 7 or 21; the
false claim "every cycle vertex is doubling" is rejected (verify_s21.py). The strong doubling property fails
widely: 15 has no doubling realiser of order `<= 8`, so "double every realiser" is not a route.

## 5. What this does and does not do for the proof

- It reorganises the target. Chan-Pak is equivalent to: *every line `L_c` is in the spectrum, with one exception
  on `L_3` and one on `L_5`*. A ray settles a line from one point on; what remains on each line is finite.
- It does not reduce the problem to finitely many checks: there are infinitely many lines, each needs a seed
  (HYP-9658 (A)) and its finitely many points below the seed. The points below a seed are again ordinary targets.
- It gives the first statement of this thread that is complete on an infinite set defined by a simple arithmetic
  condition: all `N` with `N - 1 = c 2^s`, `c <= 129`. The set has density zero (about `65 log2 X` members up to
  `X`), so this says nothing about density.
- It turns "find a seed" into an equation in known objects: `W(O) = H(B) - 1` for the cut matrix `W` of the base.
  The order-ten run solved it for every `c <= 129` in 22 seconds; a construction of solutions for all `c` is open.
- It removes one architecture: lifting every realiser along its own line (the 15 control), and it adds a concrete
  gate for any engine: it must miss the root of `L_3` and the second point of `L_5` and nothing else on them.
- Rigorous lower bounds for the counting function of the spectrum are unchanged (THM-4887: `X^(2/3)` up to a
  logarithm). Positive density is still unproved.

Next obligations: (1) a grammar for the unit sets `W(O) = H(B) - 1`, to produce seeds for all `c` rather than find
them; (2) order eleven for HYP-9658 (B); (3) decide whether `s0(c)` is bounded on average; (4) the points below
the rays for `129 < c <= 4097` (eight of them for `c <= 257` are above the order-nine range).

## 6. Connections, typed (source / target / map / preserved / lost / test)

1. **Planar groups (Session238, THM-5096) - EXACT, inherited.** Group closure of the physical port moves reaches
   7 and 21; the positive semigroup does not. Here: a ray is a positive word `Z^t`; its offset `c0` is what a
   formal inverse `Z^(-infinity)` would reach. Unit offset means the ray points at the 2-adic integer 1, which is
   a path count only as the empty tournament. Preserved: the affine law. Lost: positivity below `t = 1`.
2. **A000568 - EXACT identity (small).** The mean of the path-cover polynomial `pi_T(x) = sum_j P_j(T) C(x, j)`
   over all labelled tournaments on `n` vertices is `sum over permutations of odd order of 2^(psi - moved) x^cyc`
   (`psi` nontrivial cycles, `moved` non-fixed points); at `x = 1` it is Szele's `n!/2^(n-1)`. Davis's formula for
   A000568 sums `2^(orbits on pairs)` over the same permutations. Map: permutations of odd order. Preserved: the
   index set. Lost: everything else (one counts classes, the other paths; S239 says the same). FINITE-EXACT:
   `scripts/mean_pi_identity.py` checks it on every labelled tournament of orders 3-5 at `x = 1, 2, 3`. A proof
   would follow from the Grinberg-Stanley expansion of the path-cover polynomial over odd-cycle packings; it is
   not carried out here, so the identity is not claimed for all `n`.
3. **Collatz thread - EXACT restatement.** THM-5126 (D): the one-apex count is `1 + sum 2 floor((o+1)/4)` over the
   trajectory of the accelerated map `o -> (o+1)/2^v`. No dynamics is imported: this map always terminates.
4. **Hadwiger-Nelson and planar Jacobian - ANALOGY of kind only.** All three are spectrum problems for a positive
   construction semigroup: small exclusions by finite obstruction (7 and 21; the Moser-blocker floors of HN.md;
   the excluded small Keller degree pairs), the statement at infinity open. No map between the objects.
   Test that would upgrade it: a common exact criterion for "a seed exists"; none is known.
5. **Projective planes - COINCIDENCE.** `7 = 1 + 2 + 4` and `21 = 1 + 4 + 16` are `q^2 + q + 1` for `q = 2, 4`,
   the orders of `PG(2,2)` and `PG(2,4)`; also the numbers of vertices and arcs of a 7-vertex tournament. But
   `73 = 1 + 8 + 64` is a path count. Recorded so that it is not rediscovered as a bridge (see also HYP-9362).

## 7. Principles (names for this note; not promoted to META-PATTERNS)

- **Read the offset, not the slope.** Trigger: an affine family `c0 + M b^t`. Action: compute `c0` as a count
  (`d - delta'` here) and sort families by it. Counterindication: `c0` says nothing about points below the ray.
- **A gap is a missing prefix of a line.** Trigger: a sporadic exception. Action: find the line through it and
  the first seed above it; the exception should be isolated on the line (7, 21 are).
- **The mean prices a family.** Trigger: a construction with `2^m` members. Action: compute the exact mean of the
  count; if it is at least `2^m / C` the family cannot cover an interval by itself. (THM-5126 (C).)
- Used from S90's list: "Additivity is orientation, never presence"; from the cards: "Search the statement before
  the method" (three collisions found before writing), "Existence is a maximum or tail question".

## 8. Random search log (owner's instruction)

Seeded draw of 24 canon titles (seed `2026111021`). Read: THM-4518 (three thieves: mean `O(1)`, second moment
`log`, so hitting probability `1/log`), THM-1126 (half-coverage gate), THM-157 (alpha_1 linear identity).
Transfer kept: THM-4518's shape - a family whose mean is bounded can still fail to hit because hits cluster; the
honest test of any "expected multiplicity" argument for coverage is a second moment, which for these families is
the open linear-energy question (THM-4887). Nothing else transferred.

## 9. Audit

See `audit/AUDIT.md` (independent agent, fresh context).

## 10. Reproduction

    cd 05-knowledge/results/kps_s21_two_adic_lines_20261110/scripts
    python verify_s21.py                 # about 20 s; exits 1 on failure
    python verify_lines.py 129           # reads ../line_seeds.txt; about 5 s
    python mean_pi_identity.py 5         # the odd-order mean identity of section 6, orders 3-5
    # full regeneration (needs classgen.c from 05-knowledge/results/first_creators_genealogy_20261107/):
    gcc -O2 -o classgen classgen.c && ./classgen 8 cls && ./classgen 9 cls9      # cls_1..8.txt, cls9_9.txt
    gcc -O2 -o seeds seeds.c && for m in 2 3 4 5 6 7 8; do ./seeds $m < cls_$m.txt > seeds_$m.txt; done
    ./seeds 9 < cls9_9.txt > seeds_9.txt
    gcc -O2 -o dbl dbl.c && for m in 3 4 5 6 7 8; do ./dbl $m < cls_$m.txt > dbl_$m.txt; done
    gcc -O2 -o unitsets unitsets.c && ./unitsets 9 4 < cls9_9.txt > us_10.txt   # order-ten unit sets, 22 s
    python finalize.py 129               # writes line_seeds.txt, below_rays.txt and prints the line table
    python unitpairs.py ../unit_offset_roots_orders_2_to_9.txt 8                # identity (1') on 336 roots
    gcc -O2 -o modfam modfam.c && ./modfam 8 3 < cls_8.txt > mf_8.txt           # about 4 minutes
