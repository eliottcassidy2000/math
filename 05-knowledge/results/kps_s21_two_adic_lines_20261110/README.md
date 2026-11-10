# The 2-adic lines of the Hamiltonian-path spectrum: unit-offset rays, the isolation of 7 and 21, and the extension-mean law

**Session:** kind-pasteur-2026-11-10-S21. **Canon:** THM-5126, HYP-9658. **Status:** PROVED / FINITE-EXACT partial
results; two pre-registered runs (order nine: 2 of 4 predictions confirmed; order ten: 4 of 4); independent audit
PASS WITH CORRECTIONS (section 9; five majors on credit, scope and wording, all applied; MISTAKE-754). **Headline:** every odd `N > 1` whose `N - 1` has odd part at most 129 is a Hamiltonian-path count,
except 7 and 21 (65 infinite lines; a set of density zero).
**The Chan-Pak conjecture (every odd number other than 7 and 21 is a Hamiltonian-path count of a tournament) remains
OPEN.** No coverage cutoff, density bound or completeness statement is claimed.

**Owner prompt (summary).** Aim to prove that 7 and 21 are the only missing odd path counts; not more numerical
evidence but stronger reasoning, structural experiments, arguments about infinite objects; connect to planar
graphs and groups, A000568, Hadwiger-Nelson, the planar Jacobian conjecture; search at random when stuck; name
the principles.

## 0. Inheritance, collisions, concept board

- **Closest proved mechanism.** The fixed-root tail law (O/Z tails: Sessions 195-205, THM-4888 Setting, THM-5085 (3),
  Session236 offsets) and its cut form, the star-closure law of Session216 (`h - T + 2^t A`; three-core/PROOF.md
  (3)-(4), star/PROOF.md section 6; restated in Session229, THM-4960, Sessions 244-245, THM-5122 (1)): appending
  vertices that beat one root gives `c0 + M 2^t`.
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
- **Found by the audit, not by this session (MISTAKE-754).** The cut form of the offset is Session216's; the
  one-vertex case of the mean law is THM-4111 / THM-4115; the substitution law is THM-3121 (4), which this note
  first cited as THM-2195.

Concept board: (1) rays `c0 + M 2^t`; (2) lines through the 2-adic point 1; (3) crossed covers `delta'`;
(4) doubling vertices; (5) the mean of an extension family; (6) binary weight of `N` (Horner form).

## 1. Rays and their offsets (THM-5126 (A))

For a rooted tournament `(L, v)` let `d = H(L - v)`, `E` the number of Hamiltonian paths ending at `v`, and
`delta'` the number of *crossed covers*: ordered pairs of nonempty disjoint paths `(P1, v P2')` covering `L`, with
`P2'` nonempty and `v` beating the last vertex of `P1`. Then for `t >= 1`

    H(L o 0^t) = c0 + M 2^t,     c0 = d - delta' = H(B) - W(O),     M = H - E + delta' = U_B(O),

with `B = L - v`, `O = N+(v)`, `W(O)` the sum over the pairs `{e, f}` of `O` of the two-path-cover coupling `W_ef`
(the cut matrix of THM-4114 and THM-5122), and `U_B(O)` the sum of the responses `U_B(y)`, `y in O`.

The law and both identities are known (Session216 star closure: `h - T + 2^t A` with `h = H(B)`, `A = U_B(O)`,
`T = W(O)`; `delta'` is its coordinate `T`). What is used here is the reading: `h - T` is the *offset* of the ray,
the path count of the base minus the coupling energy of the out-neighbourhood. A ray is thus a pair (base, vertex set): its offset is a quadratic
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
hostile fact "11 cannot be doubled" is a consequence.

**All 65 lines with `c <= 129` are complete.** Every odd `c <= 129` has a seed: 54 of order `<= 9` and, for
`c = 19, 33, 71, 79, 83, 97, 101, 103, 107, 113, 115`, eleven of order ten found through the equation
`W(O) = H(B) - 1`. The 119 points below the rays are path counts of order `<= 9`, except `5313 = 3 * 1771` and
`3425 = 5 * 685`, which are order joins. Hence **every odd `N > 1` whose `N - 1` has odd part at most 129 is a
path count, except 7 and 21** (THM-5126 (B3); seeds in `line_seeds.txt`, certificates in `below_rays.txt`).
Beyond: 2,016 of the 2,049 odd `c <= 4097` have a seed among the unit-offset roots of order `<= 10` (a complete
list: the audit found no unit set of size `>= 4` over a base of order `<= 9`) and the doubling vertices of order
`<= 9` (order nine: audit); doubling vertices of order ten were not searched. The 33 without are 245, 853, 909, 1049,
1407, 1467, 1515, 1577, 1629, 1829, 1887, 2023, 2037, 2085, 2343, 2393, 2703, 3001, 3033, 3123, 3313, 3411, 3427,
3545, 3713, 3733, 3753, 3773, 3907, 3909, 3985, 3991, 4023 (HYP-9658 (B)).

**Weight three.** `7 = 111_2` and `21 = 10101_2` both have binary weight three, and weight-three numbers are the
points of the lines `c = 2^d + 1`. Rows `d <= 7` are complete except for 7 and 21 (THM-5126 (B4)). In Horner form `N = 1 + 2^(a_1)(1 + 2^(a_2)(1 + ...))`: the gaps are the two
depth-two Horner numbers with equal small steps, `1 + 2(1 + 2)` and `1 + 4(1 + 4)`; the next one,
`73 = 1 + 8(1 + 8)`, is a path count.

## 3. The extension-mean law (THM-5126 (C))

The mean of `H` over all extensions of `T` by `k` new vertices is `k! 2^(1-k) sum_j C(k+1, j) 2^(-j) P_j(T)`, and
`P_j(T) >= j! S(n, j)` with equality only for transitive `T`. For one new vertex: mean `= H + P_2/4 >= 2^(n-2) + 1/2`.

The case `k = 1` is KNOWN (THM-4111 (4), (6); THM-4115 (5), with the variance in (6)); the identity for every `k`
and the bound are new as far as a canon search shows.

Reading: a one-vertex extension family has `2^n` members and mean at least `2^(n-2)`, with the least mean for the
transitive base (the one-apex family). This is consistent with the numerical observation that one-apex counts have
density about one half and miss 13, 29, 35, ...; the mean alone implies neither statement. For `k` apexes over a
chain the same identity
gives the word-growth exponent `k / log2(k + 1)`: `1` at `k = 1`, above `1` from `k = 2` (Session246's flag
eigenvalues `(k+1)^l` are the same growth). The value 13, the first number missed by the one-apex family, is the
first point of the line `L_3` above its gap.

## 4. Experiments (universes, controls, outcomes)

| experiment | universe | result |
|---|---|---|
| unit-offset roots (`seeds.c`) | every root of every class of order `<= 9` (1,723,824 at order 9) | 2, 5, 7, 16, 31, 68, 207, 1333 roots; 636 odd `c` with a seed |
| unit sets (`unitsets.c`) | every vertex set of size 2-4 over every class of order 9 (roots of order 10 with out-degree 2-4; the audit tested all sizes: none above 3) | 19,002 unit sets, counted as pairs (base class, `O`), not class-vertex pairs (2,771 pairs, 16,231 triples, no 4-set); 5,529 odd `c` with a seed in total |
| doubling vertices (`dbl.c`) | every vertex of every class of order `<= 8` | 48 of the 319 values above 1 have a doubling realiser (the value 1 trivially); for `N = 3 mod 4`, `N < 200`: 3, 19, 23, 31, 59, 75, 95, 99, 103, 111, 155, 163, 179, 191, 199 |
| two-module families (`modfam.c`) | every class of order `<= 8`, every ordered pair of vertices (406,412 triples) | 3,477 distinct exact bilinear tuples at order 8; none is additive (`de = 0` with `be, ga != 0`), in agreement with THM-5124 |
| weight-three map (`cover3.py`) | all tuples above, `a, b <= 26` | rows `d <= 4` and `d = 7` complete from order `<= 8` cores alone |
| pre-registered order-9 run | PREREGISTRATION.md, pushed before the run | P2, P3 CONFIRMED; P1 (count in [250, 420]; actual 636) and P4 (`c = 33`) FAILED; gate G1 holds |
| pre-registered order-10 run | second part of PREREGISTRATION.md, pushed before the run | Q1-Q4 CONFIRMED (all eleven missing `c <= 129` acquire seeds; 510 of 513 odd `c <= 1025`); gate G2 holds |

Positive controls: the scan recovers the fans (`c = 1`), the family `3 - X - Y + XY` of `C_3[TT_a, TT_b, 1]`, and
the known one-apex ray `1 + 9 * 2^k` (word `1 0 1^k 0 0`). Hostile controls: no ray passes through 7 or 21; the
false claim "every cycle vertex is doubling" is rejected (verify_s21.py). Doubling vertices are rare, and 15 has no
doubling realiser of any order (audit report, A3), so "double every realiser" is not a route.

## 5. What this does and does not do for the proof

- It reorganises the target. Chan-Pak is equivalent to: *every line `L_c` is in the spectrum, with one exception
  on `L_3` and one on `L_5`*. A ray settles a line from one point on; what remains on each line is finite.
- It does not reduce the problem to finitely many checks: there are infinitely many lines, each needs a seed
  (HYP-9658 (A)) and its finitely many points below the seed. The points below a seed are again ordinary targets.
- It extends a complete statement from the classical line `c = 1` (`N - 1` a power of two; fans) to the 65 lines
  `c <= 129`: all `N` with `N - 1 = c 2^s`, `c <= 129`. The set has density zero (about `65 log2 X` members up to
  `X`), so this says nothing about density. For calibration: every odd `N <= 2^35` other than 7 and 21 was already
  a path count by THM-4875 (D) (FINITE-EXACT, another construction), so the new content is the part of each line
  above that range, reached by one explicit ray per line.
- It turns "find a seed" into an equation in known objects: `W(O) = H(B) - 1` for the cut matrix `W` of the base.
  The order-ten run solved it for every `c <= 129` in 22 seconds; a construction of solutions for all `c` is open.
- It removes one architecture: lifting every realiser along its own line (the 15 control), and it adds a concrete
  gate for any engine: it must miss the root of `L_3` and the second point of `L_5` and nothing else on them.
- Rigorous lower bounds for the counting function of the spectrum are unchanged (THM-4887: `X^(2/3)` up to a
  logarithm). Positive density is still unproved.

Next obligations: (1) a grammar for the unit sets `W(O) = H(B) - 1`, to produce seeds for all `c` rather than find
them; (2) order eleven for HYP-9658 (B), and doubling vertices of order ten; (3) decide whether `s0(c)` is bounded
on average; (4) the points below the rays for `129 < c <= 4097` (for `129 < c <= 257` six of the 108 points are not
values of order `<= 9`: 5217, 5281, 5345, 10689, 4113 and 3153; audit minor 2 gives certificates for five of them).

**The seed equation in path-deletion form (added after the audit; proof below; checked by `evenpaths.py`).** With
`h(p) = H(B - V(p))` for a directed path `p` of `B`:

    W_ef = 2 sum over directed paths p between e and f, in either direction, with an even number of vertices, of h(p),
    U_B(y) = sum over directed paths p starting at y of h(p).

So a seed for `c` is a base `B` and a set `O` with

    sum over even paths p with both ends in O of h(p) = (H(B) - 1)/2,     sum over paths p starting in O of h(p) = c 2^j.

*Proof.* By Claim A (THM-4876 (1)), `H(B + v_O) - H(B) = 2 sum over odd cycles C through v of H(B + v_O - V(C))`;
such a cycle is `v`, then a directed path of `B` from `O` to its complement with an even number of vertices. The
cut form gives `W_ef = [H(B + v_e) - H(B)] + [H(B + v_f) - H(B)] - [H(B + v_ef) - H(B)]`, and the paths that survive
are those from `e` to `f` or from `f` to `e`. The formula for `U_B(y)` is its definition, `P2` being the path. QED
The first term is the arc itself, which is THM-5122's bound `W_ef >= 2 H(B - e - f)`. This is the inherited
odd-cycle expansion in explicit form; it is recorded because it turns the seed search into a statement about
paths with both ends in `O`.

## 6. Connections, typed (source / target / map / preserved / lost / test)

1. **Planar groups (Session238, THM-5096) - EXACT (first sentence, inherited) + ANALOGY (the rest).** Group
   closure of the physical port moves reaches 7 and 21; the positive semigroup does not. Here: a ray is a positive
   word `Z^t`, and its offset `c0` is its 2-adic limit. Unit offset means the ray converges 2-adically to 1, the
   path count of the transitive tournaments; no tournament with a tail of negative length is meant. Preserved: the
   affine law. Lost: positivity below `t = 1`.
2. **A000568 - EXACT identity (small).** The mean of the path-cover polynomial `pi_T(x) = sum_j P_j(T) C(x, j)`
   over all labelled tournaments on `n` vertices is `sum over permutations of odd order of 2^(psi - moved) x^cyc`
   (`psi` nontrivial cycles, `moved` non-fixed points); at `x = 1` it is Szele's `n!/2^(n-1)`. Davis's formula for
   A000568 sums `2^(orbits on pairs)` over the same permutations. Map: permutations of odd order. Preserved: the
   index set. Lost: everything else (one counts classes, the other paths; S239 says the same). PROVED for all `n`
   by the audit (report, A4): `E[P_j] = n! C(n-1, j-1) 2^(j-n)` by linearity, so the exponential generating
   function of the mean is `((2 + t)/(2 - t))^x = exp(2x artanh(t/2))`, the generating function of the right-hand
   side. `scripts/mean_pi_identity.py` checks orders 3-5 at `x = 1, 2, 3`.
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
  count; if it is at least `2^m / C`, the family has at most `C` members per unit of its mean, so no coverage
  statement can come from counting members; whether it covers is a collision (second-moment) question.
  (THM-5126 (C); `k = 1` is THM-4111.)
- **Translate the statement into every coordinate system before typing it new.** Trigger: an identity in rooted
  coordinates `(H, d, E, U)`. Action: rewrite it in the cut coordinates `(h, A, T)` and in ear notation and grep
  those forms. Evidence: MISTAKE-754 (this note), MISTAKE-750/751 (the previous one).
- Used from S90's list: "Additivity is orientation, never presence"; from the cards: "Search the statement before
  the method" (three collisions found before writing), "Existence is a maximum or tail question".

## 8. Random search log (owner's instruction)

Seeded draw of 24 canon titles (seed `2026111021`). Read: THM-4518 (three thieves: mean `O(1)`, second moment
`log`, so hitting probability `1/log`), THM-1126 (half-coverage gate), THM-157 (alpha_1 linear identity).
Transfer kept: THM-4518's shape - a family whose mean is bounded can still fail to hit because hits cluster; the
honest test of any "expected multiplicity" argument for coverage is a second moment, which for these families is
the open linear-energy question (THM-4887). Nothing else transferred.

## 9. Audit

`audit/AUDIT.md` (independent agent, fresh context, own class generator and counters): **PASS WITH CORRECTIONS.**

Majors, all applied:
- M1: the cut form `c0 = H(B) - W(O)`, `M = U_B(O)` is Session216's star-closure law, not new.
- M2: the mean law at `k = 1` is THM-4111 (4), (6) and THM-4115 (5).
- M3: the substitution law is THM-3121 (4); THM-2195 was the wrong citation.
- M4: "no seed of order `<= 10`" lacked its search scope. The audit tested every out-set size (no unit set of size
  `>= 4` over a base of order `<= 9`) and the doubling vertices of order nine (no new `c`; `s0` lower on nine lines).
- M5: three README sentences were false or overstated (the mean "forbidding" coverage; the mean as "exact form" of a
  density observation; "1 is a path count only as the empty tournament").

Seventeen minors applied (counts 7 not 8, 48 not 49, six not eight; citations; wording; the certificates of
`below_rays.txt` now carry witness tournaments from the audit's census and are recounted by `verify_lines.py`).

Confirmed independently: all 65 seeds and all 2,016 lines of `line_seeds.txt` (6,052 literal ray counts), the 121
points below the rays, the offset identities on 1,782,476 roots, the tail law, (C) up to `k = 6` on small orders,
(D), (E), both pre-registration texts against their commits, and every census number. Not verified by anyone:
doubling vertices of order ten. The audit also supplied four short proofs (A1-A4) used above. One remark was added
after the audit and is not covered by it: the path-deletion form of the seed equation (section 5).

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
    python evenpaths.py                  # path-deletion form of W_ef and U_B(y), 1,138 tournaments
    python add_witnesses.py              # adds the audit's witness tournaments to ../below_rays.txt
    gcc -O2 -o modfam modfam.c && ./modfam 8 3 < cls_8.txt > mf_8.txt           # about 4 minutes
