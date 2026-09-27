# Hedgehogs and functional graphs: unique ergodicity on the integers is the cycle half, creeping orbits cannot exist, the trivial cycle is exactly linear with an adelic saddle multiplier, and "Collatz is the only system" holds at the level of bounded-time certification in the qn+1 family and is open pointwise

**Session:** opus, `collatz-posets-zeta5-20260927` (fourth note), 2026-09-27.
**Owner's directive:** "what can we leverage thinking along the lines of
hedgehogs? ... think about proving something along the lines that Collatz
is the only system with a single fixed point at 1 where every natural number
is represented by a node in a graph with exactly one directed edge coming
out of it", with arXiv:2609.28797 and a pasted summary as inspiration.
**Source read:** the arXiv abstract of R. Pérez-Marco, *Unique ergodicity of
hedgehogs*, arXiv:2609.28797 (23 September 2026, 21 pages): "We prove that
the dynamics on a non-linearizable hedgehog is uniquely ergodic. This solves
a conjecture of the author formulated in 1995. Its unique invariant
probability measure is the Dirac mass at the indifferent fixed point. The
proof relies on the techniques of quasi-invariant curves and the hyperbolic
form of the Denjoy–Yoccoz lemma developed by the author." The pasted
summary's displayed "distortion estimate" and its bracketed citations are
not in the abstract and are not used here.
**Inherits (cited):** Lagarias 1985 (the 2-adic extension of `T` is
conjugate to the shift and preserves Haar measure; from memory), Terras 1976
(the parity-word bijection and density-one stopping times), the S12 note
(growth families are 2-adic shadows of negative cycles), the S13 note
(Proposition 4: no prefix rank), the parallel note
[`collatz_posets_dags_20260927_spine_descent_tree.md`](collatz_posets_dags_20260927_spine_descent_tree.md)
(Proposition 4: a first-descent word is an affine bijection from a 2-adic
source class onto a 3-adic landing class; section 6.1: the rank statement),
the Kuratowski note's strategy square and drift inequality, the synthesis's
DRIFT control (the `5n+1` undecided density near `0.2`), and this session's
first note.

**Status: PROVED elementary (Propositions 1–6) + FINITE-EXACT (script) +
CITED + DIRECTION. Collatz OPEN. The hedgehog theorem transplants to a
functional graph as an exact dictionary whose Collatz side is either the
cycle half of the conjecture, or false on `Z_2`, or the excursion-rank
obstruction; the family statement the owner asked for is proved at the
density-of-certification level and is open pointwise, exactly as the
conjecture is.** Script
`04-computation/experiments/collatz_hedgehog_family_20260927.py`, output
beside it. Independent audit OWED (as for the session's other notes).

## 0. The answer in one paragraph

A hedgehog is a compact invariant set on which an indifferent,
non-linearizable fixed point is the unique statistical sink: every orbit
creeps near it without being trapped, and the Dirac mass at it is the only
invariant probability measure. For a map of the positive integers with
exactly one edge out of every node, the dictionary is exact and short. The
invariant probability measures are the convex combinations of the uniform
measures on the cycles (Proposition 1), so "uniquely ergodic with the
trivial cycle's measure" is precisely "no other cycle": the
measure-theoretic face of the conjecture is its cycle half, and divergent
orbits carry no measure at all. Creeping is impossible (Proposition 2): an
orbit that visits a finite set infinitely often is periodic, so every
integer orbit is either trapped in a cycle or spends all of its statistical
time at infinity, and the conjecture reads "the empirical measures of every
orbit converge to the trivial cycle's measure". On `Z_2` the map is
conjugate to the shift and Haar measure is invariant and ergodic, so the
2-adic Collatz is the anti-hedgehog (Proposition 3): the cycle's 2-adic
neighbourhood of radius `2^(-k)` is visited with frequency `2^(-k)`. The
trivial cycle is the opposite of an indifferent point: `T^2(x) - 1 = (3/4)(x
- 1)` on the whole class `x = 1 mod 4` (Proposition 4), an exact Koenigs
linearization with multiplier `3/4`, which is a real contraction, a 2-adic
expansion by 4 and a 3-adic contraction by 3 (an adelic saddle, product
formula), and `T^(2j)(1 + 4^j t) = 1 + 3^j t`: a deep 2-adic approach to the
cycle is converted into 3-adic structure and costs size `4^j` to make and
`3^j` to leave, so integers cannot creep. Positive cycles attract their
shadows and negative cycles repel them (Proposition 5; S12's growth families
are the negative case). The hedgehog proof's key device, a hyperbolic
distortion estimate along close returns, would transplant to a strictly
decreasing potential over completed excursions, which is the rank that
S13's Proposition 4 and the parallel note's section 6.1 show cannot be a
function of local data. On the owner's uniqueness question: in the family
`T_q(n) = n/2, (qn+1)/2`, Terras's bijection holds for every odd `q`, and
the density of integers certified by coefficient descent within `k` steps
tends to 1 for `q <= 3` (Terras) and to `1 - f_inf(q) < 1` for `q >= 5`
(`f_inf = 0.177, 0.301, 0.386` for `q = 5, 7, 9`; exact DP to `k = 200`), so
Collatz is the unique member with `q > 1` whose bounded-time certificates
exhaust the integers in density (Proposition 6); the pointwise version
("all orbits reach 1 iff `q <= 3`") is open in both directions for `q = 3`
and `q = 5`, and the class must exclude contractions such as `n -> (n+1)/2`,
which also have a single fixed point at 1.

## 1. The dictionary for a functional graph (PROVED)

Let `T : N -> N` be any map (every node has exactly one outgoing edge).

**Proposition 1 (invariant measures are carried by cycles).** Every
`T`-invariant probability measure on `N` is a convex combination of the
uniform measures on the cycles of `T`. Hence the trivial cycle's measure is
the unique invariant probability measure iff `T` has no other cycle.

*Proof.* By Poincaré recurrence `mu`-almost every point is recurrent; a
recurrent point of a functional graph on a countable set lies on a cycle
(it returns to itself); so `mu` is supported on the union of the cycles, and
invariance forces it to be uniform on each. ∎

**Proposition 2 (no creeping).** If the orbit of `n` visits some finite set
infinitely often, it is eventually periodic. Consequently a non-periodic
orbit visits every finite set finitely often, its empirical measures
`(1/N) sum_(k<N) delta_(T^k n)` converge to `delta_inf` on the one-point
compactification, and the Collatz conjecture is equivalent to: for every
`n` the empirical measures converge to the uniform measure on `{1, 2}` (the
`T`-cycle).

*Proof.* Two visits to the same point close a cycle. ∎

So the hedgehog's defining phenomenon, creeping without trapping, has no
counterpart on the integers: statistically an integer orbit is a point mass
at a cycle or at infinity, never spread. The conjecture's two halves are the
two exclusions: no mass at a second cycle (Proposition 1), no mass at
infinity (Proposition 2).

**Proposition 3 (the 2-adic side is the anti-hedgehog).** On `Z_2`, `T`
extends to a homeomorphism-conjugate of the one-sided shift (Lagarias),
Haar measure is invariant and ergodic, and for Haar-almost every `x` the
frequency of visits to the ball `1 + 2^k Z_2` around the cycle is `2^(-k)`.
So the 2-adic dynamics has a continuum of ergodic measures and the trivial
cycle is a statistical sink for no typical point; unique ergodicity fails
as badly as it can. The integers are a null, dense subset on which the
conjecture asserts the opposite behaviour.

*Proof.* Cited (Lagarias 1985); the visit frequency is the ergodic theorem
for the indicator of the ball. ∎

## 2. The trivial cycle is exactly linear (PROVED)

**Proposition 4.** For `x = 1 mod 4`, `T^2(x) - 1 = (3/4)(x - 1)`; hence for
`x = 1 mod 4^j`, `T^(2j)(x) - 1 = (3/4)^j (x - 1)`, i.e. `T^(2j)(1 + 4^j t) =
1 + 3^j t`. The multiplier `3/4` of the return map `T^2` at its fixed point
`1` has `|3/4|_oo = 3/4`, `|3/4|_2 = 4`, `|3/4|_3 = 1/3`, product 1.

*Proof.* `x = 1 mod 4` gives the word (odd, even): `T^2(x) = (3x+1)/4 = 1 +
(3/4)(x-1)`; the word of `1` is `(odd, even)^inf` and Terras's bijection
makes `x = 1 mod 4^j` follow it for `2j` steps. Checked exactly for `x <=
4·10^4` and for `j <= 7`. ∎

Three readings. (i) The cycle is Koenigs-linearizable at the real place on
its entire 2-adic neighbourhood, with no error term: the opposite of a
hedgehog's indifferent, non-linearizable point. (ii) At the three places the
multiplier is a saddle: real attractor, 2-adic repeller, 3-adic attractor.
(iii) For an integer, 2-adic closeness to the cycle is expensive and is paid
back in 3-adic coin: `x = 1 + 4^j t` with `t >= 1` has `x >= 4^j + 1`, and it
leaves the shadow at `1 + 3^j t >= 3^j + 1`. This is the special case, for
the word of the trivial cycle, of the parallel note's Proposition 4
(2-adic source class `1 mod 4^j` onto 3-adic landing class `1 mod 3^j`,
affinely); what the hedgehog reading adds is the statement that the
approach to the sink and the size of the approaching point are tied by the
place-exchange `4^j -> 3^j`, so no integer orbit can creep toward `1`: to be
`2^(-2j)`-close it must be `4^j`-large, and the approach itself is a
descent by `(3/4)^j`.

**Proposition 5 (cycles attract or repel their shadows by sign).** Let `w`
be the valuation word of a cycle of the Syracuse map with `p` odd steps and
`A` halvings, and `x_w = S_w/(2^A - 3^p)` its rational cycle point. An
integer `m = x_w mod 2^K` follows `w` for `K` steps, and its value is
multiplied over each period by `3^p/2^A` up to a bounded carry. Since
`x_w > 0` iff `2^A > 3^p`, positive cycles are real attractors of their
shadows and negative cycles real repellers.

*Proof.* The cycle equation `x_w (2^A - 3^p) = S_w > 0`. ∎

S12's growth families (shadows of `-1, -5, -17`) are the repelling case;
the classes `1 mod 4^j` are the attracting case for the only positive cycle
of the plus sheet. The reset lane's "regenerated precision around the same
negative cycle" is, in this language, a repeated approach to a real
repeller: each approach is paid in growth, and the obstruction the lane
found is that nothing local prices the approaches.

## 3. What a hedgehog-shaped proof would need, and where it is already blocked

The hedgehog proof controls the accumulated hyperbolic distortion over a
cycle of close returns (the continued-fraction denominators of the rotation
number) and shows it cannot integrate to zero against a second invariant
measure. The Collatz transplant: close returns are the completed
excursions (the parallel note's spine blocks; the return times of the
height walk), the distortion is the accumulated log-multiplier `sum log(3/2^v)`
over an excursion, and the estimate would say that over every completed
excursion the accumulated multiplier is bounded above by a strictly negative
function of the distance from a baseline. That is a strictly decreasing
potential over completed excursions, i.e. the excursion rank; S13's
Proposition 4 shows it cannot be a function of the parity prefix, and the
parallel note's Theorem 3 sharpens the witness (the same block sequence is
realized by every integer in a residue class, with every possible future).
Two further mismatches: the hedgehog's ambient metric is hyperbolic on the
complement of the invariant set, where the map is a Schwarz–Pick
contraction, while `T` on `Z_2` is uniformly expanding (multiplier 2 per
step, conjugate to the shift) and contracts only on average at the real
place; and the indifferent case corresponds to a mean multiplier of modulus
one, which in the `qn+1` family is the non-integer `q = 4`, not `q = 3`. The
Collatz cycle is attracting at the real place (Proposition 4); the
conjecture's difficulty is non-uniform hyperbolicity (unbounded excursions
of a negative-drift walk), not neutrality.

## 4. "Collatz is the only system" (PROVED at the certification level; OPEN pointwise)

Every self-map of `N` is a graph with one edge out of each node, and
contractions such as `n -> (n+1)/2` (rounded up) have a single fixed point
at 1 attracting everything, so the uniqueness must be asked inside a class.
The natural class is the family `T_q(n) = n/2` (`n` even), `(qn+1)/2` (`n`
odd), `q` odd.

**Proposition 6.** (a) For every odd `q`, `n mod 2^k -> (parity word of the
first k steps)` is a bijection (Terras's argument; checked for `q <= 11`,
`k <= 12`). (b) Let `f_k(q)` be the fraction of parity words of length `k`
with no coefficient descent (`q^(o_j) > 2^j` for all `j <= k`); it is the
natural density of the integers not certified by coefficient descent within
`k` steps. Then `f_k(1) = 0` for `k >= 1`; `f_k(3) -> 0` (Terras; here
`f_200 = 3.1·10^(-6)`, decaying like `2^((h* - 1)k)`); and for `q >= 5`,
`f_k(q)` decreases to a positive limit `f_inf(q)` (`0.177, 0.301, 0.386` for
`q = 5, 7, 9`). (c) Consequently the density of integers certified within
`k` steps tends to 1 iff `q <= 3`; for `q >= 5` the set of integers whose
orbit never has a coefficient descent has upper density at most `f_inf(q)`,
and the set certified in bounded time has density at most `1 - f_inf(q) < 1`
for every bound.

*Proof.* (a) `T_q` is affine with odd multiplier on each branch, so the
parity of `T_q^j(n)` depends only on `n mod 2^(j+1)`. (b) The log-multiplier
walk has steps `log_2 q - 1` and `-1` with probability `1/2` each, mean
`(log_2 q - 2)/2`, negative iff `q <= 3`: a negative-drift walk returns below
its start almost surely, so `f_k -> 0`; a positive-drift walk stays above
its start with positive probability, so `f_k` decreases to a positive limit
(the DP computes it exactly to `k = 200`). (c) is (a)+(b) with the densities
of the periodic sets `{n : no descent within k steps}`. ∎

So the exact form of the owner's statement is: **among the maps `T_q` with
`q > 1`, Collatz is the unique one whose bounded-time coefficient-descent
certificates exhaust the integers in density.** What remains open in both
directions is the pointwise statement: that every orbit of `T_3` reaches 1
(the conjecture), and that some orbit of `T_5` does not (believed; the
orbit of 7 exceeds `2^3110` within 20,000 steps and returns below 1000
only ten times, but no divergence is proved for any `q`). The known cycles
with minimum below 20,000 are `{1}` for `q = 1, 3, 7` and `{1, 13, 17}` for
`q = 5`; the synthesis's DRIFT control records the `5n+1` undecided density
near `0.2`, which is `f_k(5)` at the depths it used.

A rigidity theorem of the form "if a map in a natural class has all orbits
converging to a single cycle then it is Collatz" would need, for the
`T_q` family, both the conjecture for `q = 3` and divergence for every `q
>= 5`; neither is available, and the family's own dichotomy (Proposition 6)
is the strongest exact statement in that direction.

## 5. Directions (DIRECTION; none pursued)

* **D12. Price the approaches.** Proposition 5 says every deep 2-adic
  approach of an integer orbit to a negative cycle is paid in growth and to
  the positive cycle in descent. A potential `V(m) = log m + c(m)`, with
  `c(m)` a function of the 2-adic distances of `m` to the finitely many
  negative cycle points, would decrease along positive-cycle shadows and
  could absorb the growth along negative-cycle shadows only if `c` foresees
  the depth of the coming regeneration; the reset lane's construction
  regenerates arbitrary depth, so `c` must read the integer. This is the
  hedgehog "quasi-invariant curve" (a baseline against which distortion is
  measured) in Collatz clothing, and the same obstruction.
* **D13. The measure-theoretic cycle half.** Proposition 1 makes "no other
  cycle" a unique-ergodicity statement; the classical route (Steiner,
  Simons–de Weger, Hercher) proves it for cycles with few circuits by
  Diophantine approximation of `log_2 3`, which is the analogue of the
  continued-fraction bookkeeping of the rotation number in the hedgehog
  proof. Whether a unique-ergodicity argument could replace the
  cycle-by-cycle exclusion is not known to this note; on `N` the two
  statements are the same, so nothing is gained by the reformulation
  alone.

## 6. Reproduction

    cd 04-computation/experiments
    python3 collatz_hedgehog_family_20260927.py > collatz_hedgehog_family_20260927.out

Standard library only; about two minutes.
