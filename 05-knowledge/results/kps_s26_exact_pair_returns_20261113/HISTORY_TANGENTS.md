# Old threads that rhyme with the exact-pair refill

kind-pasteur-2026-11-13-S26. The owner asked for the session to follow tangents into old work, including the
problems introduced around each old session. Three read-only sweeps covered: Heegner numbers, Ramanujan's constant
and "30000"; Lonely Runner; and k-uniform hypergraphs with Collatz. Each sweep also read the sessions just before
and after its hits.

Everything below is a pointer or an analogy unless marked *checked*. A *checked* item was recomputed in this
session; its script is in `scripts/` or the line says how. No claim here enters THM-5151.

## A. Heegner numbers, Ramanujan's constant, 30000

**The 2647 coincidence.** *Checked.* Session 310's prime cyclic-frame value Q_B(2) + H(B) = 2647 is the largest
odd number that no uniquely Hamiltonian (C = 1) tournament attains. That list of 29 misses is Session 197's;
HYP-9494 (D), this machine's S14, has it to 10^9, and THM-4884 names it.

The order-6 exact pair gives the frame value 2 Q_(Y5)(2) + H(Y5) = 231, also a C = 1 miss. Among the 85 C = 1
strong tournaments of order <= 9 whose value Q_X(2) + H(X) is <= 2647 (45 distinct values), exactly two values are
misses: 231 and 2647. They come from exactly the strong members of the two detachable exact pairs.

No mechanism is known. A random 45-subset of odd numbers below 2647 hits the 29 misses about once on average.
Treat this as a coincidence until a reason appears (`scripts/match.py` and the snippet in README).

**Heegner numbers among restricted-construction misses** (sweep report; partly checked).

- 7, 43, 67, 163 are C = 1 misses; 11 and 19 are not.
- The one-arc-join family F1 has zero representations at 7, 43, 67, 163.
- Ramanujan's ternary form x^2 + y^2 + 10 z^2 misses 7, 21, 43, 67 among small values.
- HYP-9657 (B), opus-S90, gives a partial reason. Values N where N(N+1) has no small odd prime factor carry less
  mass, and are missed about 3.3 times as often.

**THM-5125 Heegner tournament.** *Checked.* On {3, 7, 11, 19, 43, 67, 163}, put p -> q iff p is a square mod q.
Quadratic reciprocity for primes ≡ 3 (mod 4) makes this a tournament, and H = 71.

**THM-4666 / THM-4751** (Oct 17 and Oct 21).

- h(-p) measures how badly the natural order ranks the Paley tournament P_p. The extreme case is exactly
  p = 7, 11, 19, 43, 67, 163.
- P_q's spectral ring is a principal ideal domain exactly at the Heegner primes.

**30000.** Several different "30000"s recur:

- 30,000 random 7-vertex tournaments with no H = 21 (S71f, March).
- A 65-chart bank attaining 14,998 of the 15,000 odd values below 30,000, missing only 7 and 21 (golden seam,
  Nov 3).
- This machine's S18 bases 4-8, exact below alpha-sum 30,000 (THM-5070 (6)).
- Berggren SYNTHESIS §5: "what 107, 73, 17 and 30000 preserve".
- *Checked:* 17, 73, 107 are the 7th, 21st and 28th primes. The indices are the two forbidden values and their sum.
- *Checked:* 740^2 + 1 and 760^2 + 1 are the only pair of primes of the form n^2 + 1 that differ by 30000.
- -j((1 + sqrt(-11))/2) = 32^3 = 32768 is also near 30000.

**Ramanujan's constant.** e^(π √163) ≈ 640320^3 + 744 (opus-S84 in March: 640320 = 2^6·3·5·23·29). The Oct 6
packet `ramanujan_heegner7_trivial_cycle_20261006.md` places the trivial Collatz cycle code {1,2,4}/7 at the
Heegner point of discriminant -7, where j = -15^3.

*Checked:* 3, 7, 19, 163 = 2·3^k + 1 for k = 0, 1, 2, 4. The k = 3 term is 55, itself a C = 1 miss. The repository
types this as numerology.

**Rabinowitsch-type clock criterion.** This is THM-5151's Proposition 7. For a clock H_k = I + c 2^(r-k) 3^k, a
prime l either never divides a term or divides exactly one residue class of k modulo ord_l(3/2). Because 3/2 is a
primitive root mod 7, a run of length 6 always meets a 7-multiple unless 7 | I.

For Session 309's clock (I = 499), the primes 13, 19, 29, 43, 47, 53 never divide a term (*checked*). The
inert-prime analogy for Q(sqrt(-163)), whose small primes are inert, is the sweep's; this note does not adopt it.

## B. Lonely Runner (sweep report)

**The fixed-bank no-go appeared four months earlier.** THM-2072 (Jul 21) shows the core {1..10, lcm(F ∪ {2..14})}
blinds every fixed clock bank F. THM-4079 (Aug 25) escapes with a clock that depends on the core. THM-4349 (Sep 2)
proves the adaptive finite clock complete: safety for every pair of odd tails holds iff one adaptive clock
certifies it, and a failure is a real bounded counterexample.

The Chan–Pak translation is Session 309's FIXED_BANKS. Its escape would be an adaptive Π modulus at which the count
depends only on the refill module's statistics, together with an absorption lemma for large refills. Not
attempted here.

**Renewal with exact tolls** (THM-4335 §7-9).

- The composed clock is H_m = (3^m H_0 - sum_i 3^(m-1-i) 2^i a_i) / 2^m.
- The sweep checked it on Session 309's run: a = 499, m = 6, from 83,405,299 to 950,033,299.
- LRC's version adds a greedy farthest-reach choice, proved minimal by exchange, and "collision-backed" handoff:
  install the refill before the last site is spent.

**THM-4025's odd-dilation semigroup** (Aug 24). φ(n) = 2n - 1 with n★m = 2nm - n - m + 1 is literally odd H under
ordinal sums. Its "survivors form a ★-upward ideal" is the prime/7p reduction in another language.

**The Jun 22 bridge** (HYP-2878/2879).

- 7 is the only permanent prime gap.
- The first strong atoms with H divisible by 7 appear at n = 7.
- Balanced 7 -> 8 ears miss exactly H = 49 and 75, and 49 is the first 7p target.

**Cycle invoice** (`clock_owned_runner_events_20261021` §4). A prescribed cyclic handoff is realizable at all
large heights iff the product of its ratio windows is below k^l. The packet names Session 309's gap: without a
re-anchoring theorem, graph cycles do not compose in real time.

**A caution from LRC's own treadmill.** The neutral-refill chain THM-4311 -> 4313 -> 4318 kept the carrier at 3,925
masks while exchanges of size 1, 43 and 9 were needed. It was never made uniform and was abandoned for direct
closures.

**Neighbours.**

- Jun 19-22: the resonance-killing game; the lcm family showing finite residue bases cannot close LRC.
- Jul 20-22: GMC(2) proved in Lean. THM-1820 calls LRC and GMC(2) "one moment nullcone".
- Sep 29 - Oct 1: the owner's "LRC local / Collatz global" thesis; THM-4532, "Tournament Clock Prime Collatz".
- Oct 20-22: Kepler, Hadwiger–Nelson and Leech-lattice reformulations.

## C. k-uniform hypergraphs and Collatz (sweep report)

**THM-4803** (Oct 26). The 3-uniform hypergraph of cyclic triples C3(T) determines the full ordered path-cover
profile, and its strong components are its connected components (Lemma C).

In that language an exact strong/decomposable pair is a connected and a disconnected cyclic-triple hypergraph with
the same profile. The fibre moves of THM-4803 cannot create a singleton component. So a site-creating rebase must
change the hypergraph, which is exactly what Session 310's trade and Proposition 2 of THM-5151 do.

**Packing grammars** (Session 124, Oct 31; Session 261, Nov 11).

- With ν <= 1 (at most one disjoint triangle), the values up to 73 form an exact list; 29, 35, 41, 43, 55, 57, 59,
  67, 71 need two disjoint triangles.
- The maximum packing over a repeatable alphabet equals the matching number of a "resource hypergraph"; deletion
  equals its cover number.
- The sweep notes that Session 309's HYPERGRAPH_SEAMS re-derives part of this without citing it.
- The conjecture there, that an unbounded prime/7p return needs ν -> ∞, is the packing-number counterpart of
  HYP-9678 (C) below.

**THM-079 Part R** (Mar 7, the H = 21 proof). A "poisoning graph" on a maximal triangle packing has out-degree at
most 1 and no cycles, so a safe vertex deletion always exists except at C3 => C3 => C3. This is a model for a safe
peel of a spent site, and possibly the reverse insertion.

**Session 163 LATER_RETURNS** (Nov 3). A strong rewiring of C3 ⊕ X overshoots by at least 2H(X - v) or
6H(X - {a,b}). This is a precursor of Session 310's contact weights.

**Collatz.**

- The anchor isometry v_p(F^j(x) - x) = t + v_p(j) (collatz_guards, Sep 21 -> Session 101, Oct 29) is the 7-adic
  isometry of Session 309, found earlier.
- The Terras clock (THM-4569/4581, Oct 7): coupled orbits differ by a simple random walk.
- THM-4594: 2-adic certificates cannot certify the class of -1.
- THM-4849 (Oct 29): circulant tournaments are rational Terras cycles with denominator 2^q - 3^k, plus a locking
  theorem. *Checked:* the circulant on Z_7 with connection set {1,2,3} has H = 175 = 2^8 - 3^4. THM-4849 itself
  warns against promoting this.
- THM-4771 (Oct 25): the Antihydra / A11 "BB(6) cryptid" machines. In x = H - a the clock H -> (3H - a)/2 is
  x -> 3x/2, so the site supply behaves like the Antihydra counter: spend one per step, gain on refill.

**Neighbours.**

- Mar 6-8: H = 21 proved for all n; the exhaustive n = 8 run; path homology.
- Oct 6-8: THM-4591, integer Terras cycles of periods {1, 2, 3, 11} read as musical intervals.
- Oct 28-29: the owner's 19, 79, 169, 101 dispatch.
- Oct 30 - Nov 1: PRIMEGAME / FRACTRAN; THM-4856's FERMATGAME halts iff a sixth Fermat prime exists.

## D. What the tangents suggest for the return rule

1. **Adaptive modulus.** In LRC, fixed banks fail (THM-2072) but core-dependent clocks complete (THM-4079,
   THM-4349). FIXED_BANKS plays the role of THM-2072 here. A core-dependent Π modulus, chosen after the right
   module, is the untried analogue.
2. **Spend before refill is wrong; refill before spending the last site.** This is LRC's collision-backed handoff.
   In THM-5151's language: before the last transfer, the left or exterior part must already hold an even component
   with a partner, because spent odd cores have no hidden capacity (Propositions 2-3).
3. **Hypergraph view.** Site creation must change C3(T) (THM-4803 Lemma C), and every such change is an R-factor
   recombination. The resource-hypergraph numbers ν, τ of Session 261 are cheap necessary filters for candidate
   partners.
