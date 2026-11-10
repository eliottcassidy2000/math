# S90 (opus, 2026-11-10): additive rigidity, the residue-tournament model, and where the mass is small

**Status.** PROVED: THM-5124 (additive rigidity of the Hamiltonian-path count on presence cubes), THM-5125
(quadratic-residue tournament model, Heegner domination, cycle-gas identities). NUMERICAL with three pre-registered
out-of-sample tests, all CONFIRMED: HYP-9657 (A)-(B). FINITE-EXACT refutations: HYP-9657 (C).
**The Chan-Pak conjecture - every odd number other than 7 and 21 is the number of Hamiltonian paths of a
tournament - remains OPEN.** That 7 and 21 are not values is THM-1370; only completeness is open. No numerical
coverage range is extended here.

Owner directive (2026-11-10): aim to prove the statement; no more numerical evidence, but structural reasoning and
experiments, connections to other areas and to past work in this repository, and names for the principles; consider
the Heegner numbers, Ramanujan's constant, the Birch and Swinnerton-Dyer conjecture, and "elliptic curves are
defined by a cubic and a linear power; odd exponents in functions relate to odd cycles in graphs", to be corrected
and refined.

## 0. Inheritance, portfolio, concept board

- *Closest proved mechanisms.* The odd-cycle formula `H = I(Omega, 2)` and its Moebius positivity on the presence
  cube (THM-4114); the universal increment floor `H(Z) = h` or `H(Z) >= (3h + 1)/2` (Session214); the exchange-cube
  selector of THM-5122 (Session246, finished the same morning on the same directive); the one-arc join law and the
  F1 mass (THM-4875 (D), THM-4885); the Monsky matrix of the reciprocity tournament (bsd_selmer_growth_20261013a).
- *Canonical hostile.* The value 13: the whole one-apex class omits it (Session 16), so no single-hub grammar is
  complete; and `7 = 1 + 2 + 2 + 2`, the additive prediction that never happens.
- *Corrected near misses.* MISTAKE-229 / 559 / 568 / 576 (Heegner or Paley language without a map); Session242's
  "elliptic halving proves descent: false"; this session's own two refuted bridges (section 4).
- *Least-used sidecar.* The presence cube itself (which induced subtournaments occur together), as opposed to the
  ear cube; and the representation count `r(N)` as a function of the arithmetic of `N` rather than of its size.
- *Portfolio.* Anchor: a proof-relevant structural theorem for Chan-Pak. Niche: the arithmetic (Redei) model the
  owner's hints point to. Wildcard: why famous numbers (7, 43, 67, 163) sit in exceptional lists of sub-families.
- *Concept board.* presence cube / exchange cube; interaction (Moebius atom); hub pair; cycle gas and its
  characters; mass versus congruence; roughness of `N(N + 1)`; fat versus thin orbit.

## 1. Additive rigidity (THM-5124)

For disjoint vertex sets `W1, W2` and `B` the rest, the interaction
`H(B+W1+W2) - H(B+W1) - H(B+W2) + H(B)` is the weight of the odd-cycle packings meeting both sets, so it is `>= 0`
(THM-4114). The new statement is the equality case:

> if both sets raise the count of `B` and the interaction vanishes, then the four counts are `c, 3c, 3c, 5c`, exactly
> one vertex in each set matters, and those two vertices are the two one-triangle vertices of a strong
> 4-tournament whose hubs lie in `B`. Three sets can never be simultaneously active and additive.

So `5 = 1 + 2 + 2` is, up to scale, the only additive identity among counts of nested induced subtournaments. A
third parallel switch between the same two hubs creates a 5-cycle through all three switches and the count is
`9` (or `15`), not `7`; `k` parallel switches give `2^k + 1`, not `2k + 1` (statement 4: the two-hub suspension has
`H(A) + D(A)` paths). Two proofs of the dimension bound are given: from the pair lemma, and in three lines from
Session214's increment floor.

*Why it matters for the target.* THM-5122 shows that in the exchange cube of one core (reversing which side of a
cut a vertex is on) additive cubes exist, because interactions are signed. THM-5124 shows the opposite for
presence. Together: **additivity of `H` is an orientation phenomenon, never a presence phenomenon.** Every proof
architecture that needs three independent additive contributions from the presence of vertices - Brown's
criterion, Erdos-Lewin sums of 3-smooth numbers over antichains, a sumset of two independent rich families,
affine presence cubes - is ruled out, not merely unfound (one_apex_join_representation_mass_pi2_over_12_20261106
section 2 had recorded that no free-staircase family was found). Remark 4 of the theorem adds a width bound for
THM-5122's orthogonal pair-exchange cubes with strong closures: an interval obtained from one such cube at baseline
`B` has length at most about `B^0.65`.

## 2. The residue-tournament model (THM-5125)

On the primes `= 3 (mod 4)`, "`p -> q` iff `p` is a quadratic residue mod `q`" is a tournament by reciprocity; by
Dirichlet and CRT it has the extension property, so it is the countable generic tournament and every finite
tournament is `T_n` for some squarefree `n`. Then `H(T_n)` is the number of *residue chains* of `n` (orderings of
the prime factors with each a residue modulo the next), Redei's theorem says this number is odd, and

> Chan-Pak  <=>  the number of residue chains of a squarefree `n` with prime factors `= 3 (mod 4)` takes every odd
> value except 7 and 21.

Least `n` with `N` chains: `a(1) = 3`, `a(3) = 231 = 3*7*11`, `a(5) = 4389 = 3*7*11*19`,
`a(9) = 188727 = 3*7*11*19*43`, `a(11) = 100947`, `a(13) = 206283`, `a(15) = 311619`, ... (`redei_model.out`).
This is a reformulation; it gives no leverage on completeness, because residue tournaments are all tournaments.

**Heegner domination (statement 2).** If `Q(sqrt(-q))` has class number one, every prime below `(q + 1)/4` is inert,
so `q` is a quadratic residue modulo every prime `p = 3 (mod 4)` below `(q + 1)/4`: in the residue tournament the
Heegner prime `q` beats all those `p`. The tournament on `3, 7, 11, 19, 43, 67, 163` therefore has 10 forced arcs
out of 21; it has `H = 71`, prefix counts `1, 1, 3, 5, 9, 33, 71`, and the six primes `= 3 (mod 8)` have `H = 29`
and Monsky nullity `4`. These numbers are data, typed FINITE-EXACT; no meaning is claimed for `71`.

## 3. The owner's themes, corrected

Each line is typed as in the protocol: source / target / map / preserved / lost / test.

1. **"Odd exponents and odd cycles."** *Correct form:* reversal of orientation acts on a directed `k`-cycle by the
   sign `(-1)^(k-1)` of a `k`-cycle permutation. In the Redei-Berge symmetric function the cycles of `T` and of its
   converse are summed with that sign, so even cycles cancel and odd cycles double: `U_T` is a polynomial in the
   odd power sums `p_1, p_3, p_5, ...` with nonnegative coefficients (Grinberg-Stanley), and `H = U_T(1, 0, 0, ...)`
   `= sum over odd packings 2^|P|`. Map: cycle of length `k` -> `p_k`. Preserved: disjointness (products). Lost:
   which cycles intersect. This is an exact object map, and it is the reason the conflict graph is built from odd
   cycles only. Earlier notes (collatz_mod6_20260917_*, seam_power_coordinates_20260925) recorded "no map" for the
   word *odd* alone; the map is the sign character, not the parity of an exponent of a function.
2. **"Cubic and linear."** The curves the repository connects to tournaments are `y^2 = x^3 - n^2 x`, an odd cubic
   in `x` (extra symmetry `x -> -x`, CM by `Z[i]`). For `n` a product of primes `= 3 (mod 8)` Monsky's formula
   (CITED) reads the pure 2-Selmer dimension off the reciprocity tournament as `nullity(M^2 + M + I)`. THM-5125 (4)
   makes the comparison exact: per(I + A), det(I + A), `H` and this nullity condition are evaluations of one
   cycle-packing enumerator. `H` uses the *parity* character with weight `2` over `Z`; the Selmer condition uses
   the *cubic* character (`x^2 + x + 1`, the element of order `3` acting on the 2-torsion) with coefficients in
   `F_2`: it holds iff the numbers of packings in the three classes of `e - o` modulo `3` have equal parity. So
   "cubic" belongs to the Selmer side and "odd" to the Hamiltonian side: two characters of one gas, not one
   statement. Test for more: HYP-9657 (C) - no congruence for `H` follows from the `F_2` data at small orders.
3. **Heegner numbers.** Two exact facts and one mechanism. (i) Domination, above. (ii) Class number one forces `q`
   and `(q + 1)/4` to be prime (Rabinowitsch, CITED). (iii) The F1 join family has representation mass of
   additive-divisor type in `N` and `N + 1` (HYP-9657 (A): classes `N = 0, -1 mod p` are favoured for every odd
   prime tested, in and out of sample; no dependence on the Kronecker symbol `(-N/p)`), so integers with
   `N(N + 1)` free of small odd primes are the lightest (mean ratio 0.72, miss-rate ratio 2.6-3.3), and the Heegner
   numbers are such integers: `r(7) = r(43) = r(67) = r(163) = 0`. This is why 7, 43, 67, 163 appear among the
   misses of F1 and of the O/Z/J readout (THM-4884's 29 values). It is a statement about sub-families. It does not
   explain 21 (`= 3 * 7`, not rough) and it is not a class-number formula: where the mass of a representation
   problem *is* a class number (sums of three squares: `r_3(d) = 24 h(-d)` for `d = 3 mod 8`, `d > 3`, Gauss, CITED;
   so `r_3 = 24` at `11, 19, 43, 67, 163`, the least possible), the Heegner numbers are extremal for a different
   reason. Type: MECHANISM (statistical), scope F1.
4. **Ramanujan.** (a) Ramanujan's constant `e^(pi sqrt 163)` is the class-number-one fact `j(tau_163)` in `Z` again
   (see ramanujan_heegner7_trivial_cycle_20261006.md and MISTAKE-568 for attribution); nothing new is claimed.
   (b) Ramanujan's ternary form `x^2 + y^2 + 10 z^2` misses the odd numbers `3, 7, 21, 31, 33, 43, 67, 79, 87, 133,
   217, 219, 223, 253, 307, 391, 679, 2719` (the last two and completeness under GRH: Ono-Soundararajan 1997;
   finiteness: Duke-Schulze-Pillot 1990). The O/Z/J readout misses 29 odd numbers, largest 2647. The lists share
   `7, 21, 43, 67`. Under a uniform null on the odd numbers below 100 (nine members each) an overlap of at least
   four has probability about 0.04; both lists favour integers of small local mass, so the null is too weak. Type:
   ANALOGY of kind (a non-regular family with a finite sporadic exceptional set), COINCIDENCE of members. The
   mechanisms differ: there, mass = class number + cusp-form coefficient; here, no Kronecker dependence is seen.
5. **Birch and Swinnerton-Dyer.** The only map is Monsky's (item 2). As a dictionary, typed ANALOGY throughout:
   the conflict polynomial `I(Omega(T), x)` is multiplicative over strong components (an Euler product,
   THM-1862/1926), `H` is its value at `x = 2`, its value at `-1` is an Euler characteristic and equals `H` mod `3`
   (THM-5125 (4a)); the spectrum problem is local-global with the single local condition at `2` (Redei), and
   `{7, 21}` plays the role of a finite obstruction set visible only globally, removed by a descent (strong floor
   plus deletion, THM-1370); the "additive rank" of `H` is at most `2` (THM-5124). None of these is a theorem about
   elliptic curves.

## 4. Bridges tested and refuted (first failed implication recorded)

- `H mod 4` (or `mod 8`) as a function of the Redei/Monsky `F_2` data: true at orders 3-4, false at order 5.
  First failed implication: equal (Laplacian corank, Monsky nullity, nullity of `I + A`) does not give equal number
  of odd cycles modulo 2. Survivor: Redei's parity and the identities of THM-5125 (4).
- "Monsky nullity `> 0` implies `3 | H`": true at orders 3-4, false at order 5 (`H = 13`, the five-vertex
  tournament with out-neighbour masks `[18, 8, 3, 5, 14]`, which is Session245's h13 cell). With the pure adjacency
  matrix it survives to order 5 and fails at order 6 (`H = 17`). Survivor: for a single triangle both hold.
- A class-number-type (Kronecker symbol) mass for F1: not seen; the bias is `N = 0, -1 mod p`.
- Sum of two independent rich families (Mann-type sumset of one-apex sets, which have density about one half):
  impossible by THM-5124 (2).

## 5. Where this leaves the proof

1. *Dead:* presence-additive architectures (THM-5124); single-realiser window inductions (THM-4876 (5), Session214);
   one-hub grammars (Session 16, value 13); finite banks as a proof pattern (their expected miss density is positive,
   HYP-9480 (C)).
2. *Alive, with a quantified limit:* exchange cubes (THM-5122) - width at most about `B^0.65` per cube with strong
   closure, so a covering needs polynomially many baselines per octave; moving cores (THM-4960, THM-5000).
3. *The analytic route.* The three-letter rooted system `S, B, F` on states `(h, R, e, u)` (it is the O/Z/J grammar;
   rules re-derived and checked against brute force in `psi_sbf.py`) satisfies `h <= R` at every state and `R` at
   most doubles with each letter, so the `2 * 3^(k-1)` words of length `k` not ending in `S` all have count at most
   `2^k`: on average at least `(3/2)^k` words per value. In this sense it is a fat orbit, not a thin one
   (word-growth exponent `log 3 / log 2 = 1.58 > 1`); its readout nevertheless misses 29 odd values (THM-4884). For
   such orbits the available technology (bilinear circle
   method with a spectral gap, as in Bourgain-Kontorovich for Zaremba) yields density one with a power-saving
   exceptional set, not all values; this typed ANALOGY is already in the repository (seam_gadget_thickness_spectrum,
   Session205) and is repeated here only to state the realistic milestone: positive density of the spectrum is
   still unproved (best: `A(X) >= X^(2/3) / (13 ln X + 102)`, THM-4887).
4. *What would finish.* Session246's two routes stand (a family of cores with exchange systems whose absolute
   intervals cover every target; or a well-founded physical predecessor rule). THM-5124 adds a constraint to the
   first (no presence switches) and Remark 4 a size requirement.

## 6. Principles, named

These are names for mechanisms that recur in this thread; each has a trigger and a counterindication. They are
recorded here, not promoted to META-PATTERNS (the evidence is from one problem family).

- **Two characters of one cycle gas.** Trigger: a count over disjoint cycles appears next to a determinant or a
  rank. Action: write both as evaluations of the cycle-packing enumerator and compare characters (parity over `Z`
  for `H`; cubic over `F_2` for Monsky). Counterindication: equal small values are not a congruence (section 4).
- **Additivity is orientation, never presence.** Trigger: a selector design that adds vertices to add to `H`.
  Action: check the interaction; by THM-5124 it vanishes only in `1, 3, 3, 5`. Use exchanges (THM-5122) or carries.
- **Mass, not congruence, decides the exceptions.** Trigger: an exceptional list contains famous numbers.
  Action: compute the class means of the representation count and its dependence on the roughness of the natural
  linear forms (`N`, `N + 1` here) before reading the list. Counterindication: families of mass below one.
  Related cards: "Separate local support from bounded-height coverage", "Existence is a maximum or tail question".
- **A height with a gap.** `log H` is superadditive over vertex partitions and multiplicative on order joins
  (Session214: `h >= p + q`), and positive increments are at least a factor `3/2`. Trigger: an induction that
  needs small increments. Counterindication: the gap says nothing about unrelated witnesses.

Cards used: "Type every analogy and every implication"; "Search the statement before the method" (the Session214
floor and THM-4114 positivity were found before the rigidity proof was written, and THM-5122 before the scope was
fixed); "Price deletion by the atoms it removes"; "Existence is a maximum or tail question".

## 7. Restatements used (no novelty claimed)

`psi(0) = 0`, `psi(2h) = psi(h) + h`, `psi(2h+1) = psi(h) + floor(h/2) + 1`; the one-apex tournament with reduced
word of value `2h` has `2 psi(h) + 1` paths; `psi(h) = h - q(h)` with `q` the weighted count of pairs of one-bits
(Session 16's weighted inversions, THM-4876 (6), THM-4885/4887 in other coordinates). Equivalently
`G(n) = sum over zero bits p of floor(n / 2^(p+1))` is the number of earlier terms of the base-2 van der Corput
sequence exceeding its `n`-th term (the exact form of the Hammersley remark in the 2026-11-06 mass note). The
splitting law `psi(m 2^j + t) = psi(m) + psi(t) + kappa_j(t) m - t_(j-1) ceil(m/2)` restates THM-4887's block law.

## 8. Reproduction

    cd 04-computation/experiments/s90_rigidity_redei_20261110
    python3 rigidity_check.py    > rigidity_check.out        # about 40 s, needs numba
    python3 redei_model.py       > redei_model.out           # about 10 min
    python3 cycle_gas.py         > cycle_gas.out
    python3 psi_sbf.py           > psi_sbf.out
    python3 f1_mass.py 21 18        > f1_mass_insample_21.out
    python3 f1_mass.py 23 21 prereg > f1_mass_prereg_23.out   # pre-registered run

Every script exits nonzero on a failed check. Universes are stated in each docstring.

## 9. Next obligations

1. Quantify the interaction: for two active vertices not in the rigid configuration, bound the interaction below
   in terms of the two increments (a second-order gap). Cheapest test: census at orders `<= 8`.
2. Singular series for F1: derive the class ratios (`about 1 + 1/p` for classes `0` and `-1`) from the joint
   distribution of `(2 psi(h) + 1, 2h + 1)` modulo `p`; explain the weak primes `13, 19, 23`.
3. Decide whether THM-5122's width bound (Remark 4 of THM-5124) survives non-strong closures and the
   unit-Lipschitz extension.
4. The mass profile of the O/Z/J readout itself (not F1), to see whether its 29 misses are the light integers of
   its own linear forms.
