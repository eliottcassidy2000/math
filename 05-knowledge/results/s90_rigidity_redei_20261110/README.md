# S90 (opus, 2026-11-10): additive rigidity, the residue-tournament model, and where the mass is small

**Status.** PROVED: THM-5124 (additive rigidity of the Hamiltonian-path count on presence cubes). THM-5125: the
Legendre (quadratic-residue) tournament is KNOWN in the repository and the literature; this session adds the
residue-chain reformulation, the Heegner data and least realisers (FINITE-EXACT) and two cubic balance criteria
(PROVED). NUMERICAL with three pre-registered out-of-sample tests, all CONFIRMED: HYP-9657 (A)-(B), restricted
after the audit. FINITE-EXACT refutations: HYP-9657 (C). Independent adversarial audit: PASS WITH CORRECTIONS
([audit/AUDIT.md](audit/AUDIT.md); 5 major and 22 minor findings, all in attribution, over-general wording and
scope, all applied; MISTAKE-752).
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
  F1 mass (THM-4875 (D), THM-4885) with its local factors (S82 note section 5); the Monsky matrix of the
  reciprocity tournament (bsd_selmer_growth_20261013a); klein's Legendre tournament and its identification with
  the random tournament (QR_TOURNAMENT.md, LEGENDRE_SIDECAR.md, HYP-9375); THM-467 (`H` and `I(Omega, t)` as
  determinant ratios: the odd cycles are the odd, arctanh-type part of `-log det(I - XA)`). The last two were
  missed at the start of the session and found by the audit (MISTAKE-752).
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

*Why it matters for the target.* For one new vertex over a fixed core there are three cubes with three signs
(Remark 2 of the theorem): on the presence cube interactions are `>= 0` and vanish only as above; on the ear cube
(single cut bits) they are strictly negative (THM-4114, THM-5122); along exchange directions (a selected and an
unselected vertex swapped) they are signed and can cancel - THM-5122 gives the criterion, a capacity bound and one
2-cube, THM-5117 gives affine exchange cubes of every dimension `k` at order `7k + 1`. Together: **additivity of
`H` beyond `1, 3, 3, 5` is an orientation (exchange) phenomenon, never a presence phenomenon.** Every proof
architecture that needs three independent additive contributions from the presence of vertices - Brown's
criterion, Erdos-Lewin sums of 3-smooth numbers over antichains, a sumset of two independent rich families,
affine presence cubes - is ruled out, not merely unfound (one_apex_join_representation_mass_pi2_over_12_20261106
section 1, route 2, had recorded that no free-staircase family was found). Remark 4 of the theorem adds a width bound for
THM-5122's orthogonal pair-exchange cubes with strong closures: an interval obtained from one such cube at baseline
`B` has length at most about `B^0.65`.

## 2. The Legendre tournament (THM-5125)

On the primes `= 3 (mod 4)`, "`p -> q` iff `p` is a quadratic residue mod `q`" is a tournament by reciprocity; by
Dirichlet and CRT it has the extension property, so it is the countable random tournament and every finite
tournament is `T_n` for some squarefree `n`. All of this is KNOWN: it is klein's Legendre tournament
(QR_TOURNAMENT.md, 2026-10-25), identified with the random tournament in LEGENDRE_SIDECAR.md (2026-10-26, citing
Cameron and Jones), with quasirandom prime prefixes (HYP-9375). What this session adds is the reading of `H`:
`H(T_n)` is the number of *residue chains* of `n` (orderings of the prime factors with each a residue modulo the
next), Redei's theorem says this number is odd, and

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
   odd power sums `p_1, p_3, p_5, ...` with nonnegative coefficients (Grinberg-Stanley), and `H` is its value at the
   one-variable alphabet `x = (1, 0, 0, ...)` (every `p_k = 1`), `H = sum over odd packings 2^|P|`. Map: cycle of
   length `k` -> `p_k`. Preserved: disjointness (products). Lost: which cycles intersect. This is KNOWN in the
   repository in a sharper form: THM-467 (June 2026) writes `H(T) = eps(det(I + XA) / det(I - XA))` and identifies
   the odd-cycle generating function with the odd, arctanh-type part of `-log det(I - XA)` ("even lengths cancel,
   odd double"; also LEM-004). So the "odd exponents" of the directive have an exact form: the odd part of the
   logarithm of a characteristic polynomial, i.e. the sign character of cycles, not the parity of an exponent of an
   arbitrary function. Earlier notes (collatz_mod6_20260917_*, seam_power_coordinates_20260925) recorded "no map"
   for the word *odd* alone, which stands.
2. **"Cubic and linear."** The curves the repository connects to tournaments are `y^2 = x^3 - n^2 x`, an odd cubic
   in `x` (extra symmetry `x -> -x`, CM by `Z[i]`). For `n` a product of primes `= 3 (mod 8)` Monsky's formula
   (CITED) reads the pure 2-Selmer dimension off the reciprocity tournament as `nullity(M^2 + M + I)`. THM-5125 (4)
   makes the comparison exact: per(I + A), det(I + A), `H` and this nullity condition are evaluations of one
   cycle-packing enumerator (the first three KNOWN: THM-467, THM-1950; the balance criteria are from this session).
   `H` uses the *parity* grading with weight `2` over `Z`; the Selmer condition uses a *cubic* grading
   (`x^2 + x + 1`) with coefficients in `F_2`: it holds iff the numbers of packings in the three classes of `e - o`
   modulo `3` have equal parity. The cubic enters at matrix level - for primes `= 3 (mod 8)` the Monsky kernel is an
   `F_4`-space on which the cube root of unity acts by `A` (repository note; stability under `(x, y) -> (y, x + y)`
   checked by the audit) - and no automorphism of the curve is involved. So "cubic" belongs to the Selmer side and
   "odd" to the Hamiltonian side: two gradings of one family of cycle packings, not one statement. Test for more:
   HYP-9657 (C) - no congruence for `H` follows from the `F_2` data at small orders.
3. **Heegner numbers.** Two exact facts and one mechanism. (i) Domination, above (a restatement of classical
   inertness). (ii) Class number one forces `q` and, for `q >= 7`, `(q + 1)/4` to be prime (Rabinowitsch, CITED).
   (iii) The F1 join family has two sheets: `N = x1 x2 + 2^s y1 y2` (`s >= 1`) and `N + 1 = x1 x2 + y1 y2` (`s = 0`).
   At the small primes `p <= 23` (and `41, 43`) the classes `N = 0` and `N = -1 mod p` carry more mass, the first
   excess coming from the first sheet and the second from the second (HYP-9657 (A); the local factor at `N = 0` and
   the reduced mass on primes were already in the S82 note, section 5; the excess is absent at
   `p = 29, 31, 37, 53, 61, 73, 89, 127` and its size is not `1/p`). So integers with `N(N + 1)` free of small odd
   primes are lighter (mean ratio 0.72, miss-rate ratio 2.6-3.3, pre-registered and confirmed), and the Heegner
   numbers `7, 43, 67, 163` are such integers (`11` and `19` are not): `r(7) = r(43) = r(67) = r(163) = 0`. This is
   a statistical tendency that makes those four likely misses of F1; it is not a derivation. It does not
   explain 21 (`= 3 * 7`, not rough) and it is not a class-number formula: where the mass of a representation
   problem *is* a class number (sums of three squares: `r_3(d) = 24 h(-d)` for `d = 3 mod 8`, `d > 3`, Gauss, CITED;
   so `r_3 = 24` at `11, 19, 43, 67, 163`, the least possible), the Heegner numbers are extremal for a different
   reason. Type: MECHANISM (statistical), scope F1.
   *Scope check on a second family* (exploratory, in-sample, `ozj_mass.out`): for the three-letter system `S, B, F`
   (the O/Z/J grammar) the number `m(N)` of rooted states with count `N` has class means `101 / 56 / 59` modulo `3`
   and `95 / 67 / 64 / 67 / 68` modulo `5` on `[3000, 12000]`: the class `N = 0` dominates (products), the class
   `N = -1` is only faintly favoured, primes carry `0.68` of the mean mass and Rabinowitsch-rough integers `0.65`.
   There the Heegner numbers are light because they are primes; the role of `N + 1` is specific to F1.
4. **Ramanujan.** (a) Ramanujan's constant `e^(pi sqrt 163)` is the class-number-one fact `j(tau_163)` in `Z` again
   (see ramanujan_heegner7_trivial_cycle_20261006.md and MISTAKE-568 for attribution); nothing new is claimed.
   (b) Ramanujan's ternary form `x^2 + y^2 + 10 z^2` misses the odd numbers `3, 7, 21, 31, 33, 43, 67, 79, 87, 133,
   217, 219, 223, 253, 307, 391, 679, 2719` (all eighteen are unconditional exceptions - the last two were added by
   Jones-Pall and Gupta; that the list is complete is Ono-Soundararajan 1997 under GRH; finiteness is
   Duke-Schulze-Pillot 1990). The O/Z/J readout misses 29 odd numbers, largest 2647 (FINITE-EXACT through 50,001,
   THM-4884). The lists share
   `7, 21, 43, 67`. Under a uniform null on the odd numbers below 100 (nine members each) an overlap of at least
   four has probability about 0.04; both lists favour integers of small local mass, so the null is too weak. Type:
   ANALOGY of kind (a non-regular family with a finite sporadic exceptional set), COINCIDENCE of members. The
   mechanisms differ: there, mass = class number + cusp-form coefficient; here the data neither establish nor
   exclude a Kronecker-symbol dependence (HYP-9657 (A)).
5. **Birch and Swinnerton-Dyer.** The only map is Monsky's (item 2). As a dictionary, typed ANALOGY throughout:
   the conflict polynomial `I(Omega(T), x)` is multiplicative over strong components (an Euler product,
   THM-1862/1926), `H` is its value at `x = 2`, its value at `-1` is an Euler characteristic and equals `H` mod `3`
   (KNOWN: THM-467, THM-4754, THM-4756); the spectrum problem is local-global with the single local condition at `2` (Redei), and
   `{7, 21}` plays the role of a finite obstruction set visible only globally, removed by a descent (strong floor
   plus deletion, THM-1370); the "additive rank" of `H` is at most `2` (THM-5124). None of these is a theorem about
   elliptic curves.

## 4. Bridges tested and refuted (first failed implication recorded)

- `H mod 4` (or `mod 8`) as a function of the Redei/Monsky `F_2` data: true at orders 3-4, false at order 5.
  First failed implication: equal (Laplacian corank, Monsky nullity, nullity of `I + A`) does not give equal number
  of odd cycles modulo 2. Survivor: Redei's parity and the identities of THM-5125 (4). The correct mod-4 statement
  was published the same day by mac-mini, THM-5127 (C): on the ear cube of a fixed core `R`, the map
  `X -> (H(R + v_X) - H(R))/2 mod 2` is a quadratic form over `F_2` with polar form `mu_R / 2 mod 2`, so the mod-4
  imbalance over the `2^m` cuts is `0` or `+-2^(m - rank/2)`: the mod-4 information lives in a quadratic form on
  cuts, not in the ranks of the adjacency matrix.
- "Monsky nullity `> 0` implies `3 | H`": true at orders 3-4, false at order 5 (`H = 13`, the five-vertex
  tournament with out-neighbour masks `[18, 8, 3, 5, 14]`, which is Session245's h13 cell). With the pure adjacency
  matrix it survives to order 5 and fails at order 6 (`H = 17`). Survivor: for a single triangle both hold.
- A general law "classes `N = 0, -1 mod p` are favoured for every odd prime, the rest equal, excess about `1/p`"
  for F1: written in the first version of HYP-9657 from eight primes and REFUTED by the audit at the first untested
  primes (`29, 31, 37, ...`; the classes are not equal at `p = 7`; the excess is not `1/p`). Survivor: the
  pre-registered statements and the sheet decomposition at small primes. A class-number-type (Kronecker symbol)
  mass is neither established nor excluded.
- Sum of two independent rich families (Mann-type sumset of one-apex sets, which have density about one half):
  impossible by THM-5124 (2).
- A second-order gap. For two active vertices with increments `d_a, d_b` over a base of count `c` and interaction
  `mu > 0`, the bound `2 mu >= max(d_a, d_b)` holds on every labelled tournament of orders 4-7 (equality on 14
  profiles, e.g. `(c, d_a, d_b, mu) = (1, 4, 2, 2)`, `(5, 4, 4, 2)`, `(11, 8, 8, 4)`) and FAILS at order 8:
  out-neighbour masks `[246, 252, 248, 177, 160, 64, 152, 32]`, `a = 1`, `b = 2`, counts `15, 29, 29, 49`, so
  `(15, 14, 14, 6)` (`interaction_gap_refuted.out`, recounted over permutations); the witnesses found by the audit
  at orders 8, 9, 10 have ratios `3/7, 5/13, 9/25 = (2^j + 1)/(3 * 2^j + 1)`, suggesting the limit `1/3`.
  First failed implication: "vanishing or large" does not extend from the zero set to a multiplicative gap.
  Survivor: THM-5124 (the zero set) and Session214's floor. A pattern exhaustive through order 7 died at the
  first untested order.

## 5. Where this leaves the proof

1. *Dead:* presence-additive architectures (THM-5124); the single-realiser induction through `[h, 2h - 1]`
   (Session214, all orders); one-hub grammars (Session 16, value 13). *Unpromising, not dead:* single-realiser
   windows `[2h + 1, 3h]` (unavailable through order 8, THM-4876 (5), FINITE-EXACT); finite banks as a proof pattern
   (conjectured to miss a positive proportion, HYP-9480 (C), OPEN).
2. *Alive, with a quantified limit:* exchange cubes (THM-5122 criterion, THM-5117 cubes of every dimension) - width
   at most about `B^0.65` per cube with strong closure, so a covering needs polynomially many baselines per octave;
   moving cores (THM-4960, THM-5000). Published the same day by kind-pasteur, THM-5126: every odd `N > 1` lies on
   one 2-adic line `1 + c 2^s` (`c` the odd part of `N - 1`), unit-offset roots put whole rays of a line into the
   spectrum, and all 65 lines with `c <= 129` lie in it except for `7` (the only miss of `c = 3`) and `21` (the
   only miss of `c = 5`): an all-order family organised by the 2-adic coordinate, the natural partner of the
   threshold count in item 5.
3. *The analytic route.* The three-letter rooted system `S, B, F` on states `(h, R, e, u)` (it is the O/Z/J grammar;
   rules re-derived and checked against brute force in `psi_sbf.py`) satisfies `h <= R` at every state and `R` at
   most doubles with each letter, so the `2 * 3^(k-1)` words of length `k` not ending in `S` all have count at most
   `2^k`: on average at least `(3/2)^k` words per value. In this sense it is a fat orbit, not a thin one
   (word-growth exponent `log 3 / log 2 = 1.58 > 1`); its readout nevertheless misses 29 odd values (THM-4884). For
   such orbits the available technology (bilinear circle method with a spectral gap, as in Bourgain-Kontorovich
   for Zaremba; a power-saving exceptional set is the Magee-Oh-Winter refinement) yields density one, not all
   values; this typed ANALOGY is already in the repository (seam_gadget_thickness_spectrum, Session205) and is
   repeated here only to state the realistic milestone: positive density of the spectrum is still unproved (best:
   `A(X) >= X^(2/3) / (60 (ln X)^(2/3))`, Session216 `density/RESULT.md`, improving THM-4887).
4. *What would finish.* Session246's two routes stand (a family of cores with exchange systems whose absolute
   intervals cover every target; or a well-founded physical predecessor rule). THM-5124 adds a constraint to the
   first (no presence switches) and Remark 4 a size requirement.
5. *A third route, with its threshold (PROGRAM; nothing proved beyond the bookkeeping).* Write a family `F` of
   tournaments with counts at most `X`, `r(N) = #{T in F : H(T) = N}` and `S(theta) = sum_{T in F} e(theta H(T))`.
   - *Threshold (heuristic).* A circle-method proof that every odd `N` in the bulk of `F` is attained needs the
     frequencies away from `0` and `1/2` to contribute `o(|F| / X)`. In every family measured the sum has
     square-root size at generic frequencies (next bullet), i.e. `r` fluctuates at least at Poisson size; then a
     pointwise bound of that strength is possible only if `|F|^(1/2)` is at most about `|F| / X`: the family must
     have more than `X^2` members. Below that threshold only mean-square information is available and the output
     is density one at best (the Bourgain-Kontorovich situation, and the reason the 2026-11-06 note called
     pointwise completeness "at least as hard").
   - *Which families pass.* One apex over a chain: `|F| = 2^n`, `X <= 2^n`: exponent `1`. F1 (two words): exponent
     `1`. The `S, B, F` system: `2 * 3^(k-1)` words, counts `<= 2^k`: exponent `1.58`. `k` apexes over a chain of `n`
     vertices: `2^(kn)` cut patterns and every count is at most `k! (k+1)^n` (a Hamiltonian path orders the apexes
     and distributes the chain vertices over the `k + 1` gaps, increasing inside each gap), so the exponent is at
     least `k / log2(k+1)`: `1.26, 1.5, 1.72, 1.93, 2.14` for `k = 2, ..., 6`. The `k`-apex family passes from `k = 6`.
   - *What the sums look like* (HEURISTIC data, `apex_sums.out`): at 16 random frequencies the median of `|S|` is
     about `2 |F|^(1/2)` for one apex (`n <= 22`, exponent `log|S| / log|F|` falling from `0.60` to `0.54`) and for two
     apexes (`n <= 11`, `0.65` to `0.54`): square-root size at generic frequencies. At special rationals the decay
     is much slower (one apex: rate `0.932` per letter at `1/7`, 2026-11-06 note); those points belong to the major
     arcs.
   - *The tool at hand.* For one new vertex over a fixed core the count is the ear quadratic
     `E(x) = c0 + L.x - (1/2) x^T W x` (THM-4114), and Weyl differencing on the cube is exact:
     `|S(theta)|^2 = sum over d in {-1,0,1}^n of e(theta (L.d - (1/2) sum_{i in supp d} (W d)_i)) prod_{i not in supp d} (1 + e(-theta (W d)_i))`
     (`ear_weyl_identity.py`, checked). So `|S|^2 <= sum_d prod_{i not in supp d} 2 |cos(pi theta (W d)_i)|`, and a
     power saving follows whenever `theta W d` has many coordinates far from the integers. The same differencing
     bound, in a decoupled form over two disjoint vertex sets, is THM-5127 (C4) (mac-mini, published the same day at
     12:11, before this note was pushed; found here independently): cite THM-5127 for the inequality.
   - *What is missing.* One cube can save at most about `2^(-n/2)`; the threshold asks for `(k+1)^(-n)`, so the
     saving must be harvested from about `2 log2(k+1)` apexes at once without the `2^(-k)` loss of iterated
     Cauchy-Schwarz; the resonant rationals (Mersenne-type denominators, where the one-apex rate is `0.89` to
     `0.93` per letter) must be handled as major arcs with vanishing local factor; and the major arcs need the
     limit law of `H / X` together with equidistribution modulo `q` (THM-4886 (C) for one apex). Type: PROGRAM.
     It is recorded because it is the one route in this thread whose exponent count allows a statement about
     *every* large `N` rather than almost all; it is not a claim that the route works.

## 6. Principles, named

These are names for mechanisms that recur in this thread; each has a trigger and a counterindication. They are
recorded here, not promoted to META-PATTERNS (the evidence is from one problem family).

- **The cycle gas** (the oldest idea in the tournament work of this repository: the odd-cycle formula, THM-002,
  March 2026; the name is Session37's, and the determinant form is THM-467). Count by families of disjoint odd
  cycles at fugacity `2`, not by path words. Every advance on this
  problem is a reading of that gas: parity is its constant term (Redei); order joins multiply because far
  components are independent (THM-1862); deletion removes the atoms through a vertex (Claim A); the presence cube
  is totally monotone because activities are positive (THM-4114); increments have a floor (Session214); the
  interaction of two sets is the weight of the packings linking them, and vanishes only in `1, 3, 3, 5`
  (THM-5124); changing the fugacity gives the other bases (THM-4940) and fugacity `-1` gives `H` mod `3`
  (THM-4754).
  Trigger: any statement about `H` of nested or glued tournaments. Counterindication: orientation changes (ear and
  exchange cubes) are signed and need the quadratic form, not the gas (THM-5122).
- **Two gradings of one cycle gas.** Trigger: a count over disjoint cycles appears next to a determinant or a
  rank. Action: write both as evaluations of the cycle-packing enumerator and compare gradings (parity over `Z`
  for `H`; modulo `3` over `F_2` for Monsky). Counterindication: equal small values are not a congruence
  (section 4).
- **Root the object.** Composition is linear in rooted states `(h, R, e, u)` and opaque in the scalar `H` (one-arc
  join law, tail law, the `S, B, F` system; Session246: "keep the arithmetic fiber"). Trigger: a scalar induction
  stalls. Action: find the root whose state closes the recursion; record which coordinate a move forgets
  (`S` forgets `e`, `F` forgets `u`). Counterindication: equal scalars with unequal root menus (THM-5117).
- **Additivity beyond 1, 3, 3, 5 is exchange, never presence.** Trigger: a selector design that adds vertices to
  add to `H`. Action: check the interaction; by THM-5124 it vanishes only in `1, 3, 3, 5`; single cut bits are
  strictly subadditive; use exchange directions (THM-5122, THM-5117) or carries.
- **Mass, not congruence, decides the exceptions** (the local factors are S82's observation; new here are the
  second sheet and the rough stratum). Trigger: an exceptional list contains famous numbers. Action: compute the
  class means of the representation count and its dependence on the roughness of the natural linear forms (`N`,
  `N + 1` here) before reading the list, and test primes that were not used to form the statement.
  Counterindication: families of mass below one; a pattern seen at eight primes failed at the ninth (section 4).
  Related cards: "Separate local support from bounded-height coverage", "Existence is a maximum or tail question".
- **A height with a gap.** `log H` is superadditive over vertex partitions and additive on order joins
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
    python3 f1_formula_check.py  > f1_formula_check.out      # the F1 formula against explicit joins
    python3 ozj_mass.py          > ozj_mass.out
    python3 interaction_gap_refuted.py > interaction_gap_refuted.out
    python3 ear_weyl_identity.py > ear_weyl_identity.out
    python3 apex_sums.py         > apex_sums.out             # about 5 min, needs numba

Every script exits nonzero on a failed check. Universes are stated in each docstring.

## 9. Next obligations

1. The interaction of two active vertices outside the rigid configuration: the naive gap `2 mu >= max(d)` is false
   from order 8 (section 4). Is `inf mu / max(d)` positive? The minima are `1/2` through order 7 and the witnesses
   of orders 8, 9, 10 have `3/7, 5/13, 9/25 = (2^j + 1)/(3 * 2^j + 1)`: is the infimum `1/3`, that is, does
   `3 mu > max(d)` hold at every order whenever `mu > 0`? OPEN; the pattern `2 mu >= max(d)` was also exhaustive
   before it failed.
2. Singular series for F1: derive the sheet excesses of the classes `0` and `-1` from the joint distribution of
   `(2 psi(h) + 1, 2h + 1)` modulo `p`, and explain which primes carry an excess. The measured sizes (HYP-9657 (A))
   are not `1/p` and not monotone: `0.56 / 0.49` at `p = 3` down to `0.18` at `11`, about `0.04` at `13, 19, 23`,
   `0.26 / 0.23` at `17`, about `0` at `29, 31, 37`, positive again at `41, 43`. (An earlier wording of this item,
   "about `1 + 1/p`", repeated the statement withdrawn by the audit, M3; corrected after close-out.)
3. Decide whether the width bound (Remark 4 of THM-5124) survives non-strong closures.
4. The program of section 5, item 5: a sup-norm bound for the ear-cube sum of one apex over a fixed core, uniform
   in the core, from the exact differencing identity; then two apexes.
