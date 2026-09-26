# Wave 13 orchestrator findings: Gersonides's four cycles, the identity 3 + 1 = 4, cheap edits under positive drift, and the Kuratowski–Tutte reading

**Status.**
- **PROVED (elementary; orchestrator).**
  - §1, the free cycles. An integer cycle of `T(x) = x/2, (3x+1)/2` on `Z` is *free* when every parity word of its shape `(p, a)` (length `p`, `a` odd letters) gives an integer periodic point. For shapes with `1 <= a <= p−1` this happens iff `|2^p − 3^a| = 1` (a shift argument); the single-word shapes `a = 0, p` give `0` and `−1`. By Levi ben Gershon's theorem (1343; CITED) the free cycles are exactly `{0}, {−1}, {1,2}, {−5,−7,−10}`, i.e. `2 − 1`, `3 − 2`, `4 − 3`, `9 − 8`. The fifth known cycle, `{−17, …}`, is *sporadic*: its shape `(11,7)` has `3^7 − 2^11 = 139`, and exactly one of the 30 necklaces of that shape is integral.
  - §1, the densities. The five densities `0, 1/2 | 1, 2/3, 7/11` are the first best lower and best upper rational approximations of `log_3 2`. `7/11` is the mediant of `2/3` (the `−5` cycle) and `5/8`, which is the critical Christoffel density of the provable strategies `sigma_k` of THM-4479 for `8 <= k <= 26`.
  - §1, the identity `3 + 1 = 2^2`. It plays three roles, each holding iff `q = 3` among odd multipliers:
    - pairing-sum preservation (THM-4470);
    - moment criticality `g_q(2) = (1+q)/4 = 1` (THM-4477);
    - the trivial cycle `1 -> 2 -> 1`.
    
    The integer roots of `g_3(s) = 1` are `s = 1, 2`, the free cycles `{0}` and `{1,2}`.
  - §2, positive drift. For odd `q >= 5` and every `W > 1`, the price of arbitrary fixed-horizon edits satisfies `eps_L(q) <= 1/W + P[Bin(L−1, 1/2) < ((L−1) + log_2 W)/log_2 q]`. So `eps_L(5) <= 2^(−0.0119 L + O(1))`, although the undecided density of `5n+1` stays positive (about `0.18`). THM-4478's capacity argument transfers verbatim and gives `eps_L(5) >= 2^(−0.0139 L − O(sqrt(L log L)))`.
- **FINITE-EXACT.**
  - The census of all words of length `p <= 20` finds exactly the five known integer cycles; the gates lane had `p <= 40`.
  - The catch-high construction of §2 was built and verified for `q = 5, 7` and `L = 8, 12, 16` on all sources `<= 4*10^5`.
- **ANALOGY / NUMEROLOGY (typed row by row).** §3, the Kuratowski–Tutte dictionary.
- **OPEN.**
  - ~~The exact exponent of `eps_L(5)`~~ SETTLED by [THM-4480](../../01-canon/theorems/THM-4480-peak-discounted-provability-price.md). Catching at the peak gives `eps_L(q) <= rho^peak_L <= (q/2) 2^(-(1-H(log_q 2))L)` for every odd `q`, so the exponent is exactly `0.013911` for `q = 5`. Proposition 4 below (catch at a fixed height `W`) is superseded; it remains a valid but weaker bound.
  - Whether periodic modifications of `5n+1` approach provability (wave-14 lanes).

Scripts: [gersonides check](../../04-computation/experiments/procgen_wave13_20260926_gersonides_check.py) → [output](procgen_wave13_20260926_gersonides_check.out); [positive-drift edits](../../04-computation/experiments/procgen_wave13_20260926_drift_arbitrary_edits.py) → [output](procgen_wave13_20260926_drift_arbitrary_edits.out). Session `collatz-procgen-20260922`, orchestrator, 2026-09-26.

## 1. Free and sporadic cycles

**Setting.** A parity word `w` of length `p` with `a` ones gives `T^p(x) = (3^a x + c_w)/2^p` on its 2-adic cylinder, where `c_(j+1) = 3 c_j + 2^j` at an odd letter `j` and `c` is unchanged at an even letter. Its unique periodic point is `x_w = c_w/(2^p − 3^a)` (Böhm–Sontacchi; THM-4471 §Banach). So `x_w` is an integer iff `(2^p − 3^a) | c_w`.

**Proposition 1 (free cycles).**
* **Statement.** The free cycles are exactly `{0}`, `{−1}`, `{1,2}` and `{−5,−7,−10}`.
  * `{0}` and `{−1}` come from the single-word shapes `(p,0)` and `(p,p)`, of which the primitive ones are `(1,0)` and `(1,1)`, i.e. `2 − 1` and `3 − 2`.
  * `{1,2}` and `{−5,−7,−10}` come from the mixed shapes with `|2^p − 3^a| = 1`, namely `(2,1)` and `(3,2)`, i.e. `4 − 3` and `9 − 8`.
* *Proof.*
  * **Single-word shapes.** For `a = 0` the only word is `0^p`, with `c = 0`, giving the fixed point `0`. For `a = p` it is `1^p`, with `c = 3^p − 2^p`, giving `−1`.
  * **Mixed shapes need `|D| = 1`.** Let `1 <= a <= p−1`. Choose a word of shape `(p,a)` whose last `1` sits at a position `j <= p−2`, and move that `1` to `j+1`. The shape is unchanged, and `c_w` changes by exactly `2^(j+1) − 2^j = 2^j`. If both words were integral, then `D = 2^p − 3^a` (odd) would divide `2^j`, forcing `|D| = 1`. Conversely, `|D| = 1` makes every word integral.
  * **Gersonides.** By Levi ben Gershon's theorem (1343), the only solutions of `|3^a − 2^p| = 1` with `1 <= a <= p−1` are `(2,1)` and `(3,2)`. Mihăilescu's theorem gives the general Catalan statement.
  * **Cycles.** Iterating the words gives `{1,2}` (from `10` and `01`) and `{−5,−7,−10}` (from the rotations of `110`).
  * FINITE-EXACT: all solutions of `|3^a − 2^p| = 1` for `p <= 400` (namely `(1,0), (1,1), (2,1), (3,2)`), and every word of these shapes. ∎
* **Sporadic.** For shape `(11, 7)`, `2^11 − 3^7 = −139`. Of the `330` words, exactly the `11` rotations of `11110111000` have `139 | c_w`, giving the cycle `−17 → −25 → −37 → −55 → −82 → −41 → −61 → −91 → −136 → −68 → −34`.
* **Census.** Among all words with `p <= 20`, the only integer cycles are these five (FINITE-EXACT; the gates lane extends this to `p <= 40`).

**Proposition 2 (Stern–Brocot position).**
* The densities of the five cycles are `0/1` and `1/2`, below `c = log_3 2 = 0.6309` (contracting, the positive side), and `1/1, 2/3, 7/11` above `c` (expanding, the negative side).
* `0/1` and `1/2` are best lower approximations of `c`. With denominators `<= 11`, the best upper approximations of `c` are exactly `1/1, 2/3, 7/11`.
* `7/11 = (2 + 5)/(3 + 8)` is the mediant of the Farey neighbours `2/3` and `5/8`, since `2·8 − 3·5 = 1`. So it is the fraction of least denominator in `(5/8, 2/3)`.
* `5/8` is `F_k` for `8 <= k <= 26`: the critical density of the provable strategies `sigma_k` of THM-4479, attained by the Christoffel orbit `319/13`.
* So the integer cycles of `3x+1` sit exactly where the continued fraction of `log_3 2` makes `|2^p − 3^a|` small:
  * `|D| = 1` gives the four free cycles;
  * `139` gives the one sporadic cycle.
  
  The next upper approximant is `12/19`, with `3^12 − 2^19 = 7153`. No cycle is known there, and none exists with `p <= 40`.

**Proposition 3 (`3 + 1 = 2^2` in three roles).** For odd `q`, the following are equivalent, and each holds iff `q = 3`:
* **(a) Pairing ladder (THM-4470).** `T(2i−1) + T(2i) = (2i−1) + 2i` for all `i`. Indeed `(q(2i−1)+1)/2 + i = (q+1)i − (q−1)/2`, and this equals `4i − 1` iff `q = 3`.
* **(b) Moment criticality (THM-4477).** `g_q(2) = (1+q)/4 = 1`.
* **(c) Trivial cycle.** `(q·1 + 1)/2 = 2`, i.e. `1 -> 2 -> 1` is a cycle.

Two further facts:
* **Continuous and discrete.** `g_3(s) = (1 + 3^(s−1))/2^s` is a smooth convex function of `s`. Its integer roots `s = 1, 2` are the solutions of `2^s − 3^(s−1) = 1`. These are exactly the free cycles whose word has a single even letter, `1^(s−1)0`: namely `{0}` (word `0`) and `{1,2}` (word `10`).
* **Other multipliers.** For odd `q >= 5` every solution of `|q^a − 2^p| = 1` has `a <= 1` (Mihăilescu; checked for `q <= 101`, `p < 400`). So `3` is the only multiplier with a free cycle of two odd steps, the `−5` cycle. `5x+1` has only the free cycles `{0}` and `{−1,−2}`, both from `5 − 4 = 1`.

**Reading.** The Pythagorean root `(3,4,5)` is `3 + 1 = 2^2 = 5 − 1`. The trivial cycle `1 -> 2 -> 1` exists for `3x+1`, where it contracts with multiplier `3/4`, and for `5x−1`, where it expands with `5/4`. So the root is the owner's "±1 sandwich" around the power of two `2^2`. The arithmetic is PROVED; the reading is ANALOGY.

## 2. Positive drift makes arbitrary edits exponentially cheap

**Setting (THM-4478).**
* `T_q(n) = n/2` or `(qn+1)/2`, with `q` odd.
* An `L`-step descent modification is any `G : N -> N` such that every `n >= 2` has `G^j(n) < n` for some `1 <= j <= L`.
* `eps_L(q)` is the infimum over such `G` of the upper density of `E(G) = {v : G(v) != T_q(v)}`.
* For `q = 3`, THM-4478 proves `eps_L(3) = 2^(−(1−h)L + o(L))`, which matches the undecided density `rho_L`.
* For `q = 5` the undecided density stays positive: `|Bad_L(5)|/2^L = 0.260, 0.234, 0.228, 0.216` at `L = 10, 14, 18, 22`, tending to about `0.176`.

**Proposition 4 (catch high).**
* **Statement.** Let `q >= 5` be odd and `W > 1`. Then
  `eps_L(q) <= 1/W + P[ Bin(L−1, 1/2) < ((L−1) + log_2 W) / log_2 q ]`.
  With `W = 2^(theta L)` and `theta = theta_q` the root of `theta = 1 − H(c_q(1+theta))`, this gives `eps_L(q) <= 2^(−theta_q L + O(1))`. The values are `theta_5 = 0.011921`, `theta_7 = 0.047071`, `theta_9 = 0.075784`.
* **Construction.**
  * Let `B` be the set of `n >= 2` with `T^j(n) >= n` for all `j <= L`.
  * For `n` in `B`, let `tau(n)` be the least `j` in `[1, L−1]` with `T^j(n) >= W n`, if there is one.
  * Put `E_high = {T^(tau(n))(n)}` and `E_low = {n in B : tau(n)` undefined`}`.
  * Define `G = 1` on `E = E_high ∪ E_low`, and `G = T_q` elsewhere.
* *Proof.*
  * **`G` is an `L`-step descent modification.** A `G`-orbit follows `T_q` until it first meets `E`, and then goes to `1 < n`.
    * A source `n` outside `B` either descends by time `j <= L` or meets `E` earlier. In the latter case it reaches `1` by time `j`.
    * A source in `E_low` goes to `1` at time `1`.
    * A source in `B` with `tau(n)` defined meets `E` by time `tau(n) <= L−1`, and so reaches `1` by time `L`.
  * **`E_high` has upper density `<= 1/W`.** Each point `v` of `E_high` has a source `n <= v/W`. So `|E_high ∩ [1,Y]| <= Y/W`.
  * **`E_low` is small.** Write `T^(L−1)(n) = w (n + h)` with word slope `w = q^e/2^(L−1)` and affine offset `h >= 0` (THM-4478 §2). Then `n` in `E_low` forces `w < W`. By the Terras bijection, `{n : w(n) < W}` is a union of residue classes mod `2^(L−1)` of density `P[q^e < W 2^(L−1)]`, where `e ~ Bin(L−1, 1/2)`.
  * **The exponent.** Chernoff, `P[Bin(m,1/2) <= xm] <= 2^(−m(1−H(x)))` for `x < 1/2`, balanced against `1/W`. ∎

**Proposition 5 (the lower bound transfers).**
* **Statement.** THM-4478 Theorem A holds for every odd `q`, with `h_k <= k/q`, `R_L = sum_(k<L)(floor(k/q) + 1)` and `M_L(K) = (floor(log_q K) + 1) R_L`. Its block construction of §4 uses only that `c_q` is irrational. Hence `eps_L(q) >= 2^(−(1 − H(c_q))L − O(sqrt(L log L)))`, and for `q = 5` the exponent is `0.013911`.
* *Proof.* The proof of THM-4478 §§2–4 is read with `1/q` in place of `1/3`:
  * the offset bound `h_k <= k/q` holds on band words (all slopes `>= 1`);
  * the number of admissible odd counts `e` for a given `(v, k)` is `floor(log_q K) + 1`;
  * the cycle-lemma rotation of blocks with `ceil(c_q b)` ones keeps every prefix slope in `[1, q^(t+b+r)]`. ∎

**Consequences.**
* For `q = 5`, `2^(−0.0139L − o(L)) <= eps_L(5) <= 2^(−0.0119L + O(1))`. The price decays exponentially, while the undecided density does not decay at all.
* **What decides the price is the critical band, not the undecided set.** These coincide exactly when the drift is negative (`q = 3`). There, conditioning on non-descent confines a source to the band, so THM-4478's price equals `rho_L`.
* **Two places.** The catch-high construction uses *height*: it edits an orbit where it is already `W` times higher than it started, and such points are sparse. A periodic edit (a union of residue classes, i.e. a 2-adic object) cannot see height. This is the same archimedean-versus-2-adic split as the S5 thin-divergence work (THM-4476, `R(d)` against `R_2(d)`).
* **Catching at the peak (refinement; digest proposal, lane `peak` auditing).** Catching at the first crossing of a fixed `W` is not optimal. Send each bad source to 1 at the highest point of its trajectory before the horizon instead. That costs `rho^peak_L(q) = 2^(-L) sum_(u in Bad_L) 1/w*(u)`, where `w*` is the peak word slope. THM-4478's capacity bound, stratified by peak height, gives the same quantity from below up to `poly(L)`.
  * For `q >= 5`, `sigma + (1 - H(c_q(1+sigma)))` is increasing in `sigma >= 0`: its derivative is at least `1 - c_q log_2((1-c_q)/c_q) > 0`. So the peak sum is dominated by the band words, and the exponent becomes exactly `1 - H(log_q 2)`: `0.013911` for `q = 5`.
  * For `q = 3` the discount changes only the subexponential factor, conjecturally `exp(-Theta(L^(1/3)))`.
* **Wave-14 lanes.** They test the periodic versions:
  * deletions via Golomb–Mykkeltveit feedback sets, where `Theta(1/k)` is expected;
  * sign flips, where the necklace bound is `>= 2/k` and whether the price tends to 0 is open.
* **Checks** (`procgen_wave13_20260926_drift_arbitrary_edits.out`).
  * For `q = 5, 7`, `L = 8, 12, 16` and `W in {1.5, 2, 4, 8, 2^(theta L)}`, `G` makes every `2 <= n <= 4*10^5` descend within `L`.
  * The measured edit densities are far below the undecided density: `0.0149` against `0.234` at `q = 5, L = 16, W = 8`.
  * The bound decays like `2^(−0.0117L)` at `L = 800`.
  * The transferred lower bound stays below the upper bound for `L <= 24`.

## 3. The Kuratowski–Tutte reading (the owner's triple)

**The triple.**
* **Kuratowski/Wagner.** Planarity is characterized by two excluded graphs, `K_5` and `K_{3,3}`.
* **Tutte.** Tutte's 4-flow conjecture excludes a third graph, the Petersen graph (cubic case proved by Robertson, Sanders, Seymour and Thomas). It is the smallest snark, and it contains both Kuratowski graphs as minors.
* **Why these three.** Each is non-planar for the counting reason `E > g(V−2)/(g−2)` at its girth `g = 3, 4, 5`: `10 > 9`, `9 > 8` and `15 > 40/3`.
* **Tutte's matroid version.** Regular matroids are characterized by excluding `{U_{2,4}, F_7, F_7*}`, a self-dual matroid and a dual pair. Graphic matroids further exclude `M*(K_5)` and `M*(K_{3,3})`.

| Kuratowski–Tutte | Collatz (this session) | Status |
|---|---|---|
| planarity characterized by excluded minors | bounded-lookahead provability characterized by excluded substructures: class (i) iff no expanding cycle in the parity graph (THM-4474 A) | both PROVED; the parallel is structural ANALOGY |
| Euler's formula: the counting reason behind the obstructions | the entropy count `h(log_3 2)`: `\|Bad_k\| <= 2^(hk)`, necklaces `>= 2^(hk)/(3k^2)`, integer capacity (THM-4478/4479) — the counting reason behind the price of provability | ANALOGY |
| `{K_5, K_{3,3}}`, forced by the counting bound, plus Petersen, sporadic, containing both | `{−1}, {−5,−7,−10}`, forced by `3 − 2 = 1`, `9 − 8 = 1` (Gersonides), plus `{−17,…}`, sporadic (`139 \| c_w`), whose word `11110111000` contains the words `1` and `110` | classification PROVED (§1); correspondence ANALOGY |
| every non-planar graph contains one of the excluded graphs | every class-(i) strategy, at every level, must break all three expanding integer cycles (each is a closed walk of every parity graph `G_0`); `−1` is forced by Proposition F, and every optimum also flips `−5` | PROVED (THM-4474, THM-4479 data) |
| Tutte's `{F_7, F_7*}` (dual pair) + `U_{2,4}` (self-dual) | level-2 cube: Collatz and 3n−1 (a `nu`-dual pair, both class (iv)) + `chi_(−4)` (self-dual, the unique provable) + `−chi_(−4)` (self-dual, maximally expanding) | classes PROVED; correspondence ANALOGY |
| Tutte's flow–colouring duality (flows vs tensions) | expanding cycles (positive circulations) vs potentials (Lemma P of THM-4479's note; Karp/Hoffman duality); lane `tension` develops the rank-function form | duality PROVED in form (LP) |
| Collatz's Schreier graph under Kohl's three involutions | Tait-coloured by construction, hence not a snark; its Kempe chains are the doubling orbits and the rising runs (wave 10) | PROVED |
| flow numbers `2, 3, 5` of `K_5, K_{3,3}`, Petersen | multipliers `2` (halving), `3` (Collatz), `5` (the drift control) | NUMEROLOGY |
| girth `3, 4, 5`; the Pythagorean root `(3,4,5)` | `3 + 1 = 2^2 = 5 − 1`: the trivial cycle of `3x+1` (contracting) and of `5x−1` (expanding) | arithmetic PROVED; reading ANALOGY |

**What the reading buys.**
* The Kuratowski pattern is "a property is characterized by excluded substructures, and a counting identity forces the smallest ones". This is exactly how bounded-lookahead provability works.
  * **Theorem A** supplies the excluded substructures.
  * **Gersonides** supplies the forced smallest ones: `−1` and `−5`, the two free expanding cycles.
  * **Entropy** supplies the count, and hence the price.
* **The Petersen analogue, `−17`, is sporadic.** It is not forced by any identity. It sits at the next Stern–Brocot approximant. It is barely expanding (`3^7/2^11 = 1.068`), and every provable strategy must still break it.
* **Where the analogy stops.** Kuratowski's list is finite, whereas the expanding cycles form an infinite family: the necklace count grows like `2^(hk)`. So there is no finite obstruction set, as the reframe lane also found. The finite part is the *integral* obstructions, and there the question "are there only three expanding integer cycles?" is the negative-side cycle problem for `3x+1`, i.e. the `3x−1` cycle conjecture.

## 4. Discrete ↔ continuous: the session's ledger as of wave 13

1. **Integer capacity versus moments.**
   * Distribution-only (Haar/moment) arguments lose a square: exponent `2(1−h)` (THM-4477).
   * Actual integer spacing, via the affine offset `h_k in [0, k/3]`, recovers it: exponent `1−h` (THM-4478).
   * Periodicity recovers it combinatorially, through disjoint necklaces (THM-4479).
2. **Christoffel words.** The critical cycles of the provable approximants `sigma_k` are upper Christoffel words, i.e. discrete lines of slope `F_k -> log_3 2`. The five integer cycles sit at the first best approximants (§1).
3. **The moment curve.** The continuous curve `g_3(s)` meets 1 at the integer points `s = 1, 2`, which are free cycles (§1).
4. **Potentials versus cycles.** A real potential certifies the absence of expanding cycles (LP duality). Bernoulli-boundary's obstruction to a periodic every-step rank for Collatz is the loop at `−1` (lane `tension`).
5. **Two places.** Periodic (2-adic) modifications cannot use height, whereas arbitrary (archimedean) modifications can. Under positive drift this is the difference between polynomial and exponential price (§2 and the wave-14 lanes).
