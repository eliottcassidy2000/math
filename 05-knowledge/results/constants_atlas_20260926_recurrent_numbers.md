# Recurrent numbers across the repository: a typed atlas (what recurs for a reason, what recurs by accident, and the one family that became a theorem)

**Status: CENSUS (lexical, FINITE-EXACT counts over 9,437 labelled files) +
TYPED VERDICTS (modelling judgements, each with its stated map or its
stated absence) + one PROVED unification (the entropy curve of `3n±1`,
[THM-4487](../../01-canon/theorems/THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md)).
A shared number is a lead, never a connection: every entry below names the
mechanism in each thread and says whether a map exists. Session
`collatz-exponent-atlas-20260926` (opus), 2026-09-26.**

Scripts: `04-computation/experiments/constants_atlas_20260926_mine.py`
(census; JSON `constants_atlas_20260926_mine.json`) and
`constants_atlas_20260926_context.py` (the sentence carrying each value,
per thread). Threads are keyword votes (collatz, lrc, jc, nc2, tournament,
hadamard, amm, pell, kakeya, erdos, additive, abc, bernoulli,
transcendence, mahler, graph, zeta, sequences); a file carries its top two
labels. Navigation and atlas files are excluded because they aggregate every
thread. Filters drop ubiquitous small numbers (decimals need four
significant digits; fractions need numerator `>= 3` or denominator `>= 9`).

## 0. What the census shows

291 constants occur in at least two files of each of at least four threads.
Read with their contexts they fall into three kinds:

1. **Structural within one problem.** The same mechanism produces the number
   in several lanes of one problem. These are real and mostly already
   canon: the entropy family of `3n±1` (section 1), the golden family of
   AMM 12592 (section 4), the rational ladder of LRC (section 5).
2. **Transported by a typed bridge that already exists in canon.** Few, and
   each has a proved map: Paley `T_7` is the normalized core of the skew
   Sylvester doubling (THM-447/448), so `189 = h(T_7)` is one number in two
   threads for one reason; the level-11 Hecke density `1/15` and the
   Collatz valuation density `2/15` were typed as an echo, not a map, in
   `level11_short` (section 6).
3. **Small-number coincidences.** Most cross-problem recurrences of
   integers (`139`, `189`, `168`, `1093`, `13`, `507`, `2187`, `1807`) are
   different mechanisms landing on the same small factorization. They are
   listed so nobody rediscovers them as connections (section 7).

The one deeper pattern that survived typing is the first: the Collatz
constants `0.95, 0.05, 0.7925, 0.2075, 1.0527, 0.488` are values, slopes and
reciprocals of one function, and that is now a theorem. Closed forms in
`alpha = log_2 3`: `h* = log_2 alpha - (1 - 1/alpha) log_2(alpha - 1) = 0.949956`,
`lambda* = log_3(1/(alpha - 1)) = 0.48808`, `theta_0 = 1 - alpha/2`.

## 1. Family A (STRUCTURAL, PROVED): the entropy curve of `3n±1`

Let `alpha = log_2 3` and `E(gamma) = h(max(1/2, gamma/alpha))`, `h` the
binary entropy. THM-4487 proves that `E(gamma)` is the exponent of
`#{n <= X : T^i(n) >= n^gamma, i <= log_2 n}` for `gamma in (log_2(3/2), 1]`,
on both sheets. The constants:

| constant | value | role | threads/lanes where it recurs |
|---|---|---|---|
| `E(1) = h(log_3 2)` | `0.949956` | Terras undecided density exponent; exceptional 2-adic dimension; thin-divergence exponent | choice ladder, THM-4476, THM-4487 |
| `1 - E(1)` | `0.050044` | sharp price of bounded-lookahead provability in the pairing family, for arbitrary edits, in the strategy cube; Chernoff rate | THM-4475, THM-4478, THM-4479, THM-4480 (`q = 3` column), synthesis 2k |
| `log_4 3 = alpha/2` | `0.792481` | Korec's exponent: `E = 1` there | barrier atlas, THM-4487 |
| `theta_0 = 1 - alpha/2` | `0.207519` | mean drift per `T`-step in bits; largest dip the counting sees | THM-4476, THM-4487 |
| `-E'(1) = -h'(log_3 2)/alpha` | `0.48808` | exactly the Chernoff tilt `lambda* = log_3(1/log_2(3/2))` of the price exponent (tilt = slope of the rate) | THM-4476 note section 1.8, THM-4480, THM-4487 |
| `1/E(1)` | `1.052681` | growth exponent below which no divergent orbit exists | THM-4476 Cor. 3 |
| `alpha - 1 = log_2(3/2)` | `0.584963` | carry exponent `X^0.585`; lower limit of the dip spectrum; `(1-rho_0)/rho_0` | THM-4476, THM-4478, THM-4487 |
| `1 - H(log_q 2)` | `0.013911 (q=5)`, `0.060510 (q=7)`, `0.100619 (q=9)` | the same curve for other multipliers (edit price); the dip spectrum is trivial for `q >= 5` because `gamma/log_2 q < 1/2` | THM-4480, THM-4487 remark |

**Map.** Terras bijection `Z/2^k -> {0,1}^k` plus the Chernoff/entropy
count of words with a prescribed odd density; the lower bounds use the
cycle-lemma rotation (THM-4478 section 4). **Preserved predicate:** the
count of residue classes with a given odd density. **Lost:** everything
archimedean beyond one window (the transversality barrier), the sheet, and
defects. **Test that separates structure from numerology here:** the same
number must appear as a value of `E` with the same `gamma`, or as `E'`, at
the point where the mechanism is a Terras count.

## 2. Family B (STRUCTURAL inside Collatz, NUMEROLOGY outside): Pillai clocks `2^K - 3^L`

| clock `K/L` | `2^K - 3^L` | Collatz role |
|---|---|---|
| `1/1` | `-1` | free cycle `{-1}` (Gersonides, THM-4484) |
| `2/1` | `+1` | free cycle `{1,2}` |
| `3/2` | `-1` | free cycle `{-5,-7,-10}`; minus-sheet expansion factor `9/8` |
| `8/5` | `+13` | first positive convergent gap; `13 = (3^3-1)/2` as well |
| `11/7` | `-139` | the sporadic cycle `{-17, ...}`: shape `(11,7)`, one integral necklace of 30; expansion `2187/2048 = 1.068` |
| `7/4` | `128 - 81 = 47` | E-SCC budget `c* = ln(128/81)/4 = 0.1144` per digit |
| `19/12` | `-7153` | Pillai clock of the plus crossings (`27` orbit), `7153 = 23 * 311` |

**Where the same integers recur elsewhere, and why it is not a link.**
`139`: LRC's `139/154 = 1 - 15/154` (seven-wall Hunter floor); Pell's
`40449 = 3*97*139`; Sylvester's `1807 = 13*139` (bernoulli lane); AMM's
`R - 3 - d_0 = 139` at `R = 512`. Four mechanisms, one small prime. `13`:
LRC(13)'s speed count against `2^8 - 3^5`. `2187`, `2048`: pure powers used
as sizes in LRC, tournament and JC tables. Verdict: STRUCTURAL inside
Collatz (the continued fraction of `log_2 3` and the cycle equation),
NUMEROLOGY across problems. No map.

## 3. Family C (an identity, no dynamics): base-3 repunits and the Wieferich prime `1093`

`(3^k - 1)/2 = 1, 4, 13, 40, 121, 364, 1093, 3280, ...` are the hubs of the
Q1-mirror climb (loops-and-escapes lane) and the multiplication points of
the trunk. `1093` is also the first Wieferich prime, the canonical hostile
of the `p`-adic row lanes (`2^(p-1) = 1 mod p^2`; all powers of two in one
row). The two Collatz roles are linked only by the identity
`1093 = 1111111_3`; Wieferich-ness is a property of `2` modulo `1093^2` and
has no known relation to `3^7` being close to `2^11`. Outside Collatz,
`1093` appears in LRC and tournament prime lists `547, 911, 1093, 2003, 2549`,
which are primes `= 1 mod 13`; that `1093 = 1 mod 13` is forced by
`3^3 = 1 mod 13` and `1093 = (3^7-1)/2`, a one-line identity with no
dynamical content. Verdict: TRANSPORTED by an identity; not a bridge.

## 4. Family D (STRUCTURAL inside AMM 12592; one untested lead): the golden family

`gamma* = log_5(phi^2) = 0.597987`, `C_* = 1 + gamma* = 1.597987`,
binding fraction `1/phi^2`, `delta = 1/phi`, base `5 = disc(phi)`; the
golden constant is optimal for separately balanced dyadic blocks (THM-3009)
and beaten uniformly: `1.377 <= C* <= 159/100 = 1.59` (THM-4467/4468),
realizable states stalling near `1.567`. **Lead, tested and dead (THM-4494, same day):** `log_2 3 = 1.58496` lay
inside the window `[1.377, 1.59]` only because THM-4468's constant was stuck
at a crude majorant threshold `203/128 = 1.5859`; the exact binomial ratio
in the same certificate gives `C* <= 197/125 = 1.576 < log_2 3`. The entry moves to
section 7. The repo's existing AMM↔Collatz bridge (Bernstein capacity,
THM-3002/3027 and the procgen bridges note) does not predict `log_2 3`.

## 5. Family E (STRUCTURAL inside LRC): the rational ladder, and the `7` cluster

`1/14, 3/41, 2/27, 3/40, 4/53, 14/183, 183 = Phi_6(14), 1/189 = (2/27 - 4/63)/2`
are exact attained or threshold values of the LRC(14) machinery
([CONSTANTS-INDEX](../../00-navigation/CONSTANTS-INDEX.md)). The integer
`189 = 27 * 7` recurs as `h(T_7) = 189` (maximum Hamiltonian paths of the
Paley tournament on 7 vertices; PROVED isomorphic to the skew Sylvester core,
THM-447/448) and as a Collatz fragile-pair index (`{377, 378}` on the orbit
of `27`, THM-4470/HYP-9135: `378 = 14 * 27`). Three factorizations `7 * 27`,
three mechanisms. The `7` cluster (Paley `T_7`, Fano plane, `|PSL(2,7)| = 168`,
LRC(14) `= 2 * 7`, the seven-wall, the sporadic cycle's seven odd steps,
`3^7 = 2187`) is a cluster of a small prime, not a structure: inside
tournaments/Hadamard the `T_7`–Sylvester link is a theorem, everywhere else
`7` is doing a different job. Verdict: STRUCTURAL inside LRC and inside
tournaments; NUMEROLOGY across. (The owner's Hamiltonian-path decomposition of tournaments was tested the same day: the insertion-slot count `1 + #{v -> x_i, x_(i+1) -> v}` is exact, the recursion it suggests is not; see the [insertion-slot note](tournament_insertion_slots_20260926.md).)

## 6. Family F: density constants that look alike

`1/15 = sum 2^(-4j)` and `2/15 = 2^(-3)/(1 - 2^(-4))` are the densities of
odd `n` with `v_2(3n+1) = 0 (mod 4)` and `= 3 (mod 4)`: geometric series in
`1/16`, structural within Collatz. The level-11 Hecke eigenform's `1/15` is a
different object; the repo typed the pair as an echo (`level11_short`). The
fair-coin entropy `8/pi^2 = 0.8106` (bits `0.70028`) belongs to the
Bernstein/AMM lane; it does not recur in Collatz. `0.1144 = ln(128/81)/4` is
Family B.

## 7. Traps recorded

| pair | values | why unrelated |
|---|---|---|
| THM-4475's lower-bound exponent vs the slope constant of Family A | `0.7737 = 0.0500 + log_2((3+sqrt13)/4)` vs `0.7736 = |log_2 log_2(3/2)|` | one involves `sqrt 13` from a two-term recurrence, the other the derivative of the entropy at `log_3 2`; they differ in the fourth decimal |
| AMM golden constant vs the Collatz drift | `1.59799` vs `1.58496` | different problems; the AMM constant is `1 + 2 log_5 phi` |
| `log_2 3` in the AMM 12592 window `[1.377, 1.59]` | `1.58496` vs the majorant threshold `203/128 = 1.5859` | the window's upper end was a proof artifact; `C* <= 197/125 = 1.576` (THM-4494) |
| `189` thrice | `h(T_7)`, `1/189` LRC reserve, fragile pair `189` | three factorizations of `7 * 27` |
| `139` four times | Collatz gap, `139/154`, `3*97*139`, `13*139 = 1807` | one small prime |
| `13` | LRC(13) vs `2^8 - 3^5` vs `(3^3-1)/2` | a small integer |
| `1093` | Wieferich vs `(3^7-1)/2` vs primes `1 mod 13` | an identity chain, no dynamics |
| `507` | `512 - 5` (THM-4475's rescue class) vs `3 * 13^2` (LRC, tournaments) | small integer |

## 8. Method note (for META-PATTERNS, if it earns a card)

The census is cheap (ninety seconds) and its top of the list is worthless
without contexts: the integers `14, 13, 24, 41` head every ranking because
they are sizes. The productive procedure was: rank by distinct threads,
read one sentence per thread, and demand a map. Only recurrences with a map
(a bijection, an isomorphism, an identity) were kept as structure; the one
family whose map is a theorem (the Terras count) turned into THM-4487.
Counterindication: a number that recurs *inside* one problem across many
lanes is usually structural and usually already known there; a number that
recurs *across* problems is usually a factorization accident. The repo's
own rule (state source, target, map, preserved predicate, loss, sidecar,
test) is the right filter, and this note applied it to 291 candidates.

## 9. Reproduction

```bash
python3 04-computation/experiments/constants_atlas_20260926_mine.py . 4 > /tmp/atlas_mine.out
PYTHONIOENCODING=utf-8 python3 04-computation/experiments/constants_atlas_20260926_context.py . i:139 i:189 i:1093 i:168 i:507 i:2187
```

The census depends on the file set at commit time; it is a snapshot, not a
frozen artifact, and no hash is claimed for it.
