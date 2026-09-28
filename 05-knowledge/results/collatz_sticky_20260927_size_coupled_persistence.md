# STICKY as a barrier: the aliquot driver's persistence is size-coupled (a theorem), Collatz's memorylessness is size-free (exact), stationary memory is not a coordinate, the barrier-atlas typing of the published mechanisms, and the nine-iteration typology re-checked

**Session:** opus, `collatz-poset-dag-20260927` (S16 continuation), 2026-09-27.
**Owner's directive (relaying the parallel session's proposal):** "provability
tracks the presence of an exact invariant, and among the invariant-free maps
the conjectured fate tracks the sign of the drift and the memory of the
driving quantity. Collatz and Juggler share the corner 'negative drift,
memoryless, no invariant'. I propose STICKY as a control for the repo's
barrier atlas: the aliquot map has Collatz's 2-adic engine with sticky
valuations and is conjecturally divergent, so any termination argument that
never uses memorylessness would prove the Catalan–Dickson conjecture. The
residue-averaging mechanisms (Terras, Korec, Tao, Krasikov–Lagarias) pass;
word-Lyapunov functions and local ranks fail."
**Source of the proposal:** the parallel opus session `collatz-posets-zeta5-20260927`,
[ninth note](collatz_aliquot_lehmer_five_20260927.md) (drivers lock in for
free, shadows are priced; transition matrices at fixed range) and
[tenth note](collatz_two_carries_typology_20260927.md) (the two carries;
the free cofactor; the typology; STICKY proposed; directions D28 "enter
STICKY into the atlas" and D29 "does the growth regime strengthen with
size?"). Both notes carry "independent audit OWED".
**Inherits (cited):** the [barrier atlas](collatz_procgen_20260922_barrier_atlas.md)
(typing rule: a result *overcomes* a control when its conclusion is false
there, and is *blind* when its conclusion still holds; controls SHEET,
DRIFT, DEFECT, INTEGRAL, UNIFORM, DIM), Terras 1976 (parity vectors are a
bijection modulo `2^k`), Erdős 1976 and Erdős–Granville–Pomerance–Spiro
1990 (aliquot ratios persist over any fixed horizon for almost all `n`;
recollection, one-sided in 1976 and two-sided in 1990), Davenport 1933
(`sigma(n)/n` has a continuous distribution), Deléglise 1998 (abundant
density in `(0.2474, 0.2480)`), Landau / Hardy–Ramanujan (the density of
integers with a bounded number of prime factors is `O((log log x)^(k-1)/log
x)`), Guy–Selfridge 1975 (drivers; the counter-conjecture), Catalan 1888 /
Dickson 1913, S13's Proposition 4 (no prefix rank), my S15 note
([spine, tight blocks, descent tree](collatz_posets_dags_20260927_spine_descent_tree.md)).
All literature is from memory and flagged as such; nothing below depends
on a literature statement except where marked CITED.

**Status: PROVED elementary (Theorem 1, Proposition 2, Propositions A and
B) + FINITE-EXACT (exact identities on `10^6` even numbers; persistence and
loss rates in seven dyadic ranges from `10^3` to `10^9`; Collatz exactness;
random-model simulations; every number of the parallel typology
re-computed) + CITED (the aliquot literature) + DIRECTION (section 6).
Collatz OPEN; Catalan–Dickson / Guy–Selfridge OPEN. Verdict on the
proposal: STICKY is a real coordinate and belongs in the atlas, but the
coordinate is *size-coupled persistence*, not memory; memory with a
stationary law changes fluctuations only (Proposition A); the aliquot map
has size-coupled persistence provably (Theorem 1), together with negative
average drift like Collatz, which is what makes the control non-vacuous;
and under the atlas's own rule Terras, Korec and Tao overcome STICKY while
Krasikov–Lagarias's predecessor count is blind to it. Independent audit:
section 8.** Scripts `04-computation/experiments/collatz_sticky_20260927.py`
and `collatz_sticky_20260927_extra.py`, outputs beside them.

## 0. The answer in one paragraph

Memory of the driving quantity is not by itself a fate-changing coordinate:
a valuation process with the Collatz marginal law (geometric of ratio
`1/2`) and any stationary dependence still has slope `log_2 3 - 2 = -0.415`
bits per step by the ergodic theorem, so its walk goes to `-infinity`
whatever its stickiness; stickiness only inflates the fluctuations
(simulated maximum heights `1.3, 3.1, 17.7, 169` at persistence `0, 0.5,
0.9, 0.99`; Proposition A). What separates the aliquot map from Collatz is
that its persistence is **size-coupled** and its growth classes carry a
**positive conditional drift**: for `n = 2^a m` the driver `2^a` is lost
exactly when `v_2(sigma(m)) <= a - 1`, which forces `m` to have at most `a -
1` prime powers to an odd exponent, so the loss probability in the class
`a` tends to zero with the size, exactly like `#odd squares` for `a = 1`
(loss `0.0200, 0.0060, 0.00184, 0.00058` in the ranges `[N, 2N)`, `N = 10^3,
..., 10^6`, matching the square count to every printed digit) and like
`1.4/ln N` for `a = 2` (`0.240, 0.170, 0.129, 0.104, 0.091, 0.077, 0.066` from
`10^3` to `10^9`), while the conditional drifts are `-0.35` bits in the class
`a = 1` and `+0.13, +0.32, +0.40, +0.45, +0.47` in the classes `a = 2, ..., 6`,
with a negative average of `-0.048` bits over even `n` (Theorem 1; the
parallel session's Proposition 2 is its algebraic core). Collatz has none of
this: `P(next valuation = 1 | valuation = 1)` is exactly `1/2` in every
interval `[N, 2N)` (it is the proportion of `n ≡ 7 (mod 8)` among `n ≡ 3
(mod 4)`), and the `0.52` measured along orbits by the parallel session is a
small-number effect (`0.4965` on orbit steps with value at least `10^4`,
`0.5334` below). Proposition B gives the dichotomy behind the fates: a
growth state left with probability `p(N)` at size `N` is held for ever with
positive probability iff `sum_k p(N_0 r^k) < infinity`; the aliquot *parity*
lock (`p ~ 0.59 N^(-1/2)`) is summable, so an even sequence stays even for
ever with probability `1 - O(N^(-1/2))`, but the *growth-driver* locks (`p
~ c/ln N`) are not summable, so in the model every growth driver is lost
infinitely often (the sequence of 276 loses `2^2·7` near step 170 and falls
from 28 to 15 digits, then regrows), and divergence, if it happens, is a
competition of phases of length `~ ln N`, not a driver held for ever. The
atlas typing (section 4): under the conclusion rule, Terras/Everett, Korec
and Tao **overcome** STICKY (their aliquot analogues fail under
Guy–Selfridge, and the fixed-horizon part fails provably by Erdős and the
abundant density), Krasikov–Lagarias's predecessor count is **blind** (all
primes reach `1`, so the aliquot analogue holds trivially), the cycle
mechanisms, conjugacies and rewriting encodings are not applicable, the
repo's thin-divergence theorems are conjecturally blind, and the word-count
theorems are not applicable. The transfer reading of the proposal is the
sharper no-go: an argument for Collatz termination whose steps use only the
engine, the marginal law of the valuations and the sign of the average
drift, and are insensitive to whether the persistence of the growth state
is size-coupled, applies verbatim to the aliquot map and proves
Catalan–Dickson; Proposition A shows that sensitivity to stationary memory
is not enough, so a valid argument must use size-free memorylessness, i.e.
Terras's bijection at every scale, for the specific integer, which is the
parallel session's cofactor identity and my spine's tight blocks read
together (section 6). Nothing here is a proof step for Collatz.

## 1. Inheritance: what is being tested, and the two readings of "pass"

The proposal has three parts: a typology of nine iterations, a new control
STICKY with the aliquot map as witness, and a typing of the published
mechanisms against it. The atlas types by *conclusions* (a mechanism
overcomes a control when its conclusion is false on the control), while the
proposal's last sentence types by *mechanism* ("residue averaging passes",
"word-Lyapunov functions fail"). The two readings agree on Terras, Korec
and Tao and disagree on Krasikov–Lagarias; and "fail" for hypothetical
arguments is only meaningful in the transfer reading (a hypothetical
termination argument concludes termination, which is false on the control
under Guy–Selfridge, so under the conclusion rule it would "overcome"). The
note keeps both readings apart and states both.

Closest proved mechanism: Terras's bijection, which is the exact form of
"memoryless at every scale" (Proposition 2). Canonical hostile: the
stationary sticky valuation chain of Proposition A, which is as sticky as
one likes and still converges. Corrected near miss: "memory" as the
coordinate; the coordinate is size-coupling. Least-used sidecar: the parity
of `sigma(p^e)`, which turns the driver lock into a statement about how many
prime factors `m` has to an odd power.

## 2. Size-coupled persistence (PROVED + FINITE-EXACT)

Notation. `s(n) = sigma(n) - n`; `n = 2^a m` with `m` odd and `a >= 1`; the
*driver class* of `n` is `a = v_2(n)`. The parallel note's Proposition 2:
`s(n) = (2^(a+1) - 1) sigma(m) - 2^a m`, hence `v_2(s(n)) = min(a, v_2
sigma(m))` when `v_2 sigma(m) != a`, and `v_2(s(n)) >= a + 1` when they are
equal (re-checked here on every even `n <= 2·10^6`).

**Theorem 1 (the aliquot driver's persistence is size-coupled).**
(a) The driver is *lost* (`v_2(s(n)) < a`) iff `v_2(sigma(m)) <= a - 1`;
*kept exactly* (`v_2(s(n)) = a`) iff `v_2(sigma(m)) > a`; *strengthened*
(`v_2(s(n)) >= a + 1`) iff `v_2(sigma(m)) = a`.
(b) `v_2(sigma(m)) = sum over p^e || m with e odd of (v_2(p + 1) + v_2(e + 1)
- 1)`; even exponents contribute nothing. In particular `v_2(sigma(m)) >=
#{p^e || m : e odd}`.
(c) For `a = 1` the driver is lost iff `m` is a perfect square; among `n in
[N, 2N)` with `v_2(n) = 1` the loss proportion is `#{odd r : N/2 <= r^2 < N}
/ #{odd m in [N/2, N)} = (2 - sqrt 2) N^(-1/2) (1 + o(1)) ≈ 0.586 N^(-1/2)`.
(d) For fixed `a >= 2` the loss proportion in the class `a` tends to zero as
`N -> infinity`: it is `O_a((log log N)^(a-2) / log N)`, and for `a = 2` it
is `≍ 1/log N` (the loss set contains every `m = p r^2` with `p ≡ 1 (mod
4)` prime and `p` not dividing `r`).
(e) Consequently, in every fixed class the persistence `P(v_2 s(n) = a |
v_2 n = a)` tends to `1` as the size grows (the strengthening event `v_2
sigma(m) = a` also has density zero for each fixed `a`), and the
conditional drifts are `E[log_2(s(n)/n) | a] = -0.348, +0.129, +0.318,
+0.404, +0.445, +0.465` for `a = 1, ..., 6` on even `n <= 2·10^6` (the
parallel note's numbers, reproduced), with the unconditional averages
`-0.048` bits over even `n` and `-2.70` bits over all `n <= 2·10^6` (odd `n`
have `s(n) < n` almost always).

*Proof.* (a) is the case split of the valuation identity. (b) For odd `p`
and even `e`, `sigma(p^e)` is a sum of `e + 1` odd terms, hence odd. For odd
`e`, `sigma(p^e) = (p + 1)(1 + p^2 + ... + p^(e-1))`, the second factor a sum
of `(e + 1)/2` odd terms, so `v_2(sigma(p^e)) = v_2(p + 1) + v_2((e + 1)/2)`.
Multiplicativity of `sigma` adds the contributions. (c) By (a),(b) the loss
for `a = 1` means `v_2(sigma(m)) = 0`, i.e. no prime divides `m` to an odd
power, i.e. `m` is a square; the count of odd squares in `[N/2, N)` is
`(sqrt N - sqrt(N/2))/2 + O(1)` against `N/4` odd numbers. (d) Loss forces
`m = q r^2` with `q` squarefree and `omega(q) <= a - 1`. Splitting `r <=
N^(1/4)` (where Landau's bound `#{q <= x : omega(q) = k} << x (log log
x)^(k-1)/log x` applies with `x = N/r^2 >= N^(1/2)`, summed over `k <= a - 1`
and over `r` with `sum 1/r^2 < infinity`) from `r > N^(1/4)` (where the trivial
bound `N/r^2` sums to `O(N^(3/4))`) gives `O_a(N (log log N)^(a-2)/log N)`
members below `N`, against `≍ N` members of the class. For the lower bound
at `a = 2`, `m = p r^2` with `p ≡ 1 (mod 4)`, `p` prime not dividing `r`, has
`v_2(sigma(m)) = v_2(p + 1) = 1`, and there are `≍ N/log N` such `m <= N` (take
`r` odd and at most `N^(1/4)`, and primes `p <= N/r^2` in the class `1 mod
4`). (e) follows from (d) and the identical argument for the event `v_2
sigma(m) = a` (which forces at most `a` odd-exponent primes). The drifts are
computed. ∎

*FINITE-EXACT (script, part 1 and extra).* Ranges `[N, 2N)`; exact by a
divisor sieve for `N <= 10^6`, sampled (20000 factorizations per range) for
`N = 10^7, 10^8, 10^9`:

| `N` | `a=1` keep | `a=1` loss (square prediction) | `a=2` keep | `a=2` loss | `loss·ln N` | `a=3` keep | `a=3` loss | `a=4` keep | `a=4` loss |
|---|---|---|---|---|---|---|---|---|---|
| `10^3` | 0.788 | 0.02000 (0.02000) | 0.576 | 0.240 | 1.66 | 0.365 | 0.444 | 0.129 | 0.742 |
| `10^4` | 0.844 | 0.00600 (0.00600) | 0.662 | 0.170 | 1.57 | 0.466 | 0.352 | 0.240 | 0.588 |
| `10^5` | 0.879 | 0.00184 (0.00184) | 0.719 | 0.129 | 1.49 | 0.527 | 0.296 | 0.333 | 0.493 |
| `10^6` | 0.902 | 0.00058 (0.00058) | 0.759 | 0.104 | 1.44 | 0.583 | 0.251 | 0.401 | 0.432 |
| `10^7` (sampled) | 0.915 | 0.0000 (0.00019) | 0.786 | 0.091 | 1.46 | 0.629 | 0.231 | 0.452 | 0.395 |
| `10^8` (sampled) | 0.925 | 0.00015 (0.00006) | 0.817 | 0.077 | 1.41 | 0.637 | 0.217 | 0.492 | 0.350 |
| `10^9` (sampled) | 0.938 | 0.0000 (0.00002) | 0.832 | 0.066 | 1.36 | 0.676 | 0.183 | 0.529 | 0.308 |

The entry rate into the growth classes from `a = 1`, `P(v_2 s(n) >= 2 | v_2
n = 1) = P(v_2 sigma(m) = 1)`, is `0.192, 0.150, 0.120, 0.098, 0.083, 0.070,
0.062` over the same ranges, i.e. `≈ 1.3/ln N`: transitions in *both*
directions between the decay class and the growth classes become rarer
like `1/ln N`. The parallel session's transition matrix (its ninth note,
even `n <= 2·10^6`) is the `N ≈ 10^6` row of this table; D29 asked whether
the growth regime strengthens with size, and the answer is yes, at the
rates above.

**Proposition 2 (Collatz is memoryless at every scale, exactly).** For odd
`n`, the first two Syracuse valuations are functions of `n mod 8`: `v_1 = 1`
iff `n ≡ 3 (mod 4)`, and then `v_2 = 1` iff `n ≡ 7 (mod 8)`. Hence in every
interval `[N, 2N)` the conditional proportion `P(v_2 = 1 | v_1 = 1)` equals
`#{n ≡ 7 (8)}/#{n ≡ 3 (4)} = 1/2` up to a boundary term `O(1/N)` (exactly
`125/250`, `12500/25000`, `1250000/2500000`, `125000000/250000000` for `N =
10^3, 10^5, 10^7, 10^9`), and Terras's bijection extends this to all finite
words: the valuations of a uniformly random residue are i.i.d. geometric of
ratio `1/2`, independently of the size. Along actual orbits of odd `n <
2·10^5` the persistence is `0.4965` on the `708316` steps whose current
value is at least `10^4` and `0.5334` on the `1313412` steps below `10^4`:
the parallel note's `0.52` is the small-number part.

*Proof.* `(3n + 1)/2 ≡ 3 (mod 4)` iff `3n + 1 ≡ 6 (mod 8)` iff `n ≡ 7 (mod
8)`. ∎

## 3. Memory is not the coordinate; size-coupling is (PROVED)

**Proposition A (stationary memory does not change the fate).** Let
`(v_l)` be a stationary ergodic sequence of positive integers with `E[v] >
log_2 3`. Then the height walk `h_l = l log_2 3 - (v_1 + ... + v_l)` satisfies
`h_l / l -> log_2 3 - E[v] < 0` almost surely, whatever the dependence
structure. In particular every stationary "sticky" version of the Collatz
random model (same geometric marginal, any persistence) has the same
negative drift `-0.415` bits per step and its walk tends to `-infinity`;
stickiness changes the fluctuations (the expected maximum height, the
entropy of the no-descent words) and not the fate.

*Proof.* The ergodic theorem for the stationary sequence `log_2 3 - v_l`. ∎

*Simulation (script, part 3).* Valuation chains that keep the previous
value with probability `p` and otherwise draw a fresh geometric: mean slope
`-0.413, -0.420, -0.427, -0.423` bits per step for `p = 0, 0.5, 0.9, 0.99`
(all within noise of `-0.415`), mean maximum height over 4000 steps `1.3,
3.1, 17.7, 169`, and paths above their start after 4000 steps `0, 0, 0, 16`
out of 200: at `p = 0.99` the fluctuations are of order `1/(1 - p)` and a
few paths have not yet felt the drift, but the slope is unchanged.

**Proposition B (size-coupled escape).** Suppose a growth state multiplies
the iterate by `r > 1` per step and is left, independently at each step,
with probability `p(N)` when the iterate has size `N`. Starting from `N_0`,
the probability of never leaving the state is `prod_k (1 - p(N_0 r^k))`,
which is positive iff `sum_k p(N_0 r^k) < infinity`. For `p(N) = c N^(-1/2)`
the sum is geometric and the state is held for ever with probability `1 -
O(N_0^(-1/2))`; for `p(N) = c / log N` the sum is harmonic and the state is
left almost surely, and, over renewed visits, infinitely often.

*Proof.* Independence and the product; `sum_k c (N_0 r^k)^(-1/2) < infinity`;
`sum_k c/(log N_0 + k log r) = infinity`. ∎

*Simulation (script, part 3).* From `N_0 = 2^20` with `r = 2^0.3`: the state
with escape `3 N^(-1/2)` is never left in 20000 steps in `97.4%` of 2000
trials; the state with escape `0.6/ln N` is always left.

*Reading for the aliquot map (with Theorem 1).* The parity lock (`a >= 1`)
is the summable case: an even number of size `N` stays even for ever with
model probability `1 - O(N^(-1/2))`, which is the "for free" of the parallel
session's ninth note, and it is exact only for parity. The growth-driver
locks (`a >= 2`) are the harmonic case: lost infinitely often in the model,
and re-entered from `a = 1` at the same rate `≈ 1.3/ln N`; phases have
expected length `≍ ln N`, the growth phases have conditional drift `+0.13`
to `+0.47` bits per step and the decay phase `-0.35`, and the fate is a
competition between them on the logarithmic scale that this model does not
decide (it is Guy–Selfridge against Catalan–Dickson). The 276 sequence
(ninth note: driver `2^2·7` to step ~170, loss, fall from 28 to 15 digits
under the driver `2`, recovery under `2^3, 2^7·3^2, ...`) is the model's
typical behaviour, not an exception to "drivers lock in".

*Consequence for the proposal.* "Memory of the driving quantity" is the
wrong name for the coordinate: what the aliquot map has and Collatz lacks
is (i) persistence that grows with the size of the iterate and (ii) a
growth class whose conditional drift is positive while the average drift is
negative. Collatz's growth class (`v = 1`, `+0.585` bits) has size-free
persistence exactly `1/2` and is the unstable one; the aliquot growth
classes (`a >= 2`) have persistence `-> 1` and are the stable ones. The
parallel note's sentence "Collatz grows only in the class that cannot be
held; the aliquot map grows only in the classes that hold themselves" is
correct, and Theorem 1 says at what rate the holding improves.

## 4. STICKY in the barrier atlas (D28 done)

**Definition (control STICKY).** The aliquot map `s(n) = sigma(n) - n` on
`n >= 2`, a map with negative average drift (`-0.048` bits per step on even
`n`, `-2.7` on all `n`), no exact invariant, and size-coupled persistence
of its growth classes (Theorem 1). Failure mode: divergence produced by a
growth state whose persistence increases with the size of the iterate,
under negative average drift; conjectural witnesses: the Lehmer five `276,
552, 564, 660, 966` (Guy–Selfridge). A mechanism *overcomes* STICKY when its
conclusion is false for the aliquot map (provably, or conjecturally under
Guy–Selfridge); it is *blind* when its conclusion still holds there.

**Typing.** Published mechanisms (rows of the atlas; only the divergence
side is typed, the cycle side has no aliquot analogue in this sense):

| mechanism | aliquot analogue of the conclusion | STICKY |
|---|---|---|
| Terras 1976 / Everett 1977 (finite stopping time on a set of density 1) | "for almost all `n` some `s^k(n) < n`": false under Guy–Selfridge; and provably, for every fixed horizon `K`, the density of `{n : s^k(n) < n for some k <= K}` is at most `1 - d(epsilon) + o(1)` with `d(epsilon) = density{sigma(n)/n > 2 + 2 epsilon} > 0` (Erdős 1976: for almost all `n`, `s_(j+1)(n)/s_j(n) > s(n)/n - epsilon` for `j < K`; Davenport: continuity of the distribution of `sigma(n)/n`; Deléglise: `d(0) ≈ 0.2476`) | **overcomes** (fixed horizon provably; full statement conjecturally) |
| Korec 1994 (`T^n(y) < y^c` on density 1) | same | **overcomes** |
| Tao 2019/22 (`Col_min(N) < f(N)` on logarithmic density 1) | "min of the aliquot sequence below `f(n)` for log-almost-all `n`": false under Guy–Selfridge | **overcomes** (conjecturally) |
| Krasikov–Lagarias 2003 (`#{n <= x : n reaches 1} >= x^0.84`) | every prime reaches `1`, so `pi(x) >> x^0.84` of the integers reach `1` | **blind** |
| Kontorovich–Lagarias stochastic models | the models assume i.i.d. valuations; STICKY is the name of that assumption | n/a (model) |
| Steiner, Simons–de Weger, Eliahou, Hercher (cycles) | cycle side | n/a |
| Bernstein–Lagarias, Monks–Yazinski (conjugacies) | no analogue | n/a |
| Yolcu–Aaronson–Heule (rewriting) | no analogue | n/a |
| Hicks–Mullen–Yucas–Zavislak (`F_2[x]`) | a different model with deterministic degree bookkeeping | n/a |
| Barina (verification below `2^71`) | finite | n/a |

Repository results: THM-4476 and THM-4499 (thin divergence: a divergent
orbit has `o(X^(h*))` elements below `X`): the aliquot analogue is expected
to hold (a divergent aliquot sequence should grow geometrically in phases),
so **blind** (conjecturally); THM-4495, THM-4487, THM-4498 (counts of
parity words): n/a; THM-4512 (certified classes): n/a; THM-4514 (spine,
descent tree): the dictionary of its Proposition 1 holds for every integer
map, so n/a; S13's Proposition 4 (no prefix rank): n/a; HYP-9161, HYP-9164:
n/a. Pattern, as for DRIFT: every mechanism that overcomes STICKY is a
residue-averaging mechanism, sign-blind (SHEET) and dimension-blind (DIM);
STICKY adds no mechanism that overcomes anything new; it names what the
averaging mechanisms use (size-free memorylessness) and what a pointwise
argument must reproduce.

**Corrections to the proposal's typing sentence.** (i) Krasikov–Lagarias's
predecessor count is blind under the atlas rule; only Terras, Korec and
Tao overcome. (ii) "Word-Lyapunov functions and local ranks fail" is the
transfer reading: such an argument, if it used only the engine, the marginal
valuation law and the average drift, would transfer to the aliquot map; it
is not a typing of an existing result. (iii) The witness is only
conjecturally divergent, so STICKY, unlike DRIFT (`5x+1`, also only
conjecturally divergent) and like it, is a conjectural control; the parallel
note's D28 asked for a Collatz-like map with provably sticky valuations and
provable divergence, and Theorem 1 supplies half of it: provable
size-coupled stickiness, without provable divergence.

**The transfer no-go, stated.** Let an argument for "every Collatz orbit
descends" use only: the multiplicative engine (`x -> (3x+1)/2^v` and its
aliquot counterpart `2^a m -> (2^(a+1)-1) sigma(m) - 2^a m`), the marginal
law of the driving quantity at each size, and the sign of the average
drift; and let it be insensitive to whether the growth state's persistence
is size-coupled. Then it applies to the aliquot map, whose average drift is
also negative, and proves that every aliquot sequence is bounded (the
Catalan–Dickson conjecture), which is expected to be false. By Proposition
A, sensitivity to stationary memory does not rescue it. So a valid
divergence-half argument must use that the persistence of the Collatz
growth class is size-free and equal to `1/2` at every scale, i.e. Terras's
bijection for the specific integer at every scale; the pointwise form of
this is the parallel session's cofactor identity (after a run of `K` ones
from `m = 2^(K+1) t - 1`, the next valuation is `1 + v_2(3^(K+1) t - 1)`;
re-checked here in 600 cases with the lifting-the-exponent law to `K <
400`) and, on the S15 side, the tight-block structure of the spine
(section 6).

## 5. The nine-iteration typology, re-checked (FINITE-EXACT)

Every number of the parallel session's table was recomputed independently
(script, parts 4–5); agreement throughout, with the size-decomposition of
the Collatz persistence as the one refinement.

| map | drift (recomputed) | memory of the driving quantity (recomputed) | invariant | fate | status |
|---|---|---|---|---|---|
| Collatz `3n+1` | `-0.415` bits per odd step (exact expectation) | memoryless, size-free: persistence exactly `1/2` in every `[N, 2N)`; orbit-measured `0.4965` above `10^4`, `0.5334` below | none | all terminate | OPEN |
| `3n-1` | `-0.415` | memoryless (same bijection) | none | three cycles, all bounded | OPEN |
| `5n+1` | `+0.32` per odd step (`+0.16` per `T`-step) | memoryless | none | most diverge | OPEN |
| aliquot | `-0.048` per step on even `n` (`-2.7` on all); `-0.35` in class 1, `+0.13 .. +0.47` in classes `2..6` | size-coupled: loss `0.586 N^(-1/2)` (class 1), `≈ 1.4/ln N` (class 2); Erdős persistence of increases `1.00, 0.81, 0.72, 0.66, 0.62, 0.56, 0.50, 0.44` over `k = 1..8` (5000 abundant starts `<= 10^6`, lower bounds since chains leaving the `2·10^7` sieve count as failures); abundant density `0.2475` | none | many diverge (Guy–Selfridge) / all bounded (Catalan–Dickson) | OPEN |
| Juggler | `-0.334` in `log_2` of the log-size per step (fair coin `-0.208`), `n <= 2000` | parity persists `0.565` | none | all terminate (all `n <= 2000`; longest 80 steps, 271 digits) | OPEN |
| reverse-and-add | `+0.403` digits per step over the first 50 steps | carries | none | `6091` non-palindromic starts below `10^5` in 300 steps (`249` below `10^4`; first `196, 295, 394, 493, 592, 689, 691, 788, 790, 879, 887, 978`) | OPEN |
| Ducci (length `2^k`) | — | — | linear over `GF(2)`, nilpotent | terminates (`4^8` vectors of length 8, worst 14 steps) | PROVED |
| Kaprekar (4 digits) | — | — | finite state space | `6174` from every non-repdigit, worst 7 steps | PROVED |
| look-and-say | `+`, ratio `1.30361` at 60 steps (Conway `1.30358`) | — | 92 atoms | diverges | PROVED |

The parallel table's Juggler drift `-0.323` (its `n <= 1500`) against
`-0.334` here (`n <= 2000`), and its reverse-and-add growth `+0.409` against
`+0.403`, are sampling differences. Its Collatz memory entry "`0.52`" should
read "exactly `1/2` at every scale" (Proposition 2).

## 6. What this says for a proof (DIRECTION)

* STICKY and my S15 spine meet at one object. A hypothetical divergent
  orbit climbs its spine by tight blocks, each beginning with a run of
  ones; the parallel session's cofactor identity says what ends a run:
  after `K` ones from `m = 2^(K+1) t - 1` the next valuation is `1 + v_2(3^(K+1)
  t - 1)`, decided by the 2-adic distance of the cofactor `t` to
  `3^(-(K+1))`. Terras's bijection is the statement that over `t` uniform
  this distance is geometric at every scale (size-free memorylessness);
  the STICKY no-go says a proof cannot avoid using it. The transversality
  statement of S13/S15 (`Z^+ ∩ E_inf = empty`) is therefore, in these
  coordinates, "the cofactors met along a spine are never, infinitely
  often, 2-adically close to the inverses of the powers of three that the
  tight blocks dictate". That is a restatement, sharpened by the block
  structure (the run lengths are the block prefixes and the block lengths
  are a Beatty sequence), not a method.
* Cheapest next probe with real content: measure, on the `5x+1` control
  (DRIFT) and on random no-descent words (the `E_inf` model), the
  distribution of `v_2(3^(K+1) t - 1)` at the ends of runs, against the
  geometric law; any deviation along actual integer orbits at large size
  would be a size-coupling in Collatz, which Proposition 2 excludes for
  the first two letters and Terras's bijection excludes for any fixed
  number of letters, but which the *orbit-coupled* law (values along one
  orbit, not residues) has never been measured at scale. The S15 chain
  probe found `O(1)` multiplicities on `5x+1`; the same script can measure
  run-end valuations.
* The aliquot side has its own open pointwise problem (the parallel
  note's D25), and Proposition B locates it: prove that a specific sequence
  spends a positive fraction of its logarithmic time in growth classes.
  Nothing transfers from one side to the other beyond the shape.

## 7. Reproduction

```text
python 04-computation/experiments/collatz_sticky_20260927.py        > 05-knowledge/results/collatz_sticky_20260927.out
python 04-computation/experiments/collatz_sticky_20260927_extra.py  > 05-knowledge/results/collatz_sticky_20260927_extra.out
```

Needs numpy (divisor sieves to `2·10^7`) and sympy (`factorint` for the
sampled ranges). Part 1 also asserts the valuation identity, the loss
criterion and the square criterion on every even `n <= 2·10^6`; the extra
script checks the cofactor identity (600 cases), the lifting-the-exponent
law (`K < 400`) and the `n mod 8` exactness (`n < 10^5`). About ten minutes
in total.

## 8. Independent audit

Pending (auditor subagent: blind re-derivation of Theorem 1, Propositions
2, A, B and the atlas typing; own script; the audit also covers the
parallel session's Propositions 1–3 of its tenth note, on which Theorem 1
rests, since those notes carry "audit OWED").
