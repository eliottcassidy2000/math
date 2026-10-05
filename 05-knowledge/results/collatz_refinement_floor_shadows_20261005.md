# Refinement floors are stopping-time bounds; descent lower bounds live on the 3-adic shadows of negative rational cycle points

2026-10-05, opus session `opus-2026-10-05-S6` (refinement-floor). Owner's seed:
"work toward a source-specific lower bound that survives refinement; pull from
recent and past repository work across topics not yet touched; search for
connections."

**Status.** PROVED (elementary): Theorem T1 (tower conjugation along an orbit
prefix), Lemma T2 (every refinement floor is a stopping-time bound), Proposition
T3 (lower bounds on the tail cost of the Codex atom criterion), Theorem T4 (the
shadow theorem: backward-descent cones are 3-adic shadows of negative rational
cycle points, forward-rising cylinders are their 2-adic shadows), Lemma T4b
(existence of primitive cones at depth `l` iff `{l log_2 3} < log_2(3/2)`,
necessity proved, sufficiency finite-exact), Proposition T5 (exact defect ledger
of every segment closure). FINITE-EXACT: all listed checks, the census of the
descent set to `2^24`, and the sampled forward/backward depths. CITED and typed:
the cross-thread inheritance of section 6. **OPEN:** a certificate-free positive
lower bound at any unresolved source, and universal Collatz. Not a canon
promotion.

## 0. Answer first

A lower bound on the Codex weight at a fixed source that survives refinement of
the residue tower is, exactly, an a priori bound on that source's odd stopping
time (T2). Nothing in the price layer, the finite representations, or any other
repository thread supplies such a bound without the orbit; every thread that
tried (Mahler `3/2`, Hensel, Rule 30, Sun's modular solubility, the LRC dyadic
tower) records the same failure boundary: positive mass in every finite class,
and the atom can still vanish (section 6). What this note adds is where lower
bounds propagate from **smaller** sources and where they cannot:

- They propagate exactly along rising words, and the set of targets reached is
  the union of the 3-adic shadows of the negative rational cycle points
  `x_w = c_w/(2^A - 3^l)` of rising words `w` (T4). Its density among odd
  integers is `0.4687` (census to `2^24`, agreeing with the cone series through
  depth 12, `0.4665`). The complement, the basin minima such as
  `7, 19, 25, 37, 43, 55, ...`, has no smaller ancestor at all: at those sources
  every lower bound is a certificate of the source itself.
- The 2-adic shadows of the same points are the cylinders where forward descent
  is delayed (Terras). Forward and backward induction are organized by one set
  of rational points; combining them halves the unpaid fraction at every depth
  but never empties it (S5: at depth 12, `5.6%` forward-unpaid, `3.1%` jointly).
- The Codex refinement criterion `p_m >= L - R(q)` is quantitatively hopeless
  for the mixtures: the certified two-step family alone forces
  `R(q) >= c/(log q)^2`, so certifying the atom at 27 (odd stopping time 41)
  would need a modulus near `10^90407` (T3). Fixed prices need `q` polynomial
  in `1/atom`, i.e. polynomial in the source, against a logarithmic orbit cost.
- Every "closure" of a source prior along a segment rule pays exactly the
  source atom minus the exit mass landing at each target (T5). Codex's C4 is
  the single-rise rule (first violation 17 -> 13, ratio 21/4) and C5 the
  refuel-block rule; the stopping rule is new and fails first at the trunk entry
  5 (ratio 1.845). The only rule without exits is the full orbit, whose mass is
  the `nu`-mean odd stopping time.

## 1. Inheritance and board

- Closest proved mechanism: Codex's exact criterion, [representation
  positivity (8)](collatz_representation_positivity_20261005.md) and [Fourier
  atom positivity (8)-(10)](collatz_fourier_atom_positivity_20261005.md):
  `w(n) > 0` iff the residue masses `mu_{5^j}(n)` stay above a positive
  `epsilon` for all large `j`, with `0 <= mu_m(n) - w(n) <= sum_{j>=m} w(j)`.
  Also the single-rise ascent closure C4 and the exact boundary debt C5 of
  [three bits and critical flow](collatz_three_bits_critical_flow_20261005.md),
  whose backward R-chain is the depth-`v_3(m+1)` shadow of `-1`.
- Canonical hostiles: Codex's `q(3)=0` measure (every residue mass positive,
  one atom zero) and the repository's MISTAKE-348 (a uniform floor was given
  to pointwise-good supports although every finite extension can be nonzero
  while the limit vanishes).
- Corrected near miss: reading the stabilization of a fixed source's residue
  (free, as THM-4072 notes for Mahler) as a floor on its mass.
- Least-used sidecar: the ordinary size of an ancestor, i.e. whether a lower
  bound comes from below or from above the source.
- Board: **refinement tower / prefix price / rising word / rational cycle point
  / basin minimum / exit mass / tail of the injection measure.**

Previous session (`opus-2026-10-05-S5`,
[Pascal boundary](collatz_pascal_boundary_leaf_section_20261005.md)): Theorem
A (counter-only flows are price mixtures), Lemma B (no power-law comparison,
witnessed by the `-1` shadows), Proposition C (Riesz decomposition), Proposition
D (leaf-section tower). This note is its continuation on the question of
lower bounds. Concurrent commit `905023c4d` (Codex, source-local floors)
audited that note and corrected four scope overclaims in place, logged in
`01-canon/MISTAKES.md` (2026-10-05, Pascal source-floor audit): the discounted
`D` obeys a weighted split and is not a Pascal trace; the lookahead theorem
excludes specified observers only; a positive potential needs no positive
injection at every vertex; summability is the integral condition
`int dmu/(1-r) < infinity`. All four are accepted. Nothing below depends on
the corrected phrasings: T2 uses only `w(L,K) <= w(L,0)`, T4 and T5 use only
the incoming inequality, and the same commit's
[localized resolvent floor](collatz_localized_resolvent_floor_20261005.md)
(monotone signed kernel minorants converging to an atom, with a conditional
finite ROOT deadline) is consistent with T2: its floors recover the atom
effectively but, by its own status line, certify no new source.

## 2. T1: the 2-adic tower of a source is the 3-adic tower of its images

Let `n` be odd with ROOT word `(a_1, ..., a_tau)`, and fix `t < tau`, `A =
a_1 + ... + a_t`, `K_pre = sum_{i<=t} floor((a_i-1)/2)`. For every `s >= 0`,
the integer `n' = n + 2^(A+1) s` has the same first `t` valuations and

    U^t(n') = U^t(n) + 2 * 3^t s.

Hence, for every fixed price `0<r<1`, whenever `U^t(n') != 1`,

    f_r(n') = (1-r)^t r^(K_pre) f_r(U^t(n')),

and the mass of the 2-adic cylinder of `n` at depth `A+1` equals the prefix
price times the fixed-price mass of the odd 3-adic class of `U^t(n)` modulo
`3^t`. Mixtures follow by integrating in `r`. (PROVED: the perturbation
`2^(A+1) s` keeps every valuation because at the last step the odd parts differ
by `2 * 3^t s`; the counters add because the word concatenates. FINITE-EXACT on
624 triples in S1.)

Reading. Refining a source 2-adically pushes the floor question forward along
its orbit and converts it into a 3-adic residue question at the image, with
one ternary digit per odd step; this is the reverse of the ternary-digit
consumption of the inverse-ray notes. The hazard-product form of the
criterion (THM-4263's `prod (1-p_j)`) becomes: the atom survives iff the orbit
reaches the root before the prefix price has decayed, because the 3-adic class
masses are bounded by the total mass.

## 3. T2: a refinement floor is a stopping-time bound

For the beta(1,2) weight, `w(L,K) <= w(L,0) = 2/(L+2)` because `w(L,K+1)/w(L,K)
= (K+1)/(L+K+3) < 1`; for beta(2,2), `6/((L+2)(L+3))`; for a fixed price,
`(1-r)^L`. So if any argument produces `W(n) >= epsilon`, it has produced
`tau(n) <= 2/epsilon - 1`. Conversely the atom itself is the floor. The
independent sweep of today's thread reached the same conclusion from the Green
solver's effective intervals: a certificate-free lower bound is an a priori
bound on the counters. There is no third thing.

## 4. T3: the tail cost of the Codex atom criterion

Let `p_m = lambda(6m+3)` and `R(q) = sum_{j>=q} p_j`. The certified two-step
family `z_k = (2 S^k(1) - 1)/3`, `k = 1 mod 9`, has `lambda(z_k) =
4/((k+1)(k+2)(k+3))` (counters `(1,k)`; checked for `k <= 28`) and `z_k ~
4^k`, so

    R(q) >= sum_{k >= log_4 q, k = 1 mod 9} 4/((k+1)(k+2)(k+3)) ~ 2/(9 (log_4 q)^2).

At `q = 10^10, 10^50, 10^200` this gives `7.6e-4, 2.8e-5, 2.0e-6`. The atom
at the leaf 27 is `lambda(27) = 1/11274451650 = 8.87e-11`, so `L - R > 0` is
impossible unless `q` exceeds about `4^(1.5e5) = 10^90407`. For a fixed price
`r` the root ray alone gives `R(q) >= c q^(-log_4(1/r))`, so the needed
modulus is at least `f_r(27)^(-1/log_4(1/r))`: `7.9e28` at `r=1/2`, `6.5e9`
at `r=1/4`, `2.4e5` at `r=1/16`. The orbit of 27 has 41 odd steps. (PROVED
lower bounds; the two-step family and the root ray are certified, so no
convergence assumption enters.)

Mechanism. The sibling rays decay polynomially in the sibling index, which is
logarithmic in size; the atoms decay exponentially in the stopping time. Any
refinement tower indexed by size therefore pays at least exponentially in the
inverse atom with heavy-tailed priors and polynomially in the source with
geometric prices, whereas the certificate costs `tau = O(log n)` steps.

## 5. T4: the shadow theorem

Call a valuation word `w = (a_1, ..., a_l)` **rising** if `2^A < 3^l`, `A =
sum a_i`. Put `c_w = sum_i 3^(l-i) 2^(a_1 + ... + a_(i-1))`, so that following
`w` forward from `m` gives `n = (3^l m + c_w)/2^A`, and let

    x_w = c_w/(2^A - 3^l)   (a negative rational, a 2-adic and a 3-adic integer).

**Theorem T4 (PROVED; 274,888 exact checks for all 68,722 rising words of
length `<= 12`).** The forward map along `w` is a bijection

    { m odd : m = x_w mod 2^(A+1) }  ->  { n odd : n = x_w mod 3^l },
    m = m_0 + 2^(A+1) s  ->  n = n_0 + 2 * 3^l s,

every `m` on the left follows `w` with exact valuations and rises to `n > m`,
and the backward chain along `w` from every `n` on the right is integral,
positive, and returns to `m`. Consequently, for every positive incoming
supersolution `v`,

    v(n) >= v(m)   with   m < n,   for all n = x_w mod 3^l.

*Proof.* The first `A+1` parity bits of the shortcut map determine and are
determined by `m mod 2^(A+1)` (Terras); the periodic point of the word has the
periodic parity vector, so the cylinder of `x_w` at depth `A+1` is the set of
sources that follow `w` with the last valuation exact. Then `n - x_w = 3^l (m -
x_w) / 2^A` after using `3^l x_w + c_w = 2^A x_w`, which gives the residue
`x_w mod 3^l` and the bijection. Backward integrality is the parity condition
at each step, which the residue class encodes. Monotonicity of `v` along
orbits is the incoming inequality. QED.

Examples: `(1) -> -1`, `(1,2) -> -5`, `(2,1) -> -7`, `(1,1,1,2,1,1,4) -> -17`
(the three integer negative cycles), `(1,1,2) -> -19/11` (the first rational
one; class `13 mod 27`, chain `7 -> 11 -> 17 -> 13`). Codex's C4 chain
`p_j = 2^j (m+1)/3^j - 1`, `j <= v_3(m+1)`, is the word `(1)^j` with
`x = -1`. Lemma B of the previous note is the family `n = 2 * 3^l - 1` on the
`-1` shadow.

**The descent set.** Let `D` be the set of odd `n` having a smaller odd
ancestor, i.e. some odd `m < n` with `n` on the orbit of `m`. Then `D` is
exactly the union over rising words of the classes `x_w mod 3^l` (up to the
finitely many small `n` where the chain would pass through zero), because a
backward chain descends below `n` iff its word is rising, up to the carries.
The census to `2^24` by the stopping-segment sieve gives

| quantity | value |
|---|---:|
| density of `D` among odd integers | 0.46867 (dyadic blocks `2^18..2^24`: 0.4685 to 0.4687) |
| among units (`3 does not divide n`) | 0.70300 |
| among `n = 1 mod 3` | 0.40600 |
| among `n = 2 mod 3` | 1 (the predecessor `(2n-1)/3`) |
| cone series through depth `l = 1..12` | 0.3333, 0.4444, 0.4444, 0.4568, 0.4568, 0.4595, 0.4632, 0.4632, 0.4646, 0.4646, 0.4653, 0.4665 |
| new (primitive) cone classes per depth | 1, 1, 0, 1, 0, 2, 8, 0, 28, 0, 124, 602 |

Leaves (`3 | n`) have no predecessors and are never in `D`. The first unit
basin minima are `1, 7, 19, 25, 37, 43, 55, 73, 79, 97, 109, 115, 127, 133,
145, ...`; all are `1` or `7 mod 9`, never `4 mod 9` (the `-5` shadow).

**Lemma T4b (new cones exist at depth `l` iff `{l log_2 3} < log_2(3/2)`).** A
primitive rising word (no rising proper suffix) must begin with `a_1 = 1`,
because the suffix of length `l-1` is not rising while the whole word is, and
`2^(a_1) < 3^l/3^(l-1) = 3`. Then `1 + ceil((l-1) log_2 3) <= floor(l log_2 3)`
is necessary, which is `{l log_2 3} < log_2(3/2) = 0.585`. This holds exactly
at `l = 1, 2, 4, 6, 7, 9, 11, 12` and fails at `3, 5, 8, 10` below 12, matching
the count sequence above (necessity PROVED; sufficiency FINITE-EXACT through
12). The depth `l = 7` is the `-17` cycle level (`3^7/2^11`), and `l = 12` is
the convergent `19/12` of `log_2 3`, where 602 new classes appear. The depths
form a Beatty-type set; this is the "tight spine block" structure of the
Collatz-thread notes of 2026-09-30 reappearing in the inverse direction, and it
is recorded here as an observation, not a theorem about all depths.

**Forward and backward together (S5, 20,000 random odd sources below `2^40`).**
Mean forward odd-step stopping depth 3.57. Backward descent within depth 12
exists for 46.6% (33.5% at depth 1, 11.1% at depth 2). The unpaid fractions:

| depth | forward-unpaid | forward-and-backward-unpaid |
|---:|---:|---:|
| 3 | 0.2562 | 0.1374 |
| 5 | 0.1523 | 0.0824 |
| 8 | 0.0917 | 0.0500 |
| 12 | 0.0559 | 0.0311 |

Backward descent roughly halves the forward-unpaid set at every depth and, by
T4, can never do better than the density of `D`; the basin minima are paid
only forward, i.e. only by their own orbit.

## 6. Cross-thread inheritance (typed)

Three sweeps were run over threads this thread had not touched. Each row gives
the mechanism, what it certifies pointwise, and the typed transfer.

| Thread / result | What it proves | Pointwise or density | Transfer (map / preserved / lost) | Status |
|---|---|---|---|---|
| Mahler `3/2`: THM-2228, THM-2352, THM-3848, THM-4072 | carry words biject with residues mod `2^m`; `A = r_m + 2^m k` gives `T^m(A) = u_m + 3^m k`; every finite reset skeleton is one residue cone; every finite local terminal-prefix test is vacuous; `A=8` vs `A=13` share state, clock and flag yet decide differently | pointwise characterization, Haar-null safe set | T1 is the same affine law (`U^t(n + 2^(A+1)s) = U^t(n) + 2*3^t s`); preserved: the 2-adic point and the cones; lost: the real tail; the open step is the Collatz step | CITED, same shape |
| Hensel: THM-3446, THM-3449, THM-3452 | free actions at every depth once first carries are independent; orbit sizes `p^(sum(a - c_i))`; Heisenberg law through `c_1 + c_2 + min` | pointwise freeness | anti-atomic: an invariant weight spreads as `c p^(-a)`; atoms need stabilizers (the dependent converse) | CITED, prunes routes |
| Rule 30: THM-4204, THM-4263 | reset words give an absorbing rank-one state stable under every extension; survivor mass is a hazard product `prod(1 - p_j)` | THM-4204 pointwise, THM-4263 density | a reset word is a ROOT certificate; the hazard product is criterion (8) rewritten as `sum p_j < infinity`, and by T1 the hazards are ratios of 3-adic class masses of consecutive images | CITED |
| LRC dyadic tower: THM-2073, THM-2075, THM-2077 | unique safe child per lift by capacity saturation; mass halves per level, harmless only because depth `<= 4` | pointwise | Collatz has no exact count forcing a unique heavy child; mass loss per level is the sibling tail | CITED, no transfer |
| Sun's modular solubility: THM-4027 and its hostile THM-4026 | every residue of every modulus represented, Hensel-stable density; one integer still has no preimage: "the missing coordinate is archimedean alignment" | density | identical to Codex's `q(3)=0`; the archimedean coordinate here is ordinary size, which T4 uses | CITED |
| AMM 12592: THM-3340, THM-3342 | pointwise floors per dyadic horizon with the horizon depending on the target; no uniform sublinear extractor | pointwise per target | the quantifier order of criterion (8): `j_0` may depend on `n` | CITED |
| Hadamard shells: THM-3417 | exchangeable laws see only Krawtchouk shells; the H8 packet is blind to every exchangeable law yet labelled responses are nonzero | density | Theorem A of the previous note is the Collatz instance: counter-only (shell) weights miss the word | CITED |
| META-PATTERNS cards: "separate local support from bounded-height coverage", "existence is a maximum or tail question", "audit the positive kernel before lifting through a flat square", "the same representation is not the same carrier"; MISTAKE-343, MISTAKE-348; reflection "the witness is pointwise, the average is blind" | the repository's own rule that measure never certifies an atom | — | this note obeys them: every lower bound is attached to an actual smaller source (T4) or declared equivalent to the orbit (T2) | policy |

The three transfers the sweep ranked highest (capacity saturation plus hazard
product; the Mahler fibre product; unit first-carry stability) all either
reformulate criterion (8) or stabilize normalized densities. None yields an
atom. That is the inheritance result: the obstruction is the same in every
thread, and it has one name in each.

## 7. T5: the exact defect ledger of segment closures

Let `nu` be any positive summable source prior (the checks use Codex's
`nu(m) = 8/(3 * 4^bitlength(m))`). A **segment rule** assigns to each odd start
`m` a finite initial orbit segment `seg(m) = {m, U(m), ..., U^k(m)}` and the
exit `e(m) = U^(k+1)(m)`. Define the closure `W_sigma = sum_m nu(m) 1_seg(m)`.

**Proposition T5 (PROVED).** For every odd target `y != 1`,

    (W_sigma - K W_sigma)(y) = nu(y) - nu{ m : e(m) = y }.

*Proof.* Each pair `(m, p)` with `p in seg(m)` and `U(p) = y` contributes
`nu(m)` to `K W_sigma(y)`; either `y in seg(m)` with `m != y` (it contributes
to `W_sigma(y)` as well) or `y = e(m)`. No nontrivial cycle passes through the
segments below the verified range, so `y = m` is impossible for `m != 1`. QED.

So `W_sigma` is a supersolution iff no target receives more exit mass than its
own atom. Instances:

- single-rise rule (exit at the first valuation `>= 2`): this is Codex's C4
  closure `B a`; the first violation is `y = 13` with exit/atom `21/4` (the
  chain `7 -> 11 -> 17 -> 13`), reproducing their (13); 540 violations among
  the 2,046 odd targets below `4096`;
- refuel-block rule: Codex's C5, `(I-K)W_r = a - C^(r+1) a`, the exit mass
  being the block boundary distribution;
- stopping rule (exit at the first value below the start): new; the first
  violation is the trunk entry `y = 5` with ratio `1.845` (starts
  `7, 11, 13, 15, 23, 35, 53, ...` exit there), then `19, 23, 37, ...`; 509
  violations below `4096`. The identity is checked exactly at all 2,046
  targets with 32,767 starts.

The only rule without exits is the full orbit, whose closure is the Green
weight with mass `E_nu[tau]`, the "stronger sufficient target" of C5. Every
truncation creates violations exactly at its exits; no truncation is free. In
T4's language, the stopping rule's exits are the landing points of the forward
descents, and the single-rise rule's exits are the targets of the first fall
after a `-1`-shadow climb.

## 8. What survives and what is owed

- A source-specific lower bound surviving refinement exists exactly when a
  smaller ancestor exists (T4) or the orbit has been run (T2). On the basin
  minima (29.7% of units) only the second option exists.
- The Codex criterion (10) is correct and, for the current weights,
  astronomically more expensive than the orbit (T3). A refinement scheme that
  could compete must index the tower by stopping time, not by size; T1 says
  that such a tower is the 3-adic class tower of the images.
- Obligations: (i) independent audit of T4's bijection and of T5; (ii) the
  density of `D` as an exact series over primitive cones, with the Beatty
  structure of T4b proved for all depths; (iii) whether the basin minima carry
  a positive share of the injection measure `lambda` (they are units, so they
  are not sources of `lambda`; their atoms are sums over leaves above them);
  (iv) the Codex obligation is unchanged: a source-dependent lower bound on
  `lambda` across a refuel boundary for a named unbounded family of leaves.
  T4 restricts where such a family can be fed from below: only through
  shadows of rising words.

## 9. Reproduction and scope

[Script](../../04-computation/experiments/collatz_refinement_floor_shadows_20261005.py),
[output](collatz_refinement_floor_shadows_20261005.out),
[JSON](collatz_refinement_floor_shadows_20261005.json):

```text
python3 04-computation/experiments/collatz_refinement_floor_shadows_20261005.py --sieve-bits 24 --cone-depth 12 --samples 20000 --json 05-knowledge/results/collatz_refinement_floor_shadows_20261005.json
python3 -O 04-computation/experiments/collatz_refinement_floor_shadows_20261005.py --sieve-bits 20 --cone-depth 8 --samples 2000
```

10,865,520 explicit checks in 3.3 s (numba for the sieve). Universe: seven
sources with prefixes `t <= 8` and shifts `s <= 11` for T1; counters
`L, K <= 60` for T2; the two-step family `k <= 28` and the tails at five
moduli for T3; all rising words of length `<= 12` with four shifts each for T4,
the sieve of all odd starts below `2^24` (maximal excursion `2.0e13`, inside
int64), and 20,000 random sources below `2^40` for S5; 32,767 starts and 2,046
targets for T5. Hostiles: leaves never in `D`; `n = 2 mod 3` always in `D`;
the shadow check with the modulus `2^A` instead of `2^(A+1)` fails at `(1)`
and is recorded as the reason for the extra bit. The sampled table is
VERIFIED-sampled, not exact; the census is exact in its range. No finite number
is treated as an unbounded statement.
