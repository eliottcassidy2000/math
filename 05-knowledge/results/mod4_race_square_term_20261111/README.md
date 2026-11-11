# The mod-4 race of Hamiltonian path counts: its triangle term, its pentagon term, its square term, and where small numbers stop being special

**Session:** mac-mini-2026-11-11-race. **Canon:** THM-5131 (triangle term, square term, the race tables), THM-5132
(the pentagon theorem). **Hypotheses:** HYP-9663.
**Audit:** [audit/AUDIT.md](audit/AUDIT.md), independent, own programs: PASS WITH CORRECTIONS. No mathematical error;
2 MAJOR findings, both "known result typed as new" (the residue mod 8 of self-converse tournaments; the closed form
of the triangle bias), and 22 minor ones; all applied in this version, nothing added afterwards. MISTAKE-760.

**Owner directive (summary).** Spend another session on the goal "7 and 21 are the only odd numbers that are not
Hamiltonian path counts of tournaments", merging in prior work on the square-sum problem, the Fermat primes,
Chebyshev's bias and other finite behaviour in small numbers that turns into the opposite behaviour in large
numbers, and prior work on dualities such as micro/macrocosm and how they relate to the MRB constant and
`pi^(n/2)/(n/2)!`.

**Outcome.** The conjecture (Chan-Pak, arXiv:2308.10214, Conjecture 6.2) is **not proved here and remains OPEN**.
The repository already explains *why 7 and 21 are small* (THM-1370, THM-4889, THM-4986: the crevasses between the
strong floors) and already has an atlas of exceptional sets (kind-pasteur S15, `root_change_local_factors_20261107/
ATLAS.md`). This session adds the piece that was missing from both: the **race** between the two classes
`H = 1` and `H = 3 (mod 4)`, treated the way the owner's 2026-10-14 prompt treats Chebyshev's race.

1. **(PROVED; the value of the mean was on record) The triangle race has a permanent, exact bias.** For every
   tournament `R` on `m` vertices the `2^m` one-vertex extensions satisfy
   `sum_X (-1)^(c_3(R + v_X)) = 2^ceil(m/2)`, the same number for every `R` (new). Hence
   `E (-1)^(c_3) = 2^(-floor((n-1)/2))` over labelled tournaments of every order `n`: an even number of cyclic
   triangles always leads, by an exactly known margin, with no oscillation (Theorem T). The closed form of the mean
   was recorded in June for the tiling ensemble and verified for `n <= 10` (THM-555; Gauss-sum mechanism in
   HYP-2707); here it is proved for every `n`.
2. **(PROVED) The pentagon race is solved too, and with it the full race up to order 6.** Let `q_1` be the number
   of vertices of odd score and `q_2` the number of vertices with an odd number of odd-score out-neighbours. Then
   `c_3 = C(q_1, 2) + floor(n/4)` and `c_3 + c_5 = C(q_2, 2) + floor(n/4) (mod 2)`, so the parity of the number of
   5-cycles is `C(q_1, 2) + C(q_2, 2)`. Consequences (Theorem P, section 6): every tournament has at least as many
   one-vertex extensions with `c_3 + c_5` even as with `c_3 + c_5` odd, with an explicit count in terms of the kernel
   of the adjacency matrix mod 2; `E (-1)^(c_3 + c_5)` and `E (-1)^(c_5)` are known in closed form at every order
   (for `n = 2 mod 4` they are `2^(-(n-2)/2)` and `2^(-(n-2))`, the triangle value and its square; otherwise about
   twice the triangle value and three times its square); and for `n <= 6`, where there are no longer odd cycles,
   `H mod 4` is decided by `q_2 mod 4`.
3. **(KNOWN + FINITE-EXACT) The full race has a square term.** Self-converse tournaments of even order all have
   `H = 1 (mod 4)`, more precisely `H = 1 + 4s (mod 8)` (KNOWN: THM-791; Session29 and Session130 for the reduction
   to an involution and for the residue mod 8); they are the "prime squares" of this race. There are at most
   `(n-1)!! 2^(n^2/4)` of them among the `2^C(n,2)` labelled tournaments (a corollary of the known count per
   involution), and exactly `48, 7040, 6737920, 31596429312` at `n = 4, 6, 8, 10` (the last two new). New here: the
   weighing. They carry `150%, 86%, 60%, 37%` of the lead `D_n = P(H = 1 mod 4) - P(H = 3 mod 4)` at these orders. At
   `n = 4` the 48 labelled self-converse tournaments (transitive, `H = 1`; strong, `H = 5`) are all on the side
   `H = 1` and the 16 others (a 3-cycle with a source or a sink) all have `H = 3`; at `n = 6` the lead of the classes
   (32 against 24) is 12 self-converse classes against a deficit of 4 among the rest.
4. **(FINITE-EXACT) Three strata.** At every order `4 <= n <= 10` the tournaments with a non-trivial automorphism
   that are not self-converse lean to `H = 3 (mod 4)`; the rigid, non-self-converse tournaments lean to `H = 1`
   from `n = 5` on, and in units of `4^(-n)` their lead is `240, 360, 1417.5, 1181.2, 1965.7, 1653.3` at
   `n = 5, ..., 10`.
5. **(FINITE-EXACT, resolves a recorded question)** HYP-9355 (opus S50) recorded the excess of `H = 1 (mod 8)` as
   present for `n <= 10` and "absent at `n = 11`" (a `10^7`-sample test). The exact values are now known: the
   relative masses of the residues `1, 3, 5, 7 (mod 8)` are `1.00330, 0.99603, 1.00155, 0.99912` at `n = 10` and
   `1.00051, 0.99949, 1.00044, 0.99956` at `n = 11`. The pull is still there at `n = 11`; it is `5 * 10^(-4)`, the
   size of one standard error of that test.
6. **(Pre-registered; one FAILED, three CONFIRMED)** Order 12, by a Monte Carlo that averages exactly over the last
   vertex (14.4 million random bases of order 11): `D_12 = (5.89 +- 0.59) * 10^(-5)`, positive at ten standard
   errors, but `D_12 4^12 = 989 +- 99`, half of the value `2 * 10^3` held for `8 <= n <= 11` (the audit's independent
   estimator, 9 million bases with other seeds: `1039 +- 125`; pooled `1008 +- 77`). The registered plateau (R1)
   FAILED; the rate of decay of the lead is OPEN. The second-order predictions held: the rank law of
   the pair form at order 11 equals the random alternating law to `3 * 10^(-4)` (R2), the deletion-pair
   correlation is `1.017` times the twin term (R3), the corank-1 bias frequency is `0.5040` (R4).
7. **(Four pre-registered predictions about the cycle-length race: two FAILED, two CONFIRMED and then PROVED.)** The
   race resolved by cycle length is exact for `n <= 10` (section 6.3). Two regularities seen at `n = 6, 8` and
   registered for `n = 10` were false (R5, R6). The two registered next (R7: the value `2^(-(n-2)/2)` at
   `n = 14, 18`; R8: equal values at `n = 11, 12` and `15, 16`) held in the sampling run, and looking for the reason
   behind a control of that run gave Theorem P. The parities of `c_3` and of `c_3 + c_5` keep biases of size
   `2^(-n/2)` for ever (that of `c_5` alone is of size `2^(-n)`); the direct continuation of the construction to the
   heptagon fails, and the longest cycles erase the bias.
8. **(Typed merge)** Section 8 places these results, and the inherited ones, in one table of small-number
   phenomena: crevasse type (Chan-Pak, square-sum), Borel-Cantelli type (Fermat primes, Sun's `2-4-6-8`), square-term
   type (Chebyshev, Polya, and the race here), fold type (ball volumes, `n^(1/n)`), permanent-bias type (the
   triangle and pentagon parities).

Nothing here enlarges the range of verified values or proves a new case of the conjecture. What is new and proved
is about the second binary digit of `H`: its exact law on the two shortest cycle lengths.

---

## 1. Inheritance pass and concurrency

Read lane by lane before deriving (MISTAKE-753 and MISTAKE-755 were both "known results typed as new"). The pass
below still missed two items, found by the audit of this packet and marked **(audit)**: the same class of error a
third time (MISTAKE-760).

- **The `H mod 4` lane.**
  - `THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts`, `THM-478-digit-one-degree-law`:
    `H = 1 + 2 c_odd (mod 4)`, `c_odd` the number of directed odd cycles; `c_3 mod 2` is a quadratic polynomial in the
    arcs.
  - `THM-552-c3-parity-dichotomy-self-converse` (kind-pasteur, 2026-06-20): a `phi`-self-converse tournament of even
    order has `c_3` even.
  - `THM-791-H-companion-laws-to-the-transitivity-flow` (kind-pasteur, 2026-07-14): a tournament with a
    fixed-point-free involutive anti-automorphism has `c_odd` even, so `H = 1 (mod 4)`.
  - **(audit)** `session29_framed_profiles_20261022/session29-converse-congruences/REPORT.md`, sections 1 and 3: a
    self-converse tournament has an anti-involution, all anti-involutions are conjugate, and
    `H = 1 + 2F + 4I (mod 8)`, at even order `H = 1 + 4I (mod 8)` with `I` the number of mirror pairs of odd cycles
    whose supports meet. `session130_blue_black_theorems_20261031/hamiltonian/PROOFS.md`, Theorem 1: the same
    congruence. The audit K1 of this machine's THM-4930 (2026-11-08) had pointed to the first report.
    `session29-halfpath/REPORT.md` section 5, THM-4743 (6), THM-4930: `2^(m^2)` (`n = 2m`) and `2^(m(m+1))`
    (`n = 2m + 1`) labelled tournaments for a fixed anti-involution; `selfk_type_ledger_opus_S312.out` prints
    `704, 7040, 380416` and the products `960, 7680, 430080`. `session29-selfconverse-asymptotics`, THM-283: the
    class counts (A002785).
  - **(audit)** `THM-555-tiling-cycle-moment-hierarchy` (kind-pasteur S20, 2026-06-20):
    `E[(-1)^c3] = 1/2^floor((n-1)/2)` in the tiling ensemble, verified for `n <= 10`;
    `HYP-2707-tournament-clifford-magic-hierarchy` (kind-pasteur S22): the mechanism, a Gauss sum of a quadratic
    form of rank `2 floor((n-1)/2)`, exact for `n = 4, ..., 8`. THM-559 (this machine, June): `c_3` as a 2-body form.
    The classical formula `c_3 = C(n,3) - sum_v C(s_v, 2)` (Kendall and Babington Smith; THM-118, THM-466).
  - THM-4922 (4) (kind-pasteur S16; proof by its audit): at every odd order the self-converse classes split evenly
    between `H = 1` and `H = 3 (mod 4)`. THM-4930: which values are ever self-converse.
  - `THM-1980-the-two-adic-edge-of-the-redei-count`: `H mod 4` is a function of the characteristic polynomial for
    `n <= 6` and not for `n = 7` (exhaustive; no formula).
  - `THM-118-trace-cycle-identity` (kind-pasteur, March): `tr A^k = k c_k` for `k = 3, 4, 5`. THM-500, THM-502: the
    odd-cycle count is spectral iff `n <= 6`; `tr A^7 = 7 (c_7 + TQ)`. Theorem P starts from these.
  - The score-parity vector: uniform on its coset for every `n`
    (`session72_even_graph_unifying_carriers_20261026/session72-even/SCORE_PARITY_QUOTIENT.md`); inside a switching
    class for odd `n` (THM-1470, klein S338; THM-1460-babai-cameron, opus S409).
  - On record before this note: `D_4, ..., D_7` and `E_5[{5}] = 7/32` (THM-466), the class race `32 : 24` at `n = 6`
    (THM-1980), a partial law for `H mod 4` at `n = 5, 6` (HYP-6935 [ii]), the matrix `A` over `F_2` and de Caen's
    rank bound (THM-4857, this machine), the integer vector `A^2 1` as the Kendall-Wei second-order score (THM-896;
    THM-555 notes that it does not determine `c_5`), `T = 1/(1 - S)` for labelled tournaments (THM-4816, THM-4817),
    the case `k = 4` of the subset double count (THM-1455).
  - HYP-9355 and `h_spectrum_multiplicities_20261021.md` section 3 (opus S50): residues of `H`; the "2-adic pull".
  - THM-5127 (this machine, yesterday): `H mod 4` on a vertex cube is a quadratic form; exact `D_n` for `n <= 11`;
    `1/4 <= P(H = 1 mod 4) <= 3/4`.
  - Today, `spectrum_multiplicity_signed_tables_20261111/WEIGHTED_ATLAS.md` (codex): `H = 4x + sigma`, the sign
    `sigma` multiplies over ordered sums and the cyclic triangle acts by `(x, sigma) -> (3x + sigma, -sigma)`.
- **The finite-exception lane.** THM-1370 (`h-spectrum-omits-7-21-all-n`), THM-4889 (exact spectra to order 11,
  strong windows, top band, zone law), THM-4986 (crevasse lemma), THM-4744 and THM-4810 (strong multiplicities),
  the ATLAS of kind-pasteur S15 (exception sets side by side, the Poisson heuristic and its failure for correlated
  families, three ways proved rows close the gap), THM-4026/4027/4028 (Sun's `2-4-6-8`).
- **The bias lane.** `prime_power_bias_collatz_collector_20261014.md` (opus S32), which records the owner's prompt
  ("... a large microcosmic region from 0 to a large number, and then a macrocosmic region ... It has to do with
  prime POWERS. Every odd prime, when squared, lands in residue 1 mod 4") and separates Chebyshev type (a fixed
  square term) from Mertens type (none). THM-4860 (klein): a square-sum arrangement of the first `k` primes needs
  Chebyshev's bias in `{0, 1, 2}`.
- **The square-sum lane.** THM-4505, THM-4510, THM-4850, THM-4851, THM-4870 (C); Session211
  (`session211/square-sum-port-bridge/PROOF.md`: no map from square-sum graphs to tournaments can send the path
  counts `p(Q_23) = 3`, `p(Q_25) = 10` to `H = 2p + 1`, because that would be 7 and 21).
- **The fold lane.** `cpq_defect_family_mirror_folds_20261104.md` section 5 (this machine): `x^(1/x)` and its
  mirror map, the MRB constant, the ball volumes `V(x) = pi^(x/2)/Gamma(x/2 + 1)`, mirror dimensions; the owner's
  prompt of 2026-11-04 ("the same values realized twice on either side of the peak ... small changes in the first
  are big changes in the latter"; "`(pi^(n/2))/(n/2)!` ... has a maximum integer value around 5, like the ratio of
  8 to 4 of isomorphism classes in tournaments at n = 5"). `session127_gsma_norm_collapse_20261031/volume/THEORY.md`
  (VOL-1).
- **The micro/macro lane.** THM-4900 (gate and doubling); THM-5127 (C) (7): the gate preserves `H mod 4`.
- **Canonical hostile example.** The owner's own "8 to 4" at `n = 5`: my note of 2026-11-04 typed its link to the ball
  volume as NUMEROLOGY. Section 8 keeps that typing for the number 5 and replaces the vague link by an exact one
  (labelled self-converse tournaments are the majority of labelled tournaments exactly for `n <= 5`).
- **Corrected near miss.** HYP-9355's "absent at `n = 11`": true at the resolution of its test, false exactly.
- **Least-used sidecar.** THM-552 / THM-791 and the Session29 congruence: parity laws for self-converse
  tournaments that had not been weighed against the size of the race.
- **Same-day sessions** (read before writing): Session249-262, codex spectrum multiplicity, odd-tables,
  decomposition16, kind-pasteur S22 (THM-5130), opus S91. None treats the distribution of `H mod 4`.

## 2. Notation

- `H(T)`: number of Hamiltonian paths; `c_k(T)`: number of directed `k`-cycles; `c_odd = c_3 + c_5 + ...`;
  `chi(T) = (-1)^(c_odd(T)) = (-1)^((H(T) - 1)/2)`.
- `D_n = E chi(T)` over uniformly random labelled tournaments on `n` vertices, so
  `D_n = P(H = 1 mod 4) - P(H = 3 mod 4)`. For a set `K` of odd lengths, `E_n[K] = E (-1)^(sum_{k in K} c_k)`.
- `R + v_X`: add a vertex `v` with out-neighbourhood `X` to `R`. `s_R(a)`: out-degree (score) of `a` in `R`.
- A tournament is *self-converse* if it is isomorphic to its converse; an *anti-automorphism* is a bijection
  `phi` of the vertices with `a -> b` iff `phi(b) -> phi(a)`.
- `SC_n`: number of labelled self-converse tournaments on `n` vertices. `k!!` is the double factorial.
- Over `F_2`: `A` the adjacency matrix (`A_ab = 1` iff `a -> b`), `1` the all-ones vector, `J = 1 1^T`, `|y|` the
  weight of a vector, `p . r` the standard pairing. `s = A1`, `t = As` (section 6).

## 3. Theorem T: the triangle race

**Theorem T.** Let `R` be a tournament on `m >= 0` vertices.

1. `c_3(R + v_X) = c_3(R) + e_R(X, V \ X)`, where `e_R(X, Y)` is the number of arcs from `X` to `Y`, and
   `e_R(X, V \ X) = sum_{a in X} s_R(a) - C(|X|, 2)`.
2. `sum_{X subset V} (-1)^(c_3(R + v_X)) = 2^ceil(m/2)`.
3. For uniformly random labelled tournaments on `n >= 1` vertices, `E (-1)^(c_3) = 2^(-floor((n-1)/2))`.

So every tournament on `m` vertices has exactly `2^(m-1) + 2^(ceil(m/2) - 1)` one-vertex extensions with an even
number of cyclic triangles, whatever its own number of triangles is.

*Proof.* (1) A cyclic triangle through `v` is `v -> a -> b -> v` with `a in X`, `b notin X`, `a -> b`. The arcs
leaving the vertices of `X` number `sum_{a in X} s_R(a)`, and `C(|X|, 2)` of them stay inside `X`.

(2) Let `sigma_a = s_R(a) mod 2`, and let `p` and `q` be the numbers of vertices of even and odd score. Since
`(-1)^C(k,2) = Re((1 - i) i^k)`,

    G(R) := sum_X (-1)^(e_R(X, V \ X)) = Re( (1 - i) prod_a (1 + (-1)^(sigma_a) i) )
          = Re( (1 - i)(1 + i)^p (1 - i)^q ) = 2^((m+1)/2) cos((p - q - 1) pi/4).

Claim: `(-1)^(c_3(R)) G(R)` is the same for all tournaments on `m` vertices. Reverse one arc `a -> b`. Then
`c_3` changes by `s(a) - 1 - s(b)` (the triangles through the arc are counted by the common structure of the two
out-neighbourhoods; the difference is `(s(a) - 1 - t) - (s(b) - t)` with `t = |N+(a) cap N+(b)|`), and the scores
of `a` and `b` change by one each. If `s(a)` and `s(b)` have the same parity, `c_3` changes parity and `p - q`
changes by `+-4`, which changes the sign of the cosine; if they have different parities, neither changes. So
`(-1)^(c_3) G` is invariant under arc reversals, hence constant. For the transitive tournament `c_3 = 0` and the
scores are `0, ..., m - 1`, so `p - q - 1` is `-1` (`m` even) or `0` (`m` odd), and `G = 2^(m/2)` or
`2^((m+1)/2)`.

(3) Average (2) over `R` on `m = n - 1` vertices. QED

*Typing* (corrected by the audit, M2 and m1). (1) is the first layer of THM-5127 (4.4). The closed form (3) is
KNOWN as a verified formula: THM-555 (kind-pasteur S20, June) records `E[(-1)^c3] = 1/2^floor((n-1)/2)` for the
tiling ensemble (tournaments containing a fixed Hamiltonian path), verified for `n <= 10`, and HYP-2707 gives the
mechanism (a Gauss sum of a quadratic form of rank `2 floor((n-1)/2)`, exact for `n = 4, ..., 8`). New here: the
proof for every `n`, in the uniform ensemble, and the pointwise statement (2), which was not found in the
repository. The two ensembles agree for `n >= 4` (the complement of a Hamiltonian path in `K_n` is connected, so the
score-parity vector stays uniform on its coset; at `n = 3` the tiling value is 0), so Theorem T also proves THM-555's
line for every `n >= 4`.

Theorem T is an elementary corollary of the classical formula `c_3 = C(n,3) - sum_v C(s_v, 2)` (Kendall and
Babington Smith; THM-118, THM-466): since `sum_v s_v^2 = q (mod 4)` for `q` the number of odd scores and
`sum_v s_v = C(n,2)`, one gets `c_3 = C(q, 2) + floor(n/4) (mod 2)` (the first formula of Theorem P (1)), and (2)
follows because `x -> (s + 1 + x, |x|)` is a bijection from the out-sets onto the admissible score-parity vectors of
`R + v_X` (the auditor's proof, `audit/AUDIT.md` m1; the proof above by arc reversal is a second one). A web search
found the tables of tournaments by number of circular triads (Kendall and Babington Smith; Harris and McShane,
J. Integer Sequences 27 (2024)) and not the parity statement; no claim about the literature is made.

*Consequence for the race.* For `n <= 4` there are no longer odd cycles, so `chi = (-1)^(c_3)` and
`D_3 = D_4 = 1/2`: the two orders at which the bound `|D_n| <= 1/2` of THM-5127 (C) (6) is attained are the
triangle race.

## 4. Theorem S: the square term

**Theorem S.**

1. (KNOWN: Session29 converse-congruences section 1 and THM-4930 for the involution; THM-791 for the parity;
   THM-552 for `c_3`.) If `T` has an anti-automorphism, it has one of order 2. An involutive anti-automorphism
   fixes at most one vertex, so none when `|T|` is even; then no odd cycle is fixed by it, `c_odd(T)` is even and
   `H(T) = 1 (mod 4)`.
2. (Mod 8. KNOWN: Session29 converse-congruences section 3, `H = 1 + 4I (mod 8)`; Session130 `hamiltonian/PROOFS.md`
   Theorem 1.) Let `psi` be a fixed-point-free involutive anti-automorphism of `T`. It acts on the directed odd
   cycles without fixed points. Let `s` be the number of orbits `{C, psi(C)}` with `V(C) cap V(psi(C))` non-empty.
   Then `H(T) = 1 + 4s (mod 8)`.
3. (Count; a corollary of the KNOWN count `2^(n^2/4)`, resp. `2^((n^2-1)/4)`, of labelled tournaments with a fixed
   anti-involution: Session29 half-path section 5, THM-4743 (6), THM-4930.) For even `n`:
   `2^(n^2/4) <= SC_n <= (n-1)!! 2^(n^2/4)`. For odd `n`: `SC_n <= n!! 2^((n^2-1)/4)`.
4. (Race; the weighing is new.) For even `n`, `D_n = SC_n / 2^C(n,2) + D_n^ns`, where `D_n^ns` is the signed fraction
   of the non-self-converse tournaments, and `0 < SC_n / 2^C(n,2) <= (n-1)!! 2^(-n(n-2)/4)`.

The first draft of this note typed (2) and (3) as new; the audit (M1) found both in the repository, (2) with the
same proof. The proofs are kept below for completeness.

*Proof.* (1) Let `phi` be an anti-automorphism. `phi^2` is an automorphism, and automorphism groups of tournaments
have odd order, so `phi^2` has odd order `k`; `psi = phi^k` is an anti-automorphism and `psi^2 = (phi^2)^k = 1`.
If `psi` fixes `u != w`, then `u -> w` iff `psi(w) -> psi(u)`, that is `w -> u`: impossible. The number of fixed
points has the parity of `|T|`. If `psi(C) = C` for a directed odd cycle `C`, then `psi` acts on `C` as a
reflection of an odd polygon and fixes a vertex. So the odd cycles pair up.

(2) By THM-466, `H = 1 + 2 alpha_1 + 4 alpha_2 (mod 8)`, `alpha_1 = c_odd`, `alpha_2` the number of unordered pairs
of vertex-disjoint odd cycles. `alpha_1 = 2t`, `t` the number of orbits. `psi` acts on the pairs; a pair
`{C, C'}` is fixed only if `C' = psi(C)` (since `psi(C) = C` is impossible), which needs `C` and `psi(C)` disjoint.
So `alpha_2 = t_0 (mod 2)`, `t_0` the number of orbits with `C`, `psi(C)` disjoint, and
`H = 1 + 4t + 4t_0 = 1 + 4(t - t_0) (mod 8)`.

(3) Fix an involution `psi` with at most one fixed point. The tournaments with anti-automorphism `psi` are given by
one free arc for each orbit of `psi` on pairs `{a, b}`: the `floor(n/2)` pairs `{a, psi(a)}` are fixed, the pairs
through the fixed point (odd `n`) form `(n-1)/2` orbits, and the others form orbits of size 2. The number of
orbits is `n^2/4` (`n` even) and `(n^2 - 1)/4` (`n` odd). There are `(n-1)!!` such involutions for even `n` and
`n!!` for odd `n`, and by (1) every self-converse tournament is counted at least once.

(4) (1) and (3). QED

**Exact values** (`scrace2.c`, `race_tables.py`, `out/race_tables.out`; every class of orders 3-10, labelled weights
`n!/|Aut|`; the weights sum to `2^C(n,2)` at every order). `SC_5, SC_6, SC_7` and the bounds at `n = 5, 6, 7` were on
file (`selfk_type_ledger_opus_S312.out`); `SC_8, SC_9, SC_10`, the shares and the splits below are new. (2) is
checked on all self-converse classes of orders 4, 6, 8 by `check_selfconverse_mod8.py`.

| `n` | `SC_n` | bound of (3) | `SC_n / 2^C(n,2)` | `D_n` | share of `D_n` carried by self-converse |
|---|---|---|---|---|---|
| 4 | 48 | 48 | 3/4 | 1/2 | 150% |
| 6 | 7,040 | 7,680 | 55/256 | 1/4 | 85.9% |
| 8 | 6,737,920 | 6,881,280 | 1645/65536 | 1361/32768 | 60.4% |
| 10 | 31,596,429,312 | 31,708,938,240 | 1928493/2147483648 | 650471/268435456 | 37.1% |

At odd orders `SC_n = 8, 704, 380416, 958578688` (`n = 3, 5, 7, 9`; bounds `12, 960, 430080, 990904320`). The
classes split evenly (THM-4922 (4): `1/1, 4/4, 44/44, 1376/1376`), but their labelled weights do not cancel: the
signed sums are `4, 256, 37824, 20782592`, because the classes with `H = 3 (mod 4)` have more automorphisms. The
self-converse classes of even order split `1/1, 6/6, 82/94, 4448/4336` between `H = 1` and `H = 5 (mod 8)` at
`n = 4, 6, 8, 10`.

**Labelled self-converse tournaments are the majority exactly for `n <= 5`** (labelled fractions
`1, 1, 1, 3/4, 11/16, 55/256, ...`; among classes `2/2, 2/4, 8/12, 12/56, 88/456, ...`, a tie at `n = 4`). This is
the owner's "8 to 4" at `n = 5`: the last order at which the squares outnumber the rest.

## 5. The race, exactly

All numbers are exact; sources: THM-5127 for `D_n` (`n <= 11`; `D_4, ..., D_7` are already in THM-466),
`out/race_tables.out` for the splits. The audit recomputed every entry of this section with its own programs.

**Labelled race in three strata** (in units of `4^(-n)`; `sc` = self-converse, `sym` = non-trivial automorphism
group and not self-converse, `rigid` = neither).

| `n` | `D_n` | `D_n 4^n` | sc | sym | rigid |
|---|---|---|---|---|---|
| 3 | 1/2 | 32 | 32 | 0 | 0 |
| 4 | 1/2 | 128 | 192 | -64 | 0 |
| 5 | 13/32 | 416 | 256 | -80 | 240 |
| 6 | 1/4 | 1024 | 880 | -216 | 360 |
| 7 | 853/8192 | 1706 | 295.5 | -7.0 | 1417.5 |
| 8 | 1361/32768 | 2722 | 1645 | -104.2 | 1181.2 |
| 9 | 63103/8388608 | 1972.0 | 79.3 | -73.0 | 1965.7 |
| 10 | 650471/268435456 | 2540.9 | 941.6 | -54.1 | 1653.3 |
| 11 | 260546929/549755813888 | 1987.8 | not split | not split | not split |

**Class race.** Classes with `H = 1` minus classes with `H = 3 (mod 4)`: `0, 0, 0, 8, 36, 264, 1268, 22600` at
`n = 3, ..., 10`. For even `n` this is (number of self-converse classes) + 2 (lead among converse pairs):
`0 = 2 - 2`, `8 = 12 - 4`, `264 = 176 + 88`, `22600 = 8784 + 13816`. At `n = 4` and `n = 6` the classes that are not
self-converse favour `H = 3`, and the self-converse classes reverse (`n = 6`) or cancel (`n = 4`) the result. (The
split `32 : 24` at `n = 6` is in THM-1980.)

**Strong race** (Proposition 5.1; the unsigned identity `T = 1/(1 - S)` is KNOWN, THM-4816, THM-4817, and the sign
is multiplicative, so this is a one-line consequence). Every tournament is uniquely an ordered sum of strong
tournaments and `H` multiplies over the summands, so with exponential generating functions over labelled structures
`T_chi(x) = 1/(1 - S_chi(x))`, where `T_chi`, `S_chi` have coefficients `sum chi` over all and over strong
labelled tournaments. Inverting with the exact `2^C(n,2) D_n` gives the signed counts of strong tournaments
`-2, 24, 176, 4080, 123952, 8135680, 332488064, 75937800960, 15259291222528` at `n = 3, ..., 11`; the same numbers
come out of the census flags for `n <= 10` (`race_tables.py` asserts it). The strong-only lead is positive for
`4 <= n <= 11` (`0.0344, 0.00519, 0.00224, 0.000433` at `n = 8, ..., 11`), and the decomposable tournaments carry
`27%, 36%, 11%, 11%` of the signed sum at these orders. So the lead is not an effect of sources and sinks.

**Residues mod 8 and HYP-9355.** Labelled relative masses of `H = 1, 3, 5, 7 (mod 8)`:

| `n` | 1 | 3 | 5 | 7 |
|---|---|---|---|---|
| 7 | 1.17920 | 0.98480 | 1.02905 | 0.80695 |
| 8 | 1.06395 | 0.97866 | 1.01912 | 0.93827 |
| 9 | 1.01529 | 0.99192 | 0.99975 | 0.99303 |
| 10 | 1.00330 | 0.99603 | 1.00155 | 0.99912 |
| 11 | 1.00051 | 0.99949 | 1.00044 | 0.99956 |

Rows 7-9 reproduce HYP-9355's table; rows 10 and 11 are exact where HYP-9355 had samples (`1.0030 +- 0.0009` at
`n = 10`; `0.9996` pooled at `n = 11`). Row 10: the census of this packet and the first audit of THM-5127. Row 11:
the census of the audit of this packet over all 9,733,056 bases of order 10 (`audit/audit_cube.c`); the labelled
counts are `9011754468980736, 9002583394936832, 9011181642270720, 9003277512775680` of `2^55`. (The first draft had
`1.00052, 0.99948` in row 11, from five-digit fractions of the second audit of THM-5127; audit m2.)

**Two deck remarks (PROVED, double counting).** Each odd cycle of length `l` lies in `n - l` vertex-deleted
subtournaments. So for even `n` the number of vertices `v` with `H(T - v) = 3 (mod 4)` has the parity of
`c_odd(T)`: `H(T) mod 4` is determined by the classes mod 4 of the deck. For odd `n` that number is even. More
generally `Y_k(T) := sum_{|S| = k} c_odd(T[S]) = sum_l C(n - l, k - l) c_l (mod 2)`, and
`Y_k = (p_k(T) - C(n,k))/2 (mod 2)` with `p_k` the number of directed paths on `k` vertices: the races below are
the races of the path counts `p_k mod 4`.

## 6. The race by cycle length: the pentagon theorem

`E_n[K] = E (-1)^(sum_{k in K} c_k)` over labelled tournaments. Section 3 is the column `K = {3}`. This section solves
the columns `{5}` and `{3,5}` at every order, shows that every vertex cube has a non-negative triangle-plus-pentagon
sum, and records what happens from the heptagon on.

### 6.1 Theorem P

For a tournament `T` on `n` vertices with adjacency matrix `A` over `F_2` (`A_ab = 1` iff `a -> b`) put

- `s = A1`: the indicator of the vertices of odd score, `q_1 = |s|`;
- `t = As`: the indicator of the vertices with an odd number of out-neighbours of odd score, `q_2 = |t|`.

Then `q_1 = C(n,2)` and `q_2 = n C(n,2) (mod 2)`.

**Theorem P.**

1. **(Score formulas.)** Modulo 2,

       c_3 = C(q_1, 2) + floor(n/4),    c_3 + c_5 = C(q_2, 2) + floor(n/4),    c_5 = C(q_1, 2) + C(q_2, 2).

   So the parity of the number of 5-cycles is read off from `q_1 mod 4` and `q_2 mod 4`.
2. **(Vertex cube.)** Let `R` have `m` vertices, `x` the indicator of `X`, `sigma = |X| mod 2`, `y = Ax + sigma 1`,
   `u_1 = As`, `u_0 = As + s + 1`. Then modulo 2

       (c_3 + c_5)(R + v_X) - (c_3 + c_5)(R) = x^T (A + A^3)(1 + x) = C(|y|, 2) + u_sigma . y.

3. **(Cube sum.)** Let `G(R) = sum_X (-1)^((c_3 + c_5)(R + v_X))`, `K_0 = ker A^T`, `K_1 = {u : A^T u = 1}` (over
   `F_2`). Every element of `K_0` has weight `0` or `1 (mod 4)`, every element of `K_1` weight `0` or `3 (mod 4)`.
   Write `n_ij` for the number of elements of `K_i` of weight `j (mod 4)`. Then

   | `m mod 4` | `G(R)` |
   |---|---|
   | 0 | `2^(m/2) |K_0|` |
   | 1 | `2^((m-1)/2) (|K_0| + n_13 - n_10)` (which is `2^((m-1)/2) |K_0|` if `K_1` is empty and `2^((m+1)/2) n_13` otherwise) |
   | 2 | `2^(m/2) (n_00 + n_13)` |
   | 3 | `2^((m+1)/2) (n_00 + n_10)` |

   In particular **`G(R) >= 0` for every tournament `R`**: no tournament has more one-vertex extensions with
   `c_3 + c_5` odd than with `c_3 + c_5` even. More precisely `G(R) = 0` or `G(R) >= 2^ceil(m/2)`, and `G(R) = 0`
   only for `m = 1 (mod 4)` (there, if `K_1` is empty then `A` is singular and `|K_0| >= 2`; the sharpening is the
   auditor's). Values at `m = 5`: `0, 8, 16`; at `m = 9`: `0, 32, 64, 128`.
4. **(Joint law.)** For a uniformly random labelled tournament on `n >= 1` vertices, `s` is uniform on
   `{|s| = C(n,2) mod 2}` (KNOWN: Session72 `SCORE_PARITY_QUOTIENT.md`; re-proved in 6.2). Given `s` with support
   `O`, `1 <= |O| <= n - 1`, the vector `t` is
   uniform on `{t : |t restricted to O| = C(|O|, 2), |t restricted to V \ O| = C(n - |O|, 2) (mod 2)}`; `t = 0` if `O`
   is empty and `t = 1` if `O = V`.
5. **(Exact races.)** With `(2/n)` the Kronecker symbol (`+1` for `n = +-1`, `-1` for `n = +-3 (mod 8)`):

   | `n` | `E_n[{3}]` | `E_n[{3,5}]` | `E_n[{5}]` |
   |---|---|---|---|
   | odd | `2^(-(n-1)/2)` | `2^(-(n-1)/2) (2 + (2/n) 2^(-(n-1)/2) - 2^(-(n-2)))` | `2^(-(n-1)) (3 - (2/n) 2^(-(n-3)/2))` |
   | `2 (mod 4)` | `2^(-(n-2)/2)` | `2^(-(n-2)/2)` | `2^(-(n-2))` |
   | `0 (mod 4)` | `2^(-(n-2)/2)` | the value at `n - 1` | the value at `n - 1` |

   So for `n = 2 (mod 4)` the bits `c_3` and `c_3 + c_5` have the same bias and are uncorrelated
   (`E[{5}] = E[{3}] E[{3,5}]`); for every other `n`, `E_n[{3,5}] / E_n[{3}] -> 2` and `E_n[{5}] / E_n[{3}]^2 -> 3`.
   Conditional form: for `1 <= |O| <= n - 1`,
   `E[(-1)^(c_3 + c_5) | s = 1_O] = +-2^(2 - ceil(|O|/2) - ceil((n - |O|)/2))`, and for `n = 2 (mod 4)` this is
   `2^(-(n-2)/2)` for every `O`: the sign `(-1)^(c_3 + c_5)` is then uncorrelated with every function of the score
   parities.

**Corollaries.**

- **The race is solved for `n <= 6`.** There `c_odd = c_3 + c_5`, so `H = 1 + 2 (C(q_2, 2) + floor(n/4)) (mod 4)`:
  the residue of `H` mod 4 of a tournament on at most six vertices is decided by the number of vertices that have an
  odd number of odd-score out-neighbours, and `D_5 = 13/32`, `D_6 = 1/4` are instances of (5). This is the formula
  behind the exhaustive fact `d(6) = 2` of THM-1980 (`H mod 4` is spectral for `n <= 6`; THM-500 has the reason
  `c_3 + c_5 = tr A^3/3 + tr A^5/5`).
- **Order 11.** `E_11[{5}] = 49/16384`, `E_11[{3,5}] = 1007/16384` (also obtained from the census of all 9,733,056
  bases of order 10).
- **Kernel reading.** For `m = 0 (mod 4)` the triangle-plus-pentagon surplus of the extensions of `R` is
  `2^(m/2)` times the number of sets `U` such that every vertex has an even number of in-neighbours in `U`.
- **Expected kernel size** (PROVED, elementary; checked exhaustively for `m <= 6`). For a uniformly random labelled
  tournament on `m` vertices, `E |ker A| = 2 + 2^((1-m)/2) cos((m-1) pi/4) - 2^(1-m)` over `F_2` (a uniformly random
  square matrix has `2 - 2^(-m)`): a non-zero `x` of weight `k` is in the kernel with probability `2^(1-m)` if
  `k = 0, 1 (mod 4)` and `0` otherwise, because the scores inside the support must all be even. With (3) this gives
  the row `n = 1 (mod 4)` of (5) a second time.

### 6.2 Proofs

*Preliminaries.* (i) `A + A^T = J + I` over `F_2`. (ii) For `y in F_2^n`, `y^T A y` is the number of arcs inside the
support of `y`, so `y^T A y = C(|y|, 2)`; polarising with (i),
`C(|p + r|, 2) = C(|p|, 2) + C(|r|, 2) + |p||r| + p . r (mod 2)`. (iii) KNOWN (THM-118): a closed walk of length 3 or 5
in a tournament is a cycle, so `c_3 = tr A^3`, `c_5 = tr A^5 (mod 2)`. (iv) A cycle of length `l` through the new
vertex of `R + v_X` is a path on `l - 1` vertices of `R` from `X` to `V \ X`, and a walk `a -> b -> c -> d` with
`a != d` is a path; so the increments of `c_3` and `c_5` are `x^T A (1 + x)` and `x^T A^3 (1 + x)`. (v) With
`pi = C(m,2) mod 2`: `|s| = pi`, `1^T A z = (m-1)|z| + s . z`, and `|As| = 1^T A s = m pi (mod 2)`.

*(2).* Put `h = I + A`. Over `F_2`, `A + A^3 = h A h` and `h^T = A + J` by (i). So
`x^T (A + A^3)(1 + x) = (h^T x)^T A (h (1 + x))`, with `h^T x = Ax + sigma 1 = y` and
`h(1 + x) = 1 + x + s + Ax = y + x + s + (1 + sigma) 1`. Expanding and using (ii) and (v):
`y^T A y = C(|y|, 2)`; `y^T A x = |Ax| + sigma 1^T A x = (1 + sigma) |y|`; `y^T A s = y . As`;
`(1 + sigma) y^T A 1 = (1 + sigma) y . s`. The sum is `C(|y|, 2) + y . (As + (1 + sigma)(s + 1)) = C(|y|, 2) + u_sigma . y`.

*(1).* Let `f_1 = C(q_1, 2)`, `f_2 = C(q_2, 2)`. It is enough to show, for `T = R + v_X`,

    f_1(T) - f_1(R) = x^T A (1 + x) + m pi,        f_2(T) - f_2(R) = x^T (A + A^3)(1 + x) + m pi     (mod 2),

because `floor((m+1)/4) - floor(m/4) = m pi (mod 2)` (both sides are 1 exactly for `m = 3 (mod 4)`), all three
quantities vanish on the one-vertex tournament, and every tournament is built by adding vertices. The score
vector of `T` is `s' = (s + 1 + x, sigma)`, and `t' = A' s' = (c_0 + y + sigma x, s . x)` with `c_0 = As + s`.

For `f_1`: `C(|s'|, 2) = C(|s + 1 + x|, 2) + sigma |s + 1 + x|`, and by (ii) twice
`C(|s + 1 + x|, 2) = C(|s|, 2) + [pi + C(|x|, 2) + m sigma + sigma] + pi (m + sigma) + pi + s . x`; adding
`sigma (pi + m + sigma)` leaves `C(|s|, 2) + C(|x|, 2) + s . x + m pi`, and `x^T A (1 + x) = s . x + C(|x|, 2)`.

For `f_2`: `C(|t'|, 2) = C(|c_0 + y + sigma x|, 2) + (s . x) |c_0 + y + sigma x|`. By (ii) and (v):
`C(|c_0|, 2) = C(|As|, 2) + m pi`; `C(|y + sigma x|, 2) = C(|y|, 2) + sigma (s . x)`;
`|c_0| = (m + 1) pi`, `|y + sigma x| = s . x (mod 2)`; `c_0 . (y + sigma x) = u_sigma . y + (1 + sigma)(s . x)`. Hence
`C(|c_0 + y + sigma x|, 2) = C(|As|, 2) + m pi + C(|y|, 2) + u_sigma . y + (s . x)(1 + (m + 1) pi)`, and the last
term cancels against `(s . x)|c_0 + y + sigma x|`. By (2) the result is `f_2(R) + m pi + x^T (A + A^3)(1 + x)`.

*(3).* By (2), `(-1)^((c_3+c_5)(R)) G(R) = sum_x (-1)^(C(|y|,2) + u_sigma . y)`, a function of the image of `x` under
`L(x) = (Ax, sigma)`. The orthogonal complement of the image of `L` is `K = {(u, tau) : A^T u = tau 1}`, so with
`Gamma(w) = sum_y (-1)^(C(|y|,2) + omega . y)` for `|omega| = w`, which equals
`Re((1 - i)(1 + i)^(m - w)(1 - i)^w) = 2^((m+1)/2) cos((m - 1 - 2w) pi/4)`,

    (-1)^((c_3+c_5)(R)) G(R) = (1/2) sum_{(u, tau) in K} [ Gamma(|u + u_0|) + (-1)^(tau + |u|) Gamma(|u + u_1|) ].

(Poisson summation: `sum_x F(Lx) = (|ker L| / |K|) sum_{k in K} sum_z (-1)^(k . z) F(z)`, `z = (y', sigma)` running
over `F_2^m x F_2`, and `|ker L| / |K| = 1/2` because `dim K = 1 + dim ker L`; the substitution `y = y' + sigma 1` gives
the sign `(-1)^(sigma (tau + |u|))`.)

For `(u, tau) in K`: `C(|u|, 2) = u^T A u = tau |u|`, which gives the weights mod 4 stated in (3);
`u . As = tau pi` and `u . s = tau m`, so `u . u_1 = tau pi` and `u . u_0 = tau pi + tau m + |u|`. With
`a = |As| mod 4`: `|u_1| = a`, `|u_0| = m - a - pi (mod 4)` (because `|As + s| = |As| + |s| - 2 s^T A s`,
`s^T A s = C(|s|, 2)`, and `k - 2 C(k,2) = k(2 - k)` is `0` or `1 (mod 4)` according to the parity of `k`), and
`|u + u_sigma| = |u| + |u_sigma| + 2 u . u_sigma (mod 4)`.
Evaluating the cosines: the element `0` contributes `Z = (1/2)(Gamma(|u_0|) + Gamma(|u_1|))`, and

| `m mod 4` | `abs(Z)` | `tau = 0`, weight 0 | `tau = 0`, weight 1 | `tau = 1`, weight 0 | `tau = 1`, weight 3 |
|---|---|---|---|---|---|
| 0 | `2^(m/2)` | `Z` | `Z` | 0 | 0 |
| 1 | `2^((m-1)/2)` | `Z` | `Z` | `-Z` | `Z` |
| 2 | `2^(m/2)` | `Z` | 0 | 0 | `Z` |
| 3 | `2^((m+1)/2)` | `Z` | 0 | `Z` | 0 |

and the sign of `Z` is `(-1)^(C(a,2) + floor(m/4))`, which is `(-1)^((c_3+c_5)(R))` by (1). This is the table of (3).

*(4).* `t_b` is the number of out-neighbours of `b` in `O`. For `b in O` it is the score of `b` in `T[O]`; for `b`
outside `O` it is the number of arcs from `b` into `O`. The condition `s = 1_O` says that the number of arcs from
`b in O` to the outside is `1 + t_b`, and that the score of `b` outside `O` in `T[V \ O]` is `t_b`. So the pairs
`(T[O], T[V \ O])` and the arcs between the two parts have prescribed out-degree parities. A tournament on `r >= 1`
vertices with a prescribed score-parity vector of weight `C(r,2) (mod 2)` exists in `2^(C(r,2) - (r-1))` ways, and a
bipartite orientation between parts of sizes `q, n - q >= 1` with prescribed out-degree parities of total
`q(n - q) (mod 2)` in `2^(q(n-q) - (n-1))` ways (reversing an arc changes two parities; the underlying graphs are
connected). The total parity condition for the arcs between the parts is `q = C(n,2) (mod 2)`, which holds. So every
admissible `t` has the same number `2^(C(n,2) - (2n - 3))` of tournaments. The two extreme cases are immediate.

*(5).* By (1) the three signs are functions of `(q_1, q_2)`, and by (4), with
`sum_{k = C(r,2) (2)} C(r,k) (-1)^C(k,2) = (-1)^floor(r/4) 2^floor(r/2)` (real and imaginary parts of `(1 + i)^r`),

    E_n[{3,5}] = (-1)^floor(n/4) 2^(1-n) [ B + 2^(2-n) sum_q C(n,q) (-1)^(C(q,2) C(n-q,2) + floor(q/4) + floor((n-q)/4)) 2^(floor(q/2) + floor((n-q)/2)) ],

the sum over `1 <= q <= n - 1`, `q = C(n,2) (mod 2)`, and `B = [C(n,2) even] + [n = C(n,2) mod 2] (-1)^C(n,2)` the
two boundary terms (`q = 0`, where `t = 0`, and `q = n`, where `t = 1`); `E_n[{5}]` has the extra sign `(-1)^C(q,2)`,
boundary terms `[C(n,2) even] + [n = C(n,2) mod 2]`, and no prefactor. The summand, divided by `C(n,q) 2^(n-2)`, is
the conditional mean given `|O| = q`; this is the conditional form. For `n = 2 (mod 4)` every admissible `q` is
odd, `q` and `n - q` are both `1` or both
`3 (mod 4)`, the exponent of `-1` in the summand is `(n-2)/4 (mod 2)` in both cases, which cancels the prefactor
`(-1)^floor(n/4)`, and `floor(q/2) + floor((n-q)/2) = (n-2)/2`: the conditional
mean is `2^(-(n-2)/2)` for every `O`, which gives the middle row of (5) directly (`(-1)^(c_5)` is `(-1)^(c_3)`, a
function of `s`, times `(-1)^(c_3 + c_5)`).

For the other rows: for `n` in a fixed class mod 8 the summand is `C(n,q)` times a function of `q mod 8` times a
fixed power of 2 (`2^(floor(q/2) + floor((n-q)/2))` depends only on `n` once the parity of `q` is fixed), and
`sum_q C(n,q) g(q mod 8)` is a combination of the powers `(1 + zeta^j)^n`, `zeta = e^(2 pi i/8)`. So, as functions of
`k` in `n = 8k + r`, both sides of each identity in (5) are combinations of at most five geometric sequences, with
ratios `2^(-4), 2^(-8), 2^(-12), -(2 + sqrt 2)^4 2^(-12), -(2 - sqrt 2)^4 2^(-12)` (from `2^(3-2n) 2^(n/2) (1 + zeta^j)^n`
for `j = 0, +-2, +-1, +-3`, from the boundary terms and from the closed forms). Two such combinations that agree at
five consecutive values of `k` agree for all `k`. `pentagon_closed_forms.py` checks the identities in exact rational
arithmetic for all `n <= 130`, at least sixteen members of each class; so they hold for every `n`. (The audit
checked that both sides satisfy the recurrence with these five roots, that only the three rational ratios occur,
and the identities for `n <= 240`.) QED

*Typing.* (1)-(5) PROVED. KNOWN and used: THM-118 (`tr A^k = k c_k`, `k <= 5`; also in the literature), THM-500,
the uniformity of the score-parity vector for every `n` (Session72 `SCORE_PARITY_QUOTIENT.md`; for odd `n` inside a
switching class THM-1470), the classical score formula for `c_3` (which gives the first formula of (1) directly),
the layers of the vertex-cube form (THM-5127 (4.4)). RELATED: over the integers `A^2 1` is the Kendall-Wei
second-order score (THM-896; THM-555 notes that it does not determine `c_5`); the matrix `A` over `F_2` and de Caen's
rank bound are in THM-4857; HYP-6935 [ii] has a partial law for `H mod 4` at `n = 5, 6`. Not found in the repository
(the audit searched as well): the reduction mod 2 of the second-order score and its weight mod 4, the formulas of
(1) for `c_3 + c_5` and `c_5`, the second equality of (2), the cube formula (3) and its positivity, the joint law
(4), the closed forms (5), the corollary for `n <= 6`. No claim about the literature (the audit found nothing on
the parity of `c_5` or on the kernel formula; `audit/literature_notes.md`).

*Machine checks* (`k35cube.c`, `k35_tables.py`, `pentagon_closed_forms.py`; `out/k35_summary.out`,
`out/pentagon_closed_forms.out`): (1) and (4) on every labelled tournament of order at most 6; (1), (3), the trace
formulas and the weights of `K_0`, `K_1` on all 9,932,003 classes of orders 2-10 and on 120,000,000 random labelled
bases of orders 11-40, without a failure; the exact labelled sums of the censuses equal the closed forms of (5) at
every order `n <= 11`.

### 6.3 The table, and what the heptagon does

Exact values for every set `K` and `n <= 10` (`cycrace.c`, `cycle_tables.py`, `out/cyc_n.txt`,
`out/cycle_tables.out`; control: `c_3` equals the score formula on every class, the triangle column equals Theorem T,
the full column equals `D_n`); row 11 from Theorem P and THM-5127.

| `n` | `{3}` | `{5}` | `{3,5}` | `{7}` | `{3,7}` | `{5,7}` | `{3,5,7}` | all |
|---|---|---|---|---|---|---|---|---|
| 5 | 1/4 | 7/32 | 13/32 | | | | | 13/32 |
| 6 | 1/4 | 1/16 | 1/4 | | | | | 1/4 |
| 7 | 1/8 | 11/256 | 67/256 | 811/8192 | 1553/8192 | 461/8192 | 853/8192 | 853/8192 |
| 8 | 1/8 | 11/256 | 67/256 | 1277/32768 | 1361/32768 | 437/32768 | 1361/32768 | 1361/32768 |
| 9 | 1/16 | 23/2048 | 263/2048 | 16343/2^21 | -1453/2^21 | 9019/2^18 | 35881/2^20 | 63103/2^23 |
| 10 | 1/16 | 1/256 | 1/16 | 35881/2^24 | 9019/2^22 | 9019/2^18 | 35881/2^20 | 650471/2^28 |
| 11 | 1/32 | 49/16384 | 1007/16384 | | | | | 260546929/2^39 |

- **Short cycles keep the bias, the longest cycles erase it.** At `n = 9`: triangles `0.0625`, triangles and
  pentagons `0.128`, lengths up to 7 `0.034`, all lengths `0.0075`. At `n = 10`: `0.0625, 0.0625, 0.034, 0.0024`.
  The columns `{3}` and `{3,5}` are permanent laws of size `2^(-n/2)` (Theorems T and P; the column `{5}` is of
  size `2^(-n)`); the full race is about `10^3 * 4^(-12)` at `n = 12` (section 7).
- **Why the first two layers are exact.** `c_3` and `c_3 + c_5` are weights mod 4 of two vectors that depend
  `F_2`-linearly on the arcs in a controlled way (`s`, then `t` given `s`), and a weight mod 4 has a Gauss sum of
  modulus `2^(n/2)`. **The direct continuation fails at the heptagon** (`tower_test.py`): the next candidate,
  `tr(A^3 + A^5 + A^7) = C(|A^3 1|, 2) + floor(n/4)`, fails for 42% of random tournaments of orders 3-14, and
  `c_7` is not a trace (THM-500, THM-502). Whether another vector does the job is open (HYP-9663 (H) (3)). From
  the heptagon on the increment of the count on a vertex cube is
  still a quadratic form (THM-5127), but its base sign is no longer a weight.
- **FINITE-EXACT, not explained.** At `n = 10`: `E[{7}] = E[{3,5}] E[{3,5,7}]` and `E[{3,7}] = E[{3,5}] E[{5,7}]`
  (the sign of `c_3 + c_5` is uncorrelated with those of `c_3 + c_5 + c_7` and `c_5 + c_7`, and with no other sign
  involving `c_7` or `c_9`); `E[{5,7}]`, `E[{3,5,7}]` have the same values at `n = 9` and `n = 10`. One order of
  evidence; recorded as a question in HYP-9663 (H). For `m = 1 (mod 4)` the cube sum vanishes
  for a fraction `0.2206 +- 0.0002` of random bases at `m = 17, 21, ..., 37` (`0.2177` at `m = 13`).
- **Subset form** (PROVED, double counting). A cycle of length `l` lies in `C(n-l, k-l)` of the `k`-subsets, so
  `Y_k = sum_l C(n-l, k-l) c_l (mod 2)` is the residue `(p_k - C(n,k))/2` of the number of `k`-vertex paths:
  `Y_5 = c_5 + C(n-3, 2) c_3` is `c_3 + c_5` for `n = 1, 2 (mod 4)` and `c_5` for `n = 0, 3 (mod 4)`. So Theorem P
  (5) is the exact law of `p_5 mod 4` at every order. (For the record, `anf5.py`: in the arc variables of a labelled
  5-set the normal form of `c_odd mod 2` has 120 monomials - 3 arcs, 25 cherries, 2 pairs of disjoint arcs, 20 paths
  on four vertices, 10 cherries with a disjoint arc and all 60 Hamiltonian paths of `K_5`.)

### 6.4 The four pre-registered predictions about this table

| | prediction | registered | outcome |
|---|---|---|---|
| R5 | "the pentagon is invisible at even order": `E_10[{3,7,9}] = D_10` | before Run B | **FAILED** (`533081/2^28` against `650471/2^28`) |
| R6 | `E_10[{5}] = 23/2048`, `E_10[{3,5}] = 263/2048` (the values at `n = 9`) | before Run B | **FAILED** (`1/256`, `1/16`) |
| R7 | `E_14[{3,5}] = 2^(-6)`, `E_18[{3,5}] = 2^(-8)` | after R5-R6 failed, before Run C | **CONFIRMED** (Run C: `0.01580`, `0.00382`, `+1.1` and `-0.5` s.e.), then PROVED (Theorem P (5)) |
| R8 | `E_12[{3,5}] = E_11[{3,5}]`, `E_16[{3,5}] = E_15[{3,5}]` | same | **CONFIRMED** (differences `+0.7`, `+0.6` s.e.), then PROVED |

R5 and R6 extrapolated two coincident orders (6 and 8) and were wrong; R7 and R8 extrapolated the residue of `n`
mod 4 and were right. Run C is plain sampling (`c35mc.c`, 40 million tournaments per order 11-18,
`out/c35_summary.out`; all 16 comparisons with the exact values of Theorem P are within 2.3 standard errors). The
proof was found afterwards, by looking for the reason behind a control: the cube sum of `c_3 + c_5` was never
negative.

## 7. The generic lead (HYP-9663 (G)) and the order-12 run

**Run A** (pre-registered R1-R4; `biasmc.c`, `merge_mc.py`, `out/mc11_*.out`, `out/mc11_summary.out`). For a
uniformly random labelled base `R` of order 11 the program computes `eps_R b_R` exactly (the Gauss sum over the
`2^11` out-sets, THM-5127 (C)); its mean is `D_12`, and the per-sample variance is `E[b_R^2]`, about `2^(-11)`, instead
of 1 for plain sampling. `N = 14,400,000` (six seeds).

| quantity | registered | observed | verdict |
|---|---|---|---|
| R1: `D_12 4^12` | in `[1500, 2500]`, positive | `988.5 +- 98.5` (`D_12 = 5.892 * 10^(-5) +- 5.9 * 10^(-6)`, ten standard errors above 0; six runs: 1017, 616, 1101, 684, 1236, 1276) | **FAILED** (5.2 standard errors below the interval); the sign holds |
| R2: rank law at order 11 | rank 10 within `0.006` of `0.83898`; rank at most 6 in `[0.0012, 0.0030]` | `0.838685 +- 0.000097`; `0.001398` (random alternating with an odd number of edges: `0.83898 / 0.15965 / 0.00137`; observed `0.83869 / 0.15992 / 0.00140`) | **CONFIRMED** |
| R3: `rho_11` | in `[1.015, 1.045]` | `1.0172` | **CONFIRMED** |
| R4: corank-1 bias frequency | in `[0.500, 0.512]` | `0.50399 +- 0.00014` | **CONFIRMED** |
| control | no failure | no failure of Theorem C (1)-(3) or of the parity law (4.5) in 14,400,000 samples | holds |
| estimator check (not registered; run after Run A) | - | the same program one order lower, 4,000,000 bases of order 10: `D_11 4^11 = 1981.5 +- 66.8` (exact `1987.8`), `rho_10 = 1.0399` (exact `1.0406`), `E abs(b_R) = 0.022552` (exact `0.02256`); `out/mc10_control_summary.out` | agrees |

Also from the run: `E_R abs(b_R) = 0.014543` (so `|D_12| <= 0.0146` rigorously up to sampling error), and the
vertex-deletion correlation `E_R[b_R] = 1.095 * 10^(-3)`, of which `2^(1-11) = 9.77 * 10^(-4)` is forced by sources and
sinks.

**Reading.** In units of `4^(-n)` the lead is `1706, 2722, 1972, 2541, 1988` at `n = 7, ..., 11` (exact) and
`989 +- 99` at `n = 12`. The four values near `2 * 10^3` were a plateau, not a law: between `n = 11` and `n = 12` the
lead falls by a factor 8, against factors `5.5, 3.1, 5.1` before. After removing the square term (at most
`10395 * 2^(-30)`, i.e. 162 in these units at `n = 12`) the non-self-converse lead is about `830` at `n = 12` against
`1599` at `n = 10`. What survives of HYP-9663 (G) is the sign: `D_n > 0` at every order computed. The rate of decay
is OPEN: one order shows that the plateau does not continue, not what replaces it (the values at `n = 7, ..., 11`
already move by a quarter from one order to the next).

**Independent estimate (audit).** The auditor wrote a separate estimator (own path DP, residues of
`H(R + v_X) mod 8` for all `2^11` out-sets, other generator and seeds; `audit/audit_cube.c`): 9,000,000 bases of order
11 give `D_12 4^12 = 1039 +- 125`, and the same program one order lower gives `2038.9 +- 66.8` against the exact
`1987.8`. Pooled with Run A: `D_12 4^12 = 1008 +- 77`. R1 FAILED stands.

**On R2.** By its registered rule (tolerance `0.006`) R2 is confirmed, and the agreement is to `3 * 10^(-4)`; but
that tolerance was sixty standard errors wide, and the observed fractions differ from the reference law by about
three standard errors (`-3.1`, `+2.8`). The rank law at order 11 is close to, and probably not equal to, the law of
a random alternating matrix with an odd number of edges.

The contrast between R1 and R2-R4 is the content: the quantities that passed are second moments and rank
statistics, which are governed by forced coincidences (twins, sources and sinks) and by a random-matrix law; the
quantity that failed is a mean, where a few small orders had drawn a straight line.


## 8. The merge: where small numbers stop being special

Each row states the small-number behaviour, the large-number behaviour, what carries the small-number regime (the
"microcosm" of the owner's 2026-10-14 prompt), where it ends, and a type. Types: **crevasse** (a structural gap
closed by overlap), **Borel-Cantelli** (a probability whose sum converges), **square term** (a duality whose fixed
points all fall in one class), **fold** (a unimodal sequence: slower growth over faster growth), **permanent bias**
(an exact Gauss sum that neither decays faster nor changes sign), **threshold** (a monotone ratio that crosses
1/2). Rows 1-7 are inherited; rows 8-13 are from this session or re-typed by it.

| # | phenomenon | small numbers | large numbers | carried by | ends at | type | status and source |
|---|---|---|---|---|---|---|---|
| 1 | odd path counts of tournaments | 7 and 21 are missed | every other odd number (conjecture) | 7: the bands `[f(n), max H(n)]` of strong tournaments, `[3,3], [5,5], [9,15], [15,45]` at orders 3-6, leave the gap between 5 and 9; 21: inside the band of order 6 but not a strong value there, with `f(7) = 25` and `21 = 3 * 7` (THM-1370 (4)-(5)) | value 21, orders 6-7 | crevasse | exclusion PROVED (THM-1370, THM-4986); coverage OPEN |
| 2 | square-sum paths on `1..n` | none for `n in {2..14, 18..22, 24}` | all `n >= 25` | leaves and the degree law | 24 | crevasse, closed by a self-similar induction | PROVED (Gerbicz; THM-4505, THM-4851) |
| 3 | Fermat numbers | prime for `k <= 4` | composite for `5 <= k <= 32`; finitely many primes expected | prime density `1/ln F_k` not yet negligible | `k = 4` | Borel-Cantelli | OPEN (THM-4847, ATLAS) |
| 4 | Sun's `2-4-6-8` | true below `8.96 * 10^14` | false | representation count `V X^(25/24)`, critical | `896315812331399` | Borel-Cantelli, critical | REFUTED (THM-4026, THM-4028) |
| 5 | primes mod 4 | class 3 leads until 26861 | class 3 leads with logarithmic density 0.9959 (conditional: GRH and linear independence of the zeros, Rubinstein-Sarnak) | prime squares | never ends | square term | CITED (opus S32 note) |
| 6 | Polya, `pi < li` | hold to `9.06 * 10^8`, `10^19+` | fail | squares | first failure | square term | CITED (same) |
| 7 | Mertens | holds to `10^16` | fails | nothing (`mu` kills squares) | - | no square term | CITED (same) |
| 8 | triangle parity of tournaments | even leads by `2^(-floor((n-1)/2))` | the same law | every base has the same surplus `2^ceil(m/2)` | never ends | permanent structural bias | PROVED (Theorem T) |
| 8b | triangle-plus-pentagon parity | even leads: `1/2, 1/2, 13/32, 1/4, 67/256, ...` | `2^(-(n-2)/2)` for `n = 2 (mod 4)`, about `2 * 2^(-floor((n-1)/2))` otherwise | every base has a non-negative surplus, a multiple of the kernel size of its adjacency matrix mod 2 | never ends | permanent structural bias | PROVED (Theorem P) |
| 9 | `H mod 4` of tournaments, labelled | `D_3 = D_4 = 1/2`, the extreme allowed; `D_5 = 13/32`, `D_6 = 1/4` | positive lead, `2 * 10^3 * 4^(-n)` for `8 <= n <= 11`, then `10^3 * 4^(-12)` | triangles (`n <= 4`) and pentagons (`n <= 6`), where the race is a weight mod 4; then self-converse tournaments (even `n <= 8`) | 6 for the exact law, about 10 for the square term | weight mod 4 and square term in the microcosm, generic lead after | PROVED for `n <= 6` (Theorems T, P); exact to 11, estimated at 12 (sections 5, 7); rate OPEN |
| 10 | `H mod 4`, classes | lead `0, 0, 0, 8` at `n = 3..6`; at 4 and 6 the non-self-converse classes favour `H = 3` | lead about `2^C(n,2) D_n / n!` | self-converse classes | `n = 8` | square term | FINITE-EXACT (section 5) |
| 11 | labelled self-converse share | majority for `n <= 5` (`11/16` at `n = 5`; 8 of 12 classes) | at most `(n-1)!! 2^(-n(n-2)/4)` | `(2k-1)!!` against `2^(k^2)`; the share is monotone and passes 1/2 between `n = 5` and `6` | 5 | threshold | PROVED bound (Theorem S (3), a corollary of known counts) |
| 12 | "the pentagon is invisible at even order" | true at 6 and 8 | false at 10 | at 6 it is the case `n = 2 (mod 4)` of Theorem P; at 8 unexplained | 10 | coincidence of two cases | pre-registered, FAILED (section 6.4) |
| 12b | "`D_n 4^n` is about `2 * 10^3`" | true at 8, 9, 10, 11 | false at 12 | unknown | 12 | plateau of four cases | pre-registered, FAILED (section 7) |
| 13 | ball volumes, `n^(1/n)` | increase to `n = 5`, to `n = 3` | decrease to 0, to 1 | `pi^(n/2)` against `Gamma(n/2 + 1)` | 5.26, `e` | fold | elementary (fold note of 2026-11-04) |

**Folds and bands (typed ANALOGY with one exact ingredient).** The owner's fold has a value realised twice, once on
each side of a peak. Two sequences in this problem have that shape for the same reason, a slower growth divided by
a faster one:

- `V_2k = pi^k/k!` (ball volumes in even dimension; the weights of a Poisson law of mean `pi`; sum `e^pi`), with its
  peak at `k = 3`; among all integer dimensions the volume peaks at 5;
- `f(n+1)/E[H_n] = (kappa/2) lambda^n/n!` with `lambda = 2 * 5^(1/3) = 3.42`, `kappa in {1, 3/5^(2/3), 9/5^(4/3)}`,
  the ratio of the least strong count at order `n+1` (Busch, CITED in THM-1370) to the mean count `n!/2^(n-1)` at
  order `n`: `3, 3.33, 3, 2, 1.11, 0.57, 0.24` at `n = 2, ..., 8`, with its peak at `n = 3` (with the true maxima
  in place of the means: `3, 5/3, 9/5, 1, 5/9`).

The self-converse share bound `(n-1)!! 2^(-n(n-2)/4)` (`3/4, 15/64, 105/4096` at `n = 4, 6, 8`) is of the same
exponential-over-factorial kind but monotone: it has no peak, it crosses 1/2 (row 11, type threshold).

EXACT: a crevasse of the path spectrum can only sit where the second ratio exceeds 1, and `f(n+1) > max H(n)`
holds exactly for `n <= 4`; this is THM-1370's reason for 7, written as a comparison of an exponential with a
factorial (21 needs more: row 1). ANALOGY: both peaks sit at small arguments because the bases `pi` and `3.42` are
small numbers; the repository rule (MISTAKE-559) applies and no map between the two is claimed. The second half of
the owner's sentence has an exact form in the path spectrum: a value `h` is realised by strong tournaments exactly
at the orders of its window `W(h) = [n(h), e(h)]` (THM-4889 (3)); at the least order it lies in or above the bulk of
that order (many cyclic triangles), at the largest order it lies at the Busch floor (few). The two ends are the two
realisations; the window has 1 to 4 orders for `h < 225`.

**The period-8 phase (IDENTITY) and the Gaussian integral (ANALOGY).** Every closed form of Theorems T and P is a
real or imaginary part of a power of `1 + i`: `sum_{k even} C(n,k) (-1)^(k/2) = Re (1 + i)^n`, and for odd `n` this
is `(2/n) 2^((n-1)/2)`. It is the identity behind the second supplement to quadratic reciprocity (reduce
`(1 + i)^p` modulo `p`: 2 is a square mod `p` iff `p = +-1 (mod 8)`). So the sign `(2/n)` in the pentagon race is the
character whose race among the primes has a Chebyshev bias of the same kind as the race mod 4 (prime squares are
`1 (mod 8)`); in the pentagon race it is exact and has no fluctuation. ANALOGY: the ball volume
`pi^(n/2)/Gamma(n/2 + 1)` is the modulus side of the `n`-th power of a one-dimensional Gaussian integral
(`int exp(-pi x^2) dx = 1`, unfolded in polar coordinates), and the race values are the phase side of the `n`-th
power of a one-dimensional Gauss sum (`1 + i`). The first has a fold, because the unfolding puts a factorial in the
denominator; the second cannot have one (modulus exactly `2^(n/2)`, phase of period 8). The MRB constant is the
alternating sum of the fold `n^(1/n) - 1` (fold note of 2026-11-04); the solvable layers of the race have no
alternating-sum constant of that kind, because they do not fold. No map between these objects is claimed.

**Micro and macro in the race (IDENTITY where stated).**

- The gate `a -> b -> R -> a` preserves `chi` (THM-5127 (C) (7)): a gate tower remembers one bit of its core
  forever.
- For even `n` the deck determines `chi`: `chi(T) = prod_v chi(T - v)` (section 5).
- Ordered sums multiply `chi` (codex signed tables), which gives the strong race of Proposition 5.1.
- A vertex cube carries a quadratic form whose triangle part is `sigma . x + e_2(x)`: the quadratic part `e_2` and
  the Gauss sum `2^ceil(m/2)` are the same for every base (Theorem T), the linear part `sigma` is the score-parity
  vector of the base; the remainder is the parity of the longer even-order paths (THM-5127 (4.4)).
- In the coordinates `y = Ax + |X| 1` the triangle-plus-pentagon part is again a triangle form,
  `C(|y|, 2) + u . y` (Theorem P (2)): the pentagon layer of a vertex cube is the triangle layer seen through the
  adjacency matrix. The step cannot be repeated in the same form (the direct continuation fails at the heptagon).
- For `n <= 6` the macrocosm is the microcosm: `H mod 4` is the weight mod 4 of one vector, `t = A^2 1`.

**What transfers to the conjecture, and what does not.**

- Chan-Pak's exceptions are of crevasse type, not of Borel-Cantelli or square-term type: no residue class mod 4
  is short of realisers (THM-5127 Corollary 4.6: every base has both kinds of extensions), and the race gives no
  obstruction at any order.
- The square-sum problem closed its crevasses by a self-similar induction from a finite base. The repository
  already records why that template does not transfer directly (ATLAS section 4, Session211). The race does not
  change this.
- The squares explain *who is special at small orders*: self-converse and symmetric tournaments. That agrees with
  THM-4922 (first creators are closed under converse) and THM-4930, and it is where a proof by induction on
  orders has to treat small bases separately.

**Principles (recorded here, not promoted).**

- **Squares own the microcosm.** Trigger: an involution on the objects whose fixed points all fall in one class
  of the statistic. Action: count the fixed points against the whole (`(n-1)!! 2^(n^2/4)` against `2^(n^2/2)`
  here) and subtract them before reading a bias. Counterindication: in the prime race the squares also own the
  macrocosm; here they do not.
- **Short cycles carry the bias; the longest cycles erase it.** Trigger: a parity that is a sum over cycle lengths.
  Action: tabulate it by length before modelling it.
- **A weight mod 4 is a permanent bias.** Trigger: a parity of the form `C(|l(T)|, 2)` for a vector `l(T)` whose law
  is uniform on cosets (scores, then second-order scores). Action: its mean is a Gauss sum,
  `Re((1 - i)(1 + i)^a (1 - i)^b)`, of exact size `2^(-n/2)`; it does not decay faster and it does not oscillate
  beyond the period 8 of `(1 + i)^n`.
- **Look for the reason behind a control.** The proof of Theorem P came from a column printed only as a sanity
  check (the number of bases with negative cube sum, which was 0), not from the registered quantities.
- **A coincidence across two orders is not a law until the binomial periods are checked.** Trigger: a quantity
  summed over `k`-subsets. Action: reduce `C(n - j, k - j)` mod 2 first; the form of the identity depends on
  `n mod 2^r`.

## 9. Routes tried in this session and why they stopped

- *A closed form for the generic lead.* Forced coincidences (sources, sinks, gates, twins) explain the second-order
  correlations of `chi` (THM-5127 (C) (5)) and are too small by a factor `2^(-n)` to explain `D_n`. No sign-reversing
  involution with a fixed-point set of size `2^(C(n,2) - 2n)` was found.
- *A proof that the rank of the pair form grows.* Not attempted again; see THM-5127 section 8.
- *The pentagon law.* First form ("invisible at even order") stated, pre-registered, refuted at `n = 10`; second
  form registered, confirmed and then proved (Theorem P).
- *A third exact layer (the heptagon).* The candidate `tr(A^3 + A^5 + A^7) = C(|A^3 1|, 2) + floor(n/4)` is false
  (`tower_test.py`); the increment of `tr A^7` on a vertex cube is quartic in the out-set (a closed 7-walk can pass
  through the new vertex twice), and `c_7` is not a trace. No substitute was found.
- *A pointwise Theorem T for the pentagon.* False: the cube sum of `c_3 + c_5` is not the same for every base
  (`0, 8, 16` at order 5). What is true is its sign and its closed form (Theorem P (3)).
- *A Borel-Cantelli count of missing path counts.* The ATLAS already has the heuristic and the reason it cannot be
  calibrated for correlated families; the strong multiplicities needed are THM-4744 / THM-4810 / THM-4889. Not
  repeated.

## 10. Reproduction

Requires `gentourng` (nauty) and a C compiler. From this directory:

    cc -O2 -o /tmp/scrace2 scrace2.c
    for n in 3 4 5 6 7 8 9; do gentourng -q $n | /tmp/scrace2 $n > out/scrace2_$n.txt; done
    for r in 0 1 2 3; do gentourng -q 10 $r/4 | /tmp/scrace2 10 > /tmp/sc10_$r.out & done; wait
    python3 merge_scrace.py out/scrace2_10.txt /tmp/sc10_*.out      # sum the four parts
    python3 race_tables.py > out/race_tables.out                   # sections 4, 5
    cc -O2 -o /tmp/cycrace cycrace.c
    for n in 3 4 5 6 7 8 9; do gentourng -q $n | /tmp/cycrace $n > out/cyc_$n.txt; done
    for r in 0 1 2 3 4 5; do gentourng -q 10 $r/6 | /tmp/cycrace 10 > /tmp/cyc10_$r.out & done; wait
    python3 merge_cyc.py out/cyc_10.txt /tmp/cyc10_*.out           # sum the six parts (about 50 CPU-minutes)
    python3 cycle_tables.py > out/cycle_tables.out                 # section 6
    python3 c3_bias.py 16 > out/c3_bias.out                        # Theorem T (3) by reduction of the quadratic form
    python3 strong_signed.py > out/strong_signed.out               # Proposition 5.1
    python3 anf5.py 5 > out/anf5.out                               # normal form of c_odd mod 2 on 5 vertices
    cc -O2 -o /tmp/biasmc biasmc.c
    for s in 101 102 103 104 105 106; do /tmp/biasmc 11 2400000 $s > out/mc11_$s.out & done; wait   # Run A
    python3 merge_mc.py out/mc11_*.out > out/mc11_summary.out      # section 7

    for s in 201 202 203 204; do /tmp/biasmc 10 1000000 $s > out/mc10_$s.out & done; wait             # estimator check
    python3 merge_mc.py out/mc10_*.out > out/mc10_control_summary.out
    cc -O2 -o /tmp/c35mc c35mc.c
    for n in 11 12 13 14 15 16 17 18; do /tmp/c35mc $n 40000000 $((1000+n)) > out/c35_$n.out & done; wait   # Run C
    python3 c35_tables.py out/c35_*.out > out/c35_summary.out      # R7, R8 and the comparison with Theorem P
    cc -O2 -o /tmp/k35cube k35cube.c                               # Theorem P: exact cube sums by elimination over F_2
    for m in 1 2 3 4 5 6 7 8 9; do gentourng -q $m | /tmp/k35cube census $m > out/k35_census_$m.txt 2> out/k35_census_$m.err; done
    for r in 0 1 2 3; do gentourng -q 10 $r/4 | /tmp/k35cube census 10 > out/k35_census_10_part$r.txt 2> out/k35_census_10_part$r.err & done; wait
    for m in $(seq 11 40); do /tmp/k35cube mc $m 4000000 $((3000+m)) > out/k35mc_$m.out 2>> out/k35mc_checks.err; done
    python3 k35_tables.py > out/k35_summary.out                    # censuses and samples against the closed forms
    python3 pentagon_closed_forms.py 6 34 130 > out/pentagon_closed_forms.out   # Theorem P (1), (4), (5)
    python3 tower_test.py > out/tower_test.out                     # the direct continuation fails at the heptagon
    python3 kernel_size_check.py 6 > out/kernel_size_check.out     # expected kernel size of A over F_2
    python3 check_selfconverse_mod8.py > out/check_selfconverse_mod8.out   # Theorem S (2) on all classes of orders 4, 6, 8

`scrace.c` and `stronghist.c` are earlier versions kept for the record (two-way split; strong histograms).
`k35cube.c` checks, on every base it reads or draws: the trace formulas against a direct count of 5-cycles (census
mode), the elimination against a direct sum over the `2^m` out-sets (`m <= 9`), Theorem T, the two score formulas,
the closed formula for the cube sum and the weights of the kernel elements; its last line reports the failures.
