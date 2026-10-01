---
id: THM-4522
title: "The multiplicative lonely runner, where LRC (local) meets Collatz (global). Every box of 3-smooth speeds {2^j 3^k} is lonely at 1/5 because 5 divides no speed, so LRC holds there trivially. The x2x3 lonely spectrum I(t) = inf ||2^j 3^k t|| has discrete top {1/5, 1/7, 1/10, 1/11, 1/13, 1/14}, reached by exactly 88 rationals. Near runners at a time h/D form disjoint 'triangles' with 6-free apexes (a coprime-to-6 Tao triangle lemma). No crowded-time dichotomy can exclude a Collatz cycle"
status: >
  PROVED + INDEPENDENTLY AUDITED (Proposition K; Theorem S, computer-
  assisted in exact rational arithmetic; Corollary S'; Theorem L; Theorem T,
  the triangle lemma; Proposition C; Proposition B, the barrier);
  FINITE-EXACT (exception lists for the discrete multiplicative LRC; the gate
  spectra for p <= 40); EMPIRICAL (Conjecture D on gate sums; the logarithmic
  crowded-time law); CITED (Furstenberg 1967; Tao 2019, Lemma 7.4;
  Perarnau-Serra survey arXiv:2409.20160; BLMV via Z. Wang
  arXiv:1004.0035); REFUTED (a bounded-numerator dichotomy for runner
  boxes).
  Setting: boxes B(J,K) = {2^j 3^k : j < J, k < K}; kappa(V) =
  max_t min_(v in V) ||v t||; I(t) = inf_(j,k >= 0) ||2^j 3^k t||.
  (K) kappa(B(J,K)) is 1/2 (J = 1), 1/3 (K = 1), 1/4 (J = 2, K >= 2) and
  1/5 (J >= 3, K >= 2). More generally kappa(V) = 1/5 for every finite set
  of 3-smooth speeds containing a dilate of {1,2,3,4}. So LRC holds for
  multiplicative speeds with unbounded margin (n+1)/5. The mechanism is a
  mod-5 obstruction; lacunarity is irrelevant.
  (S) I(t) > 0 forces t rational (Furstenberg). Moreover:
  - I(t) <= 1/5, with equality iff t is in {1,2,3,4}/5;
  - the values of I in [1/14, 1/2] are exactly 1/5, 1/7, 1/10, 1/11, 1/13
    and 1/14;
  - exactly 88 points of (0,1) have I >= 1/14, with denominators
    5, 7, 10, 11, 13, 14, 26, 28, 33, 52, 56.
  (S') For D coprime to 6 and boxes containing B(16,10) and every
  3-smooth number <= 219 D: the best time h/D has every runner >= 1/14
  iff 5, 7, 11 or 13 divides D; it is then 1/q for the least such q.
  (L) Discrete LRC by 6-free counting: L(D,J,K) >= 1/(JK+1) for all D
  coprime to 6 beyond explicit thresholds, with complete exception lists
  below them. For the box (3,2) the exceptions are D = 11, 17, 37.
  (T) For D coprime to 6 and R <= D/4, the runners near 0 at a time h/D
  form disjoint triangles with exact residues, apex numerators coprime to
  6, and disjoint shadows. The threshold D/4 is sharp.
  (B) Barrier. Write N for the number of words with c_w = 0 mod |G| at the
  gate G = 2^p - 3^a. If |S(h)| <= eps C off a structured set, Parseval
  forces minor-arc L1 mass >= lambda (1 - C/M)/eps. So bounding the
  minor-arc terms one at a time can never prove N = 0, however sharp the
  crowded-time dichotomy. Excluding a cycle needs signed cancellation,
  which is N = 0 itself.
  Reading:
  - The owner's local (LRC) / global (Collatz) duality, made concrete. For
    multiplicative speeds the local question is trivial (a mod-5
    obstruction). The global question has a discrete, prime-indexed
    spectrum (5, 7, 11, 13) and a rigid near-structure (triangles). It
    still cannot be closed by size bounds.
  - Both lanes rank the hybrid LOW as a proof route and MEDIUM as a
    diagnostic.
  Collatz and LRC are OPEN; see the external claim on LRC(14) in
  OPEN-QUESTIONS.
source: collatz-procgen-20260922 session, mlr lane (2026-09-30), answering the owner's thesis that LRC (local) and Collatz (global) are complementary and tuned to the primes; audited and promoted by the session orchestrator 2026-10-01
depends_on:
  - 05-knowledge/results/procgen_gates_20260925_gate_equidistribution.md (the gate sums S(h) and their structured spectrum)
related:
  - 05-knowledge/results/procgen_localglobal_20260930_lrc_collatz_primes.md (the companion lane: prime alignment, the duality, kappa(2,3) = 1/5)
  - 05-knowledge/results/procgen_bridges_20260923_lrc_amm_collatz.md (the earlier LRC-Collatz verdict: ANALOGY)
  - 01-canon/theorems/THM-4484-free-and-sporadic-cycles.md
note: 05-knowledge/results/procgen_mlr_20260930_multiplicative_lonely_runner.md
scripts: 04-computation/experiments/procgen_mlr_20260930_{run,core,lrc,crowd,gates}.py
script_audit: 04-computation/experiments/procgen_mlr_20260930_orchestrator_check.py
output: 05-knowledge/results/procgen_mlr_20260930.out
output_sha256: 7e29e9ee562f3b4961da0a0eb80c4b31f59ceb6504eb2bcbee5126216488c0d1
output_audit: 05-knowledge/results/procgen_mlr_20260930_orchestrator_check.out
output_audit_sha256: f79c593731e246ff82174f1b385e301b47e6fadd63d5b8524a24cdfd260241bc
hash_basis: raw bytes
audit: >
  The orchestrator read the proofs of Proposition K, Theorem S (the exact
  covering of X_delta by the box B(16,10) at five levels, then the local
  stabiliser lemma near each rational centre) and Proposition B, and read
  the statements of Theorems L and T.
  Independent code (procgen_mlr_20260930_orchestrator_check.py, the lane's
  code not read) confirms:
  - Proposition K exactly for all J <= 5, K <= 3;
  - the census of Theorem S over all rationals r/N with N <= 700: exactly
    88 points with I >= 1/14, the stated denominators and the six values;
  - the exception list {11, 17, 37} for the box (3,2) over all D coprime
    to 6 up to 3000.
  The lane's runner was re-run (269 s, ALL CHECKS PASSED). It is identical
  up to timing fields.
  Scope notes:
  - Theorem S's completeness beyond the census rests on the lane's
    computer-assisted covering, which was read but not re-implemented.
  - BLMV is not on arXiv and was quoted via Wang's paper.
---

# THM-4522 — the multiplicative lonely runner

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_mlr_20260930_multiplicative_lonely_runner](../../05-knowledge/results/procgen_mlr_20260930_multiplicative_lonely_runner.md).

## 1. Why this object

The owner sees the Lonely Runner Conjecture as *local* and Collatz as *global*:
- LRC asks for one good time for each speed set;
- Collatz must exclude every large cycle.

The Collatz cycle count at a gate `2^p - 3^a` is a sum over the "multiplicative runners" `2^j 3^k` (the gates lane). So the runners whose speeds are 3-smooth numbers carry both questions at once:
- the LRC question: is there a time at which all are far from the start?
- the Collatz question: at which times are they crowded near the start?

## 2. What holds

**LRC for multiplicative speeds is trivial.** No 3-smooth number is divisible by 5, so at `t = 1/5` every runner is at distance at least `1/5`, and nothing does better. The lonely spectrum of the whole `×2×3` semigroup is discrete at the top: `1/5, 1/7, 1/10, 1/11, 1/13, 1/14`, indexed by the small primes 5, 7, 11, 13 that are coprime to 6.

**Crowding is rigid.** Near a time `h/D`, the runners close to the start group into triangles whose apexes have numerators coprime to 6, Tao's 2019 triangle picture for moduli coprime to 6.

**But crowding cannot exclude cycles.** Parseval forces the minor-arc mass to be large whenever the crowded set is small. A cycle bound obtained by summing absolute values of the Fourier terms therefore fails by about `sqrt(C)`, however sharp the description of the crowded times.

## 3. The answer to the thesis

- The local and global questions do meet here: one product set, one prime-indexed spectrum, one triangle structure.
- What the local side knows (a mod-5 obstruction) is useless to the global side.
- What the global side needs (signed cancellation in the word sum) is invisible to every local test.
- This is the owner's "each lacks information about the other's type", now as a pair of theorems: Theorem S on one side and Proposition B on the other.
