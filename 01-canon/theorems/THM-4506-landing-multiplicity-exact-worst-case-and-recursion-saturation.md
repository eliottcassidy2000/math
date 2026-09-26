---
id: THM-4506
title: "Landing multiplicity in the thin-divergence recursion. All dippers of a landing point lie in one dyadic shell, with an odd letter between consecutive ones, so the worst case is exactly ceil((k-D)/log_2 3), about 0.631(k-D), not k. The recursion is saturated at a*(mu) = mu lambda*/h* - 3/2. Averaging over depths, orbit-blind residue-class splits and 2-adic data all leave mu = 1, so only an orbit-coupled local-time bound can lower THM-4499's exponent"
status: >
  PROVED + INDEPENDENTLY AUDITED (Lemma S, Lemma O, Theorem 2,
  Theorem 1(i) as a conditional implication, Theorem 1(ii), Proposition C,
  Lemma A, Proposition H); FINITE-EXACT (exhaustive worst case to k = 25;
  constructions to k = 400; record orbits to 10^12); CITED (Baker's
  theorem for the O(1) sharpness of the lower bound at every large k;
  OEIS A006877 and A006884 for the record lists).
  Setting: T_b(x) = x/2 or (3x+b)/2, b odd. A segment of distinct
  positive integers, scale X, k = floor(log_2 X), depth D >= 1. A dipper
  is an index i with y_i <= X and y_(i+s) < 2^(-D) y_i for some
  1 <= s <= k; its landing index j is the least such i+s. m(j) is the
  number of dippers landing at j.
  (S) If b > 0, or b < 0 and y_j >= |b|: the step into j is a halving,
  and every dipper of j lies in the single shell (2^D y_j, 2^(D+1) y_j].
  (O) Consecutive dippers of j are separated by at least one odd letter.
  (2) Hence m(j) <= ceil((k-D)/log_2 3) for b > 0, and
  m(j) <= ceil((k-D+log_2(3/2))/log_2 3) for b < 0 when y_j >= k|b|
  (smaller y_j contribute O(k^2 |b|) dippers in total). This sharpens the
  one-bit-band bound 2 ceil(k/3) of the parallel S7 note. The worst case
  over all integer segments equals the bound in every exhaustively
  scanned cell (b = 1, -1, 5, k <= 25). The actual orbit of 13255 attains
  12 = ceil(18/log_2 3) at L = 20, D = 2.
  (1) If an orbit's averaged multiplicity satisfies
  #Dip(X,D) <= C L^mu N(X 2^(-D)) + C L^2, then N(X) <= K X^(h*) L^a for
  every a > a*(mu) = mu lambda*/h* - 3/2; for mu = 0, also a = -3/2.
  Conversely, psi = X^(h*) L^(a*(mu)) satisfies the recursion inequality
  at every depth. So a*(mu) is exactly what the one-window method yields.
  Values: a*(1) = -0.9862 (THM-4499), a*(1/2) = -1.2431, a*(0) = -3/2.
  (C) An odd run gives at most 2 dippers to one landing point. A W-bit
  hover spreads its dippers over at most ceil(W)+1 landing points.
  (A) A pair (i, j) serves at most one integer depth; averaging over
  depths improves the effective multiplicity by a factor of at most 0.966.
  (H) For every M >= 2 there is a residue class mod 2^l, with
  l <= 2 log_2 3 (M-1) + D + 2, whose large members are all heavy dippers.
  So any split "light landing points <= M, heavy dippers counted as
  integers" needs M = Theta(L). A 2-adic word can be hostile at every
  scale. No orbit-blind or parity-word-only argument gives mu < 1.
  OPEN: the orbit-coupled averaged multiplicity (HYP-9161 in its
  endpoint-corrected form), a local-time statement about one orbit.
  Collatz is OPEN. Nothing here bears on the existence of divergent
  orbits beyond THM-4499.
source: collatz-procgen-20260922 session, landing lane (retry; 2026-09-26), answering the owner's "attack landing multiplicity"; audited and promoted by the session orchestrator 2026-09-26
depends_on:
  - 01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md
  - 01-canon/theorems/THM-4499-thin-divergence-is-little-o-of-x-to-the-h-star.md
  - 01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md (the ballot floor -3/2)
  - 01-canon/theorems/THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md (lambda*)
related:
  - 05-knowledge/hypotheses/HYP-9161-landing-multiplicity-polylog.md
  - 05-knowledge/results/collatz_landing_20260926_multiplicity_reassessment.md (parallel S7 note: the formula a*(beta), the one-bit band lemma 2 ceil(k/3))
  - 05-knowledge/results/crossroads_poset_20260926_integer.md (refutation of the literal HYP-9161)
  - 01-canon/theorems/THM-4478-critical-tube-affine-capacity-sharp-provability-price.md (at most floor(k/3)+1 integer ancestors per selected word: the inverse-tree counterpart of Lemma O)
note: 05-knowledge/results/procgen_landing_20260926_landing_multiplicity.md
scripts: 04-computation/experiments/procgen_landing_20260926_{lib,run}.py and procgen_landing_20260926_scan.c
script_audit: 04-computation/experiments/procgen_landing_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_landing_20260926.out
output_sha256: 31bb595b698ab29bcb4f504d139c7ef48bef9b09ad89bf7043c5adaea7962723
output_audit: 05-knowledge/results/procgen_landing_20260926_orchestrator_check.out
output_audit_sha256: 9feb829f16d68da59967db25c307255654662b519725c502ef1438205cde2d1e
hash_basis: raw bytes
audit: >
  The orchestrator read the proofs of Lemma S (the first-dip minimality
  forces y_(j-1) > y_j, hence a halving), Lemma O, Theorem 2 (ratio
  bound 3^O < 2^(n-D), with the Weierstrass factor 2/3 for b < 0),
  Theorem 1(ii) (the two cases t <= 1 and t > 1), Lemma A and
  Proposition H, and found them sound. Theorem 1(i) repeats THM-4499's
  induction with k replaced by C L^mu.
  Independent code (procgen_landing_20260926_orchestrator_check.py,
  written from the note's statements without reading the lane's scripts)
  confirms:
  - Lemmas S and O and the Theorem 2 bound on every first dipper
    y <= 2^(k+1)-1 for b = 1, -1, 5, k = 10, 12, 14 and D = 1, 1.5, 2, 3,
    with exact comparisons. For b = 1 the maximum equals the bound in all
    12 cells.
  - The orbit of 13255 has a landing point with 12 dippers at L = 20,
    D = 2.
  - A Proposition H class (b = 1, D = 3, O = 10, l = 19): in 40 random
    members every visit is a dipper landing at index l.
  - Theorem 1(ii) on a grid (mu in [0,1], C in [0.01, 100],
    L = 2^4..2^17, all D in [1, L/2]): min log_2(RHS/psi) = 0.925.
  - Lemma A on 60684 pairs along six record orbits.
  The lane's full pipeline was re-run (147 s; 130 checks, ALL CHECKS
  PASSED). Its output is identical, after normalizing timing and memory
  lines, to the lane's final run log, whose raw bytes are the committed
  .out (sha256 above).
  Scope notes:
  - The worst case is over segments of distinct positive integers. The
    hostile integers and records reach 1.
  - The sharp lower bound (k-D)/log_2 3 - O(1) at every large k rests on
    Baker (CITED). The proved lower bound is Proposition H's
    floor(O/2)+1 visits.
  - The S7 note's climb-then-drop claim was already refuted there. This
    theorem supplies the exact constant.
---

# THM-4506 — landing multiplicity: one shell, one odd letter, and the saturation of the thin-divergence recursion

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_landing_20260926_landing_multiplicity](../../05-knowledge/results/procgen_landing_20260926_landing_multiplicity.md).

## 1. The exact worst case

THM-4476's recursion charges each dipper to its landing point and bounds the multiplicity by `k`. Two facts fix the true constant:
- **One shell.** Every dipper of a landing point `y_j` lies in the single dyadic shell `(2^D y_j, 2^(D+1) y_j]`. The step into `j` must be a halving, because the orbit was still above the threshold one step earlier.
- **One odd letter each.** Between two dippers there is an odd step. Otherwise the second would have fallen out of the shell by halving.

Combined with the size ratio `3^O < 2^(n-D)`, these give `m(j) <= ceil((k-D)/log_2 3) ≈ 0.631(k - D)`. The bound is attained:
- exhaustively for every scanned `(b, k, D)`;
- by Sturmian hover-then-halve integers up to `k = 400`;
- by actual Collatz orbits, e.g. 13255 at `L = 20`.

This uses exactly two pieces of information: the real place (the size ratio) and the 2-adic place (the letters). Sharpness shows that no further two-place information improves the worst case.

## 2. What the multiplicity is worth, and why the obvious averages fail

The recursion only needs the averaged multiplicity. With an average of `L^mu`, the exponent is `a*(mu) = mu lambda*/h* - 3/2`, and the recursion cannot do better (saturation). The whole prize is the factor `(log X)^0.514`.

Three natural averages all leave `mu = 1`:
- **Over depths.** At most a factor of 0.966.
- **Over integers by residue classes.** Whole classes consist of heavy dippers, so the threshold must be `Theta(L)`.
- **Over the parity word.** A 2-adic word can be hostile at every scale.

What remains is the orbit-coupled question: how often can one orbit revisit the shell lying `D` bits above its later landing points? This is HYP-9161 in its corrected form, a local-time statement about a single orbit.
