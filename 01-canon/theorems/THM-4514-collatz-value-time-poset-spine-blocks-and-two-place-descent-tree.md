---
id: THM-4514
title: "Posets and DAGs of a Collatz orbit: the value-time poset has dimension at most two with leaders minimal and strict future minima maximal, a positive orbit has infinitely many strict future minima iff it diverges (at least one iff its global minimum is attained only once); the excursion forest has the height lower records as roots and, on a divergent orbit, one infinite tree whose unique infinite branch is the spine = the orbit's intersection with the no-descent fractal E_inf; spine blocks (THM-4495's positive min-ending words) are tight: a block of length l >= 2 starts with an odd step, has exactly ceil(l log_3 2) odd steps and height in (0, log_2 3 - 1), and blocks exist only at the lengths 1 + floor(b log_2 3/(log_2 3 - 1)) (density 1 - log_3 2), so each spine step multiplies the value by less than 3/2 plus a carry of at most 3l/4; the descent tree m -> D(m) has D(m) in [m/2, m), on orbits reaching 1 at least floor(log_2 m) + 1 lower records with equality iff m is a power of two, and each first-descent word is an affine bijection from a 2-adic source class onto a 3-adic landing class, giving the in-degree formula (exact modulo Terras's equality) and the limiting sources-per-landing constant c_D in [1.6720, 1.6983]"
status: >
  PROVED (elementary) + FINITE-EXACT + INDEPENDENTLY AUDITED (SOUND WITH
  CORRECTIONS, 2026-09-27, all applied; audit files listed below). Setting:
  T(x) = x/2 (x even), (3x+1)/2 (x odd); parity word w, o_j ones among the
  first j letters, height h_j = o_j log_2 3 - j, carry C_w, T^j n =
  (3^(o_j) n + C_(w,j))/2^j = n 2^(h_j) C_j with C_j >= 1 non-decreasing;
  E_inf = {x in Z_2 : h_j(x) >= 0 for all j >= 1} (S13); sigma(m) = first k
  with T^k m < m, kappa(m) = first k with 3^(o_k) < 2^k (Terras; kappa <=
  sigma); D(m) = T^(sigma(m)) m.
  (P1) Value-time poset (i below j iff i < j and x_i > x_j): dimension <= 2;
  minimal = strict upper records; maximal = strict future minima. Excursion
  order (i ancestor of j iff h_k > h_i on (i, j]) is a forest order with
  roots = strict height lower records (= the descent chain n > D(n) > ...
  exactly when Terras's equality holds along it); a node has infinitely
  many descendants iff it is a strict coefficient future minimum; if h ->
  infinity the forest is locally finite and has exactly one infinite
  branch, the spine (all strict coefficient future minima, rooted at the
  global argmin). For a positive orbit under T or T_q: a strict value future
  minimum exists iff the global minimum is attained only once (divergent,
  or eventually periodic with the minimum in the pre-period; the 5x+1
  orbit of 5 has one), infinitely many iff divergent, none if the orbit
  reaches 1; coefficient future minima are value future minima; value
  lower records are height lower records; the step out of a value future
  minimum is odd.
  (P2) x_j in E_inf iff j is a strict coefficient future minimum; n in E_inf
  implies the orbit is injective, x_j -> infinity, sum 1/x_j < infinity
  (THM-4476), h_j -> infinity, sum 2^(-h_j) < infinity, infinite spine;
  Z^+ ∩ E_inf = {kappa = infinity} is contained in {minima of divergent
  orbits} (sigma = infinity, orbit not eventually periodic), with equality
  iff Terras's equality holds at those minima, and is closed under the
  next-spine-point map s with s(m) > m; eventually periodic positive orbits
  meet E_inf nowhere.
  (T3) Spine blocks (positive words ending at their strict minimum over
  positive times = THM-4495's P_l, Step 1 "minimum decomposition"; their
  reversals are its Step 3 ladder blocks) generate the no-descent words
  freely; sum b_l t^l = 1 - 1/W(t) = 1 - exp(-sum B_n t^n/n). A block of
  length l >= 2 starts with 1 and has 0 < H(w) < log_2 3 - 1, hence exactly
  a(l) = ceil(l log_3 2) ones and height (1 - {l log_3 2}) log_2 3; blocks
  of length l >= 2 exist iff {l log_3 2} > log_3 2 iff l - 1 is not a
  Beatty number floor(a log_2 3) iff l = 1 + floor(b beta), beta = log_2
  3/(log_2 3 - 1) = 2.70951 (Rayleigh); 1^a 0^(l-a) is always such a
  block; along a divergent spine 2^H s_i < s_(i+1) < 2^H s_i (1 + l_i/(2
  s_i)), so s_(i+1) < (3/2) s_i + (3/4) l_i (a one-letter block gives
  exactly (3 s_i + 1)/2).
  (P4) D(m) in [ceil(m/2), m - 1]; if the orbit of m reaches 1, every dyadic
  shell below m contains a lower record, so #records >= floor(log_2 m) + 1
  with equality iff m is a power of two (unconditionally >= floor(log_2 m)
  - floor(log_2 r) + 1 with r the last record); for a first-descent word w
  (length k, o ones, carry C, residue r = -C 3^(-o) mod 2^k, threshold N(w)
  = C/(2^k - 3^o), m_0 = (3^o r + C)/2^k) the map m' -> (3^o m' + C)/2^k is
  a bijection {r + t 2^k > N(w)} -> {m_0 + t 3^o}, and these are the
  descent-tree edges from sources with sigma = kappa; indeg_D(m) = #{w : m
  in the landing class of w, source above threshold} + (Terras-defect
  count, empty for every m <= 10^6); #{m' : D(m') <= X}/X -> c_D = sum_w
  3^(-o(w)) with 3/2 < c_D < 2 and c_D in [1.6720, 1.6983] (audit
  enumeration to |w| <= 28; the session's |w| <= 26 bracket is [1.6695,
  1.7005]).
  (P5) THM-4503's cell (a, b), a >= m_b, is the interval [empty, g] of
  Young's lattice in the partition coordinates e_i = zeros before the i-th
  one, and the carry C_w is a strict order embedding of it into Z.
  FINITE-EXACT: descent tree to 2*10^6 (landing range, Terras equality,
  record bound with equality set = the 20 powers of two, in-degrees mean
  1.6903 max 17 at 293501), two-place in-degree formula exact for all
  landings <= 10^5 at sigma <= 26 (session) and for all landings <= 20000
  with sources of every sigma up to 135 (audit; N(w) < 2^|w| for all 190069
  words of length <= 26, the only least residue at or below its threshold
  being r = 1 for w = 10), block enumeration to length 16 and both
  generating-function forms to length 60, cells to a + b <= 12, the
  Sturmian element of E_inf (spine {0}, infinite non-locally-finite tree,
  divergent real series), the 5x+1 orbit of 7 (512 spine points in 3000
  steps; 165 of the first 300 spine points acquire a Higman-good partner
  within the window, after at most 141 blocks; whole-orbit landing
  multiplicity O(1) at every depth over 20000 steps).
  Collatz OPEN. Nothing here is a rank, a wqo, or a transversality
  statement; section 6 of the note types every DAG finishing move and its
  obstruction.
source: opus session collatz-poset-dag-20260927
depends_on:
  - 01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md
  - 01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md
related:
  - 01-canon/theorems/THM-4503-collatz-poset-bridge-height-selection.md
  - 01-canon/theorems/THM-4506-landing-multiplicity-exact-worst-case-and-recursion-saturation.md
  - 01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md
  - 05-knowledge/results/collatz_posets_dags_20260927_spine_descent_tree.md
  - 05-knowledge/results/collatz_directions_20260926.md
  - 05-knowledge/results/collatz_crossings_20260926_potential_and_seeds.md
  - 05-knowledge/results/forest_20260926_excursions.md
scripts:
  - 04-computation/experiments/collatz_posets_dags_20260927.py
  - 04-computation/experiments/collatz_posets_dags_20260927_chains.py
  - 04-computation/experiments/collatz_posets_dags_20260927_audit.py
outputs:
  - 05-knowledge/results/collatz_posets_dags_20260927.out
  - 05-knowledge/results/collatz_posets_dags_20260927_chains.out
  - 05-knowledge/results/collatz_posets_dags_20260927_audit.out
audit: "Independent auditor subagent, 2026-09-27: blind re-derivation from the statements plus an independent exact-integer script (84 checks); verdict SOUND WITH CORRECTIONS (C1-C12, all applied: the false 'strict future minimum iff divergent' direction, the reach-1 hypothesis of the record bound, the c_D lower end, attributions); report 05-knowledge/results/collatz_posets_dags_20260927_audit.md; session note section 9."
---

# THM-4514 — Posets and DAGs of a Collatz orbit

Full statements, proofs, numerical faces and the typed finishing-move
shapes are in the session note
[`collatz_posets_dags_20260927_spine_descent_tree.md`](../../05-knowledge/results/collatz_posets_dags_20260927_spine_descent_tree.md),
sections 2–6, with the audit in its section 9. The one genuinely new exact
object is Theorem 3 (tight spine blocks, the Beatty law of block lengths);
Propositions 1, 2, 4, 5 are elementary dictionaries that make the excursion
forest, the descent tree and THM-4503's cells exact objects with both
places visible. The statement "the spine of a divergent orbit is the
orbit's intersection with `E_inf`" reduces the board's "rank across
completed excursions" to a function on spine points, which by S13's
Proposition 4 cannot be a function of the word prefix; Theorem 3 sharpens
the witness (the prefix at a spine point is a concatenation of tight blocks
realized by a whole residue class). The two conditional clauses (equality
in (P2), vanishing of the defect term in (P4)) are exactly Terras's
equality `kappa = sigma`, verified to `10^7` (THM-4512) and here to
`2 * 10^6`, and are marked as such.
