---
id: THM-4470
title: "The pairing ladder: 3x+1 is the unique arithmetic-mean-fair consecutive pairing map, its graph is two perfect difference systems, and no pairing statistic decides tree-ness"
status: >
  PROVED + INDEPENDENTLY AUDITED (items 1-5); FINITE-EXACT (item 6).
  (1) The shortcut map T preserves the sum of every pair {2i-1, 2i}. Among
  the maps n/2 (n even), (qn+r)/2 (n odd), with q, r odd, only 3n+1 on
  {2i-1, 2i} and 3n-1 on {2i, 2i+1} preserve the sums of a consecutive
  pairing, and only q = 3 admits any sum-preserving perfect matching.
  (2) The halving edges {i, 2i} and the up edges {2i-1, 3i-1} each realise
  every difference 1, 2, 3, ... exactly once, but no Collatz subtree with
  m >= 2 edges is graceful under its own labels.
  (3) In the pairing family (a bit per pair chooses which member goes up),
  every member has properties (1) and (2). Collatz and 3n-1 are antipodal
  corners, the all-zero word and the all-one word shifted by one.
  (4) Flipping a density-zero set of pairs of Collatz yields a divergent
  orbit.
  (5) No periodic pairing admits a bounded-lookahead descent certificate.
  (6) Flipping a single pair i <= 10^7 creates a new cycle for exactly 24
  values of i, all at most 2308.
source: collatz-procgen-20260922 session, brackets/pairings lane (2026-09-24), from the owner's "graceful tree : 3N+1" analogy and the coordinator's pair-sum targets; audited and promoted by the session orchestrator 2026-09-24
depends_on: []
related:
  - 05-knowledge/results/collatz_procgen_20260922_choice_ladder.md
  - 05-knowledge/results/collatz_procgen_20260924_inverse_tree_mod192.md
  - 05-knowledge/hypotheses/HYP-9135-finitely-many-fragile-pairs.md
  - 05-knowledge/hypotheses/HYP-9136-provable-pairing-price-tends-to-zero.md
script: 04-computation/experiments/procgen_brackets_20260924_pairings.py
script_audit: 04-computation/experiments/procgen_brackets_20260924_orchestrator_check.py
output: 05-knowledge/results/procgen_brackets_20260924.out
output_audit: 05-knowledge/results/procgen_brackets_20260924_orchestrator_check.out
script_sha256: 92f8c850fc84d1c94b1350da4d86c8b32f4ba2d48ddf3678c1ad65be0bc70c42
script_audit_sha256: 3c7136e435640267f989de441191ca0b525e09a1c0db98706791ce39d372240c
output_sha256: 89951c7953634c77ea3a9c4dbb10c2f5b7d95af0a86512b002426c2d0dcd347a
output_audit_sha256: fbbfd4812b95b5c1443de04b1c434472598e28a0a29ed8b43367a34cf876b9be
hash_basis: raw LF bytes
audit: >
  The orchestrator checked each proof by hand: the q-uniqueness via
  displacements, the literal-graceful bound m <= 2, the pair-0 obstruction
  and the falsifying chain.
  Independent code (procgen_brackets_20260924_orchestrator_check.py, written
  without reading the lane's scripts) reproduces:
  * the pair-sum identity to i = 10^6;
  * the (q, r, offset) census over odd q <= 39, |r| <= 39;
  * the exact lists of 24 fragile pairs of T and 12 of 3n-1, to i = 20000;
  * the bracket-internal move lists and the complete (k, p) list;
  * the x49 square-sum partition of [25, 1249];
  * the falsifying chain (169 terms to 10^30, 82 flips).
  All five lane scripts were re-run. Their output is identical to the
  committed .out except for timing lines, the CP-SAT designs included.
  Peak memory was 594 MB.
---

# THM-4470 -- the pairing ladder

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_brackets_20260924_pairings_transitions](../../05-knowledge/results/procgen_brackets_20260924_pairings_transitions.md) §§0, 2.

`T(n) = n/2` (n even), `(3n+1)/2` (n odd), and `|T(n) - n| = ceil(n/2)`.

## 1. Arithmetic-mean fairness

* `T(2i-1) + T(2i) = (3i-1) + i = 4i-1 = (2i-1) + 2i`.
* `T(2i-1) T(2i) / ((2i-1) 2i) = 3/4 + 1/(8i-4)`.
* Each pair keeps its arithmetic mean while its geometric mean shrinks
  toward `sqrt3/2` of its value. This is the exact pairwise form of the
  drift `log(2/sqrt3)`, the AM–GM gap of the step factors `3/2` and `1/2`.
* `{1,2}` is the only pair mapped onto itself (the trivial cycle).
* `sum_(n <= 2M) T(n) = sum_(n <= 2M) n`.

**Uniqueness.** For `f(n) = n/2` (even), `(qn+r)/2` (odd):
* `f(2i-1) + f(2i) = (q+1) i + (r-q)/2`, which equals `4i-1` for all `i` iff `(q,r) = (3,1)`;
* `f(2i) + f(2i+1) = (q+1) i + (q+r)/2`, which equals `4i+1` for all `i` iff `(q,r) = (3,-1)`.

So the two sheets are exactly the two consecutive pairings. More generally,
the displacement `f(n) - n` is `-n/2 < 0` on evens and `((q-2)n + r)/2` on
odds. So a sum-preserving pair must be `{odd n, even (q-2)n + r}`. For `q = 1`
the partner is not positive. For `q >= 5` the evens `m != r (mod q-2)` are
never partners, so no sum-preserving perfect matching exists. A zero-sum partition
into blocks of bounded diameter `D` is also impossible for `q != 3`: the
displacement over `[1, 2M]` is `(q-3)M^2/2`, against `O(DM)`.

## 2. The graceful form

* The down-edges `{i, 2i}` and the up-edges `{2i-1, 3i-1}` each use every
  difference `1, 2, 3, ...` exactly once, and they share only `{1,2}`.
* Collatz is the statement that their union is connected, and then it is
  a spanning tree of `Z_(>0)`.
* The halving forest on `[1, 2n]` is Skolem graceful (Lee–Shee).
* **Literal gracefulness fails.** A subtree with `m` edges labelled by its
  own integers inside a window of length `m` needs difference 1, which is
  only the edge `{1,2}`, and difference `m`, which needs a vertex `>= 2m-1`.
  Hence `m <= 2`.
* **Near-gracefulness is bounded.** Labels in a window of length `m+k` with
  distinct differences force `m <= 3k+3`.

## 3. The pairing family and its corners

* For bits `eps_i`, `n` in pair `i` goes up by its length iff (`n` odd) XOR `eps_i`. Every member satisfies §1's sum identity and §2's difference systems.
* Shifting by one and complementing every bit exchanges the two offsets: `F^0_eps(n) + 1 = F^1_(1-eps)(n+1)`. Hence `3n-1` is the all-flipped Collatz pairing, shifted by one.

## 4. Density-zero flips falsify

* Put `c_0 = 3` and `c_(t+1) = c_t + ceil(c_t/2)`, and flip the pair of every even `c_t`.
* The chain's odd members never share a flipped pair, since consecutive terms differ by at least 2.
* So the modified map sends `c_t -> c_(t+1)` for every `t`, a divergent orbit, and the flip set has at most `log_(3/2) X + 1` elements below `X`. ∎
* **Consequence.** No property of the pairing that is insensitive to density-zero changes can decide the conjecture: pair sums, difference systems, flip density, or any residue or density statistic. This is the DEFECT control, made concrete inside the one family where the owner's balance and gracefulness hold exactly.

## 5. No periodic pairing is provable by bounded lookahead

* Let `eps_i = c(i mod 2^K)`. Near the 2-adic point `-1` (offset 0), the pair index is `0 mod 2^(j-1)`, so the bit is `c(0)`.
* If `c(0) = 0`, the odd member rises to `3i - 1`, again near `-1` with one digit fewer. If `c(0) = 1`, the even member rises to `3i`, near `0`.
* Either way one 2-adic neighbourhood rises for about `j` steps, at every level `j`. So no fixed lookahead `L` certifies descent. This is Applegate–Lagarias's "`-1` resists elimination", forced by the pair of index 0. ∎

## 6. Fragility (FINITE-EXACT)

Flipping a single pair `i` of Collatz creates a new cycle exactly for

`i = 1, 4, 5, 10, 11, 13, 20, 22, 40, 61, 84, 122, 126, 167, 189, 217, 244, 325, 334, 433, 445, 577, 1154, 2308`,

and for no other `i <= 10^7`.
* The mechanism is the cycle gate: a new cycle is a path `3i ~> 2i` (or `i-1 ~> 2i-1`) of `T`, i.e. `i (2^(K+1) - 3^(a+1)) = B_w`.
* The cycles sit at upper approximations of `log_2 3`: `(a,K) = (3,5), (5,8), (10,16), (17,27), (29,46), (34,54), (46,73)`.
* For `3n-1` there are 12 fragile pairs up to `10^7`. The last is `i = 12029`, at the convergent `84/53`.
* Whether the list for `T` is finite is HYP-9135.
