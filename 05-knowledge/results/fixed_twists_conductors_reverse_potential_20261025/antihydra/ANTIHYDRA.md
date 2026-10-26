# Antihydra's odd orbit as a 3-twisted Syracuse map: shift conjugacy, a rational-slope valuation barrier, and the golden-ratio model

2026-10-25, klein. **PROVED** structural statements (A1–A4, A6); **DERIVED under an explicit random
model** (A5); **FINITE-EXACT** orbit data (`antihydra_sim.py`, `antihydra_walk.py`; receipts). The
rule below is the public description of the BB(6) cryptid Antihydra; the repository's own cryptid
transcription (`04-computation/bb6_cryptids_s56/`) contains the eleven Dos Santos et al. systems
A1–A11, not Antihydra, so the bookkeeping constants are taken from the public description and the
structural statements are insensitive to them. **No theorem ID is taken:** opus-2026-10-25-S56 has
reserved THM-4770/4771 and HYP-9380/9381 for odd-orbit Antihydra results; this note was sent to opus.

**Rule.** `a_0 = 8`, `b_0 = 0`; if `a` is even: `(a, b) <- (3a/2, b+2)`; if `a` is odd: halt when
`b = 0`, else `(a, b) <- ((3a-1)/2, b-1)`. In both branches `a <- floor(3a/2)`.

## A1. The map is the negation mirror of the Mahler carry map (THM-2228)

`floor(3a/2) = -ceil(3(-a)/2)`, and `ceil(3a/2) = (3a + (a mod 2))/2 = T(a)` is exactly the map of
[THM-2228](../../../01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-strict-real-tails-and-ordinary-stabilization.md),
whose binary carry words are homeomorphic to `Z_2` by a map `Phi` with `Phi ∘ shift = T ∘ Phi`.
Hence `f(a) = floor(3a/2)` extended to `Z_2` is conjugate to the full 2-shift by `nu ∘ Phi`,
`nu(a) = -a`. The same negation conjugates the `3n+1` and `3n-1` sheets
([THM-4474 (E)](../../../01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md),
THM-4610's `F_p(x) = -C_p(-x)`). Antihydra's state space is therefore the Mahler `3/2` carrier read
on the negative sheet.

## A2. The odd first-return map

For odd `a` let `v = v_2(3a-1) >= 1`. Then `(3a-1)/2` has valuation `v-1`, the next `v-1` steps are
even steps multiplying by `3/2`, and the next odd term is

    A(a) = 3^{v-1} (3a-1) / 2^v = 3^{v-1} · S^-(a),      S^-(a) = (3a-1)/2^{v_2(3a-1)} = -S(-a),

with `S` the Syracuse map. So Antihydra's odd orbit is the twisted (`3n-1`) Syracuse map with the
discarded powers of two replaced by powers of three. Example: `27 -> 135 -> 303 -> 681 -> 1021`
(`v = 4, 2, 2, 1`). `A` is strictly increasing on odd `a >= 3`, so it has no cycles; the `3n-1`
Syracuse map has the cycles `{1}, {5,7}, {17,...}` (THM-4484, THM-4610).

## A3. Valuation sequences are a full shift, and cylinders are 2-adic balls

**Theorem.** `Phi_A : a -> (v_1, v_2, ...)`, `v_{k+1} = v_2(3 A^k(a) - 1)`, is a homeomorphism from
the odd 2-adic integers onto `(Z_{>=1})^N` (product topology) with `Phi_A ∘ A = shift ∘ Phi_A`. The
set of `a` with prescribed `(v_1,...,v_k)` is one residue class modulo `2^{1 + v_1 + ... + v_k}`.
*Proof.* `v_2(3a-1) = v` fixes `a` modulo `2^{v+1}`. Writing `a = a_0 + 2^{v+1} s`, `s in Z_2`,
`A(a) = A(a_0) + 2·3^v s`; as `s` runs over `Z_2`, `A(a)` runs over all odd 2-adic integers (`3^v` is a
unit), bijectively and isometrically up to the factor `2`. Induct: `(v_1,...,v_k)` determines `a`
modulo `2^{1+v_1+...+v_k}` and every sequence occurs; a bijective continuous map of compact
Hausdorff spaces is a homeomorphism. QED. The same proof gives the classical statement for `S` and
`S^-` (the repo has the cylinder formula for `S` in
[collatz_renewal_source R1](../../collatz_renewal_source_20261011.md)).

**LCA reading** (merge with [Session56's valuation-LCA carrier](../../session56_rmq_lca_connections_20261025/valuation-lca/REPORT.md)):
two odd 2-adic integers share the first `k` valuations iff their lowest-digit-first trie has
least-common-ancestor depth at least `1 + v_1 + ... + v_k`. Valuation partial sums are LCA depths.

## A4. Halting is a rational-slope ballot barrier (PROVED given the rule)

Let `b_k` be the counter after the `k`-th odd term (`b_1 = 5` at `a = 27`). Between consecutive odd
terms there are `v - 1` even steps, so `b_{k+1} = b_k + 2(v_k - 1) - 1` and

    b_k = 5 + W_{k-1},     W_k = sum_{i<=k} (2 v_i - 3).

Antihydra halts iff some `b_k` would become `-1`, i.e. iff `W_k = -6` for some `k`, i.e.

    sum_{i<=k} v_i <= (3k - 6)/2   for some k.                                        (barrier)

In trie language: the `k`-cylinder of `27` has LCA depth at most `(3k-4)/2`. Compare Collatz: an orbit
of the Syracuse map stays bounded only if `sum_{i<=k} v_i > k·log_2 3 + O(1)` eventually (THM-4499
Lemma M's moving barrier; the coefficient-descent codes `2^{A_j} < 3^j` of renewal R1). Both are
statements about the same i.i.d.-looking valuation sequence staying above a line: slope `3/2`,
exact and prefix-wise, with offset `-3`, for Antihydra; slope `log_2 3`, irrational and asymptotic,
for Collatz. The drift `E[v] = 2` exceeds both slopes; the size drift differs in sign
(`log a` grows by `v log(3/2) > 0` per odd step for Antihydra, by `log 3 - v log 2`, mean negative,
for Syracuse). This is the precise content of "differing in the barrier condition and the drift
direction": one barrier is a Diophantine approximation of `log_2 3` (THM-4633), the other is a
lattice-path problem with rational slope and an exact generating function.

## A5. The golden-ratio halting probability under the Bernoulli model (heuristic)

Under odd 2-adic Haar measure the valuations are i.i.d. with `P(v = j) = 2^{-j}` (renewal R1), so
`W` is a random walk with steps `-1` (probability `1/2`) and `2j-3` (probability `2^{-j}`, `j >= 2`),
mean `+1`. The probability `r` that it ever descends by one satisfies
`r = 1/2 + sum_{j>=2} 2^{-j} r^{2j-2}`, i.e. `r^3 - 2r + 1 = 0`, `r = (sqrt 5 - 1)/2 = 1/phi`. Hence

    P(halt | model) = r^6 = phi^{-6} = 0.05572...

This is a model statement, not a statement about the actual orbit of `27` (THM-4507 and THM-4594
show why finite-prefix certificates cannot decide such questions). The golden ratio enters through
the step law `2v - 3`, not through any Fibonacci structure of the orbit.

## A6. Exact orbit data (2,000,000 steps; receipts `antihydra_2e6.json`, `walk-normal.json`)

1,001,065 odd and 998,935 even steps; `b = 996,805` at the end; the counter identity
`b_k = 5 + W_{k-1}` holds at every odd term; the minimum of `b` after the first odd term is `5` and
`W` never drops below `0`; valuation frequencies `0.5008, 0.2509, 0.1251, 0.0619, 0.0310, 0.0156, ...`
match `2^{-j}`; `a` has 1,169,928 bits. None of this is evidence beyond the model: a halting
excursion needs six net descents and has model probability `phi^{-6}` from the start.

## Boundaries

A1–A4 and A6 are exact; A5 is conditional on the i.i.d. model; the rule constants are from the
public description. No claim that Antihydra halts or does not halt, and no claim about BB(6).
