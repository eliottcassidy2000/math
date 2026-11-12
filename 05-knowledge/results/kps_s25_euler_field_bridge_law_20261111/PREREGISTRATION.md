# Pre-registration (kind-pasteur-2026-11-11-S25)

Owner directive (2026-11-11): continue the session on Stein modules of Keller pairs (S24), and extend it into a
session that aims to prove the Chan-Pak conjecture. Both conjectures remain OPEN at the time of writing.

## Already done before this file (disclosed; not predictions)

- **Stein thread.** `scripts/euler_field.py` was written and run (`runs/euler_field.txt`): the Euler-field lemma
  `W h = (xi + delta) h` modulo `Im(d/dx - W_x) + Im(d/dy - W_y)` for a vector field with `xi(W) = W` and constant
  divergence `delta`, and its consequence that the image for Broughton's polynomial `W = x + x^2 y` is not a
  Mathieu-Zhao subspace (`W^m` in the image for all `m >= 1`, `x W^m = m! x` not in it): HYP-9668 (E') is false.
- **Chan-Pak, sweep.** A sub-agent digest of the repo (2026-11-11, after commit f90b22ab5d) found that the pieces I
  had re-derived from the odd-cycle formula are KNOWN: the vertex-addition identity (THM-4876), the join law of two
  rooted tournaments `H = x_h y_h + x_R y_R - x_e y_e` (THM-4866 (B); Session263, Theorem RJ2), the one-apex function
  `g` (THM-4875, THM-4890), the two-parameter family `(1+u)(1+v) + g(u) g(v) - e(u) e(v)` with its mass
  `(pi^2/24) X ln X` (THM-4885 (B)), and completeness of "bank x one-apex" joins up to `2^35` (THM-4885 (D)).
  A scratch run of mine (not in the repo) had recomputed the second family below `2^20` (10,895 odd numbers missed,
  the last dyadic block at 1.5 per cent); it is a re-computation of a known family, not a result.
- Concurrent: opus-2026-11-11-S93 (counting exponent 3/4 for one-apex counts, two-apex words complete below `2^23`),
  opus-2026-11-11-S94 (proof-oriented, just started), Sessions 275 and 276 (digit banks, corridors).

## What is tested next (not yet run)

**Object.** For `p >= 0` let `L(p)` be the bit length of `p` (`L(0) = 0`), `R(p)` the reversal of the complement of
`p` on `L(p)` bits (`R(0) = 0`), and `g` the one-apex function. For `i >= 1` put

    M(p, i) = [[ g(p), 2^i (R(p) + 1) ], [ p + 1, 2^(i + L(p)) ]].

Working claim (to be proved in the note from the KNOWN join law): for one-apex tournaments `Y_0, ..., Y_k` with apex
words of values `p_0, ..., p_k`, the tournament `Y_0 => C_1 => Y_1 => ... => C_k => Y_k` (`C_j` a transitive chain of
`i_j >= 1` vertices) with the arcs between consecutive apexes reversed has

    H = top-left entry of M(p_k, i_k) ... M(p_1, i_1) M(p_0, 1).

So the path-count spectrum contains the set `Sigma` of top-left entries of all products of the matrices `M(p, i)`.

**Predictions.**

- **C1.** Brute force (dynamic programming over subsets, at most 13 vertices): the formula above agrees with the
  number of Hamiltonian paths for every chain of depth 2 and 3 with words of length at most 2 and chain lengths 1,
  2 (0 mismatches).
- **C2.** Depth at most 2 (products of two matrices) misses more than 2 per cent of the odd numbers of
  `[2^17, 2^18)`.
- **C3.** Depth at most 3 misses fewer than 0.5 per cent of the odd numbers of `[2^17, 2^18)`, and fewer than 600
  odd numbers below `2^18` in all.
- **C4.** Full closure (all depths), computed exactly for all odd `N <= 2^14`: the only odd numbers missed are 7 and
  21. (If this fails, the list of misses is the result.)
- **C5.** The number of products of depth at most 3 with top-left entry at most `2^18` (counted with multiplicity
  as tuples) is at least `40 * 2^17`: mean multiplicity at least 40 per odd number.

**Reading.** C4 true would give a bank-free sharpening of the known finite completeness: the one-apex chain
semigroup alone, with no bank of small tournaments, reaches every odd number in the range. It would not prove
anything about all `N`. C4 false would show that the pure chains need the bank (or the joins without an inserted
chain), and the misses would say which.
