# Independent audit of the clone obstruction

**Mathematical result: PASS.** The author applied the requested singleton qualification in §3 and the finite-check scope correction in §6; both changes were verified. The bound applies the induced-triangle argument only to nontrivial strongly connected quotient components, whose orders are at least three. A singleton quotient is also conventionally strong, but its transitive inflation has H=1 for arbitrarily large module size, so the intermediate parameter bound must not be stated for it. The general SCC argument correctly treats singleton components separately. This qualification does not change the theorem or its conclusion.

The induced-wrapper count is exact: in cyclic orientation X→a→b→X, every path with consecutive a,b corresponds to one subset of the transitive block X; the sole remaining path is b,X,a. Thus H=2^r+1. Every vertex of a nontrivial strong tournament lies in a directed triangle, and insertion monotonicity carries the induced lower bound to the entire inflation. These facts imply r_i<=floor(log_2(N-1)) whenever its strong inflation has H<=N.

Ordered sums multiply H, so a prime output has exactly one nonunit factor. Nontrivial strong components have H>=3, singleton components contribute 1, and arbitrary further ordered sums and reversals cannot introduce new prime values. A finite library with at most Kd relevant components therefore supplies at most Kd L^d prime values. The proof actually overcounts, harmlessly, components of templates that could never themselves yield a prime.

The elementary prime-density argument is valid. In the central binomial coefficient each summand in the p-adic valuation is 0 or 1, yielding p^valuation<=2m. Its prime-factor product is therefore at most (2m)^pi(2m), while the largest-binomial-coefficient bound gives at least 4^m/(2m+1). Taking m=floor(N/2) extends the displayed even-bound estimate to all sufficiently large real or integer N. Hence pi(N) is bounded below by a positive constant times N/logN, which dominates every fixed power of logN.

This is a limitation on a finite constructor class. It implies neither global nonattainment of any particular integer nor failure of growing-interface or internally nontransitive constructions. The stated restriction on finite library size and bounded quotient order is essential.

## Reproducible evidence

- Replayed the author's checker under `python -O`; its receipt exactly matches `certificate.json`, including 1998 strong inflation tests.
- `reviewer_audit.py` imports no author counting code. It uses independent Boolean adjacency matrices, transitive-closure strongness, and recursive endpoint counting.
- 400 independent blowups with all module sizes 1 or 2 over all nontrivial strong labeled quotients of orders 3 and 4.
- 134 additional literal-permutation count comparisons.
- 16 ordered-sum multiplicativity checks.
- Inflated triangle formula through r=9.
- 100 exact central-binomial prime-bound checks.
- Explicit singleton boundary: transitive inflation of size 10 has H=1.

Normal and optimized reviewer receipts are byte-identical. The author checker source SHA-256 is pinned in those receipts. The finite checks supplement the universal proof and are not substituted for it.
