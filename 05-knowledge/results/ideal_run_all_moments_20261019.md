# All factorial moments of an ideal run

2026-10-19. A reusable consequence of Stanley's rank log-concavity for boundary interleaving counts. This note extends the quadratic lemma used in `rooted_five_logconcave_balance_20261019.md`; it does not assert a new solved rooted class by itself. No novelty claim is made.

## The bound

Let $K$ be an ideal in a finite poset $B$, and let $R$ be the initial $K$-run length in a uniform linear extension of $B$. Put $m=\mathbb E R$. For every integer $r\ge0$,

\[
\boxed{\mathbb E\binom Rr\le m^r},\qquad
\boxed{\mathbb E\binom{R+r}r\le(1+m)^r}.              \tag{1}
\]

The conventions for $r=0$ are one. The order-dual bounds hold for terminal runs of upsets.

**Proof.** The rooted-five paper proves, by adjoining a minimal marker below $B\setminus K$, that $q_j=\Pr(R\ge j)$ is a log-concave tail sequence with $q_0=1$ and initial positive support. It uses exactly the single-element specialization of Stanley's Theorem 3.1, [original primary source](https://math.mit.edu/~rstan/pubs/pubfiles/48.pdf). Hence

\[
q_{i_1+\cdots+i_r}\le q_{i_1}\cdots q_{i_r}
\quad(i_1,\ldots,i_r\ge1).
\]

For $r\ge1$, the tail-sum and composition identities give

\[
\mathbb E\binom Rr
=\sum_{s\ge r}\binom{s-1}{r-1}q_s
=\sum_{i_1,\ldots,i_r\ge1}q_{i_1+\cdots+i_r}
\le\left(\sum_{j\ge1}q_j\right)^r=m^r.
\]

Vandermonde's identity now yields

\[
\mathbb E\binom{R+r}r
=\sum_{j=0}^r\binom rj\mathbb E\binom Rj
\le\sum_{j=0}^r\binom rj m^j=(1+m)^r.
\]

All sums are finite. The empty ideal has $R=0$, so no exceptional case is omitted. $\square$

## Sharpness as a uniform bound

Take $B$ to be the disjoint union of a $k$-chain $K$ and an $\ell$-chain. A uniform linear extension is a uniform shuffle. The $K$ letters are divided among the $\ell+1$ gaps around the other chain, uniformly over weak compositions of $k$ into $\ell+1$ parts. The first gap is $R$. Thus

\[
\mathbb E R=\frac{k}{\ell+1},\qquad
\mathbb E\binom Rr
=\frac{k(k-1)\cdots(k-r+1)}{(\ell+1)(\ell+2)\cdots(\ell+r)}. \tag{2}
\]

One direct derivation of the second identity uses the weak-composition generating function. Marking $\binom Rr$ replaces the first gap's series $(1-t)^{-1}$ by $t^r(1-t)^{-r-1}$; division by the total composition count gives

\[
\frac{\binom{k+\ell}{\ell+r}}{\binom{k+\ell}{\ell}},
\]

which is (2). If $k/(\ell+1)\to m>0$ while $k,\ell\to\infty$, (2) tends to $m^r$. By Vandermonde, the second bound in (1) tends to equality as well. Hence neither constant in (1) can be improved uniformly over finite posets at a fixed limiting mean.

For an individual finite poset with nonempty $K$ and $r\ge2$, the first inequality is strict: the product expansion includes positive terms with all indices at the end of the positive support, whose sum is beyond that support and whose corresponding tail term is zero. The asymptotic sharpness is therefore not a claim of finite equality.

## Use and limitation

If an ordered $r$-point core prefix can interleave with an exterior run of length $R$ in at most $\binom{R+r}{r}$ ways, its average interleaving count is at most $(1+m)^r$. This turns an unbounded run into a degree-$r$ algebraic bound while preserving a necessary constraint of actual linear-extension laws.

The statement does not bound arbitrary joint features of two non-nested exception ideals. Their runs are measured in the same internal word and need not be independent. Replacing their mixed moments by products of means would require a separate theorem. The seven alternate-root sum certificates use exact feature counting instead of making that unsupported substitution.

## Reproduction

`verify_ideal_run_all_moments_20261019.py` independently enumerates the weak-composition shuffle laws for 54 choices of chain lengths and moment orders zero through six, checking (1), (2), the strict finite inequality, and the shifted-moment bound with exact rational arithmetic. Its receipt is `ideal_run_all_moments_checks_20261019.json`. These are formula controls, not the proof of the general result.

From the repository root run `python 05-knowledge/experiments/verify_ideal_run_all_moments_20261019.py`. All certificate guards remain active under `python -O`; the program writes its receipt to standard output. The receipt is [ideal_run_all_moments_checks_20261019.json](ideal_run_all_moments_checks_20261019.json).
