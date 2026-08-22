# Core papers: homotopy groups of spheres (opened 2026-08-02)

Overflow entry for [`CORE-PAPERS.md`](CORE-PAPERS.md).  That file is at its
bounded startup-surface byte budget and could not take even a pointer line, so
this sits beside it under the deliberately adjacent name.  Same contract: what
the repository imports, where it is consumed, and what the source does **not**
establish.  A future compaction pass on the startup surface should add the
pointer.

## Ivanov--Mikhailov--Wu — *On nontriviality of homotopy groups of spheres*

- **Primary / freshness:** [arXiv:1506.00952v1](https://arxiv.org/abs/1506.00952),
  2015-06-02, the only arXiv version.  Published as *Homology, Homotopy and
  Applications* **18** (2016), no. 2, 337--344,
  [DOI 10.4310/HHA.2016.v18.n2.a18](https://doi.org/10.4310/HHA.2016.v18.n2.a18).
  **PUBLISHED / stable.**  Record checked 2026-08-02.
- **Theorem:** `pi_n(S^2) != 0` for all `n>=2`, hence `pi_n(S^3) != 0` for
  `n>=3`.  Curtis (Bull. AMS **75** (1969) 541--544) had left the residue class
  `n = 1 mod 8` open; the paper fills it by proving
  `Z/p subset pi_(2(p-1)k+1)(S^3)` for every odd prime `p` and every `k>=1`,
  giving in particular `Z/3 subset pi_(4n+1)(S^3)` and
  `Z/15 subset pi_(8n+1)(S^3)`.
- **Attribution note the authors make themselves:** Brayton Gray,
  *Unstable families related to the image of J*, Math. Proc. Cambridge Philos.
  Soc. **96** (1984) 95--113, Theorem 12(e), had already covered the same
  dimensions by a different method.  Do not present the dimension coverage as
  first-of-its-kind; the paper's contribution is the route and the assembled
  all-`n` statement.
- **Imported role.** The repository imports section 2 only: the presentation of
  the odd-primary lambda algebra (generators `lambda_i` for `i>=1` and `mu_j`
  for `j>=0`, admissibility, the four Adem-type relation families, and the
  differential), together with Lemma 3.  Consumed by
  [THM-3205](../../01-canon/theorems/THM-3205-odd-primary-lambda-algebra-engine-and-toda-gate-family.md),
  which replays all eleven displayed identities and strengthens Lemma 3 to an
  exact `E^2` statement with an explicit generator.  The `tridiag(1,-2,1)`
  determinant in that lemma is the motivating instance of
  [THM-3204](../../01-canon/theorems/THM-3204-parabolic-continuant-single-gate-and-jacobi-smith-obstruction.md).
- **Two convention warnings before extending it.**  Neither is an error in the
  paper; both are needed by anyone who computes beyond its displayed lines.
  1. In each relation the summand `j=k` returns the left-hand side itself with
     coefficient `a(k,k)=-1`.  It must be omitted.  Reading the displayed range
     `0<=j<=N(k)` literally turns the paper's own
     `mu_0 lambda_2 = -mu_1 lambda_1` into `-(1/2) mu_1 lambda_1`.
  2. Leibniz carries the Koszul sign `d(xy)=d(x)y+(-1)^|x| x d(y)`, so only
     `lambda` factors flip it.  Without the sign `d^2 != 0` already on `mu_2`,
     where `d^2(mu_2)=4 lambda_1^2 mu_0`.  Every computation displayed in the
     paper is unaffected, because in each one `d` is applied to the right of an
     even-degree (`mu`) prefix only.
- **Does not prove:** any 2-primary statement.  Its Conjecture 1 (the
  2-component of `pi_n(S^3)` is nontrivial for `n>10`) is open, and the
  repository engine does not address `p=2`, where `Lambda` has no `mu`
  generators and the excess/stratum count of THM-3205 section 3 does not apply.
  The paper also states no `E^2` dimension, no explicit cycle, and no
  integrality result; those are THM-3205's additions.  Neither repository
  theorem proves a new homotopy group or any convergence statement.

## Supporting sources cited through the above

- A. K. Bousfield, E. B. Curtis, D. M. Kan, D. G. Quillen, D. L. Rector,
  J. W. Schlesinger, *The mod-p lower central series and the Adams spectral
  sequence*, Topology **5** (1966) 331--342 — the spectral sequence whose
  integral/unstable form is `E^1(n)=Lambda*lambda(n) => (p) pi_*(S^(2n+1))`.
  **CITED, not re-derived here.**  THM-3205 computes only `E^2` strata; it
  proves nothing about convergence or higher differentials.
- M. C. Tangora, *Computing the homology of the Lambda algebra*,
  Mem. Amer. Math. Soc. **337** (1985) — source of the `alpha_k = mu_1^(k-1)
  lambda_1` identification and of the `h_p` short exact sequence used in the
  paper's Proposition 2.  **CITED-ABSTRACT only**; not independently checked in
  this repository.
