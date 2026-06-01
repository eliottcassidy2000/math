# LRC Proof Route Landscape S505

The user asked for a long search for proofs of the LRC.  I treated LRC as the
Lonely Runner Conjecture, matching the repo's notation, and split the search
into accepted-looking proof routes, new preprints, variant results, and
low-trust complete-proof claims.

## Search Result

I did not find an accepted complete proof of the full Lonely Runner Conjecture.
The credible frontier has moved sharply, though: as of June 1, 2026, the newest
arXiv fixed-case claim I found is Sungkawichai-Trakulthongchai's April 26, 2026
preprint proving LRC(k) for `k <= 12`, i.e. up to thirteen total runners in the
runner formulation.  The general conjecture is still treated as open in the
normal literature.

I found several "complete proof" style claims outside the standard arXiv/journal
pipeline.  I did not admit them into the atlas because the current high-trust
literature still presents LRC as open and because those claims did not offer a
reproducible route comparable to the finite-sieve papers.

## Main Proof Line: Product Sieve

The modern fixed-case route is now clear:

```text
Tao finite checking
  -> Malikiosis-Santos-Schymura linearly-exponential finite bound
  -> Rosenfeld product-divisibility proof for eight
  -> Trakulthongchai fiber sieve for nine and ten
  -> Sungkawichai-Trakulthongchai intermediate sieves + polynomial method
```

The core engine is a contradiction on a minimal counterexample.  Finite-checking
gives an upper bound on the product of speeds.  For many primes `p`, a finite
ansatz check proves that every counterexample must have some speed divisible by
`p`.  If the forced product of primes exceeds the finite-checking product bound,
there is no counterexample.

Sungkawichai-Trakulthongchai sharpen this by replacing one giant lift with
intermediate sieves and by proving eventual properness of the tight residue
tuple `(1,2,...,k)` with a polynomial-method argument when `k+1` and `p>k^2+k`
are odd primes.

## General Tools

Tao's paper remains the conceptual hinge: it reduces fixed `n` to bounded
integer speeds and improves the trivial asymptotic lower bound.

Malikiosis-Santos-Schymura is the current finite-checking hinge: their
zonotope route replaces Tao's `n^{O(n^2)}` range with
`binomial(n+1,2)^(n-1) <= n^(2n)` and is now the bound used by the
fixed-case papers.

Bedert's Riesz-product paper does not prove LRC, but it gives the most relevant
general spectral pressure tool: a polynomial improvement of the global
loneliness lower bound beyond Tao's logarithmic improvement.

Jensen's mixed-threshold paper is also not a full LRC proof, but it is highly
relevant to the repo.  It supplies safe/unsafe indicator functions, an
arithmetic-progression summation formula, an unequal-threshold pair integral,
and an exact polygon for `MLPS_2`.

## Variant Warning

Shifted LRC and zonotopal covering-radius work is useful, but it is not a safe
substitute for the original conjecture.  Alcantara-Criado-Santos prove the
shifted conjecture for five total runners via dyadic fundamental-domain
certificates.  Then Blanco-Criado-Santos show explicit counterexamples to the
shifted LRC starting at `n=5` and to the Lonely Vector Property starting at
`n=12`.  The lesson for this repo: preserve the original unshifted endpoint
incidence object before importing shifted geometry.

## Tournament Analysis

I added `04-computation/lrc_proof_landscape_s505.py` and stored
`05-knowledge/results/lrc_proof_landscape_s505.out`.

The script treats proof routes as vertices.  The pairwise observable is a
feature-wise comparison of two routes.  The two gauges are:

```text
actionability = majority(repo_signal, boundary_overlap, proof_engine,
                         reproducible, frontier_signal)
credence      = majority(trust, direct_lrc, scope_strength,
                         reproducible, generality)
```

The tie Hamiltonian path is chronological.

Fingerprints:

```text
Actionability:
  SCC sizes [1,...,1], directed 3-cycles 0, Hamiltonian paths 1
  top routes ST13, Repo26, T9T10, MSS25, R8, R9

Credence:
  SCC sizes [8,1,1,1,1,1], directed 3-cycles 9, Hamiltonian paths 99
  top routes T9T10, ST13, MSS25, BS07, R9, R8

Gauge flips:
  29 of 78 edges flip between credence and actionability.
```

Interpretation: peer-reviewed older foundations win many credence edges, while
the newest product-sieve and repo endpoint-pressure routes win actionability
edges.  This is exactly the split we should expect in an active proof frontier.

## New Hypothesis

I added HYP-1969:

```text
LRC proof routes split into finite-sieve, kernel, zonotope, and endpoint-pressure currencies.
```

The concrete next move is to translate the product-sieve object `J(k,p)` or
`I(k,p,l)` into the repo's endpoint row-cover language.  The finite-sieve papers
prove that improper residue tuples disappear under enough ansatz refinement.
The repo should ask whether those disappearing tuples are exactly pressure-DAG
leaf peelings or endpoint-private row certificates in disguise.

## Sources

- Tao, "Some remarks on the lonely runner conjecture":
  https://doi.org/10.55016/ojs/cdm.v13i2.62728
- Perarnau-Serra, "The Lonely Runner Conjecture turns 60":
  https://arxiv.org/abs/2409.20160
- Malikiosis-Santos-Schymura, "Linearly-exponential checking is enough...":
  https://arxiv.org/abs/2411.06903
- Rosenfeld, "The lonely runner conjecture holds for eight runners":
  https://arxiv.org/abs/2509.14111
- Trakulthongchai, "Nine and ten lonely runners":
  https://arxiv.org/abs/2511.22427
- Rosenfeld, "The lonely runner conjecture holds for nine runners":
  https://arxiv.org/abs/2512.01912
- Sungkawichai-Trakulthongchai, "Eleven, twelve, and thirteen lonely runners":
  https://arxiv.org/abs/2604.23906
- Bedert, "Riesz products and the Lonely Runner Conjecture":
  https://arxiv.org/abs/2511.16636
- Jensen, "Mixed thresholds in the Lonely Runner Conjecture":
  https://arxiv.org/abs/2605.27941
- Alcantara-Criado-Santos, "Covering radii of 3-zonotopes...":
  https://arxiv.org/abs/2506.13379
- Blanco-Criado-Santos, "Coloopless zonotopes and counterexamples...":
  https://arxiv.org/abs/2603.24784
- Quanta overview of the 2025-2026 proof surge:
  https://www.quantamagazine.org/new-strides-made-on-deceptively-simple-lonely-runner-problem-20260306/
