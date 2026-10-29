---
id: THM-4845
title: "A000568 as a Gaussian average, and the meromorphic normalized series. (A) Let X be the centred Gaussian process on the odd integers with E[X_p X_q] = ln2 gcd(p,q) (realised by X_p = sqrt(ln2) sum_{d|p} sqrt(phi(d)) xi_d, xi_d i.i.d. N(0,1), because sum_{d|g} phi(d) = g) and W_p = 2^(-1/2) e^(X_p). Then, coefficientwise, sum_n a_n x^n = E[exp(sum_{p odd} W_p x^p / p)], i.e. a_n = E[Z_n(W)] with the Newton recursion n Z_n = sum_{p odd <= n} W_p Z_{n-p}, Z_0 = 1: the Davis/Burnside quadratic exponent is diagonalised by one independent Gaussian per odd divisor, and the core kernel B[m,t] of THM-514/516 is E[exp(v sum_{p>=3 odd} W_p x^p/p)], answering OPEN-Q-105(a) / HYP-2538 in this form. (B) The same field on all integers gives graphs (A000088) with W_j = 2^(-[j odd]/2) e^(X_j); tournaments are graphs with the even power sums switched off; doubling the covariance gives digraphs (A000273). (C) The Cameron-Martin formula E[e^(X_p) Phi(X)] = 2^(p/2) E[Phi(X + ln2 gcd(p,.))] is the inherited cycle-insertion PDE. (D) [concurrent with session98's ALL_POLES_CORE_EXPANSION, which has priority for (D)-(E)] f_n = n! a_n / 2^C(n,2) has F(y) = sum f_n y^n = sum over permutations mu of {3,5,7,...}-cycles of (|mu|!/z_mu) 2^(e(mu) - C(|mu|,2)) y^|mu| / (1 - y 4^(-k(mu)))^(|mu|+1), k(mu) = (|mu| - #cycles)/2: F is meromorphic on C with poles exactly at 4^k (k >= 0) of order exactly 3k+1, and f_n = sum_k 4^(-kn) P_k(n) with explicit polynomials P_k of degree 3k (P_1 = 32 C(n,3), P_2 = 98304 C(n,5) + 655360 C(n,6)); for graphs the poles are at every 2^d, of order 2d+1 (Q_1 = 4 C(n,2), Polya's correction). (E) Neither (a_n) nor (f_n) is P-recursive; F is not D-finite. So no linear recurrence with polynomial coefficients in n alone exists, and the Gaussian-conditional Newton recursion (A) is the recursion that does give the generating function."
status: >
  PROVED (elementary; mac-mini-2026-10-29-growth; pre-audit). (A)-(C) are exact identities between formal power
  series (each coefficient a finite expectation of lognormal monomials); (D) is THM-514's core-kernel factorization
  resummed, with a convergence argument for the meromorphic continuation; (E) follows from growth (a_n) and from the
  infinitely many poles of F (f_n). FINITE-EXACT verification (gf_gaussian.out): Davis a_n through n = 80; the
  Gaussian identity coefficientwise through n = 30; graphs through 12, digraphs through 8 and oriented graphs
  through 9 against OEIS; 1066
  coefficient identities of the PDE through degree 21; the exponential-polynomial expansion exact through n = 80 and
  for graphs through 30. CITED: Davis 1953/54 and Polya (Burnside counts), Smith 1875 (gcd matrix = Z^T Phi Z),
  the Cameron-Martin formula, the classical facts that P-recursive sequences grow at most like (n!)^r K^n and that a
  D-finite function has finitely many singular points. INHERITED: THM-514 (unit tail, core kernel K[m,t]), THM-516
  (e = C(t,2) + 1/2 sum phi(d) M_d^2), the cycle-insertion PDE and normalized recursion recorded by session98
  (05-knowledge/results/session98_tournament_growth_20261029/recurrence/PROJECTIVE_RESPONSE.md, concurrent).
  CONCURRENT, WITH PRIORITY ELSEWHERE: (D) and (E) were derived and published first by session98
  (05-knowledge/results/session98_tournament_growth_20261029/limits/ALL_POLES_CORE_EXPANSION.md, pushed 2026-10-29
  08:43, read only after this session's checkpoint). Their (1) is our (5), with the same exact pole set and orders and
  the same two non-P-recursiveness arguments. Their PALETTE_POLE_CLASSIFICATION (q pair colours, t fixed by reversal)
  contains our graph statement as q = t = 2. The derivation here is independent and agrees exactly; priority for (D)
  and (E) is theirs.
  NEW HERE: (A)-(C), and the Gaussian form of session98's whole palette family (B). Not found in the repository; the
  one-variable labelled analogue (2^C(n,2) as a lognormal moment) is folklore, and no literature search for the
  unlabelled gcd-field form was made.
session: mac-mini-2026-10-29-growth
source: 05-knowledge/results/tournament_growth_gaussian_gf_runs_fermat_20261029.md
depends_on:
  - 01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md
  - 01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md
related:
  - 05-knowledge/results/session98_tournament_growth_20261029/recurrence/PROJECTIVE_RESPONSE.md
  - 05-knowledge/results/tournament_future_boundary_20261020/README.md
  - 05-knowledge/results/weighted_core_transport_20261020/README.md
  - 01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md
  - 05-knowledge/hypotheses/HYP-9445-run-observer-torsion-and-hamiltonian-cycles-star-law-fermat-routes-and-gaussian-gf-questions.md
scripts:
  - 04-computation/experiments/growth_gf_runs_fermat_20261029/gf_gaussian.py
outputs:
  - 04-computation/experiments/growth_gf_runs_fermat_20261029/out/gf_gaussian.out
---

# THM-4845 — A000568 as a Gaussian average, and the meromorphic normalized series

## 0. Setting

For a multiplicity vector `m = (m_1, m_3, m_5, ...)` (finitely many nonzero entries, odd indices) write
`|m| = Σ m_p` for the number of cycles, `‖m‖ = Σ p m_p` for the size, `z_m = Π p^(m_p) m_p!`, and let `G` be the gcd
matrix, `G_pq = gcd(p,q)`. The Davis/Burnside count of tournaments (A000568) is

    a_n = Σ_{‖m‖ = n} 2^(e(m)) / z_m,   e(m) = (mᵀ G m − |m|)/2 = C(|m|,2) + ½ Σ_{d ≥ 3} φ(d) M_d²,   (0)

where `M_d = Σ_{d | p} m_p`. The second form is THM-516 (Smith: `G = Zᵀ Φ Z` with `Z` the divisibility matrix and
`Φ = diag φ(d)`). The first form follows from `e(m) = Σ (p−1)m_p/2 + Σ p C(m_p,2) + Σ_{p<q} gcd(p,q) m_p m_q` by
collecting the diagonal `gcd(p,p) = p`.

## 1. Statements

**(A) Gaussian representation.** Let `(ξ_d)_{d odd}` be independent standard normals and put

    X_p = √(ln 2) Σ_{d | p} √φ(d) ξ_d,   W_p = 2^(−1/2) e^(X_p)   (p odd).

Then `E[X_p X_q] = ln 2 · gcd(p,q)`, and as formal power series, coefficient by coefficient,

    Σ_n a_n x^n = E[ exp( Σ_{p odd} W_p x^p / p ) ].                                              (1)

Equivalently `a_n = E[Z_n(W)]`, where `Z_n` is the cycle index of `S_n` with the even power sums set to zero:

    Z_0 = 1,   n Z_n = Σ_{p odd ≤ n} W_p Z_{n−p}.                                                   (2)

Conditionally on the Gaussian field, (2) is the classical Newton recursion and its generating function is the
exponential in (1). The tournament count is the average.

*Core kernel.* With THM-514's core kernel `B[m,t]` (cores = permutations all of whose cycles are odd and ≥ 3),

    Σ_{m,t} B[m,t] x^m v^t = E[ exp( v Σ_{p ≥ 3 odd} W_p x^p / p ) ].                              (3)

This is the "Euler-product-like form over part sizes, modulated by the cross-terms" that OPEN-Q-105(a) and HYP-2538
asked for. The product over part sizes `p` sits inside the expectation. The coupling `½ Σ φ(d) M_d²` is carried
entirely by the shared field: one independent Gaussian per odd divisor `d`, with variance `φ(d) ln 2`. In other words
`log W` is the divisor sum of independent Gaussians, `log W = 1 ∗ g` (Dirichlet convolution, `g(d) = √(φ(d) ln2) ξ_d`),
shifted by `−½ ln 2`.

**(B) One field for graphs, tournaments and digraphs.** Extend the field to all `j ≥ 1` by the same formula
(`d` over all divisors) and put `W_j = 2^(−[j odd]/2) e^(X_j)`. Then

    graphs (A000088):      Σ g_n x^n = E[ exp( Σ_{j ≥ 1}   W_j x^j / j ) ],
    tournaments (A000568): Σ a_n x^n = E[ exp( Σ_{j odd}   W_j x^j / j ) ]   (the same W_j),

so **tournaments are graphs with the even power sums switched off**. Two more members of the family:
- *Digraphs (A000273).* Covariance `2 ln 2 · gcd` and `W_j = e^(X_j)/2`.
- *Oriented graphs (A001174).* The base-3 field: covariance `ln 3 · gcd`, with `W_j = 3^(−1/2) e^(X_j)` for odd `j`
  and `W_j = 3^(−1) e^(X_j)` for even `j`. A non-reversing pair orbit has 3 states and the one antipodal orbit of each
  even cycle has 1 state.
- *Session98's palette family.* This has `q` pair colours, `t` of them fixed by endpoint reversal. The field has
  covariance `ln q · gcd`, with `W_j = q^(−1/2) e^(X_j)` for odd `j` and `W_j = (t/q) e^(X_j)` for even `j`. Then
  `E[Π W_j^(m_j)] = q^(u(m)) t^(v(m))`, where `u` counts the non-antipodal pair orbits and `v` the even cycles, which is
  their Burnside count.
  - The special cases are tournaments `(2,0)`, graphs `(2,2)`, oriented and semicomplete graphs `(3,1)` and
    digraphs `(4,2)`.
  - `t = 0` switches the even power sums off.
  - Checked for nine `(q,t)` through order 9.
- *Self-converse tournaments (A002785, the repo's SC(n)).* Added after the audit was launched. A permutation `σ` maps
  a labelled tournament to its converse, `σ·T = T^op`, only if `σ` has at most one fixed point and every other cycle
  has length `≡ 2 (mod 4)`. In that case there are `2^((mᵀGm − m_1)/2)` such tournaments:
  - each pair orbit must flip orientation at every step;
  - a non-reversing orbit therefore needs even length;
  - the antipodal orbit of a `2r`-cycle needs odd length `r`.

  Twisted Burnside then gives, with the same field,

      Σ SC(n) x^n = E[ (1 + W_1 x) · exp( Σ_{j ≡ 2 (mod 4)} e^(X_j) x^j / j ) ].

  This reproduces OEIS A002785 for `n = 1..16`: `1, 1, 2, 2, 8, 12, 88, 176, 2752, …`.

**(C) Cameron–Martin is the inherited PDE.** Put `F(y) = E[exp(Σ_p y_p W_p / p)]`, the multivariate form whose
specialisation `y_p = x^p` is (1). Since `E[e^(X_p) Φ(X)] = e^(Var X_p / 2) E[Φ(X + Cov(X_p, X))]`,

    ∂F/∂y_p = (2^((p−1)/2) / p) · F( (2^(gcd(p,j)) y_j)_j ).                                        (4)

This is exactly the cycle-insertion PDE (1) of session98's PROJECTIVE_RESPONSE (inherited there from
weighted_core_transport). The Gaussian field is therefore an explicit solution of that hierarchy. Each inserted
`p`-cycle is a Cameron–Martin shift of the field by `ln 2 · gcd(p, ·)`. The uniform dilation `x → 2^t x` is the
deterministic shift `X_j → X_j + t j ln 2`. Session98's projective state `s_j = h_j − j h_1` is the accumulated
shift modulo this dilation direction (interpretation, not used in any proof).

**(D) The meromorphic normalized series.** Let `f_n = n! a_n / 2^(C(n,2))`, so `f_n → 1`. For a core `μ`, write
`c = |μ|` for its size, `ℓ` for its number of cycles, and `k(μ) = (c − ℓ)/2` (an integer, since all cycles are odd).
Put `κ_μ = (c!/z_μ) 2^(e(μ) − C(c,2))`. Then

    F(y) := Σ_n f_n y^n = Σ_{cores μ} κ_μ y^c / (1 − y 4^(−k(μ)))^(c+1),                            (5)

    f_n = Σ_{k ≥ 0} 4^(−kn) P_k(n),   P_0 = 1,   P_k(n) = Σ_{m = 2k+1}^{3k} K[m, m−2k] 2^(2km − C(m,2)) C(n, m),   (6)

with `K[m,t] = Σ_{|μ| = m, ℓ = t} (m!/z_μ) 2^(e(μ))` the integer core kernel of THM-514. For each `n` only finitely
many terms are nonzero. Moreover:
- `F` extends to a meromorphic function on `C`.
- Its poles are exactly the points `4^k`, `k ≥ 0`, and the pole at `4^k` has order exactly `3k + 1`. The top order
  comes only from `μ = (3^k)`, whose coefficient of `C(n,3k)` in `P_k` is `((3k)!/(3^k k!)) 2^(3k² + k)`.

Explicitly,

    P_1(n) = 32 C(n,3),
    P_2(n) = 98304 C(n,5) + 655360 C(n,6),
    P_3(n) = 12079595520 C(n,7) + 45097156608 C(n,8) + 2405181685760 C(n,9),

so that

    a_n = (2^(C(n,2)) / n!) · (1 + 32 C(n,3)/4^n + (98304 C(n,5) + 655360 C(n,6))/16^n + ...).

The expansion is exact, not asymptotic. The remainder after the poles `k ≤ K` is `4^(−(K+1)n) P_{K+1}(n)(1 + o(1))`.

*Graphs.* The same resummation gives `n! g_n / 2^(C(n,2)) = Σ_d 2^(−dn) Q_d(n)`, with cores the fixed-point-free
permutations. The poles are at every `2^d`, of order `2d + 1`, and `Q_1(n) = 4 C(n,2)`. This is Pólya's correction
`1 + n(n−1)/2^(n−1) + ...`. Tournaments see only the even `d`, because odd cycles have even `ℓ − 1`.

**(E) No P-recursive recurrence.** The sequence `(a_n)` is not P-recursive: a P-recursive sequence satisfies
`|a_n| ≤ C (n!)^r K^n` for some `r`, while `a_n ≥ 2^(C(n,2))/n!`. The sequence `(f_n)` is not P-recursive either,
because `F` has infinitely many poles and a D-finite function has finitely many singular points. So no linear
recurrence in `n` with polynomial coefficients, for `a_n` or for its normalization, can produce the generating
function. What does exist:
- the infinite-state recurrences (2) and (4) and session98's normalized (4)/(5);
- the exponential-polynomial sum (6), which says that truncating to `K` poles gives a constant-coefficient linear
  recurrence `Π_{k ≤ K} (E − 4^(−k))^(3k+1) f ≡ 0` up to `O(4^(−(K+1)n) n^(3K+3))`.

**(F) Remark: the divergence is a large deviation of the field.** For almost every realisation, `|X_p| = O(√(p log p))`
by Borel–Cantelli. So `W_p^(1/p) → 1` and the random series `exp(Σ W_p x^p/p)` has radius of convergence exactly 1,
while its coefficientwise mean has radius 0. The dominant term `2^(C(n,2))/n!` is the contribution of atypically
large `ξ_1`, the direction that counts the cycles.

## 2. Proofs

**(A).** `E[X_p X_q] = ln 2 Σ_{d | p, d | q} φ(d) = ln 2 · gcd(p,q)` by Gauss's `Σ_{d|g} φ(d) = g`. Every finite
linear combination of the `X_p` is Gaussian, so

    E[Π W_p^(m_p)] = 2^(−|m|/2) E[exp(Σ m_p X_p)] = 2^(−|m|/2) exp(½ Var Σ m_p X_p) = 2^((mᵀGm − |m|)/2) = 2^(e(m)).

The coefficient of `x^n` in `exp(Σ W_p x^p/p)` is `Σ_{‖m‖ = n} Π W_p^(m_p)/z_m`, a finite sum, so taking expectations
termwise gives (0). Differentiating the exponential gives (2). (3) is the same computation restricted to parts ≥ 3,
with `v` marking cycles.

**(B).** A permutation of type `m` (all parts) has `Σ_{i<j} gcd(i,j) m_i m_j + Σ_i (i C(m_i,2) + m_i ⌊i/2⌋)` orbits on
unordered pairs. An `i`-cycle has `⌊i/2⌋` internal orbits: `(i−1)/2` of size `i` when `i` is odd, and for even `i`,
`i/2 − 1` of size `i` plus one antipodal orbit of size `i/2`. Collecting terms, this equals
`(mᵀ G m − Σ_{i odd} m_i)/2`, and `E[Π W_j^(m_j)]` with the stated `W_j` equals 2 to that power.

For tournaments an even cycle reverses its antipodal pairs, so it fixes no tournament. Deleting the even power sums
gives Davis's formula, and on odd indices the two fields coincide.

For digraphs, the number of orbits on ordered pairs of distinct points is `mᵀ G m − |m|`, and the doubled covariance
gives `E[Π (e^(X_j)/2)^(m_j)] = 2^(mᵀGm − |m|)`. For oriented graphs, the unordered pair orbits number
`(mᵀGm − Σ_{j odd} m_j)/2`. Of these, `Σ_{j even} m_j` are reversing (the antipodal orbit of each even cycle) and
admit only "no arc". The count is therefore `3` to the power `(mᵀGm − Σ_{j odd} m_j)/2 − Σ_{j even} m_j`, again
quadratic in `m`, and the base-3 field reproduces it.

**(C).** For a centred Gaussian vector, `E[e^(⟨a,X⟩) Φ(X)] = e^(½⟨a,Ca⟩) E[Φ(X + Ca)]` (Cameron–Martin, or completing
the square). Take `a = e_p`, so `Var X_p = p ln 2` and `(Ca)_j = ln 2 · gcd(p,j)`. Then
`∂F/∂y_p = E[(W_p/p) exp(Σ y_j W_j/j)] = (2^(−1/2)/p) 2^(p/2) E[exp(Σ_j y_j 2^(gcd(p,j)) W_j / j)]`, which is (4).

**(D).** THM-514 gives `n! a_n = Σ_{m,t} C(n,m) K[m,t] 2^(C(n−m,2) + (n−m)t)`. Using
`C(n,2) − C(n−m,2) = C(m,2) + m(n−m)`,

    f_n = Σ_{m,t} K[m,t] 2^(−C(m,2)) C(n,m) 2^(−(m−t)(n−m)),

and `Σ_n C(n,m) z^(n−m) y^n = y^m / (1 − zy)^(m+1)` gives (5) termwise. Grouping by `k = (m−t)/2` gives (6).

*The range of m for fixed k.* Each cycle of length `ℓ ≥ 3` contributes `ℓ` to `m` and `ℓ − 1 ≥ (2/3)ℓ` to `m − t`.
Hence `m ≤ 3k`, with equality only for `μ = (3^k)`, and `t ≥ 1` gives `m ≥ 2k + 1`.

*Convergence.* A core of size `m` has all pair orbits of size ≥ 3, so `e(μ) ≤ C(m,2)/3` and

    Σ_{|μ| = m} κ_μ ≤ m! 2^(−2C(m,2)/3).

On a compact set `K` avoiding the points `4^j`, all but finitely many cores have `4^(k(μ)) ≥ 2 max_K |y|`. For those,
`|1 − y 4^(−k)| ≥ 1/2`, so their terms are bounded by `κ_μ (2R)^m · 2`, which is summable. The finitely many
remaining terms are rational functions with poles only at the excluded points. Hence (5) converges uniformly on `K`,
and `F` is meromorphic.

*The poles.* The pole at `4^k` collects the cores with `k(μ) = k`. Its top order `3k + 1` comes from `(3^k)` alone,
with a positive coefficient, so it does not cancel. For `μ = (3^k)`, `e = (3k² − k)/2` by (0), and
`2km − C(m,2) = (3k² + 3k)/2` at `m = 3k`, which gives the stated leading coefficient.

*Graphs.* Repeat the argument with THM-514's peel for graphs: fixed points interact with every cycle through
`gcd(1,ℓ) = 1`.
- Cores are now fixed-point-free permutations, and `d = m − t` can be odd (even cycles).
- Each cycle of length `ℓ ≥ 2` contributes `ℓ − 1 ≥ ℓ/2` to `d`, so `m ≤ 2d`. Only all-2-cycle cores reach
  `m = 2d`, and they give the top order `2d + 1`.
- For convergence: a pair orbit of a fixed-point-free permutation has size at least 2, except the fixed pair of each
  2-cycle. So `e(μ) ≤ (C(m,2) + m)/2`, and `Σ_{|μ| = m} κ_μ ≤ m! 2^(−C(m,2)/2 + m/2)`, which is still summable against
  `(2R)^m`.

**(E).** For the growth bound on P-recursive sequences, see Stanley, *Enumerative Combinatorics II*, Ch. 6: the
solutions of a linear recurrence of order `r` with polynomial coefficients of degree at most `s` satisfy
`|a_n| ≤ C (n!)^s K^n`. Then `log(2^(C(n,2))/n!) ~ n² (ln 2)/2` beats `s n log n`. For `(f_n)`: a D-finite power
series satisfies a linear ODE with polynomial coefficients, whose singular points are among the finitely many zeros
of the leading coefficient. By (D), `F` has infinitely many poles.

**(F).** `X_p ~ N(0, p ln 2)`. By Borel–Cantelli, almost surely `|X_p| ≤ 2√(p ln 2 · ln p)` eventually, so
`W_p^(1/p) → 1`.

## 3. Verification (FINITE-EXACT, `out/gf_gaussian.out`)

- Davis `a_n` computed through `n = 80`; the first 16 terms agree with OEIS A000568.
- `Σ_{d|p, d|q} φ(d) = gcd(p,q)` for `p, q < 60`.
- `a_n = E[Z_n(W)]` through `n = 30`, computed as polynomials in `W` from (2), then the moment rule.
- Graphs `g_n` through 12, digraphs through 8 and oriented graphs through 9 from the same functional, agreeing with
  OEIS A000088, A000273 and A001174.
- The palette family: nine `(q,t)` through order 9 against the direct Burnside sum `Σ q^u t^v / z`.
- (4) checked as 1066 coefficient identities through degree 21.
- (6) exact through `n = 80`. After the poles `k ≤ K` (`K = 0, 1, 2`), the remainder divided by the next term is
  `1.0000007`, `1.000014`, `1.0003` at `n = 20`, and `1` to twelve digits at `n = 40`.
- The graph expansion is exact through `n = 30`.

## 4. Inheritance and limits

- **Inherited.** THM-514 and THM-516 supply the core kernel and the quadratic form. Session98 (concurrent, 2026-10-29)
  supplies the projective fixed-degree recurrence and the formal PDE/hierarchy, inherited from weighted_core_transport
  and tournament_future_boundary. This theorem does not replace those recurrences. (A) is their explicit solution,
  (C) their mechanism, and (E) explains why nothing finite-dimensional can replace them.
- **Concurrent, with priority to session98.** The meromorphic continuation, its exact poles and the
  non-P-recursiveness, (D)-(E), are session98's ALL_POLES_CORE_EXPANSION, found independently here. Their palette
  theorem generalizes our graph remark.
- **Not claimed.**
  - A closed form for `B[m,t]` beyond the Gaussian average.
  - A faster algorithm. Evaluating (1) still sums over cycle types, and the 200-term computations of
    tournament_future_boundary remain the computational state of the art in the repo.
  - Any asymptotic statement beyond Moon's classical leading term, except as the exact expansion (6).
- **Novelty.** The formulas are elementary rearrangements of Davis and Pólya. What appears new to the repository is
  the gcd-covariance field, its identification with the PDE, the shared field for graphs and tournaments, and the
  meromorphic continuation with its exact pole orders. No literature search beyond OEIS was made.
