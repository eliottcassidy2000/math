# Recovering the missing prime converse-action sidecar

**Status:** PROVED deduction from cited modular decomposition and equivariant species calculus; FINITE-EXACT computation through weight 15; independent fixed-permutation and prime checks through 7. No literature-novelty claim. This packet does not reconstruct individual prime isotypes or their automorphism groups.

## Inheritance, distinction, and boundary

The [previous session](../session79_tournament_enumeration_20261027/PRIME_CYCLE_INDEX.md) recovered the ordinary aggregate prime cycle index, explicitly leaving anti-automorphisms unrecovered. Session67's klein packet supplied ordinary and anti cycle-type data through nine and self-converse prime counts through ten. We now recover the aggregate anti-index from the *full tournament twisted Burnside index*, without taking any supplied prime counts as input. The supplied data are only retrospective checks. This is not recovery from the ordinary prime index alone.

Concept board: canonical strong modules; pair-orbit transport; converse as an external C2 action; parity-sensitive substitution; reversal of SCC order. Hostile controls: even ordinary automorphism cycles, forbidden odd anti cycles, two fixed anti vertices, and treating SCC reversal as componentwise reversal.

For a tournament T, an automorphism preserves every arc, whereas an anti-automorphism sends every arc to its reverse. Write τ for converse. It commutes with relabeling, preserves primeness and strong connectivity, and acts on the canonical modular tree.

For any converse-stable species F define

    Z_F^g = Σ_n (1/n!) Σ_{σ∈S_n} |Fix(gσ on F[n])| p_type(σ),  g=e,τ.

Let B_e,B_τ be these indices for all tournaments including the empty one; A_g=B_g−1. Let S_g be the corresponding indices for nonempty strongly connected tournaments, including the singleton. Let P_g be the prime indices, excluding orders 1 and 2. A prime has no proper nonsingleton module.

The species decomposition T=SEQ(S), S=X+P∘A is the inherited Gallai mechanism. The canonical partition and its linear-or-prime quotient are cited in [Boussaïri–Chergui–Ille–Zaidi, Theorem 16](https://arxiv.org/html/1805.04607v1). Its compatibility with converse follows because reversing all arcs preserves exactly the modules and strong modules.

## 1. Exact twisted Burnside input

For a partition λ=(l₁,…,l_k), write z_λ=∏_d d^{m_d}m_d! and

    E(λ)=Σ_i floor(l_i/2)+Σ_{i<j} gcd(l_i,l_j).

A permutation of type λ is an anti-automorphism of some tournament **if and only if** every non-singleton cycle has length 2 modulo 4 and there is at most one singleton. When admissible, the number of tournaments on the fixed labeled set reversed by that specified permutation is 2^{E(λ)}. This includes the empty and singleton types.

**Proof.** On each orbit of unordered pairs choose one orientation, then propagate it using the anti condition. If an orbit first returns after L steps, transport must agree with reversal L times: the endpoint permutation at return must preserve the endpoints for even L and exchange them for odd L. Within an odd vertex cycle of length greater than one, any pair orbit has odd length and returns endpoints individually, so it is forbidden. Within a vertex cycle of length 2m, antipodal pairs have orbit length m and exchange endpoints at return; this is compatible exactly when m is odd. Other within-cycle pair orbits have even length and return endpoints individually. Between two distinct vertex cycles the pair-orbit length is their lcm and endpoints cannot exchange; hence it must be even. This rules out two singleton cycles, and is automatically satisfied between the other admissible cycles. These necessary conditions are sufficient on every pair orbit, independently. There are floor(l/2) internal pair orbits and gcd(l,m) between two cycles, proving the exponent.

Thus, with λ restricted to the indicated support,

    B_τ = Σ_λ 2^{E(λ)} p_λ/z_λ.

The ordinary input B_e uses all-odd partitions with the same exponent. Its even-cycle obstruction is endpoint reversal by a power of an ordinary automorphism.

Checks: fixed counts for λ=(2,1),(2,2) are 4 and 16, giving two self-converse isotypes at each of orders 3 and 4. Types (3),(4),(1,1) have zero anti-fixed tournaments.

## 2. Reversal of the SCC sequence is indispensable

The C2 action on SEQ reverses the linear order. Its two cycle indices are

    Z_SEQ^e = 1/(1−p1),        Z_SEQ^τ = (1+p1)/(1−p2).

For length 2k, the order-reversing permutation has k transpositions; for length 2k+1 it additionally has one fixed point. Summing over all labeled orders with the cycle-index normalization gives one monomial for each length, proving the second formula directly.

Write ψ_d F=F(p_d,p_{2d},…). The equivariant substitution law replaces an outer p_d at group element g by ψ_d Z_inner^{g^d}. The reason is that a child transported once around a d-cycle experiences g exactly d times. The general formula is cited in [Gainer-Dewar, *Species with an equivariant group action*, §3.2, equation (14), arXiv v4](https://arxiv.org/pdf/1401.6202). Applying it to the SCC sequence gives

    B_e = 1/(1−S_e),
    B_τ = (1+S_τ)/(1−ψ₂ S_e).

Therefore the full strong indices are explicitly known:

    S_e = 1−B_e^{-1},
    S_τ = B_τ(1−ψ₂ S_e)−1.                         (1)

In particular S_τ is NOT 1−B_τ^{-1}. At the scalar isotype level, b_τ(x)=(1+s_τ(x))/(1−s_e(x²)): a symmetric pair of positions can contain any strong tournament and its converse, not only self-converse children.

## 3. Unique recovery of the prime anti-index

Define the mixed substitution M by

    M(U)=U(q₁,q₂,…),    q_d=ψ_d A_τ if d odd,
                                  ψ_d A_e if d even.

The prime-root equations are

    P_e ∘ A_e = S_e−p1,      M(P_τ)=S_τ−p1.        (2)

Every A_g starts with p1, hence q_d starts with p_d. In the completed weighted power-series ring (weight(p_d)=d), M preserves the leading homogeneous part of every series. Consequently (2) has a unique recursively computable solution:

    (P_τ)_n = [S_τ−p1−M((P_τ)_{<n})]_n, n≥3.    (3)

The exact same argument handles the ordinary index. This is a multivariate, parity-sensitive triangular inversion, not an ordinary inverse of an isotype generating function.

**Interpretation.** P_τ=Σ_[Q self-converse] (1/|Aut Q|) Σ_{a∈Anti(Q)} p_type(a). The anti maps form a coset of Aut Q, so their number is |Aut Q|. Non-self-converse prime types contribute nothing. Therefore the sum of coefficients of degree n is exactly the number of self-converse prime isotypes. It is not half that number. A sum of weighted coset cycle types does not identify the individual cosets or groups.

Under the separately cited theorem that a prime tournament's C3 structure determines it up to converse, the aggregate prime C3 root index is (P_e+P_τ)/2. This conclusion uses the C3 theorem in addition to Gallai; the present fixed-permutation proof does not establish that theorem.

## 4. Exact output and independent checks

| n | self-converse tournaments | self-converse primes |
|---|---:|---:|
| 3 | 2 | 1 |
| 4 | 2 | 0 |
| 5 | 8 | 3 |
| 6 | 12 | 3 |
| 7 | 88 | 55 |
| 8 | 176 | 72 |
| 9 | 2752 | 2235 |
| 10 | 8784 | 5525 |
| 11 | 279968 | 259307 |
| 12 | 1492288 | 1170949 |
| 13 | 95458560 | 92996193 |
| 14 | 872687552 | 767626075 |
| 15 | 111698291584 | 110751820327 |

Examples:

    (P_τ)₃ = p1 p2,
    (P_τ)₅ = 3p1 p2²,
    (P_τ)₆ = (7/3)p2³+(2/3)p6,
    (P_τ)₇ = (161/3)p1 p2³+(4/3)p1 p6,
    (P_τ)₁₀ = (82759/15)p2⁵+(16/3)p2²p6+(12/5)p10.

`recover_converse.py` performs exact Fraction arithmetic, checks both equations coefficientwise through 15, nonnegativity, support, and integrality of all index totals. It requires no external packages or prime census input.

`check_equivariant.py` independently constructs one permutation of **every partition through seven**, solves orientation transport by XOR constraints on actual unordered pairs, exhausts every compatible orientation, and tests **every proper nonsingleton subset** for the module condition. Its fixed counts agree with the support formula and its prime counts agree with z_λ[p_λ]P_τ, including every zero-support control. This checks labeled fixed sets, not just scalar totals. A separate historical comparison matched the entire recovered anti-index through nine with the supplied table. That comparison is credited consistency evidence, not a rerun of the source census or a dependency of the published verifier.

Reproduce:

    python3 recover_converse.py 15
    python3 check_equivariant.py
    python3 -O recover_converse.py 15
    python3 -O check_equivariant.py

Normal and optimized receipts agree byte-for-byte. No Python assert is used. The order 8–15 outputs are exact consequences of (1)–(3), not independent exhaustive censuses.

## 5. Specified full even cycle: a scalar Dirichlet boundary

Let m run over positive odd integers, f_m=2^{(m−1)/2}, and let h_m be the number of prime tournaments reversed by the **specified** cycle (0 1 … 2m−1). Set h₁=0 under the convention excluding order two. Then

    (h*f)_m = 2^m−2f_m,
    h_m = 2^m−2f_m−Σ_{d|m,d<m} h_d f_{m/d}.       (4)

Here * is Dirichlet convolution on odd positive integers. The fixed-permutation theorem gives 2^m total anti-fixed tournaments; it also gives h_m=(2m)[p_{2m}]P_τ.

**Proof.** The coefficient of the single-factor monomial p_m in S_e equals [p_m]B_e=f_m/m: every nonlinear term in 1−B_e^{-1} has at least two factors. This includes m=1. In S_τ=B_τ(1−ψ₂ S_e)−1, the lone p_{2m} coefficient is therefore

    2^m/(2m)−f_m/m = (2^m−2f_m)/(2m).

A product of two nonconstant monomials cannot contribute a lone power sum. In M(P_τ), the only contributions to p_{2m} come from an outer single even cycle p_{2d} and an inner odd single cycle p_{m/d}, with d|m. Even outer cycles use A_e, hence their contributions are

    Σ_{d|m} (h_d/(2d)) (f_{m/d}/(m/d)).

Multiplication by 2m proves (4), including m=1 where both sides vanish. No prime coefficient of order two has been inserted artificially. Full cycles whose length is divisible by four have no anti-fixed tournaments at all, by the antipodal-pair obstruction; odd full cycles greater than one likewise have none.

`check_equivariant.py` directly enumerates all XOR-compatible orientations for the specified full cycle, independently of (4). It tests primeness by the following exact alternative to an exponential subset scan: for each vertex pair, repeatedly add every exterior vertex that distinguishes the current set. Every module containing the starting pair must contain every subsequently added vertex. At termination the set is itself a module, so a proper closure witnesses nonprimeness. Conversely every proper nonsingleton module contains a pair whose closure remains proper. Thus a tournament of order at least three is prime iff every pair closure fills the vertex set. The implementation cross-checks this criterion against the exhaustive all-subsets method through order ten.

| 2m | all specified-cycle anti-fixed | prime specified-cycle anti-fixed |
|---|---:|---:|
| 2 | 2 | 0 |
| 6 | 8 | 4 |
| 10 | 32 | 24 |
| 14 | 128 | 112 |
| 18 | 512 | 472 |

Each row agrees with both (4) and the independently recovered full cycle-index coefficient through weight 18. The XOR implementation also checks the divisible-by-four no-go at 4,8,12,16. These counts are labeled fixed points of one prescribed permutation, **not** self-converse isotypes, converse pairs, all anti-transitive isotypes, or a quotient by rotations. No division by 2m or by a totient alone produces those class counts.


## Independent audit

A separately written symmetric-polynomial implementation reproduced all ordinary
and twisted recovery coefficients through 15. It independently solved pair-orbit
orientation constraints and tested modules for every permutation cycle type
through seven, including forbidden types. Its independent full-cycle checks
through 18 agree with the displayed Dirichlet boundary. The combined public
checker imports only the two companion recovery modules and the standard library;
it needs no saved receipt or external classification data.
