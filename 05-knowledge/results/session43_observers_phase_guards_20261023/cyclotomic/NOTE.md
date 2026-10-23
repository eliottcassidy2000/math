# Session 43: cyclotomic divisibility from directed path covers

Status: all-size proofs below, with exact finite checks. This packet does not claim a full metagraph wiring theorem, a Hamiltonian-count spectrum theorem, or a Gaussian moment/Jacobian implication. The binary exponent in Section 4 is proved optimal for every odd order; no corresponding universal odd-prime optimality is claimed.

## 1. Carrier and inherited identities

For a tournament T on n >= 1 vertices, let F_T(x) sum x to the number of forward consecutive arcs over all n! vertex orderings. Let u_m(T) count unordered covers of its vertices by m nonempty directed paths, including one-vertex paths. A cover retains each path's direction but has no ordering of the paths.

THM-4756 proves

    F_T(x) = sum_{m=1}^n m! u_m(T) (x-1)^(n-m).                (1)

A direct proof marks a subset of the forward adjacency positions in each vertex ordering. The marked runs form an ordered path cover. In particular all factorial divisibilities here concern coefficients of this polynomial, not Hamiltonian-path counts alone.

Reversing a vertex ordering proves the reciprocal identity

    F_T(x) = x^(n-1) F_T(1/x).                               (2)

This reversal does not require T to be self-converse. Converse T itself has the same F polynomial.

## 2. All-prime cyclotomic digit bound

Fix a prime p, a primitive root zeta=zeta_p, and pi=zeta-1. Normalize v_pi(pi)=1. In Z[zeta], v_pi(p)=p-1. Define

    M_p(n) = max_{1 <= m <= n} s_p(m),
    E_p(n) = n - M_p(n),

where s_p is base-p digit sum. Legendre's factorial valuation gives

    v_pi(m! pi^(n-m))
      = (p-1) v_p(m!) + n-m
      = m-s_p(m)+n-m = n-s_p(m).

Every summand of (1) belongs to pi^E_p(n) Z[zeta], hence so does F_T(zeta). This proof is independent of any enumeration and holds for all tournaments and all n.

### Odd-prime refinement from reflection

For every odd prime p, a nonzero F_T(zeta) has EVEN pi-valuation. Indeed choose a with 2a congruent to n-1 modulo p. Then beta=zeta^(-a)F_T(zeta) is fixed by complex conjugation by (2). Write beta=pi^v u in the local ring, with u a unit. Conjugation acts trivially on the residue field F_p and bar(pi)/pi=-zeta^(-1); therefore reduction of bar(beta)/beta gives 1=(-1)^v in F_p. Since p is odd, v is even.

Consequently

    F_T(zeta_p) in pi^[2 ceil(E_p(n)/2)] Z[zeta_p]   (p odd). (3)

The zero value satisfies all such divisibilities. Equivalently this parity follows from the ramification-index-two extension of the maximal real cyclotomic subfield. The elementary argument above avoids requiring that formulation.

For p=2, zeta=-1 and (1) gives 2^E_2(n) | F_T(-1). For even n, (2) gives F_T(-1)=0 exactly. A stronger odd-order binary bound is proved in Section 4.

### Ternary scalar specialization

Let S_j=sum_{k congruent j (mod 3)} [x^k]F_T(x). Put r=(2(n-1)) mod 3, so 2r congruent n-1 modulo 3. The two other S_j are equal by (2). Thus

    F_T(omega) = omega^r a_T,
    a_T = S_r - S_j for either j != r,
    3^ceil(E_3(n)/2) divides a_T.                            (4)

Here a_T is a genuine rational integer. The analogous real-subfield coefficient for p>=5 is not generally a rational integer; only the valuation conclusion (3) transfers.

This recovers universal 3 divisibility at orders 3,4,5; 9 at orders 6,7,8; 27 at orders 9,10; and 81 at order 11. THM-085 already proved the initial mod9 result; the digit-sum bound packages an all-size hierarchy rather than replacing that prior theorem.

## 3. Recovered odd-cycle Eulerian expansion (THM-114)

The stronger input already in THM-4756 is the full principal-specialization identity

    U_T(1^k) = sum_R (2k)^r k^(n-|V(R)|),                    (5)

where R runs over collections of vertex-disjoint directed odd cycles of length at least 3, r=|R|, and the empty collection is included. Cycles are counted modulo cyclic rotation, not modulo reversal. On the path-cover side,

    U_T(1^k) = sum_m u_m(T) (k)_m,

because an unordered m-path cover can be placed in distinct members of k ordered slots in (k)_m ways. For the usual Eulerian polynomial A_d, the Stirling identity is

    A_d(x) = sum_m m! S(d,m) (x-1)^(d-m).

Expand each k^d in (5) in falling factorials and apply (1). With d=n-|V(R)|+r, this proves the exact all-size identity

    F_T(x) = sum_R 2^r (x-1)^(n-d) A_d(x).                  (6)

For n>=1, always d>=1. Moreover n-d is even, since every cycle has odd length. THM-114 already states (6) in its two-parameter master-polynomial notation; the derivation here independently recovers it from THM-4756. Formula (6) is the precise structural bridge from odd-cycle collections to forward-step observables. It is a direct reformulation of the inherited principal specialization, not a claim of new discovery of the Redei-Berge identity.

## 4. Stronger binary bound via tangent valuations

For odd d, A_d(-1) is the signed tangent number and

    v_2(A_d(-1)) = d - v_2(d+1).                            (7)

For completeness, write d=2j-1. The absolute tangent number is

    2^(2j)(2^(2j)-1) |B_(2j)| /(2j).

The von Staudt-Clausen denominator theorem gives v_2(B_(2j))=-1, proving (7). These classical facts are recorded by NIST DLMF, Sections 24.15 and 24.10:
https://dlmf.nist.gov/24.15
https://dlmf.nist.gov/24.10

For odd n every d in (6) is odd. Its R-summand has 2-adic valuation

    r+n-d + d-v_2(d+1) = n+r-v_2(d+1).

The empty collection gives n-v_2(n+1). For a nonempty collection, r>=1 and d<=n-2r<=n-2, whence v_2(d+1)<=floor(log_2(n-1)). Therefore, for odd n>=3,

    2^B(n) divides F_T(-1),
    B(n)=min(n-v_2(n+1), n+1-floor(log_2(n-1))).             (8)

At n=1 the value is 1. At n=3,5,7,9,11,13,15 the exponents B are 1,4,4,7,9,11,11. This improves the generic digit-sum bound, for example from 8 to 16 at n=5 and from 64 to 128 at n=9.

At n=5, (6) gives the especially transparent identity

    F_T(-1) = 16 - 16 c_3(T) + 32 c_5(T),                  (9)

where c_l counts directed l-cycles. This explains the finite scan's universal divisibility by 16 without enumeration. The normalized identity 1-c_3+2c_5 was already recorded in THM-061, via THM-059; equation (9) is not claimed as new.

### All-size optimality, with an explicit extremizer

In fact B(n) is the EXACT minimum finite 2-adic valuation of F_T(-1) among all n-vertex tournaments, for every odd n>=3. This means the 2-primary part of the universal gcd, not that its odd part is determined.

Write a=floor(log_2(n-1)). If n-v_2(n+1)<=n+1-a, the transitive tournament attains B(n) by (7).

Otherwise put q=n-2^a+2. This is an odd integer at least 3. On q ordered vertices orient every arc increasingly except reverse the arc between first and last. Append n-q vertices as an ordinal sum, all cross-arcs pointing forward. All directed cycles lie in the q-block and contain the reversed arc; consequently there are no two vertex-disjoint directed cycles. The unique q-cycle visits its vertices in increasing order and returns via the reversed arc. For that cycle, d=n-q+1=2^a-1, so its summand in (6) has valuation n+1-a=B(n). Every shorter cycle has 2^a<d+1<=n-1<2^(a+1), so its summand has strictly larger valuation. The empty collection also has strictly larger valuation by the case assumption. A unique lowest-valuation summand cannot cancel, proving attainment at every odd n.

Equivalently, the sharp profile has a binary-repunit boundary:

    B(n)=n-log_2(n+1)                 if n+1 is a power of 2,
    B(n)=n+1-floor(log_2(n-1))        otherwise (odd n>=3).

For the second branch, n lies strictly between consecutive binary repunits, so v_2(n+1)<=floor(log_2(n-1))-1. The repunit orders 3,7,15,31,63,127,... therefore have a precise role in this valuation observable. This does not identify them with forbidden Hamiltonian-path counts. The transitive tournament attains B(n) exactly at n=2^r-1 or n=3*2^j-1 (j>=1), besides the separately handled n=1 boundary. To see the non-repunit classification, write a=floor(log_2(n-1)); equality requires v_2(n+1)=a-1. In the interval 2^a<n+1<2^(a+1), the only integer with this valuation is 3*2^(a-1). This is a typed occurrence of the 3*2^j-1 family in the anti-evaluation problem, not a Collatz descent assertion. This refinement was supplied independently by the root reviewer.

More explicitly, the reversed-arc block has exactly binomial(q-2,l-2) directed l-cycles: choose the intermediate vertices between its first and last. This gives an independent closed scalar formula for each extremizer. The checker verifies it through order 201 and compares it to actual tournament subset-DP values through order 9. These checks support, but are not the basis of, the all-size proof.

## 5. Exact Gaussian radial transform, and its boundary

Let Z have density exp(-|z|^2)/pi on C; then |Z|^2 has exponential density exp(-t) on t>=0 and E|Z|^(2m)=m!. Consequently

    F_T(x) = E[sum_m u_m(T)|Z|^(2m)(x-1)^(n-m)]
           = integral_0^infinity exp(-t)
                sum_m u_m(T)t^m(x-1)^(n-m) dt.             (10)

The finite sum justifies interchange immediately. This is a positive radial factorial-moment realization of (1). THM-4757 already supplies a related gamma-integral representation of the signed linear-forest polynomial; THM-4001 and THM-3300 already develop factorial/Gaussian machinery on other carriers. We claim a transparent path-cover specialization, not novelty of the factorial integral itself.

This is NOT the vanishing of all powers of one fixed polynomial under a Gaussian functional, and gives no Gaussian moment conjecture or Jacobian conjecture conclusion. The input polynomial in (10) depends on T and x and changes when n changes.

## 6. Verification and priority boundaries

`check.py` uses exact subset dynamic programming for F and explicit directed odd-cycle collection enumeration for (6). It checks all 1+2+8+64+1024=1099 labelled tournaments of orders 1 through 5, plus 102 deterministic distinct samples at each order 6 through 10. Its odd-cycle equality checks cover all orders <=5 and all 102 order-6 samples. Larger orders use exact polynomial evaluation and divisibility checks, not cycle enumeration. All checks use explicit exception gates and remain enabled with Python -O.

At p=2,3,5,7 it reduces F exactly in the basis 1,zeta,...,zeta^(p-2) and divides by zeta-1 using integer arithmetic. For a coefficient vector a, division is possible iff sum(a) is divisible by p. Writing t=-sum(a)/p, quotient coefficient b_j equals -sum_(i<=j)a_i-(j+1)t. The checker separately reconstructs this formula on 500 arbitrary exact products for primes 2,3,5,7,11. Zero is recorded separately from finite valuation. The explicit sharp binary family is checked by its exact scalar formula at all 100 odd orders from 3 through 201 and against actual tournament DP through order 9.

Normal and optimized output must match byte-for-byte. The checked lower bounds pass. The ternary rounded bound is attained in every tested order 1 through 10; this is an exhaustive statement only for n<=5 and a witnessed existence statement for n>=6. Binary E_2 is not uniformly sharp even at order5. Do not transfer any of these statements to the native Hamiltonian-port/electrical mod8 guard, to the 8/9 spectral coupling gap, or to H(T) itself.

Relevant source inspection: THM-059, THM-061, THM-062, THM-085, THM-086, THM-114, THM-4756, THM-4757, THM-4001, and the older `cyclotomic_tournament_deep.py` (which studies a different independence-polynomial evaluation carrier). A theorem-directory digit-sum/cyclotomic-valuation search found other factorial valuation contexts but not this tournament statement. This is a bounded repository priority check, not a claim of external literature novelty.
