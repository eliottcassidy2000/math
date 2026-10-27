# Finite-section transfer operators of the Syracuse map: rank one plus nilpotent, second eigenvalue exactly 0 for every K and J, singular gap 2^-(K-1), exact Jordan types, blind to the multiplier

**opus-2026-10-27-S65.**

- Script: `04-computation/experiments/syracuse_transfer_finite_section_20261027.py` (+ `.out`, ALL CHECKS PASSED, about
  75 seconds).
- Canon: THM-4821 (the 2-adic carrier) and THM-4822 (the 3-adic / p-adic carrier), both filed this session.
- Related canon:
  - THM-4629 proof 2(c)–(d): the 2-adic de Bruijn part of the inverse-tree operator. 2(c): the halving-branch averaging
    E_j is rank one plus nilpotent. 2(d): the full de Bruijn walk has the single nonzero eigenvalue (u + v)/2.
  - THM-4631: the twisted level factorisation; the "fast modes" are the doubling cycles.
  - THM-4520 and THM-4519: the Fourier-side level operator and the Syracuse frequency recursion.
- Small cases already in the repository (§6):
  - `procgen_repunit_20260924_digits_rotation_repunits.md` §2.1 (Theorem H): the mod-8 and mod-9 Syracuse chains and the
    T-map window chain.
  - The TCPC note `procgen_tcpc_20261001_tournament_clock_prime_collatz.md`: §5.2 has Tao's mod-9 law and the mod-9
    chain mixed exactly in two steps; §5.4 has the square-clock IFS.

**Status.**
- **The owner's question has a sharp answer. The eigenvalue spectral gap is exactly 1 for every K and every J.** The
  second-largest eigenvalue modulus (SLEM) of both finite sections is exactly 0, so the gap is uniform and does not decay.
  - PROVED: P_K = Π + N_K with N_K nilpotent of index exactly K − 1 (2-adic carrier, THM-4821).
  - PROVED: P_J = Π + N_J with N_J nilpotent of index exactly J (Tao's chain mod 3^J, THM-4822).
  - The certificate is the exact identity (P − Π)^(K−1) = 0, equivalently Dobrushin's coefficient δ(P^(K−1)) = 0.
  - It is verified exactly for K = 8..20.
  - For Tao's chain it is verified exactly for J = 1..4 from every start, and for J = 5, 6 from five starts (from every
    start modulo a prime here, and exactly in the independent audit). For J = 7, 8 it is verified modulo a prime.
- **What does grow or decay with K is non-spectral** (PROVED, 2-adic).
  - The nilpotent part has a Jordan block of size K − 1.
  - Exact equilibrium is reached after h(a) steps, with h − 1 ~ Binomial(K − 2, 1/2) over uniform starts. So the mean is
    K/2, the worst case is K − 1, and the worst-case total variation is 1 − 2^-(K−1−n).
  - The second singular value is √(1 − 2^-(K−2)), so the singular-value gap is ≈ 2^-(K−1).
  - More precisely, ‖(P − Π)^n‖² = P(Bin(K − 2, 1/2) ≥ n) in ℓ²(π). The full singular spectrum of P^n is 1 together
    with √P(Bin(K − j, 1/2) ≥ n) of multiplicity 2^(j−2), for j = 2..K.
  - P_K is a maximally non-normal "rank one plus nilpotent" operator. The eigenvalue gap overstates its mixing, and the
    singular gap understates it.
- **Exact Jordan types (PROVED).**
  - 2-adic: one block of size K − 1 and 2^(K−2−s) blocks of size s (1 ≤ s ≤ K − 2).
  - p-adic (Tao's chain for p = 3): p − 1 blocks of size J and (p − 1)² p^(J−1−s) blocks of size s. This counts all of
    Z/p^J, transient states included. On the recurrent class the type is in §3, and the index is J in both.
  - These come from rank P_K^n = 2^(K−1−n) and rank P_J^n = p^(J−n).
- **Blind to the multiplier (PROVED; the 2-adic case is classical).**
  - For every odd p, the Ulam matrix of x ↦ (px + 1)/2^v is permutation-conjugate to the one for p = 3. The relabelling is
    the classical parity-vector conjugacy (Terras 1976; Lagarias 1985), checked entry by entry for eight multipliers.
  - The chain x ↦ (px + 1)2^-a mod p^J satisfies the same statements with p in the counts (checked for p = 5, 7).
  - So these finite sections cannot tell 3x + 1 from 5x + 1. The multiplier enters the integer problem through the
    archimedean drift log p − 2 log 2, which no 2-adic Haar section sees.
  - The arithmetic of Collatz lives elsewhere: in the stationary law (Syrac, whose Fourier decay is Tao's input and
    THM-4519's H1) and in the twisted operators of THM-4629/4631.
- Contrast (PROVED, = canon): in the T-map convention (x/2 or (3x+1)/2 with fair coins) the 3-adic chain has SLEM
  exactly 1/2 for every J. These are the doubling cycles of THM-4631 with j = 0 on the line u = v = z, where
  det(I − zU) = (1 − z) Π_d (1 − (z/2)^(L_d)). That gap is also uniform.
- Floating point, as a guide only: numpy reports |λ₂| = 0.0029 to 0.0179 for K = 8..12, growing with K. These are
  pseudo-spectral artefacts of the Jordan block of size K − 1, and the exact value is 0.
- **Independent audit: SOUND WITH CORRECTIONS (§8). All corrections are applied.** There was one error, in this note's
  display of the law of Π (Theorem B(1); MISTAKE-686). The audit also found in-repository prior art for the smallest
  cases (§6).

## 0. The request

"Build the finite-section transfer operator of the Syracuse map S(x) = (3x+1)/2^v on odd residues modulo 2^K (and, as a
second carrier, Tao's Syracuse random variables modulo 3^J), for K = 8..20 (J as far as feasible), as an exact rational
stochastic matrix (or the exact integer counting matrix), and compute its second-largest eigenvalue modulus or an exact
certified bound on it (Gershgorin/Collatz–Wielandt style bounds with exact rationals are acceptable; floating point only
as a guide, never as the result). Report whether the spectral gap is bounded away from 1 uniformly in K or decays with K."

## 1. The two carriers, exactly

**The 2-adic carrier (the Ulam matrix of S).**
- Let Ω_K be the 2^(K−1) odd residues modulo 2^K, and π the uniform law on Ω_K. Define

  P_K(a, b) = Haar{x ∈ a + 2^K Z₂ : S(x) ∈ b + 2^K Z₂} / Haar(a + 2^K Z₂).

  This is Ulam's finite section of the transfer operator of S on the odd 2-adic integers with Haar measure (which S
  preserves).
- Concretely, put v = v₂(3a + 1).
  - If v < K, row a is uniform on the 2^v residues b ≡ (3a + 1)/2^v (mod 2^(K−v)).
  - The one residue a* = −1/3 mod 2^K (v ≥ K) has the uniform row.
- The **integer counting matrix** is C_K = 2^(K−1) P_K. For every row except a*, C_K(a, b) = #{t mod 2^(K−1) : S(a + 2^K t)
  ≡ b (mod 2^K)} exactly, because only t mod 2^v matters. The row a* is the Haar limit.
- Checks:
  - The rows agree with direct Haar lift counts (lifts modulo 2^(K+12)) for K = 3..7 and p = 3, 5, 7, 1, −1 (script
    section I). For a*, the 2^(K−1) lifts with v > 12 are set aside, and the remaining lifts are exactly uniform. The
    independent audit reproduced this for K = 2..7 and seven multipliers.
  - Row sums and column sums are 2^(K−1), so P_K is doubly stochastic.
  - The number of nonzero entries is exactly K·2^(K−1), which is 10,485,760 at K = 20.
- A matrix built from integer representatives (b = S(a) mod 2^K for a < 2^K) is not a transfer operator: it depends on
  the choice of representative and is not a Markov compression of the dynamics. It is not used.

**The 3-adic carrier (Tao's chain).**
- On Z/3^J take x ↦ (3x + 1)·2^-a with P(a) = 2^-a, a ≥ 1.
- With L = 2·3^(J−1) = ord(2 mod 3^J), the transition is P_J(x, y) = 2^(L−a₀)/(2^L − 1), where 2^(a₀) ≡ (3x + 1)/y and
  a₀ ∈ [1, L] (y a unit). The exact integer matrix is (2^L − 1)P_J.
- Its stationary law is the law of Tao's Syrac(Z/3^J) (Tao, Forum Math. Pi 10 (2022), e12, arXiv:1909.03562 (2019),
  Lemma 1.12 and Remark 1.13; CITED, with the numbering checked against arXiv v7 by the audit).
- The script reproduces Tao's mod-9 law (0, 8, 16, 0, 11, 4, 0, 2, 22)/63 exactly.
- It also reproduces THM-4520 (5)'s table of the largest primitive Fourier coefficient: 0.37792, 0.25224, 0.17700,
  0.12927, 0.09611 at J = 2..6. So this is the repository's Syracuse measure.

## 2. Theorem A (THM-4821): the 2-adic carrier

Let K ≥ 2, p odd, S_p(x) = (px + 1)/2^(v₂(px+1)), and P = P_K^(p) its Ulam matrix on Ω_K (p = 3 is S).

1. **Class lemma.** For a class C = c + 2^m Z (c odd, 1 ≤ m ≤ K) let μ_C be the uniform law on C ∩ Ω_K. Then
   μ_C P = μ_(Φ(C)), where Φ(c, m) = ((pc + 1)/2^v mod 2^(m−v), m − v) if v = v₂(pc + 1) < m, and Φ(C) = Ω (the root,
   precision 1) otherwise.
2. **Reachable classes and ranks.**
   - For 1 ≤ n ≤ K − 1, the rows of P^n are exactly the laws μ_C over all classes of precision 1..K − n.
   - Each class of precision m ≥ 2 is the row of exactly C(K − m − 1, n − 1) residues: one for each choice of the n − 1
     intermediate precisions strictly between m and K. The number depends only on m and n.
   - rank P^n = 2^(K−1−n).
3. **Spectrum and Jordan type.**
   - P = Π + N with ΠN = NΠ = 0 and N^(K−1) = 0 ≠ N^(K−2).
   - The spectrum is {1} together with the eigenvalue 0 of multiplicity 2^(K−1) − 1, so **SLEM = 0**.
   - On the mean-zero space, N has one Jordan block of size K − 1 and 2^(K−2−s) blocks of size s, for 1 ≤ s ≤ K − 2.
4. **Exact mixing.** Let h(a) be the first n with P^n(a, ·) = π. Then:
   - #{a : h(a) = h} = 2·C(K−2, h−1) for 1 ≤ h ≤ K − 1, so E h = K/2;
   - the worst-case total variation is max_a ‖P^n(a, ·) − π‖_TV = 1 − 2^-(K−1−n) for 0 ≤ n ≤ K − 1;
   - Dobrushin's coefficient δ(P^n) is 1 for n ≤ K − 2 and 0 for n = K − 1.
5. **ℓ² profile and singular spectrum.**
   - ‖(P − Π)^n‖²_(ℓ²(π)) = P(Bin(K−2, 1/2) ≥ n) for 0 ≤ n ≤ K − 1.
   - In particular the second singular value is √(1 − 2^-(K−2)), and the singular-value gap is
     1 − √(1 − 2^-(K−2)) ≈ 2^-(K−1).
   - The singular values of P^n are 1 together with √P(Bin(K − j, 1/2) ≥ n) of multiplicity 2^(j−2), j = 2..K.
6. **The multiplier.** For every odd p, P_K^(p) = Q P_K^(3) Q^(−1) for a permutation Q of the odd residues. Q is induced
   by ψ = Q_3^(−1) ∘ Q_p, where Q_p is the parity-vector map of x ↦ x/2, (px + 1)/2 (a tree automorphism of Z/2^K;
   Terras 1976, Lagarias 1985). Hence 1–5, and every property invariant under relabelling, do not depend on p.

*Proof.*
- **(1)** x ↦ px + 1 maps c + 2^m Z₂ onto pc + 1 + 2^m Z₂ and preserves normalised Haar measure, since p is a unit.
  - If v < m, every element has valuation v, and division by 2^v maps the class onto r + 2^(m−v) Z₂ uniformly.
  - Otherwise the image is all of 2^m Z₂, uniformly. The odd part of a Haar-random element of 2^m Z₂ is Haar-random on Z₂^×:
    condition on the valuation.
  - Projecting modulo 2^K gives μ_(Φ(C)). The script checks the lemma for every class at every level, K = 8..20. At the
    finest level it compares each row with Φ. At coarser levels it checks the identity
    μ_(c,m) = (μ_(c,m+1) + μ_(c+2^m,m+1))/2, which carries the lemma up the tree.
- **(2)** c ↦ pc + 1 is a bijection from the odd residues modulo 2^m onto the even ones.
  - For 1 ≤ v ≤ m − 1, the residues of valuation exactly v divide by 2^v onto all odd residues modulo 2^(m−v). Hence Φ
    maps the 2^(m−1) classes of precision m onto the classes of precision 2..m − 1, hitting each exactly once. The root
    is hit twice: by v = m − 1, and by pc + 1 ≡ 0.
  - Induction on n gives the reachable sets and the level-only multiplicities: mult_(n+1)(m′) = Σ_(m=m′+1..K−n)
    mult_n(m) for m′ ≥ 2, whose solution is C(K − m − 1, n − 1).
  - The indicators of all classes of precision ≤ K − n span the same space as the 2^(K−n−1) disjoint classes of
    precision exactly K − n.
- **(3)** P^n = Π ⊕ N^n on span(1) ⊕ W, where W is the mean-zero space. So rank N^n|_W = 2^(K−1−n) − 1.
  - The number of Jordan blocks of size ≥ s is the difference of consecutive ranks, 2^(K−1−s).
  - rank N^(K−2) = 1, so the index is exactly K − 1.
- **(4)** h(a) − 1 is the number of intermediate precisions visited, a subset T of [2, K − 1].
  - By the bijection in (2) and induction on the top level, each subset T arises from exactly two residues. This gives the
    binomial count.
  - TV(μ_C, π) = 1 − 2^(1−m) for a class of precision m, and the largest precision among the rows of P^n is K − n.
  - Two rows have disjoint supports iff two incomparable classes occur, iff K − n ≥ 2.
- **(5)** For f ∈ W,
  - ‖P^n f‖² = Σ_(m≥2) w_m(n) ‖E_m f‖². Here w_m(n) is the fraction of residues whose row of P^n has precision m, and
    E_m is the conditional expectation onto the classes of precision m. This uses the level-only multiplicities of (2).
  - Each ‖E_m f‖ ≤ ‖f‖, with equality for every m ≥ 2 when f = ±1 on the classes 1 and 3 mod 4. So
    ‖P^n|_W‖² = Σ_(m≥2) w_m(n) = P(h > n), and ‖(P − Π)^n‖ = ‖P^n|_W‖.
- **(5′)** The identity in (5) holds for every f, so (P^n)*P^n = Σ_m w_m(n) E_m as operators. The E_m are nested
  projections, E_1 = Π ≤ E_2 ≤ … ≤ E_K = I.
  - On D_j = range E_j ⊖ range E_(j−1), which has dimension 2^(j−2), the operator acts as the scalar
    τ_j(n) = Σ_(m≥j) w_m(n).
  - τ_j(n) is the fraction of residues still at precision ≥ j after n steps. The set T of intermediate precisions is
    uniform on the subsets of [2, K − 1], as in (4). The precision after n steps is ≥ j iff |T ∩ [j, K − 1]| ≥ n, so
    τ_j(n) = P(Bin(K − j, 1/2) ≥ n).
  - Hence the singular values of P^n are 1 and √P(Bin(K − j, 1/2) ≥ n) with multiplicity 2^(j−2), j = 2..K.
  - Those with K − j < n vanish: 2^(K−1) − 2^(K−1−n) of them, matching rank P^n = 2^(K−1−n). ∎
- **(6)** Let Q_p : Z/2^K → {0, 1}^K send x to its first K parities under T_p(x) = x/2, (px + 1)/2. It is a bijection,
  compatible across K, for every odd p (Terras's argument).
  - So ψ = Q_3^(−1) ∘ Q_p is a Haar-preserving tree automorphism that conjugates T_p to T_3.
  - S_p is the first return of T_p to the odd residues: in parity coordinates, "shift past the next 1", which is the same
    map for every p. Hence ψ conjugates S_p to S_3.
  - It maps each class mod 2^K onto a class mod 2^K, so P_K^(p)(a, b) = P_K^(3)(ψa, ψb).
  - The special residue −1/p goes to −1/3: both have the code 1 followed by K − 1 zeros. ∎

*Checks (FINITE-EXACT, K = 8..20).*
- The class lemma at every level.
- The reachable sets: 2^(K−n) − 1 classes, all precisions 1..K − n, with multiplicities C(K − m − 1, n − 1).
- The ranks, from an exact bottom-up span count on the binary tree of classes.
- The Jordan type, the hitting-time counts, E h = K/2, the worst-case total variation and the Dobrushin sequence.
- Every tail weight τ_j(n) = P(Bin(K − j, 1/2) ≥ n) of (5′), in exact rationals. τ_2 is ‖(P − Π)^n‖² at the maximiser.
- The conjugacy of (6), entry by entry on the counting matrices:
  - p = 1, 5, 7, 9, 15, −1, −3, 1023 at K = 8..14;
  - p = 5 and −1 at K = 15..20.

Independent of the class machinery:
- For K = 8..12, C^(K−1) = 2^((K−1)(K−2)) J. It is proved by congruence modulo 3 to 7 primes whose product exceeds the
  entry bound 2^((K−1)²+1); K = 8 is also done in exact integers.
- C^(K−2) differs from 2^((K−1)(K−3)) J modulo the prime 1048573. That is a rigorous certificate that P^(K−2) ≠ Π.
- The float SVD of every power agrees with the full singular spectrum of (5) for K ≤ 11, and with its norm at K = 12.
- The reachable classes, ranks and hitting times for p = 1, 5, 7, −1 are checked at K = 8..14 (section H).
- The independent audit recomputed the following with exact dense integer matrices (python-flint) for K ≤ 11:
  - the minimal polynomial (x − 2^(K−1)) x^(K−1) of C, and (C − J)^(K−1) = 0 ≠ (C − J)^(K−2);
  - the exact Gram identity (P^n)ᵀP^n = Σ_m w_m(n) E_m for K ≤ 8.

| K | states | nnz | index | Jordan blocks of sizes 1, 2, 3, 4 | ‖P − Π‖₂² | singular gap | numpy \|λ₂\| (exact: 0) |
|---|---|---|---|---|---|---|---|
| 8 | 128 | 1,024 | 7 | 32, 16, 8, 4 | 63/64 | 7.84e-3 | 0.0029 |
| 12 | 2,048 | 24,576 | 11 | 512, 256, 128, 64 | 1023/1024 | 4.88e-4 | 0.0179 |
| 16 | 32,768 | 524,288 | 15 | 8192, 4096, 2048, 1024 | 16383/16384 | 3.05e-5 | — |
| 20 | 524,288 | 10,485,760 | 19 | 131072, 65536, 32768, 16384 | 262143/262144 | 1.91e-6 | — |

(Every K from 8 to 20 is in the `.out`.) For K = 12, the exact counts of h = 1..11 are 2, 20, 90, 240, 420, 504, 420, 240,
90, 20, 2. The average total variation profile is 0.9995, 0.9941, 0.9673, 0.8867, 0.7256, 0.5, 0.2744, 0.1133, 0.0327,
0.0059, 0.0005, 0. This is a binomial cutoff. The average total variation is exactly P(B > n) + P(B = n)/2 with
B ~ Bin(K − 2, 1/2), which equals 1/2 at n = (K − 2)/2, with a window of order √K. The mean hitting time is K/2.

## 3. Theorem B (THM-4822): the 3-adic carrier and its p-adic analogues

Let p ≥ 3 be odd and J ≥ 1, and let P = P_(p,J) be the chain x ↦ (px + 1)·2^-a on Z/p^J, a ~ Geom(1/2). For p = 3 this
is Tao's chain.
1. **Forgetting.** x_n = p^n x₀ 2^-(A_n) + Σ_(j=1..n) p^(j−1) 2^-(a_(n−j+1) + ... + a_n). So P^n(x₀, ·) depends only on
   x₀ mod p^(J−n), and **P^J = Π**. The rows of Π are the law of Σ_(j=1..J) p^(j−1) 2^-(a_(J−j+1) + ... + a_J).
   Reversing the i.i.d. a's, this is the law of Σ_(j=1..J) p^(j−1) 2^-(a_1 + ... + a_j). For p = 3 it is Syrac(Z/3^J),
   in the form of Tao's (1.29).
2. **Ranks.** rank P^n = p^(J−n) for 0 ≤ n ≤ J.
3. **Spectrum and Jordan type.**
   - P = Π + N, N nilpotent of index exactly J. The spectrum is {1} with 0 otherwise, so **SLEM = 0**.
   - On the mean-zero space, N has p − 1 blocks of size J and (p − 1)² p^(J−1−s) blocks of size s, for 1 ≤ s ≤ J − 1.
   - For p = 3 that is 2 blocks of size J and 4·3^(J−1−s) blocks of size s.
   - Every start needs exactly J steps: no row of P^(J−1) equals π.
4. **The recurrent class.** The type in 3 counts the chain on all of Z/p^J, including the p^J − r p^(J−1) transient
   states. Let R = supp π = {y : y mod p ∈ ⟨2⟩}, where r = |⟨2 mod p⟩| ≥ 2 (for p = 3, R is the units).
   - rank P_R^n = r p^(J−1−n) for 0 ≤ n ≤ J − 1, and 1 at n = J.
   - For J ≥ 2, the Jordan type of N on the π-mean-zero functions on R is r − 1 blocks of size J, r(p − 2) + 1 blocks of
     size J − 1, and r(p − 1)² p^(J−2−s) blocks of size s ≤ J − 2.
   - For p = 3 this is 1, 3 and 8·3^(J−2−s) blocks. The index is J here too.

*Proof.* (1) is induction on n.
- It is Tao's (1.21)/(1.24): Syr^n(N) mod 3^k depends only on the valuations for k ≤ n.
- The case J = 2 (P² = Π mod 9) is in the TCPC note §5.2 and in procgen_repunit_20260924 §2.1.

For (2), let F_d be the functions of y mod p^d and U the Koopman operator.
- **U F_d ⊆ F_(d−1).** (px + 1)2^-a mod p^d depends only on x mod p^(d−1).
- **U is onto F_(d−1).** Write the map F_d → F_(d−1) as a matrix whose rows are indexed by x̄ ∈ Z/p^(d−1).
  - Row x̄ is supported on the coset (px + 1)⟨2⟩ of the units mod p^d.
  - Along that cycle (length L = ord(2 mod p^d)) it carries the weights w(j) = 2^-j/(1 − 2^-L), j = 1..L, shifted to
    start at px + 1.
  - Rows on different cosets have disjoint supports.
  - Rows on the same coset are distinct shifts of w, because x̄ ↦ px + 1 is injective. They are therefore distinct rows
    of the circulant Circ(w).
  - That circulant is invertible: its eigenvalues are 1/(2ζ − 1) with ζ^L = 1, never 0 (the circle |z − 1/3| = 2/3 of
    THM-4520 (2)).
  - So the rows are independent, and U is onto F_(d−1).
- Hence U^n maps F_J onto F_(J−n), and rank P^n = rank U^n = p^(J−n).

(3) follows from the ranks as in Theorem A. The index is exactly J because rank N^(J−1)|_W = p − 1 > 0.
- The p rows of P^(J−1) (indexed by x mod p) are independent, and π = πP^(J−1) = Σ_c π_c·row_c, where
  π_c = π{x ≡ c mod p}.
- π_c > 0 for each of the r ≥ 2 classes c ∈ ⟨2⟩, since x mod p = 2^(−a₁) at stationarity. So π is no single row.

(4) Every image (px + 1)2^-a is ≡ 2^-a mod p, so R is closed and receives every step. P_R^n is the restriction of P^n to
R × R.
- Let F_d(R) be the functions on R of y mod p^d. Its dimension is r p^(d−1) for d ≥ 1, and 1 for d = 0.
- U maps F_d(R) onto F_(d−1)(R). Given g ∈ F_(d−1)(R), extend it arbitrarily off R, take a preimage f ∈ F_d under the
  surjection of (2), and restrict f to R. This works because U f(x) only reads f on R.
- So rank P_R^n = dim F_(J−n)(R).
- The number of blocks of size ≥ s is rank N_R^(s−1) − rank N_R^s. This is r − 1 for s = J, and r(p − 1)p^(J−1−s) for
  s ≤ J − 1. Consecutive differences give the type. ∎

*Checks.*
- p = 3, exact rational arithmetic:
  - P^J has equal rows from every start for J ≤ 4, and from five starts for J = 5, 6;
  - P^(J−1) does not, for J ≤ 6;
  - the stationary laws, Tao's mod-9 law, and THM-4520's Fourier table at J = 2..6.
- p = 3, modulo the prime q = 1048573:
  - P^J has equal rows from every start and P^(J−1) differs, for J ≤ 8. The difference mod q is a rigorous certificate
    that P^(J−1) ≠ Π.
  - rank_q P^n = 3^(J−n) for J ≤ 6, and 2·3^(J−1−n) on the units. A modular rank is a lower bound for the rational rank,
    and it meets the structural upper bound, so the ranks are exact.
- p = 5 and p = 7 (J ≤ 3): the same checks, including the ranks on R (r = 4 and r = 3). For p = 7, 2 is not a primitive
  root (ord 3, 21, 147), and the cosets matter.
- The independent audit, in exact integers (python-flint), also covered:
  - J = 5 and 6 from every start;
  - p = 9, 11, 13, 15, 17, 21, 23, 31 at J = 2 (primality of p is not used);
  - the type of (4) in 12 cases, including p = 7, 17, 31, 9, 15;
  - THM-4520's Fourier value at J = 7 (0.07587).
- Float guide: the ℓ²(π) norm of P_J − Π on the units is 0.674, 0.971, 0.996, 0.9995 and 0.99995 at J = 2..6. The same
  non-normality is visible, but it is not proved here.

**Fourier dictionary.** In the character basis graded by 3-adic depth, U lowers depth by one at each step, which makes
the nilpotency visible. The blocks from depth d to depth d − 1 carry the Gauss phases e(2^k/3^d) of THM-4520. The
characteristic function of Syrac(Z/3^J) at a primitive frequency is the product of the J level blocks applied to e_ξ:
the cocycle of THM-4519/4520 (their hypothesis H1). The spectral data of P_J are trivial. Everything arithmetic sits in
that product, that is in the entries of the rank-one limit Π, and not in any eigenvalue.

## 4. Answer, and what these operators can and cannot see

| operator | carrier | SLEM (exact) | uniform in the level? | exact time to equilibrium | second singular value |
|---|---|---|---|---|---|
| Syracuse Ulam P_K | odd residues mod 2^K | 0 | yes | ≤ K − 1 (mean K/2) | √(1 − 2^-(K−2)) → 1 (proved) |
| Tao's chain P_J | Z/3^J | 0 | yes | J, from every start | → 1 numerically (J ≤ 6) |
| T-map Ulam (de Bruijn) | Z/2^K | 0 | yes | K | 1 exactly: P Pᵀ is a sum of 2 × 2 averaging blocks. Spectrum: THM-4629 proof 2(d) at u = v = 1. Index K: classical (de Bruijn A^K = J), in procgen_repunit_20260924 §2.1 |
| T-map random chain | Z/3^J | 1/2 | yes | never exact | — (THM-4631 with j = 0 on the line u = v = z) |

- **The answer.** The spectral gap is bounded away from 0 uniformly, and in the Syracuse normalisation it is exactly 1
  for every K and J. It does not decay.
- **What decays** is the singular-value gap: 2^-(K−1) for the 2-adic carrier, proved. **What grows** is the Jordan index,
  and with it the time to equilibrium, linearly in K. That growth is the true mixing scale.
- **The finite sections are blind to the multiplier.**
  - The 2-adic statements hold verbatim for (px + 1)/2^v with any odd p. The p-adic chain satisfies the same statements
    with p in the counts. Its state space and its stationary law (where the arithmetic lives) do depend on p.
  - In the 2-adic case the finite sections for different p are the same Markov chain up to relabelling (Theorem A (6)).
    This is the finite-section shadow of the classical fact that every x ↦ x/2, (px + 1)/2 is the same Bernoulli shift
    under Haar measure.
  - So no property of these untwisted operators, spectral or not, can separate 3x + 1 (convergent heuristics) from
    5x + 1 (divergent heuristics). The multiplier enters the integer problem through the archimedean drift
    E log(p 2^-v) = log p − 2 log 2, which is negative for p = 3 and positive for p = 5, and which no 2-adic Haar section
    sees.
  - This is the forward-operator counterpart of THM-4629 statement 2: every finite deterministic truncation of the
    inverse-tree hierarchy is silent. THM-4631 already notes the T-map analogue: p enters only through the odd branch.
- **Where the arithmetic lives.**
  - In the stationary law: the Fourier decay of Syrac(Z/3^J), Tao's fine-scale mixing input, and THM-4519's H1.
  - In twisted operators: THM-4629/4631, where eigenvalues of modulus 1/2 appear on the 3-adic levels and the twist
    detects near-returns.
  - In the passage from Haar-random to integer starting points, which the Haar finite sections do not see at all.
- **Why the T-map 3-adic chain has a nonzero spectrum.** Its halving branch preserves 3-adic depth: it permutes the units
  by ×2^-1. The Syracuse step bundles every halving with one ×3, so every step lowers depth.

## 5. Floating point

numpy's `eigvals` on P_K (K = 8..12) returns |λ₂| = 0.0029, 0.0046, 0.0089, 0.0125, 0.0179.
- These are for the natural ordering. The values depend on the ordering and on the LAPACK build: under random
  relabellings the audit saw 0.0023 to 0.0030 at K = 8 and 0.016 to 0.022 at K = 12.
- On Tao's chain (J = 2..6) it returns up to 0.0013.
- The exact value is 0 in every case.

A nilpotent Jordan block of size m perturbed at relative size ε has eigenvalues of modulus about ε^(1/m). That is 0.005
to 0.035 for m = K − 1 = 7..11 and ε = 10^-16; the observed values are about half of that, the right order and growth.
The estimate grows toward 1 as m grows. Read naively, floating point would suggest a gap that decays with K. Here it is
exactly the artefact the request warned against.

## 6. Prior art and typing

- **CITED (classical).**
  - Under Haar measure the 2-adic Collatz map is conjugate to a one-sided Bernoulli shift by the parity-vector map
    (Terras 1976; Lagarias 1985, Amer. Math. Monthly 92, 3–23). The same argument works for every odd multiplier. From
    memory, not re-read; the audit confirmed the substance.
  - The Syracuse valuations are i.i.d. geometric (Tao, Forum Math. Pi 10 (2022), e12 [arXiv 2019], §1.3,
    Proposition 1.9).
  - Tao's (1.21)/(1.24) is Theorem B(1).
  - Theorem A is the finite-section shadow of these facts. P_K is exactly Ulam's matrix
    P(a, b) = m(A_a ∩ S^(−1)A_b)/m(A_a).
- **Canon.** THM-4629 2(c)–(d) and THM-4631 already prove the rank-one-plus-nilpotent structure of the 2-adic de Bruijn
  averaging in the inverse-tree (T-map) operator. THM-4631 also finds the doubling-cycle modes of the 3-adic levels. The
  T-map rows of §4 are those results restated (checked again here).
- **Small cases already in the repository** (found by the audit; missed in the pre-audit draft, MISTAKE-686):
  - `procgen_repunit_20260924_digits_rotation_repunits.md` §2.1 (Theorem H) has three of them:
    - the mod-8 Syracuse chain with P² uniform, which is K = 3 of Theorem A;
    - the mod-9 chain with P² of identical rows and Tao's law, which is J = 2 of Theorem B;
    - the T-map window chain mod 2^k, mixed exactly in k steps.
  - The TCPC note §5.2 has the mod-9 rank 2 (on the units), P² = 1π, and "A_(n+1) mod 3^k is a function of the last k
    valuations" (Theorem B(1)).
  - `graph_theory_as_common_carrier_20261013.md` §5 and HYP-9261 (5): the T-map class graph mod 2^j is the de Bruijn
    graph B(2, j), and E_j is rank one plus nilpotent of index j.
- **New here, as far as the repository and these sources go:**
  - the forward Syracuse Ulam operator's exact Jordan type, hitting-time law, total-variation and ℓ² profiles, and its
    full singular spectrum, for general K;
  - the exact ranks and Jordan types of Tao's chain for general J, on Z/p^J and on its recurrent class, via the
    invertible circulant Circ(w).

  The universality in p is a corollary of the classical parity-vector conjugacy: the finite sections for different p are
  the same Markov chain up to relabelling. None of this bears on the Collatz conjecture itself.
- **Possibly relevant, not read** (the audit's pointers, from memory): Matthews and Watts (Acta Arith. 1984, 1985) on
  Markov chains of generalised Syracuse maps modulo m. For the de Bruijn A^K = J: the 0–1 matrices with A^k = J (Hoffman
  1967; Lam and van Lint 1978). The audit's web search found no source for the exact Jordan types, so novelty beyond the
  repository is unestablished.
- NOT CLAIMED:
  - any statement about integer orbits;
  - an exact ℓ² profile for Tao's chain;
  - a literature search beyond the cited sources.

## 7. Reproduction

`python 04-computation/experiments/syracuse_transfer_finite_section_20261027.py` runs in about 75 seconds with numpy,
sympy and Python 3.11+. It prints the checks of sections A–I in the order they appear in this note:
- A–C, C2: the 2-adic theorem, its class multiplicities, ℓ² profile and singular spectrum;
- D: the dense certificates;
- E: the float guide;
- F: Tao's chain, including the ranks on the recurrent class;
- G: the T-map contrasts;
- H: the multiplier p (the parity-vector conjugacy, and the p-adic chains for p = 5, 7);
- I: the rows against direct Haar lift counts.

## 8. Audit record

An independent audit of the pre-audit checkpoint (commit cac59265c5) returned **SOUND WITH CORRECTIONS**.
- Report and scripts: `04-computation/experiments/syracuse_transfer_finite_section_20261027/audit_S65/AUDIT.md`.
- Scope: this note, THM-4821, THM-4822 and the script.
- Method: the auditor's own code, using python-flint exact matrices, sparse class aggregation, mpmath and negative
  controls.
- Both theorems are true as stated, every proof step checks, and every recomputed number agrees. A rerun of the script
  was byte-identical to the committed output after the timing fields were stripped.

Corrections applied:
1. **Error (this note only).** The pre-audit Theorem B(1) displayed the law of Π as Σ_j p^(j−1) 2^-(a_j + ... + a_J).
   - That indexing puts p^(j−1) on the last J − j + 1 valuations instead of the last j.
   - At p = 3, J = 2 its mod-3 marginal is (0, 5/9, 4/9), against Syrac's (0, 1/3, 2/3).
   - It is now Σ_j p^(j−1) 2^-(a_(J−j+1) + ... + a_J), equivalently Tao's (1.29).
   - The closed form for x_n displayed beside it, THM-4822, and the script (which iterates the chain) were right
     throughout. MISTAKE-686.
2. **Prior art.** Credit is now given to:
   - procgen_repunit_20260924 §2.1 (the K = 3, J = 2 and de Bruijn cases);
   - the TCPC note §5.2;
   - graph_theory_as_common_carrier_20261013 §5;
   - Tao's (1.21)/(1.24) for B(1).
3. **Citations.**
   - Tao's journal year is 2022.
   - The de Bruijn spectrum is THM-4629 proof 2(d), not 2(c). The index K is classical.
   - The T-map 3-adic contrast is THM-4631 with j = 0 on the line u = v = z. At u = v = 1 the determinant is just 0.
4. **The multiplier.** "Statements 1–5 do not depend on p" is now the permutation conjugacy P_K^(p) = Q P_K^(3) Q^(−1),
   inherited from the classical parity-vector conjugacy. It is checked entry by entry in the script (section H).
5. **Scope.**
   - The p-adic Jordan type counts transient states. The type on the recurrent class is now Theorem B(4), with a proof
     and a modular-rank check.
   - The average total variation crosses 1/2 at n = (K − 2)/2, not K/2.
   - The numpy |λ₂| values depend on the ordering.
6. **Reproducibility.**
   - The lift-count check is now script section I. It was previously an uncommitted scratch check.
   - The Status line on J = 5, 6 now says five starts exactly, every start modulo q.
   - The docstring "(2P)^K = 2^(K−1) J?" is now "(2P)^K = J".
   - The `.out` line endings are normalised.

The audit's verified strengthenings are added (Theorem A (2), (5), (5′); §4 table): the multiplicities
C(K − m − 1, n − 1), the full singular spectrum of P^n, and σ₂ = 1 for the de Bruijn matrix.
