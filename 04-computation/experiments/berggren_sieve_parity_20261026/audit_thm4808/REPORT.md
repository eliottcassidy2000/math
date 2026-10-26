# Independent audit of THM-4808 and HYP-9418 (2)–(5)

Independent adversarial audit, 2026-10-26, run in the worktree `math-wt-pyth-cuboids`. The lane saved this text from the
auditor's final message, because the auditor could not write it here.
- **Independence.** No lane script was imported; every number below comes from the audit's own code. One lane script
  was re-run, unmodified, for timing only.
- **Text audited.** The audit used the text at commit `4950589d05`, including the lane's two announced sentence edits.
- **Edits.** They are in `required_edits.json`, 34 exact-string replacements validated by `make_required_edits.py`.
  The lane applied all 34.
- **Tools.** PARI/GP 2.17.4, python-flint 0.9.0, numpy 2.2.6, sympy 1.14.0, mpmath, and C.

**Overall verdict.**
- Every proved part of THM-4808 checks out mathematically.
- Two framing claims are false and need rewording: the "4-regular tree" over Q, and the abc non-connection.
- (d) and (e) are already proved in the repo and must be cited, not claimed.
- HYP-9418's numbers reproduce, except (5)'s "real root 2", which is wrong.
- The pre-registration timing is honest.

## 1. Verdict per claim

| Claim | Verdict | Notes |
|---|---|---|
| (G) M(t) = D(V·t), each element counted 4/\|V·t\| times; τ compatibility | CONFIRMED | Exact for every dart, all primes q ≤ 401 and 12 odd composite moduli |
| (G) B = P_U + P_A + P_D is the Hashimoto operator | CONFIRMED | Holds when each dart f carries index 4/\|tail f\|; the THM should state this |
| (G) Ihara–Bass λ² − μλ + 3 = 0 | CONFIRMED; upgrade to PROVED | Proof in §3.1; exact characteristic-polynomial identity for all primes q ≤ 401 |
| (G) "over Q this graph is a 4-regular tree …" | CORRECTED (major) | Two components, swapped by C(t) = (t+1)/(t−1); the Euclid one is an edge-indexed tree with a half-edge |
| (G) exact spectra for q ≤ 211 | CONFIRMED | Exact Sturm counts; true for all primes q < 331 |
| (G) Ramanujan below 1000 except 331, 569, 719, 809, 887 | CONFIRMED, now exact | 166 primes in range, not 165 |
| (G) at 331 a degree-82 factor has roots off the circle | CONFIRMED | Real roots −1.878429 and −1.597079; μ = −3.475508 |
| (G) ⟨U, A, D⟩ has index 3 in PGL₂(Z) | CONFIRMED | ⟨U,A,D⟩ = ⟨D,S,ι⟩; also AU⁻¹ = [−1 4; 0 1] = D²ι |
| (G) uniform expansion (cited) | CONFIRMED; fix attribution | The transfer is Brooks–Burger; comparison step in §3.6 |
| (a) χ = αn − γm | CONFIRMED | Classical determinant identity for consecutive convergents of the even (theta) continued fraction |
| (a) box 0 ≤ α ≤ m, 0 ≤ γ ≤ n | CONFIRMED; upgrade to PROVED | Cone argument (§3.2); in fact 1 ≤ α ≤ m−1 and γ ≤ α |
| (b) the continued-fraction automaton counts #A | CONFIRMED | Proof correct and complete |
| (b) "Not a character" paragraph | CORRECTED (minor) | χ = (−1)^(m + c₀+…+c_L) is an XOR of two simple features; the 18-feature list is misdescribed |
| (c) (16M+2, 4M+1): word A D^(M−1) U D, χ = −1, triple ≡ (3,4,5) | CONFIRMED | Proved for all M ≥ 1 (§3.4); computed for M ≤ 5000, including M = 1, 2 |
| (c) corollary on Dirichlet/Hecke characters | CONFIRMED | Infinite-order Hecke characters need one extra step (§3.5) |
| (d) Q-recursions, sign law, Pell sign on the A ray | CONFIRMED, NOT NEW | THM-3357 (19)–(20) and THM-3341 (24); must be cited |
| (e) level identity | CONFIRMED, NOT NEW | Scalar form of THM-3357 (37)/(40), via U − A + D = I |
| (f) pairing, interleaving, node partition, S = 1 + N_U + N_D + E | CONFIRMED | Independent code; the identity needs X ≥ 5 |
| (f) "parity bias … on these two rays and on no others" | CORRECTED (moderate) | Unproved; this is the content of HYP-9418 (2) |
| (f) numbers | CONFIRMED | S(X) to 10¹⁰ matches; sampled \|E\|/√X ≤ 0.21 from 10³ holds only at the samples (true maximum 0.368) |
| 7/21 non-connection | CONFIRMED | Cite THM-1370 by file name (two files share the ID) |
| abc non-connection | CORRECTED (moderate) | ADU gives an abc hit with χ = −1 |
| F1 non-connection | CONFIRMED | The section reference should be growing_sieve field §§4–5 |
| HYP (2) numbers and sign-change count | CONFIRMED | Add full-resolution data; suggest the Rubinstein–Sarnak form |
| HYP (3) all numbers | CONFIRMED | π₂ reproduced by an independent sieve; model values reproduced |
| HYP (4) 1000–4000 numbers | CONFIRMED | 42/382; worst 1.01047 at 3539 (negative edge); largest \|μ\| 3.5004 |
| HYP (4) "165 primes" and "as the Tracy–Widom picture predicts" | CORRECTED (minor) | 166; the random-graph model predicts about 31% exceptions, not a rise |
| HYP (5) exact for q ≤ 211, multiplicities "about q/4 + 3, q/4 + 1" | CONFIRMED | Exact formulas in §3.1 |
| HYP (5) "real root 2" at q = 251 | REFUTED | Real roots −1.999903 and −1.500073, irrational; more exceptions at 479, 863, 947 |
| Pre-registration timing | CONFIRMED honest | §6 |

## 2. Errors by severity, with repairs

### Major

**M1. Over Q the graph is not one tree.**
- **Two components.** The graph is defined on all of P¹(Q). U, A, D, τ, ι and S all preserve the parity type of
  m/n (m + n odd, or m and n both odd), so there are two components (`q_structure.out`, checks Q1–Q2).
- **The odd/odd component.** Its orbifold vertex is {±1}, the degenerate triangle (0,1,1). Its out-multiset is
  [1,1,3,3], t = 1 is a self-reverse dart, and its other vertices are the leg-swapped triangles.
- **The swap.** C(t) = (t+1)/(t−1) exchanges the two components. It satisfies CDC = U, CUC = D, CAC = A, CτC = τ,
  CιC = Sι and CSC = S; this is THM-3357's intertwiner H. Mod q, C is an automorphism commuting with B and B_χ
  (check T7).
- **The Euclid component** is an edge-indexed tree whose universal cover is the 4-regular tree. At (1,0,1) = {0,∞},
  the dart 2 to (3,4,5) and the self-reverse dart ∞ (a half-edge, the THM's "doubled loop") both have index 2. All
  other indices are 1.
- **Likely cause.** The THM's phrase "up to leg swap" suggests the two components were merged by mistake.
- **Repair:** edits 1, 8, 9. The results note §2 needs the same fix.

**M2. Prior art is not cited.**
- **(d)** is THM-3357 (19)–(20) (2026-08-14, independently audited). There L, M, R are U, A, D, and
  σ = sign(b − a) satisfies σ(Lu) = +1, σ(Mu) = −σ(u), σ(Ru) = −1. The middle-ray fact |a − b| = 1 is THM-3341 (24).
- **(e)** is the scalar form of THM-3357 (37)/(40). U − A + D = I, so the signed sum Σ χ(w) M_w over a level is the
  identity matrix.
- **χ** is THM-3357's (−1)^#M Walsh character.
- **Naming hazard.** The "A-Walsh character" of THM-4057/4059 is (−1)^#U in THM-4808's letters, because their A is
  U here.
- **The non-backtracking reading over Q** is Alperin, *The modular tree of Pythagoras*, Amer. Math. Monthly 112
  (2005) 807–816 (arXiv math/0010281).
  - Alperin enumerates the tree by reduced words in the free group Γ(2) = ⟨U², L²⟩, where U² = D and
    L² = S D⁻¹ S⁻¹.
  - He quotients by the cusp stabiliser ⟨U²⟩ and an involution δ (ι here), and grows the tree by adding one of three
    letters that keep a word reduced.
  - The THM's Euclid component is exactly this Schreier graph, folded by δ. The neighbours V(t ± 2) and
    V(t/(1 ± 2t)) were checked.
- **(a), (b).** The Berggren word of m/n is its even-integer continued fraction, t = 2a₁ + ε₁/(2a₂ + …). A D-run
  followed by U or A gives an even digit with sign −1 or +1, so #A is the number of + signs.
  - (a) is therefore the classical determinant identity for consecutive convergents.
  - References: Kraaikamp–Lopes, Geom. Dedicata 59 (1996) 293–333; Romik, Trans. AMS 360 (2008) 6045–6064.
- **Repair:** edits 3, 7, 9, 17, 18.

### Moderate

**D1. The (f) bias claim is unproved.** Only the main term (1/2 + 1/√2)√X comes from the two root rays. E(X) is not
shown to be small, and its full-resolution maximum on [10², 10⁹] is 0.368√X. Repair: edit 22.

**D2. The abc non-connection is false as a general statement.**
- Use the squared packet (a², b², c²) of growing_sieve field §4.
- The word ADU maps (2,1) → (5,2) → (9,2) → (16,9). Its triple is (175, 288, 337), with
  c² = 113569 > rad = 2·3·5·7·337 = 70770, and χ = −1.
- The field packet itself shows that the words empty, UU, A and ADU realise all four (χ, hit) combinations.
- Repair: edit 24. The results note's abc row has the same problem.

**D3. HYP (5) "real root 2" is wrong.**
- At q = 251 the degree-62 irreducible factor of charpoly(B_χ) has exactly two real roots, −1.9999032793 and
  −1.5000725440. Their product is 3, and they are irrational.
- The twisted-adjacency eigenvalue is μ = −3.4999758233, from a degree-31 factor.
- The lane's value 0.6666 was |λ|/3 with the sign dropped.
- The exact census below 1000 adds exceptions at 479, 863 and 947.
- Repair: edits 26, 33.

### Minor

1. **Ihara–Bass.** Labelled "checked numerically", but provable (§3.1). The lane's (K2) check is tautological in code,
   because both sides are built from the same list, and (K3) uses numerical tolerances. Repair: edits 4, 10.
2. **"Exact" spectra.**
   - The lane's `berggren_ihara.gp` checks moduli with 38-digit `polroots`, which is not exact.
   - `berggren_exact_331.out` is truncated at 600 characters per line, so it never shows the root.
   - The audit's census is exact: Sturm counts of R(y) on [12, ∞), where R(x²) = ±F₁(x)F₁(−x) and
     F₁ = charpoly(A)/(x − 4).
   - Repair: edits 4, 11.
3. **(a) box.** Labelled FINITE-EXACT, but provable (§3.2). Repair: edits 5, 13.
4. **(b) closed form.** χ = (−1)^(m + Σcᵢ) (§3.3) makes the automaton unnecessary for χ.
   - The "Not a character" paragraph misdescribes its features. It lists Jacobi symbols of c and classes of a, b, c
     mod 16; the script used others, including the last-letter tests m > 2n and m > 3n.
   - It omitted the Stern–Brocot depth parity, with which χ is an XOR of two features.
   - Repair: edits 2, 14, 15.
5. **(c) Hecke wording.** The inference needs finite order, or the extra step in §3.5. Repair: edit 16.
6. **Expansion attribution.**
   - "Lubotzky's transfer" should be the Brooks–Burger transfer, as presented in Lubotzky 1994, Ch. 4.
   - The step from theta-group Schreier graphs to G_q should be stated (§3.6).
   - The sieve paragraph's equidistribution claim also needs the bottom of the spectrum bounded away from −4. The same
     citation gives this, applied to the finite-index subgroup of two-step words.
   - Repair: edits 6, 12, 34.
7. **(f) details.**
   - The identity needs X ≥ 5: below 5, S = 0 while the right side is 1.
   - The sentence about rays "based at nodes ending in U, D or A" should say exactly which node starts which ray.
   - Repair: edits 19, 20.
8. **HYP (2) bounds.**
   - "|E|/√X ≤ 0.30, and ≤ 0.21 from 10³" holds only at the 22 sampled heights. Over all integers in [10³, 10⁹] the
     maximum is 0.368, so the pre-registered ±0.21 band came from samples.
   - By analogy with Chebyshev's bias (Rubinstein–Sarnak), the natural conjecture is a limiting logarithmic
     distribution for E(X)/√X. A uniform O(√X) bound is stronger than that analogy supports.
   - Repair: edits 21, 28, 29.
9. **HYP (4) counts and model.**
   - There are 166 primes in 5 ≤ q < 1000, not 165. PREREGISTRATION entry (5) has the same slip.
   - "As the random-graph Tracy–Widom picture predicts" is unsupported. Huang–McKenzie–Yau (arXiv:2412.20263) give
     about 69% Ramanujan for random d-regular graphs, hence about 31% exceptions.
   - The observed 3% and 11% are far below that, and the model predicts no rise.
   - "The exceptional primes need not stay exceptional" is unclear.
   - Repair: edits 31, 32.
10. **References.** Cite THM-1370 by file name (`THM-1370-elliptic-jacobian-conjecture-all-dimensions.md` also
    exists). The radical non-descent is growing_sieve field §4 (the witness) together with §5 (the monoid argument).
    Repair: edits 23, 25.
11. **PREREGISTRATION.md (not edited).** It says "Filed under HYP-9418 (3)-(4)", but its entries (3), (4) and (5) are
    HYP (3), (2) and (4). It also says "5 of 165 primes" instead of 166.
12. **HYP (3).** The predicted quantities are fixed by published data (OEIS A007508) and an elementary product. The
    test checks arithmetic, not the model. Repair: edit 30.
13. **Results note (outside the edit scope).** It needs the same fixes as M1, D1 and D2: its §2 "4-regular tree", its
    "bias sits on the two root rays", and its abc row. Its new log-spaced statistics (maximum 0.346, log-mean −0.0043)
    agree with the full-resolution data here.

HYP-9418's scripts list also gets a pointer to the audit directory (edit 27).

## 3. Proofs supplied by the audit

### 3.1 Ihara–Bass with multiplicities

**Setup.** Let head(e) = V·e and tail(e) = V·D⁻¹e. Define (Lg)(e) = g(head e), (Tg)(e) = g(tail e), and
(Ag)(O) = Σ_{v∈V} g(V·D v s_O), which does not depend on the choice of s_O.

**Two facts.**
- As lists, [Ue, Ae, De, τe] = [Dve] over v ∈ V.
- tail(Xe) = head(e) for X ∈ {U, A, D}, and V·τe = tail(e).

**Identities.** These facts give BT = 3L and BL = LA − T.

**Eigenvectors with λ² ≠ 1.**
- Write B = K − J, where J is reversal and K[e,f] = (4/|head e|)·[head e = tail f].
- If Bφ = λφ, write the eigen-equation at e and at τe and eliminate. This shows φ = Lg₁ + Tg₂.
- Comparing terms gives g₁ = −λg₂ and Ag₂ = (λ + 3/λ)g₂, which is λ² − μλ + 3 = 0.

**Reality of μ.** A is self-adjoint for weights |O|: |O|·A[O,O′] = 4·#{s ∈ O : Ds ∈ O′}, and s ↦ ιDs is a bijection
onto the reverse set.

**The eigenvalues ±1.**
- On the joint kernel of the two dart-sum maps, B = −J, which contributes only ±1.
- The vertex part has dimension 2n_V − 1, because L1 = T1.
- Counting the eigenspaces of J gives the multiplicities below.

| Operator | q ≡ 3 (mod 4): mult(1), mult(−1) | q ≡ 1 (mod 4): mult(1), mult(−1) |
|---|---|---|
| B | (q−3)/4, (q+1)/4 | (q−5)/4, (q−1)/4 |
| B_χ | (q+9)/4, (q+1)/4 | (q+7)/4, (q−1)/4 |

These match every computed prime ≤ 401.

**Twisted case.** For χ-equivariant g, B_χ = K_χ + J, with B_χ g = A_χ g + g∘τ and B_χ(g∘τ) = −3g. This gives the same
quadratic. The orbits {0,∞} and {±1} drop out because χ(ι) = χ(Sι) = −1. The orbit {±i} survives because χ(S) = 1.

**Exact identities.** Both of the following hold for all 78 primes 3 ≤ q ≤ 401 and for N = 9, 15, 21, 25, 27, 33, 35,
45, 49, 63, 75, 105:
- charpoly(B)(x)·(x−1) = (x−1)^α (x+1)^β · x^{n_V} F_A(x + 3/x);
- charpoly(B_χ)(x) = (x−1)^α′ (x+1)^β′ · x^{n′} F_{A_χ}(x + 3/x).

### 3.2 Cone argument for (a)

- U, A and D map the cone K = {x ≥ y ≥ 0} into itself, and (1,0), (1,1) ∈ K.
- (m,n) = M_w(1,0) + M_w(1,1), so (α,γ) = M_w(1,0) ∈ K and (m−α, n−γ) = M_w(1,1) ∈ K.
- Hence 0 ≤ γ ≤ α and 0 ≤ n−γ ≤ m−α.
- α = 0 or α = m would force αn − γm = 0. So 1 ≤ α ≤ m−1 and 0 ≤ γ ≤ n.

### 3.3 Closed form for χ

**Claim.** χ = (−1)^(m + c₀ + … + c_L), where m/n = [c₀; c₁, …, c_L].

**Proof.**
1. U and A swap the parities of (m, n), and D keeps them. The root has m = 2, so m ≡ #U + #A (mod 2).
2. The digit sum c₀ + … + c_L goes up by 2 under D and A and by 1 under U, including t = 2 and integer t. This is
   THM-4057 Theorem 4.1 (whose A is U here). So Σcᵢ = 2 + #U + 2#A + 2#D ≡ #U.
3. Subtracting gives #A ≡ m + Σcᵢ (mod 2). ∎

This was checked on all 456,137 pairs with m ≤ 1500, together with both lemmas. Equivalently,
χ = (−1)^(L + [w starts with U]).

### 3.4 (c) for every M

- **Descent.** (16M+2, 4M+1) goes to (8M, 4M+1) by D⁻¹, then to (4M+1, 2) by U⁻¹, then to (5, 2) by M−1 steps of D⁻¹
  (each at slope ≥ 9/2), then to (2,1) by A⁻¹.
- **Validity.** gcd = gcd(2, 4M+1) = 1, and the parities are opposite.
- **Congruence.** D^M ≡ I (mod M), so M_w ≡ DUD⁻¹A = R = [−1 4; 0 1] (mod M). R is the reflection t ↦ 4 − t, and
  AU⁻¹ = R exactly.
- **Cross-check.** The continued fraction is [3; 1, 2M−1, 2]. Its digit sum is odd and m is even, so by §3.3, χ = −1.

### 3.5 Hecke characters

- Suppose a Hecke character ψ of Q(i) agrees with the ±1-valued χ on Euclid α = m + ni.
- On α ≡ 1 modulo the conductor, ψ(α) = (α/|α|)^k.
- The arguments of such α are dense, so k = 0. Then ψ has finite order, and (c) applies.

### 3.6 Expansion transfer to G_q

- For V-invariant f, E_A(f) = ⟨f, (4−A)f⟩ = ½ Σ_x Σ_v (f(x) − f(Dvx))² ≥ ½ Σ_x (f(x) − f(Dx))².
- Since f(Sx) = f(x), the last sum is the Dirichlet form of the Schreier graph of ⟨D, S⟩ with generators {D^±1, S}.
- Its spectral gap is uniform in q: Γ_θ ∩ Γ(q) contains Γ(2q), so Selberg's 3/16 theorem applies, transferred by
  Brooks–Burger.
- So G_q inherits half of that gap.

## 4. Numbers reproduced with independent code

### Spectra

- **Exact census of the adjacency, 5 ≤ q < 1000 (166 primes).** |μ| > 2√3 occurs exactly at:
  - 331 (μ = −3.475508);
  - 569 (−3.465441);
  - 719 (+3.510074);
  - 809 (−3.475155);
  - 887 (+3.468820).

  No eigenvalue equals ±2√3. The eigenvalue 4 is always simple. The worst ratio is 1.01327, at q = 719.
- **q = 331.** The degree-82 irreducible factor of charpoly(B) has exactly two real roots, −1.8784287477 and
  −1.5970794759.
- **1000 ≤ q < 4000 (382 primes, numerical).**
  - 42 exceptions, the same list as the lane's (0.1099).
  - Worst ratio 1.01047 at q = 3539 (μ = −3.500364); largest |μ| 3.5004.
  - Exact spot checks: 1069 (excess 2.4·10⁻⁵) and 2909 are exceptions; 1061, 1063 and 2917 are Ramanujan.
- **B_χ.**
  - Exact exceptions below 1000: 251, 479, 863, 947.
  - 1000–4000: 42 numerical exceptions, worst ratio 1.00979 at q = 1567.
- **Multiplicities.** The ±1 multiplicities of B and B_χ match `berggren_ihara.out` for all 46 primes ≤ 211.

### Tree

- **(a), (b), (d).** All hold on all 456,137 pairs with m ≤ 1500: (a), the box and cone argument, the automaton of
  (b), both closed forms, and the sign law of (d).
- **(d) symbolically.** The (d) recursions hold symbolically. The A ray has |a−b| = 1 and χ = −(a−b) for k ≤ 39.
- **(c).** Verified for M = 1 to 5000.
- **(e).** Holds symbolically for k ≤ 7.
- **(f) identity.**

  | X | 10³ | 10⁴ | 10⁵ | 10⁶ | 10⁷ |
  |---|---|---|---|---|---|
  | S | 30 | 111 | 409 | 995 | 3817 |
  | E (from local terms) | −5 | −7 | 29 | −209 | 2 |

  The identity holds at every one of these heights. Every δ is in {0, 1}, and there are no interleaving violations.
- **(f) node partition at 10⁷.**
  - 1 root;
  - 287,548 A-children;
  - 652,012 D-ray interior nodes;
  - 652,018 U-ray interior nodes.
- **Root rays.** |N_U + N_D − (1/2 + 1/√2)√X| ≤ 4.41 for X ≤ 10²⁴.

### Height sums

The audit traversed the tree depth-first, visiting 1.59·10⁹ nodes.

| X | 10⁸ | 2·10⁸ | 5·10⁸ | 10⁹ | 2·10⁹ | 5·10⁹ | 10¹⁰ |
|---|---|---|---|---|---|---|---|
| S(X) | 12158 | 18270 | 26750 | 36468 | 55546 | 86222 | 129405 |
| E(X)/√X | 0.00900 | 0.08492 | −0.01069 | −0.05382 | 0.03499 | 0.01229 | 0.08697 |

These, the triple counts N(X), and all 22 sampled E values from 100 to 10⁹ match the lane exactly.

**Full resolution, every integer X ≤ 10⁹.**
- **Maximum.** The largest |E|/√X is 0.36789, at X = 4,222,837, where E = −756. A second traversal confirmed it.
- **Per-decade maxima.** 0.319, 0.290, 0.339, 0.329, 0.368, 0.284, 0.301.
- **Sign changes.** 13,538 on [100, 10⁹].
- **Logarithmic mean of E/√X.** −0.0248 over [5, 10⁹], and −0.0077 over [10⁴, 10⁹].

### Twin sieve

- π₂(10⁸) = 440312, π₂(10⁹) = 3424506, π₂(10¹⁰) = 27412679.
- ρ = 1.11158, 1.13123, 1.14480.
- The HL–Mertens model gives 1.11417, 1.13178, 1.14555.
- 4e^(−2γ) = 1.26095.

## 5. Literature and attribution

| Item | Finding |
|---|---|
| Berggren 1934, Barning 1963, Hall 1970 | Correct (Tidskr. Elem. Mat. 17; Math. Centrum ZW-011; Math. Gazette 54, 377–379) |
| Non-backtracking reading over Q | In substance Alperin 2005; the THM's "no source found" should cite it |
| Theta-group / Γ(2) form; even continued fraction | Classical: Alperin 2005; Romik 2008; Kraaikamp–Lopes 1996; Short–Walker (arXiv 1508.01373) |
| (a) pole frame | Classical convergent determinant identity |
| (b) automaton | Regular-to-even continued fraction conversion; the parity closed form follows from THM-4057 |
| Mod-q Hashimoto / Ihara–Bass for Berggren; Ramanujan census | Not found; no priority claimed |
| Branch-character height sum S(X) | Not found; THM-4059 studies a different character at a different height |
| Selberg; Lubotzky 1994 | Correct (Progress in Math. 125, Ch. 4); the transfer is Brooks–Burger |
| Huang–McKenzie–Yau | Correct (arXiv:2412.20263); its consequence is misread in HYP (4) |
| Concurrent packets | Cited accurately in content; the field packet's own §4 data contradicts the THM's abc sentence |

## 6. Pre-registration timing

**Entries (3) and (4).**
- Committed in `2fe36c8088` at 16:07:04.
- The binaries that produce the predicted values are in the shared session scratchpad. They were compiled afterwards:
  `twin_sieve_ratio` at 16:07:50 and `parity_height_big` at 16:08:25.
- The earlier binaries (15:54 and 16:04) stop at 10⁹ and have no twin data.
- The .out files carried a uniform copy time (16:21:58) until commit `4950589d05` (16:30:52). The compile times are
  the binding evidence.

**Entry (5).**
- The scan script with the lower-bound argument was saved at 16:24:15.
- The pre-registration was committed in `a0a2bbd805` at 16:24:16.
- The output appeared at 16:25:05.
- Re-run unmodified, the lane's scan takes 43.4 s here and reproduces the stored output byte for byte.
- So the predictions were committed before the output existed.

**Verdict:** honest. Caveats: entry (3) predicts quantities already fixed by published data, and entry (4)'s ±0.21 band
was drawn from sampled heights only.

## 7. Deliverables

All files are in `04-computation/experiments/berggren_sieve_parity_20261026/audit_thm4808/`. Binaries were built in the
session scratchpad, not in the repository.

| Script | Output | Content | Runtime |
|---|---|---|---|
| `g_structure.py` | `g_structure.out` | (G) checks T1–T7; exact Ihara–Bass identities | ~19 s |
| `spectra_exact.gp`, `spectra_exact_run.gp` | `spectra_exact.out` | Exact Sturm census below 1000; the 331 and 251 factors | ~16 s |
| `spectra_exact_spot.gp` | `spectra_exact_spot.out` | Exact spot checks above 1000 | ~1 min |
| `ramanujan_scan.py` | `ramanujan_scan_5_1000.out`, `ramanujan_scan_1000_4000.out` | Numerical scans of A and A_χ | ~22 s |
| `tree_checks.py` | `tree_checks.out` | (a)–(f) and the closed forms | ~10 s |
| `cf_closed_form_search.py` | `cf_closed_form_search.out` | The search that found the closed form | <1 s |
| `q_structure.py` | `q_structure.out` | Components over Q; the involution C | ~1 s |
| `height_sum.c` | `height_sum_mode1_1e10.out`, `height_sum_mode2_1e9.out` | S(X) to 10¹⁰; full resolution to 10⁹ (250 MB window) | ~11 s, ~12 s |
| `twin_rho.c` | `twin_rho.out` | Independent twin-prime sieve to 10¹⁰ | ~10 s |
| `hl_model_check.py` | `hl_model_check.out` | HL–Mertens model values | <1 s |
| — | `prereg_timing.out` | Timing evidence for §6 | — |
| `make_required_edits.py` | `required_edits.json` | 34 validated edits (1–25 THM-4808, 26–34 HYP-9418) | <1 s |

**Sources.**
- Alperin, *The Modular Tree of Pythagoras* (arXiv math/0010281).
- Romik, *The dynamics of Pythagorean triples* (arXiv math/0406512).
- Kim–Lee–Liao, *Odd-odd continued fraction algorithm* (arXiv 2002.06310).
- Short–Walker, *Even-integer continued fractions and the Farey tree* (arXiv 1508.01373).
- Kraaikamp–Lopes, cited via arXiv 2403.12341.
- Huang–McKenzie–Yau (arXiv 2412.20263).
- Kontorovich, *From Apollonius to Zaremba* (arXiv 1208.5460).
- Kontorovich–Oh, *Almost prime Pythagorean triples in thin orbits* (arXiv 1001.0370).
- Hall, *Genealogy of Pythagorean triads* (Math. Gazette 54).
- Lubotzky 1994 (Springer).
- Cha–Kim, *Number theoretical properties of Romik's dynamical system* (arXiv 1902.03571).
