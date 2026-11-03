# Reversal-symmetric hidden storage and exact palindromic clocks

Session 174. **PROVED constructions; finite independent controls.** Reconstructed in the shared workspace after the temporary workspace reset. No prime-producing family or solution of the missing Hamilton-count candidates 7 and 21 is claimed.

## 1. Concrete reversal and resilient storage

For a tournament X let X* reverse all arcs, p_X(k) count unordered covers by k nonempty directed paths, and F_X(u)=Σ p_X(k)(u)_k. Reversing paths proves F_X*=F_X. Write ⊕ for ordered sum. Then

    (X_1⊕...⊕X_M)* ≅ X_M*⊕...⊕X_1*,
    F_(X⊕Y)=F_X F_Y.

The second identity follows by placing an X segment followed by a Y segment in each labelled path slot, either segment possibly empty.

Use five-vertex seeds with outgoing-neighbor masks

    B_0=(16,1,3,7,14), B_1=(20,1,2,7,14),
    f(u)=F_B0(u)=F_B1(u)=u^5+6u^3+2u.

The checker verifies strongness, nonisomorphism over all 120 permutations, and the full-profile equality. The seeds differ on one labelled arc. For b∈{0,1}^M put A_b=B_b1⊕...⊕B_bM. Construct T_b(r) from A_b, its disjoint dual copy A_b*, and a transitive r-vertex group Z_r with cyclic uniform orientation

    A_b → A_b* → Z_r → A_b.

**Theorem 1.** If M,r≥1 and 5M>r, these 10M+r vertex tournaments are self-dual, have vertex connectivity exactly r and unique separating r-set Z_r, and form 2^M pairwise nonisomorphic tournaments with the same complete profile. Labelled arc distance is twice Hamming distance.

**Proof.** Swap corresponding vertices of the two old sides and reverse Z_r. This explicit involution σ reverses all arcs. It is an anti-automorphism, not an automorphism. After deletion of fewer than r vertices, all three groups remain nonempty and cyclically oriented, hence strong. After exactly r deletions this remains true except when deleting precisely Z_r. Deleting Z_r leaves A_b⊕A_b*, which is not strong. Thus Z_r is the unique minimum cut.

Every isomorphism preserves this cut. Its removal reveals the uniquely ordered strong-component sequence

    B_b1,...,B_bM,B_bM*,...,B_b1*.

The first M components recover b. For equal complete profiles, collapse maximal runs inside uniform modules in every labelled cover slot. A skeleton containing k runs in X has k!p_X(k) fillings. Summing over legal skeletons depends only on the base orientation and individual profiles. Here both old sides have f^M. Finally each changed seed arc appears once in each half. ∎

This is an actual word-reversal/tournament-duality bridge. Vertex resilience concerns connectivity, not survival of stored bits after vertex erasure.

## 2. General fixed-feedback Hamilton recurrence

Write q_A(k)=k!p_A(k). In a Hamilton path of the cyclic wrapper, collapse A-runs and L-runs, and keep the r individually labelled feedback vertices. Each old side occurs at most r+1 times, because a return to that side requires a feedback vertex. Let c_r(i,j) count the legal skeletons. Then for arbitrary nonempty A,

    H(T_A(r))=Σ_(i,j=1)^(r+1) c_r(i,j)q_A(i)q_A(j).       (1)

The coefficients are nonnegative integers. Since F_A(u)=Σ q_A(k)binom(u,k), inversion gives

    q_A(i)=Σ_(s=1)^i (-1)^(i-s)binom(i,s)F_A(s).

Thus an effectively computable symmetric integer matrix K_r satisfies

    H(T_b(r))=Σ_(s,t=1)^(r+1) K_r(s,t)[f(s)f(t)]^M.       (2)

Reversing skeletons, swapping A and L, and reversing feedback indices proves symmetry. The positive formula (1) supplies the combinatorial meaning of the signed evaluation formula (2).

**Corollary.** For fixed r, H_M has a rational generating function and an integer linear recurrence of order at most (r+1)(r+2)/2. An annihilating polynomial is

    ∏_(1≤s≤t≤r+1)(E−f(s)f(t)).

Repeated or zero-coefficient bases may be omitted. This follows termwise for geometric sequences. Their generating functions for M≥1 are λt/(1−λt). Modulo any m the sequence is eventually periodic; if gcd(m,∏_(s=1)^(r+1)f(s))=1, all bases are units and the sequence is purely periodic. This concerns the displayed family, not an automaton deciding unrestricted tournament spectra.

## 3. Exact clocks and a misleading last digit

Finite skeleton enumeration gives

    K_1=[[3,−1],[−1,1]],
    K_2=[[15,−11,2],[−11,13,−4],[2,−4,2]].

Using f(1)=9, f(2)=84, f(3)=411,

    H_M(1)=3·81^M−2·756^M+7056^M,
    H_M(2)=15·81^M−22·756^M+4·3699^M
                       +13·7056^M−8·34524^M+2·168921^M.

At M=1 these are 5787 and 152757. For all M≥1,

    v_3(H_M(1))=v_3(H_M(2))=2M,
    H_M(1)≡7 (mod10), H_M(1)≡3 (mod8),
    H_M(2)≡7 (mod10) if M is odd, and ≡9 if M is even.

For r=1 the three summands have respective 3-valuations 4M+1,3M,2M. For r=2 divide by 9^M; the first three summands vanish modulo 3 and the remainder is

    13·784^M−8·3836^M+2·18769^M ≡ −8·2^M ≠0 (mod3).

The final-digit statements follow by reducing each base modulo10, and the mod8 statement follows similarly. Thus counts may end in 7 while all being divisible by 9. They cannot realize 7 or 21 and are not prime.

## 4. Reversal-fixed Hamilton paths have smaller clocks

The reversal-conjugation involution and odd fixed-count principle are inherited from `01-canon/theorems/THM-582-palindromic-hampath-odd-index.md`. The unique-center, odd-order half-data construction is inherited from `01-canon/theorems/THM-583-palindromic-count-is-half-system-hampath-count.md`. Here they specialize to σ from Theorem 1; the even-order seam condition below is a separate elementary extension. Its action on a Hamilton path reverses vertex order and applies σ. A fixed path is determined by its first half, choosing exactly one vertex from each two-element σ-orbit. If r is odd, its center is the unique fixed vertex, the middle of Z_r; if r is even there is no center.

The first half uses every original A vertex exactly once, either as its A copy or as its A* copy. It chooses one vertex from each mirrored feedback pair. Its final vertex v must beat σ(v) in the even case, or the central feedback vertex in the odd case. These conditions are sufficient: the second half is forced and its arcs follow by anti-symmetry.

Collapse old-side runs in the half path. A skeleton with k runs has q_A(k) fillings: use an ordered path cover of A and reverse the internal order of runs assigned to A*. Let d_r(k) count legal half skeletons including the feedback choices and final-vertex condition. Then

    Fix_σ(T_A(r))=Σ_(k=1)^(r+1) d_r(k)q_A(k).             (3)

For the bound, q=floor(r/2) feedback representatives allow at most two old-side runs in each of q+1 intervals. Odd r gives 2q+2=r+1. For even r, a half path must end in A or a lower-index feedback representative; these alternatives reduce the maximum to 2q+1 or 2q, respectively. Thus k≤r+1.

By binomial inversion, for A_b the fixed count has at most r+1 geometric bases f(j)^M, versus at most (r+1)(r+2)/2 product bases for all paths. This proves that the chosen reversal-fixed observable is also blind to the stored word.

The first half-skeleton vectors are d_1=(1,1) and d_2=(3,3,2). Therefore

    Fix_σ(T_b(1))=84^M−9^M,
    Fix_σ(T_b(2))=2·411^M−3·84^M+3·9^M.

For M=1 these are 75 and 597. The r=3,4 counts are 1059 and 7455 for either seed. The involution partitions Hamilton paths into fixed points and pairs, so fixed and total counts have equal parity. The general oddness of tournament Hamilton-path counts is inherited, not proved anew here. No equality of all possible equivariant invariants or other choices of anti-involution is claimed.

## 5. Decimal palindromes as finite-state encoders

Starting with a central digit c, add digit d to both ends of a length 2k+1 palindrome of value P_k. Then

    P_(k+1)=10P_k+d(10^(2k+2)+1).

Modulo any m, initialize (p,a)=(c,100) and update

    (p,a)↦(10p+d(a+1),100a).                             (4)

This two-register transducer computes residues exactly, including moduli sharing factors with 10. For constant d and gcd(m,10)=1 it is a permutation on m² states: first invert100 to recover a, then invert10 to recover p. Its trajectories are purely periodic. Nonzero outer digits avoid leading-zero ambiguity.

Taking m to be a product of any finite collection of primes and marking states nonzero at each prime gives an exact finite sieve. Passing that sieve does not certify primality.

A concrete injective bridge uses d_i=1+2b_i and the palindrome

    d_1...d_M 7 d_M...d_1.

It and T_b(r) encode the same bit word. The decimal digit sum is 7+2M+4wt(b), so avoiding 3 is exactly 1+2M+wt(b)≠0 modulo3. Any finite sieve therefore selects a well-defined subfamily of distinct hidden tournaments without changing their common complete profile or chosen reversal-fixed count. No prime infinitude or density follows.

Taking the direct product of (4) and the registers λ^M modm from (2) gives an explicit synchronized finite-state observer of both arithmetic values. Under the invertibility conditions and constant d it is a permutation clock. It preserves separate values, rather than equating a decimal palindrome with a Hamilton count.

## 6. Large expression and elementary controls

For every k≥1,

    10^(2k+1)−10^k−1 = decimal [k nines] 8 [k nines].

Subtracting one first gives 2k+1 nines, then subtracting 10^k changes the middle nine to eight. Thus 10^1888529−10^944264−1 is an odd palindrome of 1,888,529 digits. This proves its shape, not its primality. Every member is 2 modulo 3 and −2−(−1)^k modulo 11, so neither 3 nor 11 divides it.

Every even-length decimal palindrome is divisible by 11: symmetric pairs cancel in the alternating digit sum. In particular 1001=7·11·13 and 2222=2·11·101. These do not imply that odd palindromes are composite. The decimal reversal of 7417 is 7147, and 240 is not a palindrome. Those numbers may have other intended roles. No primality assertion about Belphegor's number or the enormous expression is needed, and no Banks–Hart–Sakata theorem is quoted or strengthened here.

## 7. Reproducibility, attribution, and open interfaces

Run `python3 check.py > certificate.json` and `python3 -O check.py > certificate.optimized.json`, then compare. Do the same with `fixed.py` and `fixed_certificate*.json`. These dependency-free programs use no Python assert. `check.py` verifies seed claims, literal subset-DP Hamilton counts and explicit anti-involutions for 20 graphs through 14 vertices, full profiles of both encodings for r=1,2, every potential separating set through size r for r=1,...,4, derived kernels through r=4, arithmetic clocks through M=50, and 31,590 literal decimal/modular transition comparisons. `fixed.py` implements a separate pair-mask recurrence, compares to full permutation enumeration for orders≤8, and compares with independently enumerated half-skeleton coefficients through r=4. General claims are proved above rather than extrapolated from tests.

The uniform-module substitution lemma and seed data are inherited from Session 173's spectrum note and checker. The reversal-conjugation and half-data mechanisms are inherited as THM-582/583, not novelty claims. New constructions here are the mirrored cyclic storage family, its bounded geometric recurrences, the fixed-path specialization, and the synchronized palindrome observer. No repository-wide novelty claim is made. No prohibited ledger was read and no software was installed.

Open targets: a concrete reversal-compatible orientation rule whose H tracks prescribed prime palindromes; a constructor escaping the seeds' forced 3-divisibility; and an actual explanation or resolution of missing candidates 7 and 21. A solution would need a proved Hamilton-path bijection/count formula, not merely numerical phase matching.
