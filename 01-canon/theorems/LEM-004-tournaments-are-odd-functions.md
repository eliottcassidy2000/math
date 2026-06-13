# LEM-004: Tournaments Are Odd Functions (the odd-function dictionary)

**Type:** Lemma (a four-entry dictionary)
**Certainty:** 5 — PROVED (entries (a), (b), (d), and (c1)–(c3); elementary arguments, no Weil/Gauss input)
**Status:** PROVED + VERIFIED (exhaustive n = 3..11 for (a); all p < 60, exhaustive over all sign functions p ≤ 19 for (b); 13 families exactly for (c); exact Fraction series to deg 21 for (d)) — with one honestly-labeled OPEN boundary at (c), recorded below
**Added by:** kind-pasteur-2026-06-10-S2 (Thread A, HYP-2378)
**Disputes:** none (two SCOPE FLAGS on prior statements, see "Reconciliation")
**Tags:** #odd-functions #circulant #Paley #Legendre #cherry #cluster-expansion #formal-group #arctanh #complementation #OCF-context

---

## Slogan

**The tournament condition IS the oddness condition.** A tournament chooses exactly one of each arc pair — algebraically, a ±1-function that anticommutes with negation. Every appearance of "p ≡ 3 mod 4", of the cherry weight +1, and of arctanh in this repo is the same fact read in four languages.

## Statement

**(a) The dictionary (group level).** Let n ≥ 3 be odd. The map S ↦ σ_S, σ_S(d) = +1 if d ∈ S and −1 otherwise, is a **bijection** from connection sets of circulant tournaments on Z_n (arc i→j iff j−i ∈ S) onto the **odd ±1-functions** on Z_n \ {0} (σ(−d) = −σ(d)). Both sides have exactly **2^((n−1)/2)** elements. Under it:
- T^op ↔ −σ (complementation = global sign flip);
- relabeling x ↦ ux (u ∈ Z_n^×) acts by σ(d) ↦ σ(u⁻¹d);
- for **even n** both sides are **empty** (d = n/2 is its own negative: σ(n/2) = −σ(n/2) is impossible, and a connection set cannot contain "exactly one" of {n/2, n/2});
- boundary: a general **labeled tournament** on a vertex set V is exactly an odd ±1-function on the off-diagonal **pairs** V×V∖Δ (σ(j,i) = −σ(i,j), the skew sign matrix A−Aᵀ); "circulant" is precisely the case where σ descends to the difference group Z_n.

**(b) The multiplicative classification (MISTAKE-011b as a one-line theorem).** Let p be an odd prime. The completely multiplicative ±1-functions on F_p^× are exactly the two real characters: the **trivial character** (which is EVEN, hence never a tournament) and the **Legendre symbol χ**. χ is odd ⟺ χ(−1) = −1 ⟺ **p ≡ 3 (mod 4)**. Hence a multiplicatively-structured circulant tournament on F_p exists **iff p ≡ 3 mod 4, and it is unique: the Paley tournament** (S = QR set). For p ≡ 1 mod 4, χ is even — the QR set is closed under negation and gives bidirectional arcs, not a tournament.

**(c) Cherry localization (the oddness content of the e-mechanism, HYP-2307/THM-438).** Let σ be ANY odd ±1-function on Z_n \ {0} (n odd), M the circulant skew matrix M[a,b] = σ(b−a), and
A_L = Σ_{x_0,…,x_L distinct} ∏ σ(x_{i+1}−x_i) the single-run cluster integral. Then, with NO multiplicativity, NO quasirandomness:

- **(c1)** A_L = 0 for every **odd** L (negation symmetry — the analogue of χ(−1) = −1 needs no multiplicativity);
- **(c2)** A_2 = −tr(M²) = **n(n−1) exactly**: the per-pair cherry weight is −σ(d)σ(−d) = σ(d)² = **+1, the oddness identity itself**. Hence the cherry generator a_2 = 1 for EVERY circulant tournament; the cherry ingredient of e = exp(−χ(−1)) localizes to oddness alone.
- **(c3) (single-merge reduction)** A_L = −tr(M^L) + E_L(n) with |E_L(n)| ≤ C_L·n^(L−1), C_L depending only on L. Consequently a_L := lim A_L/n^L = −lim tr(M^L)/n^L whenever the spectral limit exists. In particular suppression of the higher even generators is a **spectral** (quasirandomness) property, NOT an oddness property — see the honest boundary below.

**(d) The tournament formal group is odd.** F(x,y) = (x+y)/(1+xy) ∈ Z[[x,y]] satisfies **[−1]_F(x) = −x exactly** (F(x,−x) = 0), equivalently its logarithm is **arctanh(x) = Σ x^(2k+1)/(2k+1)**, an odd power series. General lemma: for a (commutative, 1-dimensional) formal group law F over a ring R with no additive torsion, the following are equivalent:
  (i) [−1]_F(x) = −x;  (ii) F(−x,−y) = −F(x,y);  (iii) log_F is an odd series.
Interpretation: a global sign flip of the coordinate is complementation (σ ↦ −σ is T ↦ T^op in (a); the spectral reflection λ ↦ −λ, i.e. k ↦ m−k in THM-251, is the flip-all-tiles/complement-tiling involution — NOT T^op, per MISTAKE-033). The FGL's oddness is exactly the statement that the linearizer (rapidity, THM-252) **commutes with complementation**; THM-251's functional equation Q(λ_k)·Q(λ_{m−k}) = 1 is Q = exp(2·arctanh) applied to arctanh's oddness.

## Proofs

**(a)** For ±1-values, "exactly one of d, −d in S" ⟺ σ_S(d) ≠ σ_S(−d) ⟺ σ_S(−d) = −σ_S(d). For odd n, d ↦ −d is a fixed-point-free involution on the n−1 nonzero residues, so σ is freely determined on (n−1)/2 pair-representatives: 2^((n−1)/2) functions, and S ↦ σ_S, σ ↦ {d : σ(d)=+1} are mutually inverse. T^op has arc i→j iff j→i in T iff −(j−i) ∈ S iff σ(j−i) = −1: σ_{T^op} = −σ_T. Relabeling by u: arc x→y in uT iff u⁻¹(y−x) ∈ S: σ_{uT}(d) = σ(u⁻¹d). Even n: d = n/2 = −d kills both sides. The pairs-boundary statement is the definition of a tournament read off the skew sign matrix. ∎

**(b)** A completely multiplicative ±1-function on F_p^× is precisely a group homomorphism F_p^× → {±1}. F_p^× is cyclic of even order p−1, so there are exactly two such homomorphisms (the value on a generator g must be ±1 and determines everything): the trivial one and the one with kernel the unique index-2 subgroup = the squares, i.e. χ. Trivial is even. χ is odd ⟺ χ(−1) = −1 ⟺ −1 is a non-residue ⟺ p ≡ 3 (mod 4) (Euler's criterion). Combining with (a): multiplicative ∩ tournament = {χ} if p ≡ 3 mod 4 (Paley), ∅ otherwise. ∎

**(c1)** x_i ↦ −x_i is a bijection on distinct tuples sending each factor σ(d) ↦ σ(−d) = −σ(d), so A_L = (−1)^L A_L. ∎

**(c2)** Oddness gives Σ_{d≠0} σ(d) = 0 (pair cancellation). Integrate out the endpoint: Σ_{x_2 ∉ {x_0,x_1}} σ(x_2−x_1) = 0 − σ(x_0−x_1), so
A_2 = −Σ_{x_0≠x_1} σ(x_1−x_0)σ(x_0−x_1) = −Σ σ(d)σ(−d)·n(n−1)/(n(n−1)) ... per ordered pair the weight is −σ(d)σ(−d) = +σ(d)² = +1, hence A_2 = n(n−1). Identically this is −tr(M²) = ‖M‖²_F (M is skew, so −tr(M²) is the Frobenius norm — "the cherry is Parseval for the tournament symbol"). ∎

**(c3)** Möbius inversion over the partition lattice of the L+1 walk positions: A_L = Σ_π μ(0̂,π) S_π, where S_π identifies positions within blocks and sums the rest freely. (i) S_0̂ = 𝟙ᵀM^L𝟙 = 0 since M𝟙 = 0 (zero row sums = oddness). (ii) Any π merging adjacent positions has a σ(0) = 0 factor: S_π = 0. (iii) Among single-pair merges {i,j}, j ≥ i+2: unless {i,j} = {0,L}, position 0 or L stays free and its endpoint sum Σ_x σ(x−v) = 0 kills S_π; the survivor {0,L} gives S = tr(M^L) with Möbius sign −1. (iv) Every remaining π has ≥ 2 merges, hence ≤ L−1 blocks, so |S_π| ≤ n^(L−1) (weights bounded by 1); there are at most Bell(L+1) such π with |μ| ≤ L!. ∎

**(d)** Instance: F(x,y) = (x+y)·Σ_k (−xy)^k ∈ Z[[x,y]]; F(x,−x) = 0 directly, and uniqueness of the formal inverse gives [−1]_F = −x. The addition law tanh(a+b) = (tanh a + tanh b)/(1 + tanh a·tanh b) shows log_F = arctanh, which has only odd powers.
General lemma (R torsion-free, F commutative; recall 1-dim FGLs over torsion-free rings are automatically commutative, Lazard):
- (iii)⇒(ii): over R⊗Q, F = ℓ⁻¹(ℓ(x)+ℓ(y)); the compositional inverse of an odd series with unit linear term is odd (apply ℓ to −ℓ⁻¹(y)), so F(−x,−y) = ℓ⁻¹(−ℓ(x)−ℓ(y)) = −F(x,y); the identity descends to R[[x,y]] since R → R⊗Q is injective.
- (ii)⇒(i): h(x) := F(x,−x) satisfies h(−x) = F(−x,x) = −F(x,−x) = −h(x) by (ii), while commutativity gives h(−x) = F(x,−x)... precisely: F(−x,x) = F(x,−x) = h(x). So h = −h, and torsion-freeness forces h ≡ 0, i.e. [−1](x) = −x.
- (i)⇒(iii): apply ℓ to F(x,[−1](x)) = 0: ℓ(x) + ℓ(−x) = 0. ∎

## The honest boundary: what oddness does NOT give (recorded, OPEN)

Oddness localizes the **cherry**; it does NOT localize the **limit**. Three exact findings (`odd_function_dictionary_kpo2.py`, `..._limits_kpo2.py`, `circulant_H_census_kpo2.py`):

1. **Higher even generators are spectral, not odd.** For the rotation (carousel) family S = {1,…,(n−1)/2}, by (c3) the generators are sine-spectrum power sums, and EXACTLY (Fraction identity, k ≤ 6, via ζ(2k)/π^(2k) and Bernoulli numbers):
   **a_{2k}^rot = (−1)^(k−1)·2·(2/π)^(2k)·(1−2^(−2k))·ζ(2k) = [x^(2k−1)] tanh(x)**.
   So a_4 = −1/3, a_6 = +2/15: the cluster generators of the canonical non-quasirandom circulant family are the **Taylor coefficients of tanh — the compositional inverse of the formal-group logarithm of entry (d)**. Σ_k a_{2k} = tanh(1). (Numerics: −tr(M⁴)/n⁴ → −0.33330 at n = 401, Richardson −0.333334; A_4 exact = −225680 at n = 31.)
2. **The linked-cluster EXPONENTIAL formula fails outside the quasirandom class.** The naive prediction R_rot → exp(tanh 1) ≈ 2.14169 is contradicted by exact counts: R_rot(n) = 2.00000, 2.22222, 2.30476, 2.38646, 2.44113, 2.48658, 2.52227, 2.55197 at n = 5,…,19 — already past 2.14 at n = 9 and still climbing; SIS estimates (validated against the exact values) continue 2.61(3), 2.66(5), 2.71(9), 2.86(21) at n = 21, 31, 51, 101. Multi-run connected clusters are suppressed ONLY when the symbol's Fourier mass is spread (quasirandom/DRT — exactly THM-438's proved scope). **The limit of R_rot is OPEN** (it tracks the Paley sequence; a limit = e for the rotation family is consistent with the data and would mean e survives by a mechanism beyond the single-run expansion).
3. **The cherry is universal; the crown is not.** Exact H-census of ALL circulant tournaments (reduced by the dictionary's unit-relabeling × negation symmetries, which preserve H): the circulant H-maximizer is Paley at n = 7, 11 but the **ROTATION** tournament at n = 13, 15, 17, 19; at n = 19 Paley ranks only **3rd**:
   H_rot(19) = 1,184,212,824,763 > H_2nd(19) = 1,178,609,421,219 > H_paley(19) = 1,172,695,746,915.
   (Census cross-checks canon exactly: H_paley(11) = 95095 max; n = 9 maximizer = circulant {1,2,3,5}, H = 3357, as in THM-064.)

Scope notes. The bijection (a) is at the level of labeled circulant structures (connection sets), not isomorphism classes (units identify some). (b) needs prime p (composite Z_n has zero divisors; "completely multiplicative" is not the right notion there). The exact identity (c2) holds for all odd n including composite (verified n = 9, 13, 15, 17, 21 incl. random and block symbols).

## Reconciliation with prior repo statements

- **MISTAKE-011b** ("Paley at p ≡ 1 mod 4 is not a tournament") = entry (b): trivial/even vs odd character. The mistake is now a one-line classification theorem.
- **THM-061** already established the W-polynomial's odd/even parity in n (W(−r) = ±W(r)); cited, not rediscovered. LEM-004 (d) is the formal-group face of the same complementation symmetry.
- **THM-252** stated "arctanh is the formal group logarithm of the Cayley formal group F = (x+y)/(1+xy)"; **THM-414** used the tangent analogue (x+y)/(1−xy). Neither stated oddness; (d) adds [−1] = −x ⟺ log odd ⟺ F(−x,−y) = −F and verifies both FGLs are odd (logs arctanh, arctan).
- **THM-251's functional equation** Q(λ_k)Q(λ_{m−k}) = 1 is exp(2·arctanh) of arctanh's oddness; the spectral reflection k ↦ m−k is the complement-TILING involution (MISTAKE-033 caution: not T^op).
- **SCOPE FLAG 1 (THM-438 §Significance pt. 3).** The sentence "A_{2k} = C_k p^{k+1} and R(p) → e hold for EVERY circulant tournament" **overclaims**: the rotation family has A_4 = −n⁴/3 + O(n³) ≠ C_2·n³ (exact integers above). THM-438's own proof engine (M² = J − pI, ADDENDUM-3) is a DRT property; the correct scope is quasirandom/DRT circulant families — which is what HYP-2307's INDEX entry already says ("EVERY quasirandom circulant tournament"). Core theorem untouched. Flagged for an orchestrator MISTAKE entry; not edited here.
- **SCOPE FLAG 2 (reflection `why-the-paley-path-ratio-is-e...`, pt. 3).** "Paley higher only because it is the H-maximizer (A038375)" is false for p = 19 among circulants (boundary item 3). True at p = 7 (global max 189 = A038375(7)) and p = 11.
- **LEM-003** (Aut acts freely on Ham paths) is used by the census/limits scripts: H = n·#(paths starting at 0) for any circulant tournament.

## Verification

`04-computation/odd_function_dictionary_kpo2.py` → `05-knowledge/results/odd_function_dictionary_kpo2.out` (91 checks, 0 failures, exact int/Fraction):

| Entry | Check | Result |
|---|---|---|
| (a) | counts both ways = 2^((n−1)/2), n = 3,5,7,9,11; roundtrip; A+Aᵀ=J−I; T^op = −σ (ALL T); relabel action (ALL u, ALL T); n = 4,6,8 empty | all PASS |
| (b) | all p < 60: CM = {triv, χ} (generator argument), triv even, χ odd ⟺ p ≡ 3(4), QR-set tournament ⟺ p ≡ 3(4); exhaustive over all 2^(p−2) sign functions at p ≤ 19: exactly 2 CM | all PASS |
| (c) | 13 families (Paley 7–23, rotation 9–21, block, random): Σσ = 0, σ(d)σ(−d) = −1, A_1 = A_3 = 0, A_2 = n(n−1) = −tr(M²) exact; A_4 exact vs −tr(M⁴), E/n³ bounded (≈ 2.6–2.7, both Paley and rotation); Paley tr(M⁴) = p²(p−1) exact p ≤ 43; rotation −tr(M⁴)/n⁴ → −1/3, −tr(M⁶)/n⁶ → 2/15 (n ≤ 401 + Richardson); tanh identity exact k ≤ 6; H(paley) = 3, 189, 95095 and H(rot_13) = 3711175 = MISTAKE-011b's canon value; block_n = unit relabeling of rot_n (u = 5 at n=9, u = 7 at n=13) | all PASS |
| (d) | [−1] = −x to deg 21; only odd total degrees in F; log = arctanh exactly; round trip log↔FGL; integer coefficients; tangent FGL all-3-hold (log arctan); multiplicative FGL all-3-fail ([−1] geometric); 3 random odd logs hold, 3 random non-odd logs fail together | all PASS |

`04-computation/odd_function_dictionary_limits_kpo2.py` → `05-knowledge/results/odd_function_dictionary_limits_kpo2.out` (5 checks, 0 failures): fixed-start HK validated (n = 11, 13); exact H_rot(17) = 13,689,269,499, H_rot(19) = 1,184,212,824,763; SIS validated against 4 exact values, then n ≤ 101 estimates.

`04-computation/circulant_H_census_kpo2.py` → `05-knowledge/results/circulant_H_census_kpo2.out` (5 checks, 0 failures): full exact census n = 7,…,19 (2,4,4,6,16,16,30 classes), Paley max at 7 & 11, rotation max at 13–19, Paley 3rd at 19.

## References

- HYP-2378 (05-knowledge/hypotheses/INDEX.md) — the claim this lemma proves.
- THM-002/CONJ-001 (OCF, session context), THM-061, THM-251, THM-252, THM-414, THM-438 (+HYP-2307), THM-064, LEM-003; MISTAKE-011b, MISTAKE-033.
- 07-reflections/why-the-paley-path-ratio-is-e-the-cherry-is-the-unique-cluster.md; 04-computation/cluster_universality_monad.py (+.out).
