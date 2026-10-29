# Independent audit (audit2): THM-4846 (D′), THM-4845 (B) self-converse paragraph, HYP-9445 (1)–(2)

*Second independent adversarial audit by a subagent of mac-mini-2026-10-29-growth, 2026-10-29, covering the parts added
after the first audit was launched. The producer saved this report verbatim from the auditor's final message. The
auditor's scripts and outputs are in this directory.*

**Scope.** Worktree `/Users/e/Documents/GitHub/math-wt-pyth-cuboids`, branch `codex/session-pyth-cuboids` at `92514b16e2`. I edited no existing file and committed or pushed nothing. All scripts and outputs are new files in `04-computation/experiments/growth_gf_runs_fermat_20261029/audit2/`.

**Method.** All code is written from scratch and never imports `runs.py`:
- **Path counts:** a vertex-prepending path recurrence by length. The producer appends vertices, by subset.
- **Two- and three-run observers:** subset convolution. The producer cuts permutations.
- **Smith forms and determinants:** FLINT (python-flint) instead of PARI.
- **Self-converse census:** nauty `labelg` canonical forms.

## Summary verdict

I found no mathematical error and no counterexample. Every number I checked reproduced exactly.

- **Item 1, THM-4846 (D′): PASS, with three minor GAPs.** The algebra re-derives. The formula holds on every class through order 7 (the producer stopped at order 6), and the T5 claims are true.
- **Item 2, THM-4845 (B) self-converse: PASS on the derivation and numbers; CITE-NEEDED.** The condition and count are classical and already in the repo as THM-283. The identification of A002785 with self-converse tournaments is a documented fact, not just numerical agreement.
- **Item 3, HYP-9445 (1)–(2): PASS.** Everything in (a) and (b) reproduces on all classes of orders 3–8, with certificates that do not rely on Smith forms. There are two bookkeeping slips.

## Item 1 — THM-4846 (D′), mod-3 ghost lemma: PASS (minor GAPs)

Script `audit2/ghost3_audit.py`; outputs `out/ghost3_audit_n6.out`, `out/ghost3_audit_n7.out`.

**(a) 3-local image.**
- F rebuilt from the definition of ρ equals (R2). det F = 54, Smith form 1,1,1,3,3,6, and (R4)·F = 6I.
- The (R3) lattice has index 54 (counted mod 6), so the image of F is exactly the (R3) lattice.
- Q⁻¹WQ = [[c+d, c], [a+b−c−d, a−c]]. "Upper triangular mod 3 with diagonal (u, v)" therefore reads c+d ≡ u, a+b ≡ c+d, a−c ≡ v, which is (R3) mod 3. I checked this exhaustively on (Z/3)⁶.
- The reason: ρ(g)(1,1)ᵀ ≡ (1,1)ᵀ mod 3, and det ρ = sign.

**(b) Basis and reduction; the producer's conventions.**
- I checked the conventions against the data rather than assuming them:
  - the triple (b, g) is (b[g0], b[g1], b[g2]);
  - run permutations act by g ↦ g∘s, and M_3 is invariant under them;
  - L_ij = Σ_g k_ij(g) ρ(g⁻¹).
- N = F M F⁻¹ has exactly the predicted left-multiplication block form on all 528 classes with n = 4..7. Both columns of W = [[a,b],[c,d]] (coordinates (2,4) and (3,5)) are transformed by the same L. The right-multiplication alternative fails on all 528 classes.
- `ghost3_lemma.py` reads L from the (a, c) column and conjugates as Q⁻¹(·)Q with Q⁻¹ = [[0,1],[1,−1]]. That is correct.
- Exact check: let B be the stated Z_(3)-basis, ordered (u, y21, y11 | v, w12, y22). Then B⁻¹NB is 3-integral and congruent entrywise mod 3 to diag(I, II) on all 528 classes. L21 ≡ 0, L11 ≡ P and L22 ≡ S hold throughout.

**(c) The spur_3 formula.**
- On all 528 classes:
  - rank_F3 M = rank I + rank II;
  - rank_F3 N = rank P̄ + rank S̄ + 2 rank [[P̄, L̄12], [0, S̄]];
  - the stated spur_3 formula equals the direct rank difference.
- For n ≤ 6 (72 classes), the Smith-form excess of invariant factors divisible by 3 (N minus M) also equals spur_3.
- Classes with spur_3 &gt; 0: 0, 3, 21 at n = 4, 5, 6, as claimed. At n = 7 there are 271, with distribution {0:185, 1:64, 2:79, 3:51, 4:39, 5:17, 6:11, 7:4, 8:4, 9:1, 10:1}. This is identical to the producer's `census_n7_k3.out`.
- At n = 5: T5 has spur_3 = 4, 12223/H15 has 5, 02224 (source→C3→sink) has 3, and the rest have 0, matching the (E) table.

**(d) S̄ invertible, and T5.**
- Re-derived by elimination, with c = rank(Ē : ker P̄ → coker P̄):
  - rank I = r + 2 rank P̄ + c;
  - rank II = 2r + rank P̄;
  - rank [[P̄, L̄12], [0, S̄]] = rank P̄ + r.
- Hence spur_3 = c, with E = D1 − L12 S̄⁻¹ D2. This holds on all 122 classes with S̄ invertible mod 3 (1, 5, 20, 96 at n = 4..7).
- T5 (arcs i → i+1, i+2 mod 5):
  - M_3 has Smith form 1^52 3^2 6^2 0^4 and the unglued N has 1^48 3^2 6^6 0^4, both as in (R7).
  - The sectors are P: 1^6 0^4; S: 1^6 2^4 (det 16, so invertible mod 3); standard: 1^16 3^3 6.
  - P v_w = 0 over Z for all five stars.
  - The stars have rank 4 over Q and over F_3 (the relation is Σ(−1)^w v_w = 0), which equals dim ker P̄.
  - Ē maps their span onto coker P̄ with rank 4, so it is an isomorphism, and spur_3 = 4.

**GAPs (presentation only; none affects truth):**
1. The theorem lists the basis as (u, v, y11, w12, y21, y22). The displayed lower-triangular I and II need the orders (u, y21, y11) and (v, w12, y22). In the stated order, L12 sits above the diagonal of I. The correct order appears only in the script's docstring.
2. `ghost3_lemma.py` checks only the general rank identity. Neither the S̄-invertible reduction nor the T5 isomorphism claim was machine-checked by the producer; both are confirmed here.
3. The THM-4846 status block (its PROVED and FINITE-EXACT lists) does not mention (D′).

## Item 2 — THM-4845 (B), self-converse tournaments: PASS (math and numbers) + CITE-NEEDED

Script `audit2/sc_audit.py`; output `out/sc_audit.out`.

**Derivation: correct.**
- σ·T = Tᵒᵖ means x→y forces σ(y)→σ(x). The orientation therefore flips at each step around a pair orbit.
- A non-reversing orbit of length ℓ is consistent only if ℓ is even; a reversing (antipodal) orbit only if ℓ is odd.
- This excludes two fixed points, every odd cycle of length ≥ 3, and every cycle of length ≡ 0 mod 4. Orbits between two cycles have length lcm, which is even here.
- Each surviving orbit has 2 states, and the number of orbits is (mᵀGm − Σ_{i odd} m_i)/2 = (mᵀGm − m_1)/2.
- The text says "only if", but combined with "in that case 2^(…) such tournaments" it is in fact an iff.
- Twisted Burnside plus the moment rule E[W_1^{m_1} Π e^{m_j X_j}] = 2^{(mᵀGm − m_1)/2} gives the stated Gaussian form.

**Computations:**
- **Labelled brute force:** every σ ∈ S_n for n ≤ 6, and one σ per cycle type against all 2^21 labelled tournaments at n = 7. No violation of the condition or the count.
- **Pair-orbit walk:** every cycle type for n ≤ 14, no violation.
- **Unlabelled census:** all gentourng classes, comparing labelg canonical forms of T and Tᵒᵖ, n = 2..9. Self-converse classes: 1, 2, 2, 8, 12, 88, 176, 2752.
- **Cycle-type sum, n = 1..30:** equals the entire OEIS A002785 data field (n = 1..23, read 2026-10-29), which includes the claimed n = 1..16.
- **Gaussian expansion:** with the covariance computed from the explicit divisor realisation (Σ_{d|a, d|b} φ(d), no gcd used), n ≤ 23, it equals OEIS.
- **The producer's check** (`gf_gaussian.py` (3d)) evaluates the same gcd sum against typed-in OEIS terms. Those terms are correct, but the σ-condition is never tested.

**Identification of A002785 with self-converse tournaments: a known fact, not just numerical agreement.**
- OEIS titles A002785 as self-complementary oriented graphs. A 2020 comment by B. McKay says it also counts self-converse tournaments. Its a(2n)/a(2n+1) formula (Pab Ter, 2005) is algebraically the same twisted sum.
- A. Farrugia, *Self-complementary graphs and generalizations: a comprehensive reference* (M.Sc. thesis, Univ. Malta, 1999):
  - §5.2, p. 133, after noting that complement and converse coincide for tournaments: "They are also the self-complementary oriented graphs, under another name."
  - §§4.12–4.14 (pp. 106–107) and §5.3 Theorem C give the antimorphism criterion (even cycles of length ≢ 0 mod 4, at most one fixed vertex). It is credited to Zelinka (1970) and Salvi-Zagaglia (1979).
  - §7.32, p. 201, says that M. R. Sridharan, *Self-complementary and self-converse oriented graphs*, Indag. Math. 32 (1970) 441–447 (Zbl 0214.23505), counted the self-complementary oriented graphs, which are the self-converse tournaments, by exactly this cycle-type formula.

**CITE-NEEDED.**
- The paragraph states the condition and count without credit. They are classical (Zelinka, Salvi-Zagaglia, Sridharan, via Farrugia).
- They are also already repo canon:
  - THM-283 (opus-2026-04-03-S27, "Burnside Partition Formula for SC(n)") states exactly this condition and sum;
  - THM-849 (17) has the twisted Burnside sum;
  - `04-computation/a002785_enum.py` uses the OEIS/Howroyd formula;
  - THM-793 already records SC(n) = A002785(n).
- Only the Gaussian repackaging is new, and it is a one-line consequence of the moment rule. This is the same failure mode as MISTAKE-694.
- Suggested fix: add these credits and mark the paragraph as "classical count; Gaussian form new".

## Item 3 — HYP-9445 (1)–(2) data: PASS (two bookkeeping slips); no counterexample

Scripts `audit2/m2_torsion_audit.py`, `certify_torsion.py`, `k3_unique_hc.py`, `k3_hc2.py`. Outputs in `out/`: `m2_torsion_summary.out`, `m2_classes_n3.tsv` … `m2_classes_n8.tsv`, `certify_torsion.out`, `k3_unique_hc.out`, `k3_hc2.out`.

**Census, orders 3–8 (all classes):**

| order | classes | strong | unique-HC | HC ≥ 2 |
|---|---|---|---|---|
| 3 | 2 | 1 | 1 | 0 |
| 4 | 4 | 1 | 1 | 0 |
| 5 | 12 | 6 | 3 | 3 |
| 6 | 56 | 35 | 8 | 27 |
| 7 | 456 | 353 | 21 | 332 |
| 8 | 6880 | 6008 | 55 | 5953 |

The unique-HC counts equal F(2n−6) for n ≥ 4.

**(a) The signed sector S_2.**
- S_2 has torsion for 3/3, 27/27, 332/332 and 5953/5953 strong classes with HC ≥ 2.
- On unique-HC classes it is torsion-free for 1/1, 1/1, 3/3, 8/8, 21/21 (orders 3–7) and 54/55 at order 8.
- The single exception is exactly `1010111111111111111101111101`, with S_2 = 1^21 2^4 0^3.

**(b) M_2 on unique-HC classes.**

| order | det M_2 values (count) | nonsingular | torsion-free |
|---|---|---|---|
| 3 | −1 (1) | 1 | 1 |
| 4 | −1 (1) | 1 | 1 |
| 5 | +1 (3) | 3 | 3 |
| 6 | +1 (7), 0 (1) | 7 | 8 |
| 7 | −1 (16), 0 (5) | 16 | 16 |
| 8 | −1 (36), 0 (19) | 36 | 36 |

- So det M_2 ∈ {0, ±1}, the nonsingular counts are 1, 1, 3, 7, 16, 36, and the torsion-free counts are 1, 1, 3, 8, 16, 36.
- Torsion ⟺ singular holds at orders 7 and 8 with no violations.
- The one singular class at order 6, `101011111111101`, is torsion-free (1^26 0^4).

**Other (1)–(2) statements, all reproduced:**
- M_2 has torsion for all 6315 classes with HC ≥ 2.
- "Only odd invariant factors above 1" occurs for 8 classes at order 7 and 225 at order 8.
- P_2 torsion also occurs on unique-HC classes (1, 5, 19 at orders 6, 7, 8).
- Every strong class with unimodular M_2 has a unique Hamiltonian cycle (through order 8).
- k = 3:
  - unique-HC torsion-free counts are 1, 1, 3, 8, 16 for orders 3–7;
  - all 3 + 27 + 332 = 362 classes with HC ≥ 2 through order 7 have torsion (certified).
- A045891 (first differences of A045623) is 1, 1, 3, 7, 16, 36, 80, …, so the HYP remark is accurate.

**Independence and certificates.**
- Hamiltonian-cycle counts and Smith(M_2) agree class by class with the producer's PARI outputs: 456/456 in `fast_m2_n7.out` and 6880/6880 in `fast_m2_n8.out`.
- Every verdict is also certified without trusting Smith forms:
  - torsion by |det| &gt; 1 (nonsingular), or by an F_p rank below the rational rank;
  - torsion-freeness by |det| = 1, or by maximal minors with gcd 1.
- All 6315 HC ≥ 2 classes and all 89 unique-HC classes are certified, with no failures.

**Bookkeeping slips:**
1. "torsion-free through order 7 (32 of 32)" counts orders 5–7 only, because `p2_vs_hc.py` loops over n = 5, 6, 7. Through order 7 it is 34 of 34, since C3 and S4 are also torsion-free.
2. The "invariant factors such as 3·12³ and 2·6²" in HYP (2), "Where it breaks", are the k = 3 Smith forms (1^185 3 12^3 0^21 and 1^175 2 6^2 0^32). At k = 2 the same five order-7 classes have 4², 4², 3³·6, 2², 2². The sentence should state which k.

**Not checked (outside the listed sub-items):**
- the k = n−1 exceptions;
- k = 4;
- the order-8 k = 3 delicate classes;
- the order-9 and order-10 samples;
- the pre-registered order-9 prediction. Its value 80 coincides with the extrapolation (n+1)·2^(n−6) at n = 9.

## Files

All in `/Users/e/Documents/GitHub/math-wt-pyth-cuboids/04-computation/experiments/growth_gf_runs_fermat_20261029/audit2/`:
- `alib.py`: the independent library (path recurrence, M_2 and M_3, ranks mod p, FLINT Smith forms). Its self-test `test_alib.py` checks it against brute force for n ≤ 6–7; output `out/test_alib.out`.
- Item 1: `ghost3_audit.py` → `out/ghost3_audit_n6.out`, `out/ghost3_audit_n7.out`.
- Item 2: `sc_audit.py` → `out/sc_audit.out`.
- Item 3: `m2_torsion_audit.py`, `certify_torsion.py`, `k3_unique_hc.py`, `k3_hc2.py` → the `out/` files listed under Item 3.

Sources:
- [OEIS A002785](https://oeis.org/A002785)
- [Farrugia 1999 survey (PDF)](http://www.alastairfarrugia.net/sc-graph/sc-graph-survey.pdf)
- [zbMATH, Zbl 0214.23505](https://zbmath.org/3340122)
- [Sridharan 1978, Canad. Math. Bull., citing the 1970 paper](https://www.cambridge.org/core/journals/canadian-mathematical-bulletin/article/note-on-an-asymptotic-formula-for-a-class-of-digraphs/8A52DEEF5C24C6DEB684B76786BC0FAD)
