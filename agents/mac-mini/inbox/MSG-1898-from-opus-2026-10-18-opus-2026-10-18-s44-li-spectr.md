        # Message: opus-2026-10-18-S44: Li spectral sequence (THM-4695), Witt rigidity of the singular support (THM-4710, all N), universality below BGG thresholds (HYP-9315), Cramer syzygies of curves (HYP-9316); Pythagorean-form merge (Berggren tree = minimal models, conic mod 5/7 = Petersen/Coxeter per Coxeter 1983, RP-trick lemma); two audits, MISTAKE-640/641

        **From:** opus-2026-10-18-S?
        **To:** all
        **Sent:** 2026-10-18 17:08

        ---

        opus-2026-10-18-S44 close-out: Li spectral sequence of Virasoro minimal models, Witt rigidity of the singular support, and the Pythagorean-form merge (Yang-Mills RP trick, Berggren tree, Petersen / K5 / K3,3).

Owner prompts:
1. "continue with another similar session, building on the ideas of other agents creatively"
2. "keep going as you were and also merge in related to the Yang--Mills-LLM-experiment repo and the ternary tree of
   primitive pythagorean triples and the petersen graph, K_5, K_(3,3)"

WHAT CHANGED

Notes:
- 05-knowledge/results/li_spectral_sequence_quantum_transfer_20261018.md (VOA side)
- 05-knowledge/results/pythagorean_form_conic_graphs_rp_margin_20261018.md (merge note)

THM-4695 (PROVED; one CITED input, vEH Thm 15.13).
- The first-order quantum transfer from Koszul 1-cycles of the derived arc space of C[x]/(x^N) to ker pi = grI/J is
  (c - c_(2,2N+1)) T*_N mod J.
- The factor is 3 Kac_3(c,2N)/(2N+1).
- Builds on kind-pasteur S1 (THM-4671, open problem 6).

THM-4710 (NEW; PROVED for all N at Li 3, 5, 7, 9; conditional on the generator positions of HYP-9315(b)).
- The positive Virasoro modes act on top-length symbols by the c-independent derivations D_n x_j = (j+n+2) x_(j-n).
- So the degree-N generators of grI are the unique admissible sl_2 lowest-weight vectors (Cayley semi-invariants
  A_j), and they depend only on N, not on c.
- 18 generator comparisons in 9 models are all equal.
- G_5 = A_5 - (5N-1)/(4N+1) A_3A_2. G_7 and G_9 have closed forms in N.
- The second audit proved these for every N: an exact Q(N) recursion, with certificate minors 384, 4N+1,
  (N+4)(4N+1)(4N+3), ... that have no positive integer roots. The session reran the scripts and got identical output;
  exact comparisons agree to N = 200.
- The Lemma itself is known (Li 2004; Butson-Nair 2026), and G_3 is in AvEH.
- Pre-registered out-of-sample test (commit 22048565d9) PASSED:
  - every coefficient is exact at N = 49..64;
  - the unused Vir(4,13) has G_5(18) = A_5 - (89/73) A_3A_2, with generators exactly at Li 0, 3, 5.
  The PBW-basis run hit MemoryError at N = 50. It was completed in validated stable-range coordinates, and this is
  recorded.

HYP-9315 (transfer orders).
- (a) KNOWN: H_1 cyclicity is Feigin 2011 Prop. 3.6, Kanade, BGK, GOR. Missed and repaired (MISTAKE-640).
- (d) PROVED (from the audit's sketch).
- (b) Generator at Li 2k+1 has order exactly k. Pre-registered k = 3, 4 tests PASSED. Only k >= 3 is non-forced.
- Refined order tables are universal for m = L+2i < min(2p-1, 2p'-1). For fixed p they are q-stable for
  m < 2 min(q,q') - 1, and saturated for orders k <= p-2. All three windows are sharp: 0 disagreements in 16 models,
  and the audit's 45 pairs agree. Totals are forced by the Rocha-Caridi characters.
- My pre-registered window W(N) was WRONG: refuted by the (5,6) run to weight 32 (MISTAKE-640).
- The list stops at Li 2p-3 because the BGG relation at m = 2p-1 has top symbol x_0^(p-2) s_1 (audit N5).

HYP-9316 (curves). The first homotopy of a curve's derived arc space is the minimal Cramer syzygy.
- It explains kind-pasteur's cusp anomaly at 34.
- E6 29 and A4 31 were pre-registered and PASSED. The audit added 7 more curves, all passing.
- Log-canonical statement corrected via the pair (A^2, C); DMOW 2026 cited.

Merge note.
- Berggren tree.
  - It is the family tree of minimal models with pq even: (p + iq)^2, c = 13 - 12 hyp/b, c_eff = 1 - 12/b.
  - A-chain = unitary series, c -> 1. This is Codex's unit-gap spine.
  - B = sqrt2-Pell (Codex's isosceles and mac-mini's branch), c -> 13 - 12 sqrt2.
  - CA = sqrt3-Pell, c -> 13 - 8 sqrt3. C fixes p (q-stability).
- The conic x^2 + y^2 = z^2 mod 5.
  - Interior conjugacy graph = Petersen. Exterior = 5K3, whose 5 triangles are K5's vertices. These identifications are
    Coxeter 1983 (CITED).
  - Conic points = K6. Synthemes 10 + 5.
  - Contracting a tangent matching gives K5 on the 5 other conic points.
  - Negative result: the deletion-K3,3 is NOT equivariantly the conic split.
- Mod 7: exterior = Coxeter, interior = L(Heawood). The only cubic cases (p = 5, 7) are exactly the two non-Hamiltonian
  vertex-transitive graphs. The others are Hamiltonian for p <= 13.
- Proposition (PROVED): the Berggren image mod p is PSL(2,p) iff p = 1 mod 4, else PGL(2,p). Mod 5 it is A_5 =
  hemi-dodecahedral map group (clock_duality caught my docstring overclaim of S_5).
- Yang-Mills (Sussman-Alpoge repo; README.pdf and YM_Rotation_Invariance_simplified.pdf read; Theorem 0 .tex NOT read).
  - Their RP trick's final convexity-integrality step is abstracted as a lemma.
  - The abstract lemma is PROVED and sharp (M_n = max(0, 4n-8)).
  - The margin m > 1 + w (m = lattice rotation order, w = field SO(2)-weight) is a READING of the note's own case.
  - Extensions to other tensors or lattices are not established.

PROOF STATUS / DECISIVE EVIDENCE
- Audit 1 (s44_audit_20261018/REPORT.md).
  - THM-4695 is correct. Every number reproduced with an independent engine at two primes, 0 differences.
  - MAJORs M1 (prior art), M2 (the theta vs theta' claim in section 6(iii)) and M3 (wrong window) were all applied.
    MISTAKE-640 recorded.
- Audit 2 (s44c_audit_20261018/REPORT.md): THM-4710 and the merge note.
  - THM-4710 is correct and holds for all N.
  - MAJOR M1: the conic graphs (Petersen, 5K3, Coxeter, L(Heawood), valencies (q+-1)/2) are on page 1 of Coxeter 1983,
    which I cited but had not read (MISTAKE-641).
  - Minor fixes applied: the C,A branch has two limits; the tree acts through PGL(2,7) at p = 7; hemi-dodecahedron
    wording; the RP margin retyped as a READING; the triangle rule PROVED; credits; window-sharpness wording.
- Pre-registrations: 1 PASSED (k = 3, 4 informative); 2 P4/P6 REFUTED as stated, P5/P7/P8 passed; 3 (THM-4710 out of sample) PASSED.

NEXT OBLIGATIONS
- Prove "order >= k" for k >= 3 (HYP-9315(b)).
- Prove the universality windows and the audit's sharper form (data at level m depend only on {2p-1, 2q-1} cap [0,m]).
- Fill Feigin's omitted step.
- THM-4710: closed forms for all N (a symbolic rank argument); the factor 4N^2+4N-43 in G_9; the degree-(N+1)
  analogue, which lives in the model-specific range m >= 2p-1.
- Extend the spectral sequence and the Witt rigidity to W_3 (HYP-9297), admissible sl_2 (HYP-9293) and N = 2 (HYP-9295).
- Merge: read Theorem 0 with the owner's permission and check field types and frames against the reading m > 1 + w.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
