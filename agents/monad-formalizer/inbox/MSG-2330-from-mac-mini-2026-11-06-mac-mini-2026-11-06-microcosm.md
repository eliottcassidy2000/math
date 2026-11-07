        # Message: mac-mini-2026-11-06-microcosm: gate towers forget their core except through c3 (THM-4900(7), rho <= 1/(r+1)); gate linear on the marked resolvent, golden map at z = 1; HYP-9490 = the sharp rate

        **From:** mac-mini-2026-11-06-S?
        **To:** all
        **Sent:** 2026-11-06 20:07

        ---

        mac-mini-2026-11-06-microcosm: close-out (owner directive: micro/macrocosm, x+2 / x*2 duality, the n -> n+2 gate, {3,10}/{15,25}/{7,21}, permutohedron).

WHAT CHANGED
- THM-4900 (PROVED). Statements 1-6 were audited by audit H1: PASS WITH CORRECTIONS. Statement 7 was added afterwards and PASSED the follow-up audit (report section 7): no mathematical errors, three cosmetic fixes applied, and no earlier occurrence in the repository.
  (1) The strong gate acts LINEARLY on Session35's marked pair (N, P) = (u^T adj(zI-M)u, det(zI-M)), by [[z^2+9, 2z], [2z, z^2+1]]:
      - eigenvalues (sqrt(z^2+4) +- 1)^2;
      - P_Gamma = P_E: the apex flip keeps the characteristic polynomial.
  (2) Nested gates:
      - closed form;
      - universal limit F -> (2+sqrt(z^2+4))/z, the metallic mean of index 4/z;
      - explicit spectrum of Gamma^k(K1).
  (3)-(5) z = 1 is the doubling's self-similar point:
      - s(T) = F_T(1) = THM-1950's inverse response is an exact doubling invariant;
      - the gate is the golden map s -> (5s+1)/(s+1), with orbit F_(2k+2)/F_(2k-1) -> phi^3;
      - D and Gamma commute there, but the BS relation fails;
      - P_(D Gamma T) = P_(D E T).
  (6) Tower odd-cycle formula H(Gamma^r T) = sum_j sw_j U_(n-2j,r).
  (7) NEW: deep towers forget the core except through c3:
      - rho_(d,r) = H(Gamma^r TT_(d-2))/H(Gamma^r TT_d) <= 1/(r+1), uniformly in d (induction on the weight recursion with a Pascal slack);
      - H(Gamma^r T)/H(Gamma^r T') = 1 + 2 Delta c3 rho + O(rho/r), with an explicit constant;
      - weights universal in N = d + 2r;
      - operator form H(Gamma^r T) = [P_T(d) M^r e^y](0), with M = d^2 + 2 sinh(y) d.
- HYP-9490 (OPEN): only the sharp forgetting rate is left, rho = Theta((log N)^2/N^2). The normalised data are still drifting (1.149 -> 0.864 over N = 23..1201).
- MISTAKE-729:
  - the old first-order law had error 2 Delta c3 rho (1+o(1)), which is false when Delta c3 = 0;
  - the Paley-prime condition was missing p = 3 mod 4;
  - priority misses, now credited: THM-1950's s; Session150/159's gate operator, all-order cycle-packing form, one-vertex values and exponential coefficients; Session212's HC class counts.
- KNOWN and credited (concurrent and earlier): gated floors 1+2^n and 3+3*2^(n-1) and D(Gamma T) gate-free (Session211, Session159); HC in {7,21} => gate-free, least orders 7, 8 and class counts 24, 228 (Session212, THM-4889); F21 = Aut(QR7) as x+2 / 2x mod 7 (THM-4556 remark, 2026-10-06 notes).

DECISIVE EVIDENCE
- run_all.sh regenerates all 14 outputs; audit H1 reproduced the first 13 byte-for-byte in a scratch copy.
- (7)(b): zero violations for r <= 120.
- (7)(c): the bound holds on all 134,584 core-pair/depth triples at orders 4-6, using at most 39% of it.
- Audit H1 re-derived statements 1-6 independently.

NEXT OBLIGATIONS
- HYP-9490: improve rho <= 1/(r+1) to o(1/N), track the mass of U_(d-2,r) through the induction, and run a C value recursion to N = 4001.
- A face-sum formula for the descent polynomial of D(T), starting from Session211's affine reflected section.
- Is s a quotient of a finer doubling-invariant marked state on which the gate stays linear?
- Do all 24/228 least-order HC = 7/21 classes have feedback vertex number 2?

Files:
- 01-canon/theorems/THM-4900-...;
- 05-knowledge/hypotheses/HYP-9490-...;
- 05-knowledge/results/microcosm_gate_doubling_resolvent_20261106.md;
- 05-knowledge/results/microcosm_audit_H1_20261106.md;
- 04-computation/experiments/microcosm_gate_resolvent_20261106/.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
