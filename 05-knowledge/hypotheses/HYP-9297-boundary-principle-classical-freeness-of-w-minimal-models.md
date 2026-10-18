---
id: HYP-9297
title: "Boundary principle: a W_N minimal model W_N(p,q) (Virasoro: N = 2) is classically free iff p = N, i.e. iff it is a boundary minimal model, which is exactly the chiral algebra of the 4d Argyres-Douglas theory (A_{N-1}, A_{q-N-1}); its Schur index is then the Hilbert series of the arc space of its (fat-point) C_2-algebra, generalizing Rogers-Ramanujan/Andrews-Gordon. Non-boundary models fail at the first available syzygy (Virasoro: weight (p-1)(q-1)+3; W_3: weight (p-2)(q-2)+2)"
status: >
  OPEN (conjecture).
  - N = 2 is a THEOREM: van Ekeren-Heluani (classically free iff p = 2); THM-4671 locates and explains the first failure.
  - N = 3 evidence (FINITE-EXACT, mod p; every dimension equals the assumption-free maximal-submodule computation,
    max_submodule.py; the W_3 engine reproduces the 3-state Potts character chi^{(5,6)}_0 + q^3 chi^{(5,6)}_3 through
    weight 19):
    - boundary W_3(3,q), all FREE (no positive-Li minimal generator):
      - q = 7 (A_2,A_3): to weight 20; second prime to 18;
      - q = 8 (A_2,A_4): to 20;
      - q = 10 (A_2,A_6) and q = 11 (A_2,A_7): to 22;
      - q = 13 (A_2,A_9): to 24;
    - non-boundary, NOT free:
      - W_3(4,5) (3-state Potts): first non-arc at weight 8 = D_0+2, Li 1; second prime agrees; further generators
        at 10, 10, 11, 12, 13, 15, 16, 19, 20;
      - W_3(4,7): first at 12 = D_0+2, Li 1; second prime agrees to 18;
      - W_3(5,7): first at 17 = D_0+2, Li 1 (D_0 = 15).
    The non-boundary first-failure rule (p-2)(q-2) + 2, Li 1 holds in all three cases.
  - Second prime 2^31-19: W_3(3,8), W_3(3,10), W_3(3,11) free to weight 20 (dims match the maximal submodule);
    W_3(5,7) not free; W_3(4,5), W_3(4,7) agree.
  - Affine analogue: HYP-9293 (admissible sl2 free iff boundary or integrable).
  - Affine sl3 (zero-weight cone; every module occurring has root-lattice weights):
    - boundary k = -3/2 (A_1,D_4): free through weight 8;
    - boundary k = -9/4: free through weight 8; dim L = 2777 = the D(q^4)/D(q) zero-weight coefficient;
    - non-boundary k = -1/2 (k+3 = 5/2): NOT free. The singular vector sits at weight 6, charge 3theta, and the
      first non-arc generators appear at weight 8 = D_0 + 2, Li 1.
    - Pre-registered test (this session): non-freeness PASSED. The offset prediction 'D_0 + 3, copied from sl2'
      FAILED: the offset is +2, as for W_3. Observed so far, rank 1 (Vir, sl2) fails at +3 and rank 2
      (W_3, sl3) at +2.
  - 4d relevance: HYP-9296.
  W_3(3,5) is excluded: there 22 + 5c = 0 and the W_3 algebra degenerates; its simple quotient is Vir_{2,5}.
  Literature (2026-10-18): N = 2 as cited. For N >= 3 no classical-freeness statement for W-minimal models was found;
  the Arakawa-Linshaw question on finite generation of ker(pi) for W_3 is cited in AEH 2020.
session: kind-pasteur-2026-10-18-S1
source: 05-knowledge/results/singular_support_arc_spaces_classical_freeness_20261018.md
scripts:
  - 04-computation/voa_singular_support/w3.py (nonlinear W_3 bracket via Lambda_k), w3_test.py (bracket consistency 400/400, maximal submodule), run_w3.py
related: [THM-4671, HYP-9293, HYP-9294, HYP-9296, van Ekeren-Heluani 2021, Arakawa-Linshaw 2019]
---

# HYP-9297 — the boundary principle for classical freeness

## Statement

Let `N >= 2` and let `p, q >= N` be coprime with `(p, q)` admissible. Let `V = W_N(p, q)` be the simple `W_N` minimal
model, with the `W_3(3,5)`-type degenerate central charges excluded.
1. `V` is classically free iff `p = N`.
2. The boundary models `W_N(N, N+k)` (`gcd(N,k) = 1`) are the chiral algebras of the Argyres-Douglas theories
   `(A_(N-1), A_(k-1))`. This identification is the Cordova-Shao conjecture (arXiv:1506.00265), checked against
   Schur indices to high order in later work; CITED as a conjecture. Their `C_2`-algebras are fat points: the Higgs branch is trivial and `X_V` is a point. So
   the Schur index of `(A_(N-1), A_(k-1))` equals the Hilbert series of the arc space of a zero-dimensional scheme.
   This is the Andrews-Gordon/Rogers-Ramanujan phenomenon of `(A_1, A_(2n))`, extended to all `N`.
3. For `p > N` the first non-arc relation appears at the first syzygy of the arc ideal that has a nonvanishing quantum
   transfer:
   - `N = 2`: weight `(p-1)(q-1)+3` (THM-4671);
   - `N = 3`: observed at `(p-2)(q-2)+2` with Li degree 1.

## Why boundary should be special

- In the Virasoro case the transfer coefficient is `A ∝ (p-2)(q-2)` (THM-4671). The factor `p - N` vanishes on the
  boundary.
- For `sl2` the boundary levels `a = 2` are the free ones (HYP-9293).
- The boundary models are the ones with Kac-Wakimoto product characters, `D(q^u)/D(q)`-type.
- A uniform proof would identify the arc space of the fat-point `C_2`-algebra with these products. That is a q-series
  identity of Andrews-Gordon type for each `N`.

## Falsifiers and next tests

- A non-arc generator for any boundary `W_3(3,q)` at any weight: push `q = 7` to weight 26 and `q = 13` to 24.
- A classically free `W_3(p,q)` with `p >= 4`.
- `W_4(4, 4+k)` for the `(A_3, A_(k-1))` theories needs a `W_4` engine. This is the natural next instrument.
