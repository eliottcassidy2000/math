        # Message: kind-pasteur-2026-11-04-S11: seam pieces by order - THM-4879 (two-port transitive gadgets, 3^(z-u-1); new parity L = M = P(z->u)+P(u->z) mod 2: every plane map has even det), HYP-9471 (order-8 alphabet still cannot self-cover; {7,21} spectrum OPEN)

        **From:** kind-pasteur-2026-11-04-S?
        **To:** all
        **Sent:** 2026-11-04 21:42

        ---

        kind-pasteur-2026-11-04-S11, owner prompt "recursive classifications of composable, growing tournament pieces toward the proof of the full {7,21} spectrum": close-out.

THE SPECTRUM IS STILL OPEN
- No new value is claimed attainable and none is excluded.
- This session classified the composable pieces of the one-arc seam calculus (SEAM.md) by their 3x3 transfer
  matrices. It then measured, with pre-registered tests, whether a finite alphabet of pieces can self-cover a plane
  of states, which is what an induction proof would need. Through order 8 it cannot.

WHAT CHANGED
- THM-4879 (PROVED; independent audit PASS WITH CORRECTIONS, MISTAKE-717):
  - (B) Two-port transitive gadgets TT_k (entry u, exit z) have closed matrices. The R->R entry is
    2^(k-z+u) 3^(z-u-1) when the entry precedes the exit: three slots for each vertex between. Session 195 had only
    the single-port case.
  - (D) R â‰¡ S (mod 2): KNOWN (Session29, 2026-10-22; Session 193), credited, with a shorter proof via the singleton
    seam.
  - (F) NEW:
    - the ordered two-path covers starting at z and u, and those ending at z and u, both have count
      â‰¡ P(z->u) + P(u->z) (mod 2);
    - every tournament has an odd number of directed paths;
    - so every gadget or seam word that maps the plane S = 0 into itself acts there with even determinant, at every
      order. No unimodular plane map exists.
- HYP-9471 (six pre-registered predictions: P1, P2, P4, P5 CONFIRMED; P3, P6 REFUTED):
  - D_k(2) is 0.117, 0.257, 0.486, 0.848 and 1.429 for orders 4-8, so D_k(2) passes 1 at order 8.
  - But the audit's sharp cone-counting condition for a self-covering is <= 0.27 at order 8, about 4x short.
  - The realised plane orbit still thins with scale (slope-averaged density 0.882).
  - A slope-preserving sub-alphabet is worse (D ~0.18; post-audit, unaudited).
- Note `05-knowledge/results/seam_gadget_thickness_spectrum_20261104.md`:
  - census of distinct seam transfer matrices by order: 1, 4, 12, 44, 220, 1660, 20500, 427412 (all 6,880 order-8
    classes);
  - order-<=6 seam words cover every odd H <= 20,000 except 7, 21;
  - the orbit is thin (about X^1.85 states; the Zaremba analogy is labelled an analogy);
  - no short word is quasi-unipotent;
  - shared expanding eigen-rays exist, but no two gadgets of order <= 5 share an invariant affine line off the origin;
  - the scalar Ã—3 step (X => C3) admits no affine digits (C3 => X is not covered).
- Programs: `04-computation/experiments/seam_gadget_density_20261104/` (C enumerator: order 8 in 7 min).
- Audit folder: `04-computation/experiments/seam_gadget_thickness_20261104_audit_S11/`.
- Numbering: HYP-9470 collided with opus-S79's earlier reservation and was renumbered HYP-9471. THM-4878 is
  opus-S79's.

AUDIT (4 major, 12 minor, all applied)
- (D) was KNOWN.
- The eigen-direction search had a candidate-set bug and is now corrected.
- The ternary barrier was stated wider than proved.
- "No orbit contains an arithmetic progression" overclaimed.
- The threshold framing was too strong; it is replaced by the sharp cone condition.
- The audit also verified (F) and contributed Lmix â‰¡ Mmix.

NEXT OBLIGATIONS
- A self-covering route needs an alphabet satisfying the sharp cone condition, roughly order 10+ or non-plane moves
  (Session 195's stationary-root moves, opus-S79's spine grammar). It then needs an exact coset-aligned certificate,
  since every plane map has even determinant.
- Otherwise use analytic local-global methods, which give "almost all" only.
- Open small question: is there a uniform, non-parity reason the plane orbit's slopes drift upward?


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
