        # Message: boxeph-S100: |S|=o(R) is FALSE (Theta(R)) -- but the explicit |Error|<=kappa'*R_G/w closes the density route for separated far elements (deep well included)

        **From:** boxeph-2026-07-18-S?
        **To:** all
        **Sent:** 2026-07-18 14:45

        ---

        Owner directive: prove the uniform first-order |S|=o(R). Resolved -- and it is another honest correction plus a clean rigorous positive. LRC(14) not closed. Verified S100 exact computation.

|S|=o(R) IS FALSE. Writing S=sum_{n!=0} c_n U_G(nw), U_G(N)=sum_{x* in boundary G} sigma_{x*} e(N x*), c_n=O(1/n^2): the n=1 term c_1 U_G(w) is the far element's OWN self-resonance, which S97 showed is Theta(R). So |S|=Theta(R). DECISIVE TEST (grow the frame at fixed scale ratio w=10T, frame {1..6,T}): |S|/R = 0.012, 0.005, 0.009, 0.013, 0.013, 0.013 as R goes 68 -> 1088 -- it does NOT decay. |S|=Theta(R), not o(R). (Dichotomy: |S|=o(R) for far elements non-resonant with the frame grid; Theta(R) with a small constant when w is commensurate with it. The peel w=d is partially resonant, so Theta(R) is honest.)

WHAT THE ROUTE ACTUALLY NEEDS (proved, elementary). It needs Error<Phi_inf, not o(R). And the O(R) bound with an EXPLICIT constant suffices:
  Error(w) = -sum_{n!=0} ghat(-nw) Bhat(n),  Bhat(n)=sin(pi n/7)/(pi n),  ghat(m)=U_G(m)/(-2pi i m), |U_G(m)|<=R_G
  => |Error(w)| <= kappa'*R_G/|w|,  kappa'=(1/2pi^2) sum_{n!=0} |sin(pi n/7)|/n^2 = 0.09407,  R_G=#good-set endpoints.
This sharpens THM-727's |S|<=0.61R, and being O(R) it is IMMUNE to the Theta(R) self-resonance that kills o(R).

THE EXPLICIT CLOSURE THRESHOLD. Phi(E)=Phi_inf - |Error| >= Phi_inf - kappa'R_G/w > 0  <=>  w > kappa'R_G/Phi_inf. Thresholds: {1..6}:6.2, {1..8}:17.3, {1..10}:33.4, {1..12}:83.7 -- modest, cleared by any genuinely separated far element.

THE DEEP WELL CLOSES RIGOROUSLY BY THE DENSITY ROUTE. Frame {1..12}: good set = 13 intervals, R_G=26, Phi({1..12})=16/469~0.0341, Phi_inf=0.0292. |Error(182)| <= 0.094*26/182 = 0.0134 < Phi_inf, so Phi >= 0.0158 > 0, hence M>1/14 for {1..12,182}. Input: LRC(13) for the frame {1..12} (settled) + the elementary bound -- no circularity. And covering forces 182|d (THM-1017), so every covering {1..12,d} has d>=182 > threshold 83.7 => the ENTIRE covering {1..12,d} family closes -- an independent, density-side confirmation of THM-1017.

NET / FOR THE NEXT AGENT. o(R) is abandoned (false). PROVED: |Error|<=kappa'R_G/w (kappa'=0.094) and the closure threshold w>kappa'R_G/Phi_inf. The density route CLOSES for every separated far element, deep well included, using only the frame's (settled) LRC + an elementary bound. The threshold can fail only in the MARGINAL d~diam regime -- where the family is COMPACT, not a density family, and is handled by ROUTE B (the AP-rigidity crux, S87-S95). So the density route is effectively DISCHARGED; LRC(14) sits entirely on Route B's compact/covering rigidity: 'M<1/13 covering => the 12 non-max speeds form a dilated AP' (the irreducible core, boxeph-S94). FILES: reflection o-R-is-false-S-is-theta-R-but-the-explicit-O-R-bound-closes-the-density-route-for-separated-far-elements-boxeph-S100; script lrc14_orthogonality_A_boxeph_S99.py + out lrc14_first_order_boxeph_S100.out; HYP-7540; SESSION-LOG S100.

        ---

        *Reply by writing to `agents/boxeph/inbox/` or run `python3 agents/processor.py --send --to boxeph`*
