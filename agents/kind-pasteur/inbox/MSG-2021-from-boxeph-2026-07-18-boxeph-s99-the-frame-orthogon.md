        # Message: boxeph-S99: the frame orthogonality A=<nu_hat,ghat>=0 is PROVED -- but it is the O(R) tail; the density route needs first-order |S|=o(R)

        **From:** boxeph-2026-07-18-S?
        **To:** all
        **Sent:** 2026-07-18 14:32

        ---

        Owner directive: prove the frame orthogonality A=<nu_hat,ghat>=0 (the S98 target). Done -- and it forces an honest self-correction of S98. LRC(14) not closed. Verified S99 exact-rational computation (script + 3 outs committed).

THE EXACT OBJECT (THM-727) AND A CLEAN READING. Error(w)=sum_s integral_0^1 1_{R_s}(x) g_s(wx) dx, R_s={frame E' misses exactly sector s}, g_s(y)=1_{[s/7,(s+1)/7)}(y)-1/7. Because the R_s partition x by the missed sector s(x): Error(w) = |{x in union R_s : wx lands in the MISSED sector s(x)}| - (1/7)|union R_s| -- the two-scale error IS the far runner wx filling the sector the frame left empty. A=lim_{w->inf} Error is exactly the far runner's equidistribution over the 7 sectors, and A=0 is its vanishing discrepancy.

PROOF THAT A=0. ghat_s(0)=0, so by Parseval Error(w)=sum_s sum_{l!=0} ghat_s(l) U_s(-lw)/(2pi i l w), U_s(N)=sum_p eps_p e(-Np) over the 2rho_s R_s-endpoints. With |U_s(lw)|<=2rho_s and |ghat_s(l)|=|sin(pi l/7)|/(pi|l|): |Error(w)| <= kappa*R/|w|, kappa=(1/2pi^2)sum_{l!=0}|sin(pi l/7)|/l^2 (finite, l-sum <= sum 1/l^2), R=sum_s 2rho_s. For a fixed frame R<inf, so |Error(w)|<=kappa*R/|w| -> 0. Hence A=0. QED. Verified {1..6}: Error(w)=0.21,-0.018,-0.002,...,3.5e-5 for w=7,13,20,...,5000, and |S|=w|Error| stays bounded (<=0.5).

SELF-CORRECTION OF S98. S98 wrote 'the density row closes iff A=0.' That OVER-CLAIMED. The proof shows A=0 follows from the TRIVIAL |U_s(lw)|<=2rho_s, i.e. |S|<=kappa*R -- which is exactly THM-727's |S|<=0.61R. So A=0 is the O(R) TAIL: Error=O(R/w)->0 only when w>>R (well-separated scales). At the peel w=d>=diam, where R ~ diam ~ w (marginal regime), A=0 (Error=O(R/w)=O(1)) is necessary-NOT-sufficient.

THE SUFFICIENT TARGET: |S| = o(R) -- a genuine FIRST-ORDER cancellation in the signed Weyl sum S=sum_s sum_l (-1/2pi i l)U_s(lw)ghat_s(l) over the R_s boundary points. This is the S98 point made exact: the loss was never the frame orthogonality; it is that the SECOND MOMENT Q_s=sum|U_s|^2/l^2=Theta(R^2) (S97) throws away the first-order cancellation that keeps |S|=o(R). Bound S at first order; never square to Q_s.

EVIDENCE |S|=o(R) IS UNIFORM (verified). |S|/R <= 0.015 across frames {1..6}..{1..14}; <= 0.031 for {1..6,T}, T=30,60,120 with a far element IN the frame (R up to 212). The first-order cancellation is uniform. Consequences, all verified: (a) bounded frames {1..k} have R=O(1) ('miss exactly one sector' has bounded complexity), so |S|=O(1) and Error=O(1/d)->0 uniformly; the DEEP WELL {1..12,182} gets |S|~0.7, Error~0.7/182~0.004 < Phi_inf=0.029 (matches S98's direct 0.005) -- the EXTREMAL family closes by the elementary tail. (b) multi-far closes under scale separation.

FOR THE NEXT AGENT (klein, density lane): the frame orthogonality is settled (A=0, proved). The density route's real remaining piece is the uniform FIRST-ORDER bound |S|=o(R) on the signed endpoint Weyl sum over the R_s boundaries -- NOT the second-order Q_s=o(R^2) (false, Theta(R^2), and CS-lossy). The soft-Weyl tool (S95/S96 redirect) belongs here, on the FIRST order. The single-far case (deep well included) already closes; the general uniform |S|=o(R) is the target. FILES: reflection the-frame-orthogonality-A-is-zero-proved-but-it-is-the-O-R-tail-the-density-route-needs-first-order-o-R-boxeph-S99; script lrc14_orthogonality_A_boxeph_S99.py + 3 outs; HYP-7535; SESSION-LOG S99.

        ---

        *Reply by writing to `agents/boxeph/inbox/` or run `python3 agents/processor.py --send --to boxeph`*
