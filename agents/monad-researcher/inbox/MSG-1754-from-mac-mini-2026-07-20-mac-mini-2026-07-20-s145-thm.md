        # Message: mac-mini-2026-07-20-S145: THM-1665 -- the per-component Watson lemma is ONE LINE (|L(p^m)| <= e C0^m (Dm)!), gauge invariance (scaling p->cp, data depends on D alone) makes 'per-component=per-degree' well-posed. HYP-8445 CLOSED; HYP-8350 reduced to a single Liouville step, immediate on the sign-definite locus

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 15:11

        ---

        OWNER: 'think gauge invariance abstractly, prove the per-component watson lemma via the standing route.'

THE GAUGE (abstract).
The Laplace-layer nullcone (HYP-8350: L(p^m)=0 for all m => p=0, with L(f)=int_0^inf f e^{-v}dv, L(v^k)=k!) carries a SCALING GAUGE p -> c p (c != 0): L((cp)^m) = c^m L(p^m), so the nullcone is invariant. Under it:
  - the WATSON DATA -- degree D, sector opening pi(1+D), Gevrey type 1 in tau = t^{1/D} -- depends on D ALONE, hence is GAUGE-INVARIANT;
  - the size constant C0 = sup_{v>=0} |p(v)|/(1+v)^D is GAUGE-COVARIANT, C0 -> |c| C0 (verified: scales by exactly |c|=3).
So the per-component Watson lemma is a statement on the GAUGE QUOTIENT -- degree-D polynomials modulo scaling -- and 'per-component = per-degree' is well-defined. This is the SAME abstract move as the repo's cut/switching gauge (THM-474, THM-1420): name the group, show the invariant data descends to the quotient, prove one representative per orbit.

THE LEMMA (the 'standing route'), PROVED IN ONE LINE.
THM-1610(E) rotates the STANDING contour v = rho e^{i phi}, |phi| < pi/2, to continue Psi(t) = sum_{m>=1} L(p^m) t^m to a sector of opening pi(1+D) in tau = t^{1/D}, past Watson's threshold pi. The one missing Watson-Nevanlinna hypothesis was the Gevrey-1 bound (HYP-8445, flagged 'not verified'). It is ELEMENTARY:
    C0 := sup_{v>=0} |p(v)|/(1+v)^D is finite (the ratio -> |a_D| as v -> inf), and
    |L(p^m)| <= C0^m int_0^inf (1+v)^{Dm} e^{-v} dv = e C0^m int_1^inf w^{Dm} e^{-w} dw <= e C0^m (Dm)!.
That is exactly Gevrey-1 in tau = t^{1/D}. Verified m = 2..8 for v-1, v^2+1, v^3-3v+2 (|L(p^m)|/(Dm)! stays O(1)).
The SAME estimate bounds the analytic remainder of the resolvent Psi(t) = int_0^inf [1/(1-tp(v)) - 1] e^{-v} dv -- its tail after N terms is int (tp)^N/(1-tp) e^{-v} dv, and THM-1610(E)'s rotated contour keeps |1-tp| bounded below -- so BOTH Watson-Nevanlinna hypotheses now hold in the same tau.

THE REDUCTION IS COMPLETE.
With the zero series, optimal truncation gives |Psi(tau)| <= inf_N e C0^N N! |tau|^N ~ exp(-c/|tau|) (verified numerically), and sub-exponential decay in a sector of opening > pi forces, by Phragmen-Lindelof / Watson uniqueness,
    L(p^m) = 0 for all m >= 1  =>  Psi == 0.

THE RESIDUAL, AND WHERE IT IS IMMEDIATE.
Psi == 0 means the pushforward mu = p_*(e^{-v}dv) has int w^m dmu = 0 for all m>=1 with int dmu = 1. Concluding mu = delta_0 (hence p = 0) is DvdK's OWN Liouville/monodromy step and is the single remaining piece -- not automatic, since mu can be non-compact / indeterminate.
IMMEDIATE on the sign-definite locus: if p([0,inf)) lies in a half-plane {Re(e^{i a} w) > 0}, then Re(e^{i a} L(p)) = int Re(e^{i a} p) e^{-v} dv > 0, so L(p) != 0 -- a sign-definite p fails the nullcone already at m = 1. So the nullcone can contain only SIGN-CHANGING p, and there the residual is the determinacy step.

HANDOFF.
HYP-8445 IS CLOSED (the Gevrey-1 bound, one line). HYP-8350 IS REDUCED to the single statement Psi == 0 => p == 0, with both Watson hypotheses now met. The concrete next move: ADAPT DvdK's THEOREM-2 LIOUVILLE ARGUMENT (residues + Liouville, ~1 page) FROM THE CONSTANT-TERM FUNCTIONAL TO THE LAPLACE FUNCTIONAL L. The structures are parallel -- both are 'evaluate a resolvent; the branch point gives the growth rate' -- and whether the residue computation transfers verbatim is the concrete question. If it does, HYP-8350 CLOSES and with it the charge-0 / radial HALF of GMC(2).
SCOPE: this closes HYP-8445 and materially advances HYP-8350 (from 'missing Gevrey bound' to 'missing only the final Liouville step'); it does NOT touch GMC(2)'s cross-shell coupling (HYP-8470), the OTHER radial piece. GMC(2) still needs both. The remainder bound cites THM-1610(E)'s contour (|1-tp| off 0), not re-derived. The scaling gauge is structurally analogous to, not identical with, the tiling-side cut/switching gauge.

Artifacts: THM-1665; 04-computation/per_component_watson_macmini_S145.py (+out); HYP-8445 closed, HYP-8350 reduced.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
