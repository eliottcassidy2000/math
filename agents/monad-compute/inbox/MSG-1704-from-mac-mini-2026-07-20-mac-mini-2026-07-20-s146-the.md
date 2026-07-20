        # Message: mac-mini-2026-07-20-S146: THE TWO GMC(2) RESIDUALS -- radial Liouville step CLOSED for real p (charge-0 layer done for Hermitian P) via DvdK's jump argument transplanted to L; cross-shell residual shown ALL-INTEGER (the half-shell worry is a red herring by a parity lemma). THM-1675

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 15:41

        ---

        OWNER: 'work the 2 GMC(2) residual pieces.'

Setup: for 2 real Gaussians, Z = r w, W = Zbar = r w^{-1} (w = e^{i theta}, s = r^2), Gaussian measure (1/2pi) e^{-s} ds d theta. Every P in C[Z,W] is a Laurent polynomial in w with s-coefficients, Lambda_s(w) = sum_k w^k s^{|k|/2} lambda_k(s), and E[P^m] = L(CT_w(Lambda_s^m)) with L(f) = int_0^inf f e^{-s} ds, L(s^j) = j!. THM-1665 left two residuals.

PIECE 1 (HYP-8350, the radial Liouville step) -- CLOSED FOR REAL p.
THM-1665 reduced the charge-0/radial layer to: Psi(t) = int_0^inf e^{-v}[1/(1-tp(v)) - 1] dv == 0  =>  p == 0. For REAL p of degree D >= 1 this is now PROVED, by the DISCONTINUITY (jump) argument -- DvdK's Theorem-2 mechanism transplanted from the constant-term functional to L:
  Psi == 0 means h(t) := int_0^inf e^{-v}/(1 - t p(v)) dv == 1 on the Watson open set. h is analytic on C minus the curve {1/p(v) : v >= 0} (a subset of R for real p), which is CONNECTED, so by the identity theorem h == 1 everywhere off the curve. The boundary-value jump across the curve at t = 1/p(v0) is
     h(t+i0) - h(t-i0) = (2 pi i / t) * sum_{v: p(v)=1/t} e^{-v}/|p'(v)|,
  and h == 1 (analytic across) forces this to vanish. But every summand is STRICTLY POSITIVE, and a nonconstant polynomial is unbounded on [0,inf) -- contradiction. So p is constant, and L(p) = p = 0 gives p == 0. QED.
Exhaustive check: 0 nullcone members among 16800 nonzero real p (deg <= 4, coeffs [-3,3], m = 1..12); the jump density sum e^{-v}/|p'(v)| computed strictly positive on samples.
CONSEQUENCE: the charge-0/radial layer of GMC(2) closes whenever the effective radial polynomial is real -- i.e. for HERMITIAN P (P(Z,Zbar) real, lambda_0 a real polynomial). The COMPLEX-p case is DvdK's actual Theorem 2 (same jump argument, curve in C), cited not re-proved -- the one remaining radial sub-case.

PIECE 2 (HYP-8470, cross-shell coupling) -- ALL-INTEGER, framed precisely.
E[P^m] = L(CT_w(Lambda_s^m)) = sum_j L(s^j c_j(m)), where j = (sum|k_i|)/2 over charge-balanced m-tuples (sum k_i = 0) and c_j collects prod lambda_{k_i}.
PARITY LEMMA (proved, one line): sum|k_i| == sum k_i == 0 (mod 2), so j is ALWAYS an integer. Hence E[P^m] NEVER sees a half-integer moment Gamma(j+1/2) -- the 'half-shell' worry (which comes from the s^{|k|/2} in Lambda_s itself) is a RED HERRING: it cancels in the constant term by parity. Verified: zero odd-sum|k| balanced tuples for charge sets {-1,0,1}, {-2..2}, {-3,0,3}, m <= 5.
Span-2 (charges {-1,0,1}) decouples -- E[P^m] = sum_j m!/(j!^2 (m-2j)!) L(s^j (lam1 lam-1)^j lam0^{m-2j}), which THM-1600 eliminated (m=1 => lam0=0, then lam1 lam-1 = 0).
THE GENUINE COUPLING: at FIXED m, E[P^m] is a SINGLE number mixing all charge-degrees j -- there is no per-shell separation at fixed m; separation must come from varying m. The TOP shell (all charges +-k_max) has the fastest-growing moment (m k_max/2)!, dominating for large m -- the shell-descent handle.

HANDOFF.
HYP-8350 is CLOSED for real p -- the radial charge-0 layer of GMC(2) is done for Hermitian P. HYP-8470 is FRAMED, not closed. Two things GMC(2) still needs, both now sharp:
(i) THE COMPLEX-p RADIAL CASE. Transplant DvdK's Theorem-2 jump argument with the curve {1/p(v)} in C rather than R. The real case here is the exact template; the complex case is literally DvdK's own theorem, so this may be pure citation once written carefully.
(ii) THE CROSS-SHELL / CHARGE DESCENT for span >= 3 (and non-constant coefficients). At large m the top shell dominates, so E[P^m]/(m k_max/2)! -> (leading coeff of lambda_{k_max} lambda_{-k_max}); forcing that limit to 0 kills the top charge, then induct downward on k_max -- a CHARGE-descent parallel to the TNC coefficient ladder (THM-1610) but on charges rather than radial degree. UNTRIED, and the natural next attack.
SCOPE: settled now = span-2-constant (THM-1600) + the real radial layer (THM-1675); open = complex radial (DvdK-style) + the cross-shell descent. GMC(2) as a whole remains open.

Artifacts: THM-1675; 04-computation/gmc2_two_residuals_macmini_S146.py (+out); HYP-8350 closed (real), HYP-8470 framed.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
