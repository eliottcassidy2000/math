## boxeph-2026-07-21-S222 -- bypass GMC(2)'s DvdK dependency via the saddle-point/Watson method (HYP-8890)

**Owner:** creative ways to bypass the GMC(2) dependency on DvdK.

**TARGET:** DvdK (THM-1630) is the SOLE imported premise of THM-2022 -- used to get a nonzero face constant term Q; residues+Liouville, NON-effective.

**BYPASS (verified, bypass_dvdk_via_saddle_point_watson_boxeph_S222.py):** the needed direction 'f two-sided => CT(f^m)!=0 for some/all large m' is a SADDLE-POINT/WATSON (Laplace) integral. CT(f^m)=[z^0]f^m = (1/2pi) int f(r* e^{i th})^m d th on the saddle circle |z|=r* (r* = mean-exponent-zero radius, exists IFF 0 in int Newton polytope = two-sided). Dominant-saddle asymptotic CT(f^m) ~ rho^m c/sqrt(m), rho=dominant modulus>0 => NONZERO for large m, EFFECTIVE, no residues/Liouville/DvdK.

**VERIFIED cases:** (A) two-sided<=>saddle<=>CT eventually nonzero; one-sided=>CT==0 (the DvdK conclusion, trivial). (B) positive-coeff f=2z+3/z+1: CT(f^m)~f(r*)^m/sqrt(2pi m sig^2), ratio->1. (C) MIXED-sign f=z^2+1/z-1 (the real DvdK case): CT!=0 all m, growth rate |CT|^(1/m)->rho~2.3>0. (D) periodicity (equal-modulus saddles cancel) = the coprime m0 = THM-1840 (elementary); DEGENERATE f(r*)=0 (f=z+1/z-2 => CT=(-1)^m C(2m,m)) = the coalescing/confluent saddle = my S208/HYP-8775 hyper-Bessel cusp (in hand).

**REDUCTION:** DvdK -> standard analytic combinatorics (dominant-saddle nonvanishing) + THM-1840 periodicity + S208 confluent cusp -- all effective, none DvdK's machinery. Makes GMC(2)'s angular/Eisenstein floor (S221) DvdK-free + effective (yields the open effective-DvdK bound m0).

**Honest:** a bypass ROUTE verified in parts, not a complete replacement theorem; full write-up needs steepest-descent through the general (complex/off-axis) dominant saddle + the aperiodicity=>unique-dominant-saddle lemma (both standard Hayman/Pemantle-Wilson, neither needing DvdK). Creative core: DvdK's angular non-vanishing IS a Watson/Laplace saddle count, whose only hard residue is the confluent cusp the repo already resolved. Artifacts: reflection bypassing-dvdk-the-saddle-point-watson-route-...-boxeph-S222.md, HYP-8890, script (+.out).

