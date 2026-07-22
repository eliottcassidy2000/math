## boxeph-2026-07-21-S225 -- working JC(2): the obstruction is a descent termination; coprime intervals are its tool (HYP-8905)

**Owner:** long session working to prove the planar Jacobian Conjecture; pull past threads creatively.

**HONEST:** JC(2) OPEN (JC(n>=3) FALSE, Keller THM-1300). Worked it via 3 reduction PROGRAMS (none complete): (A) de Bondt+Zhao => VC(4) (JC(2)<=>a dim-4 Laplacian moment nullcone, doubling n->2n); (B) klein-S329 Euler-Zariski cover-degree-3 bootstrap (cuspidal Jelonek curve, ramification-parabola escape, NO CF); (C) mac-mini-S137 golden-corner/Lame (subtractive Euclid on Newton slopes).

**VERIFIED (jc2_via_vanishing_conjecture_and_the_cf_termination_boxeph_S225.py):** 2D symmetric case EASY -- nilpotent Hessian <=> P prop (x+iy)^d (harmonic), Delta^m(P^m)=0, invertible. Lame worst-case = Fibonacci (longest Euclid chain <200 = (144,89)). [Restricted JC(2) proved elsewhere: THM-1345 equivariant, THM-1370 elliptic, THM-1365 poly-Galois, geom-degree<=2, THM-2063.]

**CREATIVE UNIFICATION:** the JC(2) obstruction has THREE equivalent forms -- (A) VC(4) both-signs radial nonvanishing; (B) leading-form descent {P_A,Q_B}=0 propagating to a coordinate (THM-1345 s5); (C) Lame-for-polygons CF bound -- and ALL THREE are DESCENT/RETURN-TERMINATION problems = exactly what my coprime-interval/numerical-semigroup/Frobenius engine (S223 DvdK, S224 Wall A) handles; Lame-Fibonacci = the effective bound. codex THM-2045 (smooth factorized R has NO planar Jacobian mate, an exponent-semigroup/Newton-edge obstruction) already uses this engine.

**BOUNDED:** VC(4) GMC-like (E=L o CT) but GMC(2)NOT=>JC(2) (doubling to rank>=2, S205); JC<->LRC 'shared n=12' WITHDRAWN (S137's 12 = Fibonacci proxy). POSITIVE PRIOR: Keller collision min = dim 3, so n=2 is BELOW threshold -> expect JC(2) TRUE + provable by low-rank/coprime-interval means.

**Corrections adopted (mining):** CF thread = mac-mini-S137 (not klein-S329); reify-ladder = deathstar-S75 (not defunct THM-1750); apolarity/Fischer != THM-1685/1710/1735 (those are TNC). NOT a proof of JC(2) -- an assembled route map + verified pieces + the unified descent-termination obstruction. Artifacts: reflection working-jc2-the-obstruction-is-a-descent-termination-...-boxeph-S225.md, HYP-8905, script (+.out).

