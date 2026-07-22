## boxeph-2026-07-21-S218 -- arithmetic entropy is a repo-wide invariant; the rigid extremum = the zero-entropy point (HYP-8875)

**Owner:** extend the arithmetic-entropy idea (S217) and apply it to as many repo pieces as possible.

**ONE INVARIANT:** H_arith(X|L)=log2|{X': L(X')=L(X)}| = the GLOBAL bits of X HIDDEN from a LOCAL invariant L. Zero=local-determines-global (RIGID); positive=hidden global object. FOUR instances (verified arithmetic_entropy_across_the_repo_boxeph_S218.py):
1. BINARY FORMS|genus (refines S217): the truly-hidden part is the DEEP within-genus class group (genus is congruence-detectable). h=genera x deep; Heegner -3,-7,-11 h=1 zero; -15=-3*5 pure GENUS (visible, 0 deep); -23,-47 pure DEEP (1.58,2.32 hidden bits = Hilbert class field, invisible to congruences).
2. TOURNAMENTS|score sequence: transitive (0..n-1)=UNIQUE realization (Landau) => H=0 rigid = the AP/nullcone/rank-11 vertex (S214); near-regular scores carry the hidden fiber (n=5: (1,2,2,2,3) 3 classes = kps reconstruction wall).
3. REALS|CF prefix: golden [0;1,1,..] geo-mean 1 << Khinchin = worst-approx = LRC FOIL (S206); t*=14/183=[0;13,14] geo-mean 13.5 = well-approx = extremal.
4. NULLCONE|moment depth = certificate entropy: LRC finite (bounded alphabet, Bonferroni depth ~5) vs GMC infinite (unbounded degree, Watson S211).

**DUAL entropies:** score-DISTRIBUTION entropy (transitive MAX spread) vs RECONSTRUCTION entropy (transitive 0). The AP = max-order + zero-hidden-info = rigidity; the regular/Paley = min spread + max hidden info.

**UNIFYING:** every repo RIGID extremum (AP/transitive/Heegner h=1/reify-ladder vertex) = a zero-arithmetic-entropy point (local determines global); every DIFFICULTY = its positive-entropy hidden object (deep class group / cospectral fiber / CF tail / deep moment). Rigidity = why the extremal is unique; hidden entropy = where the proof still must go.

**Honest:** genus/class-group, Landau, Khinchin/Levy, detection-depth facts are classical/verified; the contribution is the UNIFICATION (one info deficit across binary forms / tournaments / CF / nullcones) + the rigid=zero-entropy observation. Organizing lens, gate-independent (survives S217 MISTAKE-225), not a proof step. Artifacts: reflection arithmetic-entropy-is-a-repo-wide-invariant-...-boxeph-S218.md, HYP-8875, script (+.out).

