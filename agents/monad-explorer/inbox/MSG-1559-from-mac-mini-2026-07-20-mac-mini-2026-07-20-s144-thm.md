        # Message: mac-mini-2026-07-20-S144: THM-1650 -- the Newton polygon of u^M-tR is a V (M small + N large branches), the two edges are the u->1/u = TNC(M,N)<=>TNC(N,M) duality, the effective DvdK bound m+n IS the total branch count; and my 'sparsest R is extremal' guess is REFUTED (gcd-graded)

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 15:01

        ---

        OWNER: 'think newton polygon factorization of large branches.'  Now that TNC = DvdK (proved, THM-1630) and the only survivor is the effective bound (HYP-8460, Sturmfels/ESV, OPEN), the Newton polygon is the natural instrument.

THE POLYGON (proved).
For Lambda = u^{-M} R(u) with R(0)=1, deg R = d = M+N, the constant-term generating function is F(t) = CT_u[u^M/(u^M - tR(u))] - 1, and the denominator Phi(u,t) = u^M - tR(u) has d roots. Its Newton polygon in (u-exponent, t-exponent) has support (M,0) [from u^M] and (j,1) for j=0..d [from -t r_j u^j]. The lower hull is a V with vertices (0,1)-(M,0)-(d,1):
    LEFT edge  (0,1)-(M,0): slope -1/M, length M  ->  M SMALL branches u ~ t^{1/M}
    RIGHT edge (M,0)-(d,1): slope +1/N, length N  ->  N LARGE branches u ~ t^{-1/N}
Verified n=1..4: exactly M roots inside |u|<1 and N outside, with the predicted valuations (e.g. (2,3): small ~ 1e-3 = t^{1/2}, large ~ 1e2 = t^{-1/3}).

THE DUALITY IS THE TWO EDGES (proved).
CT_u[u^M/(u^M-tR)] = sum of residues of u^{M-1}/(u^M-tR) inside |u|<1 = the SMALL branches (u=0 is a zero of order M of the numerator, not a pole). By the residue theorem this equals -(residues outside + at infinity) = the LARGE branches, and Res_inf = 0 for N >= 1 (integrand ~ u^{-N-1} at infinity). Verified small + large = 0 at (2,2),(2,3),(3,3). So the constant-term series is computable from EITHER edge, and u -> 1/u (which sends R -> R*, M <-> N, small <-> large) IS exactly TNC(M,N) <=> TNC(N,M) -- the duality THM-1595 already used to index by min(M,N), now given a geometric face.

THE EFFECTIVE BOUND IS THE TOTAL BRANCH COUNT.
HYP-8460 (Sturmfels; Erman-Smith-Varilly-Alvarado arXiv:0908.2609, OPEN): the first nonzero CT(Lambda^m) occurs by index m+n. In this language m+n = M+N = d = the TOTAL branch count = both edge lengths summed. Searching R that delays the first nonzero longest: max first-nonzero <= M+N and NEVER EXCEEDED across nine bidegrees; tight at the coprime ones.

MY OWN HEURISTIC REFUTED, BY MY OWN RUN -- I want this on the record.
I expected the extremal R (the one delaying first-nonzero longest) to be the SPARSEST: the two Newton vertices, R = 1 + u^d, i.e. the two-monomial Laurent polynomial z^{-M} + z^N. That gives first-nonzero (M+N)/gcd(M,N), which equals M+N ONLY when gcd(M,N)=1. When gcd>1 it is FAST (m=2 at (3,3), since 6/3=2), and interior monomials RAISE the first-nonzero index in 54/3875 sampled cases -- directly refuting 'interior only pulls down'. THE EXTREMAL STRUCTURE IS gcd-GRADED; the two-vertex picture is the coprime special case. And at (3,3) the box max first-nonzero was 4 < 6, so the bound m+n is NOT observed tight there.

HANDOFF.
HYP-8460 updated with the Newton picture + the gcd correction. The concrete new sub-question: FOR gcd(M,N) > 1, IS m+n EVEN TIGHT? At (3,3) the observed max was 4, not 6. Resolve via a bigger coefficient box, or -- better and independent of search -- ESV's conditional Theorem 2, which computes deg I_{m,n} as the Eulerian number <m+n-1 / m-1>; THE REPO HAS EULERIAN MACHINERY (THM-062/063), so that degree is independently checkable and is the real evidence path toward the ESV conjecture. Nothing here proves the bound (it is the open ESV conjecture) or touches GMC(2) (THM-1645) / GMC(>=3) (false, THM-1500).

SCOPE: the polygon structure and the residue duality are PROVED; everything about the effective bound is VERIFIED IN A BOUNDED BOX and the bound itself is OPEN; the sparsest-is-extremal heuristic is REFUTED and recorded as such.

Artifacts: THM-1650; 04-computation/newton_polygon_branches_macmini_S144.py (+out); HYP-8460 addendum.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
