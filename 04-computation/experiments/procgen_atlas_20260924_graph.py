#!/usr/bin/env python3
"""
procgen_atlas_20260924_graph.py

Implication atlas: Collatz-family statements versus famous open problems
(session collatz-procgen-20260922, lane 'atlas', 2026-09-24).

This file is the single source of truth for the nodes and edges of the atlas.  It writes
  05-knowledge/results/procgen_atlas_20260924_implications.dot   (Graphviz)
and prints to stdout (section G*) the node table, the edge table, the reachability analysis and
the highlighted chains.

Edge labels (the task's vocabulary):
  PROVED         the implication/equivalence/relation is proved (here or in a repository note)
  CITED          proved in the literature (primary source read [P] or named secondary source [R])
  CONDITIONAL    a proved conditional theorem "open hypothesis X (abc, LW, generalized Mahler, U(1/3)) => Y";
                 the source says whether it is new here or cited
  FALSE          the bridge statement or route is refuted (counterexample or proved comparison)
  OPEN           the implication is expected/heuristic but unproved
  NO KNOWN LINK  nothing is known; the note gives the barrier
Edge kinds: '=>' implication, '<=>' equivalence, '~' structural relation (no implication), 'x' barrier.
"""
import os, sys, itertools
from collections import defaultdict, deque

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DOT = os.path.join(ROOT, "05-knowledge", "results", "procgen_atlas_20260924_implications.dot")

# ---------------------------------------------------------------------------------------------
# NODES: id -> (group, status, short label, statement)
# groups: C = Collatz family (task list), H = helper/intermediate Collatz statement,
#         F = famous problem or theorem, N = new node requested by the coordinator (2026-09-24)
NODES = {
 # --- Collatz family (task list)
 "COL":   ("C", "OPEN", "Collatz conjecture", "every n >= 1 reaches 1 under T(x)=x/2, (3x+1)/2; verified for n < 2^71 (Barina 2025)"),
 "T1":    ("C", "OPEN", "T1: no divergence", "no positive integer has an unbounded T-orbit; = PC on Z^+ (Bernstein 1994: Z^+ meets Phi(non-EP words) nowhere)"),
 "NC":    ("C", "OPEN", "NC: no nontrivial cycle", "no positive T-cycle other than {1,2}; a Pi^0_1 statement"),
 "PC":    ("C", "OPEN", "PC (Lagarias 1985)", "Periodicity Conjecture Q_inf(Q_2)=Q_2: every rational in Z_2 has an eventually periodic parity vector"),
 "C2":    ("C", "OPEN", "C2 (HYP-9123)", "no rational has a non-EP parity vector of bounded discrepancy around a supercritical slope"),
 "Y3":    ("C", "OPEN", "cube-swap (HYP-9127)", "sum_k (2^10/3^9)^(k^3) is irrational in Q_2"),
 "Q1":    ("C", "OPEN", "E-SCC Q1", "every n reaches 1 in the relaxed graph E (Le-Smith loosened graph)"),
 "Q2":    ("C", "OPEN", "E-SCC Q2", "1 reaches every m prime to 3 in E; verified below 10^18"),
 "XMIN":  ("C", "OPEN", "X_min (HYP-9124)", "no integer m >= 2 lies in the backward exceptional set Bad_inf"),
 "CST":   ("C", "OPEN", "Terras sigma = tau", "coefficient stopping time conjecture: sigma(n) = tau(n) for n >= 2; verified to 2.8*10^19 (Rozier-Terracol)"),
 # --- helpers
 "TAU":   ("H", "OPEN", "tau finite on Z^+", "every n >= 1 has 3^(a_j(n)) < 2^j for some j (finite coefficient stopping time)"),
 "PC3K":  ("H", "OPEN", "no divergence for all 3x+k", "no map (3x+k)/2, x/2 with k = +-1 mod 6 has a divergent integer orbit"),
 "SP":    ("H", "OPEN", "swap principle SP(2^L/3^a)", "for 2^L < 3^a, 1 <= a < L, non-EP S: sum_{s in S} (2^L/3^a)^s is irrational in Q_2 (PROVED for 3^a < 2^L; FALSE at a = L, base 2/3)"),
 "T1A":   ("H", "OPEN", "T1 on block classes A^N", "no positive integer has a T-parity vector eventually in A^N, A = >=2 blocks of length L and weight a, 3^a > 2^L (supercritical, positive entropy)"),
 "T1AS":  ("H", "OPEN", "T1 on an adjacent pair class", "the case A = {B,B'}, R_B' = R_B + 1, 3^a > 2^(L+1): e.g. L=10, a=7, A = {0111101110, 1101100111}"),
 "GMAS":  ("F", "OPEN", "Z(3^a/2^L, [R,R+1]/(3^a-2^L))", "e.g. no xi > 0 with xi (2187/1024)^j in Z + [4726,4727]/1163 for all j; an interval 1.88 times the FLP length 1/p"),
 "H9126": ("H", "OPEN", "HYP-9126 (3/2 wall)", "in-house hypothesis; implies X_min"),
 "LS2":   ("H", "OPEN", "Le-Smith Conj. 2", "every nontrivial E-cycle uses an E-only arrow"),
 "LBH":   ("H", "OPEN", "Rozier's LBH", "Lower Bound Hypothesis: n >= j^(-C) 2^((1-H(q/j)) j) for the first j iterates with q odd terms"),
 "LBHS":  ("H", "OPEN", "LBH slices", "n >= c_eps 2^((1-eps) j) for parity prefixes with one even term (Rozier) or of shape 1^k 0^l (Prop. R2)"),
 "MCYC":  ("H", "PROVED", "m-cycle / length bounds", "no m-cycle for m <= 91 (Hercher 2023); with Barina's 2^71 and Hercher Cor. 29 every nontrivial cycle has >= 137,528,045,312 odd terms"),
 "CYCGR": ("H", "OPEN", "L_min(N) >= N^(1/2-eps)", "least cycle length grows like the square root of the verification bound"),
 "TERRAS":("H", "PROVED", "Terras/Everett density 1", "almost every n (natural density) has finite stopping time"),
 "DDOOB": ("N", "PROVED", "deterministic Doob", "d{n : max_{j<=k} T^j(n) >= lam n} = P_k(lam) <= 1/lam, with a Pareto(1) tail (Lundberg exponent 1)"),
 "PAIR":  ("N", "PROVED", "pair-sum fairness", "T(2i-1)+T(2i) = (2i-1)+2i; (3,1) and (3,-1) are the only (q,r) preserving a consecutive pairing"),
 "PERM8": ("N", "OPEN", "orbit of 8 (Collatz 1932)", "under g(3m)=2m, g(3m+1)=4m+1, g(3m+2)=4m+3 the orbit of 8 is infinite"),
 "PERMF": ("N", "OPEN", "g has finitely many cycles", "Lagarias's (C4) for Collatz's permutation; 4 positive cycles known (min 1, 2, 4, 44); finitely many of each length (Keller) and of each m (Simons 2022)"),
 "KOHLC": ("H", "OPEN", "Kohl: T contracting", "every integer trajectory enters the 18-point centre {-136,...,2} of the five known cycles (Kohl thesis 3.6)"),
 "GTRANS":("H", "OPEN", "G_C transitive", "<tau_1(2),4(6), tau_1(3),2(6), tau_2(3),4(6)> acts transitively on N minus 0(6) (Kohl 2017)"),
 "GATE":  ("N", "PROVED", "shared gate 2^K vs 3^L", "T-cycles: n(2^K-3^L)=B(w)>0; g-cycles: n(3^k-2^K)=D_w signed; both on the convergent clocks of log_2 3"),
 "DIV5":  ("H", "OPEN", "5x+1 diverges somewhere", "control: a divergent 5x+1 orbit exists (expected)"),
 "U13":   ("H", "OPEN", "U(1/3) S-unit gap", "in-house uniform S-unit gap hypothesis (cube-theta note)"),
 "GCOLL": ("F", "PROVED", "generalized-Collatz encodings", "each Pi^0_1 / Pi^0_2 statement is equivalent to a statement about an explicit generalized Collatz map, a different map for each problem (not 3x+1)"),
 # --- famous problems and theorems
 "ABC":   ("F", "OPEN", "abc conjecture", "Oesterle-Masser; explicit forms (A. Baker 2004) and the n-term conjecture"),
 "PILLAI":("F", "OPEN", "Pillai's conjecture", "|x^p - y^q| grows; fixed-c finiteness; for bases 2,3 proved (Stroeker-Tijdeman, Bennett)"),
 "BAKER": ("F", "PROVED", "Baker / effective linear forms", "Baker's theorem; explicit two-log bounds (Laurent-Mignotte-Nesterenko, Rhin); Ellison 1971"),
 "LW":    ("F", "OPEN", "Lang-Waldschmidt", "sharp lower bound for linear forms in logarithms; gives irrationality exponent 2 for log_2 3"),
 "PSCH":  ("F", "OPEN", "p-adic Schanuel", "p-adic analogue of Schanuel's conjecture"),
 "ERDOS": ("F", "OPEN", "Erdos ternary problem", "2^n has a ternary digit 2 for every n > 8"),
 "MAHLER":("F", "OPEN", "Mahler 3/2 Z-numbers", "no xi > 0 with 0 <= {xi (3/2)^n} < 1/2 for all n"),
 "GMAH":  ("F", "OPEN", "generalized Z(p/q, I)", "no xi > 0 with xi (p/q)^n in Z + I for all n (Mahler/FLP/Bugeaud/AFS family)"),
 "FLP":   ("F", "PROVED", "FLP / Dubickas / Bugeaud", "no xi > 0 keeps {xi (p/q)^n} in an interval of length < 1/p (FLP), nor in a closed [s,s+1/p] when q<p<q^2 (Dubickas 2009)"),
 "FURST": ("F", "OPEN", "Furstenberg x2x3", "x2,x3-invariant ergodic measures are Lebesgue or finitely supported"),
 "NORMAL":("F", "OPEN", "normality of log_2 3, sqrt 2", "digit normality"),
 "LITTLE":("F", "OPEN", "(mixed) Littlewood", "liminf n ||n a|| ||n b|| = 0; inf n |n|_p ||n a|| = 0"),
 "TWIN":  ("F", "OPEN", "twin primes", "infinitely many p with p+2 prime (Pi^0_2)"),
 "GOLD":  ("F", "OPEN", "Goldbach", "every even n > 2 is a sum of two primes (Pi^0_1)"),
 "LEG":   ("F", "OPEN", "Legendre", "a prime between n^2 and (n+1)^2 (Pi^0_1)"),
 "RH":    ("F", "OPEN", "RH / GRH", "Riemann hypothesis (Pi^0_1 via an elementary criterion)"),
 "LEHMER":("F", "OPEN", "Lehmer's problem", "Mahler measures of non-cyclotomic integer polynomials are bounded away from 1"),
 "GRACE": ("F", "OPEN", "graceful tree", "every tree is graceful (Pi^0_1)"),
 "LRC":   ("F", "OPEN", "lonely runner", "LRC (Pi^0_1 for each n, finitely checkable per n by Tao 2018)"),
 "KL":    ("F", "PROVED", "Kontorovich-Lagarias models", "theorems about stochastic models; as predictions for T they are heuristic"),
 "UNDEC": ("F", "PROVED", "Conway / Kurtz-Simon", "generalized Collatz problems are undecidable (Conway 1972), Pi^0_2-complete (Kurtz-Simon 2007)"),
 "CATAL": ("F", "PROVED", "Catalan-Mihailescu", "3^2 - 2^3 = 1 is the only consecutive perfect powers; for 2,3 already Levi ben Gershon"),
 "RCWA":  ("N", "PROVED", "Kohl RCWA theory", "T is wild (Mod(T^k) = 2^k); surjective non-injective => wild; unbalanced permutations are wild, so Collatz's permutation is wild; tame orbits are rigid"),
 "DOOB":  ("N", "PROVED", "Doob martingale theorems", "convergence and maximal inequality for nonnegative martingales"),
 "ERG2":  ("N", "PROVED", "2-adic ergodic theory of T", "T on Z_2 measure-preserving, strongly mixing, conjugate to the shift (Lagarias 1985 Thm K, L; Bernstein-Lagarias 1996; Akin 2004)"),
}

# ---------------------------------------------------------------------------------------------
# EDGES: (src, dst, kind, label, short text, source)
EDGES = [
 # Collatz-internal skeleton
 ("COL", "NC", "=>", "PROVED", "a cycle other than {1,2} contains no 1", "elementary"),
 ("COL", "T1", "=>", "PROVED", "orbits reaching 1 are bounded", "elementary"),
 ("NC+T1", "COL", "=>", "PROVED", "a bounded orbit is eventually periodic; the only positive cycle is {1,2}", "elementary"),
 ("PC", "PC3K", "<=>", "CITED", "Periodicity Conjecture iff no 3x+k (k = +-1 mod 6) has a divergent integer orbit", "Bernstein-Lagarias 1996 [P], via Lagarias 1990 Cor. 2.1b [R]"),
 ("PC3K", "T1", "=>", "PROVED", "the instance k = 1 on the positive integers", "Lagarias 1985 sec. 2.8 [P]"),
 ("PC", "C2", "=>", "PROVED", "C2 is the bounded-discrepancy supercritical slice of PC", "hypothesis sweep sec. 6"),
 ("C2", "Y3", "=>", "PROVED", "an eventually-Y3 vector has width 1 at slope 9/10", "hypothesis sweep sec. 6; hard-class note"),
 ("PC", "SP", "=>", "PROVED", "block identity: Phi_T(Y_S) = -(1/3^a)[R_B/(1-rho) + (R_B'-R_B) Theta_S(rho)]", "this note, Prop. SP"),
 ("SP", "T1A", "=>", "PROVED", "for two-block alphabets A = {B,B'}", "this note, Prop. SP"),
 ("C2", "T1A", "=>", "PROVED", "A^N words have discrepancy <= L around the slope a/L", "this note"),
 ("T1", "T1A", "=>", "PROVED", "restriction", "trivial"),
 ("GMAH", "T1A", "=>", "CONDITIONAL", "Theorem M: xi = x + (1/3^a) sum R_{c_i} rho^i satisfies xi (3^a/2^L)^j in Z + I_A", "this note, Theorem M"),
 ("FLP", "T1A", "x", "FALSE", "as a route: |I_A| >= 1/(3^a-2^L) > 1/3^a, so the FLP spread bound never applies", "this note, Lemma M2"),
 ("T1A", "T1AS", "=>", "PROVED", "specialisation to adjacent pairs", "trivial"),
 ("GMAS", "T1AS", "<=>", "PROVED", "Theorem M': for R_B' = R_B + 1 and 3^a > 2^(L+1) the confinement xi alpha^j in Z + I forces the digits R_B, R_B+1, i.e. the Collatz blocks", "this note, Theorem M'"),
 ("FLP", "GMAS", "x", "FALSE", "as a route: the needed exclusion length is (alpha/(alpha-1))/p = 1.88/p, 1.59/p, 1.49/p, 1.28/p at L = 10, 16, 19, 20; proved exclusions stop at 1/p", "this note B8; FLP 1995 [R]; Dubickas 2009 [R]; Bugeaud 2004 [P]"),
 ("MAHLER", "T1", "x", "NO KNOWN LINK", "different map (all steps x3/2); SP(2/3) is FALSE; two-place orthogonality (2-adic free vs real free)", "this note; THM-2228; hard-class note"),
 ("MAHLER", "GMAH", "~", "PROVED", "Mahler's problem is the instance (3/2, [0,1/2)); Theorem M needs 3^a/2^L with a < L", "this note"),
 ("ABC", "MAHLER", "~", "CONDITIONAL", "abc gives only radical floors on carry steps; no Z-number obstruction", "THM-3848 sec. 7"),
 ("COL", "Q1", "=>", "PROVED", "an orbit to 1 is a path in E", "synthesis sec. 1"),
 ("XMIN", "Q2", "=>", "PROVED", "endgame Theorem 6.1", "Q2 endgame note"),
 ("H9126", "XMIN", "=>", "PROVED", "Proposition C1", "hypothesis sweep sec. 4"),
 ("CST", "NC", "=>", "CITED", "the minimum m of a nontrivial cycle has tau(m) <= K < sigma(m) = infinity", "Lagarias 1985 sec. 2.3 [P]; re-proved here"),
 ("T1", "TAU", "=>", "PROVED", "tau(n) = infinity forces T^j(n) >= n for all j and no periodic tail", "this note, Prop. CST"),
 ("COL", "TAU", "=>", "PROVED", "the tail (1,2)^inf has coefficients (3/4)^(j/2) -> 0", "this note"),
 ("CST+TAU", "COL", "=>", "PROVED", "sigma(n) = tau(n) < infinity for every n >= 2, then strong induction", "this note, Prop. CST"),
 ("NC", "LS2", "<=>", "PROVED", "Le-Smith Conj. 2 is equivalent to no positive Collatz cycle", "synthesis sec. 1"),
 ("LBH", "COL", "=>", "CITED", "LBH implies the Collatz conjecture (least term of a trajectory avoiding 1 forces q/j >= 0.575)", "Rozier 2017 FACM, Lemma 3.1 [P]"),
 ("LBH", "LBHS", "=>", "PROVED", "specialisation (the slices are weaker than LBH at the same j)", "trivial"),
 ("ABC", "LBHS", "=>", "CONDITIONAL", "abc => n >= K(eps) 2^((1-eps) j) for prefixes with exactly one even term; two-block 1^k 0^l version", "Rozier 2025 Thm 2.1 [P]; Prop. R2 here (PROVED)"),
 ("ABC", "NC", "x", "NO KNOWN LINK", "gate route: abc gives |2^K-3^L| >= 2^(K(1-eps))/C, exponentially weaker than proved Baker bounds; unbounded number of carry terms", "this note C4 (PROVED comparison); Rozier 2025 sec. 6 [P]"),
 ("ABC", "Y3", "x", "NO KNOWN LINK", "the n-term abc conjecture is too weak (radical exponent 2n-5 grows)", "cube-theta note sec. 4"),
 ("U13", "Y3", "=>", "CONDITIONAL", "Theorem U", "cube-theta note sec. 4"),
 ("BAKER", "MCYC", "=>", "CITED", "linear forms + continued fraction + verification", "Simons-de Weger 2005; Hercher 2023 [P]"),
 ("MCYC", "NC", "~", "PROVED", "proved partial results toward NC (bounded complexity only)", "Hercher 2023; recomputed here"),
 ("LW", "CYCGR", "=>", "CONDITIONAL", "irrationality exponent 2 of log_2 3 => L > (3 c N ln2)^(1/(2+eps))", "this note, Prop. LW"),
 ("BAKER", "CYCGR", "=>", "PROVED", "weaker form L >> N^(1/mu_eff) with the effective irrationality exponent", "this note"),
 ("CYCGR", "NC", "~", "PROVED", "quantitative NC relative to verification; never NC itself", "this note"),
 ("PILLAI", "NC", "x", "NO KNOWN LINK", "fixed-c finiteness fixes the clock of a gate value but not the carry divisibility D | B(w)", "gates note sec. 7; this note C6"),
 ("CATAL", "GATE", "=>", "PROVED", "unit gaps (1,1),(2,1),(3,2) classify the cycles of T+-, and the g-cycles on 2/1, 3/2", "catalan lane; this note P6"),
 ("GATE", "NC", "~", "PROVED", "same gate and product identity as a 3x+1 cycle", "this note P6"),
 ("GATE", "PERMF", "~", "PROVED", "g-cycles sit on 2/1, 3/2, 8/5, 19/12; Eliahou transfer: a g-cycle with elements > 10^6 has length >= 665", "this note P6"),
 ("PERM8", "PERMF", "~", "NO KNOWN LINK", "different questions (infinite orbit vs number of cycles)", "-"),
 ("PERM8", "DIV5", "~", "NO KNOWN LINK", "same type: existence of an infinite orbit against positive drift (reverse transfer)", "this note P2, P5"),
 ("RCWA", "PERM8", "x", "NO KNOWN LINK", "g is wild (unbalanced), but wildness does not force an infinite cycle: Kohl's thesis B.4 gives a wild permutation with only finite cycles; no wildness => infinite-cycle conjecture exists", "Kohl thesis 2.5.12, B.4 [P]; this note sec. 8"),
 ("RCWA", "T1", "x", "NO KNOWN LINK", "wildness is DRIFT- and SHEET-blind: every surjective non-injective rcwa map (3x-1, 5x+1) is wild", "Kohl AAM 2007 Thm 2.9 [P]; this note sec. 8"),
 ("KOHLC", "COL", "=>", "CITED", "a contraction centre meeting every trajectory gives the 3n+1 conjecture (and more: Collatz on Z)", "Kohl thesis 3.6 [P]"),
 ("GTRANS", "COL", "<=>", "CITED", "Collatz holds iff G_C acts transitively on N minus 0(6)", "Kohl, J. Group Theory 20 (2017) Prop. 1.2 [P]"),
 ("RCWA", "GTRANS", "~", "CITED", "RCWA reformulations of 3x+1 (class-transposition groups; sigma_T on Z^2)", "Kohl 2017; thesis 3.13 [P]"),
 ("BAKER", "PERMF", "~", "CITED", "finitely many m-cycles of Collatz-type permutations; none for m <= 2 besides the known ones", "Simons 2022 [P]; Keller [P via Kohl]"),
 ("ERDOS", "Q2", "x", "NO KNOWN LINK", "Q2 needs low 3-adic digits along dynamically chosen K; every 3-adic-window Erdos statement is FALSE", "transversality foundry sec. 3"),
 ("ERDOS", "T1", "x", "NO KNOWN LINK", "different map and place", "transversality foundry sec. 5"),
 ("FURST", "T1", "x", "NO KNOWN LINK", "a measure statement: DEFECT-blind for every every-orbit target", "transversality foundry sec. 5"),
 ("NORMAL", "NC", "x", "NO KNOWN LINK", "cycle bounds use the continued fraction / irrationality exponent of log_2 3, not its digits", "this note"),
 ("LITTLE", "Q2", "x", "NO KNOWN LINK", "no statement of either form controls the carry/landing conditions", "this note"),
 ("PSCH", "PC", "x", "NO KNOWN LINK", "Bernstein numbers are not exponential-type values", "cube-theta note sec. 0"),
 ("TWIN", "GCOLL", "~", "CITED", "Pi^0_2 statement; Kurtz-Simon reduction to an 'all orbits reach 1' statement of a generalized Collatz map", "Kurtz-Simon 2007 [P]"),
 ("GOLD", "GCOLL", "~", "CITED", "Pi^0_1: non-halting of an explicit program; Conway's FRACTRAN encoding", "Conway 1972/1987 [R]"),
 ("LEG", "GCOLL", "~", "CITED", "Pi^0_1, as for Goldbach", "Conway [R]"),
 ("RH", "GCOLL", "~", "CITED", "Pi^0_1 via Lagarias's elementary criterion (sigma(n) <= H_n + exp(H_n) log H_n)", "Lagarias 2002 Thm 1.1 [P]; Conway [R]"),
 ("GRACE", "GCOLL", "~", "CITED", "Pi^0_1 (finite search per tree)", "Conway [R]"),
 ("LRC", "GCOLL", "~", "CITED", "Pi^0_1 (Tao 2018: finite check for each n)", "Tao 2018 Thm 1.3, Cor. 1.4 [P]; Conway [R]"),
 ("UNDEC", "GCOLL", "~", "CITED", "universality of generalized Collatz maps", "Conway 1972; Kurtz-Simon 2007"),
 ("GCOLL", "COL", "x", "NO KNOWN LINK", "the encodings use other maps; the 3x+1 map itself is not known to be universal or undecidable", "Kurtz-Simon; GGM 2025 remark"),
 ("RH", "COL", "x", "NO KNOWN LINK", "trunk/Iwasawa coordinate only; zeta_3 zero-free; no bridge", "trunk/RH note"),
 ("LEHMER", "COL", "x", "NO KNOWN LINK", "no shared object", "this note"),
 ("LRC", "COL", "x", "NO KNOWN LINK", "analogy only (S596 two-block question)", "bridges note"),
 ("KL", "COL", "=>", "OPEN", "stochastic models predict Collatz (and 5x+1 divergence); transfer to the map is unproved", "Kontorovich-Lagarias 2010 [P]"),
 ("KL", "DIV5", "=>", "OPEN", "the 5x+1 model has positive drift", "Kontorovich-Lagarias 2010 [P]"),
 ("PAIR", "KL", "~", "PROVED", "AM factor (q+1)/4 = 1 iff q = 3; the model drift log(sqrt3/2) is the AM-GM gap", "this note F2; bridges note"),
 ("ERG2", "DOOB", "~", "CITED", "Haar-independence of parity bits makes M_j = 3^a_j/2^j a martingale", "Lagarias 1985 Thm K [P]"),
 ("DOOB", "TERRAS", "=>", "PROVED", "M_j -> 0 a.s., so tau < infinity a.s., so the stopping density is 1", "this note F4 (Terras 1976 [R] first)"),
 ("DOOB", "DDOOB", "=>", "PROVED", "Doob's maximal inequality on the exact finite martingale", "this note F3, F5"),
 ("PAIR", "DDOOB", "=>", "PROVED", "fairness = martingale = Lundberg exponent 1 (Pareto(1) excursion tail)", "this note F2, F3"),
 ("TERRAS", "COL", "x", "NO KNOWN LINK", "density statements are DEFECT-blind", "barrier atlas"),
 ("DDOOB", "T1", "x", "NO KNOWN LINK", "density statements are DEFECT-blind", "barrier atlas"),
 ("UNDEC", "COL", "x", "NO KNOWN LINK", "UNIFORM barrier binds uniform methods only", "barrier atlas sec. 3"),
 ("CATAL", "PILLAI", "~", "CITED", "Catalan is the case c = 1 of Pillai's Conj. 1.3", "Mihailescu 2004 [R]; Waldschmidt 2004 Thm 1.2 [P]"),
 ("ABC", "PILLAI", "=>", "CONDITIONAL", "abc implies Pillai's Conj. 2.6: |x^p - y^q| >= C max^(1-1/p-1/q-eps)", "Waldschmidt 2004, p. 13 [P]"),
 ("LW", "PILLAI", "=>", "CONDITIONAL", "Lang-Waldschmidt (Conj. 2.5) implies Conj. 2.6; for 2, 3: |2^a - 3^b| >= 2^a a^(-1-eps)", "Waldschmidt 2004, p. 13 [P]"),
 ("BAKER", "PILLAI", "=>", "CITED", "for bases 2, 3: |2^x - 3^y| > 2^x e^(-x/10) (x outside a finite list); 3^x - 2^y = c has at most one solution for |c| > 13; at most two for general a, b", "Ellison 1971 Thm 3 [P]; Stroeker-Tijdeman 1982 [R]; Bennett 2001 Thm 1.1 [P]"),
 ("SP", "MAHLER", "x", "FALSE", "the unrestricted swap principle at Mahler's base 2/3 is FALSE (M-orbit of 1)", "this note B5"),
 ("PC", "Q2", "x", "NO KNOWN LINK", "different quantifier structure (backward reachability, low 3-adic digits)", "foundry v4"),
]

HIGHLIGHT = {
 ("GMAS", "T1AS"), ("T1A", "T1AS"), ("GTRANS", "COL"),
 ("PC", "PC3K"), ("PC3K", "T1"), ("PC", "C2"), ("C2", "Y3"), ("PC", "SP"), ("SP", "T1A"),
 ("GMAH", "T1A"), ("T1", "T1A"),
 ("LW", "CYCGR"), ("ABC", "LBHS"), ("LBH", "COL"), ("CST", "NC"), ("CST+TAU", "COL"), ("T1", "TAU"),
 ("BAKER", "MCYC"), ("CATAL", "GATE"), ("GATE", "PERMF"), ("DOOB", "TERRAS"), ("PAIR", "DDOOB"),
 ("NC+T1", "COL"), ("XMIN", "Q2"), ("H9126", "XMIN"), ("COL", "Q1"),
}

NEW_HERE = {("GMAS", "T1AS"), ("FLP", "GMAS"), ("T1A", "T1AS"), ("PC", "SP"), ("SP", "T1A"), ("C2", "T1A"), ("GMAH", "T1A"), ("FLP", "T1A"), ("T1", "TAU"),
            ("CST+TAU", "COL"), ("LW", "CYCGR"), ("BAKER", "CYCGR"), ("CATAL", "GATE"), ("GATE", "NC"),
            ("GATE", "PERMF"), ("DOOB", "TERRAS"), ("DOOB", "DDOOB"), ("PAIR", "DDOOB"), ("SP", "MAHLER"),
            ("MAHLER", "GMAH"), ("COL", "TAU")}

def check():
    ids = set(NODES)
    for e in EDGES:
        for end in (e[0], e[1]):
            for part in end.split("+"):
                assert part in ids, (part, e)
        assert e[2] in ("=>", "<=>", "~", "x"), e
        assert e[3] in ("PROVED", "CITED", "CONDITIONAL", "FALSE", "OPEN", "NO KNOWN LINK"), e

STATUS_FILL = {"PROVED": "#d9f2d9", "OPEN": "#ffffff", "FALSE": "#f8d0d0"}
GROUP_SHAPE = {"C": "box", "H": "box", "F": "ellipse", "N": "hexagon"}
GROUP_STYLE = {"C": "filled,bold", "H": "filled,rounded,dashed", "F": "filled", "N": "filled"}
GROUP_FILL = {"C": "#dbe9ff", "H": "#eef4ff", "F": "#fff7d6", "N": "#f0e6ff"}
EDGE_STYLE = {
 "PROVED": ('solid', '#000000'), "CITED": ('solid', '#1f4fbf'), "CONDITIONAL": ('dashed', '#d07000'),
 "FALSE": ('dotted', '#c00000'), "OPEN": ('dashed', '#808080'), "NO KNOWN LINK": ('dotted', '#b0b0b0'),
}

def wrap(s, n=26):
    words, lines, cur = s.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > n and cur:
            lines.append(cur); cur = w
        else:
            cur = (cur + " " + w).strip()
    if cur:
        lines.append(cur)
    return "\\n".join(lines)

def write_dot():
    L = []
    L.append("// Implication atlas: Collatz-family statements vs famous problems (procgen_atlas_20260924).")
    L.append("// Generated by 04-computation/experiments/procgen_atlas_20260924_graph.py -- do not edit by hand.")
    L.append("// Render: dot -Tsvg procgen_atlas_20260924_implications.dot -o atlas.svg")
    L.append("digraph atlas {")
    L.append('  graph [rankdir=LR, fontname="Helvetica", fontsize=11, nodesep=0.25, ranksep=0.9, splines=true, overlap=false,')
    L.append('         label="Collatz implication atlas (2026-09-24). Boxes: Collatz family (bold) and helpers (dashed); ellipses: famous problems/theorems; hexagons: nodes added 2026-09-24.\\nEdges: black PROVED, blue CITED, orange dashed CONDITIONAL, red dotted FALSE (route refuted), grey dashed OPEN (heuristic), light dotted NO KNOWN LINK (barrier). Thick green: highlighted chains. Diamond: AND.", labelloc=b];')
    L.append('  node [fontname="Helvetica", fontsize=10];')
    L.append('  edge [fontname="Helvetica", fontsize=8];')
    for nid, (grp, st, short, stmt) in NODES.items():
        fill = GROUP_FILL[grp] if st == "OPEN" else STATUS_FILL.get(st, "#ffffff")
        lab = wrap(short) + "\\n[" + st + "]"
        L.append(f'  {nid} [label="{lab}", shape={GROUP_SHAPE[grp]}, style="{GROUP_STYLE[grp]}", fillcolor="{fill}", tooltip="{stmt}"];')
    # AND junctions
    ands = sorted({e[0] for e in EDGES if "+" in e[0]})
    for a in ands:
        aid = a.replace("+", "_AND_")
        L.append(f'  {aid} [label="AND", shape=diamond, fontsize=8, width=0.3, height=0.3, style=filled, fillcolor="#eeeeee"];')
        for part in a.split("+"):
            L.append(f'  {part} -> {aid} [arrowhead=none, color="#000000"];')
    for (s, d, kind, lab, txt, src) in EDGES:
        s_id = s.replace("+", "_AND_")
        style, color = EDGE_STYLE[lab]
        attrs = [f'style={style}', f'color="{color}"', f'tooltip="{lab}: {txt} ({src})"']
        short = {"=>": "", "<=>": "<=>", "~": "~", "x": "x"}[kind]
        elab = (short + " " if short else "") + ("new " if (s, d) in NEW_HERE else "") + lab.lower()
        attrs.append(f'label="{elab}"')
        if kind == "<=>":
            attrs.append("dir=both")
        if kind in ("~", "x"):
            attrs.append("arrowhead=none")
            attrs.append("constraint=false")
        if (s, d) in HIGHLIGHT:
            attrs.append("penwidth=3.2")
            if lab in ("PROVED", "CITED"):
                attrs[1] = 'color="#1a8a1a"'
        L.append(f'  {s_id} -> {d} [{", ".join(attrs)}];')
    L.append('  { rank=same; COL; }')
    L.append("}")
    with open(DOT, "w") as f:
        f.write("\n".join(L) + "\n")
    return DOT

def validate_dot(path):
    """Minimal syntax check (no Graphviz on this machine): balanced braces/brackets outside quotes,
    closed quotes, and every edge endpoint declared as a node."""
    import re
    txt = open(path).read()
    depth_b = depth_s = 0
    inq = False
    prev = ""
    for ch in txt:
        if ch == '"' and prev != "\\":
            inq = not inq
        elif not inq:
            if ch == "{": depth_b += 1
            if ch == "}": depth_b -= 1
            if ch == "[": depth_s += 1
            if ch == "]": depth_s -= 1
            assert depth_b >= 0 and depth_s >= 0
        prev = ch
    assert not inq and depth_b == 0 and depth_s == 0, (inq, depth_b, depth_s)
    declared = set(re.findall(r"^  (\w+) \[label=", txt, flags=re.M))
    for a, b in re.findall(r"^  (\w+) -> (\w+) \[", txt, flags=re.M):
        assert a in declared and b in declared, (a, b)
    return len(declared), len(re.findall(r" -> ", txt))

def reach():
    # directed implication graph over PROVED/CITED edges (+ CONDITIONAL marked)
    adj = defaultdict(list)
    for (s, d, kind, lab, txt, src) in EDGES:
        if kind in ("=>", "<=>") and lab in ("PROVED", "CITED", "CONDITIONAL"):
            if "+" in s:
                continue  # AND-edges handled separately
            adj[s].append((d, lab))
            if kind == "<=>":
                adj[d].append((s, lab))
    res = {}
    for src in NODES:
        seen = {src: []}
        dq = deque([src])
        while dq:
            u = dq.popleft()
            for v, lab in adj[u]:
                if v not in seen:
                    seen[v] = seen[u] + [(u, v, lab)]
                    dq.append(v)
        res[src] = seen
    return res

def main():
    check()
    print()
    print("=" * 78)
    print("G1  nodes")
    print("=" * 78)
    for nid, (grp, st, short, stmt) in NODES.items():
        print(f"  {nid:7s} [{grp}] {st:7s} {short}: {stmt}")
    print()
    print("=" * 78)
    print("G2  edges (kind: => implication, <=> equivalence, ~ relation, x barrier)")
    print("=" * 78)
    for i, (s, d, kind, lab, txt, src) in enumerate(EDGES, 1):
        flag = " *NEW*" if (s, d) in NEW_HERE else ""
        hl = " [HIGHLIGHT]" if (s, d) in HIGHLIGHT else ""
        print(f"  E{i:02d} {s:>8s} {kind:3s} {d:7s} {lab:13s}{flag}{hl}: {txt}  <{src}>")
    counts = defaultdict(int)
    for e in EDGES:
        counts[e[3]] += 1
    print("  label counts:", dict(counts), " total:", len(EDGES))
    print()
    print("=" * 78)
    print("G3  reachability through PROVED/CITED/CONDITIONAL implications (single-premise edges)")
    print("=" * 78)
    R = reach()
    fam = [n for n, v in NODES.items() if v[0] in ("F", "N")]
    col = [n for n, v in NODES.items() if v[0] in ("C", "H")]
    for f in fam:
        tgts = [c for c in R[f] if c in col and c != f]
        if tgts:
            print(f"  {f:7s} => " + ", ".join(sorted(tgts)))
    print("  Collatz-family nodes that imply a famous-problem node: " +
          (", ".join(f"{c}=>{f}" for c in col for f in fam if f in R[c]) or "none"))
    print()
    print("=" * 78)
    print("G4  highlighted chains")
    print("=" * 78)
    chains = [
     ["PC", "PC3K", "T1", "T1A", "T1AS", "GMAS"], ["PC", "C2", "Y3"], ["PC", "SP", "T1A"], ["GMAH", "T1A"],
     ["GTRANS", "COL"], ["KOHLC", "COL"],
     ["LW", "CYCGR"], ["BAKER", "MCYC"], ["ABC", "LBHS"], ["LBH", "COL"],
     ["CST", "NC"], ["T1", "TAU"], ["H9126", "XMIN", "Q2"], ["COL", "Q1"],
     ["CATAL", "GATE"], ["DOOB", "TERRAS"], ["DOOB", "DDOOB"],
    ]
    E = {(e[0], e[1]): e for e in EDGES}
    for ch in chains:
        parts = []
        for a, b in zip(ch, ch[1:]):
            e = E.get((a, b)) or E.get((b, a))
            parts.append(f"{a} ={e[3]}=> {b}" if e[2] in ("=>", "<=>") else f"{a} ~{e[3]}~ {b}")
        print("  " + " ; ".join(parts))
    print("  AND-chains: NC & T1 =PROVED=> COL ; CST & TAU =PROVED=> COL (so CST & T1 => COL).")
    p = write_dot()
    nd, ne = validate_dot(p)
    print(f"\n  wrote {os.path.relpath(p, ROOT)} ({len(NODES)} nodes, {len(EDGES)} edges); "
          f"syntax check PASS ({nd} node statements, {ne} edge statements incl. AND junctions)")

if __name__ == "__main__":
    main()
    print("\nDONE graph")
