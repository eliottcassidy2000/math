#!/usr/bin/env python3
"""Catalogue of PROVED 2-versus-3 transversality theorems, typed with the foundry
barriers, and the cross-problem matrix (collatz-procgen-20260922, transversality lane).

Access marks: [P] primary text read by this lane's literature sub-lanes; [R] statement
taken from the named secondary source because the primary was paywalled or bot-blocked
(nothing was bypassed); [repo] a result proved in this repository.  The typing is this
lane's modelling judgment under the rule printed by print_catalogue().
"""

TYPING_RULE = [
    "SHEET    B if the conclusion holds verbatim for the negated / minus-sheet analogue; O if it fails there; - no analogue.",
    "DRIFT    B if the theorem holds for every multiplicatively independent pair (e.g. (2,5), or 5/2, where the Collatz",
    "         drift 3<2^2, resp. Mahler's threshold p<q^2, has the other sign); O if it uses 3<2^2 (p<q^2) essentially.",
    "DEFECT   B for counting/density/measure/dimension statements (a planted density-zero set is invisible);",
    "         O for every-element statements (every n, every orbit, every point beyond a threshold).",
    "INTEGRAL O if the statement isolates the rational/periodic exceptional points; B if it tolerates extra ones; - n/a.",
    "UNIFORM  B if the proof runs unchanged over an explicit infinite family (all pairs, all maps of a class);",
    "         O if per-instance (explicit finite computation / bounded complexity, decidable for each instance).",
    "kind     COUNT (i) = counting/density/dimension;  EVERY (ii) = every element/orbit/point;  FINITE = finite range.",
]

# key: (short name, access, statement, kind, effective, mechanism, typing, serves)
C = {}


def add(key, name, access, statement, kind, eff, mech, SHEET, DRIFT, DEFECT, INTEGRAL, UNIFORM, serves, note=""):
    C[key] = dict(name=name, access=access, statement=statement, kind=kind, eff=eff, mech=mech,
                  typing=dict(SHEET=SHEET, DRIFT=DRIFT, DEFECT=DEFECT, INTEGRAL=INTEGRAL, UNIFORM=UNIFORM),
                  serves=serves, note=note)


# ---------------------------------------------------------------- the items asked for
add("Furstenberg1967", "Furstenberg 1967, Math. Systems Theory 1, 1-49",
    "[R] via Rudolph 1990 intro, Lindenstrauss 2005 Thm 1.1, zbMATH Zbl 0146.28502 (primary paywalled)",
    "p,q >= 2 multiplicatively independent, X subset R/Z closed with pX, qX subset X => X finite (rationals) or X = R/Z;"
    " equivalently {p^m q^n alpha} is dense mod 1 for every irrational alpha. The measure form (atomless ergodic"
    " x2,x3-invariant mu is Lebesgue) is Furstenberg's conjecture: OPEN.",
    "EVERY", "no (BLMV 2009 effective)", "topological dynamics of non-lacunary semigroups",
    "B", "B", "O", "O", "B", "none directly (orbit of a 2-parameter semigroup, not of {2^n} or of Z)")
add("RudolphJohnson", "Rudolph 1990 (ETDS 10, 395-406); Johnson 1992 (Israel J. Math. 77)",
    "Rudolph [P] Thm 4.9; Johnson [R] via Lindenstrauss survey",
    "an ergodic x p, x q-invariant probability measure (p,q coprime: Rudolph; multiplicatively independent: Johnson)"
    " with positive entropy is Lebesgue measure",
    "COUNT", "-", "entropy, symbolic cover, Z^2 inverse limit", "B", "B", "B", "B", "B", "none")
add("Host1995", "Host 1995 (Israel J. Math. 91); Hochman-Shmerkin 2015 (Invent. Math. 202)",
    "Host [R] via HS2015 Thm 1.9; HS2015 [P] arXiv 1302.5792",
    "mu x n-invariant ergodic with positive entropy, gcd(m,n)=1 => mu-a.e. point is normal in base m (Host);"
    " HS2015 Thm 1.10 extends to multiplicatively independent bases",
    "COUNT", "-", "entropy / local entropy averages", "B", "B", "B", "B", "B", "none")
add("CasselsSchmidt", "Cassels 1959 (Colloq. Math. 7); Schmidt 1960 (Pacific J. Math. 10)",
    "Cassels [R] via Khoshnevisan, HS2015; Schmidt [P] Thm 1-2",
    "Cantor-Lebesgue-a.e. point of the middle-third Cantor set is normal in every base that is not a power of 3;"
    " Schmidt: r, s multiplicatively independent => continuum many numbers normal to r and not simply normal to s",
    "COUNT", "-", "Fourier / variance estimates", "B", "B", "B", "B", "B", "none")
add("HochmanShmerkin2012", "Hochman-Shmerkin 2012, Ann. of Math. 175, 1001-1059",
    "[P] arXiv 0910.1956v2 Thm 1.3, Sec. 10.4",
    "mu, nu invariant under T_m, T_n (m,n not powers of one integer) => dim pi(mu x nu) = min(1, dim(mu x nu)) for"
    " every non-principal linear projection; sets: dim(X + sY) = min(1, dim X + dim Y) for closed x2-, x3-invariant"
    " X, Y and s != 0",
    "COUNT", "-", "local entropy averages", "B", "B", "B", "B", "B", "none")
add("ShmerkinWu2019", "Shmerkin 2019 (Ann. of Math. 189, 319-391); Wu 2019 (Ann. of Math. 189, 707-751)",
    "[P] arXiv 1609.07802v3 Thm 1.2; arXiv 1609.08053v3 Thm 1.4",
    "A, B subset [0,1] closed, x p- resp. x q-invariant, p,q multiplicatively independent: for every invertible affine"
    " g (Shmerkin) / all real u,v (Wu), upper box dim (A cap g(B)) <= max(0, dim_H A + dim_H B - 1)"
    " (Furstenberg's intersection conjecture)",
    "COUNT", "-", "L^q norms of convolutions (Shmerkin); CP-chains / ergodic theory (Wu)",
    "B", "B", "B", "B", "B", "none (a dimension-0 intersection can be nonempty: never excludes a point)")
add("Tao2022", "Tao 2022, Forum Math. Pi 10, e12",
    "[P] arXiv 1909.03562v7 and published version",
    "Prop 1.17: for n >= 1 and xi in Z/3^n not divisible by 3, |E e(-2 pi i xi Syrac(Z/3^n Z)/3^n)| <<_A n^(-A)"
    " uniformly, Syrac(Z/3^n Z) = sum_{i=1..n} 3^(i-1) 2^(-(a_1+..+a_i)) mod 3^n (eq. 1.29), a_i iid Geom(2); equivalent to"
    " fine-scale mixing (Prop 1.14). Thm 1.3: Col_min(N) < f(N) for log-almost all N, any f -> infinity",
    "COUNT", "no (A arbitrary, constants not explicit)", "2-D renewal process, 'triangles' of 2^a 3^b",
    "B", "B for Prop 1.17 (the drift 3<2^2 enters only through first-passage stabilisation; GGM's q<p^(p/(p-1)) for Thm 1.3)",
    "B", "B", "B", "none (log density; DEFECT-blind by construction)")
add("SengeStraus1973", "Senge-Straus 1973, Period. Math. Hungar. 3, 93-100",
    "[R] via Stewart 1980 intro, Spiegelhofer 2023, Dimitrov-Howe (cite 'Thm 3, p. 100'); Springer bot-check",
    "a, b >= 2: {n : s_a(n) <= c and s_b(n) <= c} is finite for every c iff log a/log b is irrational"
    " (equivalently with numbers of nonzero digits); e.g. 2^n has unboundedly many nonzero ternary digits",
    "EVERY", "no (Ridout-type p-adic Roth)", "Mahler's p-adic Thue-Siegel-Roth", "B", "B", "O", "-", "B",
    "T3/T1 only in the zero-dimensional form (few nonzero digits)")
add("Stewart1980", "Stewart 1980, J. reine angew. Math. 319, 63-72",
    "[P] GDZ scan, Thm 1, Thm 2",
    "Thm 1: log a/log b irrational => L_(alpha,a)(n) + L_(beta,b)(n) > loglog n/(logloglog n + C) - 1 for n > 25,"
    " C effective in a,b (L counts digits different from alpha, beta). Thm 2: 3^n has more than"
    " log n/(loglog n + C0) - 1 nonzero binary digits for n > 4 (same for 2^n in base 3)",
    "EVERY", "yes", "Baker's linear forms + Loxton-van der Poorten + digit-gap pigeonhole",
    "B", "B", "O", "-", "B", "T3 only in the zero-dimensional form")
add("Mahler1957", "Mahler 1957, Mathematika 4, 122-124 (with Corvaja-Zannier 2004, Zudilin 2007)",
    "[R] via Corvaja-Zannier (arXiv math/0403522 [P]), Bilu, Zudilin JTNB 19 (2007) [P]",
    "alpha > 1 rational non-integer, 0 < theta < 1 => ||alpha^n|| < theta^n for only finitely many n; so"
    " ||(3/2)^n|| > (3/4)^n for n >= n0 (ineffective, Ridout). Effective: ||(3/2)^k|| > 0.5803^k for k >= K"
    " (Zudilin 2007, K effective). Corvaja-Zannier Thm 1: ||alpha^n|| < l^n infinitely often => a power of"
    " alpha is Pisot",
    "EVERY", "no (Mahler); yes with 0.5803 (Zudilin)", "Ridout / Subspace theorem; Pade (effective)",
    "B", "B", "O", "O (Pisot/integer alpha are exactly the exceptions)", "B", "T4 only for xi = 1 (not Z-numbers)")
add("Mahler1968", "Mahler 1968, J. Austral. Math. Soc. 8, 313-321",
    "[P] reprint in Mahler Selecta (Documenta Math. 2019); Flatto 1992 [R] via Lagarias bibliography",
    "each [g, g+1) contains at most one Z-number, and it lies in [g, g+1/2) (so Z-numbers are countable); at most"
    " x^0.7 Z-numbers in [0,x] for large x (Flatto 1992: exponent log_2(3/2))",
    "COUNT", "yes", "the same series read in R and in Z_2; Fibonacci count of admissible classes (phi < 2^0.7)",
    "-", "O (for p/q = 5/2 > q^2 Z_{p/q}-numbers exist: Tijdeman 1972, Flatto 1992 [R via Andrieu-Eliahou-Vivion])",
    "B", "-", "O", "T4 (countability; partial)")
add("FLP1995", "Flatto-Lagarias-Pollington 1995, Acta Arith. 70, 125-147",
    "[R] via Lu-Zheng arXiv 2603.16794 [P], Farhi arXiv math/0611622 [P] (Acta Arith. PDF bot-blocked)",
    "Thm 1.4: for every xi > 0 and coprime p > q > 1, limsup {xi (p/q)^n} - liminf {xi (p/q)^n} >= 1/p"
    " (spread >= 1/3 for 3/2). Negative ratio -3/2: spread >= 11/27 (Lu-Zheng 2026)",
    "EVERY", "yes", "symbolic dynamics of the carry words", "B (both signs)", "B (uniform in p/q)", "O", "-", "B",
    "T4 partial (needs spread > 1/2)")
add("DubickasMossinghoff2009", "Dubickas-Mossinghoff 2009, Math. Comp. 78, 1837-1851",
    "[R] via Andrieu-Eliahou-Vivion arXiv 2510.11723v2 [P] (AMS bot-check)",
    "Thm 5.1: no Z_{p/q}-number below 2^57, 10^32, 3^42 for p/q = 3/2, 4/3, 5/3 (computation); they conjecture"
    " non-existence for q < p < q(q-1)",
    "FINITE", "yes", "computation over carry words", "-", "-", "O (below the bound)", "-", "O", "T4 finite range")
add("Narkiewicz1980", "Narkiewicz 1980 (Univ. Beograd Publ. Elektrotehn. Fak. 678-715, 173-174); HKR 2015",
    "[R] via Lagarias 2009 [P]; Holdum-Klausen-Rasmussen INTEGERS 15 (2015) A43 [P]",
    "#{n <= X : (2^n)_3 omits the digit 2} <= 1.62 X^(log_3 2); HKR Thm 2.13: <= 1.3 a^(log_3 2) for"
    " 0 <= s < a; Thm 2.14: any prime p, digits >= p/2 avoided at most 4 a^(log_p((p+1)/2)) times",
    "COUNT", "yes", "2 is a primitive root mod 3^k: the last k digits run over all units", "-", "B", "B", "B", "B",
    "T3 counting only")
add("Lagarias2009", "Lagarias 2009, J. London Math. Soc. 79, 562-588",
    "[P] arXiv math/0512006v4",
    "Thm 1.1: for lambda > 0, #{n <= X : floor(lambda 2^n)_3 omits 2} <= 25 X^0.9725 (X >= n0(lambda)); Thm 1.2:"
    " uncountably many lambda omit 2 along a sparse sequence; Thm 1.4: for nonzero lambda in Z_3,"
    " #{n <= X : (lambda 2^n)_3 omits 2} <= 2 X^(log_3 2); Thm 1.5: dim E^(1) = log_3 2,"
    " (1/2)log_3 2 <= dim E^(2) <= 1/2, (1/6)log_3 2 <= dim E^(3) <= dim E^(2). Conj. A, B: dim 0; Erdos <=> 1 not in E(Z_3)",
    "COUNT", "yes", "3-adic Cantor sets, irrationality measure of log_3 2 (Rhin)", "-", "B", "B", "B", "B",
    "T3 counting and dimension")
add("AbramLagarias", "Abram-Lagarias 2014 (J. Fractal Geom. 1); Abram-Bolshakov-Lagarias 2017 (Exp. Math. 26)",
    "[P] arXiv 1308.3133, 1508.05967",
    "path-set fractals: dim E(Z_3) <= log_3 phi = 0.4380 (ABL Thms 2.5-2.6); dim E^(2) >= log_3 phi,"
    " dim E^(3) >= 0.2284 (AL Thm 5.2)",
    "COUNT", "yes", "automata (path sets), Perron eigenvalues", "-", "B", "B", "B", "B", "T3 dimension")
add("DupuyWeirich2016", "Dupuy-Weirich 2016, J. Number Theory 158, 268-280",
    "[R] via Li-Zhao arXiv 2601.12753 [P] (ScienceDirect captcha)",
    "Thm 3: for primes p != q and a digit b, the Cesaro average over n of the frequency of b among the m least"
    " significant base-q digits of p^n tends to 1/q as m -> infinity; Thm 6: no 'higher Wieferich primes'",
    "COUNT", "yes", "orders of p mod q^m", "-", "B", "B", "B", "B", "none (low digits, averaged)")
add("Saye2022", "Saye 2022, J. Integer Seq. 25, Art. 22.3.4",
    "[P] arXiv 2202.13256v2",
    "for every 16 <= n <= 2*3^45 = 5.9e21, 2^n contains each ternary digit 0, 1, 2 (so Erdos holds for 8 < n <= 2*3^45)",
    "FINITE", "yes", "trailing-digit recursion with u_k = 2*3^(k-1)", "-", "-", "O (below the bound)", "-", "O",
    "T3 finite range")
add("DimitrovHowe2021", "Dimitrov-Howe, Rocky Mountain J. Math. 55 (2025) 45-61",
    "[P] arXiv 2105.06440v4",
    "Thm 1.1: the only powers of 3 that are sums of <= 22 distinct powers of 2 are 3^x, 0 <= x <= 25; Thm 1.2: the"
    " only powers of 2 that are sums of <= 25 distinct powers of 3 are 2^0, 2^2, 2^8 (so for x not in {0,2,8},"
    " 2^x has a ternary digit 2 or at least 26 ternary 1s)",
    "EVERY", "yes (complete, explicit)", "tower of moduli + computation", "-", "B", "O", "-", "O",
    "T3 partial (every-n, bounded complexity)")
add("Yu", "Yu 1989-2007 (Forum Math. 19 (2007) 187-280); Bugeaud-Laurent 1996 (J. Number Theory 61)",
    "[R] via Palojarvi-Seppala Thm A.2 (Yu, rational case) and Pink-Ziegler Thm 2 (Bugeaud-Laurent)",
    "explicit upper bounds for v_p(prod (x_i/y_i)^(b_i) - 1); consequence (lit_B derivation, standard, not found"
    " stated): v_3(2^n - r) <= C(r) log n for rational r not a power of 2, v_2(3^n - r) <= C'(r) log n (effective)",
    "EVERY", "yes", "p-adic linear forms in logarithms", "B", "B", "O", "O (r a power of the base: exact LTE)", "B",
    "T2 partial (terminal landing depth <= C log K)")
add("Ellison1971", "Ellison 1971, Sem. Theorie Nombres Bordeaux 1970-71, exp. 10",
    "[P] Numdam, Thm 3",
    "|2^x - 3^y| > 2^x e^(-x/10) for all positive x, y with x not in {1,...,11,13,14,16,19,27}",
    "EVERY", "yes", "Baker + continued fractions", "B", "B", "O", "O (bounds cycle denominators 2^K-3^L)", "B",
    "cycle half (denominators of rational cycles)")
add("RhinSdW", "Rhin 1987; Simons-de Weger 2005 (Lemma 12)",
    "[R] via Simons-de Weger [P preprint]",
    "|(K+L) log 2 - K log 3| > exp(-13.3 (0.46057 + log K)); irrationality measure of log 3/log 2 at most 8.616",
    "EVERY", "yes", "Pade / linear forms", "B", "B", "O", "O", "B", "cycle half; Lagarias 2009 Lemma 2.2")
add("BLMV2009", "Bourgain-Lindenstrauss-Michel-Venkatesh 2009, ETDS 29, 1705-1722",
    "[P] per lit_B; [R] via Gayfulin-Moshchevitin arXiv 2301.08212 per lit_C",
    "a, b multiplicatively independent, alpha Diophantine-generic (|alpha - p/q| >= k1 q^(-k2)) => {a^u b^v alpha :"
    " u, v <= M} is (log log M)^(-kappa)-dense for M >= M0 (effective Furstenberg)",
    "EVERY", "yes", "additive combinatorics / exponential sums", "B", "B", "O", "O", "B", "none directly")
add("RenRoettger2025", "Ren-Roettger 2025, arXiv 2511.03861",
    "[P] per lit_B (Corollary 1)",
    "the run of zeros following the leading ternary digit of 2^n has length O(log n) (by Baker's theorem)",
    "EVERY", "yes", "archimedean linear forms in logarithms", "-", "B", "O", "-", "B", "T3 partial (top-digit runs)")
add("Spiegelhofer2023", "Spiegelhofer 2023, Israel J. Math. 258, 475-502",
    "[P] arXiv 2105.11173",
    "s_2(n) = s_3(n) for infinitely many n; the count below N is >> N^(log 3/log 4 - delta)",
    "COUNT", "yes", "Fourier analysis of digital functions", "-", "B", "B", "-", "B", "none")
add("BCZ2003", "Bugeaud-Corvaja-Zannier 2003, Math. Z. 243, 79-84",
    "[R] via Cohen-Sonn JTNB 27 (2015)",
    "log gcd(2^n - 1, 3^n - 1) < eps n for all large n (ineffective); infinitely often > exp(c log n/loglog n)",
    "EVERY", "no", "Subspace theorem", "-", "B", "O", "-", "B", "none directly")
add("Koksma1935", "Koksma 1935 (classical metric theorem)",
    "[R] classical; not re-read (citation details UNVERIFIED)",
    "for every theta > 1, (xi theta^n) is uniformly distributed mod 1 for almost every xi",
    "COUNT", "-", "metric (Weyl criterion + second moment)", "B", "B", "B", "B", "B", "none (T4 blocked: DEFECT)")
add("Calegari2005", "Calegari 2005, Int. Math. Res. Not. 2005, no. 20, 1235-1249",
    "[P] arXiv math/0408214; publication checked via zbMATH and Crossref (lit_C)",
    "zeta_p(3) (Kubota-Leopoldt) is irrational for p = 2 and p = 3, and L_2(2, chi_4) is irrational (Thms 3.3, 3.4,"
    " 4.2), via the p-adic criterion Lemma 2.2 and approximations from overconvergent Eisenstein families",
    "EVERY", "yes (explicit exponents theta > 1)", "Beukers' modular proof of Apery transplanted p-adically",
    "-", "-", "O", "-", "O", "none (specific periods; see the note's holonomy remark)")
add("InHouseBCD", "repo: collatz_guards_20260921_discrepancy sec. 2a, 3, 4",
    "[repo] PROVED there; audited by the formalization lane (not in Lean)",
    "no positive-integer 3n+1 orbit has bounded critical discrepancy K_j - j log_2 3 (even after a finite prefix);"
    " extends to 3n+b, b > 0 odd (sec. 4); for b < 0 every positive orbit has discrepancy -> -infinity ((D8));"
    " corollary: no critical-slope mechanical halving word is the word or tail of a positive integer",
    "EVERY", "yes", "capacity count (integers are discrete, orbit grows at most linearly) + density-zero stopping-time lemma",
    "B (true on 3n-1 by (D8))", "O in the proof (needs mean halving 2 > log_2 3); 5n+1 statement OPEN (only (D7))",
    "O", "O (integer capacity; rationals via Proposition B)", "B (all 3n+b, b>0)",
    "T1 and PC on the positive-entropy class BCD (critical slope only)")
add("ThisLaneThmS", "this note: Theorems R, S, Lemma L, Proposition B",
    "[repo] PROVED in collatz_procgen_20260922_transversality_foundry.md",
    "no rational with odd denominator has an eventually Sturmian parity vector under 3x+r (any slope, any intercept);"
    " 5x+1 for slopes < 0.804 (< 0.861 via ADQZ squares); no Z-number has an eventually Sturmian carry word; general"
    " 2-adic repetition criterion (Thm R); PC on bounded-critical-discrepancy words (Prop B, from InHouseBCD)",
    "EVERY", "yes (explicit height bounds)", "2-adic Liouville gap + periodic stretches of Sturmian words",
    "B (Phi_{3,-1} = -Phi)", "B (5x+1 analogue PROVED too)", "O", "O (EP words give x = x')",
    "B (affine 2-adic shift maps with mu < 1.867)", "T1, PC, T4 on zero-entropy classes")

add("MonksYazinski2004", "Monks-Yazinski 2004, Discrete Math. 275, 219-236",
    "[P] author preprint AutoConjV13 (numbering matches the published version)",
    "Thm 2.1: PC <=> Autoconjugacy Conjecture Omega(Q_odd) subset Q_odd <=> no rational 2-adic integer has a divergent"
    " T-orbit; Omega(Z^+) subset Q_odd <=> no divergent positive orbit. Thm 2.7(b): a rational x with a divergent orbit"
    " has liminf kappa_n(x)/n >= ln 2/ln 3 (Lemma L of the note). Thm 2.5: {1,2} is the only self-conjugate integer cycle",
    "EVERY", "yes", "autoconjugacy (complemented parity vector); orbit counting in (1/b)Z", "B", "B", "O", "O", "B",
    "PC/T1 on the subcritical class (Thm 2.7(b))")
add("Akin2004", "Akin 2004, Contemp. Math. 356, 1-20",
    "[R] via zbMATH review (Berthe), Lagarias bibliography II, Rozier 2019 (AMS 403)",
    "Rationality Conjecture for tau_a(x) = x/2, (ax+1)/2 (a odd rational): the parity map Q_a sends Q_odd to Q;"
    " true for a = +-1, false for every non-integer odd rational a; open for odd integers |a| >= 3; heuristic: true for"
    " a = +-3, false for odd |a| >= 5",
    "EVERY", "-", "conjugacy to the shift; genericity of rationals", "B", "O (the heuristic threshold is the drift)", "O", "O", "B",
    "PC framing; the drift barrier in Akin's heuristic")
add("LopezStoll2009", "Lopez-Stoll 2009, Integers 9, A13, 141-162",
    "[P] (Integers online)",
    "Thm 1: a 2-adically convergent series for Phi(1c_alpha) (characteristic Sturmian word) in terms of the convergents"
    " of alpha; Cor. 2: a generalized continued fraction for -1/Phi(1c_alpha); sec. 4: the real sum Phi_R(m_x) is a"
    " strictly increasing devil's staircase on (ln2/ln3, 1]. They state that it is unknown whether any aperiodic word"
    " has an eventually periodic Phi (i.e. PC is open even for Sturmian words). Lopez-Stoll 2021 (arXiv 2101.12747,"
    " unrefereed) claims liminf = ln2/ln3 for divergent rationals; lit_A finds an apparent gap",
    "EVERY", "-", "continued-fraction structure of Sturmian words", "B", "B", "O", "-", "B",
    "prior study of Phi on Sturmian words (irrationality not settled there)")
add("Knight2026", "Knight 2026, Discrete Math. 349, Art. 114812 ('Collatz high cycles do not exist')",
    "[P] HAL preprint hal-04261183; published form via zbMATH",
    "Thm 5.4: no positive integer cycle other than {1,2} has a Christoffel (upper mechanical, rational slope) parity"
    " word; two cycle members differ by 2^(k-2)/(2^k - 3^x), an integer only if |2^k - 3^x| = 1. On the negative"
    " integers the Christoffel cycle {-5,-7,-10} (word 110, 2^3 - 3^2 = -1) exists",
    "EVERY", "yes", "exact difference of two cycle members; Catalan unit gaps 2^k - 3^x = +-1",
    "O (fails for 3n-1: {5,7,10} has word 110)", "-", "O", "O", "O",
    "Collatz cycles on the zero-entropy (Christoffel) class: the periodic twin of Theorem S")
add("ADQZ2001", "Allouche-Davison-Queffelec-Zamboni 2001 (J. Number Theory 91); Berthe-Holton-Zamboni 2006 (Acta Arith. 122)",
    "ADQZ [R] via BHZ intro; BHZ [P] (IMPAN CC-BY)",
    "every Sturmian word begins with infinitely many squares (initial critical exponent >= 2); BHZ Thm 1.1: some word"
    " of slope alpha has initial critical exponent exactly 2 iff a partial-quotient pattern condition holds (forces"
    " unbounded partial quotients); Prop 4.1: every slope has a word with ice <= 1 + golden ratio",
    "EVERY", "-", "Ostrowski numeration / S-adic representation", "-", "-", "O", "-", "B",
    "gives Theorem S directly for mu < 2 (u = empty, v^2 prefixes)")

ORDER_ASKED = ["Furstenberg1967", "SengeStraus1973", "Stewart1980", "Mahler1957", "Mahler1968", "FLP1995",
               "Lagarias2009", "Narkiewicz1980", "DupuyWeirich2016", "HochmanShmerkin2012", "ShmerkinWu2019",
               "Tao2022", "Yu", "Saye2022"]
ORDER_EXTRA = ["MonksYazinski2004", "Akin2004", "LopezStoll2009", "Knight2026", "ADQZ2001", "RudolphJohnson", "Host1995", "CasselsSchmidt", "AbramLagarias", "DimitrovHowe2021",
               "DubickasMossinghoff2009", "Ellison1971", "RhinSdW", "BLMV2009", "RenRoettger2025", "Spiegelhofer2023",
               "BCZ2003", "Koksma1935", "Calegari2005", "InHouseBCD", "ThisLaneThmS"]


def ref(key):
    c = C[key]
    return f"{c['name']} {c['access'].split(' ')[0]}"


def print_catalogue(P):
    P("  Typing rule:")
    for line in TYPING_RULE:
        P("    " + line)
    P("  cols: key | kind | effective | SHEET DRIFT DEFECT INTEGRAL UNIFORM | serves")
    for part, keys in (("items asked for", ORDER_ASKED), ("further items", ORDER_EXTRA)):
        P(f"  -- {part} --")
        for k in keys:
            c = C[k]
            t = c["typing"]
            short = lambda s: s.split(" ")[0]  # noqa: E731
            P(f"  {k:24s} | {c['kind']:6s} | {c['eff'][:22]:22s} | {short(t['SHEET']):2s} {short(t['DRIFT']):2s}"
              f" {short(t['DEFECT']):2s} {short(t['INTEGRAL']):2s} {short(t['UNIFORM']):2s} | {c['serves']}")
            P(f"  {'':24s}   {c['name']}  {c['access']}")
            P(f"  {'':24s}   statement: {c['statement']}")
            long_t = [f"{b}: {v}" for b, v in t.items() if len(v) > 3]
            if long_t:
                P(f"  {'':24s}   notes: " + "; ".join(long_t))
    # pattern
    every = [k for k in ORDER_ASKED + ORDER_EXTRA if C[k]["kind"] == "EVERY"]
    count = [k for k in ORDER_ASKED + ORDER_EXTRA if C[k]["kind"] == "COUNT"]
    drift_o = [k for k in ORDER_ASKED + ORDER_EXTRA if C[k]["typing"]["DRIFT"].startswith("O")]
    P(f"  pattern: EVERY items {len(every)}, COUNT items {len(count)}; items whose DRIFT entry is O: {drift_o}")
    P("  every COUNT item is DEFECT-blind.  Among the EVERY items about integers or single orbits, the Diophantine /")
    P("  transversality ones (Senge-Straus, Stewart, Dimitrov-Howe, Yu, Ren-Roettger, Knight, Theorem S) exclude only")
    P("  zero-entropy classes (few nonzero digits, bounded p-adic closeness, short top runs, Christoffel or Sturmian")
    P("  words).  The two that constrain positive-entropy classes of parity words (Monks-Yazinski 2.7(b): subcritical")
    P("  density; InHouseBCD: bounded critical discrepancy) use growth and capacity, not 2-versus-3 transversality, and")
    P("  work only where the orbit would shrink or grow at most linearly.  Supercritical positive-entropy words: nothing.")


CROSS = [
    ("Collatz divergence (T1)", "OPEN",
     "every-orbit 2-adic non-integrality of Phi on the class HARD (supercritical, unbounded discrepancy, no strong repetitions)",
     "SUB (Monks-Yazinski Thm 2.7(b)), BCD (in-house), STURM (Thm S), REP (Thm R)", "G[HARD].notZ+ (Bernstein family)"),
    ("Collatz cycles", "OPEN",
     "exclusion of positive integer cycle points R_v/(2^|v| - 3^a(v)): a divisibility/carry statement at the convergent clocks",
     "m-cycles m <= 91 (Hercher), circuits (Steiner), Christoffel cycles (Knight 2026, sheet-aware), Ellison/Rhin bounds on 2^K - 3^L",
     "outside the grammar (Theorem R is vacuous on EP words)"),
    ("E-SCC Q1", "OPEN (verified with Collatz below 2^71)",
     "a sheet-sensitive escape family for the hostile points near -1 (-1 - 2^i c/3^j): 2-adic distance of 3-power-denominator rationals",
     "dimension evidence favours 0 (dimension lane); A-L escape is sheet-blind", "mirror of D/F families; not closed"),
    ("E-SCC Q2 (T2)", "OPEN (reduced mod 27; verified below 7.87e17)",
     "chain control plus the terminal landing statement E3: 2^K w_final avoids the hostile 3-adic balls for the dynamically chosen K",
     "terminal depth <= C log K (Yu, D2); kappa formula v_3(2^K - r) = 1 + v_3(K - kappa(r))", "E3 (steering)"),
    ("Erdos ternary (T3)", "OPEN (verified 8 < n <= 2*3^45)",
     "an every-n statement about the TOP ternary digits of 2^n (archimedean), i.e. 1 not in E(Z_3); no 3-adic window statement suffices (B family FALSE)",
     "Narkiewicz/HKR counts, dim E(Z_3) <= log_3 phi, Dimitrov-Howe (>= 26 ones), Saye", "A1"),
    ("Mahler 3/2 (T4)", "OPEN (no Z-number below 2^57)",
     "2-adic non-integrality of Phi_M on the positive-entropy beta-shift S (THM-2228): the T1-type statement for the map ceil(3a/2)",
     "Mahler countability, FLP spread >= 1/3, Sturmian carry words excluded (Thm S)", "M1 (same Bernstein type as T1)"),
    ("divergent 5x+1 orbit exists", "OPEN (expected TRUE)",
     "the reverse transfer: ONE rational (e.g. 7) whose parity vector is not eventually periodic, i.e. irrationality of a parity vector",
     "none (Theorem S shows its divergent orbits cannot be Sturmian for slopes < 0.804)", "outside the grammar (opposite direction)"),
    ("Periodicity Conjecture (PC)", "OPEN",
     "Phi(w) irrational for every non-EP w (= no divergent orbit for any 3x+k, k = +-1 mod 6)",
     "SUB, BCD (Proposition B), STURM, REP", "G[ALL].irr"),
    ("Furstenberg x2x3 measure conjecture", "OPEN",
     "measure rigidity at zero entropy (Rudolph-Johnson settle positive entropy)",
     "Rudolph, Johnson, Host, Hochman-Shmerkin, Shmerkin-Wu", "none: even a proof is COUNT-type, DEFECT-blind for T1-T4"),
]


def print_cross_problem(P):
    P("S5 Cross-problem matrix: problem | status | single missing ingredient | proved partial ingredients | grammar statement")
    for row in CROSS:
        P("  " + " | ".join(row))
    P("  shared statements: the Bernstein family (Phi_T not in Z^+ / not in Q on a word class) serves T1, PC and T4 in the")
    P("  same form (different maps and classes), and Theorem S proves all three on the Sturmian class at once. The 2^K-in-Z_3")
    P("  family serves T2 and T3 with different windows (T2: low digits along dynamically chosen K; T3: top digits). The 5x+1")
    P("  question needs the opposite transfer, and Furstenberg's conjecture serves none of the every-orbit targets.")
