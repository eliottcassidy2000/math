#!/usr/bin/env python3
"""Transversality foundry (collatz-procgen-20260922).

A small grammar of 2-versus-3 statements S = (object, property, quantifier):

  objects     2^K (ternary digits / Z_3), 3^a (binary digits / Z_2),
              2^K w along the loop clocks K0(s) (E-SCC Q2 landings),
              Bernstein numbers Phi_T(w) of parity words w in a class W,
              lambda*2^n and xi*(3/2)^n (real systems of Lagarias and Mahler);
  properties  digit avoidance, number of nonzero digits, p-adic distance to a
              rational, run lengths, rationality / integrality;
  quantifiers all, almost all, infinitely many (over the natural index set).

Each generated statement carries
  status      PROVED (catalogue item or a proof in the note), FINITE-EXACT
              (bounded version checked here), OPEN, FALSE (with counterexample);
  barriers    the same statement for the 3n-1 sheet and for 5n+1 where meaningful,
              plus DEFECT (every-element vs counting), INTEGRAL, UNIFORM typing;
  targets     T1 no divergent positive Collatz orbit, T2 E-SCC Q2 endgame,
              T3 Erdos ternary, T4 Mahler Z-numbers (and PC, Lagarias's
              Periodicity Conjecture, as an extra column).
It prints the catalogue, the statements, the statement x target matrix, the
minimal unblocked statements per target, the Sturmian test and the cross-problem
matrix.  Run with --quick for a fast smoke run.  All arithmetic is exact.
"""
import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import collatz_procgen_20260922_transversality_catalogue as CAT  # noqa: E402
import collatz_procgen_20260922_transversality_engines as EN  # noqa: E402
import collatz_procgen_20260922_transversality_sturmian as ST  # noqa: E402

TARGETS = [
    ("T1", "no divergent positive 3x+1 orbit (Periodicity for Z^+, k=1): Phi(w) not in Z^+ for every non-EP word w"),
    ("T2", "E-SCC Q2 endgame (loops_and_escapes sec. 8): every hostile see-saw chain ends at 2^K w_final in a descending class"),
    ("T3", "Erdos: 2^n has a ternary digit 2 for every n >= 9"),
    ("T4", "Mahler: there is no Z-number (Phi_M(S) contains no positive integer, THM-2228)"),
    ("PC", "Lagarias's Periodicity Conjecture: Phi(w) irrational for every non-EP word w (all 3x+k, k = +-1 mod 6)"),
]
REL_SYM = {"equiv": "<=>", "implies": "=>", "cond": "=>c", "part": "part", "none": "."}


class S:
    """One generated statement."""

    def __init__(self, sid, obj, prop, quant, text, status, evidence, kind,
                 minus=None, five=None, integral="n/a", uniform="n/a", targets=None, blocked=None, finite=None):
        self.sid, self.obj, self.prop, self.quant = sid, obj, prop, quant
        self.text, self.status, self.evidence, self.kind = text, status, evidence, kind
        self.minus, self.five = minus, five
        self.integral, self.uniform = integral, uniform
        self.targets = targets or {}
        self.blocked = blocked or {}
        self.finite = finite


# ------------------------------------------------------------------ grammar
OBJECTS = {
    "2^K": "2^K in base 3 (full expansion) / in Z_3",
    "2^K|low": "2^K mod 3^D (lowest D ternary digits, a 3-adic window)",
    "2^K w|clock": "2^(K0(s)) w, K0(s)=floor((s+1)log2 3) (E-SCC Q2 landing classes)",
    "3^a": "3^a in base 2 / in Z_2",
    "Phi_T(W)": "Bernstein number Phi_T(w) = -sum 2^(d_l) 3^(-l) of parity words w in class W",
    "Phi_M(W)": "Mahler-map Bernstein number Phi_M(c) = -sum c_j 2^j 3^(-j-1) (THM-2228) of carry words c in W",
    "lambda 2^n": "lambda*2^n in R (base-3 digits of the integer part; Lagarias's real system)",
    "xi(3/2)^n": "xi*(3/2)^n mod 1 (Mahler)",
}
PROPERTIES = {
    "digit": "digit avoidance / presence",
    "nnz": "number of nonzero digits",
    "dist": "p-adic (or archimedean) distance to a rational",
    "runs": "run lengths of equal digits",
    "rat": "rationality / integrality",
}
QUANTS = {"all": "for all indices (beyond a threshold)", "almost": "for almost all indices (density 1 / Haar a.e.)",
          "inf": "for infinitely many indices"}
W_CLASSES = {
    "EP": "eventually periodic words",
    "STURM": "Sturmian words (all irrational slopes, all intercepts, lower and upper)",
    "SUB": "non-EP words with liminf density of 1s < log_3 2",
    "BCD": "non-EP words with bounded critical discrepancy sup_s |a_s - s log_3 2| < inf",
    "REP": "non-EP words with strong initial repetitions (gain lambda - log2(|P|+|Q|) unbounded; Theorem R)",
    "LIN": "non-EP words of linear factor complexity",
    "HARD": "non-EP words outside SUB, BCD, STURM and REP",
    "ALL": "all non-EP words",
}


def build_statements(F):
    """F: dict of finite-exact results. Returns the statement list."""
    out = []
    er = F["erdos"]
    tp = F["tprofile"]
    b5 = F["base5"]
    ks = F["v3"]
    ld = F["lowdig"]
    p3 = F["pow3"]
    sp = F["smallpop"]
    land = F["landing"]
    cyc_minus = "3n-1 has the positive cycles {1}, {5,7}, {17,...,91}"
    tij = "FALSE: Z_{5/2}-numbers exist (p = 5 > q^2 = 4; Tijdeman 1972, Flatto 1992 [R via Andrieu-Eliahou-Vivion])"

    # ---- A: Erdos (full expansion) ----
    out.append(S("A1", "2^K", "digit", "all",
                 "for all K >= 9, 2^K has a ternary digit 2",
                 "OPEN", f"Erdos's conjecture; FINITE-EXACT here for 9<=K<={er['Kmax']} (the K<={er['Kmax']} without a 2 are {er['no2']});"
                 f" {CAT.ref('Saye2022')}", "EVERY",
                 minus=None, five=("for all K >= 1, 2^K has a base-5 digit in {2,3,4}",
                                   f"OPEN; FINITE-EXACT for 1<=K<={b5['Kmax']} (K with only digits 0,1: {b5['only01']})"),
                 integral="n/a", uniform="per-instance (one sequence)", targets={"T3": "equiv"},
                 finite=f"FINITE-EXACT 9<=K<={er['Kmax']}"))
    out.append(S("A2", "2^K", "digit", "almost",
                 "for almost all K, 2^K has a ternary digit 2 (exceptions <= 1.62 X^(log_3 2) up to X)",
                 "PROVED", CAT.ref("Narkiewicz1980") + "; " + CAT.ref("Lagarias2009"), "COUNT",
                 five=("for almost all K, 2^K has a base-5 digit in {2,3,4}", "PROVED (same counting, log_5 2 < 1)"),
                 uniform="uniform over bases", targets={"T3": "none"},
                 blocked={"T3": "DEFECT: a counting statement allows a sparse exceptional set"}))
    out.append(S("A3", "2^K", "digit", "inf",
                 "for infinitely many K, 2^K has a ternary digit 2", "PROVED",
                 "odd K: 2^K = 2 mod 3, last digit 2 (elementary)", "COUNT", targets={"T3": "none"},
                 blocked={"T3": "DEFECT"}))
    out.append(S("A4", "2^K", "nnz", "all",
                 "for all K not in {0,2,8}, 2^K has a ternary digit 2 or at least 26 ternary digits 1",
                 "PROVED", CAT.ref("DimitrovHowe2021") + " Thm 1.2", "EVERY", uniform="per-instance (explicit computation)",
                 five=("same shape in base 5", "OPEN (not studied)"), targets={"T3": "part"},
                 blocked={"T3": "bounded complexity only: excludes counterexamples with <= 25 ternary ones"}))
    # ---- B: 3-adic window ----
    for D in (1, 13):
        Kc = 2 * 3 ** (D - 1)
        out.append(S(f"B1[D={D}]", "2^K|low", "digit", "all",
                     f"for all K, the lowest {D} ternary digits of 2^K contain a 2", "FALSE",
                     f"counterexample K = 2*3^{D - 1} = {Kc}: 2^K = 1 mod 3^{D} (LTE), low digits 0..01", "EVERY",
                     minus=("same for -2^K (digits of the negation: avoid 0 above the last digit)", "FALSE"),
                     five=(f"lowest {D} base-5 digits of 2^K contain a digit >= 2", "FALSE (K = 4*5^(D-1): 2^K = 1 mod 5^D)"),
                     integral="B (3-adic limit 2^K -> 1 along K = 2*3^j)", uniform="per-instance",
                     targets={"T3": "none", "T2": "none"},
                     blocked={"T3": "INTEGRAL: periodic in K mod 2*3^(D-1); the 3-adic window cannot see the size of 2^K",
                              "T2": "INTEGRAL: same periodicity"}))
        out.append(S(f"B2[D={D}]", "2^K|low", "digit", "almost",
                     f"for almost all K, the lowest {D} ternary digits of 2^K contain a 2", "FALSE",
                     f"exceptional density exactly 2^({D}-1)/(2*3^({D}-1)) = {2 ** (D - 1) / (2 * 3 ** (D - 1)):.3g}"
                     f" (PROVED: N_j = 2^(j-1) classes mod 2*3^(j-1); FINITE-EXACT j<={max(ld)})", "COUNT",
                     targets={"T3": "none", "T2": "none"}, blocked={"T3": "INTEGRAL", "T2": "INTEGRAL"}))
    # ---- C: nonzero digits ----
    out.append(S("C1", "2^K", "nnz", "all",
                  "the number of nonzero ternary digits of 2^K tends to infinity (>= c log K / log log K)",
                  "PROVED", CAT.ref("Stewart1980") + "; " + CAT.ref("SengeStraus1973"), "EVERY",
                  minus=("same for -2^K", "PROVED"), five=("nonzero base-5 digits of 2^K -> infinity", "PROVED (Stewart)"),
                  uniform="uniform over multiplicatively independent bases", targets={"T3": "none"},
                  blocked={"T3": "DIM: excludes only a zero-dimensional exceptional set; T3's set has dimension log_3 2"},
                  finite=f"FINITE-EXACT: min #nonzero digits over {tp['K0']}<=K<={tp['Kmax']} is {tp['worst_nnz_abs'][0]} (K={tp['worst_nnz_abs'][1]})"))
    out.append(S("C2", "2^K", "nnz", "all",
                  "for all K >= 50, at least 10% of the ternary digits of 2^K are nonzero",
                  "OPEN", "no every-K linear bound is known", "EVERY", targets={"T3": "none"},
                  blocked={"T3": "does not force the digit 2"},
                  finite=f"FINITE-EXACT {tp['K0']}<=K<={tp['Kmax']}: min proportion {tp['worst_nnz'][0]:.3f} (K={tp['worst_nnz'][1]})"))
    out.append(S("C3", "2^K", "nnz", "almost",
                  "for almost all K, 2^K has >= c log K nonzero ternary digits (among its lowest log_3 K digits)",
                  "PROVED", "elementary: 2^K mod 3^D is equidistributed over the units as K runs over a period"
                  " (proof in the note)", "COUNT", targets={"T3": "none"}, blocked={"T3": "DEFECT"}))
    # ---- R: run lengths in the ternary expansion of 2^K ----
    out.append(S("R1", "2^K", "runs", "all",
                  "the run of 0s right after the leading ternary digit of 2^K has length <= C log K (top digits, archimedean)",
                  "PROVED", CAT.ref("RenRoettger2025") + " (Baker)", "EVERY", five=("same in base 5", "PROVED (same method)"),
                  targets={"T3": "part"}, blocked={"T3": "a short top run does not force a digit 2"},
                  finite=f"FINITE-EXACT {tp['K0']}<=K<={tp['Kmax']}: longest such run {tp['top0'][0]} (K={tp['top0'][1]})"))
    out.append(S("R2", "2^K|low", "runs", "all",
                  "the run above the lowest ternary digit of 2^K has length exactly v_3(K) (0s, K even) or 1+v_3(K) (2s, K odd)",
                  "PROVED", "lifting the exponent (elementary); checked for every K in the profile range", "EVERY",
                  five=("same in base 5 with v_5", "PROVED (LTE)"), integral="B (limit points 1, -1 in Z_3)",
                  targets={"T3": "none", "T2": "part"}, blocked={"T3": "3-adic window (INTEGRAL)", "T2": "the landing needs all hostile balls, not only 1 and -1"}))
    out.append(S("R3", "2^K", "runs", "all",
                  "every run of equal ternary digits of 2^K has length <= C log K", "OPEN",
                  "PROVED only at the two ends (R1, R2); the middle digits are the Erdos-type unknown", "EVERY",
                  five=("same in base 5", "OPEN"), targets={"T3": "none"}, blocked={"T3": "does not force a digit 2"},
                  finite=f"FINITE-EXACT {tp['K0']}<=K<={tp['Kmax']}: longest run {tp['run_any'][0]} (K={tp['run_any'][1]}), longest 0-run {tp['run0'][0]} (K={tp['run0'][1]})"))
    # ---- D: 3-adic distance of 2^K to a rational ----
    kid = F["kappa"]
    out.append(S("D1", "2^K", "dist", "all",
                  "v_3(2^K - 1) is bounded", "FALSE", "K = 2*3^j gives v_3 = j+1 (LTE)", "EVERY",
                  five=("v_5(2^K - 1) bounded", "FALSE (K = 4*5^j)"), integral="B", targets={"T2": "none"},
                  blocked={"T2": "INTEGRAL: 1 is a 3-adic limit point of 2^K"}))
    out.append(S("D2", "2^K", "dist", "all",
                  "for every rational r not of the form 2^k, v_3(2^K - r) <= C_r log K (effective)",
                  "PROVED", CAT.ref("Yu") + "; exact form v_3(2^K - r) = [K=e mod 2](1+v_3(K-kappa(r))), "
                  "kappa(r)=log_3(+-r)/log_3(-2) in Z_3 (checked here for r in {1,1/2,5,7,1/5,43/32}, K<=2000)", "EVERY",
                  minus=("v_3(-2^K - r) <= C log K", "PROVED (same theorem)"),
                  five=("v_5(2^K - r) <= C_r log K", "PROVED (same theorem)"), integral="O (fails exactly at r = 2^k)",
                  uniform="uniform over primes and rationals",
                  targets={"T2": "part"}, blocked={"T2": "controls one terminal landing (fixed w_final); not the adversarial chain"}))
    worst_ex = max(v["max_excess_over_log3K"][0] for v in ks.values())
    out.append(S("D3", "2^K", "dist", "all",
                  "for every rational r (3-adic unit, not +-2^k) and eps > 0, v_3(2^K - r) <= (1+eps) log_3 K + C(r,eps)",
                  "OPEN", "equivalent to |K - kappa(r)|_3 >> K^(-1-eps) (a p-adic Lang-Waldschmidt-type bound for kappa(r));"
                  " with eps = 0 it is expected FALSE for generic r (Borel-Cantelli: sum 1/K diverges); for r = +-2^k it holds"
                  " with eps = 0 (LTE)",
                  "EVERY", five=("v_5 analogue", "OPEN"), targets={"T2": "part"},
                  blocked={"T2": "as D2"},
                  finite=f"FINITE-EXACT K<={F['v3_Kmax']}: max of v_3(2^K-r) - log_3 K over the sample r is {worst_ex:.2f}"))
    out.append(S("D4", "2^K", "dist", "almost",
                  "for almost all K, v_3(2^K - r) <= g(K) for any g -> infinity", "PROVED",
                  "density of {K : v_3(2^K - r) >= D} is at most 3^(1-D)/2 (kappa formula)", "COUNT",
                  targets={"T2": "none"}, blocked={"T2": "DEFECT"}))
    # ---- E: Q2 landing along the loop clocks ----
    hit_txt = ", ".join(f"w={a}/{b}: 1-ball {v['hits']['1']} (first s={v['first']['1']}), 1/2-ball {v['hits']['1/2']} (first s={v['first']['1/2']})"
                        for (a, b), v in land["res"].items())
    out.append(S("E1", "2^K w|clock", "dist", "all",
                  f"for all s, 2^(K0(s)) w lies outside the 3-adic balls of radius 3^-{land['D']} around 1 and 1/2",
                  "FALSE", f"K0(s)=floor((s+1)log2 3) is equidistributed mod 2*3^{land['D'] - 1} (Weyl), so each ball is hit"
                  f" with density 1/(2*3^{land['D'] - 1}); FINITE-EXACT s<={land['smax']}: {hit_txt}", "EVERY",
                  minus=("same for -w", "FALSE"), five=("2^(K) w in Z_5 along floor((s+1)log2 5)", "FALSE (same argument)"),
                  integral="B", targets={"T2": "none"}, blocked={"T2": "equidistribution forces hostile landings; T2 needs chain steering"}))
    out.append(S("E2", "2^K w|clock", "dist", "inf",
                  "for infinitely many s, 2^(K0(s)) w lands outside both balls", "PROVED", "Weyl equidistribution (positive density)",
                  "COUNT", targets={"T2": "none"}, blocked={"T2": "DEFECT: the adversary picks the chain"}))
    out.append(S("E3", "2^K w|clock", "dist", "all",
                  "every hostile see-saw chain ends at a landing 2^K w_final in a descending class (steering statement)",
                  "OPEN", "loops_and_escapes sec. 8 (the chain terminal is controlled by the ternary digits of 2^K)", "EVERY",
                  minus=("same statement on the minus sheet (E is sheet-symmetric)", "OPEN"),
                  five=("E_5 backward analogue", "OPEN (untested)"), targets={"T2": "equiv"}))
    # ---- F: 3^a in binary ----
    out.append(S("F1", "3^a", "nnz", "all", "popcount(3^a) -> infinity (>= c log a / log log a)", "PROVED",
                  CAT.ref("Stewart1980"), "EVERY", five=("popcount(5^a) -> infinity", "PROVED"),
                  uniform="uniform", targets={"T1": "none"}, blocked={"T1": "DIM"},
                  finite=f"FINITE-EXACT {p3['a0']}<=a<={p3['amax']}: min popcount {p3['worst_abs'][0]} (a={p3['worst_abs'][1]})"))
    out.append(S("F2", "3^a", "nnz", "all", "popcount(3^a) >= 23 for all a beyond the explicit list", "PROVED",
                  CAT.ref("DimitrovHowe2021") + f"; FINITE-EXACT a<={F['smallpop_amax']}: a with popcount<=22: {sp}",
                  "EVERY", targets={"T1": "none"}, blocked={"T1": "DIM"}))
    out.append(S("F3", "3^a", "dist", "all", "v_2(3^a - 1) is bounded", "FALSE",
                  "v_2(3^a - 1) = 2 + v_2(a) for even a (LTE; checked a < 3000)", "EVERY",
                  five=("v_2(5^a - 1) bounded", "FALSE (= 2 + v_2(a))"), integral="B", targets={"T1": "none"}))
    out.append(S("F4", "3^a", "runs", "all", "the longest run of 0s in the binary expansion of 3^a is <= C log a",
                  "OPEN", "PROVED only for runs just above the lowest bit (LTE: length v_2(a)+1)", "EVERY",
                  five=("same for 5^a", "OPEN"), targets={"T1": "none"},
                  finite=f"FINITE-EXACT {p3['a0']}<=a<={p3['amax']}: longest zero run {p3['run0'][0]} (a={p3['run0'][1]})"))
    # ---- G: Bernstein numbers (Collatz map) ----
    G = []
    G.append(("EP", "irr", "all", "FALSE", "every EP word gives a rational (cycle formula); e.g. (10)^inf -> 1",
              ("same", "FALSE"), ("same", "FALSE"), {"T1": "none", "PC": "none"}, {}))
    G.append(("EP", "notZ+", "all", "FALSE", "(10)^inf -> 1", ("-Phi(w) not in Z^+", "FALSE: " + cyc_minus),
              ("same", "FALSE (13 -> 33 -> 83 -> 13)"), {"T1": "none"}, {}))
    G.append(("STURM", "irr", "all", "PROVED", "Theorem S (this note; elementary; mu=max(1,alpha log2 3) < 1.867; < 2 via ADQZ squares);"
              " literature: open for slopes >= log_3 2 (lit_A; Lopez-Stoll 2009 leave it open)",
              ("same (Phi_{3,-1} = -Phi)", "PROVED"), ("5x+1, slopes alpha < 0.804 (mu<1.867)", "PROVED (Theorem S); OPEN above"),
              {"T1": "part", "PC": "part", "T4": "none"},
              {"T1": "DRIFT: the 5n+1 analogue is PROVED too", "PC": "DRIFT"}))
    G.append(("SUB", "irr", "all", "PROVED", "Monks-Yazinski 2004 Thm 2.7(b) [P] (Lemma L of the note re-proves it; Lagarias 1985 (2.31) for integers)",
              ("same", "PROVED"), ("threshold log_5 2", "PROVED"), {"T1": "part", "PC": "part"},
              {"T1": "DRIFT: the 5n+1 analogue is PROVED too", "PC": "DRIFT"}))
    G.append(("BCD", "notZ+", "all", "PROVED", "in-house discrepancy note sec. 2a (capacity count + density-zero stopping-time lemma)",
              ("minus sheet", "PROVED ((D8): q_j -> 0 on every positive 3n-1 orbit)"),
              ("critical slope log_2 5", "OPEN (in-house sec. 4: the density step has the wrong drift; only (D7) width bound)"),
              {"T1": "part"}, {"T1": "critical slope only: orbits with bounded discrepancy grow at most linearly"}))
    G.append(("BCD", "irr", "all", "PROVED", "assembly (this note): in-house sec. 2a for 3n+k, k>0 (sec. 4) on positive tails, (D8) on negative tails",
              ("same", "PROVED"), ("critical slope log_2 5", "OPEN"), {"PC": "part"}, {"PC": "critical slope only"}))
    G.append(("REP", "irr", "all", "PROVED", "Theorem R (this note; 2-adic Liouville argument)",
              ("same", "PROVED"), ("same criterion with q=5 heights", "PROVED"), {"T1": "part", "PC": "part"},
              {"T1": "DRIFT", "PC": "DRIFT"}))
    G.append(("LIN", "irr", "all", "OPEN", "Theorem R needs repetition exponent > mu; linear complexity only gives > 1",
              ("same", "OPEN"), ("same", "OPEN"), {"PC": "part", "T1": "part"}, {"T1": "DRIFT (expected)", "PC": "DRIFT (expected)"}))
    G.append(("HARD", "notZ+", "all", "OPEN", "T1 restricted to the complement of the proved classes (equivalent to T1)",
              ("no divergent positive 3n-1 orbit with such a word", "OPEN"),
              ("no divergent 5n+1 orbit with such a word", "OPEN (expected FALSE: K-L)"), {"T1": "equiv"}, {}))
    G.append(("ALL", "notZ+", "all", "OPEN", "= T1 (Lagarias 1985; Bernstein 1994)", ("no divergent positive 3n-1 orbit", "OPEN"),
              ("no divergent positive 5n+1 orbit", "OPEN (expected FALSE: Kontorovich-Lagarias)"), {"T1": "equiv", "PC": "none"}, {}))
    G.append(("ALL", "irr", "all", "OPEN", "= PC (Lagarias 1985; Bernstein-Lagarias 1996)", ("same statement", "OPEN"),
              ("PC for 5x+1", "OPEN (expected FALSE)"), {"T1": "implies", "PC": "equiv"}, {}))
    G.append(("ALL", "irr", "almost", "PROVED", "Q is countable and Phi preserves Haar measure", ("same", "PROVED"),
              ("same", "PROVED"), {"T1": "none", "PC": "none"}, {"T1": "DEFECT", "PC": "DEFECT"}))
    for (W, prop, quant, status, ev, minus, five, targets, blocked) in G:
        ptxt = "is irrational" if prop == "irr" else "is not a positive integer"
        qtxt = "for every w in " if quant == "all" else "for Haar-almost every w in "
        out.append(S(f"G[{W}].{prop}.{quant}", "Phi_T(W)", "rat", quant,
                     f"{qtxt}{W} ({W_CLASSES[W]}), Phi(w) {ptxt}", status, ev,
                     "COUNT" if quant == "almost" else "EVERY", minus=minus, five=five,
                     integral=("B" if W == "EP" else "O (EP words excluded; x != x')"),
                     uniform="uniform over affine 2-adic shift maps (sound, not complete)",
                     targets=targets, blocked=blocked))
    # ---- Mahler map ----
    out.append(S("M1", "Phi_M(W)", "rat", "all", "for every carry word c in S (safe: all tails Y_n(c) < 1), Phi_M(c) is not a positive integer",
                  "OPEN", "equivalent to T4 by THM-2228; FINITE-EXACT via Dubickas-Mossinghoff (no Z-number below 2^57, CITED)", "EVERY",
                  five=("ratio 5/2 (carry alphabet {0,1,2})", tij), targets={"T4": "equiv"}))
    out.append(S("M2", "Phi_M(W)", "rat", "all", "for every Sturmian carry word c (in S or not), Phi_M(c) is irrational; hence no Z-number has an eventually Sturmian carry word",
                  "PROVED", "Theorem S, Mahler version (mu = log2 3 < 1.867)", "EVERY",
                  five=("(5/2)-analogue a -> ceil(5a/2)", "OPEN by this method (mu = log2 5 = 2.32 > 1.867)"),
                  targets={"T4": "part"}, blocked={"T4": "Sturmian words are a zero-entropy class; S has entropy log(3/2)"}))
    out.append(S("M3", "Phi_M(W)", "rat", "all", "for every non-EP carry word c, Phi_M(c) is irrational (PC for the Mahler map)",
                  "FALSE", "Phi_M(parity word of 1 under ceil(3a/2)) = 1 and that word is not EP (THM-2228)", "EVERY",
                  targets={"T4": "none"}, blocked={"T4": "the Mahler map expands on both branches, so integers have non-EP words"}))
    # ---- H: real systems ----
    out.append(S("H0", "lambda 2^n", "digit", "all",
                  "for every lambda > 0, floor(lambda 2^n) has a ternary digit 2 for all large n", "FALSE",
                  CAT.ref("Lagarias2009") + " Thm 1.2: uncountably many lambda omit 2 along an infinite sparse sequence", "EVERY",
                  targets={"T3": "none"}, blocked={"T3": "UNIFORM (in lambda): any lambda-uniform argument meets Lagarias's exceptions; T3 needs lambda = 1"}))
    out.append(S("H1", "lambda 2^n", "digit", "almost",
                  "for every lambda > 0, the integer part of lambda 2^n has a ternary digit 2 for all n outside a set of size <= 25 X^0.9725",
                  "PROVED", CAT.ref("Lagarias2009"), "COUNT", targets={"T3": "none"}, blocked={"T3": "DEFECT"}))
    out.append(S("H2", "xi(3/2)^n", "dist", "all",
                  "for every xi > 0, limsup - liminf of frac(xi (3/2)^n) >= 1/3", "PROVED", CAT.ref("FLP1995"), "EVERY",
                  five=("(5/2)^n: range >= 1/5", "PROVED (same theorem, 1/p)"), targets={"T4": "part"},
                  blocked={"T4": "gap: needs > 1/2"}))
    out.append(S("H3", "xi(3/2)^n", "dist", "all",
                  "for every xi > 0, limsup - liminf of frac(xi (3/2)^n) > 1/2", "OPEN", "would imply T4", "EVERY",
                  minus=("ratio -3/2", "OPEN (Lu-Zheng 2026 prove spread >= 11/27)"),
                  five=("ratio 5/2: spread > 1/2 for every xi", tij), targets={"T4": "implies"}))
    out.append(S("H4", "xi(3/2)^n", "digit", "all",
                  "for every xi > 0 there is n with frac(xi (3/2)^n) >= 1/2", "OPEN", "Mahler 1968 (= no Z-number)", "EVERY",
                  five=("ratio 5/2", tij), targets={"T4": "equiv"}))
    out.append(S("H5", "xi(3/2)^n", "dist", "almost",
                  "for almost every xi > 0, (xi (3/2)^n) is uniformly distributed mod 1", "PROVED", CAT.ref("Koksma1935"),
                  "COUNT", targets={"T4": "none"}, blocked={"T4": "DEFECT"}))
    out.append(S("H6", "xi(3/2)^n", "dist", "all",
                  "||(3/2)^n|| > (3/4)^n for all n >= n0 (ineffective)", "PROVED", CAT.ref("Mahler1957"), "EVERY",
                  targets={"T4": "none"}, blocked={"T4": "concerns xi = 1 only"}))
    out.append(S("H7", "xi(3/2)^n", "digit", "all",
                  "each interval [m, m+1) contains at most one Z-number (so they are countable)", "PROVED",
                  CAT.ref("Mahler1968"), "EVERY", targets={"T4": "part"}, blocked={"T4": "countability only"}))
    out.append(S("H8", "xi(3/2)^n", "digit", "all", "no Z-number lies below the Dubickas-Mossinghoff bound",
                  "PROVED", CAT.ref("DubickasMossinghoff2009"), "EVERY", targets={"T4": "part"},
                  blocked={"T4": "finite range"}))
    return out


# ------------------------------------------------------------------ printing
def print_matrix(stmts, out):
    tk = [t for t, _ in TARGETS]
    out("S3.2 Statement x target matrix  (<=> equivalent, => implies, =>c conditional, part = PROVED/OPEN partial"
        " (excludes a thin class of counterexamples), . none).  Status column: PRV/FIN/OPN/FLS.")
    out(f"  {'id':22s} {'status':7s} " + " ".join(f"{t:>5s}" for t in tk) + "   blocked-by (where the relation is part/none)")
    for s in stmts:
        cells = [REL_SYM[s.targets.get(t, "none")] for t in tk]
        bl = "; ".join(f"{t}: {r}" for t, r in s.blocked.items())
        ab = {"PROVED": "PRV", "FINITE-EXACT": "FIN", "OPEN": "OPN", "FALSE": "FLS"}[s.status]
        out(f"  {s.sid:22s} {ab:7s} " + " ".join(f"{c:>5s}" for c in cells) + ("   " + bl if bl else ""))


def minimal_unlocking(stmts, out):
    out("S3.3 Minimal unblocked statements per target: statements that imply the target, are not FALSE, and are not"
        " blocked; among them the weakest (<=> before =>, conditional last).  PROVED partial results listed separately.")
    rank = {"equiv": 0, "implies": 1, "cond": 2}
    for t, desc in TARGETS:
        cands = [s for s in stmts if s.targets.get(t) in rank and s.status != "FALSE" and t not in s.blocked]
        cands.sort(key=lambda s: (rank[s.targets[t]], 0 if "HARD" in s.sid else 1, s.sid))
        parts = [s for s in stmts if s.targets.get(t) == "part" and s.status == "PROVED"]
        out(f"  {t}: {desc}")
        for s in cands:
            five = f" | 5n+1 analogue: {s.five[1]}" if s.five else ""
            minus = f" | 3n-1 analogue: {s.minus[1]}" if s.minus else ""
            out(f"     {REL_SYM[s.targets[t]]:4s} {s.sid}: {s.text}  [{s.status}]{minus}{five}")
        if parts:
            out("     PROVED partial: " + "; ".join(f"{s.sid}" for s in parts))


def print_statements(stmts, out):
    out("S3.1 Generated statements: id | object | property | quantifier | status | statement | evidence")
    out("      barrier profile: 3n-1 analogue | 5n+1 analogue | DEFECT (EVERY sees defects, COUNT is blind) | INTEGRAL | UNIFORM")
    for s in stmts:
        out(f"  {s.sid:20s} | {s.obj:12s} | {s.prop:5s} | {s.quant:6s} | {s.status:7s} | {s.text}")
        out(f"  {'':20s}   evidence: {s.evidence}")
        if s.finite:
            out(f"  {'':20s}   finite: {s.finite}")
        m = f"{s.minus[0]}: {s.minus[1]}" if s.minus else "n/a"
        f5 = f"{s.five[0]}: {s.five[1]}" if s.five else "n/a"
        out(f"  {'':20s}   3n-1: {m} | 5n+1: {f5} | DEFECT: {'sees' if s.kind == 'EVERY' else 'blind'} ({s.kind})"
            f" | INTEGRAL: {s.integral} | UNIFORM: {s.uniform}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="small bounds (smoke run)")
    ap.add_argument("--no-sturmian", action="store_true")
    args = ap.parse_args()
    q = args.quick
    t0 = time.time()
    P = print
    P("collatz_procgen_20260922_transversality: 2-versus-3 transversality foundry")
    P("=" * 100)
    # ---------------- S1 targets and grammar
    P("S1 Targets")
    for t, d in TARGETS:
        P(f"  {t}: {d}")
    P("  grammar: objects " + "; ".join(f"[{k}] {v}" for k, v in OBJECTS.items()))
    P("           properties " + "; ".join(f"[{k}] {v}" for k, v in PROPERTIES.items()))
    P("           quantifiers " + "; ".join(f"[{k}] {v}" for k, v in QUANTS.items()))
    P("           word classes " + "; ".join(f"[{k}] {v}" for k, v in W_CLASSES.items()))
    # ---------------- S2 catalogue
    P("")
    P("S2 Catalogue of PROVED 2-versus-3 transversality theorems (typing: B = blind, O = overcomes, - = not applicable)")
    CAT.print_catalogue(P)
    # ---------------- finite-exact engines
    P("")
    P("S3.0 FINITE-EXACT engines")
    F = {}
    Kmax = 20000 if q else 200000
    F["erdos"] = EN.erdos_scan(Kmax)
    P(f"  Erdos scan 0<=K<={Kmax}: K without ternary digit 2: {F['erdos']['no2']}; deepest lowest-2 position"
      f" {F['erdos']['deepest'][0]} at K={F['erdos']['deepest'][1]}; full conversions {F['erdos']['full_checks']}")
    F["base5"] = EN.base5_scan(Kmax)
    P(f"  base-5 analogue 0<=K<={Kmax}: K with all base-5 digits of 2^K in {{0,1}}: {F['base5']['only01']}")
    tmax = 1500 if q else 20000
    F["tprofile"] = EN.ternary_profile(tmax)
    tp = F["tprofile"]
    P(f"  ternary profile 50<=K<={tmax}: min nonzero proportion {tp['worst_nnz'][0]:.4f} (K={tp['worst_nnz'][1]}),"
      f" min digit-2 proportion {tp['worst_two'][0]:.4f} (K={tp['worst_two'][1]}), min #nonzero {tp['worst_nnz_abs'][0]}"
      f" (K={tp['worst_nnz_abs'][1]}), longest 0-run {tp['run0'][0]} (K={tp['run0'][1]}), longest run {tp['run_any'][0]} (K={tp['run_any'][1]})")
    F["lowdig"] = EN.low_digit_classes(12 if q else 16, brute_max=7 if q else 9)
    P(f"  low-digit classes N_j = #{{K mod 2*3^(j-1): lowest j ternary digits of 2^K in {{0,1}}}}: {F['lowdig']}  (= 2^(j-1), PROVED;"
      " brute force agrees for small j)")
    rs = [(1, 1), (1, 2), (5, 1), (7, 1), (1, 5), (43, 32), (59, 64)]
    F["kappa"] = EN.kappa_identity(rs, 600 if q else 3000)
    P("  kappa identity v_3(2^K b - a) = [K=e mod 2](1 + v_3(K - kappa(a/b))) verified for K <= " + str(600 if q else 3000) + ":")
    for (a, b), v in F["kappa"].items():
        P(f"     r={a}/{b}: e={v['e']}, kappa = ...{v['kappa_low_digits'][::-1]} (base 3, low digits last), max v_3 = {v['max_v3']}")
    F["v3_Kmax"] = 20000 if q else 300000
    F["v3"] = EN.v3_scan(rs[1:], F["v3_Kmax"])
    for (a, b), v in F["v3"].items():
        P(f"  v_3 scan r={a}/{b}, K<={F['v3_Kmax']}: max v_3 {v['max_v3'][0]} at K={v['max_v3'][1]};"
          f" max(v_3 - log_3 K) = {v['max_excess_over_log3K'][0]:.3f} at K={v['max_excess_over_log3K'][1]}")
    amax = 3000 if q else 30000
    F["pow3"] = EN.pow3_binary(amax)
    p3 = F["pow3"]
    P(f"  3^a in binary, 20<=a<={amax}: min popcount proportion {p3['worst_ratio'][0]:.4f} (a={p3['worst_ratio'][1]}),"
      f" min popcount {p3['worst_abs'][0]} (a={p3['worst_abs'][1]}), longest 0-run {p3['run0'][0]} (a={p3['run0'][1]}),"
      f" longest low run {p3['low_run'][0]} (a={p3['low_run'][1]}); LTE v_2(3^a-1) checked a<3000")
    F["smallpop_amax"] = 5000 if q else 30000
    F["smallpop"] = EN.small_popcount_pow3(F["smallpop_amax"])
    P(f"  exponents a<={F['smallpop_amax']} with popcount(3^a) <= 22: {F['smallpop']}")
    smax = 1000000 if q else 20000000
    F["landing"] = EN.beatty_landing([(1, 1), (5, 1), (7, 1)], smax, D=13)
    P(f"  Q2 landings: K0(s)=floor((s+1)log2 3), first values {F['landing']['K0_head']}; for s<={smax}, hits of the"
      f" 3-adic balls v_3(2^K0(s) w - h) >= 13, h in {{1, 1/2}}:")
    for (a, b), v in F["landing"]["res"].items():
        P(f"     w={a}/{b}: {v['hits']} first s {v['first']}   (expected per ball ~ {smax / (2 * 3 ** 12):.2f})")
    # ---------------- S3 statements
    P("")
    stmts = build_statements(F)
    P(f"S3 Transversality statements generated: {len(stmts)}")
    print_statements(stmts, P)
    P("")
    print_matrix(stmts, P)
    P("")
    minimal_unlocking(stmts, P)
    # ---------------- S4 Sturmian
    if not args.no_sturmian:
        P("")
        N = 4000 if q else 20000
        ST.sturmian_report(N=N, printer=P)
        P("")
        ST.control_report(N=N, printer=P)
        if not q:
            P("")
            ST.deep_report(N=100000, printer=P)
    # ---------------- S5 cross-problem matrix
    P("")
    CAT.print_cross_problem(P)
    P("")
    sys.stderr.write(f"[transversality] done in {time.time() - t0:.1f} s\n")


if __name__ == "__main__":
    main()
