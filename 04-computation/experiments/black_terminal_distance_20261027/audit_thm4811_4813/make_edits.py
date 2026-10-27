"""Build required_edits.json and validate every `old` string against the audited text (origin/main eca666b44c copies in latest/).

Each old string must occur exactly once in its file; the script also applies all edits to a scratch copy to make sure
they compose (no edit destroys another's anchor).
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LATEST = os.path.join(HERE, "latest")
F11 = "01-canon/theorems/THM-4811-arc-strong-connectivity-is-the-exact-complement-distance-to-the-nonstrong-locus-black-refinement-and-the-dipole-law.md"
F13 = "01-canon/theorems/THM-4813-collatz-numerator-coordinate-forcing-cocycle-gauss-sum-intertwiner-and-sharp-spectral-exception-budgets.md"
FN = "05-knowledge/results/dipole_calculus_black_distance_tc_gaps_collatz_numerator_20261027.md"
FR = "07-reflections/why-the-advances-work-ii-potentials-symmetry-budgets-positive-coordinates-macmini-20261027.md"

E = []


def add(f, old, new, reason, finding):
    E.append({"file": f, "old": old, "new": new, "reason": reason, "finding": finding})


# ------------------------------------------------------------------------------------------------ MAJOR 1 (THM-4813 A3)
add(F13,
    "on every native cycle the total forcing is exactly (T log 2 - d log 3)/2 (the Steiner-Eliahou cycle identity of THM-4633)",
    "on every native cycle the forcing cocycle sums to exactly T log 2 - d log 3, i.e. sum artanh(lambda) = (T log 2 - d log 3)/2 while the price sum lambda is strictly smaller (the Steiner-Eliahou cycle identity of THM-4633, factor by factor at the lifted Syracuse elements)",
    "Factor-2 inconsistency: the title itself defines the forcing cocycle as 2 artanh(lambda) = log(1 + 2^a/(3X)); its sum over a closed word is T log 2 - d log 3 (auditor: 61 rational cycles at 50 digits; trivial cycle 2 artanh(1/7) = log(4/3) = 2 log 2 - log 3, not half of it).",
    "M1")
add(F13,
    "the total forcing is `(T log 2 - d log 3)/2 > 0`: the native\nenergy identity pays, around a cycle, exactly half the logarithmic defect of `2^T/3^d`.",
    "the forcing cocycle sums to `sum_i 2 artanh(lambda_i) = T log 2 - d log 3 > 0`, i.e.\n`sum_i artanh(lambda_i)` is exactly half the logarithmic defect of `2^T/3^d`; the price actually paid by the energy\nidentity, `sum_i lambda_i`, is strictly smaller (`artanh(x) > x` for `x > 0`).",
    "Same factor-2 inconsistency (A2 calls 2 artanh(lambda) the forcing cocycle; the sentence's own trivial-cycle example has forcing log(4/3) = T log 2 - d log 3 per lap). 'Pays exactly half' is false for the payment sum lambda, which is strictly below half.",
    "M1")
add(F13,
    "observation, the total forcing `F(s)` of the native orbit of a primitive seed `s < 2 10^6` down to `1` never exceeds\n`0.11275`;",
    "observation, the total price `F(s) = sum lambda` of the native orbit of a primitive seed `s < 2 10^6` down to `1` never\nexceeds `0.11275` (its forcing `sum 2 artanh(lambda)` is about twice that, at most `0.225534` there);",
    "F(s) in forcing_tail.c / HYP-9421 (4) is the price sum lambda; under A2's definition of forcing (2 artanh lambda) the maximum below 2 10^6 is 0.225533682 (auditor price_sum.c), so 'total forcing never exceeds 0.11275' is false as worded.",
    "M1")
add(FN,
    "Around a\nnative cycle the energy identity pays exactly `(T log 2 - d log 3)/2`: the classical Steiner-Eliahou identity\n(THM-4633) in forcing language.",
    "Around a\nnative cycle the forcing cocycle sums to exactly `T log 2 - d log 3` (so `sum artanh(lambda) = (T log 2 - d log 3)/2`,\nand the price `sum lambda` actually paid is strictly smaller): the classical Steiner-Eliahou identity (THM-4633) in\nforcing language.",
    "Same normalisation error as THM-4813 A3; the payment sum lambda is strictly below half the defect.",
    "M1")

# ------------------------------------------------------------------------------------------------ MAJOR 2 (Session75 section 6 credit)
add(F13,
    "  (C) the Lipschitz constants from K' = -K + 2A_Q (THM-4811 (A)) and the path spectrum 2i cos(k pi/(n+1)); the",
    "  (C) the operator-norm Lipschitz constant 4 cos(pi/(n+1)) from K' = -K + 2A_Q and the path spectrum 2i cos(k pi/(n+1));\n  CITED: the matrix law K(tau_P T) = -K(T) + 2L, the trace-norm tariff 2||L||_* = 2 sum_k |2 cos(k pi/(n+1))| ~ 8n/pi\n  and the linear clock bound t >= (n(n-1) - C_b n^(3/2) log2 n)/(4(n-1)) = n/4 - O(sqrt(n) log n) are Session75's\n  audited RESULT section 6, (6.1)-(6.3) and item 5 (C), pushed (d5f10e1cf2) before this file; the",
    "Session75 RESULT section 6 (on main, audited, 18 minutes before THM-4813's first commit) already proves (6.1) K(tau_P T) = -K(T) + 2L, (6.2)-(6.3) the trace-norm tariff with the path eigenvalues 2cos(j pi/(n+1)) and the 8n/pi asymptotics, and a stronger linear clock bound (n/4) via exception tariffs; THM-4813 typed these as its own PROVED content.",
    "M2")
add(F13,
    "which is linear in `q`, while THM-4810's score potential gives `(q^2-1)/8`",
    "which is linear in `q` (precisely `q/(2 pi) - sqrt(q)/4 + O(1/q)`; `0.94` at `q = 11`) and weaker than the linear clock\nbound `q/4 - O(q^(1/2) log q)` already proved by Session75 (RESULT section 6 and item 5 (C): exception tariff `n - 1` per\nmove; its trace-norm tariff (6.2)-(6.3) is the second Lipschitz constant above), while THM-4810's score potential gives `(q^2-1)/8`",
    "Credit and calibration: the op-norm bound is the weakest of the available linear bounds; Session75 already has the trace-norm tariff and the q/4 bound; the 'about q/(2 pi)' is only the leading term.",
    "M2")
add(FN,
    "| `||K||_op`, `||K||_*` | `4 cos(pi/(n+1))`, `2 sum |2 cos(k pi/(n+1))|` | linear move bounds only (THM-4813 (C)) | far from tight |",
    "| `||K||_op`, `||K||_*` | `4 cos(pi/(n+1))`, `2 sum |2 cos(k pi/(n+1))|` (trace-norm tariff: Session75 RESULT (6.2)-(6.3)) | linear move bounds only (THM-4813 (C); Session75's exception tariff gives `n/4 - O(sqrt(n) log n)`) | far from tight |",
    "Credit the already-published Session75 section 6 results in the table.",
    "M2")
add(F11,
    "concurrent dyadic-module-spectrum session's matrix form K' = -K + 2A_P;",
    "concurrent dyadic-module-spectrum session's matrix form K' = -K + 2A_P (also Session75 RESULT (6.1), on main and audited\n  before this file, with the path spectrum and the trace-norm tariff (6.2)-(6.3));",
    "Session75 (d5f10e1cf2, 01:22) states K(tau_P T) = -K(T) + 2L as (6.1) before THM-4811's first commit (8ab6b3e358, 01:40).",
    "M2")

# ------------------------------------------------------------------------------------------------ MAJOR 3 (misattribution to Session75)
add(F13,
    "  three ROOT rows (2+omega)(a + b omega), and the trace-norm bound for bounded-arity substitution parses (T1-T6) are the\n  concurrent dyadic-module-spectrum session's (worktree dyadic_module_spectrum_20261026), published in audited form by\n  Session75 (05-knowledge/results/session75_overnight_black_root_spectral_20261027);",
    "  three ROOT rows (2+omega)(a + b omega) and the trace-norm theorem T1-T6 are the concurrent dyadic-module-spectrum\n  session's (worktree dyadic_module_spectrum_20261026, not on main, not audited); the integral ROOT rows themselves,\n  c^r = 2 Re(omega^r A) = (2a - b, -a - b, -a + 2b), are root_eisenstein_first_arrival_20261018j E1/E2 (recovered by\n  Session75's inheritance report, which states that the directed-C3 map was NOT source-recovered); the audited\n  bounded-arity nuclear-norm cap C_b n log2 n is Session75's own theorem\n  (05-knowledge/results/session75_overnight_black_root_spectral_20261027, session75-spectral/RESULT.md);",
    "False provenance: Session75's packet contains no three-chart lift, no beta = 1 + 2 omega, no (2+omega) row form and no T1-T6 (grep of the whole packet); its INHERITANCE_REPORT says 'The prompt's directed-C3 map is not yet source-recovered, although these exact integral rows are'. Its nuclear cap has its own constant C_b = 2 sqrt(alpha(1-alpha))/h2(alpha), not the dyadic all-tree constant cited by THM-4814.",
    "M3")

# ------------------------------------------------------------------------------------------------ MAJOR 4 (black-classification overclaims in note / reflection)
add(FN,
    "| arc connectivity `lambda = min_X f(X)` | at most 1 | distance to the non-strong locus | **always** (any colour); black: all but `C3`, `R5`, family `E` (section 2) |",
    "| arc connectivity `lambda = min_X f(X)` | at most 1 | distance to the non-strong locus | **always** (any colour); black: through order 10 all but `C3`, `R5`, `E10` (family `E(M; Z)` proved exceptional at every order; completeness OPEN, HYP-9421 (1)) |",
    "Stated the OPEN classification (HYP-9421 (1)) as an all-order fact; THM-4811 proves only the census through order 10 and the exceptionality of E(M;Z).",
    "M4")
add(FR,
    "| same, black moves only (THM-4811 (C)) | `lambda` | 1 | a black word into a minimum cut | all but `C3`, `R5`, `E(M; Z)` |",
    "| same, black moves only (THM-4811 (C)) | `lambda` | 1 | a black word starting in a minimum out-cut and ending outside it | through order 10 all but `C3`, `R5`, `E10` (FINITE-EXACT); `E(M; Z)` never exact (PROVED); no further exception is OPEN (HYP-9421 (1)) |",
    "Same overclaim; the reflection promises to carry the theorem's status.",
    "M4")

# ------------------------------------------------------------------------------------------------ MINOR
add(F13,
    "and the fixed points of the B-free native edges are exactly THM-4600's anchors c_a = 1/(2^a - 3).",
    "and the fixed points of the B-free native edges (a = 1, 2) are THM-4600's anchors c_1 = -1 and c_2 = 1 (every anchor c_a = 1/(2^a - 3) is the fixed point of the Syracuse letter X -> 3X/2^a + 1; the B-layered native edges, j >= 1, fix non-anchors such as -1/5, -3/13).",
    "Only c_1, c_2 arise from native edges (a in {1,2}); c_3 = 1/5, c_4 = 1/13, ... are fixed points of Syracuse letters, not of native edges (auditor A4 check).",
    "m1")
add(F13,
    "Anchored debts are the fixed points of the\n`B`-free native edges read in the numerator coordinate.",
    "In the numerator coordinate every anchor `c_a` (all `a >= 1`) is the\nfixed point `2^a/(2^a - 3)` of the Syracuse letter `X -> 3X/2^a + 1`; the `B`-free native edges are the letters `a = 1, 2`.",
    "Overgeneralisation (see previous edit).",
    "m1")
add(FN,
    "the\nfixed points of the `B`-free edges are THM-4600's anchors `c_a = 1/(2^a - 3)`.",
    "the\nfixed points of the `B`-free edges are THM-4600's anchors `c_1 = -1`, `c_2 = 1` (every `c_a = 1/(2^a - 3)` is the fixed point of\nthe Syracuse letter `X -> 3X/2^a + 1`).",
    "Same overgeneralisation.",
    "m1")
add(FR,
    "at the fixed points it is THM-4600's anchors `c_a = 1/(2^a - 3)`.",
    "at the fixed points of the Syracuse letters `X -> 3X/2^a + 1` it is THM-4600's anchors `c_a = 1/(2^a - 3)` (the `B`-free native edges give `c_1 = -1`, `c_2 = 1`).",
    "Same overgeneralisation.",
    "m1")
add(F11,
    "So every strong tournament of order at most 9 has `lambda = delta^0`",
    "So every strong tournament of order at most 9 has `lambda = delta^0` (for `delta^0 = 1` because strong means `lambda >= 1`;\nfor `delta^0 >= 2` because `n <= 9 <= 4 delta^0 + 1`)",
    "The criterion n <= 4 delta^0 + 1 does not cover delta^0 = 1 at orders 6-9; the missing step is lambda >= 1 for strong T.",
    "m2")
add(F11,
    "Cuts meeting `D` lose at least `|X n D| |D* \\ X| - 1` arcs to `D*` plus\n   at least one arc inside the strong `D`; the case analysis leaves only `X = D*`.",
    "If `X` meets `D` and `X n D != D`, then `f(X) >= e(X n D -> D \\ X) + e(X n D -> D* \\ X)\n   >= 1 + |X n D| |D* \\ X| - 1`, so `f(X) = 1` forces either `D* \\ X = {w*}`, `X n D = {u}` (then `f(X) >= s_D(u) >= 1` plus\n   one cross arc unless `(u, w*) = (x0, y0)`, and in that case `f(X) >= |M| >= 3`), or `D*` contained in `X`, `x0 in X` and exactly one arc\n   of `D` entering `W = D \\ X` (impossible by the source-of-`M` argument above, applied to `W`). If `D` is contained in\n   `X != V`, then `f(X) >= |D| |D* \\ X| - 1 >= 4`. Hence `X = D*` is the only minimum out-cut.",
    "The sketch omits the case X containing D* (where the cross-arc term is vacuous and the source-of-M argument is needed again) and the case X containing D (no internal D arc); the replacement is a complete case analysis (auditor verified uniqueness of D* on all 15 members through order 16).",
    "m3")
add(F11,
    "Every black line at a non-strong class is such an edge (by (B) its other end has\n`lambda = 1`), and by (B)-(C) every non-pure-blue strong class of order at most 10 is joined to a `lambda = 1` class by\nblack moves.",
    "Every black line at a non-strong class is such an edge (its other end is strong by THM-4729 L3, so has\n`lambda = 1` by (B)), and every non-pure-blue strong class is joined by black moves to a `lambda = 1` class that has a\nblack descent: at every order by the inherited strong black `c3`-descent (black_descent_20261020c; the last strong class\nof a black descent chain has `lambda = 1` and its last move is a black descent word), and through order 10 also by (C1).",
    "(B) alone gives only lambda <= 1 at the other end; strongness is THM-4729 L3 (no line joins two non-strong nodes). The reduction is valid at every order via the inherited all-order black c3-descent, not only for n <= 10.",
    "m4")
add(F11,
    "Black-core connectivity (HYP-9348 (b)) is not\nclaimed: (B)-(C) connect each strong class to the non-strong locus, not the non-strong basins to each other.",
    "Black-core connectivity (HYP-9348 (b)) is not\nclaimed beyond order 10: (B)-(C) connect each strong class to the non-strong locus, not the non-strong basins to each\nother; by (E) the general case reduces to the connectivity of `B_n` (a sufficient condition).",
    "(F) contradicted the new (E) heading 'Black-core connectivity through order 10'; B_n connectivity is a sufficient condition (not shown necessary).",
    "m4")
add(F11,
    "black distance exactly 2 (two-move\n  certificates) for the six members through order 14;",
    "black distance exactly 2 (two-move\n  certificates) for six of the seven members through order 14 (auditor: all fifteen members through order 16 have no\n  black descent, a unique minimum cut D* and black distance exactly 2);",
    "There are seven members through order 14 (orders 10, 12, 14: 1 + 2 + 4); the producer built six (missing E(C3 => K1; K1 => K1)). The auditor's family check covers all 15 members through order 16.",
    "m5")
add(FN,
    "Members through order 16\n  have no black descent; members through order 14 have black distance exactly 2.",
    "Members through order 16\n  have no black descent; all fifteen members through order 16 have black distance exactly 2 (producer: six of the seven\n  members through order 14; independent audit: all fifteen).",
    "Producer evidence covered six of the seven members through order 14.",
    "m5")
add(F11,
    "descent criterion of (B) verified literally on all legal words through order 8 (2,073,672 words at order 8, no failure);",
    "descent criterion of (B) verified literally on all legal words of all strong classes through order 8 (2,073,672 words at\n  order 8, no failure; independent audit: every class, strong or not, through order 9, 2,143,304 words at order 8 and\n  270,447,284 at order 9, no failure);",
    "The producer's check ran on strong classes only (gentourng -c); the count 2,073,672 is the strong-class word count.",
    "m6")
add(F11,
    "Protect a vertex\n`v` attaining `delta^0`, say `s(v) = k` (otherwise pass to the converse).",
    "Protect a vertex\n`v` attaining `delta^0`, say `s(v) = k` (otherwise pass to the converse); this is the protected-vertex route of the\nconcurrent regular basin theorem (B3), run here with the cut potential.",
    "Credit: the note attributes the protected-vertex route to the aliquot_power_124 basin theorem B3; the theorem file does not.",
    "m7")
add(F11,
    "For a regular `(2m+1)`-tournament the distance is `m`: the concurrent regular\nbasin theorem's cut-score bound is the case `lambda = m` of (B), and (B) shows the bound is attained from every source.",
    "For a regular `(2m+1)`-tournament the distance to the non-strong locus is\n`m`: the concurrent regular basin theorem's cut-score lower bound (for its smaller target, a singleton plus a strong\ncore) is the case `lambda = m` of (B), and (B) shows that the non-strong-locus bound `lambda(T)` is attained from every source.",
    "The basin theorem's target is CORE, a sub-locus; the sentence conflated the two targets.",
    "m8")
add(F13,
    "This is the precise reason for the shared\nspectrum;",
    "This is the precise reason for the shared\nsingular values (`|2 + omega| = |1 + 2 omega| = sqrt(3)`; the eigenvalues differ by the unit, `(3 +- i sqrt(3))/2` against\n`+- i sqrt(3)`);",
    "-RK and K have different eigenvalues on Lambda_3 (auditor: 1.5 +- 0.866i vs +-1.732i); only the moduli agree.",
    "m9")
add(FN,
    "(our THM-4811 (B) is the all-source version, (C) the black\n  census through order 10)",
    "(our THM-4811 (B) gives the exact distance lambda(T) to the larger non-strong locus from every source, (D1) the\n  singleton distance delta^0 whenever lambda = delta^0, (C) the black census through order 10)",
    "THM-4814 (D) is about the singleton locus; THM-4811 (B) is about the non-strong locus.",
    "m10")
add(FR,
    "Two of tonight's claims (the sharp spectral constant, the doubly-regular\n  lattice module) had concurrent published versions by the time they were filed.",
    "Three of tonight's claims (the sharp spectral constant, the doubly-regular\n  lattice module, and the trace-norm move tariff with its linear clock bound, Session75 RESULT section 6) had concurrent\n  published versions by the time they were filed.",
    "Session75 section 6 also predates THM-4813 (C).",
    "m11")

add(FR,
    "klein's THM-4814 (D), and\nTHM-4811 as the all-source case).",
    "klein's THM-4814 (D), and\nTHM-4811: (B) gives the exact distance to the larger non-strong locus from every source, while the all-source black-core\nversion is FINITE-EXACT only through order 9 (D2) and conjectured in general (HYP-9421 (6))).",
    "THM-4811 proves the all-source statement for the non-strong locus, not for the black-core target (CORE, black moves), which is D2 (orders 6-9) plus HYP-9421 (6).",
    "m12")
FH = "05-knowledge/hypotheses/HYP-9421-black-descent-exceptions-tc-spectrum-and-the-native-forcing-tail.md"
add(FH,
    "Evidence: order <= 10 complete; six family\nmembers through order 14 have black distance exactly 2.",
    "Evidence: order <= 10 complete; all fifteen family\nmembers through order 16 have black distance exactly 2 (producer: six of the seven through order 14; independent audit\nof THM-4811/4813: all fifteen).",
    "RELATED (outside the four audit targets, same miscount): there are seven members through order 14.",
    "m5-related")


def main():
    texts = {}
    for f in {e["file"] for e in E}:
        with open(os.path.join(LATEST, os.path.basename(f)), encoding="utf-8") as fh:
            texts[f] = fh.read()
    ok = True
    for e in E:
        c = texts[e["file"]].count(e["old"])
        if c != 1:
            ok = False
            print(f"BAD ({c} occurrences): {e['finding']} {e['file']}: {e['old'][:90]!r}")
    # compose
    work = dict(texts)
    for e in E:
        if work[e["file"]].count(e["old"]) != 1:
            ok = False
            print(f"COMPOSE FAIL: {e['finding']} {e['old'][:80]!r}")
        work[e["file"]] = work[e["file"]].replace(e["old"], e["new"], 1)
    print("all old strings unique and composable:", ok, "edits:", len(E))
    with open(os.path.join(HERE, "required_edits_vs_eca666b44c.json"), "w", encoding="utf-8") as fh:
        json.dump(E, fh, indent=2, ensure_ascii=False)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
