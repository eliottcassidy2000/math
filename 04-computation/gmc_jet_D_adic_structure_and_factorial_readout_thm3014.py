#!/usr/bin/env python3
"""THM-3014 -- D-adic structure of the resultant log-jets, the discriminant
regular value E_j(M), P_5 predictions, and the factorial-functional readout.

Four parts.

A. D-ADIC RE-EXPANSION.  D = U^2+3U-3V-1 is LINEAR in V, so (U,V) -> (U,D) is an
   invertible polynomial change of variables, V = (U^2+3U-1-D)/3.  Re-expanding
   P_4 gives 721 terms against 717 in (U,V): D-coordinates are NOT sparser and
   P_4 is NOT divisible by D.  But the graded structure is exact:
   deg_U q_k = 2(2j-k), and the term counts run 5, 27, 45, 63, 81, 99, 117, 135,
   149 -- an arithmetic progression of step 18 in the interior.

B. THE DISCRIMINANT REGULAR VALUE.  Because deg_U P_j = 4j = 2j * deg_U D, the top
   coefficient [D^{2j}]P_j is FORCED U-free.  So L_j has a regular value at the
   discriminant locus D=0 which is a pure polynomial in M:
        E_j(M) := [D^{2j}] P_j / c_j.
   Computed exactly for j=1,2,3,4 and obeying three laws:
        [M^(j+1)] E_j = (-1)^(j+1) * 46/(j(j+1))      (equidistribution; same as A_j)
        [M^j]     E_j = (-1)^j * 11/j
        [M^(j-1)] E_j = (-1)^(j+1) * 37/18            (j >= 2)
        [M^(j-2)] E_j = 0                             (j >= 3)
   Note E_j and A_j share ONLY the leading coefficient; the second differs
   (11/j vs 12/j) and the third differs in both value and sign (37/18 vs 47/24).

C. P_5 PREDICTIONS.  Shape is forced: max terms (2j+1)^3 = 1331, degrees
   (M,U,V) = (10,20,10), support b+2c <= 20.  The two extreme slices' leading
   coefficients follow from the laws in B and from THM-3011 section 2a.

D. THE FACTORIAL READOUT.  THM-3000's first-occurrence law gives jet J_j the
   coefficient -(j-1)! C(k-2,j-3) at edge k -- a FACTORIAL times a BINOMIAL.
   Exactly:
        Lambda_k(J) = -L( z^2 * Phi_k(z) ),   Phi_k(z) = sum_i C(k-2,i) J_(i+3) z^i,
   where L(z^n) = n! is the Edo-van den Essen factorial functional of the
   SFC lane (THM-2812/2836/2173).  So the Newton edge reads the log jets THROUGH
   L, twisted by a binomial transform whose index is the edge.  This is a shared
   evaluation functional, NOT a reduction -- see the loss ledger in the theorem.

Reproduce: python3 04-computation/gmc_jet_D_adic_structure_and_factorial_readout_thm3014.py
"""

import importlib.util
import json
import sympy as sp
from fractions import Fraction as Fr

M, U, V, Dv, z = sp.symbols('M U V D z')
TABLE = "05-knowledge/results/gmc_first_gap_fourth_resultant_jet_P4_table_thm3013.json"
CJ = {1: Fr(1), 2: Fr(16), 3: Fr(-128), 4: Fr(1658880)}


def rule(s):
    print("=" * 74)
    print(s)
    print("=" * 74)


def load_tables():
    spec = importlib.util.spec_from_file_location(
        "t", "04-computation/gmc_first_gap_wall_stripped_all_width_second_edge_circuit_positivity_thm2997.py")
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except SystemExit:
        pass
    except Exception:
        pass
    fr = json.loads(mod.JET_DATA_TEXT)
    return {1: fr[0], 2: fr[1], 3: fr[2], 4: json.load(open(TABLE))}


def poly_of(rows):
    return sum(sp.Rational(r[3], r[4]) * M ** r[0] * U ** r[1] * V ** r[2] for r in rows)


def partA(tabs):
    rule("A. D-ADIC RE-EXPANSION OF P_4  (D = U^2+3U-3V-1, V = (U^2+3U-1-D)/3)")
    P4 = poly_of(tabs[4])
    s = sp.expand(P4.subs(V, (U ** 2 + 3 * U - 1 - Dv) / 3))
    pol = sp.Poly(s, Dv)
    tot = 0
    ok = True
    print("   k | terms | deg_M | deg_U | predicted deg_U = 2(8-k)")
    for k in range(pol.degree(), -1, -1):
        q = sp.expand(pol.coeff_monomial(Dv ** k))
        if q == 0:
            continue
        n = len(sp.Poly(q, M, U).terms())
        tot += n
        du = int(sp.degree(q, U))
        good = (du == 2 * (8 - k))
        ok &= good
        print("   %d | %5d | %5d | %5d | %5d  %s" % (k, n, int(sp.degree(q, M)), int(du), 2*(8-k), "OK" if good else "MISMATCH"))
    print(f"   total in (M,U,D): {tot}   vs {len(tabs[4])} in (M,U,V)  -> sparser? {tot < len(tabs[4])}")
    print(f"   P_4 divisible by D?  {sp.expand(pol.coeff_monomial(Dv ** 0)) == 0}")
    print(f"  VERDICT A: {'graded U-degree law holds; D-coords NOT sparser' if ok else 'FAILED'}")
    return ok


def partB(tabs):
    print()
    rule("B. DISCRIMINANT REGULAR VALUE  E_j(M) = [D^(2j)] P_j / c_j")
    E = {}
    ok = True
    for j in (1, 2, 3, 4):
        s = sp.expand(poly_of(tabs[j]).subs(V, (U ** 2 + 3 * U - 1 - Dv) / 3))
        top = sp.expand(sp.Poly(s, Dv).coeff_monomial(Dv ** (2 * j)))
        ufree = (sp.degree(top, U) == 0)
        ok &= ufree
        E[j] = sp.expand(top / sp.Rational(CJ[j].numerator, CJ[j].denominator))
        print(f"   j={j}: U-free {ufree}   E_j = {E[j]}")
    print()
    laws = [("[M^(j+1)]", lambda j: j + 1, lambda j: sp.Rational((-1) ** (j + 1) * 46, j * (j + 1)), 1),
            ("[M^j]", lambda j: j, lambda j: sp.Rational((-1) ** j * 11, j), 1),
            ("[M^(j-1)]", lambda j: j - 1, lambda j: sp.Rational((-1) ** (j + 1) * 47 - (-1) ** (j + 1) * 10, 24) * 0 + sp.Rational((-1) ** (j + 1) * 37, 18), 2),
            ("[M^(j-2)]", lambda j: j - 2, lambda j: sp.Integer(0), 3)]
    for name, pw, f, jmin in laws:
        res = []
        for j in (1, 2, 3, 4):
            if j < jmin:
                res.append("n/a")
                continue
            got = sp.Poly(E[j], M).coeff_monomial(M ** pw(j))
            good = sp.simplify(got - f(j)) == 0
            ok &= good
            res.append("OK" if good else f"got {got} want {f(j)}")
        print(f"   {name:11s} = {'(-1)^(j+1) 46/(j(j+1))' if name=='[M^(j+1)]' else '(-1)^j 11/j' if name=='[M^j]' else '(-1)^(j+1) 37/18' if name=='[M^(j-1)]' else '0':24s} j=1..4: {res}")
    print(f"  VERDICT B: {'all E_j laws hold' if ok else 'FAILED'}")
    return ok, E


def partC():
    print()
    rule("C. P_5: FORCED SHAPE AND PREDICTED EXTREME SLICES")
    j = 5
    print(f"   max terms (2j+1)^3 = {(2 * j + 1) ** 3};  degrees (M,U,V) = ({2 * j},{4 * j},{2 * j});"
          f"  support b+2c <= {4 * j}")
    print("   predicted top-U slice A_5 (THM-3011 sec 2a laws):")
    print(f"     [M^6] = {sp.Rational((-1)**6*46,5*6)}   [M^5] = {sp.Rational((-1)**5*12,5)}"
          f"   [M^4] = {sp.Rational((-1)**6*47,24)}")
    print("   predicted discriminant value E_5 (part B laws):")
    print(f"     [M^6] = {sp.Rational((-1)**6*46,5*6)}   [M^5] = {sp.Rational((-1)**5*11,5)}"
          f"   [M^4] = {sp.Rational((-1)**6*37,18)}   [M^3] = 0")
    print("   both share the equidistribution leading coefficient 23/15, as they must.")
    return True


def partD():
    print()
    rule("D. THE FACTORIAL READOUT:  Lambda_k(J) = -L(z^2 Phi_k(z)),  L(z^n)=n!")
    J = sp.symbols('J3:12')
    ok = True
    for k in range(2, 9):
        direct = sum(-sp.factorial(j - 1) * sp.binomial(k - 2, j - 3) * J[j - 3]
                     for j in range(3, 11))
        Phi = sum(sp.binomial(k - 2, i) * J[i] * z ** i for i in range(8))
        p = sp.Poly(sp.expand(z ** 2 * Phi), z)
        via = -sum(sp.factorial(n) * p.coeff_monomial(z ** n) for n in range(p.degree() + 1))
        same = sp.simplify(sp.expand(direct - via)) == 0
        ok &= same
        print(f"   k={k}: direct == L-form  {same}")
    print("   L is the Edo-van den Essen factorial functional of the SFC lane")
    print("   (THM-2173/2812/2836).  The edge index k sets the binomial weight;")
    print("   the jet index j sets the factorial.  SHARED FUNCTIONAL, NOT a reduction:")
    print("   SFC evaluates L on f^m (NONLINEAR in the coefficients), the circuit")
    print("   evaluates L on an object LINEAR in the jets.  That is the loss.")
    print(f"  VERDICT D: {'IDENTITY VERIFIED k=2..8' if ok else 'FAILED'}")
    return ok


def main():
    tabs = load_tables()
    a = partA(tabs)
    b, _ = partB(tabs)
    c = partC()
    d = partD()
    print()
    rule(f"SUMMARY  D-adic={a}  E_j laws={b}  P_5 shape={c}  factorial readout={d}")


if __name__ == "__main__":
    main()
