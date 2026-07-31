#!/usr/bin/env python3
"""THM-3000 -- fixed-edge cumulant-curvature universality and bounded-jet transfer.

Universe.  N(n) = sum_{i=0}^d a_i n^i with a_d > 0 and a_{d-1},...,a_{d-k-1} > 0.
Normalized Newton coefficients and ratios (THM-2997 (1)):

    h_k = a_{d-k} / (a_d * C(d,k)),      h_0 = 1,
    R_k = h_k^2 / (h_{k-1} h_{k+1}).

Formal log jets (THM-2997 (4), extended):

    log( N(n) / (a_d n^d) ) = sum_{j>=1} (ell_j / j) n^{-j},
    ell_1 = u,   ell_j = (-1)^{j-1} p_j   (p_j = j-th power sum of the roots).

Normalized jets, the scale-invariant coordinates of the problem:

    m_j := (-1)^{j-1} ell_j / d = p_j / d      (raw moment of the root measure)
    x   := m_2 / m_1^2 = -d ell_2 / u^2        (THM-2997 (7))
    z   := m_3 / m_1^3 =  d^2 ell_3 / u^3      (THM-2997 (7))
    w   := m_4 / m_1^4 = -d^3 ell_4 / u^4      (new; the third-edge coordinate)

CLAIMS PROVED / VERIFIED HERE
  A. (universality) For every fixed k >= 2,
        log(R_k / R_{k-1}) = (3x^2 - 2z - 1) d^{-2} + O(d^{-3}),
     with leading coefficient literally independent of k -- no k-dependent
     constant.  Verified symbolically for k = 2..9.
  B. (cumulant form) 3x^2 - 2z - 1 = (3 kappa_2^2 - 2 kappa_1 kappa_3)/kappa_1^4,
     kappa_j = cumulants (= central moments for j = 2,3) of the root measure.
     So the curvature is 3*variance^2 - 2*mean*skew-moment; it vanishes exactly
     on the equal-root point and on the surface 3 kappa_2^2 = 2 kappa_1 kappa_3.
  C. (mechanism) The set-partition/Moebius proof: the "-1" is the third
     difference of the falling-factorial normalization log (d)_k; the k-dependent
     pieces 3 C(k,4) and -C(k,2)^2/2 collapse to the cubic k(k-1)(3-2k)/4 whose
     third difference is the constant -3.  Checked exactly.
  D. (grading) The coefficient of d^{-r} in log(R_k/R_{k-1}) is a polynomial in
     k of degree exactly r-2, and m_j first occurs at order d^{-(j-1)} with the
     binomial factor C(k-2, j-3).  Verified for r <= 5, j <= 6.
  E. (exact third edge) closed form, not asymptotic:
        h_1 h_3^3 - h_2^3 h_4 = G_3(d,x,z,w) * u^10 / (d^10 (d-1)^4 (d-2)^3 (d-3)),
     G_3 printed in full; its leading coefficient is d^6 (3x^2 - 2z - 1).
  F. (hostile) The bounded-jet hypothesis is INDISPENSABLE and the sharp
     threshold is graded, not uniform: an explicit positive-coefficient family
     with x, z fixed and w = -alpha d has R_2 > R_1 (curvature positive) but
     R_3 < R_2 for every large d, whenever 6 alpha > 3x^2 - 2z - 1.
     Exact rational witnesses printed.

Reproduce:  python3 04-computation/gmc_fixed_edge_curvature_universality_and_bounded_jet_transfer_thm3000.py
"""

from fractions import Fraction as Fr
import sympy as sp

KMAX = 9        # edges k = 2..KMAX
ORDER = 6       # expand log(R_k/R_{k-1}) to eps^ORDER-1, eps = 1/d


# --------------------------------------------------------------------------
# 0.  symbols
# --------------------------------------------------------------------------
d = sp.symbols('d', positive=True)
eps = sp.symbols('eps', positive=True)
m = {j: sp.symbols(f'm{j}', real=True) for j in range(1, KMAX + 4)}


def falling(k):
    return sp.prod([d - i for i in range(k)])


def elementary_from_powersums(kmax):
    """Newton's identities: k E_k = sum_{i=1}^k (-1)^{i-1} E_{k-i} p_i, p_j = d m_j."""
    E = {0: sp.Integer(1)}
    for k in range(1, kmax + 1):
        E[k] = sp.expand(sum((-1) ** (i - 1) * E[k - i] * (d * m[i])
                             for i in range(1, k + 1)) / k)
    return E


# --------------------------------------------------------------------------
# A + D.  symbolic 1/d expansion of log(R_k/R_{k-1}) for every fixed edge
# --------------------------------------------------------------------------
def edge_expansions():
    E = elementary_from_powersums(KMAX + 2)
    logh = {0: sp.Integer(0)}
    for k in range(1, KMAX + 2):
        hk = sp.cancel(sp.factorial(k) * E[k] / falling(k))   # = E_k / C(d,k)
        ratio = sp.cancel(sp.cancel(hk / m[1] ** k).subs(d, 1 / eps))
        logh[k] = sp.expand(sp.series(sp.log(ratio), eps, 0, ORDER).removeO())
    out = {}
    for k in range(2, KMAX + 1):
        # log(R_k/R_{k-1}) = -Delta^3 (log h)_{k-2}
        Dk = sp.expand(-logh[k + 1] + 3 * logh[k] - 3 * logh[k - 1] + logh[k - 2])
        out[k] = {r: sp.simplify(sp.expand(Dk).coeff(eps, r)) for r in range(0, ORDER)}
    return out


def curvature_symbol():
    return sp.simplify((3 * m[2] ** 2 - 2 * m[1] * m[3] - m[1] ** 4) / m[1] ** 4)


def report_A_D(exp):
    print("=" * 74)
    print("A. UNIVERSALITY OF THE LEADING d^-2 CURVATURE ACROSS FIXED EDGES")
    print("=" * 74)
    curv = curvature_symbol()
    print("  curvature C(x,z) = 3x^2 - 2z - 1 =", sp.simplify(curv))
    ok = True
    for k in range(2, KMAX + 1):
        c0 = sp.simplify(exp[k][0])
        c1 = sp.simplify(exp[k][1])
        c2 = sp.simplify(exp[k][2])
        same = sp.simplify(c2 - curv) == 0
        ok &= (c0 == 0) and (c1 == 0) and same
        print(f"   k={k}: [d^0]=0? {c0 == 0}   [d^-1]=0? {c1 == 0}   "
              f"[d^-2]==curvature? {same}")
    print("  VERDICT A:", "UNIVERSAL (k-independent, constant 1)" if ok else "FAILED")

    print()
    print("=" * 74)
    print("D. GRADING: k-degree of the d^-r coefficient, and first occurrence of m_j")
    print("=" * 74)
    kk = sp.symbols('kk')
    for r in range(2, ORDER):
        # interpolate the coefficient as a polynomial in k
        pts = [(k, sp.expand(exp[k][r])) for k in range(2, KMAX + 1)]
        # degree detection by finite differences on the vector of coefficients
        vals = [p[1] for p in pts]
        deg = -1
        cur = vals
        while any(sp.simplify(v) != 0 for v in cur):
            deg += 1
            cur = [sp.expand(cur[i + 1] - cur[i]) for i in range(len(cur) - 1)]
            if not cur:
                break
        print(f"   d^-{r}: k-degree = {deg}   (predicted r-2 = {r - 2})  "
              f"{'OK' if deg == r - 2 else 'MISMATCH'}")
    print()
    for j in range(3, KMAX + 2):
        first = None
        for r in range(2, ORDER):
            if any(sp.expand(exp[k][r]).has(m[j]) for k in range(2, KMAX + 1)):
                first = r
                break
        pred = j - 1
        tag = "OK" if first == pred else ("beyond window" if first is None else "MISMATCH")
        print(f"   m_{j}: first order d^-{first}   (predicted d^-{pred})  {tag}")
        if first is not None and first < ORDER:
            # binomial shape of the first-occurrence coefficient
            co = [sp.simplify(sp.expand(exp[k][first]).coeff(m[j], 1)) for k in range(2, KMAX + 1)]
            base = None
            shape_ok = True
            for idx, k in enumerate(range(2, KMAX + 1)):
                b = sp.binomial(k - 2, j - 3)
                if b == 0:
                    shape_ok &= (sp.simplify(co[idx]) == 0)
                else:
                    q = sp.simplify(co[idx] / b)
                    if base is None:
                        base = q
                    shape_ok &= (sp.simplify(q - base) == 0)
            print(f"        coefficient = {sp.simplify(base)} * C(k-2,{j - 3})  "
                  f"-> {'OK' if shape_ok else 'MISMATCH'}")
    return ok


# --------------------------------------------------------------------------
# B.  cumulant identity
# --------------------------------------------------------------------------
def report_B():
    print()
    print("=" * 74)
    print("B. CUMULANT FORM OF THE CURVATURE")
    print("=" * 74)
    k1, k2, k3 = sp.symbols('kappa1 kappa2 kappa3')
    sub = {m[1]: k1, m[2]: k2 + k1 ** 2, m[3]: k3 + 3 * k1 * k2 + k1 ** 3}
    lhs = sp.simplify(curvature_symbol().subs(sub))
    rhs = sp.simplify((3 * k2 ** 2 - 2 * k1 * k3) / k1 ** 4)
    print("  3x^2-2z-1  in cumulants  =", sp.factor(lhs))
    print("  target (3k2^2-2k1k3)/k1^4 =", sp.factor(rhs))
    ok = sp.simplify(lhs - rhs) == 0
    print("  VERDICT B:", "IDENTITY HOLDS" if ok else "FAILED")
    print("  => curvature > 0  <=>  3 kappa_2^2 > 2 kappa_1 kappa_3.")
    print("     equal roots (kappa_2=kappa_3=0) is the exact zero of the curvature;")
    print("     symmetric root measures (kappa_3=0) give curvature 3 kappa_2^2/kappa_1^4 >= 0.")
    return ok


# --------------------------------------------------------------------------
# C.  set-partition mechanism (the human-checkable proof)
# --------------------------------------------------------------------------
def report_C():
    print()
    print("=" * 74)
    print("C. MECHANISM: WHERE EACH OF 3x^2, -2z, -1 COMES FROM")
    print("=" * 74)
    k = sp.symbols('k')

    def D3(f):
        return sp.expand(f.subs(k, k + 3) - 3 * f.subs(k, k + 2)
                         + 3 * f.subs(k, k + 1) - f)

    # falling factorial (d)_k = d^k prod_{i<k}(1-i/d)
    S2 = sp.expand(sp.summation(sp.symbols('i') ** 2, (sp.symbols('i'), 0, k - 1)))
    d3S2 = sp.simplify(D3(S2))
    print("  log (d)_k = k log d - C(k,2)/d - S_2(k)/(2 d^2) - ... ,  S_2(k)=sum_{i<k} i^2")
    print("    Delta^3 [ -C(k,2) ] =", sp.simplify(D3(-sp.binomial(k, 2))), " (quadratic -> 0)")
    print("    Delta^3 S_2(k) =", d3S2, " => +Delta^3 log(d)_k contributes",
          sp.simplify(-d3S2 / 2), "* d^-2   <-- THE '-1'")

    # k! e_k = d^k [ 1 - C(k,2) m2/d + (3C(k,4) m2^2 + 2C(k,3) m3)/d^2 + ... ]
    # log of that: 1/d^2 bracket = 3C(k,4)m2^2 + 2C(k,3)m3 - C(k,2)^2 m2^2/2
    collapse = sp.expand(3 * sp.binomial(k, 4) - sp.binomial(k, 2) ** 2 / 2)
    print()
    print("  k! e_k = d^k [ 1 - C(k,2) m2/d + (3C(k,4) m2^2 + 2C(k,3) m3)/d^2 + O(d^-3) ]")
    print("    (partitions of [k]: singletons -> d^k; one doubleton -> -C(k,2) m2 d^{k-1};")
    print("     two doubletons -> +3C(k,4) m2^2 d^{k-2}; one tripleton -> +2C(k,3) m3 d^{k-2})")
    print("    log-correction merges 3C(k,4) with -C(k,2)^2/2 :")
    print("      3C(k,4) - C(k,2)^2/2 =", sp.factor(collapse), "= k(k-1)(3-2k)/4  (CUBIC)")
    print("    Delta^3 of that cubic =", sp.simplify(D3(collapse)), " => -(-3) m2^2 = +3 m2^2")
    print("    Delta^3 [ 2C(k,3) ] =", sp.simplify(D3(2 * sp.binomial(k, 3))),
          " => -(2) m3 = -2 m3")
    ok = (sp.simplify(collapse - k * (k - 1) * (3 - 2 * k) / 4) == 0
          and sp.simplify(D3(collapse) + 3) == 0
          and sp.simplify(D3(2 * sp.binomial(k, 3)) - 2) == 0
          and sp.simplify(d3S2 - 2) == 0)
    print("  => leading coefficient  3 m2^2 - 2 m3 - 1  (m1 = 1),  k CANCELS EXACTLY.")
    print("  VERDICT C:", "MECHANISM CONFIRMED" if ok else "FAILED")
    return ok


# --------------------------------------------------------------------------
# E.  exact second- and third-edge closed forms
# --------------------------------------------------------------------------
def exact_edges():
    u = sp.symbols('u', positive=True)
    x, z, w = sp.symbols('x z w')
    l2, l3, l4 = sp.symbols('ell2 ell3 ell4')
    e = {1: u,
         2: (u ** 2 + l2) / 2,
         3: (u ** 3 + 3 * u * l2 + 2 * l3) / 6,
         4: (u ** 4 + 6 * u ** 2 * l2 + 3 * l2 ** 2 + 8 * u * l3 + 6 * l4) / 24}
    h = {kk: sp.factorial(kk) * e[kk] / falling(kk) for kk in e}
    sub = {l2: -x * u ** 2 / d, l3: z * u ** 3 / d ** 2, l4: -w * u ** 4 / d ** 3}

    print()
    print("=" * 74)
    print("E. EXACT CLOSED FORMS (no asymptotics)")
    print("=" * 74)

    F2 = sp.cancel(sp.together((h[2] ** 3 - h[1] ** 3 * h[3]).subs(sub)))
    n2, d2 = sp.fraction(F2)
    G2 = sp.expand(sp.cancel(n2 / u ** 6))
    thm2997 = sp.expand((d - 2) * (d - x) ** 3 - (d - 1) ** 2 * (d ** 2 - 3 * d * x + 2 * z))
    ok2 = sp.simplify(G2 - thm2997) == 0
    print("  second edge:  h2^3 - h1^3 h3 = G_2 * u^6 / (%s)" % sp.factor(d2 / u ** 0))
    print("    G_2 =", sp.collect(G2, d))
    print("    matches THM-2997 (8)? ", ok2)

    F3 = sp.cancel(sp.together((h[1] * h[3] ** 3 - h[2] ** 3 * h[4]).subs(sub)))
    n3, d3f = sp.fraction(F3)
    G3 = sp.expand(sp.cancel(n3 / u ** 10))
    print()
    print("  third edge:   h1 h3^3 - h2^3 h4 = G_3 * u^10 / (%s)" % sp.factor(d3f))
    P = sp.Poly(G3, d)
    for pw in range(P.degree(), -1, -1):
        print(f"    [d^{pw}] {sp.factor(sp.simplify(G3.coeff(d, pw)))}")
    lead_ok = sp.simplify(G3.coeff(d, 6) - (3 * x ** 2 - 2 * z - 1)) == 0
    w_ok = sp.simplify(sp.expand(G3.coeff(d, 5)).coeff(w, 1) - 6) == 0
    print("    leading [d^6] == 3x^2-2z-1 ? ", lead_ok)
    print("    w enters first at [d^5] with coefficient 6 ? ", w_ok,
          "   (= 6*C(k-2,1) at k=3, matching D)")
    return ok2 and lead_ok and w_ok, G3, (x, z, w, u)


# --------------------------------------------------------------------------
# F.  the hostile: unbounded fourth jet flips the third edge only
# --------------------------------------------------------------------------
def exact_ratios(coeffs, dd):
    """coeffs = {0: a_d, 1: a_{d-1}, ...} as Fractions; returns R_1,R_2,R_3."""
    h = {}
    for kk, c in coeffs.items():
        h[kk] = Fr(c, 1) / (Fr(sp.binomial(dd, kk)) * coeffs[0]) if kk else Fr(1)
    R = {}
    for kk in range(1, max(coeffs) ):
        R[kk] = h[kk] ** 2 / (h[kk - 1] * h[kk + 1])
    return R, h


def hostile(alpha_num=1, alpha_den=1, xval=Fr(1), zval=Fr(1, 2), ds=(200, 400, 800, 1600, 3200)):
    """Family with m1=1 (u=d), x,z fixed, w = -alpha*d.

    Coefficients are forced by the jets:
      e1 = u,  e2=(u^2+l2)/2, e3=(u^3+3u l2+2 l3)/6,
      e4=(u^4+6u^2 l2+3 l2^2+8u l3+6 l4)/24,
      l2 = -x u^2/d,  l3 = z u^3/d^2,  l4 = -w u^4/d^3,  w = -alpha d.
    Only a_d..a_{d-4} matter for R_1,R_2,R_3; all are positive here.
    """
    alpha = Fr(alpha_num, alpha_den)
    print()
    print("=" * 74)
    print("F. HOSTILE: unbounded fourth normalized jet")
    print("=" * 74)
    curv = 3 * xval ** 2 - 2 * zval - 1
    print(f"  x = {xval}, z = {zval}  ->  curvature 3x^2-2z-1 = {curv} "
          f"({'>0' if curv > 0 else '<=0'})")
    print(f"  w_d = -alpha*d with alpha = {alpha};  predicted third-edge sign flip "
          f"iff 6*alpha > curvature  ->  6*alpha = {6 * alpha}")
    print()
    print("   d      R_1<R_2 ?         R_2<R_3 ?      d^2*log(R2/R1)   d^2*log(R3/R2)")
    import math
    all_flip = True
    for dd in ds:
        u = Fr(dd)                      # m1 = 1
        l2 = -xval * u ** 2 / dd
        l3 = zval * u ** 3 / dd ** 2
        wv = -alpha * dd
        l4 = -wv * u ** 4 / dd ** 3
        e1 = u
        e2 = (u ** 2 + l2) / 2
        e3 = (u ** 3 + 3 * u * l2 + 2 * l3) / 6
        e4 = (u ** 4 + 6 * u ** 2 * l2 + 3 * l2 ** 2 + 8 * u * l3 + 6 * l4) / 24
        assert e1 > 0 and e2 > 0 and e3 > 0 and e4 > 0, (dd, e1, e2, e3, e4)
        coeffs = {0: Fr(1), 1: e1, 2: e2, 3: e3, 4: e4}   # a_{d-k}/a_d
        R, h = exact_ratios(coeffs, dd)
        up12 = R[2] > R[1]
        up23 = R[3] > R[2]
        s12 = dd ** 2 * math.log(float(R[2] / R[1]))
        s23 = dd ** 2 * math.log(float(R[3] / R[2]))
        all_flip &= (up12 and not up23)
        print(f"  {dd:5d}   {str(up12):5s} (edge2)    {str(up23):5s} (edge3)   "
              f"{s12:14.6f}   {s23:14.6f}")
    print()
    print("  VERDICT F:", "HOSTILE CONFIRMED -- edge 2 rises, edge 3 falls, all d"
          if all_flip else "NOT A HOSTILE at these d")
    print("  => 'bounded higher log jets' is INDISPENSABLE, and the sharp condition")
    print("     is GRADED: m_j / m_1^j = o(d^{j-3}) for 3 <= j <= k+1, not uniform")
    print("     boundedness.  At j=4 the threshold w ~ d is attained exactly.")
    return all_flip


# --------------------------------------------------------------------------
# positive control: bounded jets -> universality actually predicts the sign
# --------------------------------------------------------------------------
def positive_control():
    print()
    print("=" * 74)
    print("POSITIVE CONTROL: real-rooted families with bounded jets")
    print("=" * 74)
    import math
    print("  N(n) = prod_{i=1}^{d} (n + r_i);  curvature predicted from the roots.")
    cases = {
        "two-point 1,2 (half/half)": lambda dd: [Fr(1)] * (dd // 2) + [Fr(2)] * (dd - dd // 2),
        "arithmetic 1..d rescaled": lambda dd: [Fr(i, dd) for i in range(1, dd + 1)],
        "geometric 2^{-i}, i<d": lambda dd: [Fr(1, 2 ** (i % 8)) for i in range(dd)],
        "heavy tail one big root": lambda dd: [Fr(1)] * (dd - 1) + [Fr(dd)],
    }
    ok = True
    for name, gen in cases.items():
        print(f"\n  -- {name}")
        for dd in (60, 120, 240):
            r = gen(dd)
            p1 = sum(r); p2 = sum(t ** 2 for t in r); p3 = sum(t ** 3 for t in r)
            m1 = p1 / dd; m2 = p2 / dd; m3 = p3 / dd
            curv = 3 * m2 ** 2 / m1 ** 4 - 2 * m3 / m1 ** 3 - 1
            # exact elementary symmetric e_1..e_5
            E = [Fr(1)] + [Fr(0)] * 5
            for t in r:
                for kk in range(5, 0, -1):
                    E[kk] = E[kk] + t * E[kk - 1]
            coeffs = {kk: E[kk] for kk in range(5)}
            R, _ = exact_ratios(coeffs, dd)
            s12 = dd ** 2 * math.log(float(R[2] / R[1]))
            s23 = dd ** 2 * math.log(float(R[3] / R[2]))
            match = (curv > 0) == (R[2] > R[1]) and (curv > 0) == (R[3] > R[2])
            ok &= match
            print(f"     d={dd:4d}  curvature={float(curv):11.6f}  "
                  f"d^2 log(R2/R1)={s12:11.6f}  d^2 log(R3/R2)={s23:11.6f}  "
                  f"signs match: {match}")
    print("\n  VERDICT (positive control):",
          "sign of both edges tracks the single universal curvature" if ok else "MISMATCH")
    return ok


# --------------------------------------------------------------------------
# G.  the invariant (jet) form and its weight grading
# --------------------------------------------------------------------------
def report_G(exp):
    """Rewrite the expansion in the NORMALIZED LOG JETS  J_j = ell_j / u^j.

    m_1 = u/d and m_j = (-1)^{j-1} p_j / d = (-1)^{j-1} ell_j, so with m_1 = 1
    (u = d) one has  m_j = (-1)^{j-1} d^{j-1} J_j.
    Grade by  wt(J_j) = j-1,  wt(1/d) = 1.  Then the eps^r coefficient of the
    expansion is exactly the weight-r part, and it is d-free in jet variables.
    """
    print()
    print("=" * 74)
    print("G. INVARIANT JET FORM AND WEIGHT GRADING  (wt J_j = j-1, wt 1/d = 1)")
    print("=" * 74)
    J = {j: sp.symbols(f'J{j}') for j in range(2, KMAX + 4)}
    sub = {m[1]: sp.Integer(1)}
    for j in range(2, KMAX + 4):
        sub[m[j]] = (-1) ** (j - 1) * d ** (j - 1) * J[j]
    t = sp.symbols('t', positive=True)
    ok = True
    for r in range(2, ORDER):
        forms = []
        for k in range(2, KMAX + 1):
            # eps^r * (coefficient) rewritten in jets; restore the d^-r factor
            wr = sp.expand(sp.simplify(sp.expand(exp[k][r]).subs(sub) / d ** r))
            forms.append((k, wr))
        # weight-homogeneity: J_j -> t^{j-1} J_j and 1/d -> t/d must scale by t^r
        def scaled(f):
            s = f.subs({J[j]: t ** (j - 1) * J[j] for j in J}, simultaneous=True)
            return sp.expand(s.subs(d, d / t))
        homog = all(sp.simplify(scaled(f) - t ** r * f) == 0 for _, f in forms)
        base = sp.expand(forms[0][1])
        kfree = all(sp.simplify(sp.expand(f) - base) == 0 for _, f in forms)
        print(f"   weight {r}: weight-{r}-homogeneous? {homog}   k-free? {kfree}")
        if r == 2:
            print("     W_2 =", sp.expand(base), "   <-- k-INDEPENDENT UNIVERSAL PART")
            ok &= homog and kfree
        else:
            # k-degree of the weight-r part
            vals = [sp.expand(f) for _, f in forms]
            deg = -1
            cur = vals
            while any(sp.simplify(v) != 0 for v in cur):
                deg += 1
                cur = [sp.expand(cur[i + 1] - cur[i]) for i in range(len(cur) - 1)]
                if not cur:
                    break
            print(f"     k-degree {deg} (predicted {r - 2})",
                  "OK" if deg == r - 2 else "MISMATCH")
            ok &= homog and (deg == r - 2)
            print("     W_%d(k=3) =" % r, sp.expand(vals[1]))
    print()
    print("  first-occurrence law:  [weight j-1] contains  -(j-1)! * C(k-2, j-3) * J_j")
    for j in range(3, 7):
        r = j - 1
        if r >= ORDER:
            continue
        coeffs = []
        for k in range(2, KMAX + 1):
            wr = sp.expand(sp.simplify(sp.expand(exp[k][r]).subs(sub) / d ** r))
            coeffs.append(sp.simplify(sp.expand(wr).coeff(J[j], 1)))
        pred = [sp.simplify(-sp.factorial(j - 1) * sp.binomial(k - 2, j - 3))
                for k in range(2, KMAX + 1)]
        good = all(sp.simplify(a - b) == 0 for a, b in zip(coeffs, pred))
        print(f"   J_{j}: coefficients {coeffs}  vs  -({j - 1})!*C(k-2,{j - 3}) {pred}"
              f"  -> {'OK' if good else 'MISMATCH'}")
        ok &= good
    print("  VERDICT G:", "INVARIANT FORM CONFIRMED" if ok else "FAILED")
    return ok


def main():
    exp = edge_expansions()
    a = report_A_D(exp)
    b = report_B()
    c = report_C()
    g = report_G(exp)
    e, _, _ = exact_edges()
    p = positive_control()
    # hostile: curvature +1/2 > 0 but 6*alpha = 3 > 1/2
    f = hostile(alpha_num=1, alpha_den=2, xval=Fr(1), zval=Fr(1, 4))
    print()
    print("=" * 74)
    print("SUMMARY  A=%s  B=%s  C=%s  D/G=%s  E=%s  positive-control=%s  F=%s"
          % (a, b, c, g, e, p, f))
    print("=" * 74)


if __name__ == "__main__":
    main()
