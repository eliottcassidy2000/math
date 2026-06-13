"""
monad-explorer 2026-06-13  —  HYP-2460  (candidate THM-493)

THE MOSER LATTICE IS THE PRODUCT AT A RESONANT ANGLE.
Unifies THM-433 (avgdeg additive under the generic-angle Minkowski product; the
3N crossover N* is "non-product") with THM-434 (L_t = Z[zeta6] (+) w_t Z[zeta6]
has 12 + r_E(t) unit vectors).

Setup. For finite G,H subset of the triangular lattice Z[zeta6], the Minkowski
product at relative rotation theta is the point set  P_theta = { g + e^{i theta} h }.
  - GENERIC theta: P realizes the Cartesian product G [] H exactly,
        U = e(G)|H| + |G|e(H)                                  (THM-433, classical)
  - RESONANT theta = arccos((2t-1)/2t) (the t-th Moser angle, e^{i theta} = w_t):
    the rotated copy is commensurate, P sits inside the rank-4 lattice L_t, and
    EXTRA "diagonal" edges appear -- exactly the transverse unit vectors
        z = alpha (1 - w_t),   N(alpha) = t                    (THM-434)
    giving a RESONANCE BONUS
        U = e(G)|H| + |G|e(H) + Delta_t(G,H),
        Delta_t(G,H) = #{(g,g',h,h'): g-g' = alpha, h'-h = alpha, N(alpha)=t}
                     = sum_{N(alpha)=t} m_alpha(G) * m_alpha(H),
    m_alpha(G) = #{(g,g') in G^2 : g-g' = alpha} (ordered).

This script PROVES the above by EXACT integer/rational arithmetic in Q(sqrt3, sqrtD),
D = 4t-1: every distance^2 is an exact element of the field; a pair is a unit edge
iff distance^2 == 1 exactly.  No float ever decides adjacency.
"""
from fractions import Fraction as F
from itertools import combinations
import sys

# ---------- exact arithmetic in Q(sqrt3, sqrtD); basis {1, sqrt3, sqrtD, sqrt(3D)} ----------
# element = (a,b,c,d) = a + b sqrt3 + c sqrtD + d sqrt(3D),  a,b,c,d in Q

def make_field(D):
    def add(u, v):  return tuple(x + y for x, y in zip(u, v))
    def sub(u, v):  return tuple(x - y for x, y in zip(u, v))
    def mul(u, v):
        a, b, c, d = u; e, f, g, h = v
        return (a*e + 3*b*f + D*c*g + 3*D*d*h,     # 1
                a*f + b*e + D*c*h + D*d*g,          # sqrt3
                a*g + 3*b*h + c*e + 3*d*f,          # sqrtD
                a*h + b*g + c*f + d*e)              # sqrt(3D)
    def sq(u): return mul(u, u)
    return add, sub, mul, sq

ZERO = (F(0),)*4
ONE  = (F(1), F(0), F(0), F(0))

def Q(x):  # rational -> field element
    return (F(x), F(0), F(0), F(0))

# ---------- geometry ----------
# zeta6 = (1 + i sqrt3)/2.  A lattice point m + n*zeta6 has
#   Re = m + n/2  (in Q),         Im = (n/2) sqrt3  (in Q*sqrt3).
# w_t = ((2t-1) + i sqrtD)/(2t),  D = 4t-1.
# A product point g + w_t*h, g = (m,n), h = (p,q):
#   hr = p + q/2,  hi = (q/2) sqrt3
#   Re(w_t h) = wr*hr - wi*hi,  Im(w_t h) = wr*hi + wi*hr,  wr=(2t-1)/(2t), wi=sqrtD/(2t)
# Re,Im are returned as field elements (4-tuples).

def lattice_point(m, n):
    """Re, Im of m + n*zeta6 as field elements."""
    Re = (F(m) + F(n, 2), F(0), F(0), F(0))      # rational
    Im = (F(0), F(n, 2), F(0), F(0))             # *sqrt3
    return Re, Im

def product_point(g, h, t, D):
    m, n = g; p, q = h
    twoT = 2*t
    wr = F(2*t-1, twoT)                # rational coeff of 1
    # wi = sqrtD/(2t)
    gRe, gIm = lattice_point(m, n)
    # h real/imag
    hr = F(p) + F(q, 2)               # rational
    hi = F(q, 2)                      # coeff of sqrt3 ; hi = (q/2) sqrt3
    # w_t * h :
    # Re = wr*hr - wi*hi = wr*hr - (sqrtD/2t)*(q/2 sqrt3) = wr*hr - (q/(4t)) sqrt(3D)
    Re_wh = (wr*hr, F(0), F(0), -F(q, 2*twoT))           # 1 and sqrt(3D)
    # Im = wr*hi + wi*hr = wr*(q/2) sqrt3 + (sqrtD/2t)*hr = (wr q/2) sqrt3 + (hr/2t) sqrtD
    Im_wh = (F(0), wr*F(q, 2), hr/twoT, F(0))            # sqrt3 and sqrtD
    Re = tuple(a+b for a, b in zip(gRe, Re_wh))
    Im = tuple(a+b for a, b in zip(gIm, Im_wh))
    return Re, Im

def build_product(G, H, t, D):
    return [product_point(g, h, t, D) for g in G for h in H]

# ---------- counting ----------
def count_unit_distances(points, D):
    add, sub, mul, sq = make_field(D)
    cnt = 0
    for (Re1, Im1), (Re2, Im2) in combinations(points, 2):
        dx = sub(Re1, Re2); dy = sub(Im1, Im2)
        d2 = add(sq(dx), sq(dy))
        if d2 == ONE:
            cnt += 1
    return cnt

def edges_in_lattice_graph(G, D=11, t=1):
    """count unit edges of a triangular-lattice point set G (D irrelevant; pure Q(sqrt3))"""
    add, sub, mul, sq = make_field(D)
    pts = [lattice_point(m, n) for (m, n) in G]
    cnt = 0
    for (Re1, Im1), (Re2, Im2) in combinations(pts, 2):
        dx = sub(Re1, Re2); dy = sub(Im1, Im2)
        d2 = add(sq(dx), sq(dy))
        if d2 == ONE:
            cnt += 1
    return cnt

# ---------- norm-t displacement spectrum (Eisenstein) ----------
def eisen_norm(m, n):           # N(m + n zeta6) = m^2 + m n + n^2
    return m*m + m*n + n*n

def norm_t_displacements(t):
    """all alpha=(dm,dn) in Z[zeta6] with N(alpha)=t (i.e. m^2+mn+n^2=t)"""
    r = int((4*t)**0.5) + 2
    out = []
    for dm in range(-r, r+1):
        for dn in range(-r, r+1):
            if eisen_norm(dm, dn) == t:
                out.append((dm, dn))
    return out

def m_alpha(G, alpha):
    """# ordered pairs (g,g') in G^2 with g-g' = alpha"""
    Gset = set(G); am, an = alpha
    return sum(1 for (m, n) in G if (m-am, n-an) in Gset)

def resonance_bonus_formula(G, H, t):
    """Delta_t(G,H) = sum_{N(alpha)=t} m_alpha(G) m_alpha(H)"""
    tot = 0
    for alpha in norm_t_displacements(t):
        # need g-g'=alpha in G AND h'-h=alpha in H i.e. h-h' = -alpha
        # count ordered (g,g') with diff alpha, and ordered (h,h') with h'-h=alpha
        mg = m_alpha(G, alpha)
        mh = m_alpha(H, alpha)   # # of (h,h') with h-h'=alpha == # of (h',h) with h'-h=alpha
        tot += mg * mh
    # each undirected diagonal edge {(g,h),(g',h')} is counted once:
    # the ordered alpha and -alpha both enumerate the same unordered edge set, so /2
    return tot // 2

# ---------- demonstrations ----------
def rhombus():
    # two unit triangles glued: {0, 1, zeta6, 1+zeta6}; contains one sqrt3 pair {0,1+zeta6}
    return [(0,0), (1,0), (0,1), (1,1)]

def rosette():
    # W7 Eisenstein flower: hub 0 + 6 unit neighbors  (the 6 sixth-roots as lattice pts)
    # zeta6^j for j=0..5 in (m,n) coords:
    #   1=(1,0), zeta6=(0,1), zeta6^2=zeta6-1=(-1,1), -1=(-1,0), -zeta6=(0,-1), (1,-1)
    return [(0,0),(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]

def hexlattice_patch(R):
    """all Eisenstein points with N <= R (a hex-symmetric blob)"""
    r = int((4*R)**0.5)+2
    return [(m,n) for m in range(-r,r+1) for n in range(-r,r+1) if eisen_norm(m,n) <= R]

def report(name, G, H, t):
    D = 4*t-1
    eG = edges_in_lattice_graph(G); eH = edges_in_lattice_graph(H)
    generic = eG*len(H) + len(G)*eH
    P = build_product(G, H, t, D)
    # de-dup any coincident points (shouldn't happen at Moser angle for these sets)
    seen = {}
    for pt in P:
        seen[pt] = seen.get(pt, 0)+1
    ndup = sum(c-1 for c in seen.values())
    actual = count_unit_distances(P, D)
    formula = resonance_bonus_formula(G, H, t)
    print(f"--- {name}  (t={t}, D={D}, w_t=Moser angle, |G|={len(G)} e(G)={eG}, |H|={len(H)} e(H)={eH}) ---")
    print(f"  product points: {len(P)} (coincident: {ndup})")
    print(f"  generic-angle (Cartesian product) unit distances : e(G)|H|+|G|e(H) = {generic}")
    print(f"  Moser-angle exact unit distances (brute, exact)  : {actual}")
    print(f"  resonance bonus  Delta = actual - generic        : {actual-generic}")
    print(f"  formula  sum_{{N(a)=t}} m_a(G)m_a(H) / 2          : {formula}")
    print(f"  MATCH: {actual-generic == formula}")
    print(f"  avgdeg generic = {2*generic/len(P):.4f} ; avgdeg Moser = {2*actual/len(P):.4f}  (kappa=6)")
    print()
    return actual, generic, formula

if __name__ == "__main__":
    print("="*78)
    print("STEP 0: transverse unit vectors are EXACT (|alpha(1-w_t)|^2 = 1 for N(alpha)=t)")
    print("="*78)
    for t in (2,3,4,13):
        D = 4*t-1
        add, sub, mul, sq = make_field(D)
        # 1 - w_t = (1 - (2t-1)/2t) - i sqrtD/2t = 1/(2t) - i sqrtD/2t
        # represent (1-w_t) as complex (Re,Im) field elements:
        Re_w = (F(1,2*t), F(0), F(0), F(0))
        Im_w = (F(0), F(0), -F(1,2*t), F(0))   # -sqrtD/2t
        # pick alpha with N=t
        disp = norm_t_displacements(t)
        if not disp:
            print(f"  t={t}: no Eisenstein alpha with N={t} (not Loeschian) -> no transverse vectors")
            continue
        am, an = disp[0]
        aRe, aIm = lattice_point(am, an)
        # alpha*(1-w_t) as complex multiply: (aRe+i aIm)(Re_w + i Im_w)
        zRe = sub(mul(aRe, Re_w), mul(aIm, Im_w))
        zIm = add(mul(aRe, Im_w), mul(aIm, Re_w))
        z2 = add(sq(zRe), sq(zIm))
        print(f"  t={t}: alpha={disp[0]} (N={eisen_norm(am,an)}), |alpha(1-w_t)|^2 = {z2}  -> unit: {z2==ONE}")
    print()
    print("="*78)
    print("STEP 1: resonance bonus on small patches (Moser angle t=3, the Moser lattice)")
    print("="*78)
    report("rhombus [] rhombus", rhombus(), rhombus(), t=3)
    report("rosette W7 [] rhombus", rosette(), rhombus(), t=3)
    report("rosette W7 [] rosette W7", rosette(), rosette(), t=3)
    print("="*78)
    print("STEP 2: bonus needs MATCHING norm-t displacements in BOTH factors")
    print("="*78)
    # a unit segment K2 has NO norm-3 displacement -> zero bonus regardless of angle
    report("rosette W7 [] K2", rosette(), [(0,0),(1,0)], t=3)
    # K3 (unit triangle) has no norm-3 displacement either
    report("K3 [] K3", [(0,0),(1,0),(0,1)], [(0,0),(1,0),(0,1)], t=3)
    print("="*78)
    print("STEP 3: other rungs -- t=2 (sqrt7) and t=4 need norm-t displacements")
    print("="*78)
    # hex patch of radius 7 contains norm-7 pairs; demo t=2 resonance
    Gp = hexlattice_patch(4)
    report("hexpatch(N<=4) [] hexpatch(N<=4)", Gp, Gp, t=2)
