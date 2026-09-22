#!/usr/bin/env python3
"""Lane berggren_edge_transport (session collatz-mod6-20260917 wave 2026-09-21, mac-mini).

A directed edge x -> y of the family E_{a,b,k} = {(x,y) odd, coprime, distinct:
a x + b = 2^k y} is encoded (odd_square_20260921_edges.md, THM-3756) by its root
pair (s,t) = (max(x,y), min(x,y)), i.e. the primitive triple
(xy, |x^2-y^2|/2, (x^2+y^2)/2).  This script proves and checks how the three
Berggren children of that root pair transport the guard (a,b,k), classifies
exactly which children remain Collatz edges (multiplier 3, parameter +-1),
gives the inverse (parent) table on the three THM-3756 cones, identifies the
inverse-fibre braid R(x)=4x+b as a Berggren word, and runs the finite censuses.

Every probe prints numbered claims S1, S2, ... with exactly one of
PROVED / FINITE-EXACT / CITED / HEURISTIC / REFUTED / SCOPE / OPEN.
Explicit `raise` only (survives python3 -O).  RAM < 1 GB, runtime ~1-2 min.

Reproduce:
  python3 04-computation/experiments/collatz_mod6_20260921_berggren_edge_transport.py \
      > 05-knowledge/results/collatz_mod6_20260921_berggren_edge_transport.out
"""
import math
import sys
from collections import defaultdict
from fractions import Fraction

import numpy as np
import sympy as sp


def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)


def is_pow2(q):
    return q >= 1 and (q & (q - 1)) == 0


def v2(n):
    return (n & -n).bit_length() - 1


def banner(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


# Berggren children in root coordinates (s,t), s>t>=1 odd coprime.
def B1(s, t):
    return (s + 2 * t, t)


def B2(s, t):
    return (2 * s + t, s)


def B3(s, t):
    return (2 * s - t, s)


CHILDREN = (("B1", B1), ("B2", B2), ("B3", B3))


def roots(x, y):
    return (max(x, y), min(x, y))


def edge_reading(p, q, a, b):
    """All k>=1 with a*p+b = 2^k q  (edge p -> q of E_{a,b,k}); returns k or None."""
    v = a * p + b
    if v <= 0 or v % q:
        return None
    r = v // q
    if r >= 2 and is_pow2(r):
        return v2(r)
    return None


def legal_readings(pair, a, bs):
    """Readings of the unordered pair as an E_{a,b,k} edge, either orientation."""
    s, t = pair
    out = []
    for b in bs:
        k = edge_reading(s, t, a, b)
        if k is not None:
            out.append((s, t, a, b, k))
        k = edge_reading(t, s, a, b)
        if k is not None:
            out.append((t, s, a, b, k))
    return out


banner("Lane berggren_edge_transport: Berggren children of Collatz-family edges")

# ---------------------------------------------------------------- S1
print("\n### S1  PROVED: Berggren children in root coordinates (s,t)=(m+n,m-n)")
m, n = sp.symbols("m n")
s_, t_ = m + n, m - n
euclid = {"B1": (2 * m - n, m), "B2": (2 * m + n, m), "B3": (m + 2 * n, n)}
root_maps = {"B1": (s_ + 2 * t_, t_), "B2": (2 * s_ + t_, s_), "B3": (2 * s_ - t_, s_)}
for name in ("B1", "B2", "B3"):
    mm, nn = euclid[name]
    lhs = (sp.expand(mm + nn), sp.expand(mm - nn))
    rhs = tuple(sp.expand(e) for e in root_maps[name])
    check(lhs == rhs, "root form of " + name)
    print("  %s: Euclid (m,n)->%s  ==  roots (s,t)->%s" % (name, euclid[name], root_maps[name]))
mats = {"B1": sp.Matrix([[1, 2], [0, 1]]), "B2": sp.Matrix([[2, 1], [1, 0]]),
        "B3": sp.Matrix([[2, -1], [1, 0]])}
print("  matrices on (s,t): B1=[[1,2],[0,1]] det %d, B2=[[2,1],[1,0]] det %d, B3=[[2,-1],[1,0]] det %d"
      % (mats["B1"].det(), mats["B2"].det(), mats["B3"].det()))
# numeric: all Euclid pairs m<=60 map to the same triples
cnt = 0
for mm in range(2, 61):
    for nn in range(1, mm):
        if (mm - nn) % 2 == 0 or math.gcd(mm, nn) != 1:
            continue
        s0, t0 = mm + nn, mm - nn
        for name, f in CHILDREN:
            em, en = [int(e.subs({m: mm, n: nn})) for e in euclid[name]]
            check(f(s0, t0) == (em + en, em - en), "numeric root form")
            cs, ct = f(s0, t0)
            check(cs > ct >= 1 and cs % 2 == 1 and ct % 2 == 1 and math.gcd(cs, ct) == 1,
                  "child in the odd coprime chamber")
            cnt += 1
print("  numeric check on all primitive Euclid pairs m<=60: %d children agree, all odd coprime s>t>=1"
      % cnt)
print("  (THM-3756 (25): L=B1, M=B2, R=B3 on the odd roots (q,d)=(2r-1,2s-1); CITED, not re-proved)")

# ---------------------------------------------------------------- S2
print("\n### S2  PROVED: transport table (a,b,k) -> child edge, symbolic identity in x")
a, b, x, K = sp.symbols("a b x K")  # K stands for 2^k
y = (a * x + b) / K
table_rows = []


def xtype(u, aa, bb):
    """x-type child {x,u}: identity  aa*x+bb = K*u ."""
    return sp.simplify(aa * x + bb - K * u) == 0


def ytype(u, alpha):
    """y-type child {y,u}: identity  a*u + alpha*b = (alpha*K + 2a)*y ."""
    return sp.simplify(a * u + alpha * b - (alpha * K + 2 * a) * y) == 0


# case A: y>x, (s,t)=(y,x)
cases = {
    "A (y>x): s=y,t=x": {
        "B1": ("x-type", B1(y, x), (a + 2 * K, b), None),
        "B2": ("y-type", B2(y, x), None, +1),
        "B3": ("y-type", B3(y, x), None, -1),
    },
    "B (y<x): s=x,t=y": {
        "B1": ("y-type", B1(x, y), None, +1),
        "B2": ("x-type", B2(x, y), (2 * K + a, b), None),
        "B3": ("x-type", B3(x, y), (2 * K - a, -b), None),
    },
}
for cname, d in cases.items():
    print("  case " + cname)
    for name in ("B1", "B2", "B3"):
        typ, child, fam, alpha = d[name]
        cs, ct = child
        # which coordinate is retained
        if typ == "x-type":
            if sp.simplify(ct - x) == 0:
                u = cs
            else:
                check(sp.simplify(cs - x) == 0, "retained x")
                u = ct
            aa, bb = fam
            check(xtype(u, aa, bb), "x-type identity " + cname + name)
            print("    %s: child (%s, %s) ; edge x -> u=%s ; family (a',b',k')=(%s, %s, k)  [x-type]"
                  % (name, sp.simplify(cs), sp.simplify(ct), sp.simplify(u), aa, bb))
            table_rows.append((cname, name, "x->u", str(aa), str(bb), "k"))
        else:
            if sp.simplify(ct - y) == 0:
                u = cs
            else:
                check(sp.simplify(cs - y) == 0, "retained y")
                u = ct
            check(ytype(u, alpha), "y-type identity " + cname + name)
            M = alpha * K + 2 * a
            print("    %s: child (%s, %s) ; edge u=%s -> y with  a*u %s b = (%s)*y   [y-type; E-edge (a,%sb,log2 M) iff M is a power of two]"
                  % (name, sp.simplify(cs), sp.simplify(ct), sp.simplify(u), "+" if alpha == 1 else "-", M,
                     "" if alpha == 1 else "-"))
            table_rows.append((cname, name, "u->y", "a", ("+b" if alpha == 1 else "-b"), "log2(%s)" % M))
print("  Both cases give the same three unordered child pairs {x,y+2x},{x+2y,y},{2y-x or 2x-y, y or x};")
print("  x-type children keep the SOURCE x and shift the multiplier by +-2^(k+1); y-type children keep")
print("  the TARGET y and the multiplier a, and are E-edges iff M = +-2^k + 2a is a positive power of two.")

# numeric check of the whole table on many edges of many families
def transport_children(xx, yy, aa, bb, kk):
    """Return list of (name, pair, typ, predicted) for the edge xx->yy of E_{aa,bb,kk}."""
    out = []
    if yy > xx:
        s0, t0 = yy, xx
        spec = {"B1": ("x", aa + 2 ** (kk + 1), bb), "B2": ("y", +1), "B3": ("y", -1)}
    else:
        s0, t0 = xx, yy
        spec = {"B1": ("y", +1), "B2": ("x", 2 ** (kk + 1) + aa, bb), "B3": ("x", 2 ** (kk + 1) - aa, -bb)}
    for name, f in CHILDREN:
        pr = f(s0, t0)
        sp_ = spec[name]
        if sp_[0] == "x":
            u = pr[0] if pr[1] == xx else pr[1]
            check((pr[0] == xx) != (pr[1] == xx) or pr[0] == pr[1], "x retained")
            a2, b2 = sp_[1], sp_[2]
            check(a2 * xx + b2 == 2 ** kk * u, "x-type numeric identity")
            out.append((name, pr, "x", (xx, u, a2, b2, kk)))
        else:
            alpha = sp_[1]
            u = pr[0] if pr[1] == yy else pr[1]
            M = alpha * 2 ** kk + 2 * aa
            check(aa * u + alpha * bb == M * yy, "y-type numeric identity")
            if M >= 2 and is_pow2(M):
                out.append((name, pr, "y", (u, yy, aa, alpha * bb, v2(M))))
            else:
                out.append((name, pr, "y", (u, yy, aa, alpha * bb, None, M)))
    return out


nchk = 0
for aa in range(1, 16, 2):
    for bb in range(-15, 16, 2):
        for xx in range(1, 400, 2):
            v = aa * xx + bb
            if v <= 0:
                continue
            kk = v2(v)
            yy = v >> kk
            if kk == 0 or yy == xx or math.gcd(xx, yy) != 1:
                continue
            transport_children(xx, yy, aa, bb, kk)
            nchk += 1
print("  numeric transport identities verified on %d edges (a<=15 odd, |b|<=15 odd, x<400): all pass" % nchk)

print("\n  lead's table for a=3, b=1, k=1 (x=3,7,...,31), children in root coordinates:")
print("  %6s %6s | %14s %-22s | %14s %-22s | %14s %-22s" % ("x", "y", "B1 pair", "family", "B2 pair", "family", "B3 pair", "family"))
for xx in range(3, 32, 4):
    yy = (3 * xx + 1) // 2
    ch = transport_children(xx, yy, 3, 1, 1)
    cells = []
    for name, pr, typ, pred in ch:
        if typ == "x":
            cells.append("%14s %-22s" % (str(pr), "x->%d in (%d,%d,%d)" % (pred[1], pred[2], pred[3], pred[4])))
        else:
            cells.append("%14s %-22s" % (str(pr), "%d->%d in (%d,%d,%d)" % (pred[0], pred[1], pred[2], pred[3], pred[4])))
    print("  %6d %6d | %s" % (xx, yy, " | ".join(cells)))

# ---------------------------------------------------------------- S3
print("\n### S3  PROVED: which children keep the multiplier")
print("  x-type: a' = a + 2^(k+1) > a, or a' = 2^(k+1) - a = a iff a = 2^k (impossible, a odd, k>=1).")
print("  So an x-type child NEVER keeps the multiplier.  y-type: a' = a always; E-edge iff")
print("  M = 2^k + 2a (alpha=+1) or M = 2a - 2^k (alpha=-1) is a positive power of two.")
print("  k=1: M = 2(a+1) or 2(a-1): power of two iff a = 2^j - 1 (alpha=+1) or a = 2^j + 1 (alpha=-1).")
print("  k>=2: 2^k+2a = 2(2^(k-1)+a) has odd part 2^(k-1)+a >= 3, never a power of two;")
print("        2a-2^k = 2(a-2^(k-1)) is a positive power of two iff a = 2^(k-1)+1 (then M=2, k'=1).")
print("  Exhaustive table for a=3, k<=12:")
keep = []
for kk in range(1, 13):
    for alpha in (+1, -1):
        M = alpha * 2 ** kk + 6
        if M >= 2 and is_pow2(M):
            keep.append((kk, alpha, M))
print("  (k, alpha, M) with M a positive power of two: %s" % keep)
check(keep == [(1, 1, 8), (1, -1, 4), (2, -1, 2)], "a=3 multiplier-keeping children")
print("  => for a=3: k=1: B2 -> (3, b, 3) and B3 -> (3, -b, 2) [y>x case, s=y]; k=1 with y<x (needs b<=-3):")
print("     B1 -> (3, b, 3); k=2 with y>x (needs b > x >= 1, so b>=3): B3 -> (3, -b, 1).  Nothing else.")
print("  SCOPE: negative multipliers 2^(k+1)-a<0 (x-type, y<x) need b<0 with |b|>(a-2^k)x; M=-2 (a=2^(k-1)-1)")
print("         would be multiplier -a.  Both are outside the positive-multiplier families and are not counted.")

# ---------------------------------------------------------------- S4
print("\n### S4  PROVED: for a=3, b=+-1 the orientation is the valuation: y>x iff k=1")
bad = 0
n_edges = 0
for bb in (1, -1):
    for xx in range(1, 100001, 2):
        v = 3 * xx + bb
        if v <= 0:
            continue
        kk = v2(v)
        yy = v >> kk
        if yy == xx:
            continue
        n_edges += 1
        if (yy > xx) != (kk == 1):
            bad += 1
print("  edges with odd source x<=10^5, b=+-1: %d ; violations of (y>x <=> k=1): %d" % (n_edges, bad))
check(bad == 0, "orientation = valuation")
print("  proof: k=1: y-x = (x+b)/2 > 0 for x>=3 (x=1,b=-1 is the fixed edge, excluded; x=1,b=1 gives k=2,y=1=x excluded);")
print("         k>=2: 2^k y = 3x+b <= 3x+1 < 4x <= 2^k x for x>=3, and x=1 only gives y=x.")

# ---------------------------------------------------------------- S5
print("\n### S5  PROVED+FINITE-EXACT: complete list of Berggren children of (3,+-1) edges that are (3,+-1) edges")
print("  Any reading of a child pair (p,q), p>q, as a (3,b',j) edge is either p->q with j>=2 or q->p with j=1")
print("  (3p+b' = 2^j q < 2^j p forces j>=2; 3q+b' = 2^j p > 2^j q forces j=1).  Each reading is a linear")
print("  equation c*x = d after y=(3x+b)/2^k.  Identities (c=d=0) are the generic transport; c!=0 gives sporadic")
print("  solutions.  The bound k<=4, j<=k+2 is proved in the note (cone ratios); here k,j<=40 are solved exactly.")
generic = []
sporadic = []
for bb in (1, -1):
    for kk in range(1, 41):
        for case in ("A", "B"):
            X = Fraction(1)  # coefficient of x
            # represent linear forms as (c0, c1): c0 + c1*x
            yform = (Fraction(bb, 2 ** kk), Fraction(3, 2 ** kk))
            xform = (Fraction(0), Fraction(1))
            if case == "A":
                sform, tform = yform, xform
            else:
                sform, tform = xform, yform

            def lin(cs, ct, f1, f2):
                return (cs * f1[0] + ct * f2[0], cs * f1[1] + ct * f2[1])

            kids = {"B1": (lin(1, 2, sform, tform), tform),
                    "B2": (lin(2, 1, sform, tform), sform),
                    "B3": (lin(2, -1, sform, tform), sform)}
            for name, (pf, qf) in kids.items():
                for jj in range(1, 41):
                    for bp in (1, -1):
                        for reading in ("p->q", "q->p"):
                            if reading == "p->q":
                                # 3p + bp = 2^j q
                                c = 3 * pf[1] - 2 ** jj * qf[1]
                                d = 2 ** jj * qf[0] - 3 * pf[0] - bp
                            else:
                                c = 3 * qf[1] - 2 ** jj * pf[1]
                                d = 2 ** jj * pf[0] - 3 * qf[0] - bp
                            if c == 0 and d == 0:
                                generic.append((bb, kk, case, name, reading, bp, jj))
                                continue
                            if c == 0:
                                continue
                            xs = d / c
                            if xs.denominator != 1 or xs <= 0 or xs.numerator % 2 == 0:
                                continue
                            xx = xs.numerator
                            vv = 3 * xx + bb
                            if vv <= 0 or v2(vv) != kk:
                                continue
                            yy = vv >> kk
                            if yy == xx or yy % 2 == 0:
                                continue
                            if (case == "A") != (yy > xx):
                                continue
                            p = int(pf[0] + pf[1] * xx)
                            q = int(qf[0] + qf[1] * xx)
                            src, tgt = (p, q) if reading == "p->q" else (q, p)
                            check(edge_reading(src, tgt, 3, bp) == jj, "sporadic recheck")
                            sporadic.append((xx, yy, bb, kk, name, (p, q), src, tgt, bp, jj))
print("  generic identities found (b, k, case, child, reading, b', j): %s" % generic)
check(sorted(generic) == sorted([(1, 1, "A", "B2", "p->q", 1, 3), (1, 1, "A", "B3", "p->q", -1, 2),
                                  (-1, 1, "A", "B2", "p->q", -1, 3), (-1, 1, "A", "B3", "p->q", 1, 2),
                                  (1, 1, "B", "B1", "p->q", 1, 3), (-1, 1, "B", "B1", "p->q", -1, 3),
                                  (1, 2, "A", "B3", "p->q", -1, 1), (-1, 2, "A", "B3", "p->q", 1, 1)]),
      "generic identities")
print("  (the case-B k=1 identities need y<x at k=1 (x<-b) and the case-A k=2 identities need y>x at k=2 (x<b):")
print("   both are vacuous for b=+-1 by S4 and are the S3 rows k=1,alpha=+1 and k=2,alpha=-1 of general b)")
sporadic = sorted(set(sporadic))
print("  sporadic solutions (x, y, b, k, child, child pair, reading src->tgt, b', j):")
for row in sporadic:
    print("    edge %d->%d in (3,%d,%d): %s child %s read as %d->%d in (3,%d,%d)" % (
        row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9]))
expected_sporadic = {(3, 5, 1, 1, "B3", (7, 5), 5, 7, -1, 1),
                     (7, 5, -1, 2, "B2", (19, 7), 19, 7, -1, 3),
                     (7, 5, -1, 2, "B3", (9, 7), 9, 7, 1, 2),
                     (3, 1, -1, 3, "B1", (5, 1), 5, 1, 1, 4),
                     (3, 1, -1, 3, "B3", (5, 3), 3, 5, 1, 1)}
check(set(sporadic) == expected_sporadic, "sporadic list")
print("  exactly %d sporadic readings (E1..E5), all with x<=7." % len(sporadic))

# brute-force cross-check over all (3,+-1) edges with max(x,y)<=5000
LIM5 = 5000
brute_generic = 0
brute_sporadic = set()
n_par = 0
for bb in (1, -1):
    for xx in range(1, 4 * LIM5, 2):
        v = 3 * xx + bb
        if v <= 0:
            continue
        kk = v2(v)
        yy = v >> kk
        if yy == xx or max(xx, yy) > LIM5:
            continue
        n_par += 1
        s0, t0 = roots(xx, yy)
        for name, f in CHILDREN:
            pr = f(s0, t0)
            for (src, tgt, a3, bp, jj) in legal_readings(pr, 3, (1, -1)):
                gen = (kk == 1 and ((name == "B2" and (src, tgt, bp, jj) == (4 * xx + bb, yy, bb, 3)) or
                                    (name == "B3" and (src, tgt, bp, jj) == (2 * xx + bb, yy, -bb, 2))))
                if gen:
                    brute_generic += 1
                else:
                    brute_sporadic.add((xx, yy, bb, kk, name, pr, src, tgt, bp, jj))
print("  brute force over %d edges with max(x,y)<=%d: generic legal children %d, sporadic %d"
      % (n_par, LIM5, brute_generic, len(brute_sporadic)))
check(brute_sporadic == expected_sporadic, "brute-force sporadic list agrees")
n_k1 = sum(1 for bb in (1, -1) for xx in range(1, 4 * LIM5, 2)
           if 3 * xx + bb > 0 and v2(3 * xx + bb) == 1 and (3 * xx + bb) // 2 != xx and (3 * xx + bb) // 2 <= LIM5)
print("  k=1 edges in that range: %d ; generic legal children = 2 per k=1 edge: %d" % (n_k1, 2 * n_k1))
check(brute_generic == 2 * n_k1, "two generic children per k=1 edge")

# ---------------------------------------------------------------- S6
print("\n### S6  PROVED (from S3+S5): the legal sub-forest of the Berggren tree for (3,+-1)")
print("  A (3,+-1) k=1 edge x->y has exactly the two generic legal children (4x+b,y) in (3,b,3) and")
print("  (2x+b,y) in (3,-b,2); both are k>=2 edges, hence leaves, unless sporadic.  A k>=2 edge has no")
print("  legal child except the sporadic parents 7->5 (b=-1,k=2) and 3->1 (b=-1,k=3).")
print("  The claim 'depth <= 1 above every k=1 edge' is REFUTED by exactly one witness:")
print("  3->5 (b=1,k=1) -B3-> (7,5) = 5->7 (b=-1,k=1) -B2-> (19,7) = 19->7 (b=-1,k=3) [and -B3-> (9,7) = 9->7 (b=1,k=2)].")
print("  Corrected statement: every k=1 edge has exactly two legal children; they are leaves except (7,5),")
print("  which is the k=1 edge 5->7 with its own two leaf children.  So 3->5 is the only k=1 edge of depth 2,")
print("  every other k=1 edge has depth exactly 1, and the root cluster {(3,1),(5,1),(5,3),(13,5),(7,5),(19,7),(9,7)}")
print("  (a tree of depth 3 rooted at the Berggren root (3,1)) is the unique legal component of size > 3.")
# verify the legal components in the range
LIM6 = 2000
legal_pairs = {}
for bb in (1, -1):
    for xx in range(1, 4 * LIM6, 2):
        v = 3 * xx + bb
        if v <= 0:
            continue
        kk = v2(v)
        yy = v >> kk
        if yy == xx or max(xx, yy) > LIM6:
            continue
        legal_pairs.setdefault(roots(xx, yy), []).append((xx, yy, bb, kk))
deg = defaultdict(int)
for pr, reads in legal_pairs.items():
    nlegal = sum(1 for name, f in CHILDREN if legal_readings(f(*pr), 3, (1, -1)))
    kmin = min(r[3] for r in reads)
    deg[(kmin, nlegal)] += 1
print("  legal pairs with max<=%d: %d ; (min k over readings, number of legal children, tested directly) census:" % (LIM6, len(legal_pairs)))
for key in sorted(deg):
    print("    k=%d, legal children=%d : %d pairs" % (key[0], key[1], deg[key]))
check(all((k == 1 and c == 2) or (k >= 2 and c == 0) or (k, c) == (3, 2) for (k, c) in deg) and deg[(3, 2)] == 1,
      "degree census: k=1 pairs have exactly 2 legal children, k>=2 pairs none except (3,1)")
cluster = [(3, 1), (5, 1), (5, 3), (13, 5), (7, 5), (19, 7), (9, 7)]
for pr in cluster:
    check(pr in legal_pairs, "cluster pair legal %s" % (pr,))
    print("    cluster pair %s readings %s legal children %s" % (
        pr, legal_pairs[pr], [(name, f(*pr)) for name, f in CHILDREN if legal_readings(f(*pr), 3, (1, -1))]))

# ---------------------------------------------------------------- S7
print("\n### S7  PROVED: parent transport on the three THM-3756 cones")
print("  Inverse maps: s/t>3: B1^-1 (s,t)->(s-2t,t); 2<s/t<3: B2^-1 (s,t)->(t,s-2t); 1<s/t<2: B3^-1 (s,t)->(t,2t-s).")
print("  For a (3,b) edge with b=+-1 the cone is a function of k alone:")
print("    k=1 (s=y,t=x, y/x in (1,2)): parent (x,(x-b)/2) = the halving edge x->(x-b)/2 of E_{1,-b,1};")
print("      forward: (1,-b,1) -B3 (y<x, x-type)-> (2^2-1, b, 1) = (3,b,1).")
print("      Opus orientation (x-b)/2 -> x is a (3, b'', 1) edge with the affine law b'' = (x+3b)/2.")
print("    k=2 (s=x,t=y, x/y in (1,2)): parent (y,(x+b)/2) = the k=1 edge (x+b)/2 -> y of E_{3,-b,1} (inverse of B3).")
print("    k=3 (x/y in (2,3)): parent (y,(x-b)/4) = the k=1 edge (x-b)/4 -> y of E_{3,b,1} (inverse of B2 = inverse braid).")
print("    k>=4 (x/y>3): parent (x-2y,y) satisfies 3(x-2y)+b = (2^k-6) y, odd cofactor 2^(k-1)-3 >= 5: never an E_{3,.,.} edge by identity.")
npar = defaultdict(int)
for bb in (1, -1):
    for xx in range(1, 8001, 2):
        v = 3 * xx + bb
        if v <= 0:
            continue
        kk = v2(v)
        yy = v >> kk
        if yy == xx or max(xx, yy) > 2000:
            continue
        s0, t0 = roots(xx, yy)
        r = Fraction(s0, t0)
        if r > 3:
            cone, par = "C1", (s0 - 2 * t0, t0)
        elif r > 2:
            cone, par = "C2", (t0, s0 - 2 * t0)
        else:
            cone, par = "C3", (t0, 2 * t0 - s0)
        if (s0, t0) == (3, 1):
            cone, par = "root", None
        expected_cone = {1: "C3", 2: "C3", 3: "C2"}.get(kk, "C1")
        if (s0, t0) != (3, 1):
            check(cone == expected_cone, "cone by k")
            if kk == 1:
                check(par == (xx, (xx - bb) // 2) and edge_reading(xx, (xx - bb) // 2, 1, -bb) == 1, "k=1 parent")
                bpp = (xx + 3 * bb) // 2
                check(edge_reading((xx - bb) // 2, xx, 3, bpp) == 1, "opus orientation law")
            elif kk == 2:
                check(par == (yy, (xx + bb) // 2) and edge_reading((xx + bb) // 2, yy, 3, -bb) == 1, "k=2 parent")
            elif kk == 3:
                check(par == (yy, (xx - bb) // 4) and edge_reading((xx - bb) // 4, yy, 3, bb) == 1, "k=3 parent")
            else:
                check(par == (xx - 2 * yy, yy) and 3 * (xx - 2 * yy) + bb == (2 ** kk - 6) * yy, "k>=4 parent")
        npar[(kk, cone)] += 1
print("  verified on all (3,+-1) edges with max(x,y)<=2000: (k,cone) counts %s" % dict(sorted(npar.items())))
print("  example 7->11 (b=1,k=1): parent (7,3); 7->3 in (1,-1,1): 7-1=6=2*3 ; 3->7 in (3,5,1): 9+5=14=2*7, b''=(7+3)/2=5")
check(edge_reading(7, 3, 1, -1) == 1 and edge_reading(3, 7, 3, 5) == 1, "7->11 parent")
print("  general inverse table (any odd a): the parent PAIR is unique (THM-3756 cones).  It is an identity")
print("  pre-image of (a,b,k) under S2 exactly when it reads as: (a-2^(k+1), b, k) [B1 of a y>x parent],")
print("  (a-2^(k+1), b, k) [B2 of a y<x parent], (2^(k+1)-a, -b, k) [B3 of a y<x parent] (x-type; needs the")
print("  retained coordinate to be the parent's source), or (a, +-b, k') with +-2^k' + 2a = 2^k (y-type).")
print("  Otherwise the parent is only a generalized edge a*u + b = 2^j*m*y with odd cofactor m>1 (k>=4 above).")
# symbolic check of the inverse statements
for kk in range(1, 6):
    Kk = 2 ** kk
    for (ap, bp, name) in ((a - 2 * Kk, b, "B1"), (a - 2 * Kk, b, "B2"), (2 * Kk - a, -b, "B3")):
        yp = (ap * x + bp) / Kk  # parent target
        if name == "B1":
            child = (yp + 2 * x, x)  # parent case A
        elif name == "B2":
            child = (2 * x + yp, x)
        else:
            child = (2 * x - yp, x)
        u = child[0]
        check(sp.simplify(a * x + b - Kk * u) == 0, "inverse x-type k=%d %s" % (kk, name))
print("  symbolic inverse x-type identities verified for k<=5.")

# ---------------------------------------------------------------- S8
print("\n### S8  FINITE-EXACT: census of coprime odd pairs s<=2000, all guards a<=15, |b|<=15")
SMAX = 2000
ss = np.arange(1, SMAX + 1, 2, dtype=np.int64)
S, T = np.meshgrid(ss, ss, indexing="ij")
mask = (T < S) & (np.gcd(S, T) == 1)
S = S[mask]
T = T[mask]
npairs = S.size
print("  coprime odd pairs (s,t), t<s<=%d: %d" % (SMAX, npairs))
real = defaultdict(list)  # (s,t) -> list of (src,tgt,a,b,k)
tot = 0
per_ab = {}
for aa in range(1, 16, 2):
    for bb in range(-15, 16, 2):
        c_ab = 0
        for (src, tgt) in ((S, T), (T, S)):
            v = aa * src + bb
            ok = v > 0
            q = np.where(ok, v // np.where(tgt > 0, tgt, 1), 0)
            ok &= (v % tgt == 0) & (q >= 2) & ((q & (q - 1)) == 0)
            idx = np.nonzero(ok)[0]
            for i in idx:
                kk = int(q[i]).bit_length() - 1
                real[(int(S[i]), int(T[i]))].append((int(src[i]), int(tgt[i]), aa, bb, kk))
                c_ab += 1
        per_ab[(aa, bb)] = c_ab
        tot += c_ab
print("  total directed realizations (pair, orientation, a, b, k): %d ; pairs with >=1 realization: %d (%.2f%%)"
      % (tot, len(real), 100.0 * len(real) / npairs))
hist = defaultdict(int)
for pr, L in real.items():
    hist[len(L)] += 1
print("  histogram of realization counts per pair: %s" % dict(sorted(hist.items())))
print("  realizations per (a,b), a=3 row: %s" % {bb: per_ab[(3, bb)] for bb in range(-15, 16, 2)})
print("  realizations per a (summed over b): %s" % {aa: sum(per_ab[(aa, bb)] for bb in range(-15, 16, 2)) for aa in range(1, 16, 2)})
# transport verification on all children of all realizations
n_id = 0
n_pred_in = 0
n_pred_out_of_range = 0
n_ytype_nonpow = 0
for pr, L in real.items():
    for (src, tgt, aa, bb, kk) in L:
        for name, cpair, typ, pred in transport_children(src, tgt, aa, bb, kk):
            n_id += 1
            if typ == "y" and pred[4] is None:
                n_ytype_nonpow += 1
                continue
            a2, b2, k2 = pred[2], pred[3], pred[4]
            if 1 <= a2 <= 15 and abs(b2) <= 15 and cpair[0] <= SMAX:
                check((pred[0], pred[1], a2, b2, k2) in real[cpair], "predicted child realization present")
                n_pred_in += 1
            else:
                n_pred_out_of_range += 1
print("  children examined: %d ; identities hold on all; predicted child guard found in census: %d ;"
      % (n_id, n_pred_in))
print("  predicted guard outside the census window (a'>15 or |b'|>15 or s'>2000): %d ; y-type with M not a power of two: %d"
      % (n_pred_out_of_range, n_ytype_nonpow))
check(n_id == n_pred_in + n_pred_out_of_range + n_ytype_nonpow, "bookkeeping")
# small table of the (3,+-1) readings for s<=31
print("  (3,+-1) readings for pairs with s<=31:")
for pr in sorted(real):
    if pr[0] > 31:
        break
    rr = [r for r in real[pr] if r[2] == 3 and abs(r[3]) == 1]
    if rr:
        print("    %s : %s" % (pr, ["%d->%d (3,%d,%d)" % (r[0], r[1], r[3], r[4]) for r in rr]))

# ---------------------------------------------------------------- S9
print("\n### S9  FINITE-EXACT: primitive triangles with hypotenuse <= 10^6 that are Collatz-family edges")
X = 10 ** 6
smax9 = int(math.isqrt(2 * X)) + 2
ss9 = np.arange(1, smax9 + 1, 2, dtype=np.int64)
S9, T9 = np.meshgrid(ss9, ss9, indexing="ij")
mask9 = (T9 < S9) & ((S9 * S9 + T9 * T9) <= 2 * X) & (np.gcd(S9, T9) == 1)
S9 = S9[mask9]
T9 = T9[mask9]
n_tri = S9.size
print("  primitive triangles (s>t odd coprime, (s^2+t^2)/2 <= 10^6): %d" % n_tri)


def count_b(bs):
    directed = defaultdict(int)
    tri = np.zeros(S9.size, dtype=bool)
    for bb in bs:
        for (src, tgt) in ((S9, T9), (T9, S9)):
            v = 3 * src + bb
            ok = v > 0
            q = np.where(ok, v // tgt, 0)
            ok &= (v % tgt == 0) & (q >= 2) & ((q & (q - 1)) == 0)
            directed[bb] += int(ok.sum())
            tri |= ok
    return dict(directed), int(tri.sum())


d1, t1 = count_b((1,))
dm1, tm1 = count_b((-1,))
dpm, tpm = count_b((1, -1))
d5, t5 = count_b((1, -1, 3, -3, 5, -5))
C_edge = sum(1 / math.sqrt(2 * (4 ** k + 9)) for k in range(1, 200))
print("  3x+1 : directed edges %d, distinct triangles %d ; opus law C_edge*sqrt(X) = %.4f*1000 = %.1f"
      % (d1[1], t1, C_edge, C_edge * 1000))
print("  3x-1 : directed edges %d, distinct triangles %d" % (dm1[-1], tm1))
print("  union b=+-1: distinct triangles %d (directed %d; overlap = triangles legal in both orientations/parameters: %d)"
      % (tpm, d1[1] + dm1[-1], d1[1] + dm1[-1] - tpm))
print("  union |b|<=5: directed by b %s ; distinct triangles %d = %.3f%% of all primitive triangles"
      % (dict(sorted(d5.items())), t5, 100.0 * t5 / n_tri))
check(abs(d1[1] - C_edge * 1000) < 40, "C sqrt X law within O(log X)")

# ---------------------------------------------------------------- S10
print("\n### S10 PROVED: the inverse-fibre braid is a Berggren word")
print("  Fibre of the target y: x_k = (2^k y - b)/3 over the k of one parity class (2^k y = b mod 3).")
print("  x_{k+2} = 4 x_k + b = x_k + 2^k y.  For k>=2 (x_k > y): (x_{k+2}, y) = B1^(2^(k-1)) (x_k, y)")
print("  since B1^m (s,t) = (s+2mt, t) and 2^k y = 2*2^(k-1)*y.  For k=1 (y > x_1): (x_3, y) = B2 (y, x_1).")
print("  Hence k=1 -> k=2j+1 is the word B1^((4^j-4)/3) B2 (B2 first), and k0>=2 -> k0+2j is B1^((2^(k0+2j)-2^k0)/6).")
print("  The fibre is a subset of the B1-ray {(x,y): t=y fixed} at the B1-heights m_j = 2^(k0-1)(4^j-1)/3.")
nfib = 0
for bb in (1, -1):
    for yy in range(1, 2002, 2):
        if yy % 3 == 0:
            continue
        for kk in range(1, 13):
            if (2 ** kk * yy - bb) % 3:
                continue
            xk = (2 ** kk * yy - bb) // 3
            xk2 = (2 ** (kk + 2) * yy - bb) // 3
            if xk <= 0 or xk == yy:
                continue
            check(xk2 == 4 * xk + bb, "R = 4x+b")
            if kk == 1:
                check(yy > xk or bb < 0, "k=1 orientation")
                if yy > xk:
                    check(B2(yy, xk) == (xk2, yy), "k=1 -> k=3 is B2")
            else:
                check(xk > yy, "k>=2 orientation")
                s0, t0 = xk, yy
                for _ in range(2 ** (kk - 1)):
                    s0, t0 = B1(s0, t0)
                check((s0, t0) == (xk2, yy), "B1^(2^(k-1)) step")
            nfib += 1
print("  verified on %d fibre steps (y<=2001, k<=12, b=+-1)." % nfib)
# the lead's example in Euclid coordinates
print("  lead's example: Euclid (9,2)=(s,t)(11,7): 7->11 k=1; (20,9)=(29,11): 29->11 k=3; (64,53)=(117,11): 117->11 k=5.")
check(edge_reading(7, 11, 3, 1) == 1 and edge_reading(29, 11, 3, 1) == 3 and edge_reading(117, 11, 3, 1) == 5, "lead example")
check(B2(11, 7) == (29, 11), "B2 step")
p = (29, 11)
for _ in range(4):
    p = B1(*p)
check(p == (117, 11), "B1^4 step")
kids = [f(29, 11) for _, f in CHILDREN]
check((117, 11) not in kids, "not a single tree step")
print("  (29,11) -> (117,11) is B1^4, not a child (children of (29,11): %s)." % kids)

# ---------------------------------------------------------------- S11
print("\n### S11 FINITE-EXACT+PROVED: consecutive orbit triangles are never parent/child (b=+1)")
LIM11 = 10 ** 5
npairs11 = 0
nhit = 0
for x0 in range(1, LIM11 + 1, 2):
    xx = x0
    prev = None
    while xx != 1:
        v = 3 * xx + 1
        yy = v >> v2(v)
        cur = roots(xx, yy)
        if prev is not None:
            npairs11 += 1
            if cur in (B1(*prev), B2(*prev), B3(*prev)) or prev in (B1(*cur), B2(*cur), B3(*cur)):
                nhit += 1
        prev = cur
        xx = yy
print("  consecutive orbit-triangle pairs over odd starts <= %d: %d ; parent/child incidences: %d" % (LIM11, npairs11, nhit))
check(nhit == 0, "no consecutive parent/child")
print("  proof sketch (note S11): the shared coordinate y must be the retained coordinate, so z in {x+2y, 2y-x}")
print("  or x in {z+2y, 2y+z, 2y-z}; with y=(3x+1)/2^k, z=(3y+1)/2^j each case is a linear equation whose only")
print("  positive odd solution is the fixed edge x=y=1.")
print("\nALL CHECKS PASSED")
