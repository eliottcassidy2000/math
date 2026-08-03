#!/usr/bin/env python3
"""
tv_fusion_bin2_verification_boxeph.py  (boxeph, 2026-08-03)

Exact verification of the n=2 TV-fusion / homogenization lemma
(Kontorovich, arXiv 2601.04079v3, "TV Homogenization Inequalities").

CLAIM 1 (closed form).  For x,y in [0,1],
    TV(Bin(2,x), Bin(2,y)) = |x-y| * (1 + |x+y-1|),
where TV = (1/2) sum_{k=0,1,2} |P(k)-Q(k)|.

CLAIM 2 (fusion inequality).  For p1,p2,q1,q2 in [0,1], with
pbar=(p1+p2)/2, qbar=(q1+q2)/2, delta_i=|p_i-q_i|:
    TV(Bin(2,pbar), Bin(2,qbar)) <= delta_1 + delta_2 - delta_1*delta_2
                                  =  1 - (1-delta_1)(1-delta_2).

CLAIM 3 (rigid equality face).  Equality holds exactly on
    F0 = {p1=q1, p2=q2}                          (both sides 0)
    F1 = {p1=p2=t, q1=q2=0},  F2 = {p1=p2=1, q1=q2=t},
    F3 = {p1=p2=0, q1=q2=t},  F4 = {q1=q2=1, p1=p2=t},   t in [0,1].
i.e. nontrivially rigid iff both blocks carry identical pairs
(p1=p2, q1=q2) and one of the two fused distributions is deterministic
(bias in {0,1}).

METHOD: exact only.  sympy (symbolic polynomial identities, sign-resolved
per closed region) + fractions.Fraction hostile scans.  No floats in any
decision; no numpy.
"""

import random
from fractions import Fraction as F

import sympy as sp

HALF = sp.Rational(1, 2)
PASS_ALL = True


def check(name, cond):
    global PASS_ALL
    tag = "PASS" if cond else "FAIL"
    if not cond:
        PASS_ALL = False
    print("[%s] %s" % (tag, name))


print("=" * 78)
print("PART 1: closed form TV(Bin(2,x),Bin(2,y)) = |x-y| (1+|x+y-1|)")
print("=" * 78)

x, y, s, d = sp.symbols("x y s d", real=True)
to_sd = {x: (s + d) / 2, y: (s - d) / 2}

# The three absolute-value arguments of the TV definition, factored
# through s = x+y, d = x-y (polynomial identities, hence valid everywhere).
A2 = x ** 2 - y ** 2                     # k=2 term
A1 = x * (1 - x) - y * (1 - y)           # k=1 term (before the factor 2)
A0 = (1 - x) ** 2 - (1 - y) ** 2         # k=0 term
check("factor k=2:  x^2-y^2            == d*s",
      sp.expand(A2.subs(to_sd, simultaneous=True) - d * s) == 0)
check("factor k=1:  x(1-x)-y(1-y)      == d*(1-s)",
      sp.expand(A1.subs(to_sd, simultaneous=True) - d * (1 - s)) == 0)
check("factor k=0:  (1-x)^2-(1-y)^2    == -d*(2-s)",
      sp.expand(A0.subs(to_sd, simultaneous=True) + d * (2 - s)) == 0)

# Sign resolution.  On the box x,y in [0,1] we have s in [0,2], so the
# factors s and 2-s are nonnegative THROUGHOUT; only sgn(d) and sgn(1-s)
# vary.  Resolve |.| on each of the four closed regions
#   {d>=0 or d<=0} x {s<=1 or s>=1}
# and check the resulting POLYNOMIAL identity.  Region closures cover the
# box, so the closed form holds everywhere on it.
for eps in (1, -1):                       # eps = sign of d, so |d| = eps*d
    for branch, absum in (("s<=1", 1 - s), ("s>=1", s - 1)):
        abs_A2 = eps * d * s              # |d*s|      = |d| * s
        abs_A1 = eps * d * absum          # |d*(1-s)|  = |d| * |1-s|
        abs_A0 = eps * d * (2 - s)        # |d*(2-s)|  = |d| * (2-s)
        tv_resolved = HALF * (abs_A2 + 2 * abs_A1 + abs_A0)
        target = eps * d * (1 + absum)    # |d| * (1+|s-1|)
        check("region  sgn(d)=%+d, %s :  (1/2)(|A2|+2|A1|+|A0|) == |d|(1+|s-1|)"
              % (eps, branch),
              sp.expand(tv_resolved - target) == 0)

# Independent exact-rational check of closed form vs raw TV definition.
def tv_def(a, b):
    return (abs(a * a - b * b)
            + 2 * abs(a * (1 - a) - b * (1 - b))
            + abs((1 - a) * (1 - a) - (1 - b) * (1 - b))) / 2


def tv_closed(a, b):
    return abs(a - b) * (1 + abs(a + b - 1))


grid33 = [F(k, 32) for k in range(33)]
mism = sum(1 for a in grid33 for b in grid33 if tv_def(a, b) != tv_closed(a, b))
check("exact Fraction grid 33x33 (step 1/32): definition == closed form "
      "(%d mismatches)" % mism, mism == 0)

print()
print("=" * 78)
print("PART 2: fusion inequality  TV(Bin(2,pbar),Bin(2,qbar)) <= d1+d2-d1*d2")
print("=" * 78)

p1, p2, q1, q2 = sp.symbols("p1 p2 q1 q2", real=True)


def mkL(P1, P2, Q1, Q2):
    pb = (P1 + P2) / 2
    qb = (Q1 + Q2) / 2
    return sp.Abs(sp.expand(pb - qb)) * (1 + sp.Abs(sp.expand(pb + qb - 1)))


def mkR(P1, P2, Q1, Q2):
    D1 = sp.Abs(sp.expand(P1 - Q1))
    D2 = sp.Abs(sp.expand(P2 - Q2))
    return D1 + D2 - D1 * D2


L0 = mkL(p1, p2, q1, q2)
R0 = mkR(p1, p2, q1, q2)

# --- 2a. Symmetry group G = <sigma, pi, kappa> (box automorphisms) -----
# sigma: swap blocks (1<->2);  pi: swap p<->q;  kappa: complement x->1-x.
transforms = {
    "sigma (block swap)":  ((p2, p1, q2, q1)),
    "pi    (p<->q swap)":  ((q1, q2, p1, p2)),
    "kappa (complement)":  ((1 - p1, 1 - p2, 1 - q1, 1 - q2)),
}
for name, args in transforms.items():
    check("LHS invariant under %s" % name, sp.simplify(mkL(*args) - L0) == 0)
    check("RHS invariant under %s" % name, sp.simplify(mkR(*args) - R0) == 0)

# Action of G on (d1, d2, S) with d_i = p_i - q_i, S = pbar + qbar:
d1e, d2e = p1 - q1, p2 - q2
Se = (p1 + p2) / 2 + (q1 + q2) / 2
act = {
    "sigma": ((p2, p1, q2, q1), (d2e, d1e, Se)),          # swaps d1,d2 ; S fixed
    "pi":    ((q1, q2, p1, p2), (-d1e, -d2e, Se)),        # negates d   ; S fixed
    "kappa": ((1 - p1, 1 - p2, 1 - q1, 1 - q2),
              (-d1e, -d2e, 2 - Se)),                      # negates d   ; S -> 2-S
}
for name, (args, (D1t, D2t, St)) in act.items():
    sub = dict(zip((p1, p2, q1, q2), args))
    ok = (sp.expand(d1e.subs(sub, simultaneous=True) - D1t) == 0
          and sp.expand(d2e.subs(sub, simultaneous=True) - D2t) == 0
          and sp.expand(Se.subs(sub, simultaneous=True) - St) == 0)
    check("action of %s on (d1,d2,S) as claimed" % name, ok)

# --- 2b. Coverage: every sign/branch region maps into a canonical one ---
# Label a region by (sgn d1, sgn d2, branch) with branch=+1 for S<=1.
# sigma: (s1,s2,b)->(s2,s1,b); pi: ->(-s1,-s2,b); kappa: ->(-s1,-s2,-b).
# Canonical targets: B=(+1,+1,+1) and A=(+1,-1,+1).
def apply_word(lbl, word):
    s1, s2, b = lbl
    for g in word:
        if g == "sigma":
            s1, s2 = s2, s1
        elif g == "pi":
            s1, s2 = -s1, -s2
        else:  # kappa
            s1, s2, b = -s1, -s2, -b
    return (s1, s2, b)


words = [tuple(w) for w in
         [[], ["sigma"], ["pi"], ["kappa"], ["sigma", "pi"],
          ["sigma", "kappa"], ["pi", "kappa"], ["sigma", "pi", "kappa"]]]
targets = {(1, 1, 1), (1, -1, 1)}
uncovered = []
for s1 in (1, -1):
    for s2 in (1, -1):
        for b in (1, -1):
            if not any(apply_word((s1, s2, b), w) in targets for w in words):
                uncovered.append((s1, s2, b))
check("all 8 open sign/branch regions reachable to canonical {B, A} "
      "(uncovered: %s)" % uncovered, not uncovered)
print("       (boundaries d_i=0, S=1 lie in region closures; the canonical")
print("        certificates below are weak inequalities on CLOSED regions,")
print("        so the union of closures covers the whole box.)")
print("       Inside case A, sigma*pi maps (d1,d2)=(a,-b) -> (b,-a) with S")
print("       fixed, so WLOG a >= b there (checked: sigma,pi actions above).")

# --- 2c. Canonical case B: d1,d2 >= 0, S <= 1 --------------------------
# Parametrize q_i = c_i >= 0, d_i = t_i >= 0, p_i = c_i + t_i (<=1).
c1, c2, t1, t2 = sp.symbols("c1 c2 t1 t2", nonnegative=True)
subB = {p1: c1 + t1, p2: c2 + t2, q1: c1, q2: c2}
S_B = c1 + c2 + (t1 + t2) / 2
check("case B parametrization: pbar-qbar == (t1+t2)/2",
      sp.expand(((p1 + p2) / 2 - (q1 + q2) / 2).subs(subB) - (t1 + t2) / 2) == 0)
check("case B parametrization: pbar+qbar == c1+c2+(t1+t2)/2",
      sp.expand(Se.subs(subB) - S_B) == 0)
# On this region |pbar-qbar| = (t1+t2)/2 and 1+|S-1| = 2-S (uses S<=1 only).
LHS_B = (t1 + t2) / 2 * (2 - S_B)
RHS_B = t1 + t2 - t1 * t2
cert_B = (t1 + t2) * (c1 + c2) / 2 + (t1 - t2) ** 2 / 4
check("case B certificate identity: RHS - LHS == "
      "(t1+t2)(c1+c2)/2 + (t1-t2)^2/4",
      sp.expand(RHS_B - LHS_B - cert_B) == 0)
print("       both certificate terms are products/squares of nonnegative")
print("       quantities (c_i,t_i >= 0), hence RHS-LHS >= 0 on case B.")
print("       AM-GM reading: 1+|S-1| = 2-S <= 2-(t1+t2)/2, so LHS <=")
print("       (T)(2-T/2)/2 = T - T^2/4 <= T - t1*t2 with T=t1+t2 (AM-GM).")

# --- 2d. Canonical case A: d1 = a >= 0 >= d2 = -b, a >= b, S <= 1 ------
# Parametrize q1 = u >= 0, p1 = u+a (<=1  =>  a <= 1); p2 = v >= 0, q2 = v+b.
u, v, a, b = sp.symbols("u v a b", nonnegative=True)
subA = {p1: u + a, q1: u, p2: v, q2: v + b}
S_A = u + v + (a + b) / 2
check("case A parametrization: pbar-qbar == (a-b)/2",
      sp.expand(((p1 + p2) / 2 - (q1 + q2) / 2).subs(subA) - (a - b) / 2) == 0)
check("case A parametrization: pbar+qbar == u+v+(a+b)/2",
      sp.expand(Se.subs(subA) - S_A) == 0)
# On this region |d1+d2| = a-b (a>=b) and 1+|S-1| = 2-S (S<=1).
LHS_A = (a - b) / 2 * (2 - S_A)
RHS_A = a + b - a * b
cert_A = b * (1 - a) + b + (a - b) * (u + v) / 2 + (a - b) * (a + b) / 4
check("case A certificate identity: RHS - LHS == "
      "b(1-a) + b + (a-b)(u+v)/2 + (a-b)(a+b)/4",
      sp.expand(RHS_A - LHS_A - cert_A) == 0)
print("       terms nonnegative on case A: b>=0; 1-a>=0 since a=p1-q1<=1;")
print("       a-b>=0 by the WLOG; u,v>=0; a+b>=0.  Hence RHS-LHS >= 0.")

print()
print("=" * 78)
print("PART 3: equality face -- symbolic + exhaustive exact grid scan")
print("=" * 78)
print("""Symbolic equality analysis (certificates are sums of nonnegative terms,
so equality forces every term to vanish):
  case A: b=0, then (a-b)(a+b)/4 = a^2/4 = 0 forces a=0  => only d1=d2=0.
  case B: (t1-t2)^2=0 and (t1+t2)(c1+c2)=0
          => t1=t2=0 (c1,c2 free)   [the diagonal p_i=q_i, both sides 0], or
          => t1=t2=t, c1=c2=0       [p1=p2=t, q1=q2=0].
Pulling {p1=p2=t, q1=q2=0} back through G = <sigma,pi,kappa> gives exactly
F1..F4; the diagonal is G-invariant and gives F0.""")

# Symbolic confirmation that each face achieves equality (signs resolved
# with t in [0,1]: |t| = t, |t-1| = 1-t).
t = sp.symbols("t", nonnegative=True)
faces_symbolic = [
    ("F1: p1=p2=t, q1=q2=0", t * (1 + (1 - t)), 2 * t - t * t),
    ("F2: p1=p2=1, q1=q2=t", (1 - t) * (1 + t), 2 * (1 - t) - (1 - t) ** 2),
    ("F3: p1=p2=0, q1=q2=t", t * (1 + (1 - t)), 2 * t - t * t),
    ("F4: q1=q2=1, p1=p2=t", (1 - t) * (1 + t), 2 * (1 - t) - (1 - t) ** 2),
]
for name, lhs_e, rhs_e in faces_symbolic:
    check("equality on %s  (LHS==RHS identically)" % name,
          sp.expand(lhs_e - rhs_e) == 0)


def on_face(P1, P2, Q1, Q2):
    if P1 == Q1 and P2 == Q2:
        return True
    return (P1 == P2 and Q1 == Q2
            and (P1 in (0, 1) or Q1 in (0, 1)))


# Exhaustive exact scan: all (p1,p2,q1,q2) in {0,1/16,...,1}^4  (17^4 pts).
vals = [F(k, 16) for k in range(17)]
viol = 0
closed_mism = 0
eq_pts = 0
face_pts = 0
eq_face_mismatch = []
for P1 in vals:
    for P2 in vals:
        for Q1 in vals:
            for Q2 in vals:
                pb = (P1 + P2) / 2
                qb = (Q1 + Q2) / 2
                lhs = tv_def(pb, qb)            # raw TV definition (hostile)
                if lhs != tv_closed(pb, qb):
                    closed_mism += 1
                D1 = abs(P1 - Q1)
                D2 = abs(P2 - Q2)
                rhs = D1 + D2 - D1 * D2
                if lhs > rhs:
                    viol += 1
                fc = on_face(P1, P2, Q1, Q2)
                if fc:
                    face_pts += 1
                if lhs == rhs:
                    eq_pts += 1
                    if not fc:
                        eq_face_mismatch.append((P1, P2, Q1, Q2))
                elif fc:
                    eq_face_mismatch.append((P1, P2, Q1, Q2))

check("dense grid 17^4 = %d pts (step 1/16): ZERO violations (%d found)"
      % (17 ** 4, viol), viol == 0)
check("dense grid: closed form == raw definition at every fused pair "
      "(%d mismatches)" % closed_mism, closed_mism == 0)
check("dense grid: equality points (%d) == predicted rigid face points (%d), "
      "pointwise both directions (%d mismatches)"
      % (eq_pts, face_pts, len(eq_face_mismatch)),
      eq_pts == face_pts and not eq_face_mismatch)

# Random hostile scan with ragged denominators (exact Fractions).
random.seed(12592)
rviol = 0
req_mismatch = 0
req_eq = 0
NRAND = 20000
for _ in range(NRAND):
    pts = []
    for __ in range(4):
        den = random.randint(1, 64)
        pts.append(F(random.randint(0, den), den))
    P1, P2, Q1, Q2 = pts
    lhs = tv_def((P1 + P2) / 2, (Q1 + Q2) / 2)
    D1 = abs(P1 - Q1)
    D2 = abs(P2 - Q2)
    rhs = D1 + D2 - D1 * D2
    if lhs > rhs:
        rviol += 1
    if lhs == rhs:
        req_eq += 1
        if not on_face(P1, P2, Q1, Q2):
            req_mismatch += 1
check("random hostile scan %d pts (denominators <= 64): zero violations "
      "(%d found), all %d equality hits on the rigid face (%d off-face)"
      % (NRAND, rviol, req_eq, req_mismatch),
      rviol == 0 and req_mismatch == 0)

print()
print("=" * 78)
print("PART 4: transfer probe numbers (AMM 12592 lane, HYP-9061 sec 2e)")
print("=" * 78)
qA = F(896, 2181)              # certificate bias A  (s_A = 7: 896 = 2^7*7)
qB = F(2974400, 11821757)      # certificate bias B  (s_B = 6: 2^6*5^2*11*13^2)
delta = abs(qA - qB)
tvAB = tv_closed(qA, qB)
fused_bound = 2 * delta - delta * delta
print("q_A = %s  (~%.6f)" % (qA, float(qA)))
print("q_B = %s  (~%.6f)" % (qB, float(qB)))
print("delta = |q_A - q_B| = %s  (~%.6f)" % (delta, float(delta)))
print("TV(Bin(2,q_A), Bin(2,q_B)) = %s  (~%.6f)" % (tvAB, float(tvAB)))
print("fusion RHS at p_i=q_A, q_i=q_B: 1-(1-delta)^2 = %s  (~%.6f)"
      % (fused_bound, float(fused_bound)))
check("lemma instance at the certificate biases holds (sanity, exact)",
      tvAB <= fused_bound)
print("""These are dimensionless O(1) rationals with NO n-dependence and no
valuation content: the lane's rate inequality (27) couples an archimedean
log-rate against 2-adic leading terms (s_A=7 vs s_B=6 misalignment); TV of
Bin(2,.) is invariant under outcome relabeling and blind to denominators
and 2-adic valuations of leaf probabilities.  See the note for the
six-field NO TRANSFER verdict.""")

print()
print("=" * 78)
print("OVERALL: %s" % ("ALL CHECKS PASS" if PASS_ALL else "SOME CHECKS FAILED"))
print("=" * 78)
if not PASS_ALL:
    raise SystemExit(1)
