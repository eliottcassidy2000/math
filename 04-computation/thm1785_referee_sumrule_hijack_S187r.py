#!/usr/bin/env python3
"""
thm1785_referee_sumrule_hijack_S187r.py — HOSTILE REFEREE checks for
THM-1785 sections 1-2 (PS-residue lemma, level-set sum rule, mu_g-rigidity,
L1-tube measurement). Referee session S187r. Pure python.

R1  PS-residue lemma at a NONZERO fold (the file's S1 model has PS == 0
    identically — degenerate check; the file's S2 never computed a
    collision limit). Here: cubic model, three independent computations:
      (a) PS formula 2v/(uc^2 L'') + (2/3) v L'''/(uc L''^2)
      (b) contour residue of v/(u(v-Lambda)) around uc (trapezoid)
      (c) collision limit of the pair-sum of residues of 1/(u(1-tL)),
          approached from TWO directions (real eps, imaginary eps).
R2  Level-set sum rule:
      (a) orientation DERIVED: Res_{u*} v/(u(v-L)) = -v/(u* L'(u*))
          (verified by contour) hence Sum PS = + Sum v/(u* L'(u*)).
      (b) all-3-folds check on the cubic model (file did fold[0] only).
      (c) TWO folds sharing one value (symmetric quartic) + 2 spectators.
      (d) global annulus check: sum of ALL residues = 0 (big minus small
          circle).
      (e) K- = 0 counterexample: rule as stated FAILS by exactly the u=0
          residue v/(v - Lambda(0)) — "mixed" hypothesis is load-bearing.
      (f) CUSP on the level (L'=L''=0): rule as stated is silent/violated;
          the cusp's own residue (contour) restores the identity — the
          statement needs the only-folds-and-simple-points hypothesis.
R3  mu_g-rotation rigidity:
      (a) exact two-term edge (4,-2): sharing pattern = mu_2 orbits, all
          sharers same modulus (rotation) — CONFIRM the scoped lemma.
      (b) CROSS-EDGE value sharing witness: palindromic
          Lambda = q0 + w^2 u^2 + w u + w/u + w^2/u^2. The two end-edge
          folds u+ ~ -1/(2w), u- ~ -2w satisfy u+ u- = 1 and share the
          value EXACTLY at every w, but |u+| != |u-|: NOT rotation-related.
          => "deleted hijacked stacks are (mu_g) symmetry stacks" is NOT
          proved for stacks meeting more than one edge.
R4  L1 tube on P4 = ZW + Z^9W^7 + W: delta-sweep at FIXED s0 (the file
    fixed delta = 3e-3 and swept s0 only). Graded s-grid resolving the
    spike; measures I(t) = int e^{-s}|Ghat| ds, the signed |A_win|, the
    peak, and the peak's delta-scaling exponent (pole vs sqrt).

boxeph-referee, 2026-07-20.
"""
import cmath
import math


# ---------------- generic root machinery (Durand-Kerner) ----------------

def droots(coeffs, seeds=None, iters=300):
    c0 = coeffs[0]
    cs = [complex(c) / c0 for c in coeffs]
    n = len(cs) - 1

    def p(x):
        v = complex(0)
        for c in cs:
            v = v * x + c
        return v
    if seeds is None:
        seeds = [(0.4 + 0.9j) ** (k + 1) * (1.3 + 0.4j) for k in range(n)]
    rs = list(seeds)
    for _ in range(iters):
        moved = 0.0
        for i in range(n):
            d = complex(1)
            for j in range(n):
                if j != i:
                    d *= (rs[i] - rs[j])
            if d == 0:
                d = 1e-30
            step = p(rs[i]) / d
            rs[i] -= step
            moved = max(moved, abs(step))
        if moved < 1e-14:
            break
    return rs


def contour_residue(F, center, radius, N=4096):
    # (1/2pi i) oint F du over circle center+radius e^{i th}
    tot = complex(0)
    for k in range(N):
        th = 2 * math.pi * k / N
        z = center + radius * cmath.exp(1j * th)
        tot += F(z) * radius * cmath.exp(1j * th)
    return tot / N


def PSformula(v, uc, L2, L3):
    return 2 * v / (uc * uc * L2) + (2.0 / 3) * v * L3 / (uc * L2 ** 2)


print("=" * 78)
print("R1: PS-residue lemma at a NONZERO fold (three independent computations)")
print("=" * 78)
# cubic model (same as the file's S2): Lambda = q0 s + w(p1 u + p2 u^2 + pm1/u)
q0c, p1, pm1, p2c = 0.4 + 0.1j, 1.0 + 0.3j, 0.8 - 0.2j, 0.5 + 0.25j
s = 2.0
w = math.sqrt(2 * s)
lam = lambda u: q0c * s + w * (p1 * u + p2c * u * u + pm1 / u)
lam1 = lambda u: w * (p1 + 2 * p2c * u - pm1 / (u * u))
lam2 = lambda u: w * (2 * p2c + 2 * pm1 / u ** 3)
lam3 = lambda u: -6 * w * pm1 / u ** 4
# critical points: roots of 2 p2 u^3 + p1 u^2 - pm1 = 0
crits = droots([2 * p2c, p1, 0, -pm1])
uc = crits[0]
v = lam(uc)
PSa = PSformula(v, uc, lam2(uc), lam3(uc))
print("fold uc = %.6f%+.6fj, v = %.6f%+.6fj" % (uc.real, uc.imag, v.real, v.imag))
print("(a) PS formula                    = %.9f%+.9fj" % (PSa.real, PSa.imag))
# (b) contour residue of v/(u(v-Lambda)) at uc (double pole)
lset = droots([w * p2c, w * p1, q0c * s - v, w * pm1])  # roots of v-Lambda (as cubic *u)
others = sorted(lset, key=lambda r: abs(r - uc))[2:]    # non-colliding root(s)
mind = min(abs(uc - o) for o in others)
F = lambda u: v / (u * (v - lam(u)))
for rad in (0.25 * mind, 0.125 * mind):
    Rb = contour_residue(F, uc, rad)
    print("(b) contour residue (r=%.4f)     = %.9f%+.9fj" % (rad, Rb.real, Rb.imag))
# (c) collision limit of pair-sum of residues of 1/(u(1-tLambda)), 2 directions
for direction, label in ((1.0, "real eps"), (1j, "imag eps")):
    for eps in (1e-3, 1e-5, 1e-7):
        t = (1 - direction * eps) / v
        # u(1-tLambda) = -t w p2 u^3 - t w p1 u^2 + (1-t q0 s) u - t w pm1
        rts = droots([-t * w * p2c, -t * w * p1, 1 - t * q0c * s, -t * w * pm1])
        rts = sorted(rts, key=lambda r: abs(r - uc))
        pair = rts[:2]
        S = sum(1 / (r * (-t) * lam1(r)) for r in pair)
        print("(c) pair-sum  %s eps=%.0e  = %.9f%+.9fj" %
              (label, eps, S.real, S.imag))

print()
print("=" * 78)
print("R2: level-set sum rule — orientation, multi-fold, K-=0, cusp")
print("=" * 78)
# (a) orientation: contour residue at the SIMPLE spectator root
ustar = others[0]
Rsimple = contour_residue(F, ustar, 0.2 * min(abs(ustar - r) for r in lset if abs(r - ustar) > 1e-9), 4096)
pred = -v / (ustar * lam1(ustar))
print("(a) Res at simple u*: contour = %.9f%+.9fj ; -v/(u*L') = %.9f%+.9fj"
      % (Rsimple.real, Rsimple.imag, pred.real, pred.imag))
print("    => residue theorem gives Sum PS_folds = + Sum v/(u* L'(u*)) (equality,")
print("       matching THM-1785; the .py caption 'PS + spectator = 0' is wrong).")

# (b) all three folds of the cubic model
print("(b) all 3 folds of the cubic model:")
for i, uci in enumerate(crits):
    vi = lam(uci)
    PSi = PSformula(vi, uci, lam2(uci), lam3(uci))
    lseti = droots([w * p2c, w * p1, q0c * s - vi, w * pm1])
    lseti = sorted(lseti, key=lambda r: -abs(r - uci))
    usti = lseti[0]
    RHS = vi / (usti * lam1(usti))
    print("    fold %d: PS = %.8f%+.8fj  RHS = %.8f%+.8fj  |PS-RHS| = %.2e"
          % (i, PSi.real, PSi.imag, RHS.real, RHS.imag, abs(PSi - RHS)))

# (c) TWO folds sharing one value: symmetric quartic
A4, C4, B4 = 1.0 + 0.2j, 0.7 - 0.1j, 0.9 + 0.4j
q0d, sd = 0.3 - 0.2j, 1.7
wd = math.sqrt(2 * sd)
lamQ = lambda u: q0d * sd + wd * (A4 * u * u + C4 * u ** 4 + B4 / (u * u))
lamQ1 = lambda u: wd * (2 * A4 * u + 4 * C4 * u ** 3 - 2 * B4 / u ** 3)
lamQ2 = lambda u: wd * (2 * A4 + 12 * C4 * u * u + 6 * B4 / u ** 4)
lamQ3 = lambda u: wd * (24 * C4 * u - 24 * B4 / u ** 5)
xs = droots([4 * C4, 2 * A4, 0, -2 * B4])  # 4C x^3 + 2A x^2 - 2B = 0 in x=u^2
xc = xs[0]
ucq = cmath.sqrt(xc)
vq = lamQ(ucq)
PSp = PSformula(vq, ucq, lamQ2(ucq), lamQ3(ucq))
PSm = PSformula(vq, -ucq, lamQ2(-ucq), lamQ3(-ucq))
# level set: C x^3 + A x^2 + ((q0 s - v)/w) x + B = 0 in x = u^2
lx = droots([C4, A4, (q0d * sd - vq) / wd, B4])
lx = sorted(lx, key=lambda r: -abs(r - xc))
xstar = lx[0]
ustq = cmath.sqrt(xstar)
RHS2 = vq / (ustq * lamQ1(ustq)) + vq / ((-ustq) * lamQ1(-ustq))
print("(c) 2 folds (+-uc) sharing v + 2 spectators (+-u*):")
print("    PS(+uc)+PS(-uc) = %.8f%+.8fj ; Sum v/(u*L') = %.8f%+.8fj ; diff %.2e"
      % ((PSp + PSm).real, (PSp + PSm).imag, RHS2.real, RHS2.imag,
         abs(PSp + PSm - RHS2)))

# (d) global annulus check on the quartic model: total residue = 0
FQ = lambda u: vq / (u * (vq - lamQ(u)))
big = contour_residue(FQ, 0.0, 30.0, 8192)
small = contour_residue(FQ, 0.0, 0.01, 8192)
print("(d) annulus total (big R=30 minus small r=0.01 circle) = %.2e (should be ~0)"
      % abs(big - small))

# (e) K- = 0: Lambda = q0 s + w(p1 u + p2 u^2 + p3 u^3): pole at u=0 appears
p3e = 0.6 - 0.15j
lamE = lambda u: q0c * s + w * (p1 * u + p2c * u * u + p3e * u ** 3)
lamE1 = lambda u: w * (p1 + 2 * p2c * u + 3 * p3e * u * u)
lamE2 = lambda u: w * (2 * p2c + 6 * p3e * u)
lamE3 = lambda u: 6 * w * p3e
ce = droots([3 * p3e, 2 * p2c, p1])
uce = ce[0]
ve = lamE(uce)
PSe = PSformula(ve, uce, lamE2(uce), lamE3(uce))
lroots_e = droots([w * p3e, w * p2c, w * p1, q0c * s - ve])
lroots_e = sorted(lroots_e, key=lambda r: -abs(r - uce))
uste = lroots_e[0]
RHSe = ve / (uste * lamE1(uste))
res0 = ve / (ve - q0c * s)   # residue of v/(u(v-L)) at u=0 when K-=0
print("(e) K-=0 model: PS = %.6f%+.6fj ; RHS-as-stated = %.6f%+.6fj" %
      (PSe.real, PSe.imag, RHSe.real, RHSe.imag))
print("    rule-as-stated violation |PS - RHS| = %.3e" % abs(PSe - RHSe))
print("    u=0 residue v/(v-L(0)) = %.6f%+.6fj ; |PS - (RHS - res0)| = %.3e"
      % (res0.real, res0.imag, abs(PSe - (RHSe - res0))))
print("    => 'mixed' (K-,K+ >= 1) is LOAD-BEARING; failure = exactly the u=0 residue")

# (f) cusp on the level: Lambda = q0 s + w(u^3 + 3u - (3/4)/u); L'=3w(u^2+1/2)^2/u^2
lamC = lambda u: q0c * s + w * (u ** 3 + 3 * u - 0.75 / u)
lamC1 = lambda u: w * (3 * u * u + 3 + 0.75 / (u * u))
ucc = 1j / math.sqrt(2)
vc = lamC(ucc)
print("(f) cusp model: uc = i/sqrt2, L'(uc) = %.2e, v = %.6f%+.6fj"
      % (abs(lamC1(ucc)), vc.real, vc.imag))
# level set: -w u^4 - 3w u^2 + (v - q0 s) u + (3/4) w = 0: triple root + 1 simple
lc = droots([-w, 0, -3 * w, vc - q0c * s, 0.75 * w])
lc = sorted(lc, key=lambda r: -abs(r - ucc))
ustc = lc[0]
RHSc = vc / (ustc * lamC1(ustc))
FC = lambda u: vc / (u * (vc - lamC(u)))
Rcusp = contour_residue(FC, ucc, 0.25 * abs(ustc - ucc), 8192)
print("    simple-side sum v/(u*L') = %.8f%+.8fj (rule-as-stated LHS = 0: VIOLATED)"
      % (RHSc.real, RHSc.imag))
print("    cusp contour residue = %.8f%+.8fj ; |Res_cusp - v/(u*L')| = %.2e"
      % (Rcusp.real, Rcusp.imag, abs(Rcusp - RHSc)))
print("    => sum rule needs the hypothesis: level set = simple folds + simple")
print("       points ONLY (cusp levels need the cusp residue term).")

print()
print("=" * 78)
print("R3: mu_g rigidity — scoped confirmation + CROSS-EDGE sharing witness")
print("=" * 78)
# (a) exact two-term edge d1=4, d2=-2 (g=2): mu_6 orbit, values by zeta^4
q1e, q2e, p0e = 0.8 + 0.3j, 1.1 - 0.5j, 0.25 + 0.6j
lamT = lambda u: p0e + q1e * u ** 4 + q2e / (u * u)
# crit: 4 q1 u^3 - 2 q2/u^3 = 0 -> u^6 = q2/(2 q1)
u6 = q2e / (2 * q1e)
ucT = u6 ** (1.0 / 6)
orb = [ucT * cmath.exp(1j * math.pi * k / 3) for k in range(6)]
vals = [lamT(u) for u in orb]
print("(a) edge (4,-2), g=gcd(4,2)=2: critical orbit values (should pair by u->-u):")
shared = []
for i in range(6):
    for j in range(i + 1, 6):
        if abs(vals[i] - vals[j]) < 1e-10 * max(1.0, abs(vals[i])):
            shared.append((i, j, abs(orb[i] / orb[j] + 1) < 1e-9))
print("    sharing pairs (i,j, ratio==-1?): %s" % shared)
print("    all sharers same modulus (rotation): %s"
      % all(abs(abs(orb[i]) - abs(orb[j])) < 1e-12 for i, j, _ in shared))

# (b) cross-edge witness: palindromic Lambda
print("(b) CROSS-EDGE witness: Lambda = q0 + w^2 u^2 + w u + w/u + w^2/u^2")
for wv in (0.4, 0.15, 0.05):
    q0p = 0.7 + 0.2j
    lamP = lambda u: q0p + wv ** 2 * u * u + wv * u + wv / u + wv ** 2 / (u * u)
    # critical points: (u^2-1)(2w u^2 + u + 2w) = 0 -> u = +-1 and inversion pair
    disc = cmath.sqrt(1 - 16 * wv * wv)
    up = (-1 + disc) / (4 * wv)
    um = (-1 - disc) / (4 * wv)
    # which is which: |up| small? assign big/small
    ubig, usmall = (up, um) if abs(up) > abs(um) else (um, up)
    dv = abs(lamP(ubig) - lamP(usmall))
    print("    w=%.2f: u+ = %.6f (edge {1,2} fold), u- = %.6f (edge {-1,-2} fold)"
          % (wv, ubig.real, usmall.real))
    print("           u+ * u- = %.12f ; |L(u+)-L(u-)| = %.2e ; |u+|/|u-| = %.3f"
          % ((ubig * usmall).real, dv, abs(ubig) / abs(usmall)))
print("    => EXACT value sharing between folds on DIFFERENT edges, related by")
print("       INVERSION u -> 1/u (moduli differ: NOT a mu_g rotation).")
print("       'Deleted hijacked stacks are symmetry stacks (mu_g)' is proved ONLY")
print("       for events on ONE two-term edge; cross-edge sharing is real.")

print()
print("=" * 78)
print("R4: L1 tube on P4 — delta-sweep at fixed s0 (the file fixed delta=3e-3)")
print("=" * 78)


def lamP4(u, s):
    w = math.sqrt(2 * s)
    return 2 * s + (w ** 16) * u * u + w / u


def lamP4_1(u, s):
    w = math.sqrt(2 * s)
    return 2 * (w ** 16) * u - w / (u * u)


def cubic_coeffs_P4(t, s):
    w = math.sqrt(2 * s)
    return [-t * w ** 16, 0.0, 1 - 2 * t * s, -t * w]


def u0_by_t_continuation(s, t, nsteps=200):
    # the file's method (0-group root by homotopy in t), returns the root
    t0 = 1e-6 * abs(t)
    cur = sorted(droots(cubic_coeffs_P4(t0, s)), key=abs)
    for k in range(1, nsteps + 1):
        tt = t0 + (t - t0) * (k / nsteps) ** 1.3
        new = droots(cubic_coeffs_P4(tt, s), seeds=cur, iters=60)
        matched = [None] * 3
        used = set()
        for i, p in enumerate(cur):
            bi, bd = None, None
            for j, r in enumerate(new):
                if j in used:
                    continue
                d = abs(r - p)
                if bd is None or d < bd:
                    bd, bi = d, j
            used.add(bi)
            matched[i] = new[bi]
        cur = matched
    return cur[0]


def v_of(s):
    w = math.sqrt(2 * s)
    ucp = (w / (2 * w ** 16)) ** (1.0 / 3)
    return lamP4(ucp, s)


def graded_grid(s0):
    # descending s grid: far tail 8 -> 3s0; log-graded window 3s0 -> s0(1+1e-6),
    # cross gap, s0(1-1e-6) -> s0/3; small tail s0/3 -> s0*1e-2
    pts = []
    n1 = 30
    for i in range(n1 + 1):
        pts.append(3 * s0 * (8.0 / (3 * s0)) ** (1 - i / n1))
    # window upper: s-s0 from 2s0 down to 1e-6 s0, 16/decade
    dec = math.log10(2.0 / 1e-6)
    n2 = int(dec * 16)
    for i in range(n2 + 1):
        x = 2.0 * (1e-6 / 2.0) ** (i / n2)
        pts.append(s0 * (1 + x))
    for i in range(n2 + 1):
        x = 1e-6 * (2.0 / 3.0 / 1e-6) ** (i / n2)   # up to 2/3 s0 below
        pts.append(s0 * (1 - x))
    n3 = 15
    for i in range(1, n3 + 1):
        pts.append((s0 / 3.0) * (1e-2 / (1.0 / 3.0)) ** (i / n3) * s0 / s0)
    pts = sorted(set(pts), reverse=True)
    return pts


def sweep(s0, deltas):
    print("  s0 = %.3f  (v0 = %.6e, arc point |t| = %.4e)" %
          (s0, abs(v_of(s0)), abs(1 / v_of(s0))))
    for de in deltas:
        v0 = v_of(s0)
        t = (1 + 1j * de) / v0
        grid = graded_grid(s0)
        # initialize tracking at the far end by the file's own t-continuation
        cur = None
        u0i = None
        Iwin = 0.0
        Awin = complex(0)
        Itail = 0.0
        Ilow = 0.0
        peak = 0.0
        prev_s = None
        prev_f = None
        prev_absf = None
        hop_flag = 0
        for si in grid:
            if cur is None:
                u0 = u0_by_t_continuation(si, t)
                allr = droots(cubic_coeffs_P4(t, si))
                allr = sorted(allr, key=lambda r: abs(r - u0))
                cur = allr
                u0i = 0
            else:
                new = droots(cubic_coeffs_P4(t, si), seeds=cur, iters=80)
                matched = [None] * 3
                used = set()
                for i, p in enumerate(cur):
                    bi, bd = None, None
                    for j, r in enumerate(new):
                        if j in used:
                            continue
                        d = abs(r - p)
                        if bd is None or d < bd:
                            bd, bi = d, j
                    used.add(bi)
                    matched[i] = new[bi]
                # hop diagnostic: matching step vs root separation
                sep = min(abs(matched[i] - matched[j])
                          for i in range(3) for j in range(i + 1, 3))
                stepd = max(abs(matched[i] - cur[i]) for i in range(3))
                if stepd > 0.6 * sep:
                    hop_flag += 1
                cur = matched
            u0 = cur[u0i]
            gh = -1 / (t * u0 * lamP4_1(u0, si))
            f = math.exp(-si) * gh
            af = abs(f)
            aG = abs(gh)
            if abs(si - s0) < 1e-3 * s0 and aG > peak:
                peak = aG
            if prev_s is not None:
                seg = abs(prev_s - si) * 0.5 * (af + prev_absf)
                segA = (prev_s - si) * 0.5 * (f + prev_f)
                lo, hi = min(prev_s, si), max(prev_s, si)
                if hi <= 3 * s0 + 1e-15 and lo >= s0 / 3 - 1e-15:
                    Iwin += seg
                    Awin += segA
                elif lo >= 3 * s0 - 1e-15:
                    Itail += seg
                else:
                    Ilow += seg
            prev_s, prev_f, prev_absf = si, f, af
        print("    delta=%.0e  peak|Ghat|=%.4e  I_win=%.4f  I_tail=%.4f "
              "I_low=%.4f  I_tot=%.4f  |A_win|=%.4f  hops=%d"
              % (de, peak, Iwin, Itail, Ilow, Iwin + Itail + Ilow,
                 abs(Awin), hop_flag))
    return


deltas = (3e-3, 1e-3, 3e-4, 1e-4, 3e-5)
for s0 in (0.1, 0.03):
    sweep(s0, deltas)
print("\n  interpretation aid: if I_win(delta) ~ a + b*log(1/delta) the tube is NOT")
print("  L1-uniform as delta->0; if I_win converges (sqrt-type fold) it is.")
print("  peak scaling: pole => peak ~ 1/delta; sqrt-branch => peak ~ 1/sqrt(delta).")

print("\nDONE.")
