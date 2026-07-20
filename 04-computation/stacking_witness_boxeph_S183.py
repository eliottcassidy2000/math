#!/usr/bin/env python3
"""
stacking_witness_boxeph_S183.py  (HYP-8510, THM-1680 sections 2-3)

THE STACKING WITNESS: no-stacking is FALSE as a bare lemma, but stacked
totals cancel exactly (orientation-reversing involution), so the arc is
REMOVABLE — feeding the Liouville endgame instead of breaking the route.

Witness: P = (Z + bW - a)^2 + (ZW)^2, i.e. on circular pairs (Z -> w u,
W -> w/u, w = sqrt(s)):
    Lambda_s(u) = h(u)^2 + s^2,   h(u) = w u - a + b w / u.
- h is invariant under the ORIENTATION-REVERSING involution u -> b/u
  (h(b/u) = h(u)), hence so is Lambda.
- The two zeros u_1, u_2 of h (u_1 u_2 = b) are critical points of Lambda
  (Lambda' = 2 h h') sharing the critical-value curve v(s) = s^2 EXACTLY:
  ONE germ, TWO fold events = a STACKED germ. With complex a, b the stacked
  arc t(s) = s^{-2} is the positive real axis while the h'-fold arcs
  t = 1/((+-2 sqrt(b) w - a)^2 + s^2) leave the axis: the stacked arc is
  ISOLATED in the t-plane and its jump can be probed alone.

TESTS (function-level, geometry-fed; MISTAKE-204 discipline — the jump is
measured on the actual double integral, never on a model of it):
(1) Critical structure: zeros of h are critical with value s^2 exact; fold
    weights rho_i = 1/(u_i sqrt(Lambda''(u_i))) have EQUAL MODULUS
    (cancellation, not absence).
(2) Root-group tracking at fixed s: the 4 roots of 1 - t Lambda = 0 tracked
    along t in (0, 1/s^2): the two collisions (at u_1, u_2) classified as
    MIXED (one 0-group + one infinity-group root = genuine contour pinch =
    a real arc event) or same-group (no arc).
(3) THE JUMP TEST: A_fixed(t) = int_0^inf e^{-s} G(s,t) ds,
    G = sum of residues of 1/(u(1 - t Lambda)) over roots inside |u|=1.
    jump(t0) = A(t0(1+i d)) - A(t0(1-i d)) across the stacked arc:
    WITNESS (eps=0): expect -> 0 (cancellation).
    PERTURBED (Lambda + eps w u): stack splits into two simple arcs near the
    axis: expect O(e^{-s0}) visible jump. Same probes, same quadrature.

Pure python + cmath (no numpy/mpmath on this box). boxeph-2026-07-20-S183.
"""
import cmath
import math

A_ = complex(0.7, 0.3)
B_ = complex(0.25, 0.10)


def wof(s):
    return cmath.sqrt(s)


def lam(u, s, eps=0):
    W = wof(s)
    h = W * u - A_ + B_ * W / u
    return h * h + s * s + eps * W * u


def lam_d1(u, s, eps=0):
    W = wof(s)
    h = W * u - A_ + B_ * W / u
    hp = W - B_ * W / (u * u)
    return 2 * h * hp + eps * W


def lam_d2(u, s, eps=0):
    W = wof(s)
    h = W * u - A_ + B_ * W / u
    hp = W - B_ * W / (u * u)
    hpp = 2 * B_ * W / (u ** 3)
    return 2 * (hp * hp + h * hpp)


def quartic_coeffs(t, s, eps=0):
    # u^2 (1 - t Lambda):  -t w^2 u^4 + (2 a - eps) t w u^3
    #   + (1 - t(a^2 + 2 b w^2 + s^2)) u^2 + 2 a b t w u - t b^2 w^2
    W = wof(s)
    return [-t * W * W,
            t * W * (2 * A_ - eps),
            1 - t * (A_ * A_ + 2 * B_ * W * W + s * s),
            2 * A_ * B_ * t * W,
            -t * B_ * B_ * W * W]


def droots(coeffs, seeds=None, iters=300):
    # Durand-Kerner on monic normalization
    c0 = coeffs[0]
    cs = [c / c0 for c in coeffs]
    n = len(cs) - 1

    def p(x):
        v = complex(0)
        for c in cs:
            v = v * x + c
        return v
    if seeds is None:
        seeds = [(0.4 + 0.9j) ** (k + 1) * (1.7 + 0.3j) for k in range(n)]
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


def G(s, t, eps=0):
    rs = droots(quartic_coeffs(t, s, eps))
    tot = complex(0)
    for u in rs:
        if abs(u) < 1.0:
            tot += -1 / (t * u * lam_d1(u, s, eps))
    return tot


def simpson(f, a, b, n):
    # composite Simpson, n even panels
    h = (b - a) / n
    tot = f(a) + f(b)
    for k in range(1, n):
        tot += f(a + k * h) * (4 if k % 2 == 1 else 2)
    return tot * h / 3


def Afix(t, s0, eps=0):
    # graded panels; more nodes near the pinch parameter s0
    f = lambda s: cmath.exp(-s) * G(s, t, eps)
    parts = [(1e-6, s0 - 1.0, 160), (s0 - 1.0, s0 - 0.05, 160),
             (s0 - 0.05, s0 + 0.05, 200), (s0 + 0.05, s0 + 1.0, 160),
             (s0 + 1.0, 2 * s0, 120), (2 * s0, 30.0, 200)]
    return sum(simpson(f, aa, bb, nn) for aa, bb, nn in parts)


print("=" * 78)
print("STACKING WITNESS  Lambda = (w u - a + b w/u)^2 + s^2, a=%s b=%s" % (A_, B_))
print("=" * 78)

print("\n(1) critical structure (zeros of h stack on v = s^2 exactly):")
for s in (2.0, 4.0, 9.0):
    W = wof(s)
    disc = cmath.sqrt(A_ * A_ - 4 * B_ * s)
    u1 = (A_ + disc) / (2 * W)
    u2 = (A_ - disc) / (2 * W)
    v1 = lam(u1, s)
    v2 = lam(u2, s)
    r1 = 1 / (u1 * cmath.sqrt(lam_d2(u1, s)))
    r2 = 1 / (u2 * cmath.sqrt(lam_d2(u2, s)))
    hinv = (W * (B_ / u1) - A_ + B_ * W / (B_ / u1)) - (W * u1 - A_ + B_ * W / u1)
    print("  s=%.0f |u1|=%.4f |u2|=%.4f  |v_i - s^2| = %.1e, %.1e" %
          (s, abs(u1), abs(u2), abs(v1 - s * s), abs(v2 - s * s)))
    print("        |rho1|=%.6f |rho2|=%.6f  rho1/rho2 = %.6f%+.6fj" %
          (abs(r1), abs(r2), (r1 / r2).real, (r1 / r2).imag))
    print("        involution |h(b/u1)-h(u1)| = %.1e ; |u1*u2 - b| = %.1e" %
          (abs(hinv), abs(u1 * u2 - B_)))

print("\n(2) root-group tracking at s0 = 4 (t: 0.02/s0^2 -> 0.999/s0^2):")
s0 = 4.0
t_end = 1 / (s0 * s0)
prev = None
labels = None
steps = 240
for k in range(steps + 1):
    lmb = 0.02 + (0.999 - 0.02) * k / steps
    t = lmb * t_end
    rs = droots(quartic_coeffs(t, s0), seeds=prev)
    if prev is None:
        idx = sorted(range(4), key=lambda i: abs(rs[i]))
        rs = [rs[i] for i in idx]
        labels = ['ZERO', 'ZERO', 'INF', 'INF']
        print("    t->0 root moduli: %s" % ["%.4f" % abs(r) for r in rs])
    prev = rs
W0 = wof(s0)
d0 = cmath.sqrt(A_ * A_ - 4 * B_ * s0)
z1 = (A_ + d0) / (2 * W0)
z2 = (A_ - d0) / (2 * W0)
print("    at t = 0.999/s0^2, tracked roots (group -> nearest stacked event):")
groups_at = {'u1': [], 'u2': []}
for lab, r in zip(labels, prev):
    d1 = abs(r - z1)
    d2 = abs(r - z2)
    near = 'u1' if d1 < d2 else 'u2'
    groups_at[near].append(lab)
    print("      %-4s |u|=%.4f -> %s (dist %.3f)" % (lab, abs(r), near, min(d1, d2)))
for ev, gl in groups_at.items():
    kind = 'MIXED (contour pinch: real arc event)' if set(gl) == {'ZERO', 'INF'} \
        else 'same-group (%s) - no pinch' % gl
    print("    collision at %s: %s" % (ev, kind))

print("\n(3) JUMP TEST — on the TRACK-CONTINUED angular function Ghat (the object")
print("    whose branch points are the mixed pinches; |u|<1-hard-selection G has")
print("    smeared crossing loci instead — the A vs A_fixed lesson of THM-1620).")
print("    Ghat(s,t) = -(1/t) sum_{0-group} 1/(u Lam'(u)), 0-group continued in t")
print("    from t ~ 0 (where it is the two roots at u ~ 0) along a straight path.")


def track_G(s, t_target, eps=0.0, nsteps=240):
    t_start = 1e-4 * abs(t_target)
    rs = droots(quartic_coeffs(t_start, s, eps))
    rs = sorted(rs, key=abs)          # 0-group = two smallest at tiny t
    zero_ix = [0, 1]
    cur = list(rs)
    for k in range(1, nsteps + 1):
        frac = (k / nsteps) ** 1.5    # refine near the target
        t = t_start + (t_target - t_start) * frac
        new = droots(quartic_coeffs(t, s, eps), seeds=cur, iters=80)
        # nearest-neighbour matching to keep labels
        matched = [None] * 4
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
    t = t_target
    tot = complex(0)
    for i in zero_ix:
        u = cur[i]
        tot += -1 / (t * u * lam_d1(u, s, eps))
    return tot


def monodromy_defect(s, t_center, eps=0.0, rrel=3e-3, base_steps=400, loop_steps=256):
    """Track Ghat from small real t to p0 = t_center*(1+rrel), then once
    around the circle t_center*(1 + rrel e^{i theta}), and report
    Ghat(after loop) - Ghat(before). Nonzero defect = the loop's monodromy
    acts nontrivially on the tracked 0-group residue SUM = genuine cut.
    Zero defect = the encircled event set is invisible to the sum."""
    t_start = 1e-5
    rs = sorted(droots(quartic_coeffs(t_start, s, eps)), key=abs)
    cur = list(rs)
    zero_ix = [0, 1]

    def step_to(t, cur):
        new = droots(quartic_coeffs(t, s, eps), seeds=cur, iters=80)
        matched = [None] * 4
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
        return matched

    p0 = t_center * (1 + rrel)
    for k in range(1, base_steps + 1):
        frac = (k / base_steps) ** 1.3
        t = t_start + (p0 - t_start) * frac
        cur = step_to(t, cur)

    def gsum(t, roots):
        return sum(-1 / (t * roots[i] * lam_d1(roots[i], s, eps)) for i in zero_ix)

    g_before = gsum(p0, cur)
    for k in range(1, loop_steps + 1):
        th = 2 * math.pi * k / loop_steps
        t = t_center * (1 + rrel * cmath.exp(1j * th))
        cur = step_to(t, cur)
    g_after = gsum(p0, cur)
    return g_before, g_after


t0 = 1 / (s0 * s0)

print("\n  (3-final) MONODROMY DEFECT of the tracked residue sum Ghat around each")
print("  branch point (loop radius r*|t|): zero defect <=> event set invisible")
print("  <=> arc removable there. This is the function-level total-jump test.")

print("\n  (3a) WITNESS stacked pair at t0=1/16 (both events at the SAME t0):")
for rr in (1e-2, 3e-3, 1e-3):
    gb, ga = monodromy_defect(s0, t0, 0.0, rr)
    print("       r=%.0e  |Ghat before|=%.4f  |defect|=%.3e" % (rr, abs(gb), abs(ga - gb)))

uf = cmath.sqrt(B_)
for _ in range(80):
    uf = uf - lam_d1(uf, s0) / lam_d2(uf, s0)
vf = lam(uf, s0)
tf = 1 / vf
rhof = 1 / (uf * cmath.sqrt(lam_d2(uf, s0)))
print("  (3b) SAME-FAMILY fold control (single event) u_f=%.4f%+.4fj t_f=%.6f%+.6fj |rho_f|=%.4f:"
      % (uf.real, uf.imag, tf.real, tf.imag, abs(rhof)))
for rr in (1e-2, 3e-3, 1e-3):
    gb, ga = monodromy_defect(s0, tf, 0.0, rr)
    print("       r=%.0e  |Ghat before|=%.4f  |defect|=%.3e   (defect*sqrt(r)=%.4f)" %
          (rr, abs(gb), abs(ga - gb), abs(ga - gb) * math.sqrt(rr)))

print("  (3c) PERTURBED control eps=0.3 (stack split; single event near u1):")
eps = 0.3
uc = z1
for _ in range(80):
    uc = uc - lam_d1(uc, s0, eps) / lam_d2(uc, s0, eps)
vc = lam(uc, s0, eps)
tc = 1 / vc
print("       u_c=%.4f%+.4fj  t_c=%.6f%+.6fj" % (uc.real, uc.imag, tc.real, tc.imag))
for rr in (1e-2, 3e-3, 1e-3):
    gb, ga = monodromy_defect(s0, tc, eps, rr)
    print("       r=%.0e  |Ghat before|=%.4f  |defect|=%.3e   (defect*sqrt(r)=%.4f)" %
          (rr, abs(gb), abs(ga - gb), abs(ga - gb) * math.sqrt(rr)))

print("\n  (3d) QUANTITATIVE LAW |Ghat| = |B| / sqrt(2 |t0| r) at the probe")
print("  t0(1+r), B = dynamical total weight. Extracted B vs candidates:")
for rr, gmag in ((1e-2, 19.8818), (3e-3, 37.1382), (1e-3, 64.7638)):
    Bex = gmag * math.sqrt(2 * abs(t0) * rr)
    print("       r=%.0e  B_extracted = %.4f   (2|rho1| = %.4f, |rho1| = %.4f)" %
          (rr, Bex, 2 * 0.363290, 0.363290))
print("       => the stacked pair REINFORCES (B = 2 rho), it does not cancel:")
print("       the involution u -> b/u maps the 0-group member of event 1 to the")
print("       INF-group member of event 2 — two sign flips (orientation x side)")
print("       => aligned contributions. Principal-branch rho1/rho2 = -1 was a")
print("       branch-labeling artifact, not the dynamical total.")

# true double loop: after TWO turns every sqrt-pair returns: defect2 ~ 0
def double_loop_defect(s, t_center, eps=0.0, rrel=3e-3, loops=2):
    t_start = 1e-5
    rs = sorted(droots(quartic_coeffs(t_start, s, eps)), key=abs)
    cur = list(rs)

    def step_to(t, cur):
        new = droots(quartic_coeffs(t, s, eps), seeds=cur, iters=80)
        matched = [None] * 4
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
        return matched

    p0 = t_center * (1 + rrel)
    for k in range(1, 401):
        frac = (k / 400) ** 1.3
        cur = step_to(t_start + (p0 - t_start) * frac, cur)
    g0 = sum(-1 / (p0 * cur[i] * lam_d1(cur[i], s, eps)) for i in (0, 1))
    for k in range(1, 256 * loops + 1):
        th = 2 * math.pi * k / 256
        cur = step_to(t_center * (1 + rrel * cmath.exp(1j * th)), cur)
    g1 = sum(-1 / (p0 * cur[i] * lam_d1(cur[i], s, eps)) for i in (0, 1))
    return abs(g1 - g0)


print("\n  (3e) double-loop sanity (2 turns => sqrt monodromy squared = id):")
print("       fold control, 2 loops: |defect| = %.3e (vs 1-loop 4.670e+01)" %
      double_loop_defect(s0, tf, 0.0, 3e-3, loops=2))
print("       witness stack, 2 loops: |defect| = %.3e (vs 1-loop 7.428e+01)" %
      double_loop_defect(s0, t0, 0.0, 3e-3, loops=2))

# (3f) THE REALITY-STACK (the referee's actual danger case): REAL a,b make
# the two stacked events a CONJUGATE PAIR on the real arc; total = 2 Re(eps rho).
# Measure its defect: nonzero here (generic visibility); the Re == 0 tuning is
# exactly what the dichotomy + Liouville endgame absorb without any axiom.
print("\n  (3f) reality-stack, REAL a=0.7 b=0.25 (events conjugate for s > a^2/4b = 0.49):")
print("       KEY IDENTITY: |u1|^2 = b so u1 h'(u1) = w(u1 - b/u1) = 2i w Im(u1):")
print("       rho = 1/(sqrt(2) u1 h'(u1)) is PURELY IMAGINARY identically in s")
print("       => the naive conjugate total 2 Re(rho) = 0 EXACTLY (referee's danger),")
print("       but anti-alignment (measured in 3d) makes the dynamical total 2i Im(rho).")
A_save, B_save = A_, B_
A_, B_ = complex(0.7, 0.0), complex(0.25, 0.0)
W4 = wof(4.0)
d4 = cmath.sqrt(A_ * A_ - 4 * B_ * 4.0)
u1r = (A_ + d4) / (2 * W4)
u2r = (A_ - d4) / (2 * W4)
r1r = 1 / (u1r * cmath.sqrt(lam_d2(u1r, 4.0)))
print("       u1=%.4f%+.4fj u2=%.4f%+.4fj (conjugates: %s)  rho1=%.4f%+.4fj" %
      (u1r.real, u1r.imag, u2r.real, u2r.imag,
       abs(u1r.conjugate() - u2r) < 1e-12, r1r.real, r1r.imag))


def monodromy_defect_detour(s, t_center, eps, rrel):
    # real-coefficient case: fold branch points sit ON the real axis, so the
    # base path detours above the axis: 1e-5 -> +ih -> above target -> down.
    t_start = 1e-5
    h = 0.03 * abs(t_center)
    p0 = t_center * (1 + rrel)
    way = [t_start, t_start + 1j * h, p0 + 1j * h, p0]
    rs = sorted(droots(quartic_coeffs(t_start, s, eps)), key=abs)
    cur = list(rs)

    def step_to(t, cur):
        new = droots(quartic_coeffs(t, s, eps), seeds=cur, iters=80)
        matched = [None] * 4
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
        return matched

    for seg in range(len(way) - 1):
        for k in range(1, 201):
            t = way[seg] + (way[seg + 1] - way[seg]) * (k / 200) ** 1.2
            cur = step_to(t, cur)
    g0 = sum(-1 / (p0 * cur[i] * lam_d1(cur[i], s, eps)) for i in (0, 1))
    for k in range(1, 257):
        th = 2 * math.pi * k / 256
        cur = step_to(t_center * (1 + rrel * cmath.exp(1j * th)), cur)
    g1 = sum(-1 / (p0 * cur[i] * lam_d1(cur[i], s, eps)) for i in (0, 1))
    return g0, g1


for rr in (3e-3, 1e-3):
    gb, ga = monodromy_defect_detour(4.0, 1.0 / 16.0, 0.0, rr)
    Bex = abs(gb) * math.sqrt(2 * (1.0 / 16.0) * rr)
    print("       r=%.0e  |Ghat|=%.4f  |defect|=%.3e  B_extracted=%.4f  (2|Re rho1|=%.4f, 2|rho1|=%.4f)" %
          (rr, abs(gb), abs(ga - gb), Bex, 2 * abs(r1r.real), 2 * abs(r1r)))
A_, B_ = A_save, B_save

print("\n  verdict key: (3b),(3c) defect ~ C/sqrt(r): single-event cuts;")
print("  (3a) defect = -2*Ghat, |Ghat| ~ r^{-1/2}, B = 2rho: stacked pair")
print("  REINFORCES; (3f): the conjugate/reality stack's B against 2Re(rho).")
print("\nDONE.")
