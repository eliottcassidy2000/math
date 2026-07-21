#!/usr/bin/env python3
"""
levelset_sumrule_L1tube_boxeph_S187.py  (HYP-8585, THM-1785 sections 1-2)

TARGET 1 of S186r: the value-hijacked (L1) sub-case, attacked with the
referee's global residue identity.

(S1) PS-RESIDUE LEMMA (verify): at a simple fold (Λ'(u_c)=0, Λ''(u_c)≠0,
     v = Λ(u_c)), the on-arc pair-sum collision limit equals
       PS·t|_{t=1/v} = Res_{u_c}[ v / (u (v − Λ(u))) ]
                     = 2v/(u_c²Λ'') + (2/3)·v·Λ'''/(u_cΛ''²).
     Check against direct root-residue sums approaching the collision.

(S2) LEVEL-SET SUM RULE (verify): summing the residue theorem for
     v/(u(v−Λ)) over ALL its poles (no pole at u = 0/∞ for mixed Λ):
       Σ_{folds with value v} PS_i·t = − Σ_{simple roots u* of Λ=v} v/(u*·(−Λ'(u*)))
     i.e. Σ PS_i·t = Σ_{u*} v/(u*Λ'(u*)) with the orientation fixed by the
     residue calculus (sign checked numerically here).

(S3) THE L1-TUBE MEASUREMENT on the hijack witness P₄ = ZW + Z⁹W⁷ + W
     (Λ = 2s + w¹⁶u² + w/u, w = √(2s)): march t = (1+iδ)/v(s₀) with
     s₀ → 0 along the hijacked arc and measure
       I(t) = ∫ e^{−s} |Ĝ(s,t)| ds
     over the spike window and off-spike region separately. HEIGHT ~ s₀^{−2}
     is the S186r refutation; the question is height × width: conjecture
     I(t) = O(1) (bounded) or O(log 1/s₀) — measured, not assumed.
     Ĝ computed by tracked continuation (0-group), pure python.

boxeph-2026-07-20-S187.
"""
import cmath
import math


# ---------- generic machinery (from S183, trimmed) ----------

def droots(coeffs, seeds=None, iters=300):
    c0 = coeffs[0]
    cs = [c / c0 for c in coeffs]
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


# ---------- (S1) + (S2): the pair model with charge-0, exact checks ----------
# Lambda = q0(s) + w(p1 u + pm1/u): folds at u_c = ±sqrt(pm1/p1),
# plus a value-level set with two more intersections when |level| large.
print("=" * 78)
print("(S1) PS-residue lemma at a simple fold (pair model + cubic model):")
q0c, p1, pm1 = 0.4 + 0.1j, 1.0 + 0.3j, 0.8 - 0.2j
s = 2.0
w = math.sqrt(2 * s)
lam = lambda u: q0c * s + w * (p1 * u + pm1 / u)
lam1 = lambda u: w * (p1 - pm1 / (u * u))
lam2 = lambda u: 2 * w * pm1 / u ** 3
lam3 = lambda u: -6 * w * pm1 / u ** 4
uc = cmath.sqrt(pm1 / p1)
v = lam(uc)
PSpred = 2 * v / (uc * uc * lam2(uc)) + (2.0 / 3) * v * lam3(uc) / (uc * lam2(uc) ** 2)
# direct: pair of roots of 1 - t Lambda = 0 near u_c at t = (1-eps)/v; sum residues of 1/(u(1-tLam))
for eps in (1e-3, 1e-5, 1e-7):
    t = (1 - eps) / v
    # roots: w p1 u^2 + (q0 s - 1/t) u + w pm1 = 0
    A2, B2, C2 = w * p1, q0c * s - 1 / t, w * pm1
    disc = cmath.sqrt(B2 * B2 - 4 * A2 * C2)
    r1 = (-B2 + disc) / (2 * A2)
    r2 = (-B2 - disc) / (2 * A2)
    # residue of 1/(u(1-tLam)) at root r: 1/(r * (-t) * lam1(r))
    S = sum(1 / (r * (-t) * lam1(r)) for r in (r1, r2))
    # compare with PSpred (multiplied by t at collision => PS*t; our S is the pair-sum of residues)
    print("  eps=%.0e  pair-residue-sum = %.6f%+.6fj   Res-formula = %.6f%+.6fj" %
          (eps, S.real, S.imag, PSpred.real, PSpred.imag))

print("\n(S2) level-set sum rule on a 2-fold model (cubic charged part):")
# Lambda = q0 s + w(p1 u + p2 u^2 + pm1/u): critical points: 3 of them
p2c = 0.5 + 0.25j
lamB = lambda u: q0c * s + w * (p1 * u + p2c * u * u + pm1 / u)
lamB1 = lambda u: w * (p1 + 2 * p2c * u - pm1 / (u * u))
lamB2 = lambda u: w * (2 * p2c + 2 * pm1 / u ** 3)
lamB3 = lambda u: -6 * w * pm1 / u ** 4
# critical points: roots of u^2 lamB1/w = p1 u^2 + 2 p2c u^3 - pm1 = 0
crs = droots([2 * p2c, p1, 0, -pm1])
uc0 = crs[0]
v0 = lamB(uc0)
PS0 = 2 * v0 / (uc0 * uc0 * lamB2(uc0)) + (2.0 / 3) * v0 * lamB3(uc0) / (uc0 * lamB2(uc0) ** 2)
# level set Lambda = v0: w p2 u^3 + w p1 u^2 + (q0 s - v0) u + w pm1 = 0: 3 roots,
# one double at uc0 counts as the fold; the remaining simple root u*:
lroots = droots([w * p2c, w * p1, q0c * s - v0, w * pm1])
# identify the simple root (farthest from uc0)
lroots = sorted(lroots, key=lambda r: -abs(r - uc0))
ustar = lroots[0]
RHS = v0 / (ustar * lamB1(ustar))
print("  fold u_c = %.4f%+.4fj ; simple level-set root u* = %.4f%+.4fj" %
      (uc0.real, uc0.imag, ustar.real, ustar.imag))
print("  PS(fold) = %.6f%+.6fj ;  v/(u* Lam'(u*)) = %.6f%+.6fj ;  SUM = %.2e" %
      (PS0.real, PS0.imag, RHS.real, RHS.imag, abs(PS0 + RHS)))
print("  (sum rule: PS + spectator term = 0 for this 1-fold-1-spectator level)")

# ---------- (S3): L1 tube integral on the hijack witness ----------
print("\n(S3) L1 tube measurement on P4 = ZW + Z^9W^7 + W  (hijacked arc):")
print("     Lambda = 2s + w^16 u^2 + w/u ; arc t(s) = 1/v(s), v ~ 2s + edge")


def lamP4(u, s):
    w = math.sqrt(2 * s)
    return 2 * s + (w ** 16) * u * u + w / u


def lamP4_1(u, s):
    w = math.sqrt(2 * s)
    return 2 * (w ** 16) * u - w / (u * u)


def quartic_coeffs_P4(t, s):
    # u * (1 - t Lambda) = -t w^16 u^3 + (1 - 2ts) u - t w  (cubic in u)
    w = math.sqrt(2 * s)
    return [-t * w ** 16, 0.0, 1 - 2 * t * s, -t * w]


def Ghat_P4(s, t, nsteps=200):
    # tracked 0-group: K- = 1 (one root near 0 at small t: from w/u term)
    t0 = 1e-6 * abs(t)
    rs = sorted(droots(quartic_coeffs_P4(t0, s)), key=abs)
    cur = list(rs)
    for k in range(1, nsteps + 1):
        tt = t0 + (t - t0) * (k / nsteps) ** 1.3
        new = droots(quartic_coeffs_P4(tt, s), seeds=cur, iters=60)
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
    u0 = cur[0]
    return -1 / (t * u0 * lamP4_1(u0, s))


def v_of(s):
    # fold of P4 nearest the hijacked arc: critical points: 2 w^16 u^3 = w
    w = math.sqrt(2 * s)
    uc = (w / (2 * w ** 16)) ** (1.0 / 3)
    return lamP4(uc, s)


def simpson_abs(f, a, b, n):
    h = (b - a) / n
    tot = abs(f(a)) + abs(f(b))
    for k in range(1, n):
        tot += abs(f(a + k * h)) * (4 if k % 2 == 1 else 2)
    return tot * h / 3


for s0 in (0.3, 0.1, 0.03, 0.01):
    v0 = v_of(s0)
    t = (1 + 1j * 3e-3) / v0
    # spike window and off-spike: integrate e^{-s}|Ghat| over s in [s0/3, 3 s0]
    # (window) and [3 s0, 8] (tail toward normal region)
    f = lambda s: math.exp(-s) * Ghat_P4(s, t)
    Iwin = simpson_abs(f, s0 / 3, 3 * s0, 60)
    Itail = simpson_abs(f, 3 * s0, 8.0, 80)
    gpk = abs(Ghat_P4(s0, t))
    print("  s0=%.3f  |t|=%.3e  peak|Ghat|=%.3e  I_window=%.4f  I_tail=%.4f  I_tot=%.4f" %
          (s0, abs(t), gpk, Iwin, Itail, Iwin + Itail))

print("\nDONE.")
