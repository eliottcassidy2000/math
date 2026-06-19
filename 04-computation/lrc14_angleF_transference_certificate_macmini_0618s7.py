#!/usr/bin/env python3
"""
lrc14_angleF_transference_certificate_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

ANGLE F part 4 — the decisive test: can a covolume / transference bound on the
relation-lattice 1/height sum certify  meas(S7(E)) <= cap_k  with the THM-532 margin?

CHAIN (all rigorous except the final numeric constant):
  corr(E) = sum_{0!=n in Lambda^o(E)} K(n),
  |K(n)| <= sum_{T<={1..6}} prod_{j: n_j!=0} |chat_T(n_j)|,  |chat_T(n)| <= |T| s7/(pi|n|),
           s7 = sin(pi/7) = 0.43388.   (per-sector Fourier decay, EXACT constant)
So  |corr(E)| <= H(E) := sum_{0!=n in Lambda^o} B(n),  B(n)=sum_T prod_{n_j!=0}|T|s7/(pi|n_j|).

We compute H(E) (truncated, the tail is geometric & controlled) and ask:
  is  M7(k) + H(E) <= cap_k  i.e.  H(E) <= margin_k ?
If YES for the AP (the worst shape), the certificate CLOSES with the absolute
constant.  We report H(E) for the dangerous shapes and the AP, and the smallest
covolume vs the bound, to see exactly where (if anywhere) the bound holds.

We ALSO compute the SUPPORT-3-ONLY bound (the canon's W(E) regime) to compare
against the honest gap reported in THM-532.
"""
from __future__ import annotations
import sys, itertools, math
from fractions import Fraction as F
import sympy
sys.stdout.reconfigure(line_buffering=True)

S7 = math.sin(math.pi / 7)          # 0.433884
KC = S7 / math.pi                   # per-|T| per-1/|n| Fourier constant = 0.138108


def measS7(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        secs = {int(((e * xm) % 1) * 7) for e in E}
        if len(secs) == 7:
            total += x1 - x0
    return total


def M7(k):
    return sum(F((-1) ** t * math.comb(6, t)) * F(7 - t, 7) ** (k - 1) for t in range(7))


# |K(n)| bound: sum over T<={1..6} of prod_{j: n_j!=0} |T| s7/(pi |n_j|)
#  = sum_{t=1}^{6} C(6,t) (t*KC)^{supp} / prod|n_j|   where supp = #{j: n_j!=0}
def Kbound(n):
    supp = sum(1 for x in n if x != 0)
    if supp == 0:
        return 0.0
    invprod = 1.0
    for x in n:
        if x != 0:
            invprod /= abs(x)
    s = 0.0
    for t in range(1, 7):
        s += math.comb(6, t) * (t * KC) ** supp
    return s * invprod


def lll(basis):
    if not basis:
        return []
    B = [[F(x) for x in r] for r in basis]
    n = len(B)
    def dot(u, v): return sum(a * b for a, b in zip(u, v))
    def gs(B):
        Bs, mu = [], [[F(0)] * n for _ in range(n)]
        for i in range(n):
            bi = list(B[i])
            for j in range(i):
                mu[i][j] = dot(B[i], Bs[j]) / dot(Bs[j], Bs[j])
                bi = [bi[t] - mu[i][j] * Bs[j][t] for t in range(len(bi))]
            Bs.append(bi)
        return Bs, mu
    Bs, mu = gs(B); k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [B[k][t] - q * B[j][t] for t in range(len(B[k]))]; Bs, mu = gs(B)
        if dot(Bs[k], Bs[k]) >= (F(3, 4) - mu[k][k - 1] ** 2) * dot(Bs[k - 1], Bs[k - 1]):
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]; Bs, mu = gs(B); k = max(k - 1, 1)
    return [[int(x) for x in r] for r in B]


def kernel_basis(nz):
    M = sympy.Matrix([nz]); out = []
    for v in M.nullspace():
        L = 1
        for x in v:
            fr = F(int(x.p), int(x.q)); L = L * fr.denominator // math.gcd(L, fr.denominator)
        iv = [int(x * L) for x in v]
        g = 0
        for c in iv:
            g = math.gcd(g, abs(c))
        if g:
            iv = [c // g for c in iv]
        out.append(iv)
    return out


def enum_lattice(nz, N0, mc):
    red = lll(kernel_basis(nz)); d = len(red); seen = set(); out = []
    for c in itertools.product(range(-mc, mc + 1), repeat=d):
        v = tuple(sum(c[i] * red[i][t] for i in range(d)) for t in range(len(nz)))
        if all(abs(x) <= N0 for x in v) and v not in seen:
            seen.add(v); out.append(v)
    return out


def Hbound(E, N0=30, mc=6, support_cap=None):
    """Truncated sum of |K(n)| upper bound over relation lattice. support_cap limits
    the support size (3 = canon W regime)."""
    nz = [e for e in E if e != 0]
    tot = 0.0
    for n in enum_lattice(nz, N0, mc):
        if all(x == 0 for x in n):
            continue
        if support_cap is not None and sum(1 for x in n if x != 0) > support_cap:
            continue
        tot += Kbound(n)
    return tot


def covol_sat(nz):
    g = 0
    for x in nz:
        g = math.gcd(g, x)
    return math.sqrt(sum((x // g) ** 2 for x in nz))


def main():
    print("=" * 94)
    print("ANGLE F part 4: does the covolume/transference 1/height bound certify S7<=cap?")
    print(f"per-sector Fourier const s7/pi = {KC:.6f}; absolute |K(n)| bound = sum_t C(6,t)(t*{KC:.4f})^supp/height")
    print("=" * 94)
    caps = {8: F(2243, 5880)}
    k = 8
    cap = float(caps[k]); m7 = float(M7(k)); margin = cap - m7
    print(f"k={k}: M7={m7:.5f}  cap={cap:.5f}  margin={margin:.5f}")
    print("-" * 94)
    print(f"{'shape':<30}{'corr_true':>11}{'covol':>9}{'Hbound':>10}{'H_supp3':>10}{'M7+H<=cap?':>12}")
    shapes = [
        ("AP {0..7}", list(range(8))),
        ("3AP {..6,8}", [0, 1, 2, 3, 4, 5, 6, 8]),
        ("2run40", [0, 1, 2, 3, 40, 41, 42, 43]),
        ("perf", [0, 2, 3, 4, 5, 6, 7, 9]),
        ("Sidon", [0, 1, 3, 7, 12, 20, 30, 44]),
        ("dissoc", [0, 1, 3, 7, 15, 31, 63, 127]),
    ]
    for name, E in shapes:
        nz = [e for e in E if e != 0]
        corr = float(measS7(E)) - m7
        cov = covol_sat(nz)
        H = Hbound(E, N0=24, mc=5)
        H3 = Hbound(E, N0=24, mc=5, support_cap=3)
        ok = "YES" if m7 + H <= cap else "no"
        print(f"{name:<30}{corr:>11.5f}{cov:>9.2f}{H:>10.4f}{H3:>10.4f}{ok:>12}")
    print("-" * 94)
    print("READOUT:")
    print(" * corr_true is the SIGNED truth; Hbound is the ABSOLUTE (triangle-ineq) bound.")
    print(" * If Hbound >> margin even when corr_true << margin, the absolute bound is too")
    print("   lossy (sign cancellation is essential) -> a covolume bound on |corr| canNOT")
    print("   close on its own; one needs the SIGNED structure (THM-503 vanishing alignment).")
    print(" * H_supp3 isolates the canon W(E) regime; compare to the THM-532 honest gap.")
    print()
    print("DONE.")


if __name__ == "__main__":
    main()
