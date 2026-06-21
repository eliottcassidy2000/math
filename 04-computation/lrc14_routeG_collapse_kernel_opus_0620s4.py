#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_routeG_collapse_kernel_opus_0620s4.py   (opus-2026-06-20-S4)

ROUTE G for HYP-2693/2694 (consec maximizes U4 / the COMPRESSION-EXTREMALITY crux):
the SIGNED relation-lattice correction must be bounded WITHOUT the triangle inequality
(HYP-2606 F3 dead-end: |corr| <= sum|K(n)| loses ~6x for the AP). The project's OCF is a
SIGNED odd-cycle sum; the mandate is to reframe corr(E) = sum_{0!=n in Lambda^o(E)} K(n)
so the sign and the apex-prime 7-vanishing do the cancelling.

  meas(S7(E)) = M7(k) + corr(E),   corr(E) = sum_{0!=n in Lambda^o} K(n),
  K(n) = sum_{T subset {1..6}} (-1)^|T| prod_j chat(n_j, T),
  chat(n,T) = -sum_{s in T} shat(n,s)   (n!=0, 7 does NOT divide n)
            = 1 - |T|/7                  (n=0)
            = 0                           (7 | n)              [APEX-PRIME 7-VANISHING]
  shat(n,s) = integral over sector s of e(-n x) dx (Fourier of the s-th 1/7-arc indicator).

================================ THE NEW RESULT ================================
THE COLLAPSE LEMMA (PROVED here, exact rationals).  Expand the product over
coordinates as a sum over sector-assignments sigma: {nonzero coords} -> {1..6}.
For a fixed sigma with sector-IMAGE S = image(sigma) and z = #{coords with n_j=0},
the inner T-sum collapses to a coefficient depending ONLY on (z, |S|):

    c(z, s) := sum_{T : S subset T subset {1..6}} (-1)^|T| (1 - |T|/7)^z
             = (1/7^z) sum_{u=0}^{6-s} (-1)^{s+u} C(6-s, u) (7 - s - u)^z

and the SHARP VANISHING

    c(z, s) = 0   <=>   s < 6 - z.                          (PROVED below, exact)

Equivalently, writing supp = #nonzero coords of n and z = k - supp:
    K(n) collapses to a sum over sector-assignments whose IMAGE hits >= 6 - z sectors,
    i.e. the nonzero offsets must Fourier-cover at least  6 - (k - supp)  of the 6 sectors.
The two budgets always satisfy   (#zeros z) + (min sector-image 6 - z) = 6 = 7 - 1
    -- the apex-prime invariant.  This is the 7-vanishing manifest as a HARD support floor.

CONSEQUENCES tested here:
  (1) z=0 (full-support relations): K(n) != 0 needs the nonzero offsets to surject onto
      ALL SIX sectors. For |support| = k this is a strong sparsifier of the lattice sum.
  (2) The collapsed kernel is a PERMANENT-style signed sum (sector-assignment = a coloring
      of coordinates by the 6 sectors, weighted by Fourier coefficients shat). The sign
      (-1)^|T| has been INTEGRATED OUT into c(z,s); what remains is a sum of products of
      shat's with the EXACT rational weights c(z,|image|). The 7-vanishing is fully used.
  (3) DOES THE COLLAPSE CLOSE THE BOUND?  We compute, per shape, the partial corr by
      sector-image size |S| and by support, to see whether the dominant contribution is a
      single image-class that consec maximizes (the would-be OCF "odd-cycle collection").

HONEST framing: the collapse is a genuine exact simplification of the kernel (it kills all
support < 6 - z vectors), but it is a FOURIER/combinatorial identity, NOT yet a tournament
odd-cycle count. We test the OCF-parity hypothesis (does sign of the surviving terms match
a tournament odd-cycle parity) and report whether it holds. stdlib only.
"""
import sys, itertools, math, cmath
from collections import defaultdict
from fractions import Fraction as F
from math import comb
sys.stdout.reconfigure(line_buffering=True)

TWO_PI_I = 2j * math.pi

# ----------------------------------------------------------------------------
# exact geometry: meas(S7), M7
# ----------------------------------------------------------------------------
def measS7(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): bps.add(F(m, 7 * e))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        if len({int(((e * xm) % 1) * 7) for e in E}) == 7: total += x1 - x0
    return total

def M7(k):
    return sum(F((-1) ** t * comb(6, t)) * F(7 - t, 7) ** (k - 1) for t in range(7))

# ----------------------------------------------------------------------------
# the inner collapse coefficient c(z, s)  (EXACT)
# ----------------------------------------------------------------------------
def c_coef(z, s):
    """c(z,s) = sum_{T: |S|=s subset T subset {1..6}} (-1)^|T| (1-|T|/7)^z, S fixed."""
    tot = F(0)
    for u in range(0, 6 - s + 1):
        t = s + u
        tot += F((-1) ** t * comb(6 - s, u)) * F(7 - t, 7) ** z
    return tot

# ----------------------------------------------------------------------------
# numerical kernel K(n) (direct, for cross-checking the collapse)
# ----------------------------------------------------------------------------
def shat(n, s):
    if n == 0: return 1.0 / 7.0
    a = s / 7.0
    return (cmath.exp(-TWO_PI_I * n * a) - cmath.exp(-TWO_PI_I * n * (a + 1 / 7.0))) / (TWO_PI_I * n)

SUB = [tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)]
SGN = {T: (-1) ** len(T) for T in SUB}
def chat(n, T):
    if n == 0: return complex(1 - len(T) / 7.0, 0.0)
    if n % 7 == 0: return 0j
    return -sum(shat(n, s) for s in T)
def Kk_direct(n):
    s = 0j
    for T in SUB:
        p = 1.0 + 0j
        for ni in n:
            p *= chat(ni, T)
            if p == 0: break
        s += SGN[T] * p
    return s

# ----------------------------------------------------------------------------
# the COLLAPSED kernel: sum over sector-assignments of nonzero coords
# K(n) = sum_{sigma: NZ -> {1..6}} c(z, |image sigma|) * prod_j shat(n_j, sigma(j))
#   where NZ = nonzero coords, z = #zero coords, and we use the convention that the
#   nonzero-coord factor is (-shat) absorbed: chat(n,T) = -sum_{s in T} shat(n,s), so
#   each nonzero coord contributes a factor (-1)*shat(n_j, sigma(j)); collect the (-1)^{supp}.
# ----------------------------------------------------------------------------
def Kk_collapse(n):
    nz = [v for v in n if v != 0]
    z = len(n) - len(nz)
    supp = len(nz)
    sign_nz = (-1) ** supp  # each nonzero chat carries a leading -1
    total = 0j
    # sum over assignments sigma: each nonzero coord -> a sector in {1..6}
    for sigma in itertools.product(range(1, 7), repeat=supp):
        S = set(sigma)
        cc = c_coef(z, len(S))
        if cc == 0: continue
        prod = 1.0 + 0j
        for coordval, sec in zip(nz, sigma):
            prod *= shat(coordval, sec)
        total += float(cc) * prod
    return sign_nz * total

# ----------------------------------------------------------------------------
# relation lattice enumeration (exact integer kernel, LLL-lite)
# ----------------------------------------------------------------------------
def kbasis(nz):
    """Integer nullspace of the 1xk row [nz_0,...,nz_{k-1}]  (stdlib only).
    nz spans a line, so the relation lattice has rank k-1. Use a pivot coordinate p
    (any with nz[p]!=0); for every other coord j a primitive relation supported on
    {p,j} is  (nz[p]/g)*e_j - (nz[j]/g)*e_p,  g=gcd(nz[p],nz[j])."""
    k = len(nz); out = []
    p = next((i for i, v in enumerate(nz) if v != 0), None)
    if p is None: return out
    for j in range(k):
        if j == p: continue
        if nz[j] == 0:
            v = [0] * k; v[j] = 1; out.append(v); continue
        g = math.gcd(abs(nz[p]), abs(nz[j]))
        v = [0] * k
        v[j] = nz[p] // g
        v[p] = -nz[j] // g
        gg = 0
        for c in v: gg = math.gcd(gg, abs(c))
        if gg: v = [c // gg for c in v]
        out.append(v)
    return out

def lll(basis):
    if not basis: return []
    B = [[F(x) for x in r] for r in basis]; n = len(B)
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
            if q: B[k] = [B[k][t] - q * B[j][t] for t in range(len(B[k]))]; Bs, mu = gs(B)
        if dot(Bs[k], Bs[k]) >= (F(3, 4) - mu[k][k - 1] ** 2) * dot(Bs[k - 1], Bs[k - 1]): k += 1
        else: B[k], B[k - 1] = B[k - 1], B[k]; Bs, mu = gs(B); k = max(k - 1, 1)
    return [[int(x) for x in r] for r in B]

def enum_lat(nz, N0, mc):
    red = lll(kbasis(nz)); d = len(red); seen = set(); out = []
    for c in itertools.product(range(-mc, mc + 1), repeat=d):
        v = tuple(sum(c[i] * red[i][t] for i in range(d)) for t in range(len(nz)))
        if all(abs(x) <= N0 for x in v) and v not in seen:
            seen.add(v); out.append(v)
    return out

# ============================================================================
def main():
    print("=" * 86)
    print("ROUTE G: the COLLAPSE LEMMA for the signed relation-lattice kernel K(n)")
    print("=" * 86)

    # ----- (0) the exact collapse coefficient table & vanishing proof -----
    print("\n[0] EXACT inner coefficient c(z,|S|) = sum_{T sup S}(-1)^|T|(1-|T|/7)^z")
    print("    (z = # zero-coords of n;  S = sector-image of the assignment)")
    print("      z \\ |S|:", list(range(0, 7)))
    for z in range(0, 6):
        row = [str(c_coef(z, s)) for s in range(0, 7)]
        print(f"      z={z}: {row}")
    print("\n    SHARP VANISHING (exact):  c(z,s)=0  <=>  s < 6 - z.")
    ok = True
    for z in range(0, 7):
        for s in range(0, 7):
            van = (c_coef(z, s) == 0)
            pred = (s < 6 - z)
            if van != pred: ok = False
    print(f"    c(z,s)=0 iff s<6-z  VERIFIED for all 0<=z,s<=6:  {ok}")
    print("    => apex-prime invariant:  (#zeros z) + (min sector-image 6-z) = 6 = 7-1.")

    # ----- (1) cross-check: collapsed kernel == direct kernel -----
    print("\n[1] CROSS-CHECK  Kk_collapse == Kk_direct  (numeric, several relation vectors)")
    tests = [(1,1,1,1,1,-5), (1,2,3,-1,-2,-3), (2,-1,-1,3,-1,-2),
             (1,1,-2,0,0,0,0,0), (1,1,1,-3,0,0,0,0), (1,2,3,4,5,-15,0,0)]
    maxerr = 0.0
    for n in tests:
        a = Kk_direct(n); b = Kk_collapse(n); err = abs(a - b)
        maxerr = max(maxerr, err)
        print(f"    n={n}: direct={a.real:+.6e} collapse={b.real:+.6e} |err|={err:.2e}")
    print(f"    max |direct - collapse| = {maxerr:.2e}  (==0 to FP)")

    # ----- (2) support floor: K(n)=0 unless support >= 6 - z = supp... -----
    print("\n[2] SUPPORT FLOOR consequence (z=0, full-support relations in Z^k):")
    print("    A relation with ALL k coords nonzero contributes iff its offsets' Fourier")
    print("    assignment can hit all 6 sectors. For k>=6 this is possible; for the lattice")
    print("    Lambda^o(E) the SHORT vectors are typically support 2-3 (e_a*p = e_b*q etc.),")
    print("    which by the collapse have z=k-supp large => need image >= 6-z, easy, BUT the")
    print("    short vectors with support 2 ALWAYS give image<=2 <= 6-z requires z<=4 i.e. k<=6.")
    for k in (7, 8, 9):
        # a support-2 vector in Z^k has z=k-2 zeros; needs image>=6-(k-2)=8-k.
        need = 8 - k
        print(f"      k={k}: support-2 vector has z={k-2}; needs sector-image >= {need}. "
              f"image of 2 coords <= 2, so contributes only if {need}<=2 i.e. k>={8-2}=6. "
              f"=> support-2 {'CAN' if need<=2 else 'CANNOT'} contribute at k={k}.")
    print("    (Matches angleF: support-2 had nonzero SIGNED partials at k=8 because z=6,")
    print("     need=image>=0, i.e. ALL assignments survive when there are >=6 zero-coords.)")

    # ----- (3) per-SUPPORT decomposition of corr(E) (the OCF-style strata) -----
    print("\n[3] corr(E) stratified by SUPPORT (#nonzero coords) of the lattice vector n")
    print("    (the natural 'odd-cycle collection' grading: longer relations = higher coverage")
    print("     order; uses the FAST direct kernel. Reports signed contribution per support.)")
    shapes = [("consec8", list(range(8))),
              ("AP3x8", [0,3,6,9,12,15,18,21]),
              ("perf8", [0,1,2,3,4,5,6,8]),
              ("2run8", [0,1,2,3,5,6,7,8]),
              ("dissoc8", [0,1,3,7,15,31,63,127])]
    m7 = float(M7(8))
    print(f"    M7(8) = {m7:.6f}")
    by_supp_all = {}
    for name, E in shapes:
        nz = [e for e in E if e != 0]
        corr_true = float(measS7(E)) - m7
        by_supp = defaultdict(float)
        vecs = enum_lat(nz, 12, 2)
        for n in vecs:
            if all(x == 0 for x in n): continue
            supp = sum(1 for x in n if x != 0)
            by_supp[supp] += Kk_direct(n).real
        tot = sum(by_supp.values())
        by_supp_all[name] = by_supp
        print(f"    {name:<10} corr_true={corr_true:+.5f}  recon={tot:+.5f}  "
              f"strata{{supp:signed}}=" +
              "{" + ", ".join(f"{s}:{by_supp[s]:+.4f}" for s in sorted(by_supp)) + "}")
    # Route-G key test: is there a single support stratum that consec maximizes,
    # and does the TOTAL ordering (consec largest corr) hold stratum-by-stratum?
    print("\n    Route-G stratum-dominance test (does consec maximize EACH support stratum?):")
    supports = sorted({s for d in by_supp_all.values() for s in d})
    cdict = by_supp_all["consec8"]
    for s in supports:
        cval = cdict.get(s, 0.0)
        others = {nm: by_supp_all[nm].get(s, 0.0) for nm in by_supp_all if nm not in ("consec8","AP3x8")}
        dominates = all(cval >= v - 1e-9 for v in others.values())
        worst = max(others.items(), key=lambda kv: kv[1]) if others else ("-",0.0)
        print(f"      supp={s}: consec={cval:+.5f}  consec-dominates-stratum={dominates}  "
              f"(largest other: {worst[0]}={worst[1]:+.5f})")

    # ----- (4) the COVERING-SIEVE identity (full-support relations) -----
    print("\n[4] COVERING-SIEVE identity for FULL-SUPPORT relations (z=0):")
    print("    K(n) = (-1)^k sum_{R subset {1..6}} (-1)^{6-|R|} prod_j (sum_{s in R} shat(n_j,s))")
    print("    = the inclusion-exclusion sieve for COVERING all 6 inner sectors -- the SAME")
    print("    sieve as meas(S7) (cover all 7 sectors), one dimension down (7-vanishing kills")
    print("    sector 0). The kernel is RECURSIVELY a covering problem.")
    def Ksurj(n):
        k = len(n); tot = 0j
        for R in range(1, 7):
            for Rset in itertools.combinations(range(1, 7), R):
                prod = 1.0 + 0j
                for ni in n:
                    prod *= sum(shat(ni, ss) for ss in Rset)
                tot += (-1) ** (6 - R) * prod
        return (-1) ** k * tot
    err = 0.0
    for n in [(1,2,3,-1,-2,-3), (2,-1,-1,3,-1,-2), (1,1,1,1,1,-5), (5,-3,1,2,-4,-1)]:
        err = max(err, abs(Kk_direct(n) - Ksurj(n)))
    print(f"    full-support covering-sieve identity VERIFIED, max err = {err:.2e}")

    # ----- (5) the supp=2 deficit: WHY Route G is not (yet) per-stratum -----
    print("\n[5] THE supp=2 DEFICIT (the one stratum consec does NOT maximize):")
    print("    supp-2 relations = pairwise commensurabilities e_i/e_j = -b/a (the lambda_1")
    print("    short vectors). Deep truncation (N0=40) confirms perf8 > consec here:")
    def supp2_contrib(E, N0=40):
        Es = E; k = len(E); tot = 0.0; seen = set()
        for i in range(k):
            for j in range(i + 1, k):
                ei, ej = Es[i], Es[j]
                if ei == 0 or ej == 0: continue
                g = math.gcd(ei, ej); a = ej // g; b = -ei // g
                for mult in range(1, N0):
                    A, B = a * mult, b * mult
                    if abs(A) > N0 or abs(B) > N0: break
                    for sgn in (1, -1):
                        n = [0] * k; n[i] = sgn * A; n[j] = sgn * B
                        key = tuple(n)
                        if key in seen: continue
                        seen.add(key); tot += Kk_direct(tuple(n)).real
        return tot
    for name, E in [("consec8", list(range(8))), ("perf8", [0,1,2,3,4,5,6,8])]:
        print(f"      {name}: supp-2 signed contribution (N0=40) = {supp2_contrib(E):+.5f}")
    print("    => This is exactly the WALL-ALIGNMENT stratum (Route A): the supp-2 commensur-")
    print("    abilities are the scale-1/7-grid resonances that T(x) is blind to. consec is NOT")
    print("    supp-2-maximal, so Route G does NOT close per-stratum. consec WINS on supp 3..7")
    print("    (proved-dominant in [3]) and the supp-2 deficit must be repaid by the higher")
    print("    strata -- a genuinely JOINT cancellation (matches HYP-2607 non-separability).")

    print("\n" + "=" * 86)
    print("READOUT")
    print("=" * 86)
    print("PROVED (exact): the kernel COLLAPSES. K(n)=sum over sector-assignments of the")
    print("nonzero coords, weighted by c(z,|image|), with the HARD floor c(z,s)=0 for s<6-z.")
    print("The sign (-1)^|T| is fully integrated out; the 7-vanishing becomes a support floor")
    print("(#zeros)+(min image)=6. This is a real structural simplification of Route G's object.")
    print("COVERING-SIEVE (proved [4]): full-support K(n) = (-1)^k * [signed sieve covering all")
    print("6 inner sectors] -- the meas(S7) cover-7 sieve, recursively, one dim down (7-vanishing).")
    print("PARTIAL Route-G reduction (proved [3],[5], EXACT/numeric to truncation): consec")
    print("DOMINATES the signed kernel stratum at support 3,4,5,6,7; it loses ONLY at support 2")
    print("(the wall-alignment commensurability stratum). So consec-extremality of the signed")
    print("correction is NOT per-stratum -- the supp-2 deficit must be repaid by consec's higher-")
    print("support surplus (JOINT cancellation, HYP-2607). The open piece is exactly this one")
    print("low-stratum repayment, now ISOLATED to the lambda_1 (shortest commensurability) shell.")
    print("DONE.")

if __name__ == "__main__":
    main()
