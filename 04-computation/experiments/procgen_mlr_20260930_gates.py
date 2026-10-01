#!/usr/bin/env python3
"""procgen_mlr_20260930_gates.py -- Part C (Collatz side) of the multiplicative lonely runner lane.

The cycle-gate exponential sum S(h) = sum_w e(h c_w/M), M = |2^p - 3^a|, c_w = sum_i 3^(a-1-i) 2^(s_i),
is a product over the runners h 3^m 2^s / M along the staircase path of the word.  Here:

C1  full spectra (FFT, M <= 2^20, 14 <= p <= 20): the large values vs the extended runner numerator
    u_ext(h) = min{|r| : r = 2^s 3^m h mod M, 0 <= s < p, -a <= m < a}; the dichotomy function
    U_eps = max{u_ext(h) : |S(h)| >= eps C}; the gates lane's unidentified peaks.
C2  crowding along the path vs |S|: the switching bound B(h) (gates lane, PROVED |S|/C <= B) and the
    path-weighted near fraction.
C3  major arcs vs minor arcs: structured contribution E_U, minor-arc L1 mass, the Parseval floor;
    full-spectrum clocks exactly, sparse clocks 24 <= p <= 40 via the DP on the structured set + samples.
Prints to stdout; raises on failed checks.
"""
import math
import random
import sys
import time

import numpy as np

import procgen_mlr_20260930_core as core
from procgen_mlr_20260930_core import check

T0 = time.time()


def hdr(s):
    print("\n" + "=" * 100 + "\n" + s + "\n" + "=" * 100)
    sys.stdout.flush()


def u_scan(p, a, M, chunk=1 << 19):
    """u_box (m >= 0) and u_ext (-a <= m < a) for all h mod M (chunked over h to cap memory)."""
    ubox = np.full(M, M, dtype=np.int64)
    uext = np.full(M, M, dtype=np.int64)
    G = [(core.mult(s, m, M), m >= 0) for s in range(p) for m in range(-a, a)]
    for h0 in range(0, M, chunk):
        h = np.arange(h0, min(M, h0 + chunk), dtype=np.int64)
        be = uext[h0:h0 + len(h)]
        bb = ubox[h0:h0 + len(h)]
        for g, inbox in G:
            r = (g * h) % M
            c = np.minimum(r, M - r)
            np.minimum(be, c, out=be)
            if inbox:
                np.minimum(bb, c, out=bb)
    ubox[0] = uext[0] = 0
    return ubox, uext


def full_clocks():
    gc = core.gates_core()
    out = []
    for p in range(14, 21):
        for a in range(2, p):
            M = abs(gc.gate(p, a))
            C = math.comb(p, a)
            if 2000 <= M <= (1 << 20) and C >= 1000:
                out.append((p, a))
    return out


def c1_c3_full():
    hdr("C1/C3a. full spectra (certified DP), 14 <= p <= 20, M <= 2^20, C >= 1000: structure of the large values")
    import procgen_gates_20260925_stats as gs
    rows = []
    unexpl = []
    for (p, a) in full_clocks():
        M, S, nres = core.gate_spectrum_dp(p, a)
        C = math.comb(p, a)
        N = int(round(nres[0]))
        absS = np.abs(S) / C
        ubox, uext = u_scan(p, a, M)
        hs = np.arange(1, M)
        rnd = math.sqrt(math.log(M) / C)
        top = hs[np.argsort(-absS[hs])[:20]]
        sixfree_small = sum(1 for h in top if uext[h] <= 50)
        top1 = int(top[0])
        # gates-lane decomposition of the top 4 (u minimal over h = u 2^-j 3^k, j<=p, k<=a)
        for h in top[:4]:
            uG = gs.decompose(int(h), p, a, M, 3)
            if abs(uG[0]) > 3 and uext[h] <= 3:
                unexpl.append((p, a, int(h), uG, int(uext[h])))
        # dichotomy function
        Ueps = {}
        for eps in (0.10, 0.05, 0.03):
            sel = absS[hs] >= eps
            if eps < 4 * rnd:
                Ueps[eps] = None          # eps not above 4x the random level sqrt(ln M / C): not meaningful
            else:
                Ueps[eps] = int(uext[hs][sel].max()) if sel.any() else 0
        # 6-free check for u_ext of all frequencies with |S| >= 0.05 C and u_ext < M^(1/2)
        sel = (absS[hs] >= 0.05) & (uext[hs] < math.isqrt(M))
        nsix = int(((uext[hs][sel] % 2 != 0) & (uext[hs][sel] % 3 != 0)).sum())
        # major / minor arcs
        Coll = float((nres ** 2).sum())
        E2 = M * Coll - C * C  # sum_{h != 0} |S|^2 (exact identity, float)
        arcs = {}
        for U in (1, 7, 50):
            struct = uext <= U          # includes h = 0
            EU = float(np.real(S[struct].sum())) / M
            off = ~struct
            L1 = float(np.abs(S[off]).sum()) / M
            supoff = float(np.abs(S[off]).max()) if off.any() else 0.0
            e2off = float((np.abs(S[off]) ** 2).sum())
            share = 1 - e2off / E2 if E2 > 0 else 0
            nsig = int(struct.sum())
            # PROVED lower bound (Prop. B): (1/M) sum_off |S| >= (sum_off |S|^2) / (M sup_off), with
            # sum_off |S|^2 = M Coll - sum_{Sigma} |S|^2 (exact identity; Coll = sum_r n_r^2 >= C)
            lb = e2off / (M * supoff) if supoff > 0 else 0.0
            check(L1 >= lb * (1 - 1e-9), f"Prop B lower bound violated at {(p, a, U)}")
            arcs[U] = (nsig - 1, EU, N - EU, L1, supoff / math.sqrt(C), share, nsig / M, lb)
        rows.append(dict(p=p, a=a, M=M, C=C, N=N, rnd=rnd, top1=top1, top1S=float(absS[top1]),
                         top1u=int(uext[top1]), sixfree_small=sixfree_small, Ueps=Ueps, nsix=nsix,
                         nsel=int(sel.sum()), arcs=arcs, CM=C / M,
                         topu=[int(uext[h]) for h in top[:8]]))
    # print per-clock summary (subset) and aggregates
    print("  per clock: top-1 |S|/C and u_ext; u_ext of top 8; #top-20 with u_ext <= 50; U_eps = max u_ext over |S| >= eps C")
    print("  (U_eps shown only when eps >= 4 sqrt(ln M / C), i.e. eps above the random level; '-' otherwise)")
    print("   p  a        M        C   N   top1|S|/C u_ext   u_ext(top8)                       #<=50   U_.10  U_.05  U_.03")
    for r in rows:
        if r["a"] % 3 == 0 or r["p"] == 20:
            U = r["Ueps"]
            fm = lambda v: '-' if v is None else str(v)
            print(f"  {r['p']:2d} {r['a']:2d} {r['M']:8d} {r['C']:8d} {r['N']:2d}   {r['top1S']:.3f}  {r['top1u']:5d}   "
                  f"{str(r['topu']):34s} {r['sixfree_small']:3d}   {fm(U[0.10]):>5s}  {fm(U[0.05]):>5s}  {fm(U[0.03]):>6s}")
    n = len(rows)
    n1 = sum(1 for r in rows if r["top1u"] <= 1)
    n7 = sum(1 for r in rows if r["top1u"] <= 7)
    tot_sel = sum(r["nsel"] for r in rows)
    tot_six = sum(r["nsix"] for r in rows)
    print(f"\n  {n} clocks.  top-1 frequency has u_ext = 1 on {n1}/{n}, u_ext <= 7 on {n7}/{n}.")
    print(f"  number of the top-20 frequencies with u_ext <= 50: min over the {n} clocks = {min(r['sixfree_small'] for r in rows)}")
    allU = [(r['Ueps'][0.10], r['Ueps'][0.05]) for r in rows]
    u10 = [u for u, _ in allU if u is not None]
    u05 = [u for _, u in allU if u is not None]
    print(f"  U_0.10 over the {len(u10)} clocks where 0.10 >= 4x random level: max {max(u10)}; "
          f"U_0.05 over the {len(u05)} clocks where meaningful: max {max(u05) if u05 else None}")
    print(f"  among all h with |S| >= 0.05 C and u_ext < sqrt(M): {tot_six}/{tot_sel} have u_ext coprime to 6 "
          f"(Triangle Lemma (c): apex numerators are 6-free)")
    for p in range(14, 21):
        rp = [r for r in rows if r["p"] == p]
        if not rp:
            continue
        def med(key):
            v = [r["Ueps"][key] for r in rp if r["Ueps"][key] is not None]
            return (int(np.median(v)), max(v), len(v)) if v else None
        print(f"   p={p}: clocks {len(rp):2d};  U_0.10 (median,max,#) = {med(0.10)};  U_0.05 = {med(0.05)};  U_0.03 = {med(0.03)}")
    if unexpl:
        print("\n  gates-lane top-4 frequencies whose minimal |u| over h = u 2^-j 3^k (j<=p,k<=a) exceeds 3, but u_ext <= 3:")
        for (p, a, h, uG, ue) in unexpl[:12]:
            print(f"    ({p},{a}) h={h}: gates decomposition (u,j,k)={uG}, u_ext={ue}")
        print(f"    ({len(unexpl)} such cases)")
    # C3a: arcs
    hdr("C3a. major vs minor arcs on the full-spectrum clocks (Sigma_U = {h : u_ext(h) <= U}, h = 0 included)")
    print("  E_U = (1/M) sum_{Sigma_U} S(h);  N = E_U + (1/M) sum_{off} S(h);  L1_off = (1/M) sum_{off} |S(h)|;")
    print("  Prop. B: L1_off >= LB = (sum_off |S|^2)/(M sup_off |S|), sum_off |S|^2 = M Coll - sum_Sigma |S|^2 (checked on every row)")
    print("   p  a        M        C    C/M   N  | U: |Sig|/M  Esh(Sig)   E_U       N-E_U    L1_off   LB(PropB)  sup_off/sqrtC")
    for r in rows:
        if r["a"] % 4 == 0 or r["CM"] > 1:
            for U in (1, 7, 50):
                s_, EU, res, L1, sup, share, frac, lb = r["arcs"][U]
                lead = f"  {r['p']:2d} {r['a']:2d} {r['M']:8d} {r['C']:8d} {r['CM']:6.3f} {r['N']:2d}" if U == 1 else " " * 39
                print(f"{lead}  | {U:2d}: {frac:6.3f}   {share:6.3f}  {EU:9.4f}  {res:9.4f}   {L1:8.2f}   {lb:8.2f}     {sup:6.2f}")
    # whenever Sigma_U carries at most 90% of the energy, the triangle-inequality route fails (|E_U| + L1_off >= 1)
    nreg = 0
    for r in rows:
        for U in (1, 7, 50):
            s_, EU, res, L1, sup, share, frac, lb = r["arcs"][U]
            if share <= 0.9:
                nreg += 1
                check(abs(EU) + L1 >= 1 and lb >= 1, f"triangle inequality certifies at {(r['p'], r['a'], U)}")
    med_ratio = np.median([r["arcs"][7][3] / math.sqrt(r["C"]) for r in rows])
    inc = [r["arcs"][U] for r in rows for U in (1, 7, 50) if r["arcs"][U][5] <= 0.9]
    print(f"\n  over the included cases: LB in [{min(x[7] for x in inc):.2f}, {max(x[7] for x in inc):.2f}], "
          f"L1_off in [{min(x[3] for x in inc):.1f}, {max(x[3] for x in inc):.1f}] (units of M)")
    print(f"  in all {nreg} (clock, U) cases where Sigma_U carries <= 90% of the energy: LB >= 1, so |E_U| + L1_off >= 1")
    print("  (no certificate of N = 0 by the triangle inequality, whatever the bound on |S| off Sigma_U);")
    print(f"  median L1_off/sqrt(C) at U = 7: {med_ratio:.3f}  (minor-arc L1 mass ~ sqrt(C) * M)")
    return rows


def c1b_larger():
    hdr("C1b. larger clocks (blocked DP for |S(h)|, all h): does the S-dichotomy U_eps stay bounded as p grows?")
    gc = core.gates_core()
    print("   p  a        M         C   top1 u_ext  #top20 u_ext<=50   U_.10  U_.05  U_.03  (eps shown only if >= 4 sqrt(ln M/C))")
    for (p, a) in ((21, 13), (21, 11), (22, 13), (23, 14), (24, 15)):
        t = time.time()
        M = abs(gc.gate(p, a))
        C = math.comb(p, a)
        absS = gc.S_abs_all_dp(p, a) / C          # h = 0 .. M//2
        H = len(absS)
        ubox, uext = u_scan(p, a, M)
        uh = uext[:H]
        hs = np.arange(1, H)
        rnd = math.sqrt(math.log(M) / C)
        top = hs[np.argsort(-absS[hs])[:20]]
        Ue = {}
        for eps in (0.10, 0.05, 0.03):
            sel = absS[hs] >= eps
            Ue[eps] = '-' if eps < 4 * rnd else (str(int(uh[hs][sel].max())) if sel.any() else '0')
        print(f"  {p:2d} {a:2d} {M:9d} {C:9d}   {int(uh[top[0]]):5d}      {sum(1 for h in top if uh[h] <= 50):3d}          "
              f"{Ue[0.10]:>5s}  {Ue[0.05]:>5s}  {Ue[0.03]:>5s}   [{time.time() - t:.0f}s]")
        sys.stdout.flush()
        del absS, ubox, uext, uh


def c1c_named():
    hdr("C1c. the gates lane's named frequencies, re-read through u_ext; a two-piece chain")
    gc = core.gates_core()

    def uext_one(h, p, a, M):
        best = None
        for s_ in range(p):
            for m in range(-a, a):
                r = core.centred(h * core.mult(s_, m, M), M)
                if best is None or abs(r) < abs(best[0]):
                    best = (r, s_, m)
        return best
    named = [((24, 15), 729492, "u=-1,j=24,k=10"), ((24, 15), 786952, "u=+18431,j=4,k=10"),
             ((24, 15), 486328, "u=-1,j=23,k=9"), ((24, 15), 1094238, "u=-18431,j=3,k=9"),
             ((27, 17), 1294383, "u=-37,j=27,k=5")]
    for (p, a), h, lab in named:
        M = abs(gc.gate(p, a))
        r, s_, m = uext_one(h, p, a, M)
        print(f"  ({p},{a}) h={h:8d} [gates: {lab:18s}]  u_ext = {abs(r):5d}:  h * 2^{s_} * 3^{m} = {r:+d} (mod {M})")
    for (p, a), u, j, k in (((25, 15), -95, 12, 1), ((29, 18), -175, 29, 7)):
        M = abs(gc.gate(p, a))
        h = u * pow(pow(2, j, M), -1, M) * pow(3, k, M) % M
        r, s_, m = uext_one(h, p, a, M)
        print(f"  ({p},{a}) h={h:9d} [gates: u={u},j={j},k={k}]  u_ext = {abs(r)}:  h * 2^{s_} * 3^{m} = {r:+d} (mod {M})")
        check(abs(r) == abs(u), "named partially structured frequency has a different u_ext")
    # two-piece chain at (20,8)
    p, a, h = 20, 8, 141509
    M = abs(gc.gate(p, a))
    R = M // 8
    res, near, comps = core.near_components(M, h, R, range(0, p), range(0, a))
    W = core.path_weights(p, a)
    pieces = []
    for comp in comps:
        j0 = min(P[0] for P in comp)
        k0 = min(P[1] for P in comp)
        wt = sum(W[P] for P in comp) / a
        if wt > 0.1:
            pieces.append(((j0, k0), res[(j0, k0)], len(comp), round(float(wt), 3)))
    r, s_, m = uext_one(h, p, a, M)
    print(f"  two-piece chain at (20,8), h = {h}, R = M/8: path-carrying near components (apex, numerator, size, path weight):")
    print(f"    {pieces};  u_ext = {abs(r)} (h * 2^{s_} * 3^{m} = {r:+d}): one orbit point seen on both sides of 2^20 = 3^8 (mod M)")


def c2_switching():
    hdr("C2. crowding along the path: |S|/C <= B(h) (switching bound, gates Prop. SW) and the path-near fraction")
    import procgen_gates_20260925_switching as sw
    rng = random.Random(11)
    for (p, a) in ((20, 8), (20, 12), (18, 11)):
        M, S, nres = core.gate_spectrum_dp(p, a)
        C = math.comb(p, a)
        absS = np.abs(S) / C
        W = core.path_weights(p, a)
        hs = np.arange(1, M)
        top = [int(h) for h in hs[np.argsort(-absS[hs])[:10]]]
        rnd = [rng.randrange(1, M) for _ in range(60)]
        B = sw.switching_B(p, a, 3, M, top + rnd)
        R = M // 8

        def pnear(h):
            tot = 0.0
            for s in range(p):
                for m in range(a):
                    if W[s, m] > 0:
                        r = core.centred(h * core.mult(s, m, M), M)
                        if abs(r) < R:
                            tot += W[s, m]
            return tot / a
        pt = [pnear(h) for h in top]
        pr = [pnear(h) for h in rnd]
        bt, br = B[:10], B[10:]
        for i, h in enumerate(top):
            check(absS[h] <= bt[i] * (1 + 1e-9) + 1e-12, "switching bound violated")
        for i, h in enumerate(rnd):
            check(absS[h] <= br[i] * (1 + 1e-9) + 1e-12, "switching bound violated")
        print(f"  ({p},{a}) M={M}: top-10 h: median |S|/C={np.median(absS[top]):.3f}, median B={np.median(bt):.3f}, "
              f"median P_near={np.median(pt):.2f};  60 random h: median |S|/C={np.median(absS[rnd]):.2e}, "
              f"median B={np.median(br):.2e}, median P_near={np.median(pr):.2f}")
        # monotone link: B vs P_near over all 70 h
        allp = np.array(pt + pr)
        allB = np.array(list(bt) + list(br))
        for lo, hi in ((0, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 1.01)):
            sel = (allp >= lo) & (allp < hi)
            if sel.any():
                print(f"        P_near in [{lo:.1f},{hi:.1f}): #h={int(sel.sum()):2d}  max B={allB[sel].max():.2e}")


def c3_sparse():
    hdr("C3b. sparse clocks 24 <= p <= 40 (no full spectrum): exact structured contribution, sampled minor arcs")
    gc = core.gates_core()
    rng = random.Random(5)
    print("  Sigma = {h : h 2^s 3^m = +-u (mod M), u in {1,5}, 0<=s<p, -a<=m<a};  E_Sig = (C + sum_Sig S(h))/M exactly (DP);")
    print("  minor arcs sampled at 150 random h off Sigma; RMS_abs = RMS of |S(h)| (the triangle route needs sup_off |S| < ~1)")
    print("   p  a     M (bits)   C (bits)   C/M        |Sig|  max_Sig|S|/C  E_Sig        RMS_off/sqrtC  max_off/sqrtC  RMS_abs")
    out = []
    for p in (24, 28, 32, 36, 40):
        acrit = int(p * math.log(2) / math.log(3))
        for a in sorted({p // 2, acrit, acrit + 1}):
            M = abs(gc.gate(p, a))
            C = math.comb(p, a)
            sig = set()
            for s in range(p):
                for m in range(-a, a):
                    g = core.mult(-s, -m, M)
                    for u in (1, 5):
                        for sgn in (1, -1):
                            sig.add((sgn * u * g) % M)
            sig.discard(0)
            sig = sorted(sig)
            Ssig = gc.S_dp(p, a, sig)
            ESig = (C + float(np.real(Ssig.sum()))) / M
            samp = []
            sset = set(sig)
            while len(samp) < 150:
                h = rng.randrange(1, M)
                if h not in sset:
                    samp.append(h)
            Ss = np.abs(gc.S_dp(p, a, samp))
            rms = float(np.sqrt((Ss ** 2).mean())) / math.sqrt(C)
            mx = float(Ss.max()) / math.sqrt(C)
            # requirement for the triangle-inequality route: sup_off |S| < ~1, floor ~ RMS ~ sqrt(C)
            need_floor = rms * math.sqrt(C)
            out.append((p, a, M, C, ESig))
            print(f"  {p:2d} {a:2d}   {math.log2(M):6.1f}     {math.log2(C):6.1f}   {C / M:9.2e}   {len(sig):5d}   {float(np.abs(Ssig).max()) / C:.3f}      "
                  f"{ESig:10.3e}   {rms:6.3f}         {mx:6.2f}        {need_floor:9.3g}")
            sys.stdout.flush()
            check(rms > 0.3, "minor-arc RMS unexpectedly small")
    print("  (all these clocks carry N = 0: gates-lane census, p <= 40, two independent methods)")
    print("  => N = E_Sig + (signed minor-arc sum); E_Sig is an exactly computable O(1) number, but excluding a")
    print("     cycle needs the signed minor-arc sum to be < 1 - |E_Sig| while its terms have RMS ~ sqrt(C) = 10^3..10^5.")


def main():
    core.lean_malloc()
    c1_c3_full()
    c1b_larger()
    c1c_named()
    c2_switching()
    c3_sparse()
    print(f"\n[Part C done in {time.time() - T0:.0f}s, peak RSS {core.mem_mb():.0f} MB]")
    check(core.mem_mb() < 500, "memory budget exceeded")
    print("PART C CHECKS PASSED")


if __name__ == "__main__":
    main()
