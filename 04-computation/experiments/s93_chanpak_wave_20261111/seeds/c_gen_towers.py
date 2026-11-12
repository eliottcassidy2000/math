#!/usr/bin/env python3
"""S93 lane "seeds", deliverable (c): generalized towers = towers whose levels are ANY unit neutral gadgets.

A unit neutral gadget is (X, P) with kappa_X(P) = 0 and rho = r_X(P) = +1 or -1; then r+ = rho r for every base.
Padding: X' = TT_(i-1) => X => TT_m (i - 1 sources, m sinks), P unchanged, is again unit neutral with the same rho and
    H' = H,  Pi' = 2^(i-1+m) Pi,  u' = 2^m u,  sigma' = [i = 1] sigma,  eps' = [m = 0] eps,  nu' = 2^(i-1) nu,
    u3' = (2u - sigma + u3) 2^(i-1) 3^m - 2^(m+1) u + [i = 1] sigma          (PROOFS.txt, Lemma G)
so the level map of the padded gadget is
    r+ = rho r,  h+ = H r + 2^(i-1+m) Pi M,  M+ = 2^m u (r + h) + c 2^(i-1) 3^m M - [i = 1] sigma E,
    E+ = [m = 0] eps (r - h) + 2^(i-1) nu M,        c = 2 (2u - sigma + u3).
(For the transitive pairs: H = 1, Pi = 2^(delta+1), u = alpha, sigma = eps = 1, nu = alpha, c = 2 gamma.)

G1  unit neutral gadgets of order 8 and 9 (from cls_8.txt, cls_9.txt): primitive ones (no source / sink of X outside
    P), their constants; the padding formulas checked on the padded gadgets found at order 9.
G2  exponent bound sigma_0(G): root of sum over gadgets and (i, m) of (c 2^(i-1) 3^m - [i = 1] sigma)^(-s) = 1.
G3  exact enumeration of the generalized tower states with slope <= 2^K over the transitive bases: distinct slopes
    and values, compared with the pure towers.
G4  one explicit new seed: the sign-preserving gadget of order 8 put on the 3-cycle seed, counted directly.
Usage: python3 c_gen_towers.py K     (default 26)
"""
import math
import sys

import numpy as np
from numba import njit

from seedlib import out_from_colcode, closed_out, hpaths
from c_neutral import scan as neutral_scan, NCH, MAXREC
from c_neutral_verify import gadget_constants, closing_with, state


def unit_gadgets(n):
    reps = np.atleast_1d(np.loadtxt(f"cls_{n}.txt", dtype=np.int64))
    rec = np.zeros((NCH, MAXREC, 4), np.int64)
    nrec = np.zeros(NCH, np.int64)
    signs = np.zeros((NCH, 3), np.int64)
    extremes = np.zeros((NCH, 2), np.float64)
    neutral_scan(n, reps, rec, nrec, signs, extremes)
    res = []
    for ch in range(NCH):
        for c in range(int(nrec[ch])):
            code, P, rX, H = (int(x) for x in rec[ch, c])
            if abs(rX) == 1 and H > 1:
                res.append((code, P))
    return sorted(res)


def is_padded(n, out, P):
    full = (1 << n) - 1
    for v in range(n):
        if (P >> v) & 1:
            continue
        if int(out[v]) == full ^ (1 << v) or int(out[v]) == 0:
            return True
    return False


def rho_fn(gadgets, s, terms=80):
    tot = 0.0
    for c, sg in gadgets:
        for m in range(terms):
            tot += (c * 3.0 ** m - sg) ** (-s)
        tot += c ** (-s) * 2.0 ** (-s) / (1 - 2.0 ** (-s)) / (1 - 3.0 ** (-s))
    return tot


def root(f):
    lo, hi = 0.3, 1.5
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 1:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


@njit(cache=True)
def enumerate_gen(X, G, outM, outr):
    """G[g] = (rho, H, Pi, u, sigma, c, eps, nu).  All generalized tower states with M <= X over the transitive bases."""
    ng = G.shape[0]
    MAXD = 64
    sr = np.zeros(MAXD, np.int64)
    sM = np.zeros(MAXD, np.int64)
    sh = np.zeros(MAXD, np.int64)
    sE = np.zeros(MAXD, np.int64)
    sg = np.zeros(MAXD, np.int64)
    si = np.zeros(MAXD, np.int64)
    sm = np.zeros(MAXD, np.int64)
    cap = outM.shape[0]
    cnt = 0
    nb = 0
    g0 = 0
    while (1 << g0) <= X:
        nb += 1
        g0 += 1
    for b in range(nb + 2 * 56):
        if b < nb:
            r0, M0, h0, E0 = 1, 1 << b, 1, (1 if b == 0 else 0)
        else:
            bb = b - nb
            m0 = bb >> 1
            M0 = (3 if (bb & 1) == 0 else 5) * (1 << m0)
            if M0 > X:
                continue
            r0, h0, E0 = -1, 1, (1 if m0 == 0 else 0)
        if cnt < cap:
            outM[cnt] = M0
            outr[cnt] = r0
        cnt += 1
        depth = 0
        sr[0], sM[0], sh[0], sE[0] = r0, M0, h0, E0
        sg[0], si[0], sm[0] = 0, 1, 0
        while depth >= 0:
            g = sg[depth]
            if g >= ng:
                depth -= 1
                continue
            i = si[depth]
            m = sm[depth]
            M = sM[depth]
            p3 = 1
            for _ in range(m):
                p3 *= 3
            c = G[g, 5]
            lam = c * (1 << (i - 1)) * p3 - (G[g, 4] if i == 1 else 0)
            if lam > X // M + 1 or lam * M > X:
                if m == 0:
                    sg[depth] = g + 1
                    si[depth] = 1
                else:
                    si[depth] = i + 1
                sm[depth] = 0
                continue
            sm[depth] = m + 1
            r = sr[depth]
            h = sh[depth]
            E = sE[depth]
            M2 = (1 << m) * G[g, 3] * (r + h) + c * (1 << (i - 1)) * p3 * M - (G[g, 4] * E if i == 1 else 0)
            if M2 > X:
                continue
            h2 = G[g, 1] * r + (1 << (i - 1 + m)) * G[g, 2] * M
            E2 = (G[g, 6] * (r - h) if m == 0 else 0) + (1 << (i - 1)) * G[g, 7] * M
            r2 = G[g, 0] * r
            if cnt < cap:
                outM[cnt] = M2
                outr[cnt] = r2
            cnt += 1
            depth += 1
            sr[depth], sM[depth], sh[depth], sE[depth] = r2, M2, h2, E2
            sg[depth], si[depth], sm[depth] = 0, 1, 0
    return cnt


def count_values(arr, K):
    b = 1 << K
    sub = arr[arr <= b]
    acc = [sub]
    cur = sub
    while True:
        cur = cur[cur <= b // 2] * 2
        if cur.shape[0] == 0:
            break
        acc.append(cur)
    return int(np.unique(np.concatenate(acc)).shape[0])


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 26
    lines = []

    def P(s=""):
        print(s, flush=True)
        lines.append(s)

    trans = [dict(rho=-1, H=1, Pi=4, u=3, sigma=1, c=10, eps=1, nu=3, tag="TT_2 {1,2}"),
             dict(rho=-1, H=1, Pi=8, u=5, sigma=1, c=26, eps=1, nu=5, tag="TT_3 {1,3}")]
    prim = {}
    for n in (8, 9):
        ug = unit_gadgets(n)
        consts = []
        npad = 0
        okpad = True
        for code, Pm in ug:
            out = np.zeros(n, np.int64)
            out_from_colcode(n, code, out)
            g = gadget_constants(n, out, Pm)
            g["c"] = 2 * (2 * g["u"] - g["sigma"] + g["u3"])
            g["code"], g["P"], g["n"] = code, Pm, n
            if is_padded(n, out, Pm):
                npad += 1
                # check the padding formulas against some primitive gadget of order n - 1
                match = False
                for q in prim.get(n - 1, []):
                    for (i, m) in ((2, 0), (1, 1)):
                        pred = (q["H"], (1 << (i - 1 + m)) * q["Pi"], (1 << m) * q["u"], q["sigma"] if i == 1 else 0,
                                q["eps"] if m == 0 else 0, (1 << (i - 1)) * q["nu"],
                                (2 * q["u"] - q["sigma"] + q["u3"]) * (1 << (i - 1)) * 3 ** m - (1 << (m + 1)) * q["u"] + (q["sigma"] if i == 1 else 0),
                                q["rX"])
                        if pred == (g["H"], g["Pi"], g["u"], g["sigma"], g["eps"], g["nu"], g["u3"], g["rX"]):
                            match = True
                if not match:
                    okpad = False
            else:
                consts.append(g)
        prim[n] = consts
        P(f"G1  order {n}: unit neutral gadgets over non-transitive X: {len(ug)} (labelled P of the class representatives); padded {npad}"
          + (f" (constants = padding formulas of an order-{n - 1} primitive gadget: {okpad})" if npad else "") + f"; primitive {len(consts)}")
        byrho = {}
        for g in consts:
            byrho.setdefault(g["rX"], []).append(g)
        for rho in sorted(byrho):
            gs = byrho[rho]
            P(f"      rho = {rho:+d}: {len(gs)} gadgets; |P| in {sorted(set(bin(g['P']).count('1') for g in gs))}; H(X) in [{min(g['H'] for g in gs)}, {max(g['H'] for g in gs)}]; "
              f"c in [{min(g['c'] for g in gs)}, {max(g['c'] for g in gs)}]")
        with open(f"c_gen_gadgets_{n}.txt", "w") as fh:
            fh.write("# order colcode Pmask rho H Pi u sigma u3 eps nu c\n")
            for g in consts:
                fh.write(f"{n} {g['code']} {g['P']} {g['rX']} {g['H']} {g['Pi']} {g['u']} {g['sigma']} {g['u3']} {g['eps']} {g['nu']} {g['c']}\n")

    sets = [("transitive pairs only", [(10, 1), (26, 1)])]
    sets.append(("+ order 8", sets[0][1] + [(g["c"], g["sigma"]) for g in prim[8]]))
    sets.append(("+ orders 8, 9", sets[1][1] + [(g["c"], g["sigma"]) for g in prim[9]]))
    for name, gs in sets:
        s0 = root(lambda s: rho_fn(gs, s))
        P(f"G2  gadget set '{name}' ({len(gs)} primitive gadgets): sigma_0 = {s0:.6f};  rho(1) = {rho_fn(gs, 1.0):.6f}")
    for n in (8, 9):
        tot1 = sum(rho_fn([(g["c"], g["sigma"])], 1.0) for g in prim[n])
        P(f"      contribution of the primitive gadgets of order {n} to rho(1): {tot1:.6f};  to rho(0.8): {sum(rho_fn([(g['c'], g['sigma'])], 0.8) for g in prim[n]):.6f}")

    # G3
    X = 1 << K
    for name, glist in (("transitive pairs only", trans), ("+ order 8", trans + prim[8]), ("+ orders 8, 9", trans + prim[8] + prim[9])):
        G = np.array([[g["rho"] if "rho" in g else g["rX"], g["H"], g["Pi"], g["u"], g["sigma"], g["c"], g["eps"], g["nu"]] for g in glist], np.int64)
        cap = enumerate_gen(X, G, np.zeros(0, np.int64), np.zeros(0, np.int8))
        outM = np.zeros(cap, np.int64)
        outr = np.zeros(cap, np.int8)
        n_ = enumerate_gen(X, G, outM, outr)
        Mp = np.unique(outM[outr == 1])
        Mm = np.unique(outM[outr == -1])
        P(f"G3  '{name}': X = 2^{K}: states {n_}; distinct slopes +: {Mp.shape[0]}, -: {Mm.shape[0]}; distinct values 2^t M <= X: +: {count_values(Mp, K)}, -: {count_values(Mm, K)}")

    # G4: the sign-preserving gadget on the seed (C_3, two vertices)
    keep = [g for g in prim[8] if g["rX"] == 1]
    if keep:
        g = keep[0]
        outb = np.zeros(3, np.int64)
        outb[0], outb[1], outb[2] = 1 << 1, 1 << 2, 1 << 0          # 0 -> 1 -> 2 -> 0
        O = 0b011
        st0 = state(3, outb, O)
        outx = np.zeros(8, np.int64)
        out_from_colcode(8, g["code"], outx)
        N, res = closing_with(3, outb, O, 8, outx)
        st1 = state(N, res, g["P"] << 3)
        pred = (g["rX"] * st0[0], g["u"] * (st0[0] + st0[2]) + g["c"] * st0[1] - g["sigma"] * st0[3], g["H"] * st0[0] + g["Pi"] * st0[1],
                g["eps"] * (st0[0] - st0[2]) + g["nu"] * st0[1])
        vals = []
        for t in (1, 2):
            res2 = np.zeros(N + t + 1, np.int64)
            N2 = closed_out(N, res, g["P"] << 3, t, res2)
            vals.append(int(hpaths(N2, res2)))
        P(f"G4  seed (C_3, two vertices) with state {st0}; the sign-preserving gadget of order 8 (H = {g['H']}, colcode {g['code']}, P mask {g['P']}) on top: "
          f"state of the new pair (order {N}) counted directly {st1}, gadget law {pred}; closed tournaments t = 1, 2: H = {vals} = 1 + 2^t * {st1[1]}: "
          f"{vals == [1 + 2 * st1[1], 1 + 4 * st1[1]] and st1 == pred}")
    with open(f"c_gen_towers_{K}.out", "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
