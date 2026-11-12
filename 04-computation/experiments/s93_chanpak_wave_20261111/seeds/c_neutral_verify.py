#!/usr/bin/env python3
"""S93 lane "seeds": independent verification of the neutral gadgets of order N (default 8) that violate HYP-9662 (D),
in particular those with kappa_X(P) = 0 and r_X(P) = +1 or -1 over a non-transitive X.

For each such gadget (X, P), found again here from cls_<N>.txt:
 (i)  literal enumeration over all N! vertex orders and all cut points: H(X), Pi_2(X), W_X(P), W3_X(P)
      -> kappa = Pi_2 - 2 W - 2 W3 and r_X = H - W;
 (ii) the gadget law itself: for random bases (B, O) of order 1..5 the tournament B+ = (B => X => closer,
      closer -> O) is built explicitly and the state (r, M, h, E) of the pair (B+, P) is counted directly by the subset
      recursion; compared with THM-5130 (A):  r+ = r_X r + kappa M (= r_X r),  M+ = u (r + h) + (4u - 2 sigma + 2 u3) M
      - sigma E,  h+ = H_X r + Pi_2 M,  E+ = eps (r - h) + nu M.
Writes the list of unit neutral gadgets with their constants (c_neutral_unit_<N>.txt).
"""
import sys

import numpy as np

from seedlib import dp_tables, pair_mats, out_from_colcode
from a_verify import literal
from c_neutral import gadget_tables


def tables(n, out):
    fe = np.zeros((1 << n, n), np.int64)
    gs = np.zeros((1 << n, n), np.int64)
    hS = np.zeros(1 << n, np.int64)
    dp_tables(n, out, fe, gs, hS)
    return fe, gs, hS


def state(n, out, O):
    fe, gs, hS = tables(n, out)
    Q = np.zeros((n, n), np.int64)
    U = np.zeros(n, np.int64)
    pair_mats(n, fe, gs, hS, Q, U)
    vs = [v for v in range(n) if (O >> v) & 1]
    full = (1 << n) - 1
    h = int(hS[full])
    W = sum(int(Q[e, f]) for e in vs for f in vs if e != f)
    return (h - W, sum(int(U[v]) for v in vs), h, sum(int(fe[full, v]) for v in vs))


def gadget_constants(n, out, P):
    fe, gs, hS = tables(n, out)
    full = (1 << n) - 1
    UR = np.zeros((1 << n, n), np.int64)
    K = np.zeros((n, n), np.int64)
    Wm = np.zeros((n, n), np.int64)
    pi2 = int(gadget_tables(n, fe, gs, hS, UR, K, Wm))
    vs = [v for v in range(n) if (P >> v) & 1]
    H = int(hS[full])
    w = sum(int(Wm[e, f]) for i, e in enumerate(vs) for f in vs[i + 1:])
    ksum = sum(int(K[e, f]) for i, e in enumerate(vs) for f in vs[i + 1:])
    kappa = pi2 - ksum
    sigma = sum(int(gs[full, v]) for v in vs)
    eps = sum(int(fe[full, v]) for v in vs)
    u = sum(int(UR[full, v]) for v in vs)
    # nu = V_X(P): pairs (P1 ending at o, P2 possibly empty)
    nu = sum(int(fe[S, o]) * int(hS[full ^ S]) for o in vs for S in range(1, full + 1) if (S >> o) & 1)
    # P_2(Y) for every subset Y: ordered pairs of nonempty paths covering Y
    u3 = 0
    for o in vs:
        for S in range(1, full):
            if (S >> o) & 1 and gs[S, o]:
                Y = full ^ S
                p2 = 0
                T = (Y - 1) & Y
                while T:
                    p2 += int(hS[T]) * int(hS[Y ^ T])
                    T = (T - 1) & Y
                u3 += int(gs[S, o]) * p2
    return dict(H=H, Pi=pi2, w=w, kappa=kappa, rX=H - w, u=u, sigma=sigma, eps=eps, nu=nu, u3=u3)


def closing_with(nb, outb, O, nx, outx):
    """B => X => closer, closer -> O, B - O -> closer; X on the vertices nb..nb+nx-1, closer last."""
    N = nb + nx + 1
    res = np.zeros(N, np.int64)
    xmask = ((1 << nx) - 1) << nb
    vbit = 1 << (nb + nx)
    for x in range(nb):
        res[x] = int(outb[x]) | xmask
        if not (O >> x) & 1:
            res[x] |= vbit
    for x in range(nx):
        res[nb + x] = (int(outx[x]) << nb) | vbit
    res[nb + nx] = O
    return N, res


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    rng = np.random.default_rng(8)
    reps = np.atleast_1d(np.loadtxt(f"cls_{n}.txt", dtype=np.int64))
    found = []
    allviol = 0
    for code in reps.tolist():
        out = np.zeros(n, np.int64)
        out_from_colcode(n, code, out)
        fe, gs, hS = tables(n, out)
        UR = np.zeros((1 << n, n), np.int64)
        K = np.zeros((n, n), np.int64)
        Wm = np.zeros((n, n), np.int64)
        pi2 = int(gadget_tables(n, fe, gs, hS, UR, K, Wm))
        H = int(hS[(1 << n) - 1])
        for P in range(1, 1 << n):
            vs = [v for v in range(n) if (P >> v) & 1]
            if len(vs) < 2 or len(vs) > 4:
                continue
            ks = sum(int(K[e, f]) for i, e in enumerate(vs) for f in vs[i + 1:])
            if ks == pi2:
                w = sum(int(Wm[e, f]) for i, e in enumerate(vs) for f in vs[i + 1:])
                if len(vs) != 2 or H - w != -H:
                    allviol += 1
                if abs(H - w) == 1 and H > 1:
                    found.append((code, P))
    print(f"order {n}: gadgets with kappa = 0 violating HYP-9662 (D) (|P| in 2..4): {allviol}; with |r_X| = 1 over a non-transitive X: {len(found)}")
    lines = ["# order colcode Pmask |P| H(X) Pi_2 W(P) r_X kappa u sigma u3 eps nu   (X = seedlib.out_from_colcode(order, colcode))"]
    ok_lit = True
    ok_law = True
    nlaw = 0
    for code, P in found:
        out = np.zeros(n, np.int64)
        out_from_colcode(n, code, out)
        g = gadget_constants(n, out, P)
        # (i) literal
        Q2 = np.zeros((n, n), np.int64)
        U2 = np.zeros(n, np.int64)
        au = np.zeros((n, n), np.int64)
        H3 = np.zeros((n, n), np.int64)
        H, pi2 = literal(n, out, Q2, U2, au, H3)
        vs = [v for v in range(n) if (P >> v) & 1]
        w = sum(int(Q2[e, f]) for e in vs for f in vs if e != f)
        w3 = sum(int(H3[e, f]) for e in vs for f in vs if e != f)
        if (H, pi2, w, pi2 - 2 * w - 2 * w3) != (g["H"], g["Pi"], g["w"], 0) or g["kappa"] != 0:
            ok_lit = False
        lines.append(f"{n} {code} {P} {len(vs)} {g['H']} {g['Pi']} {g['w']} {g['rX']} {g['kappa']} {g['u']} {g['sigma']} {g['u3']} {g['eps']} {g['nu']}")
        # (ii) the law on random bases
        for _ in range(6):
            nb = int(rng.integers(1, 6))
            codeb = int(rng.integers(0, 1 << (nb * (nb - 1) // 2))) if nb > 1 else 0
            outb = np.zeros(nb, np.int64)
            out_from_colcode(nb, codeb, outb)
            O = int(rng.integers(1, 1 << nb))
            r, M, h, E = state(nb, outb, O)
            N, res = closing_with(nb, outb, O, n, out)
            got = state(N, res, P << nb)
            pred = (g["rX"] * r + g["kappa"] * M,
                    g["u"] * (r + h) + (4 * g["u"] - 2 * g["sigma"] + 2 * g["u3"]) * M - g["sigma"] * E,
                    g["H"] * r + g["Pi"] * M,
                    g["eps"] * (r - h) + g["nu"] * M)
            if got != pred or got[0] != g["rX"] * r:
                ok_law = False
                print("LAW FAILS", code, P, (r, M, h, E), got, pred)
            nlaw += 1
    print(f"(i) literal enumeration of all {n}! vertex orders confirms kappa = 0 and r_X for all {len(found)}: {ok_lit}")
    print(f"(ii) gadget law on {nlaw} explicit closings (bases of order 1..5, orders up to {n + 6}): offset r+ = r_X r and the full state: {ok_law}")
    for ln in lines:
        print(ln)
    with open(f"c_neutral_unit_{n}.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    if not (ok_lit and ok_law):
        sys.exit(1)


if __name__ == "__main__":
    main()
