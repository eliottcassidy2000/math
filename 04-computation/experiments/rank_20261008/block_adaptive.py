#!/usr/bin/env python3
"""Adaptive block certificates (state-dependent block length).

A block-length rule k(h) that depends only on the current hidden state (M, e) mod d^K is a stopping rule, so the
Lyapunov argument of THM-4611 (2) applies to the walk sampled at tau_(n+1) = tau_n + k(h_(tau_n)) as soon as one form Q
balances every nonzero block matrix A_(k(h))(h mod d^(k(h))) actually used (block lengths <= KMAX, overshoot <= KMAX B).
Construction: start from all hidden states mod d^K0; a node whose block has margin < TAU under a pilot form is refined
(replaced by its lifts mod d^(s+1): M-lifts in H_(s+1), e-lifts e + d^s t), up to level KMAX.  The leaves partition the
hidden states, i.e. they define a rule k(h).  Then optimize Q over the leaf blocks and verify every leaf block exactly.
Only full tables up to level KMAX - 1 are stored; leaf blocks at higher levels come from the recursion
A_(s+1)(h) = d^s D(pi(h)) + sum_j A_s(h'_j).   Usage: python3 block_adaptive.py MAPNAME K0 KMAX TAU"""
import contextlib, io, math, sys, time
import numpy as np
from block_balance import int_coords, ratio_group, coupling_D, std_r
from block_balance2 import block_arrays, integer_forms
from block_lmi import ellipsoid_opt
from block_lmi2 import exact_check_vec, active_set_opt

MAPS = {
    'Z5_1_6_11_11_4': (5, [1, 6, 11, 11, 4]),
    'Z5_1_1_6_39_11': (5, [1, 1, 6, 39, 11]),
    'Z5_1_14_4_1_29': (5, [1, 14, 4, 1, 29]),
    'Z5_1_4_31_1_6': (5, [1, 4, 31, 1, 6]),
    'Z5_1_4_1_11_34': (5, [1, 4, 1, 11, 34]),
    'Z5_12371': (5, [1, 2, 3, 7, 1]),
    'Z5_11237': (5, [1, 1, 2, 3, 7]),
    'Z4_1357': (4, [1, 3, 5, 7]),
    'Z5_1_6_4_31_1': (5, [1, 6, 4, 31, 1]),
}

def margins_Q(Q, As):
    L = np.linalg.cholesky(Q); Li = np.linalg.inv(L)
    res = np.empty(len(As))
    for s in range(0, len(As), 400000):
        W = Li @ As[s:s + 400000].astype(float) @ Li.T
        ev = np.linalg.eigvalsh(W); tr = ev.sum(axis=1)
        safe = np.where(tr > 0, tr, 1.0)
        res[s:s + 400000] = np.where(tr > 0, (tr - 2 * ev[:, -1]) / safe, np.inf)
    return res

class Recursor:
    def __init__(self, d, m, r, v, rho):
        self.d, self.m, self.r, self.rho = d, m, r, rho
        self.D = np.zeros((d, d, rho, rho), dtype=np.int64)
        for a in range(1, d):
            if math.gcd(a, d) != 1: continue
            for b in range(d): self.D[a, b] = np.array(coupling_D(v, a, b, d, rho), dtype=np.int64)
        self.mm = np.array(m, dtype=np.int64); self.rr = np.array(r, dtype=np.int64)
    def next_level(self, Hs, As, s, Mt, Et):
        """Blocks A_(s+1) at nodes (Mt, Et) given mod d^(s+1), from the full level-s table (Hs, As)."""
        d = self.d; mod = d ** (s + 1); pmod = d ** s
        pos = -np.ones(pmod, dtype=np.int64); pos[Hs] = np.arange(len(Hs))
        inv = np.array([pow(x, -1, mod) for x in self.m], dtype=np.int64)
        a = Mt % d; b = Et % d
        A = self.D[a, b] * pmod
        for j in range(d):
            i = (a * j + b) % d
            t = (self.mm[i] * inv[j]) % mod; t = (t * Mt) % mod; t = (t * self.rr[j]) % mod
            N = (self.mm[i] * Et + self.rr[i] - t) % mod
            assert np.all(N % d == 0)
            e2 = (N // d) % pmod
            M2 = (((self.mm[i] * inv[j]) % mod) * Mt) % pmod
            idx = pos[M2]; assert np.all(idx >= 0)
            A = A + As[idx, e2]
        return A

def block_arrays_lean(d, m, r, v, rho, k, chunk=128):
    """Same table as block_balance2.block_arrays(d, m, r, v, rho, k), built level by level in int32 with chunked in-place
    accumulation (memory about half, no full-size temporaries).  Values are asserted to fit in int32."""
    R = Recursor(d, m, r, v, rho)
    Hp, Ap = None, None
    for s in range(1, k + 1):
        mod = d ** s; H = np.array(ratio_group(m, mod), dtype=np.int64); E = np.arange(mod, dtype=np.int64)
        A = np.empty((len(H), mod, rho, rho), dtype=np.int32)
        for c0 in range(0, len(H), chunk):
            Mt = H[c0:c0 + chunk][:, None] * np.ones((1, mod), dtype=np.int64); Et = np.ones((Mt.shape[0], 1), dtype=np.int64) * E[None, :]
            if s == 1:
                blk = R.D[Mt % d, Et % d]
            else:
                blk = R.next_level(Hp, Ap, s - 1, Mt, Et)
            assert np.abs(blk).max() < 2 ** 31
            A[c0:c0 + chunk] = blk
        Hp, Ap = H, A
    return Hp, Ap

def adaptive_certify(d, m, r, K0, KMAX, TAU, rounds=4, log=print, lean=False, chunk=20000):
    """Return (Qi, nleaf_blocks, counts, leaf_margin) for an exactly verified adaptive certificate, or None."""
    rho, v = int_coords(m)
    T = {s: (block_arrays_lean(d, m, r, v, rho, s) if lean else block_arrays(d, m, r, v, rho, s)) for s in range(K0, KMAX)}
    R = Recursor(d, m, r, v, rho)
    H0, A0 = T[K0]
    flat0 = A0.reshape(-1, rho, rho)
    nz0 = np.any(flat0.reshape(len(flat0), -1) != 0, axis=1)
    U0 = np.unique(flat0[nz0].reshape(-1, rho * rho), axis=0).reshape(-1, rho, rho).astype(float)
    with contextlib.redirect_stdout(io.StringIO()):
        P0, mg0 = active_set_opt(U0 / np.trace(U0, axis1=1, axis2=2)[:, None, None], rho)
    Q = np.linalg.inv(P0)
    log(f"   pilot: level-{K0} optimum {mg0:+.5f}")
    for rnd in range(rounds):
        leaves = []; counts = []
        mg = margins_Q(Q, flat0).reshape(A0.shape[0], A0.shape[1])
        bad = mg < TAU
        gb = A0[~bad]; leaves.append(gb[np.any(gb.reshape(len(gb), -1) != 0, axis=1)])
        idx = np.argwhere(bad)
        Mf = H0[idx[:, 0]]; Ef = idx[:, 1].astype(np.int64)
        counts.append(len(Mf))
        for s in range(K0, KMAX):
            if len(Mf) == 0: break
            Hs, As = T[s]
            Hn = np.array(ratio_group(m, d ** (s + 1)), dtype=np.int64)
            red = Hn % d ** s
            order = np.argsort(red); red_sorted = red[order]
            final = (s + 1 == KMAX)
            newM = []; newE = []
            # process the frontier in chunks; lifts of a chunk are built, evaluated and reduced to unique leaves at once
            for c0 in range(0, len(Mf), chunk):
                Mc, Ec = Mf[c0:c0 + chunk], Ef[c0:c0 + chunk]
                lo = np.searchsorted(red_sorted, Mc, side='left'); hi = np.searchsorted(red_sorted, Mc, side='right')
                Mt_list = []; Et_list = []
                for t in range(d):
                    for off in range(int((hi - lo).max())):
                        sel = lo + off < hi
                        Mt_list.append(Hn[order[(lo + off)[sel]]]); Et_list.append(Ec[sel] + d ** s * t)
                Mt = np.concatenate(Mt_list); Et = np.concatenate(Et_list)
                blocks = R.next_level(Hs, As, s, Mt, Et)
                mgn = margins_Q(Q, blocks)
                badn = (mgn < TAU) & (not final)
                gb = blocks[~badn]
                gb = gb[np.any(gb.reshape(len(gb), -1) != 0, axis=1)]
                if len(gb): leaves.append(np.unique(gb.reshape(-1, rho * rho), axis=0))
                newM.append(Mt[badn]); newE.append(Et[badn])
            Mf = np.concatenate(newM) if newM else np.zeros(0, dtype=np.int64)
            Ef = np.concatenate(newE) if newE else np.zeros(0, dtype=np.int64)
            counts.append(len(Mf))
        U = np.unique(np.concatenate([x.reshape(-1, rho * rho) for x in leaves]), axis=0).reshape(-1, rho, rho)
        Uf = U.astype(float); Uf /= np.trace(Uf, axis1=1, axis2=2)[:, None, None]
        with contextlib.redirect_stdout(io.StringIO()):
            P, mgl = active_set_opt(Uf, rho)
        log(f"   round {rnd}: nodes refined per level {counts}; {len(U)} distinct leaf blocks; leaf optimum {mgl:+.5f}")
        Qn = np.linalg.inv(P)
        if mgl > 0:
            for Qi in integer_forms(Qn, rho):
                ok, unsure = exact_check_vec(Qi, U)
                if ok: return Qi, len(U), counts, mgl
            log("   leaf optimum positive but rationalization failed")
        Q = Qn
    return None

if __name__ == '__main__':
    name, K0, KMAX, TAU = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
    d, m = MAPS[name]; r = std_r(m)
    rho, v = int_coords(m)
    Lam = sum(math.log(x / d) for x in m) / d
    t0 = time.time()
    print(f"{name} m = {m} r = {r} rank {rho} Lambda = {Lam:+.4f}", flush=True)
    res = adaptive_certify(d, m, r, K0, KMAX, TAU, log=lambda x: print(x, flush=True))
    if res:
        Qi, nU, counts, mgl = res
        print(f"   EXACT: integer form Q = {Qi} balances all {nU} distinct leaf blocks (block lengths {K0}..{KMAX}; nodes refined per level {counts})  [{time.time() - t0:.1f}s]", flush=True)
    else:
        print("   no certificate found", flush=True)
