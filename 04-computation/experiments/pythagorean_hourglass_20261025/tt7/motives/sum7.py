# Fast exact evaluation of normalised seven-vertex Paley graph sums
#   S_H(q) = sum over injective x : {0..6} -> F_q with x_a = 0, x_b = 1 of prod_{(i<j) in H} chi(x_j - x_i)
# (order orientation; T_H = q(q-1) S_H when |E(H)| is even).  Every pair of free vertices carries either chi(x_w - x_v)
# (an arc) or the distinctness factor 1 - delta, so the network is a complete K5 on the free variables; it is
# contracted in O(q^5) flops with one large GEMM per value of the first free variable.  float64 is exact here
# (all partial sums are integers below q^5 < 2^53).  mac-mini-2026-10-25-pyth.
import itertools
import numpy as np

def chi_table(p):
    c = np.full(p, -1.0); c[0] = 0.0
    c[[(x * x) % p for x in range(1, p)]] = 1.0
    return c

def pair_factor_matrix(p, chi, is_arc):
    """matrix F[s, t] for a pair (i < j) with x_i = s, x_j = t: chi(t - s) for an arc, 1 - delta otherwise"""
    idx = np.arange(p)
    D = (idx[None, :] - idx[:, None]) % p
    return chi[D] if is_arc else 1.0 - np.eye(p)

def S7(p, edges, fixed=(0, 1)):
    chi = chi_table(p)
    E = {tuple(sorted(e)) for e in edges}
    a, b = fixed
    val = {a: 0, b: 1}
    const = 1.0
    if tuple(sorted((a, b))) in E:
        const = chi[(val[max(a, b)] - val[min(a, b)]) % p]
    free = [v for v in range(7) if v not in fixed]
    idx = np.arange(p)
    # unary factors from the two fixed vertices
    u = {}
    for v in free:
        w = np.ones(p)
        for f in fixed:
            c = val[f]
            arc = tuple(sorted((v, f))) in E
            if v < f:   # pair (v, f): chi(x_f - x_v) = chi(c - t)
                w = w * (chi[(c - idx) % p] if arc else (idx != c))
            else:       # pair (f, v): chi(x_v - x_f) = chi(t - c)
                w = w * (chi[(idx - c) % p] if arc else (idx != c))
        u[v] = w
    # binary matrices M[(v, w)][x_v, x_w]
    Mf = {}
    for v, w in itertools.combinations(free, 2):
        F = pair_factor_matrix(p, chi, (v, w) in E)      # rows x_v, cols x_w with v < w
        Mf[(v, w)] = F; Mf[(w, v)] = F.T
    f1, f2, f3, f4, f5 = free
    M12, M13, M14, M15 = Mf[(f1, f2)], Mf[(f1, f3)], Mf[(f1, f4)], Mf[(f1, f5)]
    M23, M24, M25 = Mf[(f2, f3)], Mf[(f2, f4)], Mf[(f2, f5)]
    M34, M35, M45 = Mf[(f3, f4)], Mf[(f3, f5)], Mf[(f4, f5)]
    tot = 0.0
    for x1 in range(p):
        if u[f1][x1] == 0: continue
        u2 = u[f2] * M12[x1]; u3 = u[f3] * M13[x1]; u4 = u[f4] * M14[x1]; u5 = u[f5] * M15[x1]
        U3 = u3[None, :] * M23            # (x2, x3)
        U4 = u4[None, :] * M24            # (x2, x4)
        U5 = u5[None, :] * M25            # (x2, x5)
        # K[x2, x3, x5] = sum_x4 M34[x3, x4] U4[x2, x4] M45[x4, x5]
        K = ((M34[None, :, :] * U4[:, None, :]).reshape(p * p, p) @ M45).reshape(p, p, p)
        K *= M35[None, :, :]
        K *= U3[:, :, None]
        K *= U5[:, None, :]
        tot += u[f1][x1] * float(u2 @ K.sum(axis=(1, 2)))
    return int(round(const * tot))

def S7_brute(p, edges, fixed=(0, 1)):
    chi = [0] + [1 if pow(x, (p - 1) // 2, p) == 1 else -1 for x in range(1, p)]
    E = {tuple(sorted(e)) for e in edges}
    a, b = fixed
    free = [v for v in range(7) if v not in fixed]
    tot = 0
    for vals in itertools.permutations([t for t in range(p) if t not in (0, 1)], 5):
        x = [None] * 7; x[a] = 0; x[b] = 1
        for v, t in zip(free, vals): x[v] = t
        pr = 1
        for i, j in E:
            pr *= chi[(x[j] - x[i]) % p]
            if pr == 0: break
        tot += pr
    return tot

if __name__ == '__main__':
    import json, sys, time
    G = json.load(open('genuine37.json'))
    for p in (7, 11):
        t0 = time.time()
        ok = all(S7(p, H['edges']) == S7_brute(p, H['edges']) for H in G)
        alt = all(S7(p, H['edges'], fixed=(2, 5)) == S7(p, H['edges']) for H in G)
        print(f"p = {p}: fast = brute for all 38 survivors: {ok}; independent of the fixed pair: {alt}  [{time.time() - t0:.1f}s]")
    t0 = time.time(); S7(151, G[0]['edges']); print(f"timing: one sum at q = 151 in {time.time() - t0:.2f}s")
