"""(*) at gamma=3/5 for large R: beam search controlling SEVERAL low-order
residual coefficients per step instead of clamping each one greedily.

Recursion: sigma_{-1}=q^{R-1};  p*sigma_i = sigma_{i-1} - Delta_i.
delta_{d}   is forced (= sigma_{i-1}(0), must be +-1).
delta_{d-1} sets sigma_i(0);  delta_{d-2} sets sigma_i(1); etc.
So the first CTRL coefficients of the next residual are directly steerable.
We enumerate small targets for them and keep the beam-best states by |sigma|.
"""
import sys
from math import comb
sys.setrecursionlimit(10000)

def basis_poly(d, k):
    r = [0]*(d+1)
    for s in range(k+1):
        r[d-k+s] += comb(k, s)*(-1)**s
    return r

def qpow(n):
    r = [1]
    for _ in range(n):
        nr = [0]*(len(r)+1)
        for i, c in enumerate(r):
            nr[i] += c; nr[i+1] -= c
        r = nr
    return r

def trim(a):
    a = list(a)
    while a and a[-1] == 0: a.pop()
    return a

def step(sigma, d, targets):
    """targets: tuple of desired sigma_next[j] for j=0..len-1 (None = clamp)."""
    if not sigma or abs(sigma[0]) != 1: return None
    res = list(sigma) + [0]*(d+4)
    deltas = [0]*(d+1)
    v = sigma[0]; deltas[d] = v
    for t, c in enumerate(basis_poly(d, d)): res[t] -= v*c
    for j, tgt in enumerate(targets):
        k = d-1-j
        if k < 0: return None
        cap = comb(d, k)
        want = res[j+1] - (tgt if tgt is not None else 0)
        if tgt is None:
            want = max(-cap, min(cap, res[j+1]))
            if (want-cap) % 2: want = want-1 if want-1 >= -cap else want+1
        if abs(want) > cap or (want-cap) % 2: return None
        deltas[k] = want
        for t, c in enumerate(basis_poly(d, k)): res[t] -= want*c
    for k in range(d-1-len(targets), -1, -1):
        cap = comb(d, k); want = res[d-k]
        vv = max(-cap, min(cap, want))
        if (vv-cap) % 2: vv = vv-1 if vv-1 >= -cap else vv+1
        deltas[k] = vv
        if vv:
            for t, c in enumerate(basis_poly(d, k)): res[t] -= vv*c
    if res[0] != 0: return None
    return deltas, trim(res[1:])

def solve(R, g1=3, g2=5, D0=0, beam=250, ctrl=2, span=2):
    d = [(g1*(R+i))//g2 + D0 for i in range(R)]
    states = [([], qpow(R-1))]
    opts = [None] + [x for x in range(-span, span+1)]
    from itertools import product
    for i in range(R-1):
        nxt = []
        for acc, sig in states:
            for tg in product(opts, repeat=ctrl):
                if tg[0] is None: continue          # sigma(0) must be +-1
                if tg[0] not in (1, -1): continue
                r = step(sig, d[i], tg)
                if r is None: continue
                de, ns = r
                if not ns or abs(ns[0]) != 1: continue
                nxt.append((acc+[de], ns))
        if not nxt: return None, f"died at row {i}"
        nxt.sort(key=lambda s: (len(s[1]), sum(abs(c) for c in s[1][:6])))
        seen, uniq = set(), []
        for a, s in nxt:
            key = tuple(s[:10])
            if key in seen: continue
            seen.add(key); uniq.append((a, s))
            if len(uniq) >= beam: break
        states = uniq
    di = d[R-1]
    for acc, sig in states:
        if len(sig)-1 > di: continue
        res = list(sig)+[0]*(di+2); de = [0]*(di+1); ok = True
        for k in range(di, -1, -1):
            cap = comb(di, k); want = res[di-k]
            if abs(want) > cap or (want-cap) % 2: ok = False; break
            de[k] = want
            for t, c in enumerate(basis_poly(di, k)): res[t] -= want*c
        if ok and not trim(res): return acc+[de], "SOLVED"
    return None, f"beam {beam} exhausted at final row"

def verify(R, sol, g1=3, g2=5, D0=0):
    d = [(g1*(R+i))//g2 + D0 for i in range(R)]
    acc = [0]*(3*R+8)
    for i, de in enumerate(sol):
        for k, v in enumerate(de):
            assert abs(v) <= comb(d[i], k) and (v-comb(d[i], k)) % 2 == 0
            if v:
                for t, c in enumerate(basis_poly(d[i], k)): acc[i+t] += v*c
    return trim(acc) == trim(qpow(R-1))

if __name__ == "__main__":
    for R in [8, 16, 32, 64]:
        sol, msg = solve(R)
        ok = verify(R, sol) if sol else False
        print(f"gamma=3/5 R={R:4d}: {msg}" + (f" | VERIFY {'PASS' if ok else 'FAIL'}" if sol else ""), flush=True)
