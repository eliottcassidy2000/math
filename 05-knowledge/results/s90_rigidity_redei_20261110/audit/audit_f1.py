#!/usr/bin/env python3
"""Independent audit of HYP-9657 (A), (B): the F1 family, its representation count and the class-mean pattern.
Own code, starting from the CANON coordinates of THM-4875 (D) / THM-4885:

      F1(u, v) = (1 + u)(1 + v) + g(u) g(v) - e(u) e(v),   u, v >= 0,
      g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2),   e(v) = [v = 0 or v even],

not from the packet's psi-coordinates.  Steps:
 (T)  brute force on actual tournaments: the apex triple (H, e_apex, R_apex) of the one-apex tournament T(v) equals
      (g(v), e(v), 1 + v) for v < 64, and H(T(u) => T(v)^op with the apex arc reversed) = F1(u, v) for u, v < 16
      (depth-first path counting).
 (C)  coordinate change: g(2h) = 2 psi(h) + 1 (psi as in the packet) for h < 2^20, so that with
      u = (2 h1 + 1) 2^s1 - 1, v = (2 h2 + 1) 2^s2 - 1:  F1(u, v) = (2psi(h1)+1)(2psi(h2)+1) + 2^(s1+s2)(2h1+1)(2h2+1) - [s1+s2 = 0].
      The packet's r(N) counts triples (h1, h2, s) with s = s1 + s2, i.e. ORDERED pairs (u, v) with v even (s2 = 0);
      THM-4885 (B)/(C) counts UNORDERED pairs u <= v.  Both are computed here:
            r_tri(N)  = #{(u, v) ordered, v even : F1(u, v) = N}      (= the packet's r)
            r_pair(N) = #{u <= v : F1(u, v) = N}                      (THM-4885's multiplicity)
      and r_tri is cross-checked against a direct enumeration of (h1, h2, s) below 2^15.
 (A)  class means of r modulo p, for both multiplicities, with standard errors; the spread of the classes other than
      0 and -1; split by the Kronecker symbol (-N/p).
 (B)  r at the Heegner numbers; the N = 4q - 1 < 3000 with N, q prime and r = 0; mean and miss rate on
      "N prime and (N+1)/4 prime" and on "15 | N(N+1)".
Usage: python3 audit_f1.py [K=21] [lo=18]     (window: odd N in [2^lo, 2^K))
"""
import sys
import time

import numpy as np


def hp_tables(out):
    """cnt[mask][start] = number of directed paths with vertex set mask starting at start (depth-first)."""
    n = len(out)
    cnt = [[0] * n for _ in range(1 << n)]
    for s in range(n):
        stack = [(s, 1 << s)]
        while stack:
            v, used = stack.pop()
            cnt[used][s] += 1
            c = out[v] & ~used
            while c:
                low = c & -c
                c ^= low
                stack.append((low.bit_length() - 1, used | low))
    return cnt


def hp_count(out):
    n = len(out)
    full = (1 << n) - 1
    total = 0
    stack = [(v, 1 << v) for v in range(n)]
    while stack:
        v, used = stack.pop()
        if used == full:
            total += 1
            continue
        c = out[v] & ~used
        while c:
            low = c & -c
            c ^= low
            stack.append((low.bit_length() - 1, used | low))
    return total


def one_apex(v):
    """chain 0 -> 1 -> ... -> m-1 (i -> j for i < j), apex = m, apex -> i iff bit i of the word (leading bit first)."""
    bits = [int(b) for b in bin(v)[2:]] if v else []
    m = len(bits)
    out = [0] * (m + 1)
    for i in range(m):
        for j in range(i + 1, m):
            out[i] |= 1 << j
        if bits[i]:
            out[m] |= 1 << i
        else:
            out[i] |= 1 << m
    return out, m          # apex index m


def g_e_arrays(V):
    g = np.zeros(V, dtype=np.int64)
    g[0] = 1
    k = 0
    while (1 << k) < V:
        lo, hi = 1 << k, min(1 << (k + 1), V)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        g[lo:hi] = np.where(idx & 1 == 1, g[half], g[half] + half + (half & 1))
        k += 1
    e = ((np.arange(V) & 1) == 0).astype(np.int64)
    return g, e


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 21
    LO = int(sys.argv[2]) if len(sys.argv) > 2 else 18
    X = 1 << K
    t0 = time.time()
    ok = True
    g, e = g_e_arrays(X + 2)
    # scalar definition of g on a prefix, as a check of the vectorised recursion
    gs = [1]
    for v in range(1, 5000):
        gs.append(gs[v >> 1] if v & 1 else gs[v >> 1] + (v >> 1) + ((v >> 1) & 1))
    ok &= gs == [int(x) for x in g[:5000]]

    # ---- (T) tournaments
    bad = 0
    triples = {}
    for v in range(64):
        out, apex = one_apex(v)
        n = len(out)
        cnt = hp_tables(out)
        full = (1 << n) - 1
        h = sum(cnt[full])
        conv = [0] * n
        for i in range(n):
            for j in range(n):
                if out[i] >> j & 1:
                    conv[j] |= 1 << i
        cntc = hp_tables(conv)                      # paths of the converse starting at t = paths ending at t
        e_apex = cntc[full][apex]
        R = 0
        for Xm in range(1 << n):
            if Xm >> apex & 1:
                continue
            hx = sum(cnt[Xm]) if Xm else 1
            R += hx * cnt[full ^ Xm][apex]
        triples[v] = (h, e_apex, R)
        if (h, e_apex, R) != (int(g[v]), int(e[v]), 1 + v):
            bad += 1
    print(f"(T) apex triple (H, e, R) of T(v) = (g(v), e(v), 1 + v) for v < 64: mismatches {bad}", "OK" if bad == 0 else "FAIL")
    ok &= bad == 0
    bad = 0
    for u in range(16):
        A, a = one_apex(u)
        for v in range(16):
            B, b = one_apex(v)
            nA, nB = len(A), len(B)
            # J = A => B^op with the arc b -> a reversed
            J = [0] * (nA + nB)
            for i in range(nA):
                J[i] = A[i] | (((1 << nB) - 1) << nA)
            for i in range(nB):
                for j in range(nB):
                    if B[j] >> i & 1:              # j -> i in B  =>  i -> j in B^op
                        J[nA + i] |= 1 << (nA + j)
            J[a] &= ~(1 << (nA + b))
            J[nA + b] |= 1 << a
            want = (1 + u) * (1 + v) + int(g[u]) * int(g[v]) - int(e[u]) * int(e[v])
            if hp_count(J) != want:
                bad += 1
    print(f"(T) H(T(u) => T(v)^op, apex arc reversed) = F1(u, v) for u, v < 16 (256 joins, orders up to 10): mismatches {bad}",
          "OK" if bad == 0 else "FAIL", f"[{time.time()-t0:.0f}s]")
    ok &= bad == 0

    # ---- (C) coordinates
    H2 = 1 << 20
    psi = np.zeros(H2, dtype=np.int64)
    for h in range(1, H2):
        psi[h] = psi[h >> 1] + ((h >> 1) if h % 2 == 0 else ((h >> 1) // 2 + 1))
    good = bool(np.all(g[0:2 * H2:2][:H2] == 2 * psi + 1)) if 2 * H2 <= X + 2 else bool(
        np.all(g[0:X + 2:2] == 2 * psi[:len(g[0:X + 2:2])] + 1))
    print("(C) g(2h) = 2 psi(h) + 1 on the computed range:", "OK" if good else "FAIL")
    ok &= good

    # ---- representation counts
    r_pair = np.zeros(X + 1, dtype=np.int32)
    r_tri = np.zeros(X + 1, dtype=np.int32)
    u = 0
    even_values = 0
    while (1 + u) * (1 + u) <= X:
        vmax = X // (1 + u) - 1
        vs = np.arange(u, vmax + 1, dtype=np.int64)
        N = (1 + u) * (1 + vs) + int(g[u]) * g[u:vmax + 1] - int(e[u]) * e[u:vmax + 1]
        keep = N <= X
        N = N[keep]
        vs = vs[keep]
        even_values += int(np.count_nonzero(N % 2 == 0))
        r_pair += np.bincount(N, minlength=X + 1).astype(np.int32)
        w = ((vs & 1) == 0).astype(np.int64)                       # ordered (u, v), v even
        w2 = np.where(vs != u, 1 - (u & 1), 0)                     # ordered (v, u), u even, v != u
        r_tri += np.bincount(N, weights=w + w2, minlength=X + 1).astype(np.int32)
        u += 1
    print(f"    F1 values up to 2^{K}: even values {even_values} (expected 0)", f"[{time.time()-t0:.0f}s]")
    ok &= even_values == 0

    # direct (h1, h2, s) enumeration below 2^15 as a check of r_tri
    Xs = 1 << 15
    rd = np.zeros(Xs + 1, dtype=np.int64)
    hmax = Xs // 2 + 1
    gg = 2 * psi[:hmax] + 1
    yy = 2 * np.arange(hmax, dtype=np.int64) + 1
    for s in range(0, 16):
        for h1 in range(hmax):
            y1 = int(yy[h1]) << s
            if y1 > Xs:
                break
            Nv = int(gg[h1]) * gg + y1 * yy - (1 if s == 0 else 0)
            Nv = Nv[Nv <= Xs]
            if len(Nv):
                rd += np.bincount(Nv, minlength=Xs + 1)
    good = bool(np.all(rd == r_tri[:Xs + 1]))
    print("(C) r_tri (ordered pairs with v even) = direct count of triples (h1, h2, s), all N <= 2^15:", "OK" if good else "FAIL")
    ok &= good
    good = bool(np.all((r_pair > 0) == (r_tri > 0)))
    print("    zero sets of r_pair and r_tri coincide:", "OK" if good else "FAIL")
    ok &= good

    odd = np.arange(1, X + 1, 2)
    win = odd[odd >= (1 << LO)]
    print(f"window: odd N in [2^{LO}, 2^{K}): {len(win)} integers; first misses: {odd[r_tri[odd] == 0][:40].tolist()}")
    for name, r in (("r_tri (packet)", r_tri), ("r_pair (THM-4885)", r_pair)):
        rw = r[win].astype(np.float64)
        mean = rw.mean()
        sd = rw.std()
        miss = int((rw == 0).sum())
        print(f"--- {name}: mean {mean:.4f}, sd {sd:.3f}, misses {miss}, miss rate {miss/len(win):.5f}")
        for p in (3, 5, 7, 11, 13, 17, 19, 23):
            res = win % p
            means = np.array([rw[res == c].mean() for c in range(p)])
            ns = np.array([(res == c).sum() for c in range(p)])
            se = sd / np.sqrt(ns.mean())
            order = np.argsort(-means)
            top2 = set(int(x) for x in order[:2])
            others = [c for c in range(1, p - 1)]
            om = means[others]
            # Kronecker symbol (-N/p) on the classes other than 0 and -1
            plus = [c for c in others if pow((-c) % p, (p - 1) // 2, p) == 1]
            minus = [c for c in others if c not in plus]
            mp = np.mean(means[plus]) if plus else float("nan")
            mm = np.mean(means[minus]) if minus else float("nan")
            sep = se * np.sqrt(1 / max(len(plus), 1) + 1 / max(len(minus), 1))
            line = (f"p={p:2d}: class 0 {means[0]:.4f}, class -1 {means[p-1]:.4f}; others min {om.min():.4f} max {om.max():.4f} "
                    f"(spread {om.max()-om.min():.4f} = {(om.max()-om.min())/se:.1f} s.e. of one class mean, s.e. {se:.4f}); "
                    f"top two = {{0,-1}}: {top2 == {0, p - 1}}")
            if plus and minus:
                line += f"; (-N/p)=+1 classes {mp:.4f} vs -1 classes {mm:.4f} (diff {mp-mm:+.4f}, {(mp-mm)/sep:+.1f} s.e.)"
            print(line)
            if p == 7:
                print("      p = 7 class means c = 0..6:", [round(float(x), 4) for x in means])
        if name.startswith("r_tri"):
            sieve = np.ones(X + 2, dtype=bool)
            sieve[:2] = False
            for i in range(2, int((X + 2) ** 0.5) + 1):
                if sieve[i]:
                    sieve[i * i:: i] = False
            m4 = win % 4 == 3
            rough = m4 & sieve[win] & sieve[(win + 1) // 4]
            m15 = (win * (win + 1)) % 15 == 0
            rough29 = np.ones(len(win), dtype=bool)
            for p in (3, 5, 7, 11, 13, 17, 19, 23, 29):
                rough29 &= (win % p != 0) & ((win + 1) % p != 0)
            for label, mask in (("N prime and (N+1)/4 prime", rough), ("N(N+1) free of odd primes <= 29", rough29),
                                ("15 | N(N+1)", m15)):
                sel = rw[mask]
                print(f"    {label:34s}: n = {mask.sum():7d}, mean {sel.mean():.4f} (ratio {sel.mean()/mean:.3f}), miss rate "
                      f"{(sel == 0).mean():.5f} (ratio {(sel == 0).mean()/(miss/len(win)):.3f})")
            heeg = {d: int(r[d]) for d in (3, 7, 11, 19, 43, 67, 163)}
            print("    r at the Heegner numbers:", heeg, "; r_pair there:", {d: int(r_pair[d]) for d in heeg})
            cand = [n for n in range(7, 3000, 4) if sieve[n] and sieve[(n + 1) // 4]]
            zero = [n for n in cand if r[n] == 0]
            print(f"    N = 4q - 1 < 3000 with N, q prime: {len(cand)} integers; r = 0 at {zero}")
            ok &= heeg == {3: 1, 7: 0, 11: 2, 19: 4, 43: 0, 67: 0, 163: 0}
            ok &= zero == [7, 43, 67, 163, 211, 547, 1171, 1723] and len(cand) == 30
    print("AUDIT F1:", "ALL CHECKS PASSED" if ok else "SOME CHECK FAILED", f"[{time.time()-t0:.0f}s]")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
