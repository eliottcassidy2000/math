"""P106 generalized: every 'all-above' gadget G (G beats T-A, arbitrary arcs G<->A, hub stays A) acts on
(H, U, sigma, K) by   H' = g H + s2 U - geo sigma + s2p K,  U' = aall U + aA K,  sigma' = gst K,  K' = g K.
Continuation data = (g, s2, aall, geo); deposit column = (s2p, aA, gst).  Verify the formula on random T, then
search all G with |G| <= 5 for equal continuation data and different deposit columns."""
import itertools, random, sys
sys.path.insert(0, '.')
from w45_hk import ham_split, sub

def paths(k, out):            # all directed paths (as tuples) of a k-vertex tournament
    res = []
    def ext(p, used):
        res.append(tuple(p))
        for w in range(k):
            if not used >> w & 1 and out[p[-1]] >> w & 1:
                ext(p + [w], used | 1 << w)
    for v in range(k):
        ext([v], 1 << v)
    return res

def Hsub(k, out, rem):        # H of G restricted to vertex set rem (bitmask)
    keep = [v for v in range(k) if rem >> v & 1]
    return ham_split(*sub(out, keep))[0]

def data(k, out, X):          # X = bitmask of G-vertices that A beats
    full = (1 << k) - 1
    P = paths(k, out)
    g, st, en = ham_split(k, out)
    geo = sum(en[v] for v in range(k) if X >> v & 1)
    gst = sum(st[v] for v in range(k) if X >> v & 1)
    s2 = s2p = 0
    for gam in P:                         # split (alpha, gamma): gamma nonempty, alpha = path cover of rest
        if not X >> gam[0] & 1: continue
        rem = full
        for v in gam: rem &= ~(1 << v)
        if rem == 0:
            s2 += 1; s2p += 1; continue
        for al in P:
            m = 0
            for v in al: m |= 1 << v
            if m == rem and len(al) == bin(rem).count('1'):
                s2 += 1
                if not X >> al[-1] & 1: s2p += 1
    aall = 1 * g  # alpha empty -> H(G)
    aA = 0
    for al in P:
        m = full
        for v in al: m &= ~(1 << v)
        h = Hsub(k, out, m)
        aall += h
        if not X >> al[-1] & 1: aA += h
    return (g, s2, aall, geo), (s2p, aA, gst)

def regs(n, out, A):          # (H, U, sigma, K) of (T, A) by brute force
    H, st, en = ham_split(n, out)
    rest = [v for v in range(n) if v != A]
    K = ham_split(*sub(out, rest))[0]
    E = 0
    def ext(p, used):
        nonlocal E
        if p[-1] == A:
            E += ham_split(*sub(out, [v for v in range(n) if not used >> v & 1]))[0]
            return
        for w in range(n):
            if not used >> w & 1 and out[p[-1]] >> w & 1:
                ext(p + [w], used | 1 << w)
    for v in range(n):
        ext([v], 1 << v)
    return H, E - K, st[A], K

def randtour(n, rng):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < .5: out[i] |= 1 << j
            else: out[j] |= 1 << i
    return out

if __name__ == '__main__':
    rng = random.Random(7); bad = 0
    for trial in range(150):
        n = rng.randint(1, 6); k = rng.randint(1, 4)
        T = randtour(n, rng); G = randtour(k, rng); X = rng.randrange(1 << k); A = 0
        o = T[:] + [0] * k
        for a in range(k):
            o[n + a] = G[a] << n
            for i in range(n):
                if i != A: o[n + a] |= 1 << i
            if X >> a & 1: o[A] |= 1 << (n + a)
            else: o[n + a] |= 1 << A
        H, U, s, K = regs(n, T, A); H2, U2, s2_, K2 = regs(n + k, o, A)
        (g, s2, aall, geo), (s2p, aA, gst) = data(k, G, X)
        if (H2, U2, s2_, K2) != (g*H + s2*U - geo*s + s2p*K, aall*U + aA*K, gst*K, g*K): bad += 1
    print("gadget matrix formula: failures", bad, "of 150")
    for kmax in (3, 4, 5):
        groups = {}
        for k in range(1, kmax + 1):
            pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
            for bits in range(1 << len(pairs)):
                out = [0] * k
                for t, (i, j) in enumerate(pairs):
                    if bits >> t & 1: out[i] |= 1 << j
                    else: out[j] |= 1 << i
                for X in range(1 << k):
                    c, d = data(k, out, X)
                    groups.setdefault(c, set()).add(d)
        multi = {c: ds for c, ds in groups.items() if len(ds) > 1}
        g1 = [c for c in multi if c[0] == 1]
        print(f"|G|<={kmax}: {len(groups)} continuation classes; {len(multi)} carry 2+ deposit columns; with g=1: {len(g1)}")
        for c in sorted(multi)[:4]:
            print("   continuation (g,s2,aall,geo) =", c, "deposits (s2p,aA,gst):", sorted(multi[c])[:6])
