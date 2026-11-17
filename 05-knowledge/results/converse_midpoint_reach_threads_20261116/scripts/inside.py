#!/usr/bin/env python3
"""inside.py (Q-4): for each strong digit-0 pair (S, z) of orders 4-8 whose base S - z has a strong terminal component,
does an induced subtournament T of S with a port t and k >= 1 satisfy H(S) = h(T) + 2^k B_t(T)?  Also whether such a T
avoids z, and whether a one-vertex deletion S - u works."""
exec(open('nearjoin_check.py').read().split('# (1)')[0])
def rev(T): n = len(T); return [[T[j][i] for j in range(n)] for i in range(n)]
def Bvec(T):
    n = len(T); N = 1 << n; dp, H = ham_all(T); dpr, _ = ham_all(rev(T))
    return H[N-1], [sum(dpr[m][t] * H[(N-1) & ~m] for m in range(1, N) if (m >> t) & 1) for t in range(n)]
def term_comp(B):
    n = len(B)
    def reach(s):
        r = {s}; f = [s]
        while f:
            a = f.pop()
            for b in range(n):
                if B[a][b] and b not in r: r.add(b); f.append(b)
        return r
    for v in range(n):
        c = scc_of(B, v)
        if reach(v) == c: return sorted(c)
cases = 0; found = 0; found_noz = 0; found_del = 0; details = []
for n in range(4, 9):
    for S in classes(n):
        if not strong(S): continue
        H, B = Bvec(S); d0 = [z for z in range(n) if B[z] == H]
        for z in d0:
            rest = [v for v in range(n) if v != z]; comp = term_comp(sub(S, rest))
            if len(comp) == 1: continue
            cases += 1; ok = False; ok_noz = False; ok_del = False; wit = None
            for W in range(1, 1 << n):
                verts = [v for v in range(n) if (W >> v) & 1]
                if len(verts) == n: continue
                T = sub(S, verts); hT, BT = Bvec(T)
                for t, b in enumerate(BT):
                    k = 1
                    while hT + (b << k) <= H:
                        if hT + (b << k) == H:
                            ok = True
                            if z not in verts: ok_noz = True
                            if len(verts) == n - 1: ok_del = True
                            if wit is None: wit = (len(verts), hT, b, k, z in verts)
                        k += 1
            found += ok; found_noz += ok_noz; found_del += ok_del
            if not ok: details.append((H, n))
print(f"strong-D cases (orders 4-8): {cases}; with a near-join base induced inside S: {found}; with one avoiding z: {found_noz}; with a one-vertex deletion S - u as base: {found_del}")
print("cases with no inside base (H, order):", details)
