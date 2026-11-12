# Verifier (lane architecture): independent exact check of P4 (one-exit extension in source coordinates) and of the
# one-vertex closing law used by route R4:  H(B => TT_t => v, v -> y, every other vertex of B -> v) = H(B) + 2^t U_B(y),
# t >= 1.  Independent code: bitmask dynamic programming (the lane used permutations).
# Exhaustive for rooted labelled tournaments of orders 1..5 (5,405 cases) and a random sample at orders 6, 7.
import random, sys

def tables(n, out):
    # out[v] = bitmask of out-neighbours of v.  st[mask][v] = number of Hamiltonian paths of mask starting at v,
    # en[mask][v] = number ending at v, hp[mask] = number of Hamiltonian paths of mask (hp[0] = 1).
    full = 1 << n
    en = [[0] * n for _ in range(full)]
    for v in range(n):
        en[1 << v][v] = 1
    for mask in range(1, full):
        row = en[mask]
        for v in range(n):
            c = row[v]
            if c == 0: continue
            nxt = out[v] & ~mask
            while nxt:
                b = nxt & -nxt
                w = b.bit_length() - 1
                en[mask | b][w] += c
                nxt ^= b
    inn = [0] * n
    for v in range(n):
        for w in range(n):
            if out[v] >> w & 1:
                inn[w] |= 1 << v
    st = [[0] * n for _ in range(full)]
    for v in range(n):
        st[1 << v][v] = 1
    for mask in range(1, full):
        row = st[mask]
        for v in range(n):
            c = row[v]
            if c == 0: continue
            nxt = inn[v] & ~mask
            while nxt:
                b = nxt & -nxt
                w = b.bit_length() - 1
                st[mask | b][w] += c
                nxt ^= b
    hp = [sum(en[mask]) for mask in range(full)]
    hp[0] = 1
    return st, en, hp

def invariants(n, out, x):
    st, en, hp = tables(n, out)
    full = (1 << n) - 1
    h = hp[full]
    e = en[full][x]
    u = 0
    for S in range(1, full + 1):
        if S >> x & 1:
            u += st[S][x] * hp[full ^ S]
    return h, e, u

def random_out(n, rng):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.getrandbits(1): out[i] |= 1 << j
            else: out[j] |= 1 << i
    return out

def out_from_bits(n, bits):
    out = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> k & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
            k += 1
    return out

def extend(n, out, x):
    # B = B' + y, y -> x, every other vertex -> y ; y has index n
    o2 = out[:] + [1 << x]
    for v in range(n):
        if v != x: o2[v] |= 1 << n
    return o2

def closing(n, out, y, t):
    # B => TT_t => v with v -> y and every other vertex of B -> v.  Vertices: B (0..n-1), TT (n..n+t-1), v = n+t.
    N = n + t + 1
    o = out[:] + [0] * (t + 1)
    tt = list(range(n, n + t)); v = n + t
    for b in range(n):
        for s in tt: o[b] |= 1 << s
        if b != y: o[b] |= 1 << v
    for a in range(t):
        for c in range(a + 1, t): o[tt[a]] |= 1 << tt[c]
        o[tt[a]] |= 1 << v
    o[v] |= 1 << y
    st, en, hp = tables(N, o)
    return hp[(1 << N) - 1]

bad = 0; cnt = 0; badc = 0; cntc = 0
for n in range(1, 6):
    m = n * (n - 1) // 2
    for bits in range(1 << m):
        out = out_from_bits(n, bits)
        for x in range(n):
            h, e, u = invariants(n, out, x)
            o2 = extend(n, out, x)
            HB, EB, UB = invariants(n + 1, o2, n)
            cnt += 1
            if not (HB == h - e + u and UB == h + u and EB == h - e and u >= h):
                bad += 1
                if bad < 5: print("FAIL P4", n, bits, x, (h, e, u), (HB, EB, UB))
            if n <= 4:
                for t in (1, 2, 3):
                    cntc += 1
                    if closing(n, out, x, t) != h + (1 << t) * u:
                        badc += 1
                        if badc < 5: print("FAIL closing", n, bits, x, t)
print("P4 exhaustive, rooted labelled tournaments of orders 1..5:", cnt, "failures:", bad)
print("one-vertex closing law H = H(B) + 2^t U_B(y), orders 1..4, t = 1, 2, 3:", cntc, "failures:", badc)
rng = random.Random(20261111)
cnt2 = 0; bad2 = 0; cntc2 = 0; badc2 = 0
for n, reps in ((6, 1500), (7, 600)):
    for _ in range(reps):
        out = random_out(n, rng)
        x = rng.randrange(n)
        h, e, u = invariants(n, out, x)
        o2 = extend(n, out, x)
        HB, EB, UB = invariants(n + 1, o2, n)
        cnt2 += 1
        if not (HB == h - e + u and UB == h + u and EB == h - e and u >= h):
            bad2 += 1
        if n == 6 and cntc2 < 300:
            t = rng.choice((1, 2))
            cntc2 += 1
            if closing(n, out, x, t) != h + (1 << t) * u:
                badc2 += 1
print("P4 random sample, orders 6 and 7:", cnt2, "failures:", bad2)
print("closing law random sample, order 6, t in {1, 2}:", cntc2, "failures:", badc2)
print("CHECK", "PASS" if bad == 0 and badc == 0 and bad2 == 0 and badc2 == 0 else "FAIL")
