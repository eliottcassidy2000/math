# Hamiltonian paths and cycles of Paley tournaments P_p (p = 3 mod 4): closing / non-closing paths, orbits under Aut.
import sys
def paley(p):
    qr = {(x * x) % p for x in range(1, p)}
    out = [0] * p
    for a in range(p):
        for b in range(p):
            if a != b and ((b - a) % p) in qr: out[a] |= 1 << b
    return out, qr
def count(out, n):
    # paths starting at vertex 0 (vertex-transitive): dp over subsets containing 0
    full = (1 << n) - 1
    dp = {}
    dp[(1, 0)] = 1
    layers = [dict() for _ in range(n + 1)]
    layers[1][(1, 0)] = 1
    for size in range(1, n):
        for (S, v), c in layers[size].items():
            m = out[v] & ~S
            while m:
                w = (m & -m).bit_length() - 1
                key = (S | (1 << w), w)
                layers[size + 1][key] = layers[size + 1].get(key, 0) + c
                m &= m - 1
        layers[size] = None
    paths0 = 0; closing0 = 0
    for (S, v), c in layers[n].items():
        paths0 += c
        if (out[v] >> 0) & 1: closing0 += c
    return paths0, closing0
for p in (3, 7, 11, 19):
    out, qr = paley(p)
    p0, c0 = count(out, p)
    H = p * p0; closing = p * c0; cyc = c0  # Hamiltonian cycles = closing paths from vertex 0
    aut = p * (p - 1) // 2
    print("P_%d: H = %d = %d * %d (|Aut| = p(p-1)/2); Hamiltonian cycles c_p = %d; closing paths %d; non-closing paths %d = %d * |Aut|" % (
        p, H, aut, H // aut, cyc, closing, H - closing, (H - closing) // aut), "| H/|Aut| =", H // aut, "| cycles mod |Aut|:", cyc % aut)
