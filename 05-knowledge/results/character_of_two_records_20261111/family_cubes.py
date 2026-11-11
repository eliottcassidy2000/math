#!/usr/bin/env python3
"""Pentagon cube sums G_35(R) = sum_X (-1)^((c_3 + c_5)(R + v_X)) of Paley and Cayley-Dickson tournaments, by the
kernel formula of THM-5132 (3), cross-checked (a) by the Gauss sum of the quadratic form x^T (A + A^3)(1 + x) computed
by symplectic elimination, and (b) for m <= 11 by direct counting of 3- and 5-cycles in all 2^m extensions.

Kernel formula (K_0 = ker A^T, K_1 = {u : A^T u = 1} over F_2; n_ij = number of elements of K_i of weight j mod 4):
  m = 0 (4): 2^(m/2) |K_0|;   m = 1 (4): 2^((m-1)/2) (|K_0| + n13 - n10);
  m = 2 (4): 2^(m/2) (n00 + n13);   m = 3 (4): 2^((m+1)/2) (n00 + n10).
Usage: family_cubes.py"""
import sys, itertools
from paley_tools import paley_out
from cayley_dickson import tournament as cd_tournament

def transpose_cols(out, m):
    inn = [0] * m
    for a in range(m):
        for b in range(m):
            if (out[a] >> b) & 1: inn[b] |= 1 << a
    return inn
def nullspace(rows, ncols):
    """basis of {v : every row . v = 0} over F_2; rows are bitmasks over ncols columns."""
    red = []
    for r in rows:
        for p, rr in red:
            if (r >> p) & 1: r ^= rr
        if r:
            p = r.bit_length() - 1
            red = [(pp, (rr ^ r) if (rr >> p) & 1 else rr) for pp, rr in red] + [(p, r)]
    piv = dict(red); basis = []
    for f in range(ncols):
        if f in piv: continue
        v = 1 << f
        for p, rr in piv.items():
            if (rr >> f) & 1: v |= 1 << p
        basis.append(v)
    return basis
def kernel_data(out, m):
    """K = {(u, tau) : A^T u = tau 1}; returns dim K_0, [K_1 nonempty], counts n[tau][w mod 4] (exact, by enumeration
    if the dimension is small, by the doubly-even criterion otherwise)."""
    inn = transpose_cols(out, m)                      # (A^T u)_b = <inn[b], u>
    rows = [inn[b] | (1 << m) for b in range(m)]      # equation <inn[b], u> + tau = 0
    basis = nullspace(rows, m + 1)
    full = (1 << m) - 1
    d = len(basis)
    n = [[0] * 4, [0] * 4]
    if d <= 22:
        v = 0
        for g in range(1 << d):
            if g:
                v ^= basis[(g & -g).bit_length() - 1]
            n[(v >> m) & 1][bin(v & full).count('1') & 3] += 1
        return d, n, "enumerated"
    # large kernel: use weights of a basis of K_0 and one element of K_1
    b0 = [b & full for b in basis if not (b >> m) & 1]
    b1 = [b for b in basis if (b >> m) & 1]
    # make a basis of K_0 plus one representative w of K_1
    if b1:
        w = b1[0]; b0 += [(b ^ w) & full for b in b1[1:]]; w &= full
    doubly_even = all(bin(b).count('1') % 4 == 0 for b in b0) and all(bin(a & b).count('1') % 2 == 0 for a in b0 for b in b0)
    if not doubly_even: raise RuntimeError("kernel too large and not doubly even")
    n[0][0] = 1 << len(b0)
    if b1:
        if not all(bin(w & b).count('1') % 2 == 0 for b in b0): raise RuntimeError("coset weights not constant mod 4")
        n[1][bin(w).count('1') & 3] = 1 << len(b0)
    return d, n, "doubly even basis"
def G_formula(m, n):
    k0 = n[0][0] + n[0][1]
    if m % 4 == 0: return (1 << (m // 2)) * k0
    if m % 4 == 1: return (1 << ((m - 1) // 2)) * (k0 + n[1][3] - n[1][0])
    if m % 4 == 2: return (1 << (m // 2)) * (n[0][0] + n[1][3])
    return (1 << ((m + 1) // 2)) * (n[0][0] + n[1][0])
def mat_rows_product(out, M, m):
    P = []
    for a in range(m):
        acc = 0; t = out[a]
        while t:
            b = (t & -t).bit_length() - 1; acc ^= M[b]; t &= t - 1
        P.append(acc)
    return P
def G_gauss(out, m):
    """(-1)^(c_3 + c_5)(R) times the Gauss sum of q(x) = sum_{a != d} W_ad x_a (1 + x_d), W = A + A^3 over F_2."""
    A2 = mat_rows_product(out, out, m); A3 = mat_rows_product(out, A2, m); A4 = mat_rows_product(out, A3, m); A5 = mat_rows_product(out, A4, m)
    eps = (sum((A3[a] >> a) & 1 for a in range(m)) + sum((A5[a] >> a) & 1 for a in range(m))) & 1
    W = [(out[a] ^ A3[a]) & ~(1 << a) for a in range(m)]
    Wt = [0] * m; l = 0
    for a in range(m):
        if bin(W[a]).count('1') & 1: l |= 1 << a
        t = W[a]
        while t:
            b = (t & -t).bit_length() - 1; Wt[b] |= 1 << a; t &= t - 1
    B = [(W[a] ^ Wt[a]) & ~(1 << a) for a in range(m)]
    active = (1 << m) - 1; c = 0; r = 0
    while True:
        a = -1; t = active
        while t:
            i = (t & -t).bit_length() - 1; t &= t - 1
            if B[i] & active: a = i; break
        if a < 0: break
        d = ((B[a] & active) & -(B[a] & active)).bit_length() - 1
        active &= ~((1 << a) | (1 << d)); Aa = B[a] & active; Bd = B[d] & active; la = (l >> a) & 1; ld = (l >> d) & 1
        t = Aa
        while t:
            i = (t & -t).bit_length() - 1; B[i] ^= Bd; t &= t - 1
        t = Bd
        while t:
            i = (t & -t).bit_length() - 1; B[i] ^= Aa; t &= t - 1
        l ^= (Aa & Bd)
        if la: l ^= Bd
        if ld: l ^= Aa
        c ^= (la & ld); r += 2
    if l & active: return 0, r
    val = 1 << (m - r // 2)
    if c ^ eps: val = -val
    return val, r
def G_direct(out, m):
    tot = 0
    def c35(o, n):
        c3 = n * (n - 1) * (n - 2) // 6 - sum(bin(x).count('1') * (bin(x).count('1') - 1) // 2 for x in o)
        c5 = 0
        for v in range(n):
            hi = ((1 << n) - 1) & ~((1 << (v + 1)) - 1)
            a = o[v] & hi
            while a:
                x = (a & -a).bit_length() - 1; a &= a - 1; b = o[x] & hi & ~(1 << x)
                while b:
                    y = (b & -b).bit_length() - 1; b &= b - 1; cc = o[y] & hi & ~((1 << x) | (1 << y))
                    while cc:
                        z = (cc & -cc).bit_length() - 1; cc &= cc - 1; dd = o[z] & hi & ~((1 << x) | (1 << y) | (1 << z))
                        while dd:
                            w = (dd & -dd).bit_length() - 1; dd &= dd - 1
                            if (o[w] >> v) & 1: c5 += 1
        return c3 + c5
    for X in range(1 << m):
        o = [out[a] | ((0 if (X >> a) & 1 else 1) << m) for a in range(m)] + [X]
        tot += -1 if c35(o, m + 1) & 1 else 1
    return tot
def report(name, out, m):
    d, n, how = kernel_data(out, m)
    Gf = G_formula(m, n); Gg, r = G_gauss(out, m)
    line = "%-8s m = %3d: dim ker A^T = %2d (%s); K_1 %s; weights mod 4: K_0 %s, K_1 %s; G_35 by the kernel formula = 2^%s%s; Gauss sum by elimination = %s (rank %d)" % (
        name, m, d - (1 if sum(n[1]) else 0), how, "non-empty" if sum(n[1]) else "empty", n[0], n[1],
        (Gf.bit_length() - 1) if Gf > 0 and Gf & (Gf - 1) == 0 else "?", "" if Gf > 0 and Gf & (Gf - 1) == 0 else " (value %d)" % Gf,
        "equal" if Gg == Gf else "DIFFERENT (%d)" % Gg, r)
    if m <= 11:
        Gd = G_direct(out, m); line += "; direct count over 2^m extensions = %s" % ("equal" if Gd == Gf else "DIFFERENT (%d)" % Gd)
    print(line); sys.stdout.flush()
    return Gf
if __name__ == "__main__":
    print("Paley tournaments P_q (R4: 2^q for q = 7 (mod 8), 2^((q+1)/2) for q = 3 (mod 8))")
    for q in (3, 7, 11, 19, 23, 27, 31, 43, 47):
        out = paley_out(q); G = report("P_%d" % q, out, q)
        pred = (1 << q) if q % 8 == 7 else (1 << ((q + 1) // 2))
        print("   prediction %s" % ("HOLDS" if G == pred else "FAILS"))
    print("Cayley-Dickson sign tournaments CD_k (R5: 2^(2^k - 1) for k = 3, 4, 5; k = 2, 6, 7 exploratory)")
    for k in (2, 3, 4, 5, 6, 7):
        out, g = cd_tournament(k); m = (1 << k) - 1
        G = report("CD_%d" % k, out, m)
        print("   G_35 = 2^m: %s" % ("yes" if G == (1 << m) else "no"))
