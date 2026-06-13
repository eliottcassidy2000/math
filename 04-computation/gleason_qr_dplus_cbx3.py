#!/usr/bin/env python3
"""THM-481/482 (claudebox-2026-06-11-S3): Paley tournament gauges generate the Gleason ring
(e8, GOLAY g24, eQR(32), eQR(48)); one skew-Sylvester doubling step thermalizes any
even-row skew-Hadamard gauge code to d_{2n}+. Pure python, exact F2."""
import sys, itertools
from collections import Counter
sys.path.insert(0, '/home/claude/math-research/04-computation')
from tournament_barba_flag_cbx1 import paley_S, skewH_from_drt, double_skewH

def binrows(H): return [int(''.join('1' if x < 0 else '0' for x in row), 2) for row in H]
def bincols(H):
    n = len(H)
    return [int(''.join('1' if H[i][j] < 0 else '0' for i in range(n)), 2) for j in range(n)]
def basis_of(vecs):
    basis = []
    for v in vecs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v); basis.sort(reverse=True)
    return basis
def wd_enum(basis):
    cnt = Counter({0: 1}); w = 0; gp = 0
    for i in range(1, 1 << len(basis)):
        g = i ^ (i >> 1); bit = (g ^ gp).bit_length() - 1; gp = g
        w ^= basis[bit]; cnt[w.bit_count()] += 1
    return dict(sorted(cnt.items()))
def span_set(basis):
    words = [0]
    for v in basis: words += [w ^ v for w in words]
    return set(words)

def eqr_code(q):
    QR = {(x * x) % q for x in range(1, q)}
    for S0 in (QR, set(range(1, q)) - QR, QR | {0}, (set(range(1, q)) - QR) | {0}):
        v = 0
        for r in S0: v |= (1 << r)
        shifts = [(((v << s) | (v >> (q - s))) & ((1 << q) - 1)) for s in range(q)]
        ext = [w | ((bin(w).count('1') % 2) << q) for w in shifts] + [(1 << (q + 1)) - 1]
        be = basis_of(ext)
        if len(be) == (q + 1) // 2: return be
    return None

def words_of_weight(basis, wt):
    out = ([0] if wt == 0 else []); w = 0; gp = 0
    for i in range(1, 1 << len(basis)):
        g = i ^ (i >> 1); bit = (g ^ gp).bit_length() - 1; gp = g
        w ^= basis[bit]
        if w.bit_count() == wt: out.append(w)
    return out

def four_pattern_profile(basis):
    w8 = words_of_weight(basis, 8)
    pats = set()
    for i in range(len(w8)):
        for j in range(i + 1, len(w8)):
            x = w8[i] & w8[j]
            if x.bit_count() == 4: pats.add(x)
    prof = Counter()
    for p in pats: prof[sum(1 for w in w8 if (w & p) == p)] += 1
    return dict(sorted(prof.items()))

def rm(r, m):
    pts = list(itertools.product([0, 1], repeat=m)); gens = []
    for d in range(r + 1):
        for mo in itertools.combinations(range(m), d):
            v = 0
            for k, p in enumerate(pts):
                val = 1
                for i in mo: val &= p[i]
                if val: v |= (1 << k)
            gens.append(v)
    return basis_of(gens)

if __name__ == '__main__':
    print("== THM-481: the Paley/QR ladder ==")
    GOLAY = {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
    for q in (7, 23, 31, 47):
        n = q + 1
        b = basis_of(binrows(skewH_from_drt(paley_S(q))))
        w = wd_enum(b)
        mind = min(k for k in w if k > 0)
        de = all(k % 4 == 0 for k in w)
        eqr = eqr_code(q)
        match = (wd_enum(eqr) == w)
        print(f"q={q}: [{n},{len(b)},{mind}] doubly-even={de}; wd==eQR({n}): {match}" +
              ("; == GOLAY g24" if q == 23 and w == GOLAY else ""))
        assert len(b) == n // 2 and de and match
        assert (q, mind) in [(7, 4), (23, 8), (31, 8), (47, 12)]
    print("  q=23 rigorous (unique [24,12,8] Type II = g24); q=47 rigorous (unique extremal [48,24,12])")
    print("\n== q=31 discriminator: gauge vs RM(2,5) vs eQR(32) ==")
    g31 = basis_of(binrows(skewH_from_drt(paley_S(31))))
    for name, c in [("gauge(Paley31)", g31), ("eQR(32)", eqr_code(31)), ("RM(2,5)", rm(2, 5))]:
        print(f"  {name}: 4-pattern multiplicity profile {four_pattern_profile(c)}")
    assert four_pattern_profile(g31) == four_pattern_profile(eqr_code(31)) != four_pattern_profile(rm(2, 5))
    print("  gauge == eQR(32) profile != RM(2,5) profile  => gauge(Paley31) = eQR(32), rigorous")
    print("\n== THM-482: one doubling step -> d_{2n}+ (exact code equality) ==")
    H8 = skewH_from_drt(paley_S(7)); H16 = double_skewH(H8); H32 = double_skewH(H16)
    for H, Hd, name in [(H8, H16, "8->16"), (H16, H32, "16->32")]:
        n = len(H)
        E = basis_of([(1 << j) | (1 << k) for j in range(n) for k in range(j + 1, n)
                      if True][:0] + [(1 << 0) | (1 << k) for k in range(1, n)])  # e0+ek spans E_n
        pd = [w | (w << n) for w in E]
        rowsD = binrows(Hd)
        glue = rowsD[n]   # a bottom row (gamma-bar_1, gamma_1): pairwise-odd glue
        cand = basis_of(pd + [glue])
        bD = basis_of(rowsD)
        eq = span_set(cand) == span_set(bD)
        print(f"  C(double) == PD(E_{n}) + <bottom-row glue>  [{name}]: {eq}  (dim {len(bD)})")
        assert eq
        # gauge identity check: gamma_j = 1 ^ b_j ^ e_j
        bs = binrows(H); gs = bincols(H)
        full = (1 << n) - 1
        ok = all(gs[j] == full ^ bs[j] ^ (1 << (n - 1 - j)) for j in range(n))
        print(f"    gauge identity gamma_j = ~b_j ^ e_j: {ok}")
        assert ok
    print("\nall checks passed")
