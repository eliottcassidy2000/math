# c3_spectrum_kpc1_verify.py — ADVERSARIAL VERIFICATION (session kind-pasteur-2026-06-10-S1)
# Independent computation of the c3 spectrum A(n) for n=3..60 (claims A2/A4 and the
# computational side of A3/A5).
#
# Fresh methods, different from worker's block-based bitset DP:
#  (1) BRUTE: plain backtracking over nondecreasing score sequences with the Landau
#      prefix condition checked at EVERY prefix (no concavity/block reduction at all),
#      for n=3..14; also counts Landau sequences (cross-check vs A000571).
#  (2) DP: one-vertex-at-a-time layered DP keyed (k, S) -> bitset of achievable
#      f = sum C(s_i,2); processes score values v=0..n-1 in increasing order, and within
#      a value sweeps k upward so arbitrary multiplicities accumulate; Landau prefix
#      bound checked at every single prefix, plus the exact prefix-average upper prune
#      S_k <= floor(k(n-1)/2) (the k smallest scores average at most the global mean).
# c3 = C(n,3) - f. Gap-freeness of c3-spectrum <=> f-set is a contiguous bit run.
import time
from math import comb

def brute(n):
    """Backtracking over Landau sequences; returns (count, fset)."""
    total = comb(n, 2)
    fset = set()
    count = 0

    def rec(k, S, last, f):
        nonlocal count
        if k == n:
            if S == total:
                count += 1
                fset.add(f)
            return
        rem = n - k - 1
        for v in range(last, n):
            S2 = S + v
            if S2 < comb(k + 1, 2):
                continue
            if S2 + rem * v > total:
                break          # later values only larger
            if S2 + rem * (n - 1) < total:
                continue       # cannot reach total
            rec(k + 1, S2, v, f + comb(v, 2))

    rec(0, 0, 0, 0)
    return count, fset

def dp_fbits(n):
    """Layered bitset DP; returns int whose bit f is set iff f achievable."""
    total = comb(n, 2)
    layers = [dict() for _ in range(n + 1)]
    layers[0][0] = 1
    lo = [comb(k, 2) for k in range(n + 1)]
    hi = [(k * (n - 1)) // 2 for k in range(n + 1)]
    for v in range(n):
        cv = comb(v, 2)
        for k in range(n):
            layer = layers[k]
            if not layer:
                continue
            nxt = layers[k + 1]
            l2, h2 = lo[k + 1], hi[k + 1]
            for S, bits in list(layer.items()):
                S2 = S + v
                if S2 < l2 or S2 > h2:
                    continue
                nb = bits << cv
                if S2 in nxt:
                    nxt[S2] |= nb
                else:
                    nxt[S2] = nb
    return layers[n].get(total, 0)

def M_formula(n):
    return (n**3 - n) // 24 if n % 2 == 1 else (n**3 - 4 * n) // 24

def m_nearreg(n):
    h = n // 2
    if n % 2 == 1:
        return n * comb(h, 2)
    return h * comb(h - 1, 2) + h * comb(h, 2)

# A000571 (number of score sequences of n-tournament), n=0..16, from OEIS:
A000571 = [1, 1, 1, 2, 4, 9, 22, 59, 167, 490, 1486, 4639, 14805, 48107,
           158808, 531469, 1799659]

t0 = time.time()
print("=== Part 1: brute backtracking n=3..14 (every-prefix Landau check) ===")
brute_sets = {}
for n in range(3, 15):
    cnt, fset = brute(n)
    brute_sets[n] = fset
    c3set = sorted(comb(n, 3) - f for f in fset)
    gaps = [x for x in range(c3set[0], c3set[-1] + 1) if x not in set(c3set)]
    ok571 = (cnt == A000571[n])
    print(f"n={n}: #Landau={cnt} (A000571 match={ok571}), c3 range "
          f"[{c3set[0]},{c3set[-1]}], |A|={len(c3set)}, gaps={gaps if gaps else 'NONE'}, "
          f"max==formula:{c3set[-1] == M_formula(n)}, minf==nearreg:"
          f"{min(fset) == m_nearreg(n)}")

print()
print("=== Part 2: independent layered DP n=3..60 ===")
all_ok = True
for n in range(3, 61):
    bits = dp_fbits(n)
    if bits == 0:
        print(f"n={n}: EMPTY f-set -- FAIL")
        all_ok = False
        continue
    fmax = bits.bit_length() - 1
    fmin = (bits & -bits).bit_length() - 1
    size = bin(bits).count("1")
    contiguous = (bits >> fmin) == (1 << (fmax - fmin + 1)) - 1
    c3max = comb(n, 3) - fmin
    c3min = comb(n, 3) - fmax
    okM = (c3max == M_formula(n))
    okmin = (c3min == 0)
    oknr = (fmin == m_nearreg(n))
    oksize = (size == M_formula(n) + 1)
    line_ok = contiguous and okM and okmin and oknr and oksize
    extra = ""
    if n in brute_sets:
        match = set(range(fmin, fmax + 1)) if contiguous else None
        # exact set equality vs brute:
        dpset = {f for f in range(fmin, fmax + 1) if (bits >> f) & 1}
        extra = f", DPset==brute:{dpset == brute_sets[n]}"
        line_ok = line_ok and (dpset == brute_sets[n])
    if not line_ok:
        all_ok = False
    print(f"n={n}: maxc3={c3max} (formula:{okM}), minc3={c3min} (==0:{okmin}), "
          f"|A|={size} (==M+1:{oksize}), gapfree={contiguous}, "
          f"minf==nearreg:{oknr}{extra}")

print()
print("max-c3 terms n=3..18 (for OEIS value search):",
      [M_formula(n) for n in range(3, 19)])
print("spectrum-size terms n=3..15:", [M_formula(n) + 1 for n in range(3, 16)])
print(f"elapsed {time.time() - t0:.1f} s")
print("VERDICT(A2 spectrum n=3..60 gap-free + formulas):",
      "PASS" if all_ok else "FAIL")
