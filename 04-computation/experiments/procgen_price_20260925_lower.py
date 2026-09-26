#!/usr/bin/env python3
"""procgen_price_20260925 -- part 3: the provable lower bound for the price, and where the local rescues fail.

(a) Weighted-multiplicity bound (PROVED in the note, section 3).  Every member of P_L satisfies
        delta >= rho_L / W*_L,   W*_L = max_v  sum_{k < L} sum_{n : T^k(n) = v} v/n,
    and the backward tree of the shortcut map T gives W*_L <= (1 + o(1)) sum_{k<L} g(k) with
        g(0) = 1,  g(1) = 2,  g(k) = 1.5 g(k-1) + 0.25 g(k-2)   (growth lambda = (3 + sqrt 13)/4 = 1.65139).
    We print g, the bound, and the resulting exponent (1 - h) + log2(lambda) = 0.7737.
    Cross-check: a brute-force maximum of the weighted multiplicity over all v <= V (exact backward trees).
(b) The '111' criterion.  A low rescue point with a 2-step partner is an odd up-image y = 3 mod 4, i.e. a factor
    111 of the parity word.  Count the Bad_L classes whose word avoids 111 altogether (the rescue-free part of the
    undecided set); its growth rate is the dimension of the bad no-111 subshift (containing the 2-adic cycle -5 -> -7 -> -10).
(c) The 5x+1 (DRIFT) undecided density beta_L: its limit is positive.
Runtime < 1 minute; memory < 100 MB.
"""
import math

LOG32 = math.log(2) / math.log(3)
h = -(LOG32 * math.log2(LOG32) + (1 - LOG32) * math.log2(1 - LOG32))
lam = (3 + math.sqrt(13)) / 4

def thr(k):
    a = 0
    while 3 ** a <= 2 ** k:
        a += 1
    return a

def bad_count(L):
    T = [thr(k) for k in range(L + 1)]
    dist = {0: 1}
    for k in range(1, L + 1):
        nd = {}
        for a, c in dist.items():
            for up in (0, 1):
                b = a + up
                if b >= T[k]:
                    nd[b] = nd.get(b, 0) + c
        dist = nd
    return sum(dist.values())

def g_seq(L):
    g = [1.0, 2.0]
    while len(g) < L:
        g.append(1.5 * g[-1] + 0.25 * g[-2])
    return g[:L]

def brute_W(V, L):
    """max over v <= V of sum over backward paths (length < L) of 3^a/2^k (the leading term of v/n), exact trees."""
    best = (0.0, None)
    for v in range(1, V + 1):
        tot = 0.0
        stack = [(v, 0, 0)]
        while stack:
            x, k, a = stack.pop()
            tot += 3.0 ** a / 2.0 ** k
            if k + 1 < L:
                stack.append((2 * x, k + 1, a))
                if x % 3 == 2:
                    stack.append(((2 * x - 1) // 3, k + 1, a + 1))
        if tot > best[0]:
            best = (tot, v)
    return best

def no111_bad(L):
    """number of words of length L, prefix-bad (a_k >= thr(k)), avoiding the factor 111."""
    T = [thr(k) for k in range(L + 1)]
    dist = {(0, 0): 1}           # (a, trailing ones run length capped at 2)
    for k in range(1, L + 1):
        nd = {}
        for (a, run), c in dist.items():
            for up in (0, 1):
                r = run + 1 if up else 0
                if r >= 3:
                    continue
                b = a + up
                if b >= T[k]:
                    key = (b, r)
                    nd[key] = nd.get(key, 0) + c
        dist = nd
    return sum(dist.values())

def drift_beta(L):
    p = math.log(2) / math.log(5)
    dist = {0: 1.0}
    for k in range(1, L + 1):
        t = math.floor(k * p) + 1
        nd = {}
        for a, c in dist.items():
            for up in (0, 1):
                b = a + up
                if b >= t:
                    nd[b] = nd.get(b, 0) + 0.5 * c
        dist = nd
    return sum(dist.values())

if __name__ == "__main__":
    print("=" * 100)
    print("PART 3. Provable lower bound for the price; the 111 criterion; the DRIFT undecided density")
    print("=" * 100)
    print(f"  lambda = (3 + sqrt 13)/4 = {lam:.6f};  log2(lambda) = {math.log2(lam):.6f};  1 - h = {1 - h:.6f}")
    print(f"  lower-bound exponent (1 - h) + log2(lambda) = {1 - h + math.log2(lam):.6f}   (conjectured true exponent: 1 - h)")
    print("   L    rho_L        W*_L bound (sum g)    delta_L >= rho_L / W*_L     -log2(bound)/L")
    for L in (2, 4, 5, 8, 12, 16, 20, 24, 28, 32, 36, 40):
        rho = bad_count(L) / 2 ** L
        W = sum(g_seq(L))
        lb = rho / W
        print(f"  {L:3d}  {rho:.6f}   {W:14.2f}         {lb:.3e}                {-math.log2(lb) / L:.4f}")
    print("  brute-force check of the tree bound (max over v <= V of sum of 3^a/2^k over backward paths of length < L):")
    for L, V in ((6, 20000), (8, 20000), (10, 20000), (12, 8000)):
        bw, bv = brute_W(V, L)
        print(f"    L = {L:2d}, v <= {V}: max = {bw:.3f} at v = {bv};  bound sum g = {sum(g_seq(L)):.3f}")
    print()
    print("  (b) bad classes avoiding the factor 111 (no low rescue point with a 2-step partner):")
    print("   L    |Bad_L|         no-111 part     fraction     growth (log2 ratio per step, window 10)")
    prev = {}
    for L in (8, 12, 16, 20, 24, 30, 40, 60, 80, 100):
        b = bad_count(L) if L <= 60 else None
        n1 = no111_bad(L)
        g = ""
        if L >= 20:
            g = f"{(math.log2(n1) - math.log2(no111_bad(L - 10))) / 10:.4f}"
        bs = f"{b:>14d}" if b is not None else " " * 14
        fr = f"{n1 / b:.2e}" if b else "   -    "
        print(f"  {L:3d}  {bs}   {n1:>14d}   {fr}     {g}")
    print("   => bad words with no factor 111 at all grow like 2^(0.3 L): on them no F-rescue (low up-image = 3 mod 4 with a")
    print("      2-step partner) ever exists; they contain the 2-adic cycle (110)^inf = -5 -> -7 -> -10.  The construction")
    print("      needs its B-rescue (partner not 2-step) only when also A is blocked, which forces n = -5 mod 2^9 (note, lemma C6).")
    print()
    print("  (c) DRIFT: undecided density beta_L of 5x+1 (words with a_k > k log_5 2 for all k <= L):")
    for L in (4, 8, 12, 16, 20, 24, 30, 40, 60, 100, 200, 400):
        print(f"    L = {L:3d}: beta_L = {drift_beta(L):.5f}")
    print("   => beta_L decreases to a positive limit (about 0.17): for 5x+1 the undecided set does not thin out.")
