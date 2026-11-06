# Independent audit of THM-4887 (opus-2026-11-06-S83 audit), part 1: the two-block identity.
# Written from scratch; does not import the session's scripts.
#  (i)  g(v) = number of Hamiltonian paths of the one-apex tournament of the binary word of v
#       (brute-force subset DP on the actual tournament, all v < 2^11, leading zeros too);
#  (ii) g(h 2^m + l) = g(h) + h C_m(l) + gamma_h(l), C_m(l) = sum_{i<m} (1 - bit_i(l)) 2^(m-1-i),
#       gamma_h(l) = g(l) - 1 + (h mod 2)[l < 2^(m-1)]: exhaustive for small (h, m, l), random big-int cases;
#  (iii) 1 <= g(l) <= l + 1 and 0 <= gamma_h(l) <= 2^m - 1 (sharper than the stated <= 2^m);
#  (iv) the identity is exactly Session 193 eq. (1) H(uv) = H(u) + H(v) - 1 + A(u) C(v) + p(u) eps(v), checked
#       on words by brute force (u = word of h, v = m-letter word of l with leading zeros).
import random, sys
from functools import lru_cache

def g_rec(v):
    # g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v mod 2); iterative top-down on the bits of v
    bits = bin(v)[2:] if v > 0 else ""
    G, A = 1, 0          # (H, A) of the empty word; appending bits left to right
    for ch in bits:
        if ch == "1":
            A = 2 * A + 1
        else:
            G = G + A + (A & 1)
            A = 2 * A
    return G

def ham_paths(word):
    # one-apex tournament: chain c_1 -> c_2 -> ... -> c_L (i -> j for i < j), apex x beats c_i iff word[i] == '1'
    L = len(word)
    N = L + 1            # vertices 0..L-1 chain, L apex
    adj = [[False] * N for _ in range(N)]
    for i in range(L):
        for j in range(i + 1, L):
            adj[i][j] = True
    for i in range(L):
        if word[i] == "1":
            adj[L][i] = True
        else:
            adj[i][L] = True
    # dp[mask][last] = number of directed paths visiting mask ending at last
    full = (1 << N) - 1
    dp = [[0] * N for _ in range(1 << N)]
    for s in range(N):
        dp[1 << s][s] = 1
    for mask in range(1 << N):
        row = dp[mask]
        for last in range(N):
            c = row[last]
            if not c:
                continue
            for nxt in range(N):
                if not (mask >> nxt) & 1 and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += c
    return sum(dp[full])

def Cm(l, m):
    return sum((1 - ((l >> i) & 1)) << (m - 1 - i) for i in range(m))

# (i) brute force
bad = 0; tested = 0
for L in range(0, 11):
    for v in range(1 << L):
        w = format(v, "0%db" % L) if L > 0 else ""
        H = ham_paths(w) if L > 0 else 1
        if H != g_rec(v):
            bad += 1
        tested += 1
print(f"(i) g(v) = H(one-apex tournament of the L-letter word of v), all words of length <= 10 "
      f"(incl. leading zeros): {tested} words, {bad} mismatches")

# (ii) exhaustive small + random big
bad = 0; cnt = 0; gmax_ratio = 0.0; gam_min = 10**9; gam_max_minus_C = -10**9
for m in range(1, 9):
    for h in range(0, 1 << 9):
        for l in range(1 << m):
            lhs = g_rec((h << m) + l)
            gam = g_rec(l) - 1 + (h & 1) * (1 if l < (1 << (m - 1)) else 0)
            rhs = g_rec(h) + h * Cm(l, m) + gam
            cnt += 1
            if lhs != rhs:
                bad += 1
            gam_min = min(gam_min, gam); gam_max_minus_C = max(gam_max_minus_C, gam - (1 << m))
print(f"(ii) exhaustive m <= 8, 0 <= h < 2^9, all l: {cnt} cases, {bad} violations; "
      f"min gamma = {gam_min}, max (gamma - 2^m) = {gam_max_minus_C}")
random.seed(20261106)
bad = 0
for _ in range(100000):
    m = random.randint(1, 64)
    h = random.getrandbits(random.randint(1, 200)) | 1 if random.random() < 0.5 else random.getrandbits(random.randint(1, 200))
    l = random.getrandbits(m)
    lhs = g_rec((h << m) + l)
    rhs = g_rec(h) + h * Cm(l, m) + g_rec(l) - 1 + (h & 1) * (1 if l < (1 << (m - 1)) else 0)
    bad += lhs != rhs
print(f"(ii) random big-int cases (m <= 64, h < 2^200): 100000 cases, {bad} violations")

# (iii) bounds
worst = 0
for v in range(1 << 20):
    gv = g_rec(v)
    assert 1 <= gv <= v + 1, v
eq = [v for v in range(1 << 12) if g_rec(v) == v + 1]
print(f"(iii) 1 <= g(v) <= v + 1 for all v < 2^20; equality g(v) = v+1 for v < 2^12 at {len(eq)} values, first {eq[:12]}")
# sharper gamma bound: max over l < 2^m, both parities, of gamma
for m in range(1, 16):
    mx = max(g_rec(l) - 1 + (1 if l < (1 << (m - 1)) else 0) for l in range(1 << m))
    mx0 = max(g_rec(l) - 1 for l in range(1 << m))
    assert mx <= (1 << m) - 1
    if m <= 6 or m == 15:
        print(f"     m={m}: max gamma (h odd) = {mx}, (h even) = {mx0}, 2^m - 1 = {(1 << m) - 1}")

# (iv) Session 193 eq. (1) on words, brute force
def C_word(w):
    return sum((1 - int(w[j])) << j for j in range(len(w)))
bad = 0; cnt = 0
for a in range(1, 6):
    for b in range(1, 6):
        if a + b > 10: continue
        for hu in range(1 << (a - 1), 1 << a):
            u = format(hu, "0%db" % a)
            for lv in range(1 << b):
                v = format(lv, "0%db" % b)
                lhs = ham_paths(u + v)
                rhs = ham_paths(u) + ham_paths(v) - 1 + int(u, 2) * C_word(v) + (int(u[-1])) * (1 if v[0] == "0" else 0)
                cnt += 1; bad += lhs != rhs
                assert C_word(v) == Cm(lv, b)
print(f"(iv) Session 193 eq. (1) H(uv) = H(u)+H(v)-1+A(u)C(v)+p(u)eps(v), brute force, |u|,|v| <= 5: "
      f"{cnt} cases, {bad} violations; C(v) = C_m(l) on all of them")
