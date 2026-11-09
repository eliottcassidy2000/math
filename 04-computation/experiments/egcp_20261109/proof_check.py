"""Mechanical check of the case analysis in the S89 proof of the sharp even-graph Chan-Pak theorem.

Claim_3(n): every odd x <= 2^n - n - 3 with x not in {3, 7, 13} is d(G) for an even graph on <= n vertices.
Base n <= 9: census (E(n) = O(n), verified separately).
Step n >= 10, x odd, 2^(n-1) - n - 2 < x <= 2^n - n - 3:
  Case 1  x < 2^(n-1), x = 2^(n-1) - j:   j in T(n) (top band, leaf graphs) or x = M(n-1) (cocktail party).
  Case 2a D = 2^n - x <= 2^(n-2) + 1:       E_2^t v G, g = least g = n mod 2 with 2^g + n - g - 1 >= D,
                                            d(G) = 2^g - D + n - g in Claim_3(g) range and not in {3,7,13};
                                            else g + 2 (if <= n - 2); else x = 3*2^(n-2) + {1,5,11} via cone over L(n-1).
  Case 2b D >= 2^(n-2) + 3:                 cone over G on n-1 vertices without universal vertex, d(G) = y = 2^(n-1) - D:
                                            y in Claim_3(n-2) range (padded), y in {3,7,13} (x via L(n)),
                                            y = 2^(n-2) - j with j in T(n-1) or y = M(n-2) (padded cocktail party).
The script checks, for every x in the step range, that one of the branches applies with all side conditions,
and that the leaf-family bookkeeping (top-band j-sets and offsets c) is supplied by the base + SR + frontier families.
"""
import sys

def M(n): return 2**n - n if n % 2 else 2**n - n - 1
def Tj(m):
    """j-set of the top band at size m (values 2^(m-1) - j needing exactly m vertices)"""
    N = m - 1
    if N % 2 == 0: return set(range(1, N, 2))           # j <= N - 1
    return set(range(1, N - 1, 2)) | {N + 2}            # j <= N - 2, plus X(N)

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20

# leaf-family bookkeeping: LT[m] = top-band j's realised by leaf graphs (no universal vertex) on m vertices;
# LC[m] = offsets c with 2^(m-1) + c realised by leaf graphs on m vertices
LT = {8: {1, 3, 5, 9}, 9: {1, 3, 5, 7}}                 # census (verify_base.py)
LC = {8: {1, 3, 5, 7, 11, 13}, 9: {1, 3, 5, 7, 11, 13}}  # census
for m in range(10, NMAX + 1):
    fr = set()
    if m % 2 == 0:
        if m >= 6: fr.add(m - 3)                         # P_{rho-1} u K_1, z endpoint, rho = m/2
        if m >= 10: fr.add(m + 1)                        # spider S(l1,1,2), z end of leg 1, rho = m/2 >= 5
    else:
        if m >= 9: fr.add(m - 2)                         # P_{rho-1} u K_1, z off-centre interior, rho = (m+1)/2 >= 5
    LT[m] = LT[m - 2] | fr                               # SR preserves j
    LC[m] = set(LC[m - 2])                               # SR preserves c
    assert Tj(m) <= LT[m], (m, Tj(m) - LT[m])

def in_claim3(g, v):
    return v >= 1 and v % 2 == 1 and v <= 2**g - g - 3 and v not in (3, 7, 13) or (g <= 9 and v in base_E(g))

_base = {}
def base_E(g):
    # exact census sets for g <= 9: O(g)
    if g not in _base:
        s = set(range(1, M(g) + 1, 2)) - {3, 7, 13}
        if g % 2 and g >= 5: s.discard(2**g - g - 2)
        _base[g] = s
    return _base[g]

stats = {}
for n in range(10, NMAX + 1):
    lo = 2**(n - 1) - n - 2
    hi = 2**n - n - 3
    cnt = {'1T': 0, '1M': 0, '2a': 0, '2a+2': 0, '2aL': 0, '2bpad': 0, '2bL': 0, '2bT': 0, '2bM': 0}
    x = lo + 1 if (lo + 1) % 2 else lo + 2
    while x <= hi:
        if x < 2**(n - 1):
            j = 2**(n - 1) - x
            if j in LT[n]: cnt['1T'] += 1
            elif x == M(n - 1): cnt['1M'] += 1
            else: raise SystemExit("Case 1 fails n=%d x=%d j=%d" % (n, x, j))
        else:
            D = 2**n - x
            assert n + 3 <= D <= 2**(n - 1) - 1
            if D <= 2**(n - 2) + 1:
                g = 1 if n % 2 else 2
                while 2**g + n - g - 1 < D: g += 2
                assert g <= n - 2
                dG = 2**g - D + n - g
                assert dG >= 1 and dG % 2 == 1 and dG <= 2**g - g - 3 + (2 if g <= 9 else 0)
                if in_claim3(g, dG):
                    cnt['2a'] += 1
                elif g + 2 <= n - 2:
                    d2 = 2**(g + 2) - D + n - g - 2
                    if in_claim3(g + 2, d2): cnt['2a+2'] += 1
                    else: raise SystemExit("Case 2a(+2) fails n=%d x=%d" % (n, x))
                else:
                    c = x - 3 * 2**(n - 2)
                    if c in LC[n - 1] and x == 2**(n - 1) + 2**(n - 2) + c: cnt['2aL'] += 1
                    else: raise SystemExit("Case 2a(L) fails n=%d x=%d c=%d" % (n, x, c))
            else:
                y = 2**(n - 1) - D
                assert 1 <= y <= 2**(n - 2) - 3
                if y <= 2**(n - 2) - (n - 2) - 3 and y not in (3, 7, 13): cnt['2bpad'] += 1
                elif y in (3, 7, 13):
                    if y in LC[n]: cnt['2bL'] += 1
                    else: raise SystemExit("Case 2b(L) fails n=%d x=%d" % (n, x))
                else:
                    j = 2**(n - 2) - y
                    if j in LT[n - 1]: cnt['2bT'] += 1
                    elif y == M(n - 2): cnt['2bM'] += 1
                    else: raise SystemExit("Case 2b(T) fails n=%d x=%d j=%d" % (n, x, j))
        x += 2
    # top values beyond Claim_3: M(n) explicit; X(n) excluded
    stats[n] = cnt
    print(n, cnt, flush=True)
print("all cases verified for 10 <= n <=", NMAX)
