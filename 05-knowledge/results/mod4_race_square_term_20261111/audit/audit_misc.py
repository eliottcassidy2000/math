#!/usr/bin/env python3
"""audit_misc.py -- independent audit: Theorem S (1)-(3), the deck remarks, the subset form, zero cube sums by class.

Reads classes from `gentourng -q n` ('1' at (i,j), i<j, read as i -> j).

(S1/S2) for every class of even order n <= NS: search for an anti-automorphism (backtracking); if one exists, build an
        involutive one (psi = phi^k, k the odd part of the order of phi), check it has no fixed vertex, list ALL directed odd
        cycles explicitly, check that psi acts on them without fixed points, count the orbits {C, psi(C)} whose vertex sets
        meet (s), and check H = 1 + 4 s (mod 8), H computed by subset DP.  Also the odd-order analogue (H = 1 + 2F + 4I).
(S3)    for every involution psi with at most one fixed point (all of them for n <= 6, one representative otherwise):
        the number of labelled tournaments with anti-automorphism psi is 2^(n^2/4) resp. 2^((n^2-1)/4)  (n <= 6 by brute force).
(deck)  chi(T) = prod_v chi(T - v) for even n; prod_v chi(T - v) = 1 for odd n; Y_k = sum_l C(n-l,k-l) c_l and
        Y_k = (p_k - C(n,k))/2 (mod 2), p_k = number of directed paths on k vertices, for all k.
(zero)  number of classes of order m = 5, 9 with G(R) = 0 by the kernel criterion (K_1 non-empty, all weights 0 mod 4).
"""
import subprocess, sys, itertools
from math import comb

GENT = "/opt/homebrew/bin/gentourng"

def classes(n):
    out = subprocess.run([GENT, "-q", str(n)], capture_output=True, text=True).stdout.split("\n")
    res = []
    for line in out:
        line = line.strip()
        if len(line) != n * (n - 1) // 2:
            continue
        adj = [[0] * n for _ in range(n)]
        pos = 0
        for i in range(n):
            for j in range(i + 1, n):
                if line[pos] == "1": adj[i][j] = 1
                else: adj[j][i] = 1
                pos += 1
        res.append(adj)
    return res

def ham(adj, verts):
    """number of Hamiltonian paths of the subtournament induced on verts (subset DP)."""
    k = len(verts)
    if k == 0: return 1
    dp = [[0] * k for _ in range(1 << k)]
    for i in range(k): dp[1 << i][i] = 1
    for mask in range(1, 1 << k):
        for last in range(k):
            c = dp[mask][last]
            if not c: continue
            for w in range(k):
                if not (mask >> w) & 1 and adj[verts[last]][verts[w]]:
                    dp[mask | (1 << w)][w] += c
    return sum(dp[(1 << k) - 1])

def odd_cycles(adj):
    """all directed odd cycles, each as a tuple starting at its least vertex."""
    n = len(adj); res = []
    def dfs(start, path, vis):
        last = path[-1]
        if len(path) >= 3 and len(path) % 2 == 1 and adj[last][start]:
            res.append(tuple(path))
        for w in range(start + 1, n):
            if not (vis >> w) & 1 and adj[last][w]:
                path.append(w); dfs(start, path, vis | (1 << w)); path.pop()
    for v in range(n):
        dfs(v, [v], 1 << v)
    return res

def all_cycles_by_length(adj):
    n = len(adj); cnt = [0] * (n + 1)
    def dfs(start, last, ln, vis):
        if ln >= 3 and adj[last][start]: cnt[ln] += 1
        for w in range(start + 1, n):
            if not (vis >> w) & 1 and adj[last][w]:
                dfs(start, w, ln + 1, vis | (1 << w))
    for v in range(n): dfs(v, v, 1, 1 << v)
    return cnt

def find_anti(adj):
    n = len(adj); score = [sum(r) for r in adj]; img = [-1] * n; used = [False] * n
    def bt(v):
        if v == n: return True
        for w in range(n):
            if used[w] or score[w] != n - 1 - score[v]: continue
            if all(adj[u][v] == adj[w][img[u]] for u in range(v)):
                img[v] = w; used[w] = True
                if bt(v + 1): return True
                used[w] = False
        return False
    return img[:] if bt(0) else None

def perm_pow(p, k):
    r = list(range(len(p)))
    for _ in range(k): r = [p[x] for x in r]
    return r
def perm_order(p):
    k = 1; q = p[:]
    while q != list(range(len(p))): q = [p[x] for x in q]; k += 1
    return k

def canon_cycle(c):
    i = c.index(min(c)); return tuple(c[i:] + c[:i])

def check_S(NS):
    print("Theorem S (1), (2): all classes, orders 3 ..", NS)
    for n in range(3, NS + 1):
        cls = classes(n); nsc = 0; bad1 = 0; bad2 = 0; fixedcyc = 0; bad_odd = 0; hist = {}
        for adj in cls:
            phi = find_anti(adj)
            if phi is None: continue
            nsc += 1
            o = perm_order(phi); k = o
            while k % 2 == 0: k //= 2
            psi = perm_pow(phi, k)
            # psi is an anti-automorphism and an involution with at most one fixed vertex
            ok = all(adj[a][b] == adj[psi[b]][psi[a]] for a in range(n) for b in range(n) if a != b)
            ok = ok and perm_pow(psi, 2) == list(range(n))
            nfix = sum(1 for a in range(n) if psi[a] == a)
            if not ok or nfix != n % 2: bad1 += 1; continue
            cyc = odd_cycles(adj); cs = set(cyc)
            F = 0; I = 0; D = 0; seen = set()
            for c in cyc:
                if c in seen: continue
                img = canon_cycle([psi[x] for x in reversed(c)])      # psi(C): the image, traversed in the direction of T
                assert img in cs
                if img == c:
                    F += 1; seen.add(c)
                else:
                    seen.add(c); seen.add(img)
                    if set(c) & set(img): I += 1
                    else: D += 1
            H = ham(adj, list(range(n)))
            if n % 2 == 0:
                if F: fixedcyc += 1
                if H % 8 != (1 + 4 * I) % 8: bad2 += 1
                hist[H % 8] = hist.get(H % 8, 0) + 1
            else:
                if H % 8 != (1 + 2 * F + 4 * I) % 8: bad_odd += 1
                hist[H % 8] = hist.get(H % 8, 0) + 1
        print(f"  n={n}: classes {len(cls)}, self-converse {nsc}; failures of (1) {bad1}; classes with a psi-fixed odd cycle (even n) {fixedcyc}; "
              f"failures of H = 1 + 4s mod 8 (even n) {bad2}; failures of H = 1 + 2F + 4I mod 8 (odd n) {bad_odd}; H mod 8 of the self-converse classes {dict(sorted(hist.items()))}")

def check_S3():
    print("Theorem S (3): number of labelled tournaments with a given involutive anti-automorphism psi (all psi with at most one fixed point)")
    for n in range(2, 7):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        invs = []
        for p in itertools.permutations(range(n)):
            if all(p[p[i]] == i for i in range(n)) and sum(1 for i in range(n) if p[i] == i) == n % 2:
                invs.append(p)
        counts = set(); union = set()
        for psi in invs:
            c = 0
            for t in range(1 << len(pairs)):
                adj = [[0] * n for _ in range(n)]
                for b, (i, j) in enumerate(pairs):
                    if (t >> b) & 1: adj[i][j] = 1
                    else: adj[j][i] = 1
                if all(adj[a][b] == adj[psi[b]][psi[a]] for a, b in pairs):
                    c += 1; union.add(t)
            counts.add(c)
        e = n * n // 4 if n % 2 == 0 else (n * n - 1) // 4
        print(f"  n={n}: involutions {len(invs)}; tournaments per involution {sorted(counts)} (2^{e} = {2**e}); union over involutions = SC_n = {len(union)}")

def check_deck(NMAX):
    print("Deck remarks and subset form: all classes, orders 3 ..", NMAX)
    for n in range(3, NMAX + 1):
        cls = classes(n); bad_deck = 0; bad_Y = 0; bad_p = 0
        for adj in cls:
            V = list(range(n))
            H = ham(adj, V); chi = 1 if H % 4 == 1 else -1
            prod = 1
            for v in V:
                Hv = ham(adj, [u for u in V if u != v]); prod *= (1 if Hv % 4 == 1 else -1)
            if n % 2 == 0 and prod != chi: bad_deck += 1
            if n % 2 == 1 and prod != 1: bad_deck += 1
            cl = all_cycles_by_length(adj)
            if n <= 7:
                for k in range(1, n + 1):
                    Y = 0; p = 0
                    for S in itertools.combinations(V, k):
                        sub = [[adj[a][b] for b in S] for a in S]
                        cs = all_cycles_by_length(sub)
                        Y += sum(cs[l] for l in range(3, k + 1, 2))
                        p += ham(adj, list(S))
                    if Y != sum(comb(n - l, k - l) * cl[l] for l in range(3, k + 1, 2)): bad_Y += 1
                    if (p - comb(n, k)) % 2 != 0 or ((p - comb(n, k)) // 2 - Y) % 2 != 0: bad_p += 1
        print(f"  n={n}: classes {len(cls)}; failures deck {bad_deck}; Y_k formula {bad_Y}; Y_k = (p_k - C(n,k))/2 mod 2 {bad_p}")

def check_zero():
    print("Classes with G(R) = 0 by the kernel criterion (m = 1 mod 4): K_1 by elimination over F_2, whole coset enumerated")
    for m in (5, 9):
        cls = classes(m); z = 0; k1 = 0
        for adj in cls:
            # rows of A^T: row c = in-neighbours of c; solve A^T u = 1
            rows = [sum(1 << a for a in range(m) if adj[a][c]) | (1 << m) for c in range(m)]   # bit m = right-hand side 1
            where = [-1] * m; r = 0
            for col in range(m):
                p = next((i for i in range(r, m) if (rows[i] >> col) & 1), None)
                if p is None: continue
                rows[r], rows[p] = rows[p], rows[r]
                for i in range(m):
                    if i != r and (rows[i] >> col) & 1: rows[i] ^= rows[r]
                where[col] = r; r += 1
            if any((rows[i] >> m) & 1 for i in range(r, m)): continue
            k1 += 1
            u0 = sum(1 << col for col in range(m) if where[col] >= 0 and (rows[where[col]] >> m) & 1)
            basis = []
            for col in range(m):
                if where[col] < 0:
                    v = 1 << col
                    for c2 in range(m):
                        if where[c2] >= 0 and (rows[where[c2]] >> col) & 1: v |= 1 << c2
                    basis.append(v)
            allz = True
            for k in range(1 << len(basis)):
                u = u0
                for j, b in enumerate(basis):
                    if (k >> j) & 1: u ^= b
                # verify the solution directly
                assert all(bin(u & sum(1 << a for a in range(m) if adj[a][c])).count("1") & 1 for c in range(m))
                if bin(u).count("1") % 4: allz = False; break
            if allz: z += 1
        print(f"  m={m}: classes {len(cls)}; K_1 non-empty in {k1}; G(R) = 0 in {z}")

if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("S", "all"): check_S(int(sys.argv[2]) if len(sys.argv) > 2 else 8)
    if what in ("S3", "all"): check_S3()
    if what in ("deck", "all"): check_deck(7)
    if what in ("zero", "all"): check_zero()
