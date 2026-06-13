"""
paley_H23_monad.py  (monad-explorer)

Fast int64 Held-Karp Hamiltonian-path counter for the Paley tournament, to extend
the canonical ratio sequence r(p)=H(T_p)/|Aut(T_p)| to p=23 (the prime named in
the-tessellation.md).  H(T_23) ~ 1.6e16 < int64 max (9.2e18), and every partial dp
value is <= paths_from_0 = H/p ~ 7e14, so int64 is exact and overflow-free.

Memory trick: vertex 0 is always in the path (start fixed by vertex-transitivity,
H = p * #paths-from-0), so index dp rows by (mask >> 1) -> 2^(n-1) rows.
For n=23: 2^22 * 23 * 8 bytes = 770 MB.

Validated against H(T_7)=189, H(T_11)=95095, H(T_19)=1172695746915.
"""
import sys, time
import numpy as np
from sympy import factorint

def paley_adj(p):
    assert p % 4 == 3
    Q = {(x*x) % p for x in range(1, p)}
    A = np.zeros((p, p), dtype=np.int64)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in Q:
                A[i, j] = 1
    return A

def count_ham_paths_i64(A):
    n = A.shape[0]
    rows = 1 << (n - 1)              # masks containing bit 0, indexed by mask>>1
    dp = np.zeros((rows, n), dtype=np.int64)
    dp[0, 0] = 1                     # mask = {0}
    A = A.astype(np.int64)
    t0 = time.time()
    for r in range(rows):
        vec = dp[r]
        # quick skip of empty rows
        if not vec.any():
            continue
        mask = (r << 1) | 1
        contrib = vec @ A           # contrib[w] = sum_{v: v->w} dp[mask][v]
        for w in range(1, n):       # w != 0 (0 already in mask)
            if mask & (1 << w):
                continue
            c = contrib[w]
            if c:
                nr = (mask | (1 << w)) >> 1
                dp[nr, w] += c
        if r and r % (1 << 20) == 0:
            print(f"   ... {r}/{rows} rows  ({time.time()-t0:.0f}s)", flush=True)
    full_row = rows - 1
    paths_from_0 = int(dp[full_row].sum())
    return n * paths_from_0

def classify_mod3(q):
    if q == 3: return "ramified"
    return "split" if q % 3 == 1 else "inert"

if __name__ == "__main__":
    KNOWN = {7: 189, 11: 95095, 19: 1172695746915}
    ps = [int(x) for x in sys.argv[1:]] or [7, 11, 19, 23]
    for p in ps:
        t0 = time.time()
        print(f"[computing H(T_{p}) ...]", flush=True)
        H = count_ham_paths_i64(paley_adj(p))
        dt = time.time() - t0
        tag = ""
        if p in KNOWN:
            tag = "  PASS" if H == KNOWN[p] else f"  FAIL (expected {KNOWN[p]})"
        print(f"H(T_{p}) = {H}   [{dt:.1f}s]{tag}")
        aut = p * (p - 1) // 2
        if H % aut == 0:
            r = H // aut
            fr = factorint(r)
            cls = {q: classify_mod3(q) for q in fr}
            allsplit = all(c != "inert" for c in cls.values())
            print(f"   |Aut|={aut}  r(p)={r}")
            print(f"   factor(r)={fr}")
            print(f"   classes mod3={cls}")
            print(f"   completely split in Q(sqrt-3)? {allsplit}   smooth(<=100)? {bool(fr) and max(fr)<=100}")
            print(f"   factor(H)={factorint(H)}")
        sys.stdout.flush()
