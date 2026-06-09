"""
kind-pasteur-2026-06-09-S2 : BRANCH III (dyadic-gap hunt, HYP-2359 / Erdos-Gyarfas)
Part 3: STRUCTURED CANDIDATES.
  (A) Generalized Petersen graphs GP(n,k), 3<=n<=40, 1<=k<n/2 (2n vertices, <=80):
      girth, #C4, #C8; for C8-free ones also #C16 (and #C32 if 2n>=32, budgeted).
  (B) Cubic Cayley graphs:
      - circulants Cay(Z_n, {+-a, n/2}), n even, 6<=n<=80
      - dihedral D_m, 3 reflections {s, s r^j, s r^k}, 2m<=80
      - dihedral D_m, reflection + rotation {s, r^a, r^-a}, 2m<=80
      - Z_2 x Z_m (m even), S={(1,0),(0,+-a)}
      girth, #C4, #C8; C8-free ones get the full dyadic profile.

Output -> 05-knowledge/results/dyadic_gp_cayley_kps2.out
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dyadic_cycle_checker_kps2 import count_cycles_len, girth, edges_to_adj
from dyadic_gap_hunt_kps2 import count_cycles_capped, is_connected

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))
def save():
    with open("05-knowledge/results/dyadic_gp_cayley_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

def GP(n, k):
    edges = []
    for i in range(n):
        edges.append((i, (i + 1) % n))            # outer
        edges.append((i, n + i))                  # spoke
        edges.append((n + i, n + (i + k) % n))    # inner
    # dedupe (inner loop may double-add when k=n/2; we exclude that case anyway)
    es = set(frozenset(e) for e in edges)
    return edges_to_adj(2 * n, [tuple(e) for e in es])

def dyadic_row(adj, deep=True):
    n = len(adj)
    g = girth(adj)
    c4 = count_cycles_len(adj, 4)
    c8 = count_cycles_len(adj, 8)
    c16 = c32 = None
    # deep dyadic profile only for CONNECTED C8-free graphs (disconnected ones reduce
    # to their components; not counterexample candidates per se)
    if c8 == 0 and deep and is_connected(adj):
        c16 = count_cycles_capped(adj, 16)
        if n >= 32:
            c32 = count_cycles_capped(adj, 32, budget=10_000_000)
    return g, c4, c8, c16, c32

def main():
    t0 = time.time()
    log("=" * 100)
    log("(A) GENERALIZED PETERSEN GRAPHS GP(n,k), 3<=n<=40, 1<=k<n/2  (2n vertices)")
    log("=" * 100)
    c8free = []
    rows = 0
    for n in range(3, 41):
        line = []
        for k in range(1, (n - 1) // 2 + 1):
            adj = GP(n, k)
            if any(len(a) != 3 for a in adj):
                line.append(f"k={k}:NONCUBIC")
                continue
            g, c4, c8, c16, c32 = dyadic_row(adj)
            rows += 1
            tag = f"k={k}:g{g},c4={c4},c8={c8}"
            if c8 == 0:
                tag += f",c16={c16}" + (f",c32={c32}" if c32 is not None else "")
                c8free.append((n, k, g, c4, c8, c16, c32))
            line.append(tag)
        log(f"GP(n={n:2d}): " + "  ".join(line))
        save()
    log(f"\n[{rows} GP graphs checked]")
    log("\nGP graphs with NO 8-cycle:")
    if not c8free:
        log("  NONE -- every GP(n,k) with 3<=n<=40 contains an 8-cycle.")
    for (n, k, g, c4, c8, c16, c32) in c8free:
        log(f"  GP({n},{k}): girth={g} c4={c4} c8=0 c16={c16} c32={c32}"
            + ("  <-- ALSO C4-FREE" if c4 == 0 else ""))
    save()

    log("\n" + "=" * 100)
    log("(B) CUBIC CAYLEY GRAPHS (order <= 80)")
    log("=" * 100)
    findings = []

    # --- circulants Cay(Z_n, {a,-a,n/2})
    log("\n(B1) circulants Cay(Z_n, {+-a, n/2}), n even:")
    cnt = 0
    free_list = []
    for n in range(6, 81, 2):
        for a in range(1, n // 2):
            edges = set()
            for v in range(n):
                edges.add(frozenset((v, (v + a) % n)))
                edges.add(frozenset((v, (v + n // 2) % n)))
            adj = edges_to_adj(n, [tuple(e) for e in edges])
            if any(len(x) != 3 for x in adj):
                continue
            cnt += 1
            g, c4, c8, c16, c32 = dyadic_row(adj)
            if c8 == 0:
                free_list.append((f"Cay(Z_{n},{{+-{a},{n//2}}})", len(adj), g, c4, c16, c32,
                                  is_connected(adj)))
    log(f"  checked {cnt}; C8-free: {len(free_list)}")
    for t in free_list:
        log(f"    {t[0]}: n={t[1]} girth={t[2]} c4={t[3]} c16={t[4]} c32={t[5]} conn={t[6]}")
    findings += free_list
    save()

    # --- dihedral, three reflections {s, s r^j, s r^k}
    log("\n(B2) dihedral D_m, three reflections {s, s r^j, s r^k} (0<j<k<m), 2m<=80:")
    cnt = 0
    free_list = []
    for m in range(3, 41):
        seen = set()
        for j in range(1, m):
            for k in range(j + 1, m):
                # canonical class: multiset of cyclic gaps {j, k-j, m-k} up to rotation/reflection
                gaps = tuple(sorted((j, k - j, m - k)))
                if gaps in seen:
                    continue
                seen.add(gaps)
                # group elements: 0..m-1 = rotations r^i ; m..2m-1 = reflections s r^(i-m)
                # right-multiplication Cayley graph by involutions {s r^t} for t in {0, j, k}
                # (r^i)*(s r^t) = s r^(t-i) ; (s r^i)*(s r^t) = r^(t-i)
                edges = set()
                for t in (0, j, k):
                    for i in range(m):
                        edges.add(frozenset((i, m + (t - i) % m)))
                adj = edges_to_adj(2 * m, [tuple(e) for e in edges])
                if any(len(x) != 3 for x in adj):
                    continue
                cnt += 1
                g, c4, c8, c16, c32 = dyadic_row(adj)
                if c8 == 0:
                    free_list.append((f"D_{m} refl(0,{j},{k})", 2 * m, g, c4, c16, c32,
                                      is_connected(adj)))
        save()
    log(f"  checked {cnt}; C8-free: {len(free_list)}")
    for t in free_list:
        log(f"    {t[0]}: n={t[1]} girth={t[2]} c4={t[3]} c16={t[4]} c32={t[5]} conn={t[6]}")
    findings += free_list
    save()

    # --- dihedral, reflection + rotation {s, r^a, r^-a}
    log("\n(B3) dihedral D_m, {s, r^a, r^-a} (1<=a<m/2), 2m<=80:")
    cnt = 0
    free_list = []
    for m in range(3, 41):
        for a in range(1, (m + 1) // 2):
            edges = set()
            for i in range(m):
                edges.add(frozenset((i, (i + a) % m)))                  # rotations among rotations
                edges.add(frozenset((m + i, m + (i - a) % m)))          # (s r^i)*r^a = s r^(i+a)? careful
            # right multiplication: (r^i)*r^a=r^(i+a); (s r^i)*r^a = s r^(i+a).
            edges = set()
            for i in range(m):
                edges.add(frozenset((i, (i + a) % m)))
                edges.add(frozenset((m + i, m + (i + a) % m)))
                # by s: (r^i)*s = s r^(-i); (s r^i)*s = r^(-i)... right-mult by s: x -> x s
                # r^i s = s r^{-i}; index of s r^{-i} is m + (-i)%m
                edges.add(frozenset((i, m + (-i) % m)))
            adj = edges_to_adj(2 * m, [tuple(e) for e in edges])
            if any(len(x) != 3 for x in adj):
                continue
            cnt += 1
            g, c4, c8, c16, c32 = dyadic_row(adj)
            if c8 == 0:
                free_list.append((f"D_{m} {{s,r^+-{a}}}", 2 * m, g, c4, c16, c32,
                                  is_connected(adj)))
    log(f"  checked {cnt}; C8-free: {len(free_list)}")
    for t in free_list:
        log(f"    {t[0]}: n={t[1]} girth={t[2]} c4={t[3]} c16={t[4]} c32={t[5]} conn={t[6]}")
    findings += free_list
    save()

    # --- Z_2 x Z_m, S = {(1,0), (0,a), (0,-a)}
    log("\n(B4) Z_2 x Z_m, S={(1,0),(0,+-a)}, m even (odd m isomorphic to circulant), 2m<=80:")
    cnt = 0
    free_list = []
    for m in range(4, 41, 2):
        for a in range(1, (m + 1) // 2):
            edges = set()
            for b in range(2):
                for i in range(m):
                    v = b * m + i
                    edges.add(frozenset((v, b * m + (i + a) % m)))
                    edges.add(frozenset((v, (1 - b) * m + i)))
            adj = edges_to_adj(2 * m, [tuple(e) for e in edges])
            if any(len(x) != 3 for x in adj):
                continue
            cnt += 1
            g, c4, c8, c16, c32 = dyadic_row(adj)
            if c8 == 0:
                free_list.append((f"Z2xZ_{m} a={a}", 2 * m, g, c4, c16, c32, is_connected(adj)))
    log(f"  checked {cnt}; C8-free: {len(free_list)}")
    for t in free_list:
        log(f"    {t[0]}: n={t[1]} girth={t[2]} c4={t[3]} c16={t[4]} c32={t[5]} conn={t[6]}")
    findings += free_list

    log("\n" + "=" * 100)
    log("SUMMARY (B): C8-free cubic Cayley graphs found: " + str(len(findings)))
    log(f"total time {time.time()-t0:.0f}s")
    save()

if __name__ == "__main__":
    main()
