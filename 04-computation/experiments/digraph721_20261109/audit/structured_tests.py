"""AUDIT D1: structured symmetric test graphs (n <= 11) for cross-checking aw_verify's stabiliser-chain parity
against brute force over all permutations ("aw_verify bf").  Each graph is emitted under a random relabelling
too.  Output format: tag E n 0 | u-v ...  (k = 0 is a dummy; only the BF lines matter)."""
import itertools, random
random.seed(20261109)
out = []
def emit(tag, n, edges):
    out.append((tag, n, sorted(set(tuple(sorted(e)) for e in edges))))
    p = list(range(n)); random.shuffle(p)
    out.append((tag + "_rl", n, sorted(set(tuple(sorted((p[u], p[v]))) for u, v in edges))))
def comp(n, edges):
    E = set(tuple(sorted(e)) for e in edges)
    return [(u, v) for u in range(n) for v in range(u + 1, n) if (u, v) not in E]
for t in range(1, 6):            # cocktail party CP(2t) = E_2 v ... v E_2 (even), and K_1 v CP(2t)
    n = 2 * t; E = comp(n, [(2 * i, 2 * i + 1) for i in range(t)]); emit(f"CP{n}", n, E)
    if n + 1 <= 11: emit(f"K1vCP{n}", n + 1, E + [(n, v) for v in range(n)])
for r in range(1, 6):            # K_{r,r}: odd iff r odd
    emit(f"K{r}{r}", 2 * r, [(i, r + j) for i in range(r) for j in range(r)])
for m in range(3, 12): emit(f"C{m}", m, [(i, (i + 1) % m) for i in range(m)])
for m in range(2, 12): emit(f"P{m}", m, [(i, i + 1) for i in range(m - 1)])
for m in range(2, 11): emit(f"Star{m}", m + 1, [(0, i) for i in range(1, m + 1)])
pet = [(i, (i + 1) % 5) for i in range(5)] + [(5 + i, 5 + (i + 2) % 5) for i in range(5)] + [(i, 5 + i) for i in range(5)]
emit("Petersen", 10, pet); emit("coPetersen", 10, comp(10, pet))
q3 = [(a, b) for a in range(8) for b in range(a + 1, 8) if bin(a ^ b).count("1") == 1]
emit("Q3", 8, q3); emit("coQ3", 8, comp(8, q3))
emit("3P3", 9, [(0, 1), (1, 2), (3, 4), (4, 5), (6, 7), (7, 8)])
emit("2P3uK1", 7, [(0, 1), (1, 2), (3, 4), (4, 5)])
emit("K4-e", 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)])
emit("Wheel6", 7, [(i, (i + 1) % 6) for i in range(6)] + [(6, i) for i in range(6)])
emit("K3xK3", 9, [(a, b) for a in range(9) for b in range(a + 1, 9) if (a // 3 == b // 3) != (a % 3 == b % 3)])
emit("E_6", 6, [])
emit("2C5", 10, [(i, (i + 1) % 5) for i in range(5)] + [(5 + i, 5 + (i + 1) % 5) for i in range(5)])
emit("C4uC4uK2", 10, [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (8, 9)])
emit("P4vP4", 8, [(0, 1), (1, 2), (2, 3), (4, 5), (5, 6), (6, 7)] + [(i, j) for i in range(4) for j in range(4, 8)])
emit("E3vE3", 6, [(i, j) for i in range(3) for j in range(3, 6)])
emit("E2vE2vP3", 7, comp(7, [(0, 1), (2, 3), (4, 6)]))
# hub-move images of small even graphs (PH(G,u,l) with l = 1, 2)
def ph(n, edges, u, l):
    h = n; E = list(edges) + [(v, h) for v in range(n) if v != u] + [(h, n + 1 + i) for i in range(l)]
    return n + 1 + l, E
for (name, n, E) in [("P3uK1", 4, [(0, 1), (1, 2)]), ("P3uP3", 6, [(0, 1), (1, 2), (3, 4), (4, 5)]),
                     ("C4", 4, [(0, 1), (1, 2), (2, 3), (3, 0)])]:
    for u in range(n):
        for l in (1, 2, 3):
            m, F = ph(n, E, u, l)
            if m <= 11: emit(f"PH_{name}_{u}_{l}", m, F)
with open("structured_tests.txt", "w") as f:
    for tag, n, E in out:
        f.write(f"{tag} E {n} 0 | " + " ".join(f"{u}-{v}" for u, v in E) + "\n")
print(len(out), "graphs written")
