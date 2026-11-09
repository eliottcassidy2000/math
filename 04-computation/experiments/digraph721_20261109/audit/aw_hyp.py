"""AUDIT D1 (HYP-9625): explicit graphs for the reduction's named values and the arithmetic facts.
Writes aw_hyp_graphs.txt (aw_verify format; set E is a dummy - only d and the reported parity matter here) and
prints the arithmetic facts.  Graphs (vertex-labelled edge lists):
  25   = K_2 v E_3;  185 = cone over (K_2 v P_4) u K_1;  517 = paw u T_47;  1369 = cone over (K_2 u T_115) u K_1;
  THM-5056 seeds (light check for HYP-9625's new status): K_4 - e, K_2 u C_4, K_2 v (P_3 u E_2), K_2 u T_55 and
  their PH(w, 2) images (69, 257, 549, 1061)."""
def thr(k, m):                      # threshold graph on m vertices with d = k (THM-4996 Lemma B1)
    assert k % 2 == 1 and 1 <= k < 2**m
    if m == 1: return []
    if k < 2**(m - 1): return thr(k, m - 1)                      # pad
    return thr(k - 2**(m - 1), m - 1) + [(v, m - 1) for v in range(m - 1)]   # cone (apex m-1)
def cone(g, m): return g + [(v, m) for v in range(m)]
def k2join(g, m): return g + [(m, m + 1)] + [(v, m) for v in range(m)] + [(v, m + 1) for v in range(m)]
def union(g1, m1, g2): return g1 + [(a + m1, b + m1) for a, b in g2]
def ph(g, m, u, l): return g + [(v, m) for v in range(m) if v != u] + [(m, m + 1 + i) for i in range(l)]
G = []
G.append(("25:K2vE3", 5, 25, k2join([], 3)))
P4 = [(0, 1), (1, 2), (2, 3)]
G.append(("185:cone(K2vP4+K1)", 8, 185, cone(k2join(P4, 4), 7)))
paw = [(0, 1), (0, 2), (0, 3), (1, 2)]
G.append(("517:paw+T47", 10, 517, union(paw, 4, thr(47, 6))))
G.append(("1369:cone(K2+T115+K1)", 11, 1369, cone(union([(0, 1)], 2, thr(115, 7)), 10)))
K4e = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]                    # x, y = 0, 1; w = 2
G.append(("seed:K4-e", 4, 13, K4e)); G.append(("seed:PH(K4-e,w,2)", 7, 69, ph(K4e, 4, 2, 2)))
K2C4 = [(0, 1), (2, 3), (3, 4), (4, 5), (5, 2)]                   # w = 2
G.append(("seed:K2uC4", 6, 33, K2C4)); G.append(("seed:PH(K2uC4,w,2)", 9, 257, ph(K2C4, 6, 2, 2)))
K2vP3E2 = k2join([(0, 1), (1, 2)], 5)                             # w = 3 (in E_2)
G.append(("seed:K2v(P3uE2)", 7, 101, K2vP3E2)); G.append(("seed:PH(K2v(P3uE2),w,2)", 10, 549, ph(K2vP3E2, 7, 3, 2)))
T55 = thr(55, 6); K2T55 = union([(0, 1)], 2, T55)
deg = [0] * 8
for a, b in K2T55: deg[a] += 1; deg[b] += 1
w = [v for v in range(2, 8) if deg[v] == 2][0]
G.append(("seed:K2uT55", 8, 165, K2T55)); G.append(("seed:PH(K2uT55,w,2)", 11, 1061, ph(K2T55, 8, w, 2)))
with open("aw_hyp_graphs.txt", "w") as f:
    for tag, n, k, g in G:
        f.write(f"{tag} E {n} {k} | " + " ".join(f"{a}-{b}" for a, b in sorted(set(tuple(sorted(e)) for e in g))) + "\n")
def isprime(p):
    if p < 2: return False
    i = 2
    while i * i <= p:
        if p % i == 0: return False
        i += 1
    return True
def fac(p):
    out = []; i = 2
    while i * i <= p:
        while p % i == 0: out.append(i); p //= i
        i += 1
    if p > 1: out.append(p)
    return out
for v in (517, 1025, 1029, 1061, 2053, 4097, 4133, 4613, 549, 65537, 1369, 185, 25, 11, 23, 47, 37, 5):
    print(v, "prime" if isprime(v) else "= " + " * ".join(map(str, fac(v))))
print("37^2 =", 37**2, "; 2^10 + 3*115 =", 2**10 + 3 * 115)
