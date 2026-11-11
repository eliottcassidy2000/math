# audit-5129, step 1: brute force over ALL labelled tournaments of order <= 6 and every root.
#   - literal (permutation) counts of D, E, U_even, U_odd, Pi_2 against the subset programme
#   - E = U_odd - U_even, U >= D
#   - closing laws H(B + z) = D - E + U (= D + 2 U_even) and H(B => TT_t => z) = D + 2^t U, t = 1, 2, 3,
#     on the explicit closed tournaments (literal counts wherever affordable, subset programme otherwise)
#   - appended sink keeps D and doubles U
#   - product law (A1) on every non-strong tournament, components counted literally
#   - the sets of pairs (D, U) and triples (D, E, U)
import sys, time
import numpy as np
import acore as A

t0 = time.time()
NMAX = 6
perms = {k: A.perm_array(k) for k in range(1, 9)}
rng = np.random.default_rng(5129)

out_pairs = {}
out_triples = {}
tot_roots = 0
fail = {"lit_vs_dp": 0, "E_parity": 0, "U_ge_D": 0, "head": 0, "head_even": 0, "tail1": 0, "tail2": 0,
        "tail3": 0, "pi2": 0, "pi2_sum": 0, "sink": 0, "product_D": 0, "product_U": 0}
cnt = {"head_literal": 0, "tail1_literal": 0, "tail2_literal": 0, "tail1_dp": 0, "tail2_dp": 0, "tail3_dp": 0,
       "nonstrong": 0, "product_roots": 0, "root_trivial_component": 0, "sink_literal": 0, "sink_dp": 0,
       "U_eq_D_nonsink": 0, "U_eq_D_sink": 0}

for n in range(1, NMAX + 1):
    npairs = n * (n - 1) // 2
    adj = np.zeros(n, np.int64)
    hp = np.zeros((1 << n, n), np.int64)
    sp = np.zeros((1 << n, n), np.int64)
    Hs = np.zeros(1 << n, np.int64)
    E = np.zeros(n, np.int64)
    Ue = np.zeros(n, np.int64)
    Uo = np.zeros(n, np.int64)
    st = np.zeros(4, np.int64)
    comp = np.zeros(n, np.int64)
    cadj = np.zeros(12, np.int64)
    dps = {m: np.zeros((1 << m, m), np.int64) for m in range(2, 11)}
    pairs = set()
    triples = set()
    for code in range(1 << npairs):
        A.adj_from_code(code, n, adj)
        Hlit = A.literal_stats(adj, n, perms[n], E, Ue, Uo)
        A.tables(adj, n, hp, sp, Hs)
        p2 = A.pi2_from_tables(n, Hs)
        if p2 != A.literal_pi2(adj, n, perms[n]):
            fail["pi2"] += 1
        if p2 != int((Ue + Uo).sum()) + Hlit:
            fail["pi2_sum"] += 1
        strong = A.is_strong(adj, n)
        ncomp = A.components(adj, n, comp)
        if strong != (ncomp == 1):
            raise SystemExit("component routine inconsistent")
        # literal component data for the product law
        if ncomp > 1:
            cnt["nonstrong"] += 1
            cH = []
            cPi = []
            cU = []
            cverts = []
            for c in range(ncomp):
                vs = [v for v in range(n) if comp[v] == c]
                k = len(vs)
                sub = np.zeros(k, np.int64)
                for a_, v in enumerate(vs):
                    for b_, w in enumerate(vs):
                        if (adj[v] >> w) & 1:
                            sub[a_] |= (1 << b_)
                e_ = np.zeros(k, np.int64); ue_ = np.zeros(k, np.int64); uo_ = np.zeros(k, np.int64)
                h_ = A.literal_stats(sub, k, perms[k], e_, ue_, uo_)
                cH.append(h_)
                cPi.append(A.literal_pi2(sub, k, perms[k]))
                cU.append({v: int(ue_[a_] + uo_[a_]) for a_, v in enumerate(vs)})
                cverts.append(vs)
            # the condensation order returned must satisfy C_i => C_j for i < j
            for v in range(n):
                for w in range(n):
                    if comp[v] < comp[w] and not ((adj[v] >> w) & 1):
                        raise SystemExit("condensation order wrong")
            Dprod = 1
            for h_ in cH:
                Dprod *= h_
            if Dprod != Hlit:
                fail["product_D"] += 1
        for y in range(n):
            tot_roots += 1
            A.root_stats(n, hp, sp, Hs, y, st)
            D, Ey, ue, uo = int(st[0]), int(st[1]), int(st[2]), int(st[3])
            U = ue + uo
            if not (D == Hlit and Ey == E[y] and ue == Ue[y] and uo == Uo[y]):
                fail["lit_vs_dp"] += 1
            if Ey != uo - ue:
                fail["E_parity"] += 1
            if U < D:
                fail["U_ge_D"] += 1
            if U == D:
                if adj[y] == 0:
                    cnt["U_eq_D_sink"] += 1
                else:
                    cnt["U_eq_D_nonsink"] += 1
            pairs.add((D, U))
            triples.add((D, Ey, U))
            # head
            A.head_adj(adj, n, y, cadj)
            Hh = A.literal_H(cadj, n + 1, perms[n + 1])
            cnt["head_literal"] += 1
            if Hh != D - Ey + U:
                fail["head"] += 1
            if Hh != D + 2 * ue:
                fail["head_even"] += 1
            # tails
            for t in (1, 2, 3):
                m = n + t + 1
                A.tail_adj(adj, n, y, t, cadj)
                if m <= 8 and (m <= 7 or n <= 5):
                    Ht = A.literal_H(cadj, m, perms[m])
                    cnt["tail%d_literal" % t] = cnt.get("tail%d_literal" % t, 0) + 1
                    # also the subset programme, to validate it
                    if Ht != A.ham_count(cadj, m, dps[m]):
                        raise SystemExit("ham_count disagrees with the literal count")
                else:
                    Ht = A.ham_count(cadj, m, dps[m])
                    cnt["tail%d_dp" % t] = cnt.get("tail%d_dp" % t, 0) + 1
                if Ht != D + (1 << t) * U:
                    fail["tail%d" % t] += 1
            # appended sink
            for v in range(n):
                cadj[v] = adj[v] | (1 << n)
            cadj[n] = 0
            if n + 1 <= 6:
                e2 = np.zeros(n + 1, np.int64); ue2 = np.zeros(n + 1, np.int64); uo2 = np.zeros(n + 1, np.int64)
                h2 = A.literal_stats(cadj, n + 1, perms[n + 1], e2, ue2, uo2)
                cnt["sink_literal"] += 1
                if not (h2 == D and ue2[y] + uo2[y] == 2 * U):
                    fail["sink"] += 1
            else:
                hp2 = np.zeros((1 << (n + 1), n + 1), np.int64); sp2 = np.zeros((1 << (n + 1), n + 1), np.int64)
                Hs2 = np.zeros(1 << (n + 1), np.int64)
                A.tables(cadj, n + 1, hp2, sp2, Hs2)
                A.root_stats(n + 1, hp2, sp2, Hs2, y, st)
                cnt["sink_dp"] += 1
                if not (st[0] == D and st[2] + st[3] == 2 * U):
                    fail["sink"] += 1
            # product law
            if ncomp > 1:
                cnt["product_roots"] += 1
                k = int(comp[y])
                pred = cU[k][y]
                if len(cverts[k]) == 1:
                    cnt["root_trivial_component"] += 1
                    if cU[k][y] != 1:
                        raise SystemExit("U of a point is not 1")
                for i in range(ncomp):
                    if i < k:
                        pred *= cH[i]
                    elif i > k:
                        pred *= cPi[i]
                if pred != U:
                    fail["product_U"] += 1
    out_pairs[n] = pairs
    out_triples[n] = triples
    print("order %d: %d labelled tournaments, %d distinct pairs (D,U), %d distinct triples (D,E,U); elapsed %.1fs"
          % (n, 1 << npairs, len(pairs), len(triples), time.time() - t0), flush=True)

print()
print("rooted labelled tournaments checked:", tot_roots)
print("counters:", cnt)
print("failures:", fail)
print("ALL CHECKS PASS" if all(v == 0 for v in fail.values()) else "FAILURES PRESENT")

# extra literal tail checks at order 6 (random sample): t = 1 (order 8), t = 2 (order 9)
n = 6
adj = np.zeros(n, np.int64); cadj = np.zeros(12, np.int64)
hp = np.zeros((1 << n, n), np.int64); sp = np.zeros((1 << n, n), np.int64); Hs = np.zeros(1 << n, np.int64)
st = np.zeros(4, np.int64)
p9 = A.perm_array(9)
bad = 0
ns1 = 600
ns2 = 60
for i in range(ns1):
    code = int(rng.integers(0, 1 << 15)); y = int(rng.integers(0, n))
    A.adj_from_code(code, n, adj); A.tables(adj, n, hp, sp, Hs); A.root_stats(n, hp, sp, Hs, y, st)
    D = int(st[0]); U = int(st[2] + st[3])
    A.tail_adj(adj, n, y, 1, cadj)
    if A.literal_H(cadj, 8, perms[8]) != D + 2 * U:
        bad += 1
    if i < ns2:
        A.tail_adj(adj, n, y, 2, cadj)
        if A.literal_H(cadj, 9, p9) != D + 4 * U:
            bad += 1
print("order 6 literal tail sample: %d (t=1), %d (t=2) rooted tournaments, failures %d" % (ns1, ns2, bad))

# summary of small pair sets
for n in range(1, 6):
    print("order %d pairs (D,U):" % n, sorted(out_pairs[n]))
print("order 6: offsets", sorted({d for d, u in out_pairs[6]}))
# heads of 11, 13, 47, 67 at a non-sink root
want = {(5, 1, 7): 11, (5, 1, 9): 13, (17, 1, 31): 47, (33, 7, 41): 67}
for n in range(1, NMAX + 1):
    for tr in sorted(out_triples[n]):
        if tr in want:
            print("order %d has triple (D,E,U) = %s : head value D-E+U = %d" % (n, tr, tr[0] - tr[1] + tr[2]))
np.save("pairs_le6.npy", np.array(sorted({(n, d, u) for n in out_pairs for d, u in out_pairs[n]}), dtype=np.int64))
tri = np.array(sorted({(n,) + t for n in out_triples for t in out_triples[n]}), dtype=np.int64)
np.save("triples_le6.npy", tri)
print("elapsed %.1fs" % (time.time() - t0))
