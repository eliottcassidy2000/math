#!/usr/bin/env python3
"""
audit_analysis.py -- S72 audit of THM-4863 (independent code; reads only the auditor's own C dumps in DATA).

Checks:
  T  Theorem H on every merged-deck collision and ordinary same-deck pair (orders 3-9): OCF re-derivation by subset DP,
     equality of the spanning cycle-cover counts with c >= 2 cycles, equality of the non-spanning part, and the
     determinant analogue det(T) - det(T') = (-1)^(n+1) (hc(T) - hc(T')) (Pouzet 1979 / Bankoussou-mabiala et al. 2019)
  L  Lemma L and the classification of merged-deck collisions
  J  Lemma J on every join A => B of ordinary classes with |A| + |B| <= 8
  F  Fibonacci-valued face rows at orders 4-8 and their join structure
  P  PRIMEGAME calibration (fractions as moves between order-6 node codes)
  M  misc: tiling sums, S_4 example, order-7 pure-blue node 1 -> T5 -> 1
"""
import sys, random, itertools
from math import gcd
from functools import reduce
from fractions import Fraction

DATA = sys.argv[1] if len(sys.argv) > 1 else "data"
NMAX = 9

# ---------------------------------------------------------------- canonical form (port of audit_tower.c)
def popc(x):
    return bin(x).count("1")

def refine(n, adj, col, ncol):
    while True:
        k = ncol
        sig = []
        for v in range(n):
            s = [col[v]] + [0] * k
            m = adj[v]
            while m:
                u = (m & -m).bit_length() - 1
                m &= m - 1
                s[col[u] + 1] += 1
            sig.append(tuple(s))
        order = sorted(range(n), key=lambda v: sig[v])
        newcol = [0] * n
        nc = 0
        for i, v in enumerate(order):
            if i > 0 and sig[order[i]] != sig[order[i - 1]]:
                nc += 1
            newcol[v] = nc
        nc += 1
        col[:] = newcol
        if nc == k:
            return col, ncol
        ncol = nc

def canon(n, adj):
    score = [popc(adj[v]) for v in range(n)]
    key = []
    for v in range(n):
        inn = 0
        for u in range(n):
            if u != v and not (adj[v] >> u) & 1:
                inn |= 1 << u
        c3 = 0
        so = 0
        for u in range(n):
            if (adj[v] >> u) & 1:
                c3 += popc(adj[u] & inn)
                so += score[u]
        key.append(score[v] * 10000 + c3 * 100 + so)
    d = sorted(set(key))
    col = [d.index(key[v]) for v in range(n)]
    col, ncol = refine(n, adj, col, len(d))
    best = [None, 0]

    def search(col, ncol):
        if ncol == n:
            at = [0] * n
            for v in range(n):
                at[col[v]] = v
            code = 0
            for j in range(1, n):
                for i in range(j):
                    code = (code << 1) | ((adj[at[j]] >> at[i]) & 1)
            if best[0] is None or code > best[0]:
                best[0] = code
                best[1] = 1
            elif code == best[0]:
                best[1] += 1
            return
        cnt = [0] * n
        for v in range(n):
            cnt[col[v]] += 1
        c = 0
        while cnt[c] == 1:
            c += 1
        for v in range(n):
            if col[v] != c:
                continue
            col2 = [(x + 1 if x > c else (c + 1 if x == c else x)) for x in col]
            col2[v] = c
            c2, nc2 = refine(n, adj, col2, ncol + 1)
            search(c2, nc2)

    search(list(col), ncol)
    return best[0], best[1]

def decode(n, code):
    adj = [0] * n
    b = n * (n - 1) // 2 - 1
    for j in range(1, n):
        for i in range(j):
            if (code >> b) & 1:
                adj[j] |= 1 << i
            else:
                adj[i] |= 1 << j
            b -= 1
    return adj

def converse(n, adj):
    op = [0] * n
    for i in range(n):
        for j in range(n):
            if (adj[i] >> j) & 1:
                op[j] |= 1 << i
    return op

def delv(n, adj, x):
    res = []
    for i in range(n):
        if i == x:
            continue
        a = adj[i]
        lo = a & ((1 << x) - 1)
        hi = (a >> (x + 1)) << x
        res.append(lo | hi)
    return res

def hampaths(n, adj):
    full = (1 << n) - 1
    def ends(g):
        dp = [[0] * n for _ in range(full + 1)]
        for v in range(n):
            dp[1 << v][v] = 1
        for S in range(1, full + 1):
            row = dp[S]
            for v in range(n):
                c = row[v]
                if not c:
                    continue
                m = g[v] & ~S
                while m:
                    w = (m & -m).bit_length() - 1
                    m &= m - 1
                    dp[S | (1 << w)][w] += c
        return [dp[full][v] for v in range(n)]
    he = ends(adj)
    hs = ends(converse(n, adj))
    return hs, he

def join(A, a, B, b):
    """A => B: vertices 0..a-1 = A, a..a+b-1 = B, every A vertex beats every B vertex"""
    n = a + b
    out = [0] * n
    for i in range(a):
        out[i] = A[i] | (((1 << b) - 1) << a)
    for j in range(b):
        out[a + j] = B[j] << a
    return out

# ---------------------------------------------------------------- load the auditor's dumps
classes = {}   # n -> list of dict
nodes = {}
code2idx = {}
for n in range(1, NMAX + 1):
    cl = []
    for line in open(f"{DATA}/audit_classes_{n}.txt"):
        if line.startswith("#"):
            continue
        a, b = line.split("|")
        f = a.split()
        cl.append(dict(idx=int(f[0]), code=int(f[1]), conv=int(f[2]), aut=int(f[3]), mid=int(f[4]), H=int(f[5]), hc=int(f[6]),
                       odeck=[int(x) for x in b.split()]))
    classes[n] = cl
    code2idx[n] = {c['code']: c['idx'] for c in cl}
    if n >= 2:
        nd = []
        for line in open(f"{DATA}/audit_nodes_{n}.txt"):
            if line.startswith("#"):
                continue
            a, W, md = line.split("|")
            f = a.split()
            Wd = {int(k): int(v) for k, v in (t.split(":") for t in W.split())}
            nd.append(dict(mid=int(f[0]), rep=int(f[1]), code=int(f[2]), sc=int(f[3]), aut=int(f[4]), H=int(f[5]), hc=int(f[6]),
                           W=Wd, mdeck=[int(x) for x in md.split()]))
        for x in nd:
            den = 2 * x['aut'] if x['sc'] else x['aut']
            assert all(v % den == 0 for v in x['W'].values())
            x['R'] = {k: v // den for k, v in x['W'].items()}
            x['til'] = x['H'] // x['aut'] if x['sc'] else 2 * x['H'] // x['aut']
        nodes[n] = nd

def merged_id(n, adj):
    c, _ = canon(n, adj)
    return classes[n][code2idx[n][c]]['mid']

def class_idx(n, adj):
    c, _ = canon(n, adj)
    return code2idx[n][c]

print("=" * 110)
print("0. Port check: Python canon reproduces the C canonical codes and |Aut| (orders <= 7 all classes, order 8-9 samples)")
for n in range(2, NMAX + 1):
    cl = classes[n]
    sample = cl if len(cl) <= 500 else random.Random(n).sample(cl, 300)
    bad = 0
    for c in sample:
        adj = decode(n, c['code'])
        perm = list(range(n))
        random.Random(c['code']).shuffle(perm)
        r = [0] * n
        for a in range(n):
            for b in range(n):
                if (adj[a] >> b) & 1:
                    r[perm[a]] |= 1 << perm[b]
        cc, aut = canon(n, r)
        if cc != c['code'] or aut != c['aut']:
            bad += 1
    print(f"  n={n}: {len(sample)} classes checked, {bad} mismatches")

print("=" * 110)
print("M. Tiling sums sum_u |F(u)| = 2^binom(n-1,2); S_4 example")
for n in range(2, NMAX + 1):
    tot = sum(x['til'] for x in nodes[n])
    print(f"  n={n}: sum of node tilings = {tot}, 2^m = {2 ** ((n - 1) * (n - 2) // 2)}")
S4 = [x for x in nodes[4] if x['H'] == 5][0]
adj = decode(4, S4['code'])
hs, he = hampaths(4, adj)
cards = [merged_id(3, delv(4, adj, x)) for x in range(4)]
print(f"  S_4: h_s={hs} h_e={he} h_s+h_e={[a + b for a, b in zip(hs, he)]} cards(merged id, 0/1 = H of order-3 node {[x['H'] for x in nodes[3]]})={cards} W={S4['W']} R={S4['R']}")

# ---------------------------------------------------------------- T. Theorem H
print("=" * 110)
print("T. Theorem H: decomposition of H(T) by the odd-cycle formula on every collision pair")

def subset_hc(n, adj):
    """hc[U] = number of directed Hamiltonian cycles of T[U] (|U| >= 3), by DP rooted at min(U)"""
    hc = [0] * (1 << n)
    for U in range(1, 1 << n):
        if popc(U) < 3:
            continue
        r = (U & -U).bit_length() - 1
        verts = [v for v in range(n) if (U >> v) & 1]
        dp = {(1 << r, r): 1}
        # iterate over subsets of U containing r in increasing popcount
        layer = {(1 << r, r): 1}
        for _ in range(len(verts) - 1):
            nl = {}
            for (S, v), c in layer.items():
                m = adj[v] & U & ~S
                while m:
                    w = (m & -m).bit_length() - 1
                    m &= m - 1
                    key = (S | (1 << w), w)
                    nl[key] = nl.get(key, 0) + c
            layer = nl
        tot = 0
        for (S, v), c in layer.items():
            if (adj[v] >> r) & 1:
                tot += c
        hc[U] = tot
    return hc

def cover_poly(n, adj):
    """G[S] = polynomial (list) in y: sum over spanning collections of vertex-disjoint odd directed cycles of T[S], y^(#cycles)"""
    hc = subset_hc(n, adj)
    G = [None] * (1 << n)
    G[0] = [1]
    for S in range(1, 1 << n):
        r = (S & -S).bit_length() - 1
        rest = S & ~(1 << r)
        poly = [0]
        # U = {r} | sub, sub subset of rest, |U| odd >= 3
        sub = rest
        while True:
            U = sub | (1 << r)
            if popc(U) >= 3 and popc(U) % 2 == 1 and hc[U]:
                g2 = G[S & ~U]
                if g2 and any(g2):
                    need = len(g2) + 1
                    if len(poly) < need:
                        poly += [0] * (need - len(poly))
                    for i, c in enumerate(g2):
                        poly[i + 1] += hc[U] * c
            if sub == 0:
                break
            sub = (sub - 1) & rest
        G[S] = poly
    return G, hc

def det_int(M):
    """exact integer determinant (Bareiss)"""
    M = [row[:] for row in M]
    n = len(M)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = None
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    sw = i
                    break
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]

def adjmat(n, adj):
    return [[(adj[i] >> j) & 1 for j in range(n)] for i in range(n)]

def H_profile(n, code):
    adj = decode(n, code)
    G, hc = cover_poly(n, adj)
    full = (1 << n) - 1
    H_ocf = sum(sum(c * 2 ** i for i, c in enumerate(G[S])) for S in range(1 << n))
    span = G[full] + [0] * 6
    nonspan = H_ocf - sum(c * 2 ** i for i, c in enumerate(G[full]))
    hs, he = hampaths(n, adj)
    return dict(H=sum(he), H_ocf=H_ocf, hc=hc[full] if n >= 3 else 0, span=span[:6], nonspan=nonspan, det=det_int(adjmat(n, adj)))

def check_pair(n, c1, c2, label):
    p1, p2 = H_profile(n, c1), H_profile(n, c2)
    ocf_ok = p1['H'] == p1['H_ocf'] and p2['H'] == p2['H_ocf']
    span2_equal = p1['span'][2:] == p2['span'][2:]
    nons_equal = p1['nonspan'] == p2['nonspan']
    dH, dhc = p1['H'] - p2['H'], p1['hc'] - p2['hc']
    thmH = (dH == 0) if n % 2 == 0 else (dH == 2 * dhc)
    ddet = p1['det'] - p2['det']
    detlaw = ddet == (-1) ** (n + 1) * dhc
    print(f"  n={n} {label}: H {p1['H']}/{p2['H']} (OCF ok {ocf_ok}), hc {p1['hc']}/{p2['hc']}, spanning covers by #cycles {p1['span'][1:5]} / {p2['span'][1:5]},"
          f" c>=2 equal {span2_equal}, non-spanning part {p1['nonspan']}/{p2['nonspan']} equal {nons_equal}; Theorem H {'OK' if thmH else 'FAIL'};"
          f" det {p1['det']}/{p2['det']}, det law {'OK' if detlaw else 'FAIL'}")
    return thmH and span2_equal and nons_equal and ocf_ok

# merged-deck collisions
mcoll = {}
for n in range(3, NMAX + 1):
    groups = {}
    for x in nodes[n]:
        groups.setdefault(tuple(x['mdeck']), []).append(x)
    mcoll[n] = [g for g in groups.values() if len(g) > 1]
# ordinary same-deck classes
ocoll = {}
for n in range(3, NMAX + 1):
    groups = {}
    for c in classes[n]:
        groups.setdefault(tuple(c['odeck']), []).append(c)
    ocoll[n] = [g for g in groups.values() if len(g) > 1]
allok = True
for n in range(3, NMAX + 1):
    for g in mcoll[n]:
        for a, b in itertools.combinations(g, 2):
            allok &= check_pair(n, a['code'], b['code'], "merged-deck collision")
    for g in ocoll[n]:
        for a, b in itertools.combinations(g, 2):
            allok &= check_pair(n, a['code'], b['code'], "ordinary same-deck pair")
print(f"  all decomposition checks passed: {allok}")

# OCF sanity on all classes n <= 6 and a sample at 7
bad = 0
cnt = 0
for n in range(3, 8):
    pool = classes[n] if n <= 6 else random.Random(7).sample(classes[n], 60)
    for c in pool:
        adj = decode(n, c['code'])
        G, _ = cover_poly(n, adj)
        H_ocf = sum(sum(cc * 2 ** i for i, cc in enumerate(G[S])) for S in range(1 << n))
        cnt += 1
        bad += H_ocf != c['H']
print(f"  OCF H = sum_S G_S(2) re-derived on {cnt} classes (all n <= 6, 60 at n = 7): {bad} failures")

# ---------------------------------------------------------------- L. Lemma L and the classification
print("=" * 110)
print("L. Lemma L (source lift) and the classification of merged-deck collisions")
def strongly_connected(n, adj):
    # reachability from 0 in adj and in converse
    def reach(g):
        seen = 1
        stack = [0]
        while stack:
            v = stack.pop()
            m = g[v] & ~seen
            while m:
                w = (m & -m).bit_length() - 1
                m &= m - 1
                seen |= 1 << w
                stack.append(w)
        return seen == (1 << n) - 1
    return reach(adj) and reach(converse(n, adj))

for m in (2, 4, 6, 8):
    pairs = [g for g in ocoll.get(m, [])] if m >= 3 else []
    for g in pairs:
        a, b = g
        conv_pair = a['conv'] == b['idx']
        S = decode(m, a['code'])
        Sop = decode(m, b['code'])
        L1 = join([0], 1, S, m)
        L2 = join([0], 1, Sop, m)
        u1, u2 = merged_id(m + 1, L1), merged_id(m + 1, L2)
        d1 = nodes[m + 1][u1]['mdeck']
        d2 = nodes[m + 1][u2]['mdeck']
        print(f"  order-{m} same-deck pair (H={a['H']}, converse pair: {conv_pair}): lifts 1->S, 1->S' are merged nodes {u1}, {u2} at order {m + 1};"
              f" distinct {u1 != u2}; equal merged decks {d1 == d2}")
for n in (3, 5, 7, 9):
    for g in mcoll[n]:
        hcs = [x['hc'] for x in g]
        desc = []
        for x in g:
            adj = decode(n, x['code'])
            # source / sink
            src = [v for v in range(n) if popc(adj[v]) == n - 1]
            snk = [v for v in range(n) if popc(adj[v]) == 0]
            desc.append(f"node {x['mid']} H={x['H']} hc={x['hc']} sc={x['sc']} strong={strongly_connected(n, adj)} source={bool(src)} sink={bool(snk)}")
        # is it a same-ordinary-deck pair?
        same_ord = False
        for og in ocoll.get(n, []):
            mids = sorted(c['mid'] for c in og)
            if mids == sorted(x['mid'] for x in g):
                same_ord = True
        lift = None
        if all(h == 0 for h in hcs):
            rems = []
            for x in g:
                adj = decode(n, x['code'])
                src = [v for v in range(n) if popc(adj[v]) == n - 1]
                if src:
                    rem = delv(n, adj, src[0])
                else:
                    snk = [v for v in range(n) if popc(adj[v]) == 0]
                    rem = converse(n - 1, delv(n, adj, snk[0]))   # S -> 1 is the converse of 1 -> S^op
                rems.append(class_idx(n - 1, rem))
            pair_ok = any(sorted(rems) == sorted(c['idx'] for c in og) and og[0]['conv'] == og[1]['idx'] for og in ocoll.get(n - 1, []))
            lift = f"source-lift remainders {rems} = an order-{n - 1} converse same-deck pair: {pair_ok}"
        print(f"  n={n} collision: " + " || ".join(desc) + f" | same ordinary deck: {same_ord}" + (f" | {lift}" if lift else ""))

# ---------------------------------------------------------------- J. Lemma J on all joins |A| + |B| <= 8
print("=" * 110)
print("J. Lemma J (join rule) on every join A => B of ordinary classes with |A| + |B| <= 8")
def W_direct(n, adj):
    hs, he = hampaths(n, adj)
    W = {}
    for x in range(n):
        v = merged_id(n - 1, delv(n, adj, x))
        W[v] = W.get(v, 0) + hs[x] + he[x]
    return {k: v for k, v in W.items() if v}

def W_join(A, a, B, b):
    n = a + b
    if a > 1:
        hsA, heA = hampaths(a, A)
    else:
        hsA, heA = [1], [1]
    if b > 1:
        hsB, heB = hampaths(b, B)
    else:
        hsB, heB = [1], [1]
    HA, HB = sum(heA), sum(heB)
    W = {}
    for x in range(a):
        sub = join(delv(a, A, x), a - 1, B, b) if a > 1 else B
        v = merged_id(n - 1, sub)
        W[v] = W.get(v, 0) + HB * hsA[x]
    for y in range(b):
        sub = join(A, a, delv(b, B, y), b - 1) if b > 1 else A
        v = merged_id(n - 1, sub)
        W[v] = W.get(v, 0) + HA * heB[y]
    return {k: v for k, v in W.items() if v}

checked = fails = 0
for n in range(3, 9):
    for a in range(1, n):
        b = n - a
        for ca in classes[a]:
            A = decode(a, ca['code']) if a > 1 else [0]
            for cb in classes[b]:
                B = decode(b, cb['code']) if b > 1 else [0]
                t = join(A, a, B, b)
                checked += 1
                if W_direct(n, t) != W_join(A, a, B, b):
                    fails += 1
print(f"  joins checked: {checked}, failures: {fails}")

# ---------------------------------------------------------------- F. Fibonacci rows
print("=" * 110)
print("F. Face rows that are runs of consecutive Fibonacci numbers; strong-component structure of their nodes")
fib = [1, 2]
while fib[-1] < 10 ** 9:
    fib.append(fib[-1] + fib[-2])
fidx = {v: i for i, v in enumerate(fib)}

def components(n, adj):
    """strong components in the order of the condensation (first component beats all later ones)"""
    # reachability matrix
    reach = []
    for v in range(n):
        seen = 1 << v
        stack = [v]
        while stack:
            u = stack.pop()
            m = adj[u] & ~seen
            while m:
                w = (m & -m).bit_length() - 1
                m &= m - 1
                seen |= 1 << w
                stack.append(w)
        reach.append(seen)
    comps = []
    done = 0
    for v in range(n):
        if (done >> v) & 1:
            continue
        comp = [u for u in range(n) if (reach[v] >> u) & 1 and (reach[u] >> v) & 1]
        for u in comp:
            done |= 1 << u
        comps.append(comp)
    # order: component X before Y if X beats Y
    comps.sort(key=lambda c: -popc(adj[c[0]] & ~sum(1 << u for u in c)))
    names = []
    for c in comps:
        k = len(c)
        if k == 1:
            names.append("1")
        else:
            sub = [0] * k
            for i, u in enumerate(c):
                for j, w in enumerate(c):
                    if (adj[u] >> w) & 1:
                        sub[i] |= 1 << j
            hs, he = hampaths(k, sub)
            names.append({3: "C3", 4: "S4"}.get(k, f"strong{k}") + ("" if k <= 4 else f"(H={sum(he)})"))
    return names

for n in range(4, 9):
    multi = [x for x in nodes[n] if len(x['R']) >= 2]
    hits = []
    for x in multi:
        vals = sorted(x['R'].values())
        if len(set(vals)) == len(vals) and all(v in fidx for v in vals):
            ids = [fidx[v] for v in vals]
            if all(ids[i + 1] - ids[i] == 1 for i in range(len(ids) - 1)):
                hits.append(x)
    print(f"  order {n}: {len(hits)} of {len(multi)} multi-entry rows are consecutive-Fibonacci runs")
    for x in hits:
        adj = decode(n, x['code'])
        names = components(n, adj)
        print(f"     node {x['mid']} H={x['H']} aut={x['aut']} {'sc' if x['sc'] else 'pair'} row values {sorted(x['R'].values())}; strong components (top to bottom): {' => '.join(names)}")

# source extensions of the order-5 rows (5,8) and (2,3)
for x in nodes[5]:
    if sorted(x['R'].values()) in ([2, 3], [5, 8]):
        adj = decode(5, x['code'])
        t = join([0], 1, adj, 5)
        u = merged_id(6, t)
        print(f"  order-5 node {x['mid']} (H={x['H']}, row {sorted(x['R'].values())}) -> source extension = order-6 node {u}, row values {sorted(nodes[6][u]['R'].values())}")

# ---------------------------------------------------------------- order-7 pure-blue node with 3 tilings
print("=" * 110)
print("M2. Order-7 node with H=15, |Aut|=5, self-converse: is it 1 -> T5 -> 1 ?")
T5 = [0] * 5
for i in range(5):
    for d in (1, 2):
        T5[i] |= 1 << ((i + d) % 5)
t = join(join([0], 1, T5, 5), 6, [0], 1)
u = merged_id(7, t)
cand = [x for x in nodes[7] if x['H'] == 15 and x['aut'] == 5 and x['sc']]
print(f"  1 -> T5 -> 1 is merged node {u} (H={nodes[7][u]['H']}, aut={nodes[7][u]['aut']}, tilings={nodes[7][u]['til']}); candidates with H=15, aut=5, sc: {[c['mid'] for c in cand]}")

# ---------------------------------------------------------------- P. PRIMEGAME calibration
print("=" * 110)
print("P. PRIMEGAME fractions as moves between order-6 node codes")
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
fracs = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17), (11, 13), (13, 11), (15, 2), (1, 7), (55, 1)]
def val(x):
    v = [0] * 10
    for i, p in enumerate(primes):
        while x % p == 0:
            x //= p
            v[i] += 1
    assert x == 1
    return v
fv = {}
for p, q in fracs:
    d = tuple(a - b for a, b in zip(val(p), val(q)))
    fv[d] = fv.get(d, 0) + 1
low = nodes[5]
order5 = sorted(low, key=lambda x: (x['H'], x['til'], sorted(x['R'].items())))
canon_assign = {x['mid']: i for i, x in enumerate(order5)}   # mid -> prime index
print("  canonical assignment (by H, tilings, row): " + ", ".join(f"{primes[canon_assign[x['mid']]]}<-node{x['mid']}(H={x['H']},til={x['til']},R={sorted(x['R'].values())})" for x in order5))

def codes_and_moves(assign, use_prim=False):
    vecs = []
    for x in nodes[6]:
        v = [0] * 10
        R = x['R']
        g = reduce(gcd, R.values()) if use_prim else 1
        for k, e in R.items():
            v[assign[k]] += e // g
        vecs.append(tuple(v))
    hits = 0
    for a in range(len(vecs)):
        for b in range(len(vecs)):
            if a != b:
                d = tuple(x - y for x, y in zip(vecs[b], vecs[a]))
                hits += fv.get(d, 0)
    return vecs, hits

vecs, h0 = codes_and_moves(canon_assign)
ints = []
for v in vecs:
    N = 1
    for i, e in enumerate(v):
        N *= primes[i] ** e
    ints.append(N)
print(f"  distinct integer codes: {len(set(ints))}/34; examples: TT6 -> {ints[[i for i, x in enumerate(nodes[6]) if x['H'] == 1][0]]}")
for x, N in zip(nodes[6], ints):
    if N < 10 ** 6:
        adj = decode(6, x['code'])
        print(f"     order-6 node {x['mid']} H={x['H']} {'sc' if x['sc'] else 'pair'} components {' => '.join(components(6, adj))}: N = {N}")
print(f"  canonical assignment: {h0} (pair, fraction) incidences")
for use_prim, label in ((False, "raw rows"), (True, "primitive rows")):
    _, hc0 = codes_and_moves(canon_assign, use_prim)
    rng = random.Random(4242)
    ids = [x['mid'] for x in low]
    samp = []
    for _ in range(20000):
        perm = list(range(10))
        rng.shuffle(perm)
        samp.append(codes_and_moves({k: perm[i] for i, k in enumerate(ids)}, use_prim)[1])
    samp.sort()
    ge = sum(1 for s in samp if s >= hc0)
    print(f"  {label}: canonical {hc0}; random (20,000, my seed): min {samp[0]}, median {samp[len(samp) // 2]}, mean {sum(samp) / len(samp):.3f}, max {samp[-1]}; P(random >= canonical) = {ge / len(samp):.4f}")
