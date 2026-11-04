#!/usr/bin/env python3
"""audit_B_statement7_census.py -- audit B of THM-4871 statement 7: independent census of the N-even H-fixed graphs.

Orbit method (own code, brute-force canonical forms for the 4-coloured objects on the k <= 4 fixed points):
  for every N-even N-orbit representative (C_P, C_Q, cross, gamma) build X, compute Aut X (pynauty), its sign
  (sign-even iff all generators even), and the type certificate; report
    tau_E, W (sign-odd types), E_even, E_wit (multiplicities among N-even representatives), and for tournaments
    tau_T, E_T;  defect = tau_T - tau_E;  and the identity W + E_even + E_wit - E_T.
  Statement 7 predicts for primes p < q < 2p:  k < q-p: W = E_even = E_wit = E_T = 0;  k = q-p: W = 1 (the witness
  is K_{q,q} with parts Q and P u F), E_even = E_wit = E_T = 0.
  The proof of 7(a) (Steps 1-4) in fact only uses k < min(p, q-p): tested on (3,7) k <= 2 and the NEW family
  (5,11) k <= 4 (q > 2p, never computed before), and statement 7 itself on the NEW family (13,17) k <= 4.

Definition-level brute force (--brute): all 2^c H-fixed graphs and tournaments; N enumerated by closure; N-evenness
  computed directly (every g in N fixing X has even reversal parity), types by certificate, N-orbits by union-find.
Usage: python3 audit_B_statement7_census.py p q kmax [--brute kb]
"""
import sys, time, itertools, math
from collections import defaultdict, Counter
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from auditb_common import *
import numpy as np

def orbit_census(p, q, k):
    n = p + q + k
    skewP = dilation_classes(skew_sets(p), p); skewQ = dilation_classes(skew_sets(q), q)
    symP = [C for C in dilation_classes(symmetric_sets(p), p) if multiplier_even(C, p)]
    symQ = [C for C in dilation_classes(symmetric_sets(q), q) if multiplier_even(C, q)]
    colT = coloured_objects(k, True)
    colG_all = coloured_objects(k, False)
    colG = [(E, col, ev) for (E, col, ev) in colG_all if ev]
    bk = burnside_b_k4(k)
    assert len(colT) == bk, (len(colT), bk)
    assert len(colG) == bk, (len(colG), bk)   # ec_k(4) = b_k(4), THM-4721 Theorem R
    # tournaments
    typesT = Counter(); R = 0
    for SP in skewP:
        for SQ in skewQ:
            for cross in (0, 1):
                for (AF, col, _) in colT:
                    typesT[tournament_cert(n, build_tournament(p, q, k, SP, SQ, cross, AF, col))] += 1; R += 1
    tau_T = len(typesT); E_T = R - tau_T
    # graphs
    typesG = {}; RG = 0; witnesses = []
    for CP in symP:
        for CQ in symQ:
            for cross in (0, 1):
                for (EF, col, _) in colG:
                    edges = build_graph(p, q, k, CP, CQ, cross, EF, col)
                    cert, gens, order = graph_cert_and_gens(n, edges)
                    RG += 1
                    if cert in typesG:
                        typesG[cert][0] += 1; continue
                    ev = sign_even(gens, edges)
                    typesG[cert] = [1, ev]
                    if not ev:
                        witnesses.append((CP, CQ, cross, EF, col, len(edges), order))
    tau_E = sum(1 for m, ev in typesG.values() if ev)
    E_even = sum(m - 1 for m, ev in typesG.values() if ev)
    W = sum(1 for m, ev in typesG.values() if not ev)
    E_wit = sum(m - 1 for m, ev in typesG.values() if not ev)
    defect = tau_T - tau_E
    print("p=%d q=%d k=%d n=%d R=%d RG=%d tau_T=%d tau_E=%d defect=%d W=%d E_T=%d E_even=%d E_wit=%d identity=%s theoremR=%s  t=(%d,%d) e=(%d,%d) b_k=%d"
          % (p, q, k, n, R, RG, tau_T, tau_E, defect, W, E_T, E_even, E_wit, defect == W + E_even + E_wit - E_T, R == RG,
             len(skewP), len(skewQ), len(symP), len(symQ), bk), flush=True)
    for (CP, CQ, cross, EF, col, ne, order) in witnesses:
        isKqq = (len(CP) == 0 and len(CQ) == 0 and cross == 1 and len(EF) == 0 and all(c == 2 for c in col) and ne == q * q)
        print("   witness: C_P=%s C_Q=%s cross=%d E_F=%s col=%s edges=%d |Aut|=%d  -> K_{q,q} on Q vs P u F: %s"
              % (sorted(CP), sorted(CQ), cross, sorted(EF), col, ne, order, isKqq))
    return dict(p=p, q=q, k=k, R=R, RG=RG, tau_T=tau_T, tau_E=tau_E, defect=defect, W=W, E_T=E_T, E_even=E_even, E_wit=E_wit)

# ---------------------------------------------------------------- definition-level brute force
def closure(gens, n):
    ident = tuple(range(n)); seen = {ident}; frontier = [ident]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                x = tuple(h[g[i]] for i in range(n))
                if x not in seen:
                    seen.add(x); nxt.append(x)
        frontier = nxt
    return sorted(seen)

def brute_census(p, q, k):
    t0 = time.time()
    n = p + q + k
    P = list(range(p)); Q = list(range(p, p + q)); F = list(range(p + q, n))
    h = list(range(n))
    for i in range(p): h[i] = (i + 1) % p
    for j in range(q): h[p + j] = p + (j + 1) % q
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    orbits = []; orb_of = {}
    for pr in pairs:
        if pr in orb_of: continue
        orb = []; x = pr
        while x not in orb_of:
            orb_of[x] = len(orbits); orb.append(x); x = tuple(sorted((h[x[0]], h[x[1]])))
        orbits.append(orb)
    c = len(orbits)
    # N by closure: h, dilations on P and Q, Sym(F)
    gens = [tuple(h)]
    for (blk, m) in ((P, p), (Q, q)):
        for a in units(m):
            if len({pow(a, e, m) for e in range(1, m)}) != m - 1:
                continue   # keep only primitive roots (they generate the dilation group)
            g = list(range(n))
            for i in range(m): g[blk[i]] = blk[(a * i) % m]
            gens.append(tuple(g))
    if k >= 2:
        g = list(range(n)); g[F[0]], g[F[1]] = F[1], F[0]; gens.append(tuple(g))
        g = list(range(n))
        for i in range(k): g[F[i]] = F[(i + 1) % k]
        gens.append(tuple(g))
    N = closure(gens, n)
    assert len(N) == p * (p - 1) * q * (q - 1) * math.factorial(k)
    # every g in N normalises <h>: g h g^-1 is a power of h
    powers = set()
    x = list(range(n))
    for _ in range(p * q):
        x = [h[i] for i in x]; powers.add(tuple(x))
    for g in N:
        ginv = [0] * n
        for i in range(n): ginv[g[i]] = i
        conj = tuple(g[h[ginv[i]]] for i in range(n))
        assert conj in powers
    # induced action on orbit indices and parity vectors
    Xs = np.arange(1 << c, dtype=np.int64)
    bits = [((Xs >> o) & 1) for o in range(c)]
    # parity-of-popcount table
    par_tab = np.zeros(1 << c, dtype=np.int8)
    for o in range(c):
        par_tab ^= ((Xs >> o) & 1).astype(np.int8)
    Neven = np.ones(1 << c, dtype=bool)
    parent = np.arange(1 << c, dtype=np.int64)        # union-find for N-orbits on graphs
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    # tournaments: orientation words; orbit o oriented 'forward' = H-orbit of the arc (u, v) of its first pair
    # g in N maps the arc-orbit of (u,v) to the arc-orbit of (g u, g v); record whether that is forward or backward
    arc_orbit_fwd = {}
    for o, orb in enumerate(orbits):
        (u, v) = orb[0]; src, dst = u, v
        for _ in orb:
            arc_orbit_fwd[(src, dst)] = (o, 0); arc_orbit_fwd[(dst, src)] = (o, 1); src, dst = h[src], h[dst]
    parentT = np.arange(1 << c, dtype=np.int64)
    def findT(a):
        while parentT[a] != a:
            parentT[a] = parentT[parentT[a]]; a = parentT[a]
        return a
    gens_set = set(gens)
    for g in N:
        pi = [0] * c; par = 0; flipmask = 0; piT = [0] * c
        for o, orb in enumerate(orbits):
            (u, v) = orb[0]
            pi[o] = orb_of[tuple(sorted((g[u], g[v])))]
            rp = sum(1 for (a, b) in orb if g[a] > g[b]) & 1
            if rp: flipmask |= 1 << o
            (o2, flip) = arc_orbit_fwd[(g[u], g[v])]
            piT[o] = (o2, flip)
        gX = np.zeros(1 << c, dtype=np.int64)
        for o in range(c):
            gX |= bits[o] << pi[o]
        fixed = (gX == Xs)
        odd = par_tab[Xs & flipmask].astype(bool)
        Neven &= ~(fixed & odd)
        if g not in gens_set:
            continue
        # N-orbits = connected components under the generators (union-find on X ~ gX, generators only)
        for a, b in zip(Xs[gX != Xs], gX[gX != Xs]):
            ra, rb = find(int(a)), find(int(b))
            if ra != rb: parent[ra] = rb
        # tournaments: word bit o = 1 means the orbit is oriented backward
        gT = np.zeros(1 << c, dtype=np.int64)
        for o in range(c):
            (o2, flip) = piT[o]
            gT |= (bits[o] ^ flip) << o2
        for a, b in zip(Xs[gT != Xs], gT[gT != Xs]):
            ra, rb = findT(int(a)), findT(int(b))
            if ra != rb: parentT[ra] = rb
    R = len({find(int(x)) for x in Xs}); RT = len({findT(int(x)) for x in Xs})
    # types and signs
    typesG = {}; evenG = np.zeros(1 << c, dtype=bool); certG = [None] * (1 << c)
    for X in range(1 << c):
        edges = sorted(e for o in range(c) if (X >> o) & 1 for e in orbits[o])
        cert, gens_, order = graph_cert_and_gens(n, edges)
        certG[X] = cert; evenG[X] = sign_even(gens_, edges)
    certT = [None] * (1 << c)
    for X in range(1 << c):
        arcs = []
        for o, orb in enumerate(orbits):
            (u, v) = orb[0]
            src, dst = ((u, v) if not (X >> o) & 1 else (v, u))
            for _ in orb:
                arcs.append((src, dst)); src, dst = h[src], h[dst]
        certT[X] = tournament_cert(n, arcs)
    assert all(Neven[X] for X in range(1 << c) if evenG[X])   # sign-even => N-even
    tau_T = len(set(certT))
    tau_E = len({certG[X] for X in range(1 << c) if evenG[X]})
    W = len({certG[X] for X in range(1 << c) if Neven[X] and not evenG[X]})
    # multiplicities: N-orbits per type (N-even orbits for witnesses)
    orbrepG = {}
    for X in range(1 << c):
        orbrepG.setdefault(find(X), X)
    m_even = Counter(); m_wit = Counter()
    for r_, X in orbrepG.items():
        if evenG[X]: m_even[certG[X]] += 1
        elif Neven[X]: m_wit[certG[X]] += 1
    E_even = sum(m - 1 for m in m_even.values()); E_wit = sum(m - 1 for m in m_wit.values())
    E_T = RT - tau_T
    RGeven = int(sum(1 for r_, X in orbrepG.items() if Neven[X]))
    defect = tau_T - tau_E
    print("  BRUTE p=%d q=%d k=%d c=%d |N|=%d R=%d RGeven=%d tau_T=%d tau_E=%d defect=%d W=%d E_T=%d E_even=%d E_wit=%d identity=%s theoremR=%s (%.1fs)"
          % (p, q, k, c, len(N), RT, RGeven, tau_T, tau_E, defect, W, E_T, E_even, E_wit, defect == W + E_even + E_wit - E_T, RT == RGeven, time.time() - t0), flush=True)
    # the witnesses, explicitly
    for X in range(1 << c):
        if Neven[X] and not evenG[X]:
            edges = sorted(e for o in range(c) if (X >> o) & 1 for e in orbits[o])
            deg = Counter();
            for (u, v) in edges: deg[u] += 1; deg[v] += 1
            print("     N-even sign-odd H-fixed graph: %d edges, degree multiset %s, type %s" % (len(edges), sorted(Counter(deg[v] for v in range(n)).items()), certG[X][:8].hex() if isinstance(certG[X], bytes) else "?"))
    return dict(R=RT, RGeven=RGeven, tau_T=tau_T, tau_E=tau_E, defect=defect, W=W, E_T=E_T, E_even=E_even, E_wit=E_wit)

if __name__ == "__main__":
    p, q, kmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    kb = int(sys.argv[sys.argv.index("--brute") + 1]) if "--brute" in sys.argv else -1
    print("audit B: (%d,%d,1^k), k <= %d%s" % (p, q, kmax, ", brute force for k <= %d" % kb if kb >= 0 else ""))
    rows = []
    for k in range(kmax + 1):
        t0 = time.time()
        row = orbit_census(p, q, k)
        rows.append(row)
        if k <= kb:
            b = brute_census(p, q, k)
            agree = all(b[key] == row[key] for key in ("R", "tau_T", "tau_E", "defect", "W", "E_T", "E_even", "E_wit")) and b["RGeven"] == row["RG"]
            print("  brute vs orbit: %s" % ("AGREE" if agree else "DISAGREE"), flush=True)
    print("\nk  R  tau_T tau_E defect W E_T E_even E_wit")
    for r in rows:
        print("%d %d %d %d %d %d %d %d %d" % (r["k"], r["R"], r["tau_T"], r["tau_E"], r["defect"], r["W"], r["E_T"], r["E_even"], r["E_wit"]))
