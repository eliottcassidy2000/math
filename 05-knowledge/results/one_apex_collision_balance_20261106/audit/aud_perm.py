"""Independent audit of THM-4891 (auditor's own code; the session's class/literal encoding is NOT used).

Construction: vertices = permutations of {0..n-1} (arrangements); edge {p, p s_i} for each adjacent transposition.
Each edge e = (u, v) with u < v (vertex indices) carries a variable o_e: o_e = 1 means the edge points INTO v.
Geometry: x(p)[val] = position of val in p.  The direction vector of an oriented edge is x(head) - x(tail).
2-faces are generated directly from the Cayley graph:
  squares  p, p s_i, p s_i s_j, p s_j (|i - j| >= 2), opposite edges = the two s_i-edges and the two s_j-edges;
  hexagons c_0 = p, c_(k+1) = c_k s_i or s_(i+1) alternately (6-cycle), opposite edges E_k, E_(k+3).
"Opposite edges point the same way" = equal direction vectors (the geometric definition, translation of the side).
The coherence relations are o_e = o_e' XOR parity(e, e'), parity from the direction vectors.

Part A: union-find with parity over (squares), (hexagons), (both); class counts vs C(n,2)2^(n-2), C(n,2)(n-2)!, C(n,2);
        consistency; each class described by (pair, set before) resp. (pair, relative order of the others); (i) checked
        for n <= 7.
Part B: exhaustive sink spectra for n = 3, 4 (all class assignments), with acyclicity, via the C helper aud_perm_enum.
Part C: SAT at n = 5 (Glucose4, sequential-counter cardinality, edge variables) for the full sink spectra of both
        relaxations; acyclic witnesses for 7 and 21 by my own CEGAR loop; independent verification of every witness;
        verification of the session's two printed tables (n4_L2_7, n5_L1_7) and realisability of the session's
        printed sink SETS for n5_L1_21, n5_L2_7, n5_L2_21 as exact sink sets of acyclic coherent orientations.
Part D: n = 6 squares: SAT witness with >= 90 sinks (Glucose4 / seqcounter), verified independently.
Usage: python aud_perm.py [A|B|C|D ...]
"""
import itertools, sys, time, subprocess, os, json
from math import comb, factorial

HERE = os.path.dirname(os.path.abspath(__file__))


class Perm:
    def __init__(self, n):
        self.n = n
        self.V = list(itertools.permutations(range(n)))
        self.idx = {p: k for k, p in enumerate(self.V)}
        self.E = []          # (u, v, i) with u < v, the edge swapping positions i, i+1
        self.eid = {}
        for u, p in enumerate(self.V):
            for i in range(n - 1):
                q = list(p); q[i], q[i + 1] = q[i + 1], q[i]
                v = self.idx[tuple(q)]
                key = (min(u, v), max(u, v))
                if key not in self.eid:
                    self.eid[key] = len(self.E)
                    self.E.append((key[0], key[1], i))
        self.inc = [[] for _ in self.V]
        for e, (u, v, i) in enumerate(self.E):
            self.inc[u].append(e); self.inc[v].append(e)

    def swap(self, p, i):
        q = list(p); q[i], q[i + 1] = q[i + 1], q[i]
        return tuple(q)

    def x(self, p):
        r = [0] * self.n
        for pos, val in enumerate(p):
            r[val] = pos
        return tuple(r)

    def dvec(self, e):
        """direction vector of edge e when o_e = 1 (points into the larger index v): x(v) - x(u)."""
        u, v, _ = self.E[e]
        xu, xv = self.x(self.V[u]), self.x(self.V[v])
        return tuple(b - a for a, b in zip(xu, xv))

    def e_of(self, p, q):
        a, b = self.idx[p], self.idx[q]
        return self.eid[(min(a, b), max(a, b))]

    def faces(self, kind):
        """yield opposite edge pairs (e, e') of every square ('sq') or hexagon ('hex')."""
        n = self.n
        seen = set()
        for p in self.V:
            if kind == 'sq':
                for i in range(n - 1):
                    for j in range(i + 2, n - 1):
                        pi, pj = self.swap(p, i), self.swap(p, j)
                        pij = self.swap(pi, j)
                        assert pij == self.swap(pj, i)
                        key = frozenset([p, pi, pj, pij])
                        if key in seen:
                            continue
                        seen.add(key)
                        yield (self.e_of(p, pi), self.e_of(pj, pij))
                        yield (self.e_of(p, pj), self.e_of(pi, pij))
            else:
                for i in range(n - 2):
                    c = [p]
                    for k in range(6):
                        c.append(self.swap(c[-1], i if k % 2 == 0 else i + 1))
                    assert c[6] == p and len(set(c[:6])) == 6
                    key = frozenset(c[:6])
                    if key in seen:
                        continue
                    seen.add(key)
                    Ek = [self.e_of(c[k], c[k + 1]) for k in range(6)]
                    for k in range(3):
                        yield (Ek[k], Ek[k + 3])

    def classes(self, kinds):
        """union-find with parity; returns (comp id per edge, parity per edge, #comps, #contradictions, #nonparallel)."""
        m = len(self.E)
        par = list(range(m)); rel = [0] * m      # rel[e] = parity of e relative to par[e]

        def find(e):
            r = 0
            path = []
            while par[e] != e:
                path.append(e)
                r ^= rel[e]
                e = par[e]
            # path compression
            root = e
            acc = r
            for z in path:
                pz = rel[z]
                par[z] = root
                rel[z] = acc
                acc ^= pz
            return root, r

        contra = nonpar = 0
        D = [self.dvec(e) for e in range(m)]
        for kind in kinds:
            for e, f in self.faces(kind):
                if D[e] == D[f]:
                    want = 0
                elif D[e] == tuple(-z for z in D[f]):
                    want = 1
                else:
                    nonpar += 1
                    continue
                re, pe = find(e); rf, pf = find(f)
                if re == rf:
                    if pe ^ pf != want:
                        contra += 1
                else:
                    par[re] = rf; rel[re] = pe ^ pf ^ want
        comp = {}
        cid = [0] * m; cpar = [0] * m
        for e in range(m):
            r, p_ = find(e)
            if r not in comp:
                comp[r] = len(comp)
            cid[e] = comp[r]; cpar[e] = p_
        return cid, cpar, len(comp), contra, nonpar


def describe_edge(P, e):
    u, v, i = P.E[e]
    p = P.V[u]
    a, b = p[i], p[i + 1]
    pair = (min(a, b), max(a, b))
    before = frozenset(p[:i])
    others = tuple(z for z in p if z not in pair)
    return pair, before, others


def partA():
    out = []
    for n in range(3, 8):
        P = Perm(n)
        res = {}
        for name, kinds in (('squares', ['sq']), ('hexagons', ['hex']), ('both', ['sq', 'hex'])):
            cid, cpar, nc, contra, nonpar = P.classes(kinds)
            # describe classes
            desc = {}
            ok_desc = True
            for e in range(len(P.E)):
                pair, before, others = describe_edge(P, e)
                key = {'squares': (pair, before), 'hexagons': (pair, others), 'both': pair}[name]
                if cid[e] in desc and desc[cid[e]] != key:
                    ok_desc = False
                desc[cid[e]] = key
            injective = len(set(desc.values())) == len(desc)
            # parity = geometric coherence: within a class, all edges have the same direction (o = 1 vector sign-adjusted)
            res[name] = (nc, contra, nonpar, ok_desc and injective)
        exp = (comb(n, 2) * 2 ** (n - 2), comb(n, 2) * factorial(n - 2), comb(n, 2))
        out.append('n=%d: edges %d; classes squares %d (expected C(n,2)2^(n-2) = %d), hexagons %d (expected C(n,2)(n-2)! = %d), '
                   'both %d (expected C(n,2) = %d); contradictions %s; non-parallel opposite pairs %s; '
                   'class = (pair, set before) / (pair, order of others) / pair: %s'
                   % (n, len(P.E), res['squares'][0], exp[0], res['hexagons'][0], exp[1], res['both'][0], exp[2],
                      [res[k][1] for k in res], [res[k][2] for k in res], [res[k][3] for k in res]))
    return out


def write_struct(P, kinds, path):
    cid, cpar, nc, contra, nonpar = P.classes(kinds)
    assert contra == 0 and nonpar == 0
    with open(path, 'w') as f:
        f.write('%d %d %d\n' % (len(P.V), len(P.E), nc))
        for e, (u, v, i) in enumerate(P.E):
            f.write('%d %d %d %d\n' % (u, v, cid[e], cpar[e]))
    return nc


def partB():
    out = []
    exe = os.path.join(HERE, 'aud_perm_enum.exe')
    for n in (3, 4):
        P = Perm(n)
        for name, kinds in (('tournaments (both)', ['sq', 'hex']), ('hexagon-coherent', ['hex']), ('square-coherent', ['sq'])):
            path = os.path.join(HERE, '_struct_n%d_%s.txt' % (n, kinds[0] if len(kinds) == 1 else 'both'))
            nc = write_struct(P, kinds, path)
            r = subprocess.run([exe, path], capture_output=True, text=True)
            out.append('n=%d %-20s classes=%2d  %s' % (n, name, nc, r.stdout.strip().replace('\n', '\n' + ' ' * 36)))
            os.remove(path)
    return out


# ---------------- SAT (Glucose4, sequential counter, edge variables) ----------------

def sat_setup(P, kinds):
    from pysat.formula import IDPool
    cid, cpar, nc, contra, nonpar = P.classes(kinds)
    assert contra == 0 and nonpar == 0
    m = len(P.E)
    pool = IDPool(start_from=m + 1)
    clauses = []
    # coherence: o_e <-> o_f xor parity, via classes: o_e <-> o_rep(e) xor cpar
    rep = {}
    for e in range(m):
        if cid[e] not in rep:
            rep[cid[e]] = e
    for e in range(m):
        r = rep[cid[e]]
        if r == e:
            continue
        x = cpar[e] ^ cpar[r]
        oe, orr = e + 1, r + 1
        if x == 0:
            clauses += [[-oe, orr], [oe, -orr]]
        else:
            clauses += [[-oe, -orr], [oe, orr]]
    svars = []
    for p in range(len(P.V)):
        sv = pool.id(('s', p))
        svars.append(sv)
        lits = []
        for e in P.inc[p]:
            u, v, _ = P.E[e]
            lits.append((e + 1) if v == p else -(e + 1))     # points into p
        for l in lits:
            clauses.append([-sv, l])
        clauses.append([sv] + [-l for l in lits])
    return pool, clauses, svars


def sinks_and_acyclic(P, model_bits, kinds=None):
    """independent verification: build the oriented graph, count sinks, topological sort, coherence on each face type."""
    m = len(P.E)
    out = [[] for _ in P.V]
    head = [0] * m
    for e, (u, v, i) in enumerate(P.E):
        if model_bits[e]:
            out[u].append(v); head[e] = v
        else:
            out[v].append(u); head[e] = u
    sinks = [p for p in range(len(P.V)) if not out[p]]
    indeg = [0] * len(P.V)
    for p in range(len(P.V)):
        for q in out[p]:
            indeg[q] += 1
    st = [p for p in range(len(P.V)) if indeg[p] == 0]
    seen = 0
    while st:
        p = st.pop(); seen += 1
        for q in out[p]:
            indeg[q] -= 1
            if indeg[q] == 0:
                st.append(q)
    acyclic = seen == len(P.V)

    def odir(e):
        u, v, _ = P.E[e]
        tail, hd = (u, v) if model_bits[e] else (v, u)
        return tuple(b - a for a, b in zip(P.x(P.V[tail]), P.x(P.V[hd])))
    coh = {}
    for kind in ('sq', 'hex'):
        coh[kind] = all(odir(e) == odir(f) for e, f in P.faces(kind))
    return sinks, acyclic, coh


def find_cycle(P, model_bits):
    adj = [[] for _ in P.V]
    for e, (u, v, i) in enumerate(P.E):
        if model_bits[e]:
            adj[u].append((v, e))
        else:
            adj[v].append((u, e))
    color = [0] * len(P.V); parent = [None] * len(P.V)
    for s in range(len(P.V)):
        if color[s]:
            continue
        stack = [(s, iter(adj[s]))]; color[s] = 1
        while stack:
            x, it = stack[-1]
            nxt = next(it, None)
            if nxt is None:
                color[x] = 2; stack.pop(); continue
            y, e = nxt
            if color[y] == 0:
                color[y] = 1; parent[y] = (x, e); stack.append((y, iter(adj[y])))
            elif color[y] == 1:
                cyc = [e]; z = x
                while z != y:
                    pz, pe = parent[z]; cyc.append(pe); z = pz
                return cyc
    return None


def solve(P, kinds, k=None, mode='eq', acyclic=False, fixed_sinks=None, solver='g4', enc='seqcounter', max_iter=50000):
    from pysat.card import CardEnc, EncType
    from pysat.solvers import Glucose4, Minisat22, Lingeling
    pool, clauses, svars = sat_setup(P, kinds)
    if fixed_sinks is not None:
        for p in range(len(P.V)):
            clauses.append([svars[p]] if p in fixed_sinks else [-svars[p]])
    elif k is not None:
        et = {'seqcounter': EncType.seqcounter, 'kmtotalizer': EncType.kmtotalizer, 'sortnetwrk': EncType.sortnetwrk}[enc]
        if mode == 'eq':
            card = CardEnc.equals(lits=svars, bound=k, vpool=pool, encoding=et)
        else:
            card = CardEnc.atleast(lits=svars, bound=k, vpool=pool, encoding=et)
        clauses += card.clauses
    S = {'g4': Glucose4, 'm22': Minisat22, 'lgl': Lingeling}[solver](bootstrap_with=clauses)
    m = len(P.E)
    it = 0
    try:
        while True:
            it += 1
            if not S.solve():
                return False, it, None
            model = S.get_model()
            bits = [1 if model[e] > 0 else 0 for e in range(m)]
            if not acyclic:
                return True, it, bits
            cyc = find_cycle(P, bits)
            if cyc is None:
                return True, it, bits
            S.add_clause([-(e + 1) if bits[e] else (e + 1) for e in cyc])
            if it > max_iter:
                return None, it, None
    finally:
        S.delete()


def partC():
    out = []
    P = Perm(5)
    for name, kinds, top in (('square-coherent', ['sq'], 60), ('hexagon-coherent', ['hex'], 60)):
        t0 = time.time()
        yes, no = [], []
        for k in range(0, top + 1):
            ok, it, bits = solve(P, kinds, k)
            (yes if ok else no).append(k)
            if ok:
                sinks, ac, coh = sinks_and_acyclic(P, bits)
                assert len(sinks) == k and coh[kinds[0]], (name, k)
        out.append('n=5 %s (Glucose4, seqcounter, edge variables): achievable %s; not achievable %s [%.1fs]'
                   % (name, compress(yes), compress(no), time.time() - t0))
    # cross-check the decisive cases with a second solver/encoding
    for name, kinds, ks in (('square-coherent', ['sq'], (7, 21, 25, 26)), ('hexagon-coherent', ['hex'], (7, 21, 54, 55, 56, 57, 58, 59, 60))):
        r = {}
        for k in ks:
            ok, it, bits = solve(P, kinds, k, solver='m22', enc='kmtotalizer')
            r[k] = ok
        out.append('n=5 %s (Minisat22, kmtotalizer) decisive cases: %s' % (name, r))
    # acyclic witnesses by my own CEGAR
    wit = {}
    for name, kinds in (('square-coherent', ['sq']), ('hexagon-coherent', ['hex'])):
        for k in (7, 21):
            t0 = time.time()
            ok, it, bits = solve(P, kinds, k, acyclic=True)
            sinks, ac, coh = sinks_and_acyclic(P, bits)
            wit['%s_%d' % (kinds[0], k)] = bits
            out.append('n=5 %s acyclic witness k=%d: SAT=%s after %d CEGAR rounds [%.1fs]; independent check: sinks %d, acyclic %s, '
                       'square-coherent %s, hexagon-coherent %s; sinks %s'
                       % (name, k, ok, it, time.time() - t0, len(sinks), ac, coh['sq'], coh['hex'],
                          ' '.join(''.join(map(str, P.V[p])) for p in sinks)))
    # session witnesses: full tables printed for n4_L2_7 and n5_L1_7
    wpath = os.path.join(HERE, '..', 'permutohedron', 'witnesses.out')
    lines = open(wpath).read().splitlines()
    tab4 = json.loads(lines[-2].split(':', 1)[1].strip())
    tab5 = json.loads(lines[-1].split(':', 1)[1].strip())
    for (nn, tab, kind) in ((4, tab4, 'hex'), (5, tab5, 'sq')):
        Q = Perm(nn)
        bits = []
        for e, (u, v, i) in enumerate(Q.E):
            p = Q.V[u]
            a, b = p[i], p[i + 1]
            lo, hi = min(a, b), max(a, b)
            if kind == 'hex':
                key = '%d%d|%s' % (lo, hi, ''.join(str(z) for z in p if z not in (a, b)))
            else:
                key = '%d%d|{%s}' % (lo, hi, ''.join(str(z) for z in sorted(p[:i])))
            first = tab[key]
            # the edge points into the arrangement where `first` precedes the other element
            qv = Q.V[v]
            first_before_in_v = qv.index(first) < qv.index(b if first == a else a)
            bits.append(1 if first_before_in_v else 0)
        sinks, ac, coh = sinks_and_acyclic(Q, bits)
        out.append('session table n=%d (%s): sinks %d [%s], acyclic %s, square-coherent %s, hexagon-coherent %s'
                   % (nn, kind, len(sinks), ' '.join(''.join(map(str, Q.V[p])) for p in sinks), ac, coh['sq'], coh['hex']))
    # session sink SETS (printed without tables): realisable as exact acyclic coherent sink sets?
    for ln in lines:
        if ln.startswith('n5_L1_21') or ln.startswith('n5_L2_7') or ln.startswith('n5_L2_21'):
            nm = ln.split()[0]
            words = ln.split('sinks:')[1].split()
            kinds = ['sq'] if '_L1_' in nm else ['hex']
            fixed = set(P.idx[tuple(int(c) for c in wd)] for wd in words)
            ok, it, bits = solve(P, kinds, fixed_sinks=fixed, acyclic=True)
            if ok:
                sinks, ac, coh = sinks_and_acyclic(P, bits)
                out.append('session sink set %s (%d words): realisable exactly by an acyclic %s orientation: %s '
                           '(check: sinks %d equal to the set %s, acyclic %s, coherent %s)'
                           % (nm, len(words), kinds[0], ok, len(sinks), set(sinks) == fixed, ac, coh[kinds[0]]))
            else:
                out.append('session sink set %s: realisable: %s (CEGAR rounds %d)' % (nm, ok, it))
    return out


def compress(lst):
    if not lst:
        return '{}'
    lst = sorted(lst)
    parts = []
    s = prev = lst[0]
    for z in lst[1:]:
        if z == prev + 1:
            prev = z; continue
        parts.append('%d..%d' % (s, prev) if prev > s else '%d' % s)
        s = prev = z
    parts.append('%d..%d' % (s, prev) if prev > s else '%d' % s)
    return '{' + ', '.join(parts) + '}'


def partD():
    out = []
    P = Perm(6)
    t0 = time.time()
    ok, it, bits = solve(P, ['sq'], 90, mode='atleast', solver='g4', enc='seqcounter')
    if ok:
        sinks, ac, coh = sinks_and_acyclic(P, bits)
        out.append('n=6 square-coherent, >= 90 sinks (Glucose4, seqcounter): SAT [%.1fs]; independent check: sinks %d, '
                   'square-coherent %s, hexagon-coherent %s, acyclic %s' % (time.time() - t0, len(sinks), coh['sq'], coh['hex'], ac))
    else:
        out.append('n=6 square-coherent >= 90: %s [%.1fs]' % (ok, time.time() - t0))
    t0 = time.time()
    ok2, it2, bits2 = solve(P, ['sq'], 90, mode='atleast', acyclic=True, solver='g4', enc='seqcounter')
    if ok2:
        sinks, ac, coh = sinks_and_acyclic(P, bits2)
        out.append('n=6 square-coherent, >= 90 sinks and acyclic: SAT after %d CEGAR rounds [%.1fs]; check: sinks %d, acyclic %s, square-coherent %s'
                   % (it2, time.time() - t0, len(sinks), ac, coh['sq']))
    else:
        out.append('n=6 square-coherent, >= 90 sinks and acyclic: %s after %d rounds [%.1fs]' % (ok2, it2, time.time() - t0))
    return out


if __name__ == '__main__':
    parts = sys.argv[1:] or ['A', 'B', 'C', 'D']
    for part in parts:
        res = {'A': partA, 'B': partB, 'C': partC, 'D': partD}[part]()
        print('\n'.join(res))
        sys.stdout.flush()
