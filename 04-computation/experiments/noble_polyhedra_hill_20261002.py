"""Noble polyhedra (C. Hill, arXiv:2607.28711) seen from this repo: the fissary objects, the genus spectrum, the
number fields of the orbit locations, antipodal folding onto the Petersen graph, and an exact A5 dictionary for
5-tournaments.

Companion to 05-knowledge/results/noble_polyhedra_hill_fissary_20261002.md.

Data, fetched at run time into a temporary cache and not vendored here:
  Hill's model library (GPL-3.0): https://github.com/Plasmath/noble-tools-revised
    commit a801da7582fa927ba07e0af74d7c9c39445b7264
  Hill's appendix tables, parsed from the arXiv HTML: https://arxiv.org/html/2607.28711v1

Sections:
  A  surface invariants of all 148 models (146 noble + the 2 fissary models tI-F, rD-F): chi, orientability, genus
  B  independent check of Hill's fissary claim (which nobles have coplanar faces)
  C  abstract map automorphism groups and Petrie walks of selected models
  D  number fields: the fissary field Q(phi, sqrt(4 phi - 3)) of discriminant -5^2 19, and a census of the fields
     of all orbit locations in Hill's Tables 10-13
  E  antipodal folding of the dodecahedral nobles onto the Petersen graph / K10, with cover multiplicities
  F  5-tournaments through A5: cyclic triangles = dodecahedron vertices; regular 5-tournaments = faces of D-1, D-6;
     the gap 7 at n = 5 is Camion + Moon; the two classes of tournaments with three cyclic triangles
  F2 6-tournaments: the odd-cycle formula on all 32768, the gap 21 first appears at n = 6, the I-1 / I-4 split
  F3 7-tournaments: witnesses for H = 35 and H = 39
  G  where 7 and 21 occur
  H  rows of the paper's appendix tables that disagree with the model files, checked against the printed dual rows
  H2 the printed orbit locations of Tables 10-13 against the parameters recovered from the model files
"""
import hashlib
import html
import itertools
import json
import os
import random
import re
import tempfile
import urllib.parse
import urllib.request
from collections import Counter, defaultdict

import networkx as nx
import numpy as np
import sympy as sp

COMMIT = 'a801da7582fa927ba07e0af74d7c9c39445b7264'
REPO = 'Plasmath/noble-tools-revised'
PAPER = 'https://arxiv.org/html/2607.28711v1'
CACHE = os.path.join(tempfile.gettempdir(), 'noble-tools-revised-' + COMMIT[:7])


# ============================================================================ data
def fetch_library():
    os.makedirs(CACHE, exist_ok=True)
    marker = os.path.join(CACHE, 'index.json')
    if os.path.exists(marker):
        return json.load(open(marker))
    tree = json.load(urllib.request.urlopen(f'https://api.github.com/repos/{REPO}/git/trees/{COMMIT}?recursive=1'))
    files = {}
    for t in tree['tree']:
        p = t['path']
        if p.startswith('library/') and p.endswith('.off'):
            url = f'https://raw.githubusercontent.com/{REPO}/{COMMIT}/' + urllib.parse.quote(p)
            local = os.path.join(CACHE, os.path.basename(p))
            urllib.request.urlretrieve(url, local)
            files[os.path.basename(p)[:-4]] = local
    json.dump(files, open(marker, 'w'))
    return files


def fetch_paper_tables():
    """Hill's tables parsed from the arXiv HTML: {caption: [row cells]}.  Formulas are read from their alttext."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, 'paper.html')
    if not os.path.exists(path):
        req = urllib.request.Request(PAPER, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=120) as r, open(path, 'wb') as f:
            f.write(r.read())
    raw = open(path, 'rb').read()
    print(f'  paper HTML: {len(raw)} bytes, sha256 {hashlib.sha256(raw).hexdigest()[:16]}')
    s = raw.decode('utf-8')

    def clean(cell):
        def formula(m):
            alt = re.search(r'alttext="([^"]*)"', m.group(0))
            return html.unescape(alt.group(1)) if alt else '?'
        cell = re.sub(r'<math\b.*?</math>', formula, cell, flags=re.S)
        return re.sub(r'\s+', ' ', html.unescape(re.sub(r'<[^>]+>', ' ', cell))).strip()

    tables = {}
    for t in re.findall(r'<figure[^>]*class="ltx_table"[^>]*>(.*?)</figure>', s, flags=re.S):
        cap = re.search(r'<figcaption[^>]*>(.*?)</figcaption>', t, flags=re.S)
        rows = [[clean(c) for c in re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', r, flags=re.S)]
                for r in re.findall(r'<tr[^>]*>(.*?)</tr>', t, flags=re.S)]
        tables[clean(cap.group(1)) if cap else ''] = rows
    return tables


def table(tables, k):
    (rows,) = [v for c, v in tables.items() if c.startswith(f'Table {k}:')]
    return rows


def read_off(path):
    toks = [l.strip() for l in open(path, encoding='utf-8') if l.strip() and not l.startswith('#')]
    nv, nf, _ = map(int, toks[1].split())
    pts = [tuple(map(float, toks[2 + i].split()[:3])) for i in range(nv)]
    faces = []
    for i in range(nf):
        a = list(map(int, toks[2 + nv + i].split()))
        faces.append(tuple(a[1:1 + a[0]]))
    return pts, faces


def abstract_structure(faces):
    """Split every point whose face-link is disconnected (fissary models): one abstract vertex per link component.
    Corners are joined through a shared neighbour point, which is unambiguous because every segment carries exactly
    two face-sides (asserted)."""
    seg = Counter(frozenset((f[t], f[(t + 1) % len(f)])) for f in faces for t in range(len(f)))
    assert set(seg.values()) == {2}
    corners = defaultdict(list)
    for fi, f in enumerate(faces):
        for t in range(len(f)):
            corners[f[t]].append((fi, t))
    newid, count = {}, 0
    for v, cs in corners.items():
        G = nx.Graph(); G.add_nodes_from(range(len(cs)))
        by_w = defaultdict(list)
        for ci, (fi, t) in enumerate(cs):
            f = faces[fi]; k = len(f)
            for w in (f[(t - 1) % k], f[(t + 1) % k]):
                by_w[w].append(ci)
        for cis in by_w.values():
            G.add_edges_from(itertools.combinations(cis, 2))
        for comp in nx.connected_components(G):
            for ci in comp:
                newid[(v, cs[ci])] = count
            count += 1
    return [tuple(newid[(f[t], (fi, t))] for t in range(len(f))) for fi, f in enumerate(faces)]


def surface(nfaces):
    edges = defaultdict(list)
    for fi, f in enumerate(nfaces):
        for t in range(len(f)):
            edges[frozenset((f[t], f[(t + 1) % len(f)]))].append((fi, f[t]))
    assert all(len(l) == 2 for l in edges.values())
    V = len(set(x for f in nfaces for x in f)); E = len(edges); F = len(nfaces)
    adj = defaultdict(list)
    for (f1, a1), (f2, a2) in edges.values():
        rel = -1 if a1 == a2 else 1
        adj[f1].append((f2, rel)); adj[f2].append((f1, rel))
    sign, stack, ori = {0: 1}, [0], True
    while stack:
        f = stack.pop()
        for g, rel in adj[f]:
            if g not in sign:
                sign[g] = sign[f] * rel; stack.append(g)
            elif sign[g] != sign[f] * rel:
                ori = False
    assert len(sign) == F
    chi = V - E + F
    return V, E, F, chi, ori, ((2 - chi) // 2 if ori else 2 - chi)


# ============================================================================ A, B
def section_A(lib):
    print('== A. surface invariants of all models')
    rows = {}
    for name, path in sorted(lib.items()):
        pts, faces = read_off(path)
        V, E, F, chi, ori, g = surface(abstract_structure(faces))
        rows[name] = dict(points=len(pts), V=V, E=E, F=F, chi=chi, orientable=ori, genus=g,
                          p=sorted(set(len(f) for f in faces)))
    n_noble = sum(1 for k in rows if not k.endswith('-F'))
    print(f'  models: {len(rows)} ({n_noble} noble + fissary {sorted(k for k in rows if k.endswith("-F"))})')
    assert n_noble == 146
    ori = Counter(r['genus'] for r in rows.values() if r['orientable'])
    non = Counter(r['genus'] for r in rows.values() if not r['orientable'])
    print('  orientable genus spectrum (genus: count):', dict(sorted(ori.items())))
    print('  non-orientable (crosscap number: count):', dict(sorted(non.items())))
    for g in (7, 21):
        print(f'  genus {g}:', sorted(k for k, r in rows.items() if r['orientable'] and r['genus'] == g))
    print('  crosscap 14:', sorted(k for k, r in rows.items() if not r['orientable'] and r['genus'] == 14))
    for k in ('tI-F', 'rD-F', 'gD-19.1', 'gD-28.1', 'D-4', 'D-5'):
        r = rows[k]
        print(f'  {k}: {r["points"]} points, abstract V={r["V"]}, E={r["E"]}, F={r["F"]}, chi={r["chi"]}, '
              f'{"orientable genus" if r["orientable"] else "crosscap"} {r["genus"]}')
    assert rows['D-4']['orientable'] and rows['D-4']['genus'] == 21
    assert [k for k, r in rows.items() if r['orientable'] and r['genus'] == 21] == ['D-4']
    return rows


def section_B(lib):
    print('== B. Hill\'s fissary claim: nobles with two faces in one plane')
    found = {}
    for name, path in sorted(lib.items()):
        pts, faces = read_off(path)
        X = np.array(pts)
        poles = []
        for f in faces:
            P = X[list(f)]; c = P.mean(axis=0)
            u, s, vt = np.linalg.svd(P - c)
            assert s[-1] < 1e-6
            nrm = vt[-1]; d = float(nrm @ c)
            if d < 0: nrm, d = -nrm, -d
            poles.append(nrm / d)
        poles = np.array(poles)
        used = np.zeros(len(poles), bool); mult = Counter()
        for i in range(len(poles)):
            if used[i]: continue
            same = np.where(np.linalg.norm(poles - poles[i], axis=1) < 1e-7)[0]
            used[same] = True; mult[len(same)] += 1
        if any(m > 1 for m in mult):
            found[name] = dict(mult)
    for k, v in found.items():
        print(f'  {k}: plane multiplicities {v}')
    assert set(found) == {'D-4', 'D-5', 'gD-19.1', 'gD-28.1'}
    print('  exactly the four parents named by Hill (their duals tI-F, rD-F, D-F1, D-F2 are fissary)')


# ============================================================================ C
def flag_system(nfaces):
    flags = [(fi, t, s) for fi, f in enumerate(nfaces) for t in range(len(f)) for s in (0, 1)]
    idx = {fl: i for i, fl in enumerate(flags)}
    occ = defaultdict(list)
    for fi, f in enumerate(nfaces):
        for t in range(len(f)):
            occ[frozenset((f[t], f[(t + 1) % len(f)]))].append((fi, t))
    r0, r1, r2 = [0] * len(flags), [0] * len(flags), [0] * len(flags)
    for i, (fi, t, s) in enumerate(flags):
        f = nfaces[fi]; k = len(f); v = f[t]
        if s == 0:
            r0[i] = idx[(fi, (t + 1) % k, 1)]; r1[i] = idx[(fi, t, 1)]; e = frozenset((v, f[(t + 1) % k])); pos = t
        else:
            r0[i] = idx[(fi, (t - 1) % k, 0)]; r1[i] = idx[(fi, t, 0)]; e = frozenset((f[(t - 1) % k], v)); pos = (t - 1) % k
        (gj, u), = [o for o in occ[e] if o != (fi, pos)]
        g = nfaces[gj]
        r2[i] = idx[(gj, u, 0)] if g[u] == v else idx[(gj, (u + 1) % len(g), 1)]
    return len(flags), (r0, r1, r2)


def aut_order(nfaces):
    n, rs = flag_system(nfaces)
    cnt = 0
    for target in range(n):
        img, stack, ok = {0: target}, [0], True
        while stack and ok:
            x = stack.pop()
            for r in rs:
                y, iy = r[x], r[img[x]]
                if y in img:
                    ok = img[y] == iy
                    if not ok: break
                else:
                    img[y] = iy; stack.append(y)
        cnt += ok and len(img) == n and len(set(img.values())) == n
    return cnt, n


def section_C(lib):
    print('== C. abstract map automorphism groups (flag-orbit counts) and Petrie walks')
    for name in ['D-1', 'D-3', 'D-4', 'D-5', 'D-6', 'I-2', 'I-3', 'tC-1.1', 'sC-5.1', 'gD-19.1', 'gD-28.1', 'tI-F',
                 'rD-F']:
        nf = abstract_structure(read_off(lib[name])[1])
        a, n = aut_order(nf)
        nfl, (r0, r1, r2) = flag_system(nf)
        seen, walks = [False] * nfl, Counter()
        for i in range(nfl):
            if seen[i]: continue
            j, L = i, 0
            while True:
                seen[j] = True; j = r2[r1[r0[j]]]; L += 1
                if j == i: break
            walks[L] += 1
        print(f'  {name:8s} flags {n:5d}  |Aut| {a:4d}  flag orbits {n // a:2d}  Petrie (r0r1r2) orbit lengths {dict(walks)}')


# ============================================================================ D
S5 = sp.sqrt(5)
PHI = (1 + S5) / 2
RADICANDS = {'phi': PHI, '2': sp.Integer(2), '1+4phi': 1 + 4 * PHI, '4phi-3': 4 * PHI - 3}


def tex_poly(s):
    a, b = sp.symbols('a b')
    s = s.replace('^{', '**(').replace('}', ')')
    s = re.sub(r'(\d)\s*([ab])', r'\1*\2', s)
    return sp.Poly(sp.sympify(s, locals={'a': a, 'b': b}))


def location_rows(tables):
    """Tables 10-13: orbit -> [(parameter name, printed location, minimal polynomial)], and orbit -> table number"""
    locs, src = {}, {}
    for k in (10, 11, 12, 13):
        cur = None
        for r in table(tables, k)[1:]:
            if len(r) == 3:
                cur, r = r[0], r[1:]
                locs[cur] = []; src[cur] = k
            pol = tex_poly(r[1])
            locs[cur].append((str(pol.gens[0]), float(r[0]), pol))
    return locs, src


def is_square_Qphi(e):
    x = sp.Symbol('x')
    return any(sp.degree(f, x) == 1 for f, _ in sp.factor_list(x ** 2 - sp.expand(e), extension=S5)[1])


def field_class(pol):
    """The field generated by one root, up to isomorphism.  A quartic that splits over Q(sqrt5) generates
    Q(phi)(sqrt D), D the discriminant of a quadratic factor; that field is Q(phi, sqrt R) iff D R or conj(D) R is a
    square in Q(phi)."""
    d, g = pol.degree(), pol.gens[0]
    if d == 1:
        return 'Q'
    if d == 2:
        disc = sp.discriminant(pol)
        return 'Q(phi)' if sp.sqrt(sp.Rational(disc, 5)).is_rational else f'quadratic, disc {disc}'
    if d == 4:
        fl = sp.factor_list(pol.as_expr(), extension=S5)[1]
        if len(fl) == 2:
            c = sp.Poly(fl[0][0], g).all_coeffs()
            D = sp.expand((c[1] / c[0]) ** 2 - 4 * c[2] / c[0])
            for name, R in RADICANDS.items():
                if any(is_square_Qphi(DD * R) for DD in (D, D.subs(S5, -S5))):
                    return f'Q(phi, sqrt({name}))'
            return f'quartic over Q(phi), radicand {D}'
        return 'quartic, not over Q(phi)'
    return f'degree {d}'


def section_D(tables):
    print('== D. number fields of the orbit locations')
    a = sp.symbols('a')
    polys = {'tI-F': a**4 + 2*a**3 + 2*a**2 + a - 1, 'rD-F': a**4 - 2*a**3 + 2*a**2 - a - 1,
             'gD-19 (a)': a**4 + a**3 - 2*a**2 + 2*a - 1, 'gD-28 (a = b)': a**4 - a**3 - 2*a**2 - 2*a - 1}
    for k, p in polys.items():
        d = sp.discriminant(p, a)
        fac = sp.factor(p, extension=S5)
        print(f'  {k}: {p}   disc {d} = {sp.factorint(d)};  over Q(sqrt5): {fac}')
        assert d == -475
    x = (-1 + sp.sqrt(4 * PHI - 3)) / 2
    y = x + 1
    g28 = (PHI + sp.sqrt(5 * PHI + 1)) / 2
    g19 = (-PHI + sp.sqrt(5 * PHI + 1)) / 2
    checks = {
        'tI-F location x = (-1 + sqrt(4phi-3))/2 has Hill\'s minimal polynomial': sp.minimal_polynomial(x, a) - polys['tI-F'],
        'rD-F location y = x + 1 has Hill\'s minimal polynomial': sp.minimal_polynomial(y, a) - polys['rD-F'],
        'gD-28 a = (phi + sqrt(5phi+1))/2 has Hill\'s minimal polynomial': sp.minimal_polynomial(g28, a) - polys['gD-28 (a = b)'],
        'gD-19 a = (-phi + sqrt(5phi+1))/2 has Hill\'s minimal polynomial': sp.minimal_polynomial(g19, a) - polys['gD-19 (a)'],
        'x(x+1) = 1/phi': sp.simplify(x * (x + 1) - 1 / PHI),
        'y(y-1) = 1/phi': sp.simplify(y * (y - 1) - 1 / PHI),
        '(5phi+1)(4phi-3) = (phi+4)^2, so sqrt(5phi+1) sqrt(4phi-3) = phi + 4':
            sp.simplify(sp.expand((5 * PHI + 1) * (4 * PHI - 3)) - sp.expand((PHI + 4) ** 2)),
    }
    for k, v in checks.items():
        print(f'  {k}: {sp.expand(v) == 0}')
        assert sp.expand(v) == 0
    print('  numerical: x =', sp.N(x, 15), ' y =', sp.N(y, 15), ' gD-28 a =', sp.N(g28, 15), ' gD-19 a =', sp.N(g19, 15))
    N = lambda u, v: u * u + u * v - v * v   # norm of u + v phi
    print('  norms of the radicands: N(phi) =', N(0, 1), ', N(2) =', N(2, 0), ', N(1+4phi) =', N(1, 4),
          ', N(4phi-3) =', N(-3, 4), ', N(5phi+1) =', N(1, 5))
    norms = set(abs(N(u, v)) for u in range(-60, 61) for v in range(-60, 61))
    print('  |N(z)| for z in Z[phi] (|coefficients| <= 60) avoids 3, 7, 21:', not ({3, 7, 21} & norms),
          '(3 and 7 are inert in Q(sqrt5))')
    assert not ({3, 7, 21} & norms) and sp.legendre_symbol(5, 3) == sp.legendre_symbol(5, 7) == -1

    # census of Tables 10-13 (minimal polynomials of the orbit locations)
    locs, src = location_rows(tables)
    orbits = {o: [pol for _, _, pol in lst] for o, lst in locs.items()}
    off_root = [(o, v) for o, lst in locs.items() for v, loc, pol in lst
                if np.min(np.abs(np.roots([float(c) for c in pol.all_coeffs()]) - loc)) > 1e-6]
    print(f'  Tables 10-13: {len(orbits)} orbits, {sum(map(len, orbits.values()))} location polynomials; printed '
          f'locations that are not roots of the printed polynomial: {off_root} (copy errors, see section H)')
    assert off_root == [('sD-3', 'b'), ('sD-16', 'b'), ('gD-15', 'a')]
    fc = {(o, i): field_class(p) for o, lst in orbits.items() for i, p in enumerate(lst)}
    by_class = defaultdict(list)
    for (o, i), c in fc.items():
        by_class[c].append(f'{o}:{orbits[o][i].gens[0]}')
    for c in sorted(by_class):
        if c.startswith('Q(phi, sqrt'):
            print(f'  {c}: {by_class[c]}')
    t10 = {c: sorted(o for (o, i), cc in fc.items() if cc == c and src[o] == 10) for c in by_class}
    assert t10['Q(phi, sqrt(phi))'] == ['rD-4', 'tD-3', 'tI-2']
    assert t10['Q(phi, sqrt(2))'] == ['rD-1', 'tI-1']
    assert t10['Q(phi, sqrt(1+4phi))'] == ['rD-6', 'rD-7', 'tD-1', 'tI-4', 'tI-7']
    assert t10['Q(phi, sqrt(4phi-3))'] == ['rD-F', 'tI-F']
    print('  Table 10 quartics by field (polynomial discriminants):')
    for c in ('Q(phi, sqrt(phi))', 'Q(phi, sqrt(2))', 'Q(phi, sqrt(1+4phi))', 'Q(phi, sqrt(4phi-3))'):
        print(f'    {c}:', {o: sp.factorint(sp.discriminant(orbits[o][0])) for o in t10[c]})
    others = sorted(c for c in t10 if t10[c] and not c.startswith('Q(phi, sqrt'))
    print('  Table 10, other locations:', {c: t10[c] for c in others})
    print('    tO-1 cubic discriminant:', sp.factorint(sp.discriminant(orbits['tO-1'][0])))
    inK = sorted(o for o, lst in orbits.items()
                 if all(fc[(o, i)] in ('Q', 'Q(phi)', 'Q(phi, sqrt(4phi-3))') for i in range(len(lst)))
                 and any(fc[(o, i)] == 'Q(phi, sqrt(4phi-3))' for i in range(len(lst))))
    print('  orbits all of whose location parameters lie in the fissary field K = Q(phi, sqrt(4phi-3)):', inK)
    assert inK == ['gD-19', 'gD-28', 'rD-F', 'sD-9', 'tI-F']
    sd9 = orbits['sD-9']
    assert sp.expand(sd9[0].as_expr() - polys['tI-F'].subs(a, sd9[0].gens[0])) == 0
    print('  so K holds the four fissary-related orbits and one more: sD-9, whose a has the tI-F minimal polynomial and')
    print('  whose b is the root 1/phi of', sd9[1].as_expr(), '; its noble sD-9.1 has no coplanar faces (section B)')


# ============================================================================ E
def section_E(lib):
    print('== E. dodecahedral nobles folded by the antipodal map')
    ref = np.array(read_off(lib['D-1'])[0])
    Gd = nx.Graph()
    for f in read_off(lib['D-1'])[1]:
        Gd.add_edges_from((f[t], f[(t + 1) % len(f)]) for t in range(len(f)))
    dist = dict(nx.all_pairs_shortest_path_length(Gd))
    anti = {i: int(np.argmin(np.linalg.norm(ref + ref[i], axis=1))) for i in range(20)}
    assert all(dist[i][anti[i]] == 5 for i in range(20))
    cls = {i: min(i, anti[i]) for i in range(20)}
    pet = nx.petersen_graph()
    print('  quotient edges: a pair of antipodal classes is a Petersen edge if its lifts are at distances 1 and 4,')
    print('  and a J(5,2) edge if they are at distances 2 and 3 (each quotient edge has 4 lifts, 2 of each distance)')
    for name in ['D-1', 'D-2', 'D-3', 'D-4', 'D-5', 'D-6', 'D-7']:
        pts, faces = read_off(lib[name])
        X = np.array(pts)
        m = [int(np.argmin(np.linalg.norm(ref - x * (np.linalg.norm(ref[0]) / np.linalg.norm(x)), axis=1))) for x in X]
        F = [tuple(m[v] for v in f) for f in faces]
        E = set(frozenset((f[t], f[(t + 1) % len(f)])) for f in F for t in range(len(f)))
        dcls = Counter(dist[min(e)][max(e)] for e in E)
        Fset = set(frozenset(f) for f in F)
        sym = all(frozenset(anti[v] for v in f) in Fset for f in F)
        Q = nx.Graph(); Q.add_edges_from((cls[min(e)], cls[max(e)]) for e in E)
        assert all(cls[min(e)] != cls[max(e)] for e in E)
        kind = 'Petersen' if nx.is_isomorphic(Q, pet) else ('K10' if Q.number_of_edges() == 45 else '?')
        mult_P, mult_J = Counter(), Counter()
        for e in E:
            u, v = cls[min(e)], cls[max(e)]
            (mult_P if dist[u][v] in (1, 4) else mult_J)[frozenset((u, v))] += 1
        print(f'  {name}: edges at dodecahedral distances {dict(sorted(dcls.items()))}; centrally symmetric {sym}; '
              f'antipodal quotient = {kind}; Petersen edges covered {len(mult_P)}/15 x{sorted(set(mult_P.values()))}; '
              f'J(5,2) edges covered {len(mult_J)}/30 x{sorted(set(mult_J.values()))}')
        if name == 'D-5':
            assert len(mult_P) == 15 and set(mult_P.values()) == {4} and len(mult_J) == 30 and set(mult_J.values()) == {1}
            d3 = [e for e in E if dist[min(e)][max(e)] == 3]
            assert all(frozenset(anti[v] for v in e) not in E for e in d3)
            print('    D-5 takes exactly one of the two distance-3 lifts of every J(5,2) edge; the central inversion,')
            print('    which maps the chiral D-5 to its mirror image, takes the other one')


# ============================================================================ F
def count_hamiltonian_paths(n, beats):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1 << n):
        for v in range(n):
            if dp[S][v]:
                for w in range(n):
                    if not S >> w & 1 and beats[(v, w)]:
                        dp[S | 1 << w][w] += dp[S][v]
    return sum(dp[(1 << n) - 1])


def galois_edge_rank(coords):
    """coords(f) lists the vertices of a polyhedron as functions of f = phi.  Return the Euclidean distance class (1 =
    nearest) of the image, under f -> -1/phi, of the edges (nearest pairs) of the original."""
    phi, phibar = (1 + 5 ** 0.5) / 2, (1 - 5 ** 0.5) / 2
    X, Y = np.array(coords(phi), float), np.array(coords(phibar), float)
    dX = np.linalg.norm(X[:, None] - X[None], axis=2); dY = np.linalg.norm(Y[:, None] - Y[None], axis=2)
    levels = sorted(set(np.round(dY[dY > 1e-9], 9)))
    edge = dX[dX > 1e-9].min()
    ranks = {levels.index(round(dY[i, j], 9)) + 1 for i, j in zip(*np.where(np.abs(dX - edge) < 1e-9))}
    (r,) = ranks
    return r


def tournament(n, mask):
    beats = {}
    for t, (a, b) in enumerate(itertools.combinations(range(n), 2)):
        beats[(a, b)] = bool(mask >> t & 1); beats[(b, a)] = not beats[(a, b)]
    return beats


def section_F(lib):
    print('== F. 5-tournaments through A5')
    n = 5
    def comp(p, q): return tuple(p[q[i]] for i in range(n))
    def inv(p):
        r = [0] * n
        for i, x in enumerate(p): r[x] = i
        return tuple(r)
    def cyc(c):
        p = list(range(n))
        for i in range(len(c)): p[c[i]] = c[(i + 1) % len(c)]
        return tuple(p)
    def even(p):
        s, seen = 0, set()
        for i in range(n):
            j, L = i, 0
            while j not in seen:
                seen.add(j); j = p[j]; L += 1
            s += (L - 1) if L else 0
        return s % 2 == 0
    def order(p):
        q, k = p, 1
        while q != tuple(range(n)): q = comp(p, q); k += 1
        return k
    A5 = [p for p in itertools.permutations(range(n)) if even(p)]
    three = [p for p in A5 if order(p) == 3]; five = [p for p in A5 if order(p) == 5]
    idx = {p: i for i, p in enumerate(three)}
    conj = lambda g, p: comp(comp(g, p), inv(g))
    # A5-orbitals on pairs of 3-cycles.  Two of them are dodecahedron graphs (geometrically: distance 1 and distance
    # 4, the edges of D-1 and of D-6), swapped by the outer automorphism.  We take the first, so 'D-1' versus 'D-6'
    # below is fixed only up to the outer automorphism; every count below is symmetric under it.
    orbit_of = {}
    for a, b in itertools.combinations(range(20), 2):
        key = min(tuple(sorted((idx[conj(g, three[a])], idx[conj(g, three[b])]))) for g in A5)
        orbit_of.setdefault(key, []).append((a, b))
    dodecs = [E for E in orbit_of.values() if len(E) == 30 and nx.is_isomorphic(nx.Graph(E), nx.dodecahedral_graph())]
    assert len(dodecs) == 2
    D = nx.Graph(dodecs[0])
    dist = dict(nx.all_pairs_shortest_path_length(D))
    assert all(dist[i][idx[inv(three[i])]] == 5 for i in range(20))
    print('  the 20 three-cycles of A5 carry two A5-invariant dodecahedron graphs (swapped by the outer automorphism);')
    print('  in either, the inverse 3-cycle is the antipode')
    layers = set()
    rho_orbit = lambda rho, i: frozenset(idx[p] for p in itertools.accumulate(range(4), lambda p, _: conj(rho, p), initial=three[i]))
    for rho in five:
        for i in range(20):
            layers.add(rho_orbit(rho, i))
    shape = lambda S: tuple(sorted(dist[a][b] for a, b in itertools.combinations(S, 2)))
    print('  layers (orbits of 5-fold rotations):', len(layers), dict(Counter(shape(L) for L in layers)))
    # check the two layer types against the geometric D-1 / D-6 faces: map the abstract dodecahedron to D-1
    ref = np.array(read_off(lib['D-1'])[0])
    Gd = nx.Graph()
    for f in read_off(lib['D-1'])[1]:
        Gd.add_edges_from((f[t], f[(t + 1) % len(f)]) for t in range(len(f)))
    iso = next(nx.algorithms.isomorphism.GraphMatcher(D, Gd).isomorphisms_iter())
    X6 = np.array(read_off(lib['D-6'])[0])
    m6 = [int(np.argmin(np.linalg.norm(ref - x * (np.linalg.norm(ref[0]) / np.linalg.norm(x)), axis=1))) for x in X6]
    faces1 = set(frozenset(f) for f in read_off(lib['D-1'])[1])
    faces6 = set(frozenset(m6[v] for v in f) for f in read_off(lib['D-6'])[1])
    img = lambda L: frozenset(iso[x] for x in L)
    assert set(img(L) for L in layers) == faces1 | faces6
    print('  layers = the 12 faces of the dodecahedron D-1 + the 12 pentagram faces of the great stellated dodecahedron D-6')
    pairs = list(itertools.combinations(range(n), 2))
    def relabel(mask, g):
        out = 0
        for t, (a, b) in enumerate(pairs):
            win = (a, b) if mask >> t & 1 else (b, a)
            ga, gb = g[win[0]], g[win[1]]
            if ga < gb: out |= 1 << pairs.index((ga, gb))
        return out
    def is_strong(beats):
        for fwd in (True, False):
            seen, stack = {0}, [0]
            while stack:
                u = stack.pop()
                for w in range(n):
                    if w not in seen and (beats[(u, w)] if fwd else beats[(w, u)]):
                        seen.add(w); stack.append(w)
            if len(seen) < n:
                return False
        return True
    table_, regular, camion, strong_scores, three_tri = Counter(), [], Counter(), Counter(), []
    for mask in range(1 << 10):
        beats = tournament(n, mask)
        c3 = []
        for a, b, c in itertools.combinations(range(n), 3):
            if beats[(a, b)] and beats[(b, c)] and beats[(c, a)]: c3.append(idx[cyc((a, b, c))])
            elif beats[(a, c)] and beats[(c, b)] and beats[(b, a)]: c3.append(idx[cyc((a, c, b))])
        rhos = [cyc(q) for q in ((0,) + r for r in itertools.permutations(range(1, n)))
                if all(beats[(q[i], q[(i + 1) % n])] for i in range(n))]
        c5 = len(rhos)
        H = count_hamiltonian_paths(n, beats)
        assert H == 1 + 2 * (len(c3) + c5)
        table_[(len(c3), c5, H)] += 1
        st = is_strong(beats)
        camion[(st, len(c3) >= 3, c5 >= 1, H)] += 1
        if st:
            strong_scores[tuple(sorted(sum(beats[(v, w)] for w in range(n) if w != v) for v in range(n)))] += 1
        if len(c3) == 5:
            regular.append((mask, frozenset(c3)))
        if len(c3) == 3:
            three_tri.append((mask, beats, frozenset(c3), rhos))
    print('  (c3, c5, H) over the 1024 labelled 5-tournaments:', dict(sorted(table_.items())))
    print('  H takes the values', sorted(set(h for _, _, h in table_)), '- the first gap of the H-spectrum, 7, appears here:')
    print('    c3 <= 2 forces c5 = 0 and c3 = 3 forces c5 = 1, so c3 + c5 = 3 never occurs')
    # the gap 7 is Camion + Moon
    assert all(st == big3 == ham for st, big3, ham, _ in camion)
    Hns = sorted({h for st, _, _, h in camion if not st}); Hst = sorted({h for st, _, _, h in camion if st})
    print('  strong <=> c3 >= 3 <=> Hamiltonian cycle (c5 >= 1), on all 1024; H of non-strong:', Hns, '; of strong:', Hst)
    print('  strong score sequences:', dict(strong_scores))
    assert Hns == [1, 3, 5] and min(Hst) == 9
    print('  so H = 7 is skipped because a strong 5-tournament has c3 >= 3 (Moon: >= n - 2) and c5 >= 1 (Camion),')
    print('  hence H >= 9, while a non-strong one has c5 = 0 and c3 <= 2, hence H <= 5')
    # regular tournaments
    assert all(S in layers for _, S in regular) and len(regular) == 24
    print('  the cyclic triangles of each of the 24 regular 5-tournaments form a layer: a face of D-1 or of D-6')
    masks = {m for m, _ in regular}
    orbits, seen = [], set()
    for m in sorted(masks):
        if m in seen: continue
        orb = {relabel(m, g) for g in A5}
        seen |= orb; orbits.append(orb)
    S_of = dict(regular)
    types = [Counter(img(S_of[m]) in faces1 for m in orb) for orb in orbits]
    print('  A5-orbits of regular tournaments:', [len(o) for o in orbits], '; D-1-face membership per orbit:', types)
    assert len(orbits) == 2 and all(len(t) == 1 for t in types)
    print('  -> the two A5-orbits (swapped by odd permutations, i.e. by the Galois conjugation sqrt5 -> -sqrt5) are')
    print('     exactly the faces of the dodecahedron and the faces of the great stellated dodecahedron')
    rank = galois_edge_rank(lambda f: [(s1, s2, s3) for s1 in (1, -1) for s2 in (1, -1) for s3 in (1, -1)]
                            + [c for s1 in (1, -1) for s2 in (1, -1)
                               for c in ((0, s1 / f, s2 * f), (s1 / f, s2 * f, 0), (s2 * f, 0, s1 / f))])
    print(f'  Galois sqrt5 -> -sqrt5 on exact dodecahedron coordinates sends edges to distance class {rank} of the image')
    assert rank == 4   # the image edges are the edges of a great stellated dodecahedron
    # tournaments with exactly three cyclic triangles
    S5perms = list(itertools.permutations(range(n)))
    kinds = Counter()
    for mask, beats, S, rhos in three_tri:
        canon_form = min(relabel(mask, g) for g in S5perms)
        inl = [L for L in layers if S <= L]
        if inl:
            (L,), (rho,) = inl, rhos
            assert any(rho_orbit(rho, i) == L for i in S)
            kind = 'in a layer of its 5-cycle rho, a ' + ('D-1' if img(L) in faces1 else 'D-6') + ' face'
        else:
            common = [w for w in range(20) if all(D.has_edge(w, s) for s in S)]
            (w,) = common
            tri = [i for i in range(n) if three[w][i] != i]
            a, b, c = tri
            cyclic = (beats[(a, b)] and beats[(b, c)] and beats[(c, a)]) or (beats[(a, c)] and beats[(c, b)] and beats[(b, a)])
            assert shape(S) == (2, 2, 2) and not cyclic
            kind = 'in no layer: the 3 neighbours of a vertex whose triple is transitive in T'
        kinds[(kind, canon_form)] += 1
    print('  the 240 tournaments with c3 = 3 (score sequence (1,1,2,3,3)), by the position of their 3 cyclic triangles:')
    for (kind, cf), v in sorted(kinds.items()):
        print(f'    {v:4d}  {kind}   [isomorphism class {cf}]')
    classes = defaultdict(set)
    for (kind, cf) in kinds:
        classes[cf].add(kind.split(',')[0])
    assert len(classes) == 2 and sorted(len(v) for v in classes.values()) == [1, 1]
    print('  -> the two isomorphism classes with c3 = 3 are told apart exactly by layer / no layer')


# ============================================================================ F2
def section_F2():
    print('== F2. n = 6: A5 = PSL(2,5) on the 6 axes of the icosahedron')
    n = 6
    Hs = Counter()
    for mask in range(1 << 15):
        beats = tournament(n, mask)
        H = count_hamiltonian_paths(n, beats)
        Hs[H] += 1
        c3 = [t for t in itertools.combinations(range(n), 3)
              if (beats[(t[0], t[1])] and beats[(t[1], t[2])] and beats[(t[2], t[0])])
              or (beats[(t[0], t[2])] and beats[(t[2], t[1])] and beats[(t[1], t[0])])]
        c5 = sum(1 for S5 in itertools.combinations(range(n), 5) for r in itertools.permutations(S5[1:])
                 if all(beats[(c, d)] for c, d in zip((S5[0],) + r, r + (S5[0],))))
        d33 = sum(1 for a, b in itertools.combinations(c3, 2) if not set(a) & set(b))
        assert H == 1 + 2 * (len(c3) + c5) + 4 * d33
    missing = sorted(set(range(1, max(Hs) + 1, 2)) - set(Hs))
    print('  odd-cycle formula H = 1 + 2(c3 + c5) + 4 d33 verified on all 32768 labelled 6-tournaments')
    print('  H values at n = 6:', sorted(Hs), '; missing odd values up to the maximum 45:', missing)
    assert missing == [7, 21, 35, 39]
    print('  so 21 first becomes a gap at n = 6 (7 at n = 5); 35 and 39 are attained at n = 7 (section F3)')
    phi = (1 + 5 ** 0.5) / 2
    V = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            V += [(0, s1, s2 * phi), (s1, s2 * phi, 0), (s2 * phi, 0, s1)]
    V = np.array(V, float)
    axes = []
    for v in V:
        if not any(np.allclose(v, w) or np.allclose(v, -w) for w in axes): axes.append(v)
    ax = lambda v: next(i for i, w in enumerate(axes) if np.allclose(v, w) or np.allclose(v, -w))
    dist = np.linalg.norm(V[:, None] - V[None], axis=2)
    dv = sorted(set(np.round(dist[dist > 1e-9], 6)))
    tri = lambda e: {tuple(sorted(ax(V[i]) for i in t)) for t in itertools.combinations(range(12), 3)
                     if all(abs(dist[a, b] - e) < 1e-6 for a, b in itertools.combinations(t, 2))}
    t1, t4 = tri(dv[0]), tri(dv[1])   # icosahedron I-1 faces, great icosahedron I-4 faces
    allt = set(itertools.combinations(range(6), 3))
    assert len(t1) == len(t4) == 10 and t1 | t4 == allt
    assert all(tuple(sorted(set(range(6)) - set(t))) in t4 for t in t1)
    rank = galois_edge_rank(lambda f: [c for s1 in (1, -1) for s2 in (1, -1)
                                       for c in ((0, s1, s2 * f), (s1, s2 * f, 0), (s2 * f, 0, s1))])
    print(f'  Galois sqrt5 -> -sqrt5 on exact icosahedron coordinates sends edges to distance class {rank} of the image')
    assert rank == 2   # the image edges are the edges of a great icosahedron
    print('  the 20 triples of the 6 axes = 10 face triples of the icosahedron I-1 + 10 of the great icosahedron I-4')
    print('  (the two hemi-icosahedra), and complementation swaps the classes.  So every complementary pair of triples,')
    print('  in particular every vertex-disjoint pair of cyclic triangles counted by d33, is {I-1 triple, I-4 triple},')
    print('  for every tournament and every labelling of its vertices by the axes: a tautology, no tournament content')


def section_F3():
    print('== F3. n = 7: H = 35 and H = 39 (missing at n = 6) are attained')
    n, rng, found = 7, random.Random(20261002), {}
    for _ in range(100000):
        mask = rng.getrandbits(21)
        H = count_hamiltonian_paths(n, tournament(n, mask))
        if H in (35, 39) and H not in found:
            found[H] = mask
            if len(found) == 2:
                break
    for H, mask in sorted(found.items()):
        beats = tournament(n, mask)
        score = sorted(sum(beats[(v, w)] for w in range(n) if w != v) for v in range(n))
        print(f'  H = {H}: arc mask {mask} (bit t orients the t-th pair (a, b), a < b, as a -> b), scores {score}')
    assert set(found) == {35, 39}


# ============================================================================ G, H
def rotation_group(X):
    """The rotations mapping the point set X (an orbit, so on a sphere about 0) onto itself."""
    U = X / np.linalg.norm(X, axis=1)[:, None]; G = U @ U.T; n = len(U)
    b = next(j for j in range(n) if abs(abs(G[0, j]) - 1) > 1e-6)
    c = next(k for k in range(n) if abs(np.linalg.det(U[[0, b, k]])) > 1e-6)
    Binv = np.linalg.inv(U[[0, b, c]].T)
    out = []
    for x in range(n):
        for y in np.where(np.abs(G[x] - G[0, b]) < 1e-7)[0]:
            for z in np.where((np.abs(G[x] - G[0, c]) < 1e-7) & (np.abs(G[y] - G[b, c]) < 1e-7))[0]:
                M = U[[x, y, z]].T @ Binv
                if np.linalg.det(M) > 0 and np.allclose(M @ M.T, np.eye(3), atol=1e-7):
                    img = U @ M.T
                    if np.max(np.min(np.linalg.norm(img[:, None] - U[None], axis=2), axis=1)) < 1e-6:
                        out.append(M)
    return out


def orbit_parameters(orbit_type, X):
    """Hill's parameters of the orbit X = V_G(d1, d2, d3): d_i is the distance of an orbit point to the mirror r_i of
    its Moebius triangle, r_i being opposite the vertex u_i, where u = (3-fold, 5-fold, 2-fold axis) for *532 and
    (3-fold, 2-fold, 4-fold axis) for *432 (Hill: V(1,0,0) = D or C, V(0,1,0) = I or CO, V(0,0,1) = ID or O).
    Normalised as in Hill's Table 1: tI, tO = V(0,a,1); tD, rC = V(a,0,1); rD, tC = V(a,1,0); s*, g* = V(a,b,1)."""
    Rs = rotation_group(X)
    big = {60: 5, 24: 4}[len(Rs)]
    ax = defaultdict(list)
    for M in Rs:
        if np.allclose(M, np.eye(3)): continue
        w, V = np.linalg.eig(M)
        v = np.real(V[:, np.argmin(np.abs(w - 1))]); v /= np.linalg.norm(v)
        k = next(k for k in (2, 3, 4, 5) if np.allclose(np.linalg.matrix_power(M, k), np.eye(3), atol=1e-6))
        for s in (v, -v):
            if not any(np.allclose(s, u, atol=1e-6) for u in ax[k]):
                ax[k].append(s)
    ax[2] = [s for s in ax[2] if not any(abs(s @ q) > 1 - 1e-6 for q in ax[4])]   # 4-fold axes carry half-turns
    P, T, S = ax[big], ax[3], ax[2]
    near = lambda A, B: max(p @ q for p in A for q in B)
    cpt, cps, cts = near(P, T), near(P, S), near(T, S)
    plane = lambda p, q: np.cross(p, q) / np.linalg.norm(np.cross(p, q))
    x = X[0]
    for p, t, s in itertools.product(P, T, S):
        if p @ t > cpt - 1e-9 and p @ s > cps - 1e-9 and t @ s > cts - 1e-9 \
                and np.all(np.linalg.solve(np.array([t, p, s]).T, x) >= -1e-9):
            u2, u3 = (p, s) if big == 5 else (s, p)
            d = (abs(x @ plane(u2, u3)), abs(x @ plane(t, u3)), abs(x @ plane(t, u2)))
            break
    if orbit_type in ('tI', 'tO'): return [d[1] / d[2]]
    if orbit_type in ('tD', 'rC'): return [d[0] / d[2]]
    if orbit_type in ('rD', 'tC'): return [d[0] / d[1]]
    return [d[0] / d[2], d[1] / d[2]]


def section_G(rows):
    print('== G. where 7 and 21 occur')
    print('  canon: H(T) is odd, and never 7 (THM-343; THM-338 at n = 5) or 21 (THM-1370), for every n; every other odd')
    print('  value up to 609 occurs by n = 8 (THM-1370).  That 7 and 21 are the ONLY gaps is the open spectrum-completeness')
    print('  conjecture (THM-1370, THM-4094).  Mechanism of the two gaps: multiplicativity over strong components and')
    print('  the strong minima 3, 5, 9, 15, 25, 45, 75, ... (Busch).')
    print('  in the noble classification:')
    print('    genus 21 occurs once, for D-4 (whose dual is the fissary D-F1);')
    print('    genus 7 occurs for', sorted(k for k, r in rows.items() if r['orientable'] and r['genus'] == 7), ';')
    assert not any({r['points'], r['V'], r['E'], r['F']} & {7, 21} for r in rows.values())
    print('    no noble polyhedron has 7 or 21 points, vertices, edges or faces; the fissary field has discriminant')
    print('    -5^2 19, and no element of Z[phi] has norm +-7 or +-21.')
    print('  verdict: NUMEROLOGY (no mechanism links surface genus to Hamiltonian path counts).')


def section_H(lib, tables):
    print("== H. rows of the paper's appendix Tables 5-9 against the model files")
    rows = {}
    for k in (5, 6, 7, 8, 9):
        for r in table(tables, k)[1:]:
            if len(r) != 7: continue
            (pq,) = re.findall(r'(\d+)\s*,\s*(\d+)', r[1])
            rows[r[0]] = dict(pq=tuple(map(int, pq)), V=int(r[2]), E=int(r[3]), F=int(r[4]), dual=r[6])
    assert len(rows) == 146 and set(rows) == {k for k in lib if not k.endswith('-F')}

    def model(name):
        pts, faces = read_off(lib[name])
        V, E, F, chi, ori, g = surface(abstract_structure(faces))
        (p,) = {len(f) for f in faces}
        return ((p, 2 * E // V), V, E, F)

    def printed(name):
        r = rows[name]
        return (r['pq'], r['V'], r['E'], r['F'])

    bad = {k: (printed(k), model(k)) for k in rows if printed(k) != model(k)}
    internal = {k for k in rows if not (2 * rows[k]['E'] == rows[k]['pq'][0] * rows[k]['F'] == rows[k]['pq'][1] * rows[k]['V'])}
    print(f'  {len(rows)} rows parsed; rows disagreeing with their model: {sorted(bad)}; '
          f'internally inconsistent rows (2E = pF = qV fails): {sorted(internal)}')
    assert sorted(bad) == ['D-2', 'D-7', 'rD-5.2', 'tI-5.6'] and internal <= set(bad)
    fmt = lambda t: f'{{{t[0][0]},{t[0][1]}}} V={t[1]} E={t[2]} F={t[3]}'
    for k, (pr, m) in sorted(bad.items()):
        d = rows[k]['dual']
        dual_of_model = ((m[0][1], m[0][0]), m[3], m[2], m[1])
        ok = printed(d) == model(d) == dual_of_model and rows[d]['dual'] == k
        print(f'  {k}: printed {fmt(pr)};  model {fmt(m)};  printed dual row {d}: {fmt(printed(d))} = dual of the model: {ok}')
        assert ok
    print('  every wrong row is refuted by its own printed dual row; the count 146 is unaffected')

    print("== H2. the printed orbit locations of Tables 10-13 against the model files")
    locs, _ = location_rows(tables)
    bad = []
    for orbit, lst in locs.items():
        name = next(k for k in sorted(lib) if k == orbit or k.startswith(orbit + '.'))
        got = orbit_parameters(orbit.split('-')[0], np.array(read_off(lib[name])[0]))
        assert len(got) == len(lst)
        for (v, loc, pol), g in zip(lst, got):
            c = [float(x) for x in pol.all_coeffs()]
            root = abs(np.polyval(c, g)) < 1e-9 * np.polyval(np.abs(c), abs(g))
            if not (root and abs(g - loc) < 1e-8):
                bad.append((orbit, v, loc, g, root))
    print(f'  {len(locs)} orbits: parameters recovered from each model (distances to the mirrors of its Moebius triangle)')
    for orbit, v, loc, g, root in bad:
        print(f'  {orbit} {v}: printed {loc:.14f}, model {g:.14f}, a root of the printed minimal polynomial: {root}')
    assert [(o, v) for o, v, *_ in bad] == [('sD-3', 'b'), ('sD-16', 'b'), ('gD-15', 'a')] and all(r for *_, r in bad)
    print('  three printed decimals are copied from a neighbouring row (sD-3 b = sD-1 b, sD-16 b = sD-17 a,')
    print('  gD-15 a = gD-14 a); the polynomials and the models agree, so nothing else is affected')


if __name__ == '__main__':
    print('== data')
    lib = fetch_library()
    print(f'  model library: {len(lib)} .off files (commit {COMMIT[:7]})')
    tables = fetch_paper_tables()
    rows = section_A(lib)
    section_B(lib)
    section_C(lib)
    section_D(tables)
    section_E(lib)
    section_F(lib)
    section_F2()
    section_F3()
    section_G(rows)
    section_H(lib, tables)
    print('ALL CHECKS PASSED')
