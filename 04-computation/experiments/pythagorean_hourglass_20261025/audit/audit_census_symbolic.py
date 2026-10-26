# Audit checks C and D, written independently of census_level_forecast.py / census_interpolation.py.
#  C. Survivor census k = 4..7: own canonical forms and automorphism groups (brute force over S_k), the odd-reversal
#     test, completeness of the isomorphism classes (orbit sizes must sum to 2^C(k,2); geng is used for k = 7 only as
#     a source of representatives), signed multiplicities as the SUM OF SIGNS OVER LABELLED COPIES (not the
#     |Aut|^-1 sum_sigma formula), and the Desarmenien-Foata identity for n <= 9.
#  D. An interpolation-free derivation of the census: THM-4652's elimination calculus implemented symbolically
#     (exact integer polynomials in q, q = 3 mod 4), applied to every labelled subgraph of K_k (k <= 6; for k = 7 to
#     the class representatives with their signed multiplicities).  This yields 2^C(k,2) TT_k exactly as
#     A_k(q) + sum_H c_H(q) T_H over the dense survivors H, with no fitting.
# Polynomials are tuples of integer coefficients (low -> high) in q.  Audit 2026-10-25.
import itertools, subprocess, sys, time
from math import comb
from functools import lru_cache
import numpy as np

# ---------------------------------------------------------------- polynomials
def padd(a, b):
    n = max(len(a), len(b)); r = [0] * n
    for i, x in enumerate(a): r[i] += x
    for i, x in enumerate(b): r[i] += x
    while len(r) > 1 and r[-1] == 0: r.pop()
    return tuple(r)
def pscale(a, c):
    r = [c * x for x in a]
    while len(r) > 1 and r[-1] == 0: r.pop()
    return tuple(r) if r else (0,)
def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i + j] += x * y
    while len(r) > 1 and r[-1] == 0: r.pop()
    return tuple(r)
def peval(a, q): return sum(c * q ** i for i, c in enumerate(a))
def pstr(a):
    terms = [f'{c}*q^{i}' for i, c in enumerate(a) if c]
    return ' + '.join(terms) if terms else '0'
def combo_add(A, B, scale=(1,)):
    R = dict(A)
    for k, v in B.items():
        R[k] = padd(R.get(k, (0,)), pmul(v, scale))
        if R[k] == (0,): del R[k]
    return R

# ---------------------------------------------------------------- graphs on k labelled vertices
EDGES = {k: [(i, j) for i in range(k) for j in range(i + 1, k)] for k in range(0, 8)}
EIDX = {k: {e: n for n, e in enumerate(EDGES[k])} for k in EDGES}
PERMS = {}
def perm_tables(k):
    if k in PERMS: return PERMS[k]
    P = list(itertools.permutations(range(k)))
    E = EDGES[k]
    img = np.zeros((len(P), max(1, len(E))), dtype=np.int64)
    rev = np.zeros((len(P), max(1, len(E))), dtype=np.int64)
    for a, s in enumerate(P):
        for n, (i, j) in enumerate(E):
            x, y = s[i], s[j]
            img[a, n] = EIDX[k][(min(x, y), max(x, y))]
            rev[a, n] = 1 if x > y else 0
    PERMS[k] = (P, img, rev)
    return PERMS[k]

def bits(mask): return [n for n in range(mask.bit_length()) if mask >> n & 1]

@lru_cache(maxsize=None)
def canon(k, mask):
    """(canonical mask, sign, odd) with T(mask, order orientation) = sign * T(canonical, order orientation),
    odd = the underlying graph has an automorphism reversing an odd number of order-oriented arcs"""
    if k <= 1 or mask == 0: return mask, 1, False
    P, img, rev = perm_tables(k)
    E = bits(mask)
    pw = (np.int64(1) << img[:, E]).sum(axis=1)
    a = int(np.argmin(pw)); cm = int(pw[a])
    sign = -1 if int(rev[a, E].sum()) % 2 else 1
    auts = np.nonzero(pw == mask)[0]
    odd = bool((rev[np.ix_(auts, E)].sum(axis=1) % 2 == 1).any())
    return cm, sign, odd

def degrees(k, mask):
    d = [0] * k
    for n in bits(mask):
        i, j = EDGES[k][n]; d[i] += 1; d[j] += 1
    return d

# ---------------------------------------------------------------- the elimination calculus (q = 3 mod 4)
def delete_vertex(k, edges_set, v):
    """edges_set: set of (i,j) i<j on 0..k-1 not touching v -> mask on 0..k-2, order-preserving relabel"""
    m = 0
    for (i, j) in edges_set:
        i2 = i - (i > v); j2 = j - (j > v)
        m |= 1 << EIDX[k - 1][(i2, j2)]
    return m

def toggle(S, a, z):
    """multiply by chi(x_z - x_a): toggles the edge {a,z}; sign +1 if a < z else -1 (chi(-1) = -1)"""
    e = (min(a, z), max(a, z))
    S = set(S)
    if e in S: S.remove(e)
    else: S.add(e)
    return S, (1 if a < z else -1)

@lru_cache(maxsize=None)
def R(k, mask):
    """T(G) for the canonical mask with order orientation, as {key: poly}; key 'one' = the empty graph on 0 vertices,
    key ('dense', k, mask) = a dense survivor (left symbolic)"""
    if bin(mask).count('1') % 2: return {}
    cm, s, odd = canon(k, mask)
    assert cm == mask and s == 1
    if odd: return {}
    if k == 0: return {'one': (1,)}
    deg = degrees(k, mask)
    v = min(range(k), key=lambda x: deg[x])
    if deg[v] >= 3: return {('dense', k, mask): (1,)}
    S = {EDGES[k][n] for n in bits(mask)}
    nbrs = [u for u in range(k) if (min(u, v), max(u, v)) in S]
    S0 = {e for e in S if v not in e}
    def T_of(Sx, sign):
        m = delete_vertex(k, Sx, v)
        cm2, s2, odd2 = canon(k - 1, m)
        if odd2 or sign == 0: return {}
        return {key: pscale(val, sign * s2) for key, val in R(k - 1, cm2).items()}
    out = {}
    others = [z for z in range(k) if z != v]
    if deg[v] == 0:
        out = combo_add(out, T_of(S0, 1), (-(k - 1), 1))                      # (q - (k-1)) T(G - v)
    elif deg[v] == 1:
        u = nbrs[0]; su = 1 if u < v else -1
        for z in others:
            if z == u: continue
            S1, t = toggle(S0, u, z)
            out = combo_add(out, T_of(S1, -su * t))
    else:
        u, w = nbrs; su = 1 if u < v else -1; sw = 1 if w < v else -1
        out = combo_add(out, T_of(S0, -su * sw))                               # -1 (Jacobi sum)
        for z in others:
            if z in (u, w): continue
            S1, t1 = toggle(S0, u, z); S2, t2 = toggle(S1, w, z)
            out = combo_add(out, T_of(S2, -su * sw * t1 * t2))
    return out

def T_oriented(k, mask):
    cm, s, odd = canon(k, mask)
    if odd: return {}
    return {key: pscale(v, s) for key, v in R(k, cm).items()}

# ---------------------------------------------------------------- C: classes, survivors, multiplicities
def g6_decode(sline):
    n = ord(sline[0]) - 63
    bl = []
    for ch in sline[1:]:
        v = ord(ch) - 63; bl.extend((v >> (5 - i)) & 1 for i in range(6))
    m = 0; t = 0
    for j in range(1, n):
        for i in range(j):
            if bl[t]: m |= 1 << EIDX[n][(i, j)]
            t += 1
    return n, m

def orbit_signs(k, mask):
    """dict image-mask -> set of signs over all relabellings: the labelled copies of the graph and their signs
    T(copy, order) = sign * T(mask, order)"""
    P, img, rev = perm_tables(k)
    E = bits(mask)
    pw = (np.int64(1) << img[:, E]).sum(axis=1)
    par = rev[:, E].sum(axis=1) % 2
    d = {}
    for m, r in zip(pw.tolist(), par.tolist()):
        d.setdefault(m, set()).add(-1 if r else 1)
    return d

def class_reps(k):
    if k <= 6:
        reps = {}
        for m in range(1 << comb(k, 2)):
            cm = canon(k, m)[0]
            reps.setdefault(cm, 0); reps[cm] += 1
        return sorted(reps)
    outp = subprocess.run(['geng', '-q', str(k)], capture_output=True, text=True).stdout.split()
    return [canon(k, g6_decode(s)[1])[0] for s in outp]


@lru_cache(maxsize=None)
def class_info(k):
    """(reps, info): info[cm] = (number of labelled copies, m_H = sum of signs over labelled copies, odd)"""
    reps = class_reps(k); info = {}
    for cm in reps:
        osg = orbit_signs(k, cm)
        consistent = all(len(s) == 1 for s in osg.values())
        odd = canon(k, cm)[2]
        assert consistent == (not odd)
        info[cm] = (len(osg), sum(next(iter(s)) for s in osg.values()) if consistent else 0, odd)
    return reps, info

def mask_of(k, edges): return sum(1 << EIDX[k][(min(i, j), max(i, j))] for i, j in edges)
W4_doc = mask_of(5, [(0, 1), (0, 2), (0, 3), (0, 4), (1, 3), (1, 4), (2, 3), (2, 4)])
OCT_doc = mask_of(6, [e for e in EDGES[6] if e not in [(0, 1), (2, 3), (4, 5)]])
KP_doc = mask_of(6, [e for e in EDGES[6] if e not in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]])

def sgn(s):
    s = list(s); n = len(s); seen = [False] * n; c = 0
    for i in range(n):
        if not seen[i]:
            c += 1; j = i
            while not seen[j]: seen[j] = True; j = s[j]
    return -1 if (n - c) % 2 else 1
def eulerian(m):
    A = [0] * max(1, m)
    for s in itertools.permutations(range(m)):
        A[sum(1 for i in range(m - 1) if s[i] > s[i + 1])] += 1
    return tuple(A)

# ---------------------------------------------------------------- D: symbolic census
def census(k):
    if k <= 6:
        tot = {}
        for m in range(1 << comb(k, 2)):
            tot = combo_add(tot, T_oriented(k, m))
        return tot
    reps, info = class_info(k)
    tot = {}
    for cm in reps:
        n_copies, m_H, odd = info[cm]
        if odd or m_H == 0: continue
        tot = combo_add(tot, {key: pscale(v, m_H) for key, v in R(k, cm).items()})
    return tot

def express(res, k):
    """translate dense keys into the documents' labellings"""
    lab = {}
    for key, v in res.items():
        if key == 'one': lab['1'] = v; continue
        _, kk, m = key
        name = None
        for nm, k2, m2 in [('T_W4', 5, W4_doc), ('T_oct', 6, OCT_doc), ('T_K6-P6', 6, KP_doc)]:
            if kk == k2:
                c2, s2, _ = canon(k2, m2)
                if c2 == m: name = nm; v = pscale(v, s2)        # T(canonical) = s2 * T(doc labelling)
        lab[name or f'T[{kk}-vertex dense {m}]'] = v
    return lab

if __name__ == '__main__':
    lines = []
    def say(s=''):
        print(s, flush=True); lines.append(s)
    t0 = time.time()
    for k in range(4, 8):
        reps, info = class_info(k)
        tot_labelled = sum(info[cm][0] for cm in reps)
        dense = [cm for cm in reps if min(degrees(k, cm)) >= 3]
        even = [cm for cm in dense if bin(cm).count('1') % 2 == 0]
        killed = [cm for cm in even if info[cm][2]]
        surv = [cm for cm in even if not info[cm][2]]
        say(f'k = {k}: {len(reps)} isomorphism classes, orbit sizes sum to {tot_labelled} = 2^{comb(k, 2)}: {tot_labelled == 1 << comb(k, 2)}; '
            f'dense (min deg >= 3): {len(dense)}; even |E|: {len(even)}; killed by an odd-reversal automorphism: {len(killed)}; survivors: {len(surv)}')
        for cm in surv:
            deg = sorted(degrees(k, cm)); P, img, rev = perm_tables(k)
            E = bits(cm); pw = (np.int64(1) << img[:, E]).sum(axis=1); aut = int((pw == cm).sum())
            alleven = all(x % 2 == 0 for x in deg)
            comp = [e for e in EDGES[k] if not cm >> EIDX[k][e] & 1]
            say(f'    |E| = {len(E):2d} degrees {deg} |Aut| = {aut:4d} m (sum of signs over labelled copies, canonical labelling) = {info[cm][1]:+d}'
                f'{"  [all degrees even]" if alleven else ""}  complement {comp}')
        if k == 7:
            ms = sorted({abs(info[cm][1]) for cm in surv})
            sg = sorted({info[cm][1] for cm in surv})
            say(f'    k = 7: |m| values {ms}; signed values (canonical labelling) {sg}')
    say(f'[time {time.time() - t0:.0f}s]')
    for name, k, m in [('W4 (hub 0, non-edges {1,2},{3,4})', 5, W4_doc), ('octahedron (antipodal {1,2},{3,4},{5,6})', 6, OCT_doc),
                       ('K6-P6 (path 1-2-...-6)', 6, KP_doc)]:
        osg = orbit_signs(k, m)
        say(f'signed multiplicity of {name} in K{k} with that labelling: {sum(next(iter(s)) for s in osg.values()):+d} '
            f'({len(osg)} labelled copies)')
    okDF = True
    for n in range(1, 10):
        L = [0] * n
        for s in itertools.permutations(range(n)):
            L[sum(1 for i in range(n - 1) if s[i] > s[i + 1])] += sgn(s)
        rhs = eulerian((n + 1) // 2)
        for _ in range(n // 2): rhs = pmul(rhs, (1, -1))
        okDF &= padd(tuple(L), pscale(rhs, -1)) == (0,)
        if n == 6: val6 = sum(c * (-1) ** i for i, c in enumerate(L))
    say(f'Desarmenien-Foata sum_(S_n) sgn t^des = (1-t)^floor(n/2) A_ceil(n/2)(t) for n = 1..9: {okDF}; at n = 6, t = -1: {val6} (half = {val6 // 2})')
    bad = 0
    Ekp = [e for e in EDGES[6] if e not in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]]
    for s in itertools.permutations(range(6)):
        revc = sum(1 for i, j in Ekp if s[i] > s[j])
        inv = sum(1 for i in range(6) for j in range(i + 1, 6) if s[i] > s[j]); des = sum(1 for i in range(5) if s[i] > s[i + 1])
        bad += revc != inv - des
    say(f'rev_(K6-P6)(sigma) = inv(sigma) - des(sigma) for all 720 sigma: {bad == 0}')
    t1 = time.time()
    for k in range(2, 7):
        lab = express(census(k), k)
        say(f'2^{comb(k, 2)} TT{k} = ' + ' + '.join(f'({pstr(v)})*{nm}' for nm, v in lab.items()))
        if k == 5:
            A5 = lab['1']; tgt = pmul(pmul((0, 1), (-1, 1)), pmul((-3, 1), (105, -21, 1)))
            say(f'   A5 = q(q-1)(q-3)(q^2-21q+105): {A5 == tgt};  coefficient of T_W4 = {lab.get("T_W4")}')
        if k == 6:
            A6 = lab['1']; tgt = pmul(pmul((0, 1), (-1, 1)), pmul((-3, 1), (-3239, 682, -46, 1)))
            B = padd(lab.get('T_W4', (0,)), pmul(lab.get('T_oct', (0,)), (-5, 1)))   # substitute T_oct = (q-5) T_W4
            say(f'   A6 = q(q-1)(q-3)(q^3-46q^2+682q-3239): {A6 == tgt}')
            say(f'   coefficient of T_W4 from non-dense graphs: {pstr(lab.get("T_W4", (0,)))}; of T_oct: {pstr(lab.get("T_oct", (0,)))}; '
                f'total with T_oct = (q-5)T_W4: {pstr(B)} (= 7(q-15)? {B == (-105, 7)}); coefficient of T_K6-P6: {pstr(lab.get("T_K6-P6", (0,)))}')
    say(f'[time {time.time() - t1:.0f}s]')
    t2 = time.time()
    lab7 = express(census(7), 7)
    dense7 = {nm: v for nm, v in lab7.items() if nm.startswith('T[7')}
    say(f'2^21 TT7: A7 = {pstr(lab7["1"])} (degree {len(lab7["1"]) - 1})')
    say(f'   coefficient of T_W4: {pstr(lab7.get("T_W4", (0,)))} (degree {len(lab7.get("T_W4", (0,))) - 1}); '
        f'of T_K6-P6: {pstr(lab7.get("T_K6-P6", (0,)))}; of T_oct: {pstr(lab7.get("T_oct", (0,)))}')
    say(f'   seven-vertex dense survivors with nonzero coefficient: {len(dense7)}; all coefficients constant: '
        f'{all(len(v) == 1 for v in dense7.values())}')
    say(f'[time {time.time() - t2:.0f}s]')
    with open('audit_census_symbolic.out', 'w') as fh: fh.write('\n'.join(lines) + '\n')
