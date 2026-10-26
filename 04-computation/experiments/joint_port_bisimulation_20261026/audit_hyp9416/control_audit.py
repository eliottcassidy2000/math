# Independent audit of HYP-9416 / note section 3, check E: the 81/27 control.
#   T+ = C3[A+, A+, A+], T- = C3[A+, A+, A-], A+ = C3 -> 1 (triangle dominates the singleton), A- = 1 -> C3.
# Recomputes H (subset DP and full enumeration), |Aut| (leaf count of the audit's own canonical labelling), tc, the
# labelled cyclic triples, injectivity of the partner map on presentations, the shared partner classes, the
# entrance-phase masses and phase families of T-, and (cross-check of vitali FROBENIUS-ABC section 2) the Delta c3
# flux table, Delta c3 = score(first) - score(last) - 1.  Imports only the audit's tcanon.py.
import sys, os, time
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from tcanon import canon, free_key, ham_paths, partner

def child(kind):
    out = [0] * 4
    for u, v in ((0, 1), (1, 2), (2, 0)): out[u] |= 1 << v
    for u in range(3):
        if kind == '+': out[u] |= 1 << 3
        else: out[3] |= 1 << u
    return out
def build(dec):
    out = [0] * 12
    for b, kind in enumerate(dec):
        c = child(kind)
        for u in range(4): out[4 * b + u] |= sum(1 << (4 * b + v) for v in range(4) if (c[u] >> v) & 1)
        for u in range(4):
            for v in range(4): out[4 * b + u] |= 1 << (4 * ((b + 1) % 3) + v)
    for u in range(12): assert not (out[u] >> u) & 1
    for u in range(12):
        for v in range(u + 1, 12): assert ((out[u] >> v) & 1) ^ ((out[v] >> u) & 1)
    return tuple(out)

def hp_dp(out):
    n = len(out); dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for m in range(1, 1 << n):
        row = dp[m]
        for v in range(n):
            c = row[v]
            if not c: continue
            nxt = out[v] & ~m
            while nxt:
                w = nxt & -nxt; nxt ^= w; dp[m | w][w.bit_length() - 1] += c
    return sum(dp[(1 << n) - 1])

def c3(out):
    n = len(out); cnt = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                a, b, c = (out[i] >> j) & 1, (out[j] >> k) & 1, (out[k] >> i) & 1
                if a == b == c: cnt += 1
    return cnt

def pres_code(a, p):
    code = 0; b = 0
    for i in range(11):
        row = a[p[i]]
        for j in range(i + 1, 12):
            if (row >> p[j]) & 1: code |= 1 << b
            b += 1
    return code

t0 = time.time(); R = {}
for name, dec in (('81', '+++'), ('27', '++-')):
    T = build(dec)
    H = hp_dp(T)
    code, aut, order = canon(T)
    P = ham_paths(T); assert len(P) == H
    pres = {}
    mult = Counter()
    for p in P:
        c = pres_code(T, p); mult[c] += 1
        if c not in pres: pres[c] = p
    assert set(mult.values()) == {aut}
    part = {c: free_key(partner(T, p)) for c, p in pres.items()}
    score = [o.bit_count() for o in T]
    flux = Counter(score[p[0]] - score[p[-1]] - 1 for p in P)
    R[name] = dict(T=T, H=H, aut=aut, tc=len(pres), c3=c3(T), part=part, pres=pres, P=P, flux=flux)
    print(f"{name}: H = {H} (DP; enumeration agrees), |Aut| = {aut}, tc = {len(pres)}, cyclic triples = {c3(T)}, "
          f"distinct partner classes = {len(set(part.values()))}  [{time.time() - t0:.0f}s]", flush=True)
a, b = R['81'], R['27']
sa, sb = set(a['part'].values()), set(b['part'].values())
print(f"injective partner map: 81 -> {len(sa) == a['tc']}, 27 -> {len(sb) == b['tc']}; shared partner classes {len(sa & sb)}")
# entrance phases of T- relative to the odd child (block 2)
ph = lambda v: ((v // 4) - 2) % 3
mass = Counter(ph(p[0]) for p in b['P'])
fam = defaultdict(set)
for c, p in b['pres'].items(): fam[ph(p[0])].add(b['part'][c])
print(f"T- entrance-phase path masses {dict(sorted(mass.items()))}; phase families: sizes "
      f"{ {k: len(v) for k, v in sorted(fam.items())} }, pairwise disjoint: "
      f"{all(not (fam[i] & fam[j]) for i in range(3) for j in range(i + 1, 3))}")
# the shared classes by phase family
print(f"shared partner classes by T- entrance phase: { {k: len(v & sa) for k, v in sorted(fam.items())} }")
# exit phases and joint (entrance, exit) phases, for completeness
mass2 = Counter((ph(p[0]), ph(p[-1])) for p in b['P'])
print(f"T- joint (entrance, exit) phase masses: {dict(sorted(mass2.items()))}")
# T+ entrance-block masses (all blocks equivalent by rotation)
massp = Counter(p[0] // 4 for p in a['P'])
print(f"T+ entrance-block masses: {dict(sorted(massp.items()))}")
print(f"Delta c3 flux (labelled paths): 81 {dict(sorted(a['flux'].items()))}; 27 {dict(sorted(b['flux'].items()))}")
print(f"[{time.time() - t0:.0f}s]")
