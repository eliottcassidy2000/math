# The 81/27 hostile control.  T+ = C3[A+,A+,A+] and T- = C3[A+,A+,A-], A+ = C3->1, A- = 1->C3 (12 vertices).
# They share H and the labelled cyclic-triple data but |Aut| = 81 vs 27.  For every Hamilton path P (presentation)
# compute the line partner K_P(T) = T^op xor E(P), its free class (nauty), and the attachment phase: the blocks of
# the entrance and exit relative to the odd child (the Frobenius phase of the period-3 decoration).  Compare the
# partner histograms at several resolutions.  mac-mini-2026-10-25-pyth (2026-10-26).
import time
from collections import Counter, defaultdict
from rooted import ham_paths, relabel
from canon_nauty import canon_batch

def child(kind):
    # vertices 0,1,2 = the 3-cycle 0->1->2->0, vertex 3 = the singleton; A+ : cycle -> singleton, A- : singleton -> cycle
    adj = [0] * 4
    for u, v in ((0, 1), (1, 2), (2, 0)): adj[u] |= 1 << v
    for u in range(3):
        if kind == '+': adj[u] |= 1 << 3
        else: adj[3] |= 1 << u
    return adj
def substitute(dec):
    n = 12; adj = [0] * n
    blocks = [list(range(4 * b, 4 * b + 4)) for b in range(3)]
    for b, kind in enumerate(dec):
        c = child(kind)
        for u in range(4):
            for v in range(4):
                if (c[u] >> v) & 1: adj[4 * b + u] |= 1 << (4 * b + v)
    for b in range(3):                       # quotient C3: block b -> block b+1
        for u in blocks[b]:
            for v in blocks[(b + 1) % 3]: adj[u] |= 1 << v
    return tuple(adj)
def partner(adj, P):
    n = len(adj)
    op = [0] * n
    for u in range(n):
        for v in range(n):
            if (adj[v] >> u) & 1: op[u] |= 1 << v
    for a, b in zip(P, P[1:]):
        op[a] |= 1 << b; op[b] &= ~(1 << a)
    return tuple(op)

t0 = time.time()
res = {}
for name, dec in (('81 (+++)', '+++'), ('27 (++-)', '++-')):
    T = substitute(dec)
    P = ham_paths(T)
    # presentations and |Aut|
    pres = Counter()
    for p in P:
        perm = [0] * 12
        for i, v in enumerate(p): perm[v] = i
        pres[relabel(T, perm)] += 1
    aut = len(P) // len(pres)
    K = [partner(T, p) for p in P]
    canon = canon_batch(K)
    odd = 2 if dec == '++-' else None
    def phase(v):
        b = v // 4
        return (b - 2) % 3 if odd is not None else 'sym'
    hist_free = Counter(canon)
    hist_phase = Counter((c, phase(p[0]), phase(p[-1])) for c, p in zip(canon, P))
    ent_phase = Counter(phase(p[0]) for p in P)
    res[name] = dict(T=T, H=len(P), aut=aut, tc=len(pres), free=hist_free, phase=hist_phase, ent=ent_phase)
    print(f"{name}: H = {len(P)}, |Aut| = {aut}, tc = {len(pres)}, distinct partner classes {len(hist_free)}; "
          f"entrance-block phase counts {dict(ent_phase)}  [{time.time() - t0:.0f}s]", flush=True)
a, b = res['81 (+++)'], res['27 (++-)']
# coarse invariants agree?
def c3(adj):
    n = len(adj); return sum(1 for i in range(n) for j in range(n) for k in range(n) if i < j < k and (
        ((adj[i] >> j) & 1 and (adj[j] >> k) & 1 and (adj[k] >> i) & 1) or ((adj[i] >> k) & 1 and (adj[k] >> j) & 1 and (adj[j] >> i) & 1)))
print(f"cyclic triples: {c3(a['T'])} vs {c3(b['T'])}; H: {a['H']} vs {b['H']}")
print(f"partner-class histograms (labelled paths) equal: {a['free'] == b['free']}")
print(f"  partner classes common to both: {len(set(a['free']) & set(b['free']))}; only in 81: {len(set(a['free']) - set(b['free']))}; only in 27: {len(set(b['free']) - set(a['free']))}")
# per presentation (divide by |Aut|): the incidence histograms
pa = {c: m // a['aut'] for c, m in a['free'].items()}; pb = {c: m // b['aut'] for c, m in b['free'].items()}
print(f"presentation-level partner histograms equal: {pa == pb}; sums {sum(pa.values())} vs {sum(pb.values())}")
# Frobenius phase: in the 27 case, does the partner histogram depend on the entrance phase?
byph = defaultdict(Counter)
for (c, pe, px), m in b['phase'].items(): byph[pe][c] += m
print("27 case: partner histograms by entrance phase (0 = the odd child's block):",
      {ph: (sum(h.values()), len(h)) for ph, h in sorted(byph.items())},
      "; phases 1 and 2 equal:", byph[1] == byph[2], "; phase 0 equals phase 1:", byph[0] == byph[1])
