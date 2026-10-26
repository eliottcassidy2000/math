# Independent audit of HYP-9416 (3), check D: joint two-port response histograms of strong atoms (rooted-lines R1)
#   F_A(h) = #presentations x of A (Hamilton paths modulo Aut A) with [tau x; first, last] = h,
#   B_A(h) = the blue ones (the reversal P_i -> P_(n-1-i) is an isomorphism A -> A^op).
# Distinctness of F_A, also after forgetting ports, plus coarser invariants (to measure how much of the histogram is
# needed), totals, and an exact coefficient-by-coefficient comparison with the vitali-response strong-atlas.json.
# Imports only the audit's tcanon.py.  Usage: atoms.py n
import sys, os, json, time
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from tcanon import canon, rooted_key, onesided_key, free_key, all_tournaments, is_strong, ham_paths, partner, relabel
import subprocess
from math import factorial

n = int(sys.argv[1]); t0 = time.time()
def classes_of(n):
    if n <= 6:
        seen = {}
        for out in all_tournaments(n):
            c, aut, order = canon(out)
            if c not in seen: seen[c] = relabel(out, order)
        return list(seen.values())
    reps = []
    for s in subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split():
        out = [0] * n; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if s[k] == '1': out[i] |= 1 << j
                else: out[j] |= 1 << i
                k += 1
        reps.append(tuple(out))
    return reps
atoms = [a for a in classes_of(n) if is_strong(a)]

def pres_code(a, p):
    m = len(p); code = 0; b = 0
    for i in range(m - 1):
        row = a[p[i]]
        for j in range(i + 1, m):
            if (row >> p[j]) & 1: code |= 1 << b
            b += 1
    return code

def blue(a, p):
    m = len(p)
    return all(((a[p[i]] >> p[j]) & 1) == ((a[p[m - 1 - j]] >> p[m - 1 - i]) & 1) for i in range(m) for j in range(m) if i != j)

info = {}
for a in atoms:
    code, aut, order = canon(a)
    P = ham_paths(a)
    pres = Counter(pres_code(a, p) for p in P)
    assert set(pres.values()) == {aut}, "Aut acts freely on paths and |Aut| = H/tc"
    F = Counter(); Bc = Counter(); Ffree = Counter(); Fside = Counter(); Fscore = Counter()
    for p in P:
        K = partner(a, p)
        h = rooted_key(K, p[0], p[-1])
        F[h] += 1; Bc[h] += blue(a, p)
    for h in F: assert F[h] % aut == 0 and Bc[h] % aut == 0
    F = Counter({h: c // aut for h, c in F.items()}); Bc = Counter({h: c // aut for h, c in Bc.items() if c})
    info[code] = dict(rep=a, aut=aut, H=len(P), tc=len(pres), F=F, B=Bc)
# free classes, one-sided classes and score sequences of all targets
targets = sorted({h for d in info.values() for h in d['F']})
from tcanon import from_code
trep = {}
for h in targets:
    nn, tag, c = h
    out = from_code(nn, c)          # canonical order puts the entrance at position 0 and the exit at position 1
    trep[h] = (out, 0, 1)
tfree = {h: free_key(trep[h][0]) for h in targets}
tside = {h: (onesided_key(trep[h][0], 0), onesided_key(trep[h][0], 1)) for h in targets}
tscore = {h: tuple(sorted(o.bit_count() for o in trep[h][0])) for h in targets}
tstrong = {h: is_strong(trep[h][0]) for h in targets}

def distinct(f):
    g = defaultdict(list)
    for code, d in info.items(): g[f(d)].append(code)
    return len(g), sorted(len(v) for v in g.values() if len(v) > 1)
fz = lambda c: frozenset(c.items())
res = {
    'atoms': len(info), 'presentations': sum(d['tc'] for d in info.values()),
    'largest_fibre': max(max(d['F'].values()) for d in info.values()),
    'distinct_joint_targets': len(targets), 'strong_targets': sum(tstrong.values()),
    'atoms_with_blue': sum(1 for d in info.values() if d['B']),
    'F_distinct': distinct(lambda d: fz(d['F'])),
    'F_portsforgotten_distinct': distinct(lambda d: fz(Counter({tfree[h]: 0 for h in d['F']}) + sum((Counter({tfree[h]: m}) for h, m in d['F'].items()), Counter()))),
    'F_onesided_distinct': distinct(lambda d: fz(sum((Counter({tside[h]: m}) for h, m in d['F'].items()), Counter()))),
    'support_free_set_distinct': distinct(lambda d: frozenset(tfree[h] for h in d['F'])),
    'target_score_multiset_distinct': distinct(lambda d: fz(sum((Counter({tscore[h]: m}) for h, m in d['F'].items()), Counter()))),
    'H_tc_Aut_score_distinct': distinct(lambda d: (d['H'], d['tc'], d['aut'], tuple(sorted(o.bit_count() for o in d['rep']))))
}
print(f"n = {n}: " + json.dumps(res) + f"  [{time.time() - t0:.0f}s]", flush=True)

# exact comparison with the vitali-response atlas
atlas_path = os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '05-knowledge', 'results',
                          'strong_atom_continuation_vitali_20261026', 'strong-atlas.json')
atlas = [e for e in json.load(open(atlas_path)) if e['n'] == n]
mism = 0; matched = 0
for e in atlas:
    src = tuple(e['source']); code = canon(src)[0]
    d = info[code]
    if e['Aut'] != d['aut'] or e['H'] != d['H'] or e['presentations'] != d['tc']: mism += 1; continue
    F = Counter(); Bc = Counter()
    for r in e['response']:
        h = rooted_key(tuple(r['rows']), r['roots'][0], r['roots'][1])
        F[h] += r['total']
        if r['blue']: Bc[h] += r['blue']
    if F == d['F'] and Bc == d['B']: matched += 1
    else: mism += 1
print(f"  atlas: {len(atlas)} sources, {matched} with identical (Aut, H, tc, F_A, B_A), {mism} mismatches")
