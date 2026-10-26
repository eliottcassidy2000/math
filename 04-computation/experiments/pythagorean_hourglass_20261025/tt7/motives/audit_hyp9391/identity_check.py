# Audit of HYP-9391 (A, B): the auditor's brute-force seven-vertex sums against the lane's data, and every one of the
# 37 stated identities S_H = (known combination) + (+-Y) of exact_tt7_formula.out evaluated with the auditor's own
# S_H (brute_sums.c), Y sums (ysums.c) and PARI coefficients.  Also sum_H m_H S_H = g(q) with g from the auditor's
# TT7 counts and the reduced form of the earlier audit.
# Usage: identity_check.py brute_file [brute_file ...] -- ysums_file T5file formsfile ; writes identity_check.out
import sys, os, re, json
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
args = sys.argv[1:]
k = args.index('--')
bfiles, (Yf, T5f, Ff) = args[:k], args[k + 1:]
G = json.load(open(os.path.join(HERE, '..', 'genuine37.json')))
mine = {}
for fn in bfiles:
    for line in open(fn):
        r = list(map(int, line.split()))
        if r: mine[r[0]] = r[1:]
lane = {}
for line in open(os.path.join(HERE, '..', 'sums_7_199.jsonl')):
    d = json.loads(line); lane.setdefault(d['q'], {}).update({i: v for i, v in enumerate(d['S'])})
for fn in sorted(os.listdir(os.path.join(HERE, '..'))):
    if fn.startswith('fast_') and fn.endswith('.jsonl'):
        for line in open(os.path.join(HERE, '..', fn)):
            d = json.loads(line); lane.setdefault(d['q'], {}).update({int(i): v for i, v in d['S'].items()})
forms = {}
for line in open(Ff):
    p, x84, x86, x46, nb = map(int, line.split()); forms[p] = {'a84': x84, 'a86': x86, 'a46': x46, '|b|^2': nb}
Ys = {}
for line in open(Yf):
    r = list(map(int, line.split())); Ys[r[0]] = {'Y02': r[1], 'Y05': r[3], 'Y09': r[5], 'Y18': r[6]}
T5 = {}
for line in open(T5f):
    p, t = map(int, line.split()); T5[p] = t
ident = {}
for line in open(os.path.join(HERE, '..', 'exact_tt7_formula.out')):
    m = re.match(r'H(\d\d) \(m = ([+-]\d+), (\d+) primes\): (.*)$', line.strip())
    if not m: continue
    terms = re.findall(r'\(([-\d/]+)\)\*(\S+(?: a84)?)', m.group(4))
    assert len(terms) == m.group(4).count('*'), line
    ident[int(m.group(1))] = [(Fraction(cf), nm) for cf, nm in terms]
def col(nm, p):
    if nm == '1': return 1
    if nm == 'q': return p
    if nm == 'q^2': return p * p
    if nm == 'q^3': return p ** 3
    if nm == 'q a84': return p * forms[p]['a84']
    if nm in forms[p]: return forms[p][nm]
    return Ys[p][nm]
out = []
out.append(f"auditor's brute-force sums at {len(mine)} primes: {sorted(mine)}")
# 1. against the lane's data
cmp_lines = []
for p in sorted(mine):
    if p not in lane: cmp_lines.append(f"  q={p}: no lane data"); continue
    idx = sorted(lane[p])
    bad = [i for i in idx if lane[p][i] != mine[p][i]]
    cmp_lines.append(f"  q={p}: lane has {len(idx)} sums; equal: {len(idx) - len(bad)}{'; DIFFER at ' + str(bad) if bad else ''}")
out.append("1. auditor's S_H vs lane (sums_7_199.jsonl, fast_*.jsonl):")
out += cmp_lines
# 2. identities
out.append(f"2. stated identities S_H = known + (+-Y) (exact_tt7_formula.out, {len(ident)} graphs) at the auditor's primes:")
allok = True
for i in sorted(ident):
    res = []
    for p in sorted(mine):
        if p not in Ys or p not in forms: continue
        rhs = sum(cf * col(nm, p) for cf, nm in ident[i])
        res.append(rhs == mine[p][i])
    ok = all(res); allok &= ok
    ytag = [nm for _, nm in ident[i] if nm.startswith('Y')]
    out.append(f"  H{i:02d} m={G[i]['m']:+4d} deg={sorted(G[i]['degrees'])} Y={ytag or '-'}: holds at {sum(res)}/{len(res)} primes")
out.append(f"  ALL 37 identities hold at every auditor prime: {allok}")
# 3. sum m_H S_H = g(q) with g from the auditor's TT7 counts
gl = []
for p in sorted(mine):
    if p not in T5 or p not in forms: continue
    TT7 = p * (p - 1) // 2 * T5[p]
    A7 = p * (p - 1) * (p - 3) * (p ** 4 - 87 * p ** 3 + 2785 * p ** 2 - 38827 * p + 198842)
    B7 = 28 * p * p - 1288 * p + 13238; C7 = 1912 - 64 * p
    TW4 = p * (p - 1) * (-4 - forms[p]['a84']); TK = p * (p - 1) * (forms[p]['|b|^2'] - p * p - 7)
    num = (1 << 21) * TT7 - A7 - B7 * TW4 - C7 * TK
    assert num % (p * (p - 1)) == 0
    g = num // (p * (p - 1))
    comb = sum(H['m'] * mine[p][i] for i, H in enumerate(G) if not H['all_even'])
    # K7 - C7 (index of the all-even graph) against eps (q-6) S_(K6-P6), eps = +1 (audit)
    ia = next(i for i, H in enumerate(G) if H['all_even'])
    k7 = mine[p][ia] == (p - 6) * (forms[p]['|b|^2'] - p * p - 7)
    gl.append(f"  q={p}: sum m_H S_H = {comb}, g(q) from auditor TT7 = {g}: {comb == g}; S_(K7-C7) = (q-6) S_(K6-P6): {k7}")
out.append("3. multiplicities: sum over the 37 genuine m_H S_H = g(q), and the K7 - C7 reduction:")
out += gl
print('\n'.join(out))
open(os.path.join(HERE, 'identity_check.out'), 'w').write('\n'.join(out) + '\n')
