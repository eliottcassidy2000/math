import re, glob
from runs import tour_from_bits
from torsion_vs_hc import ham_cycles_sub, strong_components
res = {}
exc = []
for f in sorted(glob.glob('out/census_n*_k*.out')):
    n, k = map(int, re.findall(r'n(\d+)_k(\d+)', f)[0])
    if k >= n:
        continue
    for l in open(f):
        if not re.match(r'^[01]{3,}\s', l) or ' strong ' not in l:
            continue
        m = re.search(r' M=(.*?)( unglued=| P=|$)', l)
        if not m:
            continue
        bits = l.split()[0]
        A = tour_from_bits(n, bits)
        hc = ham_cycles_sub(A, range(n))
        snf = m.group(1).split()
        tf = all(t.startswith('1^') or t.startswith('0^') for t in snf)
        uni = all(t.startswith('1^') for t in snf)
        d = res.setdefault((n, k), {'HC1': 0, 'HC1_tf': 0, 'HC1_uni': 0, 'HC>=2': 0, 'HC>=2_tf': 0})
        if hc == 1:
            d['HC1'] += 1; d['HC1_tf'] += tf; d['HC1_uni'] += uni
            if not tf: exc.append(('HC1 with torsion', n, k, bits, ' '.join(snf)[:60]))
        else:
            d['HC>=2'] += 1; d['HC>=2_tf'] += tf
            if tf: exc.append(('HC>=2 torsion-free', n, k, bits, hc, ' '.join(snf)[:60]))
for key in sorted(res):
    print('n=%d k=%d strong:' % key, res[key])
for e in exc:
    print('  ', e)
