"""Certify one representative of every distinct two-clock parameter tuple found in the censuses with the rigorous
complete-table program c/certify.c, recompute the tuple from the table, spot-check literal counts, and list the
unit-offset families H = 1 + (alpha 3^j - sig2) 2^k (beta = 0, c11 = 1) with their seed multipliers."""
import gzip, re, subprocess, sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from tern_check import parse
from certify_twoclock import substitute, TT
from funnel import H_only
def lines(name):
    with gzip.open(os.path.join(ROOT, 'data', name), 'rt') as f: return f.read().splitlines()
reps = {}
def add(key, n, adj, u, v):
    if key not in reps: reps[key] = (n, adj, u, v)
def tup(l):
    g = lambda pat: int(re.search(pat, l).group(1))
    return (g(r' A=(\d+)'), g(r'sigma=(-?\d+)'), g(r' al=(\d+)'), g(r'sig2=(-?\d+)'))
for l in lines('gen_all_orders.txt.gz'):
    if not l.startswith('GEN'): continue
    p = l.split(); n = int(p[1]); code = p[2]; u = int(p[3].split('=')[1]); v = int(p[4].split('=')[1])
    add(tup(l), n, parse(code, n), u, v)
for f in ('backscan_n8_b4.txt.gz', 'backscan_n10_b5.txt.gz', 'backscan_n11_b4.txt.gz'):
    for l in lines(f):
        if not l.startswith('GEN'): continue
        n = int(re.search(r'n=(\d+)', l).group(1)); u = int(re.search(r' u=(\d+)', l).group(1)); v = int(re.search(r' v=(\d+)', l).group(1))
        bk = set(tuple(map(int, x.split('>'))) for x in re.search(r'back=([\d>,]+)', l).group(1).strip(',').split(','))
        adj = [0]*n
        for i in range(n):
            for j in range(i+1, n):
                if (j, i) in bk: adj[j] |= 1 << i
                else: adj[i] |= 1 << j
        add(tup(l), n, adj, u, v)
S60 = [2+4+8+64+128, 4+8+16+32+64+128, 8+16+32, 16+32+64+128, 1+32+64+128, 1+64+128, 4+128, 4]
def s60_ext(pat, e):
    n = 8 + e; adj = list(S60) + [0]*e
    for k in range(e):
        for j in range(8):
            if (pat >> (8*k + j)) & 1: adj[8+k] |= 1 << j
            else: adj[j] |= 1 << (8+k)
    if e == 2:
        if (pat >> 16) & 1: adj[8] |= 1 << 9
        else: adj[9] |= 1 << 8
    return adj, n
for f in ('extscan_e1.txt.gz', 'extscan_e2.txt.gz'):
    for l in lines(f):
        if not l.startswith('GEN'): continue
        e = int(re.search(r'e=(\d+)', l).group(1)); pat = int(re.search(r'pat=(\d+)', l).group(1)); v = int(re.search(r' v=(\d+)', l).group(1))
        adj, n = s60_ext(pat, e)
        add(tup(l), n, adj, 3, v)
for l in lines('extone_from_certified.txt.gz'):
    if not l.startswith('GEN'): continue
    n = int(re.search(r'n=(\d+)', l).group(1)); u = int(re.search(r' u=(\d+)', l).group(1)); v = int(re.search(r' v=(\d+)', l).group(1))
    rows = list(map(int, re.search(r'rows=([\d,]+)', l).group(1).split(',')))
    add(tup(l), n, rows, u, v)
def main(certifier):
  keys = sorted(reps)
  inp = ''.join('%d %d %d %s\n' % (reps[k][0], reps[k][2], reps[k][3], ' '.join(map(str, reps[k][1]))) for k in keys)
  res = subprocess.run([certifier], input=inp, capture_output=True, text=True, check=True).stdout.splitlines()
  if len(res) != len(keys): raise RuntimeError('certifier output mismatch')
  unit = []; certified = 0
  for k, line in zip(keys, res):
    tab = {}
    for item in line.split()[1:]:
        ij, val = item.split(':'); i, ip = map(int, ij.split(',')); tab[(i, ip)] = int(val)
    if any(c for (i, ip), c in tab.items() if i >= 4 or ip >= 3): raise RuntimeError('higher visits in %s' % (k,))
    A, al = tab.get((3, 1), 0), tab.get((3, 2), 0)
    if tab.get((2, 1)) != 3*A or tab.get((2, 2)) != 3*al: raise RuntimeError('not ternary-shaped %s' % (k,))
    t = (A, 3*A - tab[(1, 1)], al, 3*al - tab.get((1, 2), 0))
    if t != k: raise RuntimeError('tuple mismatch %s %s' % (k, t))
    n, adj, u, v = reps[k]
    for (j, kk) in ((1, 1), (2, 1), (1, 2)):
        B, tot = substitute(adj, n, u, v, TT(j), TT(kk))
        A_, s, al_, s2 = t
        pred = 3**j*(al_*2**kk + A_ - 2*al_) - s2*2**kk - (s - 2*s2)
        if H_only(B, tot) != pred: raise RuntimeError('literal mismatch %s' % (k,))
    certified += 1
    A_, s, al_, s2 = t
    if A_ == 2*al_ and s - 2*s2 == -1:
        unit.append({'alpha': al_, 'sig2': s2, 'n': n, 'u': u, 'v': v, 'rows': adj})
  def odd(x):
    while x % 2 == 0: x //= 2
    return x
  seeded = sorted({odd(f['alpha']*3**j - f['sig2']) for f in unit for j in range(1, 40) if f['alpha']*3**j - f['sig2'] > 0})
  print('certified tuples:', certified, ' unit-offset families:', len(unit))
  print('distinct odd seed multipliers c <= 10^6:', len([c for c in seeded if c <= 10**6]),
      ' first 40:', [c for c in seeded][:40])
  json.dump({'certified_tuples': certified, 'unit_families': unit}, open(os.path.join(ROOT, 'data', 'unit_two_clock_families.json'), 'w'), indent=0)
  allrep = [{'A': k[0], 'sigma': k[1], 'alpha': k[2], 'sig2': k[3], 'n': reps[k][0], 'u': reps[k][2], 'v': reps[k][3], 'rows': reps[k][1]} for k in keys]
  json.dump(allrep, open(os.path.join(ROOT, 'data', 'two_clock_representatives.json'), 'w'), indent=0)

  return certified, unit, seeded

if __name__ == '__main__':
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        exe = os.path.join(d, 'certify')
        subprocess.run(['cc', '-O2', '-o', exe, os.path.join(ROOT, 'c', 'certify.c')], check=True)
        main(exe)
