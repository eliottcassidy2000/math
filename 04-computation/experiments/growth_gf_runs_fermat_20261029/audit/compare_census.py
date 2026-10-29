"""Compare the audit census (audit/out/census_n*_k*.out) with the producer's census files class by class."""
import re, sys, glob
def parse_prod(path):
    d = {}
    for l in open(path):
        if not re.match(r'^[01]{3,}\s', l): continue
        bits = l.split()[0]
        m = re.search(r' M=(.*?)( unglued=| P=|$)', l)
        u = re.search(r' unglued=(.*?) spur2', l)
        sp = re.search(r'spur2=(\d+) spur3=(\d+)', l)
        kb = re.search(r'killsBdry=(\w+)', l)
        d[bits] = (m.group(1).strip() if m else None, u.group(1).strip() if u else None,
                   (int(sp.group(1)), int(sp.group(2))) if sp else None, kb.group(1) == 'True' if kb else None)
    return d
def parse_mine(path):
    d = {}
    for l in open(path):
        if not re.match(r'^[01]{3,}\s', l): continue
        bits = l.split()[0]
        m = re.search(r' M=(.*?)( N=|$)', l)
        u = re.search(r' N=(.*?) spur2', l)
        sp = re.search(r'spur2=(-?\d+) spur3=(-?\d+)', l)
        kb = re.search(r'kills=(\w+)', l)
        hc = int(re.search(r'HC=(\d+)', l).group(1))
        d[bits] = (m.group(1).strip(), u.group(1).strip() if u else None,
                   (int(sp.group(1)), int(sp.group(2))) if sp else None, kb.group(1) == 'True' if kb else None, hc)
    return d
for n, k in [(int(a), int(b)) for a, b in (x.split(',') for x in sys.argv[1:])]:
    P = parse_prod('../out/census_n%d_k%d.out' % (n, k))
    Mi = parse_mine('out/census_n%d_k%d.out' % (n, k))
    assert set(P) == set(Mi), 'class lists differ'
    bad = [b for b in P if (P[b][0], P[b][1], P[b][2]) != Mi[b][:3] or (P[b][3] is not None and Mi[b][3] is not None and P[b][3] != Mi[b][3])]
    print('n=%d k=%d: %d classes, %d disagreements with the producer (Smith of M, unglued, spur2/3, kills)' % (n, k, len(P), len(bad)))
    for b in bad[:10]: print('   ', b, P[b], Mi[b])
