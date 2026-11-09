"""AUDIT G1: band certificates (THM-4985 (3)) with own DP (aud_core), S(11) from census T-lines."""
import os, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.normpath(os.path.join(HERE, '..', '..', '..', '..', '04-computation', 'experiments', 'genesis_20261109'))
S11 = set(int(x) for x in open(os.path.join(HERE, 'S_11.txt')))
S10 = set(int(x) for x in open(os.path.join(HERE, 'S_10.txt')))
core = os.path.join(HERE, 'aud_core')
def run(mode, inp, *args):
    return subprocess.run([core, mode] + list(args), input=inp, capture_output=True, text=True, check=True).stdout
def props(n, code):
    line = run('props', '%d %s\n' % (n, code)).strip()
    a, b, c, d = line.split('|')
    t = a.split()
    kv = dict(x.split('=') for x in t[2:])
    return int(t[1]), {k: int(v) for k, v in kv.items()}, list(map(int, b.split())), list(map(int, c.split()[1:]))
# witnesses
wit = {}
for fn in ('wit_band3_13_H.out', 'wit_band4_16_dels.out', 'wit_band4_16_H.out'):
    for l in open(os.path.join(EXP, fn)):
        p = l.split()
        if len(p) == 3 and p[0] == 'W':
            L = len(p[2]); m = int(round((1 + (1 + 8 * L) ** 0.5) / 2)); assert m * (m - 1) // 2 == L
            H, kv, dels, sc = props(m, p[2])
            ok = (H == int(p[1]))
            print('witness %s: order %d claimed %s recomputed %d -> %s' % (fn, m, p[1], H, 'OK' if ok else 'MISMATCH'))
            if ok: wit.setdefault(H, []).append(m)
for s, n in ((3, 13), (4, 16)):
    code = run('band', '', str(s), str(n)).split()[1]
    H, kv, dels, sc = props(n, code)
    print('band_%d(%d): H=%d lambda=%d prime=%d moddef=%d twin=%d scores=%s' % (s, n, H, kv['lam'], kv['prime'], kv['moddef'], kv['twin'], sorted(sc)))
    print('   deletions:', dels)
    print('   distinct deletion values:', sorted(set(dels)), len(set(dels)))
    # membership H in S(n-1)
    cert = [m for m in wit.get(H, []) if m <= n - 1]
    print('   H in S(%d): witness orders %s -> %s' % (n - 1, cert, bool(cert)))
    if n - 2 == 11:
        print('   deletions in S(11) (census):', all(d in S11 for d in dels), ' max deletion', max(dels))
    else:
        dc = {d: [m for m in wit.get(d, []) if m <= n - 2] for d in set(dels)}
        print('   deletion certificates (witness orders <= %d):' % (n - 2), dc, ' all certified:', all(dc.values()))
