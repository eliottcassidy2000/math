# Combine the three verification runs over [957, 1.1e11), compare records with OEIS A031142/A031143,
# and compare the suffix tail with the exact 5-adic model N * Z_k / T_k.
import json, re
W = 288
runs = [('run_957_1e10.out', 'R1', 957, 10**10), ('run_1e10_6e10.out', 'R2', 10**10, 6*10**10), ('run_6e10_11e10.out', 'R3', 6*10**10, 11*10**10)]
H = [0]*(W+1); recs = []; total = 0; surv = 0
for f, tag, a, b in runs:
    txt = open(f).read().splitlines()
    done = [l for l in txt if 'DONE' in l][0]
    m = re.search(r'count=(\d+) survivors=(\d+)', done); total += int(m.group(1)); surv += int(m.group(2))
    h = list(map(int, [l for l in txt if l.startswith(tag + ' HIST')][0].split()[2:]))
    H = [x + y for x, y in zip(H, h)]
    for l in txt:
        if 'RECORD' in l:
            mm = re.search(r'n=(\d+) first_zero_from_right=(\d+)', l); recs.append((int(mm.group(1)), int(mm.group(2)) - 1))
print('exponents checked by C (n in [957, 1.1e11)):', total, ' survivors (no zero in last 288 digits):', surv)
# global records of zeroless-suffix length
recs.sort(); glob = []; best = -1
for n, s in recs:
    if s > best: glob.append((n, s)); best = s
print('global records (n, zeroless suffix length):', glob)
a142 = [int(x.split()[1]) for x in open('b031142.txt') if x.strip() and not x.startswith('#')]
a143 = [1,2,3,4,5,6,8,9,10,11,12,15,16,21,22,23,24,25,26,36,38,54,57,59,93,115,119,120,121,136,138,164,174,176,191,196,212,217,227,233,249,250,260,268,275,308]
oeis = {n: s for n, s in zip(a142, a143)}
inrange = [(n, s) for n, s in oeis.items() if 957 <= n < 11*10**10]
ours = dict(glob)
print('A031142 records in range:', inrange)
print('all reproduced exactly:', all(ours.get(n) == s for n, s in inrange if n != min(n for n, _ in inrange) or True))
Z = {int(k): int(v) for k, v in json.load(open('zk_values.json'))['Z'].items()}
tau = 1.25 * 0.887694043115148264
def trail(k): return Z[k] / (4 * 5**(k-1)) if k in Z else tau * 0.9**k
ge = [sum(H[k:]) for k in range(W+1)]
print('tail #(zeroless suffix >= k): obs vs N*Z_k/T_k')
for k in [20, 40, 60, 100, 140, 180, 200, 220, 230, 240, 245, 249, 250]:
    print(k, ge[k], round(total * trail(k), 2))
json.dump({'total': total, 'survivors': surv, 'records': glob, 'hist': H}, open('verify_combined.json', 'w'))
