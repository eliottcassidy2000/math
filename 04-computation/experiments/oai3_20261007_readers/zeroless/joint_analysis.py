# Adelic independence test: leading (archimedean) vs trailing (5-adic) zeroless lengths over n in [957, 1e10).
import json, math
from mpmath import mpf, nstr
W = 288; PMAX = 36
txt = open('run_957_1e10.out').read().splitlines()
hist = list(map(int, [l for l in txt if l.startswith('R1 HIST')][0].split()[2:]))   # index pos-1, pos=1..289
J = {}
for l in txt:
    if l.startswith('R1 J '):
        p = l.split(); J[int(p[2])] = list(map(int, p[3:]))
N = sum(hist)
print('N =', N)
Z = {int(k): int(v) for k, v in json.load(open('zk_values.json'))['Z'].items()}
mlead = {int(a): float(v) for a, v in json.load(open('leading_measure.json')).items()}
tau = 1.25 * 0.887694043115148264; lam = 1.0845342222539724319
def trail(k): return 1.0 if k == 0 else (Z[k] / (4 * 5**(k-1)) if k in Z else tau * 0.9**k)
def lead(a): return 1.0 if a == 0 else mlead.get(a, lam * 0.9**a)
# suffix tail counts: suffix length s = pos-1 ; count(s >= k) = sum_{pos >= k+1}
suf_ge = [sum(hist[k:]) for k in range(W+1)]       # suf_ge[k] = #(suffix >= k)
print('\nTrailing (5-adic) side: #(zeroless suffix >= k) vs N*Z_k/T_k')
for k in [1,2,3,5,8,10,13,16,20,30,40,60,80,100,120,140,160,180,200,210,215]:
    e = N * trail(k)
    print('k=%3d obs=%12d exp=%14.2f ratio=%.5f' % (k, suf_ge[k], e, suf_ge[k]/e if e else float('nan')))
# prefix tail counts
pre_ge = [sum(sum(J[p]) for p in range(a, PMAX+1)) for a in range(PMAX+1)]
print('\nLeading (archimedean) side: #(zeroless prefix >= a) vs N*m_a')
for a in [1,2,3,5,8,10,15,20,25,30,35]:
    e = N * lead(a)
    print('a=%2d obs=%12d exp=%14.2f ratio=%.5f' % (a, pre_ge[a], e, pre_ge[a]/e))
# joint
def joint_ge(a, k):
    return sum(sum(J[p][k:]) for p in range(a, PMAX+1))
print('\nJoint: #(prefix >= a and suffix >= k) vs N*m_a*Z_k/T_k (independence)')
worst = 0
for a in [2,5,10,15,20,25,30]:
    row = []
    for k in [2,10,20,40,60,80,100,120,140,160]:
        o = joint_ge(a, k); e = N * lead(a) * trail(k)
        row.append('%d/%.1f' % (o, e))
        if e > 50: worst = max(worst, abs(o - e) / math.sqrt(e))
    print('a=%2d ' % a + '  '.join(row))
print('max |obs-exp|/sqrt(exp) over cells with exp>50: %.2f' % worst)
# distribution of total coverage s = prefix + suffix (prefix capped at 36)
cov = {}
for p in range(PMAX+1):
    for pos in range(1, W+2):
        s = p + (pos - 1)
        cov[s] = cov.get(s, 0) + J[p][pos-1]
smax = max(s for s in cov if cov[s] > 0)
print('\nmax total zeroless coverage (prefix+suffix) over n in [957,1e10):', smax)
# expected #(prefix+suffix >= s) under independence: sum over a of P(pre = a) P(suf >= s-a)
def P_pre_eq(a): return lead(a) - lead(a+1) if a < PMAX else lead(PMAX)
for s in [100, 150, 200, 220, 230, 240]:
    obs = sum(v for t, v in cov.items() if t >= s)
    exp = N * sum(P_pre_eq(a) * trail(max(s - a, 0)) for a in range(0, PMAX+1))
    print('#(coverage >= %d): obs=%d exp=%.2f' % (s, obs, exp))
