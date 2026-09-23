#!/usr/bin/env python3
"""collatz_procgen_20260923_sweep_records.py  (HYP-9122 / HYP-9125 record objects)

Drives the sparse tight DP (collatz_procgen_20260923_sweep_tightdp.c) over the record
denominators of log_2 3, reconstructs every record object, verifies each one with an
independent exact checker, checks the hub conditions of the splicing theorems
(loops note Theorem 4.2 / mirror note Theorem 2.2), and builds + verifies explicit spliced
loops for sample lengths.

Usage:
  python3 collatz_procgen_20260923_sweep_records.py SIDE QMAX [--recon] [--samples s1,s2,...]
     SIDE = pos | neg ; QMAX = largest record to treat
     binary: env TIGHTDP (default scratch/procgen_sweep/tightdp_new), checkpoints in env CKDIR
Objects are written to env OBJDIR (default scratch/procgen_sweep/objects) as text files.
"""
import os, sys, subprocess, time, json
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
EXE = os.environ.get('TIGHTDP', os.path.join(ROOT, 'scratch/procgen_sweep/tightdp_new'))
CKDIR = os.environ.get('CKDIR', os.path.join(ROOT, 'scratch/procgen_sweep/ck'))
OBJDIR = os.environ.get('OBJDIR', os.path.join(ROOT, 'scratch/procgen_sweep/objects'))
LN2, LN3 = 0.6931471805599453, 1.0986122886681098
FB = float(os.environ.get('FB', '2.0'))   # base-loop cap = FB * (height-law estimate)
FC = float(os.environ.get('FC', '3.0'))   # cycle cap = FC * layers / G_final
RT = int(os.environ.get('RT', '64'))       # checkpoint spacing

# ---------------- exact arithmetic ----------------
def ceil_qL(q):          # least p with 2^p > 3^q
    return (3**q).bit_length()
def ceil_ntheta(n):      # least a with 3^a > 2^n
    a = (n * 6309) // 10000
    while 3**a <= 2**n: a += 1
    while a > 0 and 3**(a-1) > 2**n: a -= 1
    return a
def eps_rec(qmax):       # upper best approximations: successive minima of ceil(qL)-qL
    from mpmath import mp, mpf, log, ceil
    mp.dps = 60; L = log(3)/log(2); out = []; best = mpf(2)
    for q in range(1, qmax+1):
        e = ceil(q*L) - q*L
        if e < best: best = e; out.append((q, ceil_qL(q)))
    return out
def eta_rec(nmax):       # lower best approximations: successive minima of ceil(n th)-n th
    from mpmath import mp, mpf, log, ceil
    mp.dps = 60; th = log(2)/log(3); out = []; best = mpf(2)
    for n in range(1, nmax+1):
        e = ceil(n*th) - n*th
        if e < best: best = e; out.append((n, ceil_ntheta(n)))
    return out
def y(n): return (3**(n+1)-1)//2     # positive climb points 1,4,13,40,...
def w(t): return (3**t+1)//2         # negative climb points (in u=-v) 1,2,5,14,41,...

def G_pos_loop(q):
    p = ceil_qL(q); return Fraction(2**p - (3**q - 1), 2**p)
def G_pos_cycle(q, h):
    p = ceil_qL(q); return 1 - Fraction(3**q*(2*h+1) - 1, 2**(p+1)*h)
def G_neg_loop(n):
    a = ceil_ntheta(n); return Fraction(3**a + 1, 2**n) - 1
def G_neg_cycle(n, h):
    a = ceil_ntheta(n); return Fraction(3**a*(2*h-1) + 1, 2**(n+1)*h) - 1

# ---------------- independent verifier ----------------
def verify(sign, h, moves, K_target):
    """Exact check of a reverse-move object from h back to h.
    sign=+1: x -> (2^k x - 1)/3 ; sign=-1: u -> (2^k u + 1)/3 (u = -v, negative E-graph).
    Checks: legality of every move, return to h, total K, forward E-simulation, carry identity."""
    x = h; K = 0; pts = [h]
    for k in moves:
        t = (x << k) + (-1 if sign > 0 else 1)
        if t % 3: return False, "non-integer"
        yv = t // 3
        if yv <= 0 or yv % 3 == 0: return False, "illegal"
        x = yv; K += k; pts.append(x)
    if x != h: return False, "does not close"
    if K != K_target: return False, "K mismatch"
    s = len(moves)
    v = h if sign > 0 else -h                 # forward E-simulation (signed values)
    for idx in range(s, 0, -1):
        v = 3*v + 1
        for _ in range(moves[idx-1]):
            if v % 2: return False, "halving an odd number"
            v //= 2
        if v != (pts[idx-1] if sign > 0 else -pts[idx-1]): return False, "forward mismatch"
    Ki = 0; B = 0
    for i, k in enumerate(moves, start=1):
        Ki += k; B += 3**(i-1) << (K - Ki)
    ok = (2**K*h == 3**s*h + B) if sign > 0 else (3**s*h == 2**K*h + B)
    if not ok: return False, "carry identity fails"
    return True, dict(s=s, K=K, height=max(pts), pts=pts)

# ---------------- C driver ----------------
def run(sign, h, layers, K, G, xmax, climb=0, recon=0):
    os.makedirs(CKDIR, exist_ok=True)
    args = [EXE, str(sign), str(h), str(layers), str(K), "%.17e" % float(G), str(xmax),
            str(climb), str(recon), CKDIR]
    r = subprocess.run(args, capture_output=True, text=True)
    d = {'raw': r.stdout.strip(), 'err': r.stderr.strip()[-300:]}
    for line in r.stdout.strip().split('\n'):
        if line.startswith('moves='): d['moves'] = [int(t) for t in line[6:].split(',')]
        elif line.startswith('sign='):
            for kv in line.split():
                a, b = kv.split('='); d[a] = b
    return d

def objfile(name): os.makedirs(OBJDIR, exist_ok=True); return os.path.join(OBJDIR, name)

def treat(side, recs, recon, log):
    objs = {}
    sign = 1 if side == 'pos' else -1
    for (r, c) in recs:                        # pos: (q,p) ; neg: (n,a)
        if side == 'pos':
            q, p = r, c
            if q < 3: continue                 # C_1 = trivial loop at 1 (moves [2]) is added by hand
            e = float(p - q*Fraction(15849625007211561814537389439478, 10**31))
            layersB, KB, GB = q-1, p-1, G_pos_loop(q)
            xmaxB = int(max(10**7, FB*1.3*q/(LN2*e)))
            hubs = [y(n) for n in range(0, 20)]
            cyc = lambda h: (q, p, G_pos_cycle(q, h))
            ctrip = lambda h: FC*q
        else:
            n, a = r, c
            if n < 3: continue                 # C_1 = basic loop MH at -1 (moves [1]) is added by hand
            e = float(a - n*Fraction(6309297535714574370995271143427, 10**31))
            layersB, KB, GB = a, n-1, G_neg_loop(n)
            xmaxB = int(max(10**7, FB*1.3*a/(LN3*e)))
            hubs = [w(t) for t in range(0, 22)]
            cyc = lambda h: (a, n, G_neg_cycle(n, h))
            ctrip = lambda h: FC*a
        t0 = time.time()
        needB = (side == 'pos') or (r >= 11)   # mirror: base loops only for records n >= 11 (Theorem 2.2)
        dB = run(sign, 1, layersB, KB, GB, xmaxB, 1, RT if recon else 0) if needB else {}
        rec = dict(record=r, B_found=dB.get('found'), B_climb=dB.get('climb'), B_maxstates=dB.get('maxstates'),
                   B_xmax=xmaxB, B_time=round(time.time()-t0, 1))
        if recon and 'moves' in dB:
            ok, info = verify(sign, 1, dB['moves'], KB)
            rec['B_verified'] = ok
            if ok:
                rec['B_height'] = info['height']
                cp = [i for i, hv in enumerate(hubs) if hv in set(info['pts'])]
                rec['B_hubs_visited_max'] = max(cp) if cp else -1
                with open(objfile('%s_B_%d.txt' % (side, r)), 'w') as f: f.write(','.join(map(str, dB['moves'])))
        # cycle: lowest climb hub that works
        for i, h in enumerate(hubs):
            L, P, G = cyc(h)
            if G <= 0 or float(G) < 0.25*e*(LN2 if side == 'pos' else LN3): continue
            xmaxC = int(max(10**6, ctrip(h)/float(G)))
            t0 = time.time()
            dC = run(sign, h, L, P, G, xmaxC, 0, RT if recon else 0)
            if dC.get('found') == '1':
                rec.update(C_hub_index=i, C_hub=h, C_maxstates=dC.get('maxstates'), C_xmax=xmaxC, C_time=round(time.time()-t0, 1))
                if recon and 'moves' in dC:
                    ok, info = verify(sign, h, dC['moves'], P)
                    rec['C_verified'] = ok
                    if ok:
                        rec['C_height'] = info['height']
                        with open(objfile('%s_C_%d_%d.txt' % (side, r, h)), 'w') as f: f.write(','.join(map(str, dC['moves'])))
                break
        objs[r] = rec
        print(json.dumps(rec), flush=True); log.write(json.dumps(rec) + '\n'); log.flush()
    return objs

if __name__ == '__main__':
    side = sys.argv[1]; qmax = int(sys.argv[2]); recon = '--recon' in sys.argv
    qmin = 1
    for a in sys.argv:
        if a.startswith('--qmin='): qmin = int(a.split('=')[1])
    recs = eps_rec(qmax) if side == 'pos' else eta_rec(qmax)
    recs = [rc for rc in recs if rc[0] >= qmin]
    with open(objfile('%s_log.jsonl' % side), 'a') as log:
        treat(side, recs, recon, log)
